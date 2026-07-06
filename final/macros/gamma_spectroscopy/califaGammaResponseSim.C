// CALIFA gamma-response simulation using R3BCALIFATestGenerator
// with INCL-style vertex smearing conventions:
//   SetXYZ(x,y,z)       -> central vertex [cm]
//   SetDxDyDz(sx,sy,sz) -> x,y Gaussian sigmas [cm], z uniform in [z-sz, z+sz]
//
// -------------------------------------------------------------------------
// SIMULATION MODES for fitting the 22O spectrum:
//
//  1. califaGammaResponseSimSingle(3199., ...)
//       --> single gamma line, e.g. direct 2+ -> 0+ (3199 keV)
//
//  2. califaGammaResponseSimCascade({e1, e2, ...}, ...)
//       --> any cascade of N gammas emitted in the same event from the
//           same nucleus (same vertex, beta, beam direction).
//           All lines have BR=1; sumBR=N>1 so ALL are fired every event.
//
// -------------------------------------------------------------------------
// RECOMMENDED WORKFLOW for the multiplicity-resolved population fit
// (subtract_background.C). beta = 0.79 for the 25F(p,2p) -> 22O channel;
// use the measured dispersion of beta_frag as betaDispersion.
//
//   .L califaGammaResponseSim.C+
//
//   // >= 1e6 events each so template statistics never limit the fit
//   califaGammaResponseSimSingle (3199.,          1000000, 0.79, 0.006);
//   califaGammaResponseSimCascade({1383., 3199.}, 1000000, 0.79, 0.006);
//
//   // Multiplicity-resolved, PER-DECAY-normalisable templates
//   // (binning, cluster acceptance and multiplicity definition are taken
//   //  from GammaCfg in gammaSpectra.h — identical to the data macro):
//   buildAll22OResponseSlices(0.79);
//
// IMPORTANT — two fixes w.r.t. the old examples:
//   * The cascade template must contain the FULL per-cluster spectrum
//     (BOTH the 1383 and the 3199 lines). The old example applied an
//     energy window 2500-4000 keV, which amputated the 1383 peak and made
//     the cascade template almost identical to the direct-3199 one — that
//     is precisely the template degeneracy seen in the 1D fits.
//     buildGammaResponseSlices applies NO energy window.
//   * Templates are written with the SAME binning as the data histograms
//     (GammaCfg::NBINS_E over [E_MIN_MEV, E_MAX_MEV]) so the fit can be a
//     bin-by-bin Poisson likelihood. The old 1600-bin/0-5.76 MeV binning
//     did not match the data.
//
// The legacy builders buildGammaResponseHist / buildGammaResponseHistPerCluster
// are kept unchanged below for backwards compatibility, but the population
// fit uses ONLY the output of buildGammaResponseSlices.
// -------------------------------------------------------------------------

#include "gammaSpectra.h"   // GammaCfg: shared binning / acceptance / names

#include <TClonesArray.h>
#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TGraph.h>
#include <TMath.h>
#include <TStopwatch.h>
#include <TString.h>
#include <TSystem.h>
#include <TTree.h>
#include <TRandom3.h>
#include <TVector3.h>
#include <TROOT.h>
#include <TParameter.h>
#include <algorithm>
#include <numeric>
#include <sstream>
#include <iomanip>
#include <vector>

namespace
{
    void ConfigureEnvironment()
    {
        TString dir = gSystem->Getenv("VMCWORKDIR");
        TString r3b_geomdir = dir + "/geometry/";
        gSystem->Setenv("GEOMPATH", r3b_geomdir.Data());
        r3b_geomdir.ReplaceAll("//", "/");
        TString r3b_confdir = dir + "/gconfig/";
        gSystem->Setenv("CONFIG_DIR", r3b_confdir.Data());
        r3b_confdir.ReplaceAll("//", "/");
    }

    Double_t DopplerCorrectEnergy(const Double_t eLabMeV,
                                  const Double_t beta,
                                  const Double_t thetaLabRad)
    {
        const Double_t gamma = 1.0 / TMath::Sqrt(1.0 - beta * beta);
        return gamma * (1.0 - beta * TMath::Cos(thetaLabRad)) * eLabMeV;
    }

    // -----------------------------------------------------------------------
    // Recompute the polar angle of a CALIFA cluster as seen from the true MC
    // emission vertex, rather than the (fixed) origin used by the standard
    // R3BCalifaClusterData::GetTheta().
    //
    // The cluster 3D position is reconstructed from (theta_fixed, phi_cluster)
    // assuming the cluster centroid sits at a fixed CALIFA radius R from the
    // origin. For the CALIFA barrel R ~ 30 cm is a reasonable approximation.
    // The vertex offset is then subtracted and the new polar angle is taken.
    // -----------------------------------------------------------------------
    Double_t ClusterThetaWrtVertex(Double_t clusterTheta,
                                   Double_t clusterPhi,
                                   Double_t vx,
                                   Double_t vy,
                                   Double_t vz,
                                   Double_t R_califa_cm = 30.0)
    {
        const Double_t Cx = R_califa_cm * TMath::Sin(clusterTheta) * TMath::Cos(clusterPhi);
        const Double_t Cy = R_califa_cm * TMath::Sin(clusterTheta) * TMath::Sin(clusterPhi);
        const Double_t Cz = R_califa_cm * TMath::Cos(clusterTheta);

        const Double_t dx = Cx - vx;
        const Double_t dy = Cy - vy;
        const Double_t dz = Cz - vz;
        const Double_t r  = TMath::Sqrt(dx * dx + dy * dy + dz * dz);
        if (r <= 0.) return clusterTheta;
        return TMath::ACos(dz / r);
    }

    // -----------------------------------------------------------------------
    // Locate the primary MC vertex of the current event.
    // Returns true and fills (vx, vy, vz) if the first MotherId == -1 track
    // is found; returns false otherwise (then caller should fall back to origin).
    // -----------------------------------------------------------------------
    bool GetPrimaryMCVertex(TClonesArray* mcTracks,
                            Double_t& vx, Double_t& vy, Double_t& vz)
    {
        if (!mcTracks) return false;
        for (Int_t k = 0; k < mcTracks->GetEntriesFast(); ++k)
        {
            auto* mc = dynamic_cast<R3BMCTrack*>(mcTracks->At(k));
            if (mc && mc->GetMotherId() == -1 && mc->GetPdgCode() == 22)
            {
                vx = mc->GetStartX();
                vy = mc->GetStartY();
                vz = mc->GetStartZ();
                return true;
            }
        }
        return false;
    }

    // Build a filename tag from the cascade energies, e.g. "1383+3199" or "2354+1383+3199"
    TString CascadeTag(const std::vector<Double_t>& energies_keV)
    {
        TString tag = "";
        for (size_t i = 0; i < energies_keV.size(); ++i)
        {
            if (i > 0) tag += "+";
            tag += Form("%d", (int)TMath::Nint(energies_keV[i]));
        }
        return tag;
    }

    FairRunSim* BuildRun(const TString& outFile,
                         const TString& parFile,
                         Bool_t includeTargetArea)
    {
        FairRunSim* run = new FairRunSim();
        run->SetName("TGeant4");
        run->SetOutputFile(outFile.Data());
        FairRuntimeDb* rtdb = run->GetRuntimeDb();

        FairParAsciiFileIo* parIo1 = new FairParAsciiFileIo();
        parIo1->open("califaDigi.par", "in");
        rtdb->setFirstInput(parIo1);
        rtdb->print();
        rtdb->getContainer("califaCrystalPars4Sim");
        UInt_t runId = 1;
        rtdb->initContainers(runId);

        run->SetMaterials("media_r3b.geo");

        FairModule* cave = new R3BCave("CAVE");
        cave->SetGeometryFileName("r3b_cave.geo");
        run->AddModule(cave);

        if (includeTargetArea)
        {
            auto* foots = new R3BTra("target_area_v2025_5cm.geo.root");
            run->AddModule(foots);
        }

        auto* calsim = new R3BCalifa("califa_v2025.6.geo.root");
        calsim->SelectGeometryVersion(2025);
        run->AddModule(calsim);

        auto* califaDig = new R3BCalifaDigitizer();
        califaDig->SetNonUniformity(1.);
        califaDig->SetRealConfig(false);
        califaDig->SetExpGammaEnergyRes(12.); // 5. means 5% at 1 MeV
        califaDig->SetComponentRes(6.);
        califaDig->SetDetectionThreshold(0.000010); // in GeV!! 0.000010 means 10 keV
        run->AddTask(califaDig);

        auto* cal2Clus = new R3BCalifaCrystalCal2Cluster();
        cal2Clus->SetCrystalThreshold(0.1);
        cal2Clus->SetProtonClusterThreshold(18.);
        cal2Clus->SetGammaClusterThreshold(0.1);
        cal2Clus->SelectGeometryVersion(2025);
        cal2Clus->SetRoundWindow(0.4);
        cal2Clus->SetRandomization(kTRUE);
        TString randomization_file = "/home/e12exp/ssd/R3BParams_g249/califa/gammas_minus_4.root";
        randomization_file.ReplaceAll("//", "/");
        cal2Clus->SetRandomizationFile(randomization_file);
        run->AddTask(cal2Clus);

        FairLogger::GetLogger()->SetLogVerbosityLevel("LOW");
        FairLogger::GetLogger()->SetLogScreenLevel("INFO");

        return run;
    }

    void FinaliseRun(FairRunSim* run, const TString& parFile, Int_t nEvents, Int_t randomSeed)
    {
        run->Init();
        TVirtualMC::GetMC()->SetRandom(new TRandom3(randomSeed));
        TVirtualMC::GetMC()->SetMaxNStep(-15000);

        Bool_t kParameterMerged = kTRUE;
        auto* parOut = new FairParRootFileIo(kParameterMerged);
        FairRuntimeDb* rtdb = run->GetRuntimeDb();
        parOut->open(parFile.Data());
        rtdb->setOutput(parOut);
        rtdb->saveOutput();

        if (nEvents > 0)
            run->Run(nEvents);
    }
} // anonymous namespace

// ===========================================================================
// 1. Single gamma-line simulation
//    e.g. direct 2+ -> 0+ at 3199 keV
//    Optional beam angular divergence (1-sigma projected angles in mrad).
// ===========================================================================
void califaGammaResponseSimSingle(Double_t eGamma_keV = 3199.,
                                  Int_t nEvents = 1000000,
                                  Double_t beta = 0.79,
                                  Double_t betaDispersion = 0.0,
                                  Double_t sigmaX_cm = 0.1,
                                  Double_t sigmaY_cm = 0.1,
                                  Double_t halfThicknessZ_cm = 2.5,
                                  Double_t zCenter_cm = 4.0,
                                  Bool_t includeTargetArea = kTRUE,
                                  Int_t randomSeed = 0,
                                  TString outFile = "",
                                  TString parFile = "",
                                  Double_t beamDivX_mrad = 0.,
                                  Double_t beamDivY_mrad = 0.)
{
    ConfigureEnvironment();
    gROOT->SetBatch(kTRUE);

    if (outFile.IsNull())
        outFile = Form("./sim_gamma_%dkeV.root", (int)TMath::Nint(eGamma_keV));
    if (parFile.IsNull())
        parFile = Form("./par_gamma_%dkeV.root", (int)TMath::Nint(eGamma_keV));

    if (beta <= 0.0 || beta >= 1.0)
    {
        Error("califaGammaResponseSimSingle", "beta must be in (0,1), got %.6f", beta);
        return;
    }

    TStopwatch timer;
    timer.Start();

    FairRunSim* run = BuildRun(outFile, parFile, includeTargetArea);

    auto* gammaGen = new R3BCALIFATestGenerator(22, 1);
    gammaGen->SetPRange(eGamma_keV * 1.e-6, eGamma_keV * 1.e-6); // keV -> GeV
    gammaGen->SetThetaRange(0., 180.);
    gammaGen->SetCosTheta();
    gammaGen->SetPhiRange(0., 360.);
    gammaGen->SetXYZ(0., 0., zCenter_cm);
    gammaGen->SetDxDyDz(sigmaX_cm, sigmaY_cm, halfThicknessZ_cm);
    gammaGen->SetLorentzBoost(beta, betaDispersion);
    if (beamDivX_mrad > 0. || beamDivY_mrad > 0.)
        gammaGen->SetBeamAngularDivergence(beamDivX_mrad, beamDivY_mrad);

    auto* primGen = new FairPrimaryGenerator();
    primGen->AddGenerator(gammaGen);
    run->SetGenerator(primGen);
    run->SetStoreTraj(kFALSE);

    FinaliseRun(run, parFile, nEvents, randomSeed);

    timer.Stop();
    std::cout << "\nFinished single-line simulation: " << eGamma_keV << " keV\n"
              << "  beta      = " << beta << " +/- " << betaDispersion << "\n"
              << "  output    = " << outFile << "\n"
              << "  real time = " << timer.RealTime() << " s\n\n";
}

// ===========================================================================
// 2. Cascade simulation — works for ANY number of gammas (double, triple, ...)
//
//    Pass the energies in decay order, e.g.:
//      {1383., 3199.}          double cascade
//      {2354., 1383., 3199.}   triple cascade
//
//    All gammas share the same vertex, beta, and beam direction per event
//    (correct for sequential emission from one moving nucleus).
//    BR=1.0 for every line; sumBR=N>1 so ALL lines fire every event.
//
//    Output file is named automatically, e.g.:
//      sim_cascade_1383+3199keV.root
// ===========================================================================
void califaGammaResponseSimCascade(std::vector<Double_t> energies_keV,
                                   Int_t nEvents = 1000000,
                                   Double_t beta = 0.79,
                                   Double_t betaDispersion = 0.0,
                                   Double_t sigmaX_cm = 0.1,
                                   Double_t sigmaY_cm = 0.1,
                                   Double_t halfThicknessZ_cm = 2.5,
                                   Double_t zCenter_cm = 4.0,
                                   Bool_t includeTargetArea = kTRUE,
                                   Int_t randomSeed = 0,
                                   TString outFile = "",
                                   TString parFile = "",
                                   Double_t beamDivX_mrad = 0.,
                                   Double_t beamDivY_mrad = 0.)
{
    ConfigureEnvironment();
    gROOT->SetBatch(kTRUE);

    if (energies_keV.empty())
    {
        Error("califaGammaResponseSimCascade", "No gamma energies provided.");
        return;
    }
    if (beta <= 0.0 || beta >= 1.0)
    {
        Error("califaGammaResponseSimCascade", "beta must be in (0,1), got %.6f", beta);
        return;
    }

    const TString tag = CascadeTag(energies_keV);
    if (outFile.IsNull()) outFile = Form("./sim_cascade_%skeV.root", tag.Data());
    if (parFile.IsNull()) parFile = Form("./par_cascade_%skeV.root", tag.Data());

    TStopwatch timer;
    timer.Start();

    FairRunSim* run = BuildRun(outFile, parFile, includeTargetArea);

    // One generator, nuclear decay chain.
    // BR=1 for every line => sumBR = N > 1 => all lines fire every event.
    //
    // NOTE: SetDecayChainPoint expects energy in MeV (it internally divides
    // by 1000 to convert to GeV). Our vector is in keV, so divide by 1000.
    auto* cascadeGen = new R3BCALIFATestGenerator(22, 1);
    cascadeGen->SetNuclearDecayChain();
    for (auto e : energies_keV)
        cascadeGen->SetDecayChainPoint(e / 1000., 1.0); // keV -> MeV

    cascadeGen->SetThetaRange(0., 180.);
    cascadeGen->SetCosTheta();
    cascadeGen->SetPhiRange(0., 360.);
    cascadeGen->SetXYZ(0., 0., zCenter_cm);
    cascadeGen->SetDxDyDz(sigmaX_cm, sigmaY_cm, halfThicknessZ_cm);
    cascadeGen->SetLorentzBoost(beta, betaDispersion);
    if (beamDivX_mrad > 0. || beamDivY_mrad > 0.)
        cascadeGen->SetBeamAngularDivergence(beamDivX_mrad, beamDivY_mrad);

    auto* primGen = new FairPrimaryGenerator();
    primGen->AddGenerator(cascadeGen);
    run->SetGenerator(primGen);
    run->SetStoreTraj(kFALSE);

    FinaliseRun(run, parFile, nEvents, randomSeed);

    timer.Stop();
    std::cout << "\nFinished cascade simulation: " << tag << " keV\n";
    for (auto e : energies_keV)
        std::cout << "  gamma: " << e << " keV\n";
    std::cout << "  beta      = " << beta << " +/- " << betaDispersion << "\n"
              << "  output    = " << outFile << "\n"
              << "  real time = " << timer.RealTime() << " s\n\n";
}

// ===========================================================================
// buildGammaResponseSlices  — THE template builder for the population fit.
//
// For each simulated DECAY (tree entry) of a given 22O state:
//   * collect clusters passing GammaCfg::PassCluster on the LAB energy
//     (identical acceptance to the data macro; multiplicity M := their count)
//   * Doppler-correct each accepted cluster (nominal beta, theta w.r.t. z,
//     same approximation as the data procedure with measured beta/direction)
//   * fill each accepted cluster into the all / M==m / M<=m / M>=m spectra
//     with EXACTLY the data binning (GammaCfg::NBINS_E over [0,10] MeV)
//
// Written to outRoot:
//   hResp_<state>_all, hResp_<state>_mult1..MULT_MAX,
//   hResp_<state>_mleq1..MULT_MAX, hResp_<state>_mgeq2..MULT_MAX  [raw counts]
//   hRespMult_<state>                                             [M distribution]
//   h2Resp_E1_vs_E2_<state>, hResp_<state>_gate1383/_gate3199     [QA]
//   TParameter<double> nDecays_<state>  = number of SIMULATED DECAYS
//                                         (INCLUDING zero-cluster events!)
//
// The fit divides every slice histogram of a state by the SAME nDecays, so a
// fitted amplitude is directly "number of decays of that state in the data".
// Do NOT unit-normalise the slices individually — the relative slice
// populations carry the multiplicity information that breaks the
// direct-vs-cascade degeneracy.
//
// energyScaleToMeV: applied to R3BCalifaClusterData::GetEnergy(). In this
// setup GetEnergy() returns MeV (scale 1.0). If your branch is in keV, pass
// 1.e-3; a loud sanity check prints the mean cluster energy at start-up.
// ===========================================================================
void buildGammaResponseSlices(TString simFile,
                              TString stateTag,   // GammaCfg::StateDirect() / StateCascade()
                              Double_t betaForCorrection,
                              TString outRoot = "./gamma_response_functions.root",
                              Bool_t recreate = kFALSE,
                              Double_t energyScaleToMeV = 1.0,
                              // ---- extra broadening to match data peak width ----
                              // Applied per corrected cluster as a Gaussian smear of
                              //   sigma(E) = sqrt( (extraResFrac*E)^2 + extraSigAbs^2 )
                              // extraResFrac is the energy-PROPORTIONAL term (fixes both
                              // the 1383 and 3199 peaks with ONE number); extraSigAbs is
                              // an optional constant floor. Both default to 0 (no smear).
                              Double_t extraResFrac = 0.0,
                              Double_t extraSigAbs  = 0.0,
                              TString treeName = "evt",
                              TString branchName = "CalifaClusterData")
{
    using namespace GammaCfg;

    if (extraResFrac > 0. || extraSigAbs > 0.)
        std::cout << "[buildGammaResponseSlices] extra broadening: sigma(E) = sqrt(("
                  << extraResFrac << "*E)^2 + " << extraSigAbs << "^2) MeV\n";

    auto* inFile = TFile::Open(simFile, "READ");
    if (!inFile || inFile->IsZombie())
    {
        Error("buildGammaResponseSlices", "Cannot open %s", simFile.Data());
        return;
    }
    TTree* tree = dynamic_cast<TTree*>(inFile->Get(treeName));
    if (!tree) { Error("buildGammaResponseSlices", "Tree '%s' not found", treeName.Data());
                 inFile->Close(); return; }

    auto* clusters = new TClonesArray("R3BCalifaClusterData");
    tree->SetBranchAddress(branchName.Data(), &clusters);

    const Long64_t nDecays = tree->GetEntries();
    std::cout << "[buildGammaResponseSlices] state=" << stateTag
              << "  file=" << simFile << "  nDecays=" << nDecays << "\n";

    // ---- Book histograms (data binning!) -----------------------------------
    auto bookE = [&](const TString& slice) {
        auto* h = new TH1D(RespHist(stateTag, slice),
                           Form("CALIFA response %s, %s;E_{#gamma}^{DC} [MeV];Counts per cluster",
                                stateTag.Data(), slice.Data()),
                           NBINS_E, E_MIN_MEV, E_MAX_MEV);
        h->Sumw2();
        h->SetDirectory(nullptr);
        return h;
    };

    TH1D* hAll = bookE(SliceAll());
    TH1D* hExact[MULT_MAX];
    TH1D* hLeq  [MULT_MAX];
    for (int m = 1; m <= MULT_MAX; ++m) {
        hExact[m-1] = bookE(SliceExact(m));
        hLeq  [m-1] = bookE(SliceLeq(m));
    }
    TH1D* hGeq[MULT_MAX + 1] = { nullptr };
    for (int m = 2; m <= MULT_MAX; ++m)
        hGeq[m] = bookE(SliceGeq(m));

    auto* hMult = new TH1D(RespMultHist(stateTag),
                           Form("M_{#gamma} per decay, %s;M_{#gamma};Decays", stateTag.Data()),
                           16, -0.5, 15.5);
    hMult->Sumw2();
    hMult->SetDirectory(nullptr);

    auto* hGate1383 = new TH1D("hResp_" + stateTag + "_gate1383",
                               "Companions of 1383-gate clusters;E_{#gamma}^{DC} [MeV];Counts",
                               NBINS_E, E_MIN_MEV, E_MAX_MEV);
    auto* hGate3199 = new TH1D("hResp_" + stateTag + "_gate3199",
                               "Companions of 3199-gate clusters;E_{#gamma}^{DC} [MeV];Counts",
                               NBINS_E, E_MIN_MEV, E_MAX_MEV);
    auto* h2E1E2 = new TH2D("h2Resp_E1_vs_E2_" + stateTag,
                            "M==2;E_{low} [MeV];E_{high} [MeV]",
                            NBINS_E, E_MIN_MEV, E_MAX_MEV,
                            NBINS_E, E_MIN_MEV, E_MAX_MEV);
    hGate1383->Sumw2(); hGate1383->SetDirectory(nullptr);
    hGate3199->Sumw2(); hGate3199->SetDirectory(nullptr);
    h2E1E2->Sumw2();    h2E1E2->SetDirectory(nullptr);

    // ---- Units sanity check on the first events -----------------------------
    {
        double sumE = 0.; long nE = 0;
        const Long64_t nCheck = std::min<Long64_t>(nDecays, 2000);
        for (Long64_t i = 0; i < nCheck; ++i) {
            tree->GetEntry(i);
            for (Int_t j = 0; j < clusters->GetEntriesFast(); ++j) {
                auto* cl = dynamic_cast<R3BCalifaClusterData*>(clusters->At(j));
                if (!cl) continue;
                sumE += cl->GetEnergy() * energyScaleToMeV;
                ++nE;
            }
        }
        if (nE > 0) {
            const double meanE = sumE / nE;
            std::cout << "[buildGammaResponseSlices] mean cluster energy (scaled): "
                      << meanE << " MeV\n";
            if (meanE > 50.)
                std::cout << "  *** WARNING: mean cluster energy > 50 MeV for a few-MeV line.\n"
                          << "  *** GetEnergy() is probably in keV in your build —\n"
                          << "  *** rerun with energyScaleToMeV = 1.e-3\n";
        }
    }

    // ---- Event loop ----------------------------------------------------------
    std::vector<double> eCorr;
    for (Long64_t i = 0; i < nDecays; ++i)
    {
        tree->GetEntry(i);
        if (i % 200000 == 0)
            std::cout << "\r  event " << i << " / " << nDecays << std::flush;

        eCorr.clear();
        for (Int_t j = 0; j < clusters->GetEntriesFast(); ++j)
        {
            auto* cl = dynamic_cast<R3BCalifaClusterData*>(clusters->At(j));
            if (!cl) continue;

            const Double_t eLab = cl->GetEnergy() * energyScaleToMeV;
            if (!PassCluster(eLab)) continue;   // SAME acceptance as data

            Double_t eDC = DopplerCorrectEnergy(eLab, betaForCorrection, cl->GetTheta());

            // Optional extra broadening so the simulated peak width matches data.
            if (extraResFrac > 0. || extraSigAbs > 0.) {
                const Double_t sig = TMath::Sqrt(extraResFrac * eDC * extraResFrac * eDC
                                                 + extraSigAbs * extraSigAbs);
                if (sig > 0.) eDC += gRandom->Gaus(0., sig);
            }
            eCorr.push_back(eDC);
        }

        const int M = static_cast<int>(eCorr.size());
        hMult->Fill(static_cast<double>(M));
        if (M == 0) continue;    // still counted in nDecays (normalisation!)

        for (int ig = 0; ig < M; ++ig) {
            const double e = eCorr[ig];
            hAll->Fill(e);
            if (M <= MULT_MAX) hExact[M-1]->Fill(e);
            for (int m = 1; m <= MULT_MAX; ++m)
                if (M <= m) hLeq[m-1]->Fill(e);
            for (int m = 2; m <= MULT_MAX; ++m)
                if (M >= m) hGeq[m]->Fill(e);
        }

        if (M >= 2) {
            for (int j = 0; j < M; ++j)
                for (int k = 0; k < M; ++k) {
                    if (k == j) continue;
                    if (In1383Gate(eCorr[j])) hGate1383->Fill(eCorr[k]);
                    if (In3199Gate(eCorr[j])) hGate3199->Fill(eCorr[k]);
                }
        }
        if (M == 2) {
            const double eLo = std::min(eCorr[0], eCorr[1]);
            const double eHi = std::max(eCorr[0], eCorr[1]);
            h2E1E2->Fill(eLo, eHi);
        }
    }
    std::cout << "\n";

    // ---- Write ---------------------------------------------------------------
    auto* outFile = TFile::Open(outRoot, recreate ? "RECREATE" : "UPDATE");
    if (!outFile || outFile->IsZombie()) {
        Error("buildGammaResponseSlices", "Cannot open output %s", outRoot.Data());
        inFile->Close();
        return;
    }
    outFile->cd();

    hAll->Write(hAll->GetName(), TObject::kOverwrite);
    for (int m = 1; m <= MULT_MAX; ++m) {
        hExact[m-1]->Write(hExact[m-1]->GetName(), TObject::kOverwrite);
        hLeq  [m-1]->Write(hLeq  [m-1]->GetName(), TObject::kOverwrite);
    }
    for (int m = 2; m <= MULT_MAX; ++m)
        hGeq[m]->Write(hGeq[m]->GetName(), TObject::kOverwrite);
    hMult->Write(hMult->GetName(), TObject::kOverwrite);
    hGate1383->Write(hGate1383->GetName(), TObject::kOverwrite);
    hGate3199->Write(hGate3199->GetName(), TObject::kOverwrite);
    h2E1E2->Write(h2E1E2->GetName(), TObject::kOverwrite);

    TParameter<double> pN(NDecaysName(stateTag), static_cast<double>(nDecays));
    pN.Write(pN.GetName(), TObject::kOverwrite);

    outFile->Close();
    inFile->Close();

    std::cout << "[OK] Wrote multiplicity-sliced response '" << stateTag
              << "' (nDecays=" << nDecays << ") to " << outRoot << "\n";
}

// ===========================================================================
// Convenience: build BOTH 22O templates expected by subtract_background.C
// from the default simulation file names.
//
// extraResFrac : energy-PROPORTIONAL broadening added to the templates so that
//                the simulated peak widths match the DATA. Determine it once
//                from a Gaussian fit to the data vs template 3199 peak:
//                   extraResFrac = sqrt( (sig_data/E)^2 - (sig_tmpl/E)^2 )
//                e.g. sig_data=0.302, sig_tmpl=0.166, E=3.199 -> ~0.079.
//                Because it scales with E, ONE value fixes BOTH the 1383 and
//                3199 peaks. Default 0.079 matches the measured G249 widths;
//                pass 0 to disable and see the raw (too-narrow) templates.
// ===========================================================================
void buildAll22OResponseSlices(Double_t beta = 0.79,
                               TString simDir = ".",
                               TString outRoot = "./gamma_response_functions.root",
                               Double_t energyScaleToMeV = 1.0,
                               Double_t extraResFrac = 0.0,
                               Double_t extraSigAbs  = 0.0)
{
    buildGammaResponseSlices(simDir + "/sim_gamma_3199keV.root",
                             GammaCfg::StateDirect(),  beta, outRoot,
                             kTRUE /*recreate*/, energyScaleToMeV,
                             extraResFrac, extraSigAbs);
    buildGammaResponseSlices(simDir + "/sim_cascade_1383+3199keV.root",
                             GammaCfg::StateCascade(), beta, outRoot,
                             kFALSE /*update*/, energyScaleToMeV,
                             extraResFrac, extraSigAbs);
}

// ===========================================================================
//  ARMEL-METHOD single-line response builder.
//
//  For ONE transition (single gamma of energy eGamma_keV), build:
//    * hResp_line_<keV>  : the per-cluster Doppler-corrected response, with the
//                          SAME binning, acceptance and (optional) broadening as
//                          the data, normalised to UNIT INTEGRAL so a fit
//                          amplitude is the number of DETECTED photopeak counts.
//    * eps_<keV>         : TParameter<double> = photopeak efficiency
//                          = (counts in +-window around the line) / (thrown decays)
//
//  These are all subtract_background.C (Armel method) needs: single templates to
//  measure photopeak areas, and eps(E) to convert areas into populations. No
//  multiplicity slicing and no cascade templates are produced here.
//
//  The input simFile is the output of califaGammaResponseSimSingle(eGamma_keV,...).
// ===========================================================================
// ---- Doppler-correction geometry mode (fixes the (0,0,0)-vs-target bug) ----
//   kVtxOrigin  : legacy — angle from (0,0,0). WRONG when target is not at origin;
//                 kept only to reproduce the old (broken) behaviour.
//   kVtxTargetZ : angle from the FIXED target midpoint (0,0,targetZ). This mirrors
//                 the experiment's "target position in a par file" prescription.
//                 RECOMMENDED — corrects the systematic offset with the SAME fixed
//                 geometry the data uses, so sim and data broaden consistently.
//   kVtxTrueMC  : angle from the TRUE per-event MC emission vertex. Removes the
//                 vertex-spread broadening entirely; use ONLY if the experiment
//                 also reconstructs a per-event vertex (event-by-event tracking).
enum VtxMode { kVtxOrigin = 0, kVtxTargetZ = 1, kVtxTrueMC = 2 };

double buildSingleResponse(TString simFile,
                           int eGamma_keV,
                           Double_t betaForCorrection,
                           TString outRoot = "./gamma_response_functions.root",
                           Bool_t recreate = kFALSE,
                           Double_t energyScaleToMeV = 1.0,
                           Double_t extraResFrac = 0.079,  // match sim width to data
                           Double_t extraSigAbs  = 0.0,
                           Int_t    vtxMode      = kVtxTargetZ, // see VtxMode
                           Double_t targetZ_cm   = 4.0,    // target midpoint z [cm]
                           Double_t R_califa_cm  = 30.0,   // CALIFA barrel radius [cm]
                           TString treeName = "evt",
                           TString branchName = "CalifaClusterData",
                           TString mcBranchName = "MCTrack")
{
    using namespace GammaCfg;

    auto* inFile = TFile::Open(simFile, "READ");
    if (!inFile || inFile->IsZombie())
    { Error("buildSingleResponse", "Cannot open %s", simFile.Data()); return -1.; }
    TTree* tree = dynamic_cast<TTree*>(inFile->Get(treeName));
    if (!tree) { Error("buildSingleResponse", "Tree '%s' not found", treeName.Data());
                 inFile->Close(); return -1.; }

    auto* clusters = new TClonesArray("R3BCalifaClusterData");
    tree->SetBranchAddress(branchName.Data(), &clusters);

    // MC tracks only needed for the true-vertex mode.
    TClonesArray* mcTracks = nullptr;
    if (vtxMode == kVtxTrueMC) {
        if (tree->GetBranch(mcBranchName)) {
            mcTracks = new TClonesArray("R3BMCTrack");
            tree->SetBranchAddress(mcBranchName.Data(), &mcTracks);
        } else {
            Warning("buildSingleResponse",
                    "MCTrack branch '%s' not found — falling back to fixed target z.",
                    mcBranchName.Data());
            vtxMode = kVtxTargetZ;
        }
    }

    const char* modeName = (vtxMode == kVtxOrigin) ? "origin(0,0,0) [legacy]"
                         : (vtxMode == kVtxTargetZ) ? "fixed target z" : "true MC vertex";
    std::cout << "[buildSingleResponse] " << eGamma_keV << " keV: Doppler angle from "
              << modeName;
    if (vtxMode == kVtxTargetZ) std::cout << " (z=" << targetZ_cm << " cm, R=" << R_califa_cm << " cm)";
    std::cout << "\n";

    const Long64_t nDecays = tree->GetEntries();

    auto* hResp = new TH1D(RespLineHist(eGamma_keV),
                           Form("Response %d keV;E_{#gamma}^{DC} [MeV];Counts per cluster",
                                eGamma_keV),
                           NBINS_E, E_MIN_MEV, E_MAX_MEV);
    hResp->Sumw2();
    hResp->SetDirectory(nullptr);

    // 2D angular-QA: Doppler-corrected energy vs CALIFA cluster theta (deg).
    // A correct Doppler correction -> a FLAT horizontal band at E = line energy.
    auto* hEvsTh = new TH2D(RespEvsThetaHist(eGamma_keV),
                            Form("Response %d keV: DC energy vs #theta;"
                                 "#theta_{CALIFA} [deg];E_{#gamma}^{DC} [MeV]", eGamma_keV),
                            NBINS_TH, TH_MIN_DEG, TH_MAX_DEG,
                            NBINS_E,  E_MIN_MEV,  E_MAX_MEV);
    hEvsTh->SetDirectory(nullptr);

    for (Long64_t i = 0; i < nDecays; ++i)
    {
        tree->GetEntry(i);

        // Determine the emission vertex used for the angle correction.
        Double_t vx = 0., vy = 0., vz = 0.;
        if (vtxMode == kVtxTargetZ) {
            vz = targetZ_cm;
        } else if (vtxMode == kVtxTrueMC) {
            if (!GetPrimaryMCVertex(mcTracks, vx, vy, vz)) { vx = 0.; vy = 0.; vz = targetZ_cm; }
        }

        for (Int_t j = 0; j < clusters->GetEntriesFast(); ++j)
        {
            auto* cl = dynamic_cast<R3BCalifaClusterData*>(clusters->At(j));
            if (!cl) continue;
            const Double_t eLab = cl->GetEnergy() * energyScaleToMeV;
            if (!PassCluster(eLab)) continue;

            // Angle for the Doppler correction: recompute w.r.t. the chosen vertex
            // unless we are in legacy origin mode.
            Double_t thetaCorr = cl->GetTheta();  // legacy: origin-based
            if (vtxMode != kVtxOrigin) {
                thetaCorr = ClusterThetaWrtVertex(cl->GetTheta(), cl->GetPhi(),
                                                  vx, vy, vz, R_califa_cm);
            }

            Double_t eDC = DopplerCorrectEnergy(eLab, betaForCorrection, thetaCorr);
            if (extraResFrac > 0. || extraSigAbs > 0.) {
                const Double_t sig = TMath::Sqrt(extraResFrac * eDC * extraResFrac * eDC
                                                 + extraSigAbs * extraSigAbs);
                if (sig > 0.) eDC += gRandom->Gaus(0., sig);
            }
            hResp->Fill(eDC);
            // QA uses the SAME angle that was used for the correction (in deg).
            hEvsTh->Fill(thetaCorr * TMath::RadToDeg(), eDC);
        }
    }

    // Photopeak efficiency: counts within +- window of the line / thrown decays.
    const double Enom = eGamma_keV / 1000.0;
    const double hw   = LineHalfWindowMeV(eGamma_keV);
    const int    b1   = hResp->GetXaxis()->FindBin(Enom - hw + 1e-6);
    const int    b2   = hResp->GetXaxis()->FindBin(Enom + hw - 1e-6);
    const double ppCounts = hResp->Integral(b1, b2);
    const double eps = (nDecays > 0) ? ppCounts / double(nDecays) : 0.0;

    // Normalise template to unit integral -> amplitude = detected photopeak counts
    const double tot = hResp->Integral();
    if (tot > 0.) hResp->Scale(1.0 / tot);

    auto* outFile = TFile::Open(outRoot, recreate ? "RECREATE" : "UPDATE");
    if (!outFile || outFile->IsZombie())
    { Error("buildSingleResponse", "Cannot open output %s", outRoot.Data());
      inFile->Close(); return -1.; }
    outFile->cd();
    hResp->Write(hResp->GetName(), TObject::kOverwrite);
    hEvsTh->Write(hEvsTh->GetName(), TObject::kOverwrite);
    TParameter<double> pEps(EffName(eGamma_keV), eps);
    pEps.Write(pEps.GetName(), TObject::kOverwrite);
    outFile->Close();
    inFile->Close();

    std::cout << "[buildSingleResponse] " << eGamma_keV << " keV: "
              << "eps(photopeak) = " << 100.0 * eps << " %"
              << "  (" << ppCounts << " / " << nDecays << " thrown)\n";
    return eps;
}

// ===========================================================================
//  Convenience: build single-line responses for a list of transitions, and
//  write an EFFICIENCY-vs-ENERGY curve (in percent) for bookkeeping.
//
//  Expects the simulation files to be named sim_gamma_<keV>keV.root, i.e. run
//  califaGammaResponseSimSingle(keV, ...) for each transition first.
//  The FIRST call recreates the output file; the rest update it.
//
//  Outputs (in outRoot), in addition to per-line hResp_line_<keV>,
//  hResp_EvsTheta_<keV> and eps_<keV>:
//     gEfficiencyPercent : TGraph  efficiency [%] vs E_gamma [MeV]
//     hEfficiencyPercent : TH1D    same, one bin per transition (labelled)
// ===========================================================================
void buildAllSingleResponses(Double_t beta = 0.8025,
                             TString simDir = ".",
                             TString outRoot = "./gamma_response_functions.root",
                             Double_t energyScaleToMeV = 1.0,
                             Double_t extraResFrac = 0.079,
                             Double_t extraSigAbs  = 0.0,
                             Int_t    vtxMode      = kVtxTargetZ,
                             Double_t targetZ_cm   = 4.0,
                             Double_t R_califa_cm  = 30.0)
{
    const std::vector<int> lines = GammaCfg::AllTransitionsKeV();

    std::vector<double> eKeV, effPct;
    for (size_t i = 0; i < lines.size(); ++i) {
        const int keV = lines[i];
        const TString simFile = simDir + Form("/sim_gamma_%dkeV.root", keV);
        const double eps = buildSingleResponse(simFile, keV, beta, outRoot,
                                               (i == 0) /*recreate on first*/,
                                               energyScaleToMeV,
                                               extraResFrac, extraSigAbs,
                                               vtxMode, targetZ_cm, R_califa_cm);
        if (eps >= 0.) { eKeV.push_back(keV); effPct.push_back(100.0 * eps); }
    }

    // ---- Efficiency curve (percent) -----------------------------------------
    if (!eKeV.empty()) {
        // sort by energy for a clean curve
        std::vector<size_t> idx(eKeV.size());
        for (size_t i = 0; i < idx.size(); ++i) idx[i] = i;
        std::sort(idx.begin(), idx.end(),
                  [&](size_t a, size_t b){ return eKeV[a] < eKeV[b]; });

        auto* g = new TGraph((int)eKeV.size());
        g->SetName(GammaCfg::EffCurveGraph());
        g->SetTitle("CALIFA #gamma efficiency;E_{#gamma} [MeV];Efficiency [%]");
        g->SetMarkerStyle(20);
        g->SetMarkerSize(1.1);

        auto* h = new TH1D(GammaCfg::EffCurveHist(),
                           "CALIFA #gamma efficiency;E_{#gamma} [MeV];Efficiency [%]",
                           (int)eKeV.size(), 0.5, (int)eKeV.size() + 0.5);
        h->SetDirectory(nullptr);

        for (size_t k = 0; k < idx.size(); ++k) {
            const size_t i = idx[k];
            g->SetPoint((int)k, eKeV[i] / 1000.0, effPct[i]);
            h->SetBinContent((int)k + 1, effPct[i]);
            h->GetXaxis()->SetBinLabel((int)k + 1, Form("%d", (int)eKeV[i]));
        }

        auto* outFile = TFile::Open(outRoot, "UPDATE");
        if (outFile && !outFile->IsZombie()) {
            outFile->cd();
            g->Write(g->GetName(), TObject::kOverwrite);
            h->Write(h->GetName(), TObject::kOverwrite);
            outFile->Close();
        }
        std::cout << "[buildAllSingleResponses] efficiency curve (percent):\n";
        for (size_t k = 0; k < idx.size(); ++k) {
            const size_t i = idx[k];
            std::cout << "    " << std::setw(5) << eKeV[i] << " keV : "
                      << std::fixed << std::setprecision(2) << effPct[i] << " %\n";
        }
    }

    std::cout << "[buildAllSingleResponses] wrote " << lines.size()
              << " single-line templates + efficiencies + EvsTheta QA + eff curve to "
              << outRoot << "\n";
}

// ===========================================================================
// LEGACY builders (unchanged, kept for backwards compatibility).
// NOT used by the population fit.
// ===========================================================================

// buildGammaResponseHist
// For SINGLE gamma simulations: sum all clusters per event, fill once.
void buildGammaResponseHist(TString simFile,
                            TString histName,
                            Double_t eGamma_keV,
                            Double_t betaForCorrection,
                            TString outRoot = "./gamma_response_functions.root",
                            Bool_t recreate = kFALSE,
                            TString treeName = "evt",
                            TString branchName = "CalifaClusterData")
{
    auto* inFile = TFile::Open(simFile, "READ");
    if (!inFile || inFile->IsZombie())
    {
        Error("buildGammaResponseHist", "Cannot open %s", simFile.Data());
        return;
    }
    TTree* tree = dynamic_cast<TTree*>(inFile->Get(treeName));
    if (!tree) { inFile->Close(); return; }

    auto* clusters = new TClonesArray("R3BCalifaClusterData");
    tree->SetBranchAddress(branchName.Data(), &clusters);

    const Double_t eGammaMeV = eGamma_keV * 1.e-3;
    auto* hResp = new TH1D(histName.Data(),
                           Form("CALIFA response, E_{0}=%.0f keV", eGamma_keV),
                           1600, 0, eGammaMeV * 1.8);
    hResp->GetXaxis()->SetTitle("E_{#gamma}^{Doppler corr.} [MeV]");
    hResp->GetYaxis()->SetTitle("Counts");
    hResp->Sumw2();

    for (Long64_t i = 0; i < tree->GetEntries(); ++i)
    {
        tree->GetEntry(i);

        Double_t eRestAddback = 0.;
        for (Int_t j = 0; j < clusters->GetEntriesFast(); ++j)
        {
            auto* cl = dynamic_cast<R3BCalifaClusterData*>(clusters->At(j));
            if (!cl) continue;

            const Double_t eRest = DopplerCorrectEnergy(cl->GetEnergy(),
                                                        betaForCorrection,
                                                        cl->GetTheta());
            if (eRest > 0.) eRestAddback += eRest;
        }
        if (eRestAddback > 0.) hResp->Fill(eRestAddback);
    }

    auto* outFile = TFile::Open(outRoot, recreate ? "RECREATE" : "UPDATE");
    hResp->Write(histName.Data(), TObject::kOverwrite);
    outFile->Close();
    inFile->Close();
    std::cout << "Wrote '" << histName << "' to " << outRoot << "\n";
}

// buildGammaResponseHistPerCluster
// For CASCADE simulations: fill once per cluster (not per event sum).
// WARNING: do NOT use an energy window when producing FIT templates — the
// cascade template must keep BOTH lines. Window kept only for diagnostics.
void buildGammaResponseHistPerCluster(TString simFile,
                                      TString histName,
                                      Double_t eRef_keV,
                                      Double_t betaForCorrection,
                                      TString outRoot = "./gamma_response_functions.root",
                                      Bool_t recreate = kFALSE,
                                      Double_t eMin_keV = 0.,
                                      Double_t eMax_keV = 0.,
                                      TString treeName = "evt",
                                      TString branchName = "CalifaClusterData")
{
    auto* inFile = TFile::Open(simFile, "READ");
    if (!inFile || inFile->IsZombie())
    {
        Error("buildGammaResponseHistPerCluster", "Cannot open %s", simFile.Data());
        return;
    }
    TTree* tree = dynamic_cast<TTree*>(inFile->Get(treeName));
    if (!tree) { inFile->Close(); return; }

    auto* clusters = new TClonesArray("R3BCalifaClusterData");
    tree->SetBranchAddress(branchName.Data(), &clusters);

    const Double_t eRefMeV = eRef_keV * 1.e-3;
    auto* hResp = new TH1D(histName.Data(),
                           Form("CALIFA cascade response (per cluster), E_{ref}=%.0f keV", eRef_keV),
                           1600, 0, eRefMeV * 1.8);
    hResp->GetXaxis()->SetTitle("E_{#gamma}^{Doppler corr.} [MeV]");
    hResp->GetYaxis()->SetTitle("Counts per cluster");
    hResp->Sumw2();

    const bool useWindow = (eMax_keV > eMin_keV && eMax_keV > 0.);
    const Double_t eMinMeV = eMin_keV * 1.e-3;
    const Double_t eMaxMeV = eMax_keV * 1.e-3;

    for (Long64_t i = 0; i < tree->GetEntries(); ++i)
    {
        tree->GetEntry(i);

        for (Int_t j = 0; j < clusters->GetEntriesFast(); ++j)
        {
            auto* cl = dynamic_cast<R3BCalifaClusterData*>(clusters->At(j));
            if (!cl) continue;

            const Double_t eRest = DopplerCorrectEnergy(cl->GetEnergy(),
                                                        betaForCorrection,
                                                        cl->GetTheta());
            if (eRest <= 0.) continue;
            if (useWindow && (eRest < eMinMeV || eRest > eMaxMeV)) continue;
            hResp->Fill(eRest);
        }
    }

    auto* outFile = TFile::Open(outRoot, recreate ? "RECREATE" : "UPDATE");
    hResp->Write(histName.Data(), TObject::kOverwrite);
    outFile->Close();
    inFile->Close();
    std::cout << "Wrote '" << histName << "' to " << outRoot << "\n";
}
