#include "gammaSpectra.h"
#include <TTreeFormula.h>
#include <TLeaf.h>
#include <TBranch.h>
#include <TF1.h>
#include <TObjArray.h>
#include <algorithm>

// ─────────────────────────────────────────────────────────────────────────────
//  gammaSpectra.C — Doppler-corrected CALIFA gamma spectra for one fragment
//                   species (22O signal channel, or 23O / 24O background
//                   channels) from an already-filtered FilterDataTree file.
//
//  WORKFLOW (population analysis):
//    1) run this macro on the 22O-, 23O- and 24O-gated filtered files:
//         root -l -b -q 'gammaSpectra.C("data_22O.root","gamma_spectra_22O.root",
//                                        "","22O.txt",2.5,1.3,1.6,"22O","neuland_mult>=1")'
//         root -l -b -q 'gammaSpectra.C("data_23O.root","gamma_spectra_23O.root",
//                                        "","23O1n.txt",2.5,1.3,1.6,"23O","")'
//         root -l -b -q 'gammaSpectra.C("data_24O.root","gamma_spectra_24O.root",
//                                        "","24O.txt",2.5,1.3,1.6,"24O","")'
//       (23O and 24O have NO bound excited states -> their spectra are pure
//        background templates for the fit, with the same reaction topology.)
//    2) produce response templates with califaGammaResponseSim.C
//    3) fit populations with subtract_background.C
//
//  The multiplicity M is DEFINED here (and identically in the simulation) as
//  the number of clusters passing GammaCfg::PassCluster on the LAB energy.
//  All histogram names / binning come from GammaCfg (gammaSpectra.h).
// ─────────────────────────────────────────────────────────────────────────────

// OPA (proton-proton opening angle) cut in radians – same as fitMomdis.C
static constexpr double OPA_MIN = 1.3;
static constexpr double OPA_MAX = 1.6;

// ─────────────────────────────────────────────────────────────────────────────
//  DopplerQA: slice a (theta, E_corr) 2D histogram in theta and fit a Gaussian
//  to ONE gamma line in each slice, then save the per-slice fit results.
//
//  Purpose: a PERFECT Doppler correction makes the peak CENTROID independent of
//  theta. If the "<line>_mean vs theta" histogram is flat, the correction is
//  good and the measured width is genuine detector resolution. If it TILTS, the
//  correction is imperfect (wrong beta / vertex / angle convention) and part of
//  the 1D peak width is a correction artefact that shows up as a theta-dependent
//  centroid drift.
//
//  Restricting the fit to a Y-window around the line (rather than ROOT's blind
//  FitSlicesY over the whole axis) is essential — otherwise the fit latches onto
//  the wrong structure in slices where the line is weak.
//
//  Writes (into the current directory):
//    <base>_mean   : Gaussian mean   vs theta   (THE diagnostic — should be flat)
//    <base>_sigma  : Gaussian sigma  vs theta
//    <base>_amp    : Gaussian amp    vs theta
//    <base>_chi2ndf: reduced chi2    vs theta
// ─────────────────────────────────────────────────────────────────────────────
static void DopplerQA_FitLineVsTheta(TH2F* h2,
                                     const char* base,
                                     double lineCenterMeV,
                                     double yLoMeV,
                                     double yHiMeV,
                                     int groupTheta = 10,          // theta bins merged per fit
                                     int minEntriesPerSlice = 40)
{
    if (!h2) return;
    if (groupTheta < 1) groupTheta = 1;

    // Merge groupTheta adjacent theta bins so each slice has enough statistics.
    // With the default 1 deg/bin booking, groupTheta=10 => 10 deg slices.
    TH2F* hg = h2;
    bool  ownHg = false;
    if (groupTheta > 1) {
        hg = static_cast<TH2F*>(h2->Clone(Form("%s_grp", base)));
        hg->SetDirectory(nullptr);
        hg->RebinX(groupTheta);   // merge along theta (X); Y untouched
        ownHg = true;
    }

    const int    nx   = hg->GetNbinsX();          // theta super-bins
    const double xlo  = hg->GetXaxis()->GetXmin();
    const double xhi  = hg->GetXaxis()->GetXmax();

    auto* hMean = new TH1F(Form("%s_mean",   base),
        Form("%s: %.2f MeV peak mean vs #theta (#Delta#theta=%d bins);#theta_{lab} [deg];peak mean [MeV]",
             base, lineCenterMeV, groupTheta), nx, xlo, xhi);
    auto* hSig  = new TH1F(Form("%s_sigma",  base),
        Form("%s: %.2f MeV peak #sigma vs #theta;#theta_{lab} [deg];peak #sigma [MeV]",
             base, lineCenterMeV), nx, xlo, xhi);
    auto* hAmp  = new TH1F(Form("%s_amp",    base),
        Form("%s: %.2f MeV peak amp vs #theta;#theta_{lab} [deg];amplitude",
             base, lineCenterMeV), nx, xlo, xhi);
    auto* hChi  = new TH1F(Form("%s_chi2ndf",base),
        Form("%s: %.2f MeV peak #chi^{2}/ndf vs #theta;#theta_{lab} [deg];#chi^{2}/ndf",
             base, lineCenterMeV), nx, xlo, xhi);
    hMean->SetDirectory(nullptr); hSig->SetDirectory(nullptr);
    hAmp ->SetDirectory(nullptr); hChi->SetDirectory(nullptr);

    for (int ix = 1; ix <= nx; ++ix) {
        TH1D* proj = hg->ProjectionY(Form("%s_py_%d", base, ix), ix, ix);
        if (!proj) continue;
        proj->SetDirectory(nullptr);

        // integral within the fit window only
        const int by1 = proj->GetXaxis()->FindBin(yLoMeV + 1e-6);
        const int by2 = proj->GetXaxis()->FindBin(yHiMeV - 1e-6);
        if (proj->Integral(by1, by2) < minEntriesPerSlice) { delete proj; continue; }

        TF1 fg(Form("%s_g_%d", base, ix), "gaus", yLoMeV, yHiMeV);
        fg.SetParameters(proj->GetMaximum(), lineCenterMeV, 0.30);
        fg.SetParLimits(1, yLoMeV, yHiMeV);     // mean stays in window
        fg.SetParLimits(2, 0.05, 0.80);         // sensible sigma range
        const int st = proj->Fit(&fg, "QNR");

        if (st == 0 || st == 4000) {
            hMean->SetBinContent(ix, fg.GetParameter(1)); hMean->SetBinError(ix, fg.GetParError(1));
            hSig ->SetBinContent(ix, fg.GetParameter(2)); hSig ->SetBinError(ix, fg.GetParError(2));
            hAmp ->SetBinContent(ix, fg.GetParameter(0)); hAmp ->SetBinError(ix, fg.GetParError(0));
            hChi ->SetBinContent(ix, fg.GetNDF() > 0 ? fg.GetChisquare()/fg.GetNDF() : 0.);
        }
        delete proj;
    }

    hMean->Write(); hSig->Write(); hAmp->Write(); hChi->Write();
    if (ownHg) delete hg;
    std::cout << "[DopplerQA] " << base << ": wrote _mean/_sigma/_amp/_chi2ndf "
              << "for the " << lineCenterMeV << " MeV line, "
              << groupTheta << " theta-bins/slice "
              << "(check that _mean is FLAT vs theta).\n";
}

// k-sigma for PID ellipse (same convention as CrossSections.cpp)
static constexpr double PID_K_DEFAULT = 2.5;

// ─────────────────────────────────────────────────────────────────────────────
//  Outgoing-fragment PID window presets (seed windows for the 2D Gaussian fit;
//  the fit refines the centroids/sigmas, the ellipse cut uses the fit result).
//  23O / 24O seeds are the 22O window shifted by the nominal A/Z spacing.
//  Adjust here if your calibrated A/Z scale differs.
// ─────────────────────────────────────────────────────────────────────────────
struct PidWindow { double aoqMin, aoqMax, zMin, zMax; };

static PidWindow pidPresetFor(const TString& species)
{
    if (species == "23O") return { 2.79, 2.89, 7.6, 8.5 };
    if (species == "24O") return { 2.91, 3.01, 7.6, 8.5 };
    // default: 22O
    return { 2.67, 2.76, 7.6, 8.5 };
}

// 25F incoming PID window (used only on the unreacted file if provided)
static constexpr double AOQ_25F_MIN = 2.71;
static constexpr double AOQ_25F_MAX = 2.77;
static constexpr double Z_25F_MIN   = 8.3;
static constexpr double Z_25F_MAX   = 9.5;

// ─────────────────────────────────────────────────────────────────────────────
//  Helpers
// ─────────────────────────────────────────────────────────────────────────────

// Open a FilterDataTree file and return the TTree; caller does NOT own TFile.
static TTree* openTree(TString path, TFile*& fout)
{
    fout = TFile::Open(path, "READ");
    if (!fout || fout->IsZombie()) {
        std::cerr << "[ERROR] Cannot open file: " << path << "\n";
        fout = nullptr;
        return nullptr;
    }
    auto* t = dynamic_cast<TTree*>(fout->Get("FilterDataTree"));
    if (!t)
        std::cerr << "[ERROR] FilterDataTree not found in " << path << "\n";
    return t;
}

// ─────────────────────────────────────────────────────────────────────────────
//  PID pass: fill a 2D histogram, fit 2D Gaussian, return parameters
// ─────────────────────────────────────────────────────────────────────────────
static Fit2DParams fitPID(TTree* tree,
                           const char* hname,
                           double aoqMin, double aoqMax,
                           double zMin,   double zMax,
                           int uid)
{
    int nBins = (tree->GetEntries() > 5000) ? 500 : 100;
    auto* h = new TH2F(hname,
                        Form("%s;AoQ;Z", hname),
                        nBins, aoqMin - 0.04, aoqMax + 0.04,
                        nBins, zMin   - 0.15, zMax   + 0.15);
    h->SetDirectory(nullptr);

    double z = 0, aoq = 0;
    tree->SetBranchStatus("*",         0);
    tree->SetBranchStatus("Z_frag_est",1);
    tree->SetBranchStatus("AoQ_frag",  1);
    tree->SetBranchAddress("Z_frag_est", &z);
    tree->SetBranchAddress("AoQ_frag",   &aoq);

    for (Long64_t i = 0; i < tree->GetEntries(); ++i) {
        tree->GetEntry(i);
        if (aoq > aoqMin && aoq < aoqMax && z > zMin && z < zMax)
            h->Fill(aoq, z);
    }

    tree->SetBranchStatus("*", 1);
    tree->ResetBranchAddresses();

    std::cout << "[PID] " << hname << ": " << h->GetEntries()
              << " entries in window, fitting 2D Gaussian...\n";

    Fit2DParams p = Fit2DGaussian(h, aoqMin, aoqMax, zMin, zMax, uid);
    h->Write(); // save PID histogram to output file

    if (p.valid)
        std::cout << "[PID] muAoQ=" << p.muAoQ << " sigAoQ=" << p.sigmaAoQ
                  << "  muZ=" << p.muZ << " sigZ=" << p.sigmaZ << "\n";
    else
        std::cerr << "[WARN] 2D Gaussian fit did not converge — will use rectangular window.\n";

    delete h;
    return p;
}

// ─────────────────────────────────────────────────────────────────────────────
//  Main analysis
// ─────────────────────────────────────────────────────────────────────────────

void gammaSpectra(TString reactionFile,
                  TString outFile       = "gamma_spectra_22O.root",
                  TString unreactedFile = "",          // optional: data_25F.root
                  TString offsetFile    = "",   // settings file in $repopath/final/settings/
                  double  kPID          = PID_K_DEFAULT,
                  double  opaMin        = OPA_MIN,
                  double  opaMax        = OPA_MAX,
                  TString species       = "22O",  // "22O" | "23O" | "24O"  (PID preset)
                  TString neulandCut    = "",     // e.g. "neuland_mult>=1" ; "" = no cut
                  // ── EXACT vertex+geometry Doppler correction (default) ───────
                  //  Uses the FOOT-reconstructed vertex (vertex_x/y/z) as the
                  //  emission point and R3BCalifaGeometry to get each gamma's TRUE
                  //  mother-crystal centre from its crystal ID — NO CALIFA-radius
                  //  assumption. Falls back to origin-based (theta,phi) if the
                  //  branches or geometry are unavailable.
                  Bool_t  useVertex     = kTRUE,
                  Int_t   califaGeoVersion = 2025,          // g249
                  Double_t maxVertexDCA = 0.2,   // cm; reject poorly reconstructed vertices
                  TString vtxXBranch    = "vertex_x",
                  TString vtxYBranch    = "vertex_y",
                  TString vtxZBranch    = "vertex_z",
                  TString vtxDcaBranch  = "vertex_dca",
                  TString crystalIDBranch = "califa_gamma_crystalID")
{
    using namespace GammaCfg;

    // ── Resolve output path ───────────────────────────────────────────────────
    // Bare filenames are placed in $repopath/results/final/; absolute paths
    // (containing '/') are used as-is.
    if (!outFile.Contains("/")) {
        const char* repopath = getenv("repopath");
        if (repopath)
            outFile = TString(repopath) + "/results/final/" + outFile;
        else
            std::cerr << "[WARN] $repopath not set — writing output to current directory\n";
    }

    // ── Open output file ─────────────────────────────────────────────────────
    TFile* fout_root = new TFile(outFile, "RECREATE");
    if (!fout_root || fout_root->IsZombie()) {
        std::cerr << "[ERROR] Cannot create output file: " << outFile << "\n";
        return;
    }

    const PidWindow pw = pidPresetFor(species);

    std::cout << "\n=== gamma spectroscopy: 25F(p,2p) -> " << species << " ===\n";
    std::cout << "  Reaction file : " << reactionFile   << "\n";
    std::cout << "  Output file   : " << outFile        << "\n";
    std::cout << "  Unreacted file: " << (unreactedFile.IsNull() ? "(none)" : unreactedFile.Data()) << "\n";
    std::cout << "  Offset file   : " << offsetFile     << "\n";
    std::cout << "  kPID          : " << kPID    << "\n";
    std::cout << "  OPA cut [rad] : " << opaMin  << " - " << opaMax << "\n";
    std::cout << "  PID seed      : AoQ [" << pw.aoqMin << "," << pw.aoqMax
              << "]  Z [" << pw.zMin << "," << pw.zMax << "]\n";
    std::cout << "  NeuLAND cut   : " << (neulandCut.IsNull() ? "(none)" : neulandCut.Data()) << "\n";
    std::cout << "  Cluster cut   : lab E in (" << ECL_LAB_MIN_MEV << ", "
              << ECL_LAB_MAX_MEV << ") MeV  (defines multiplicity M)\n\n";

    // ── Load kinematic offsets ────────────────────────────────────────────────
    // Same format as 23O1n.txt / 24O.txt in $repopath/final/settings/:
    //   line 0: fragment x-direction offset (fx)
    //   line 1: fragment y-direction offset (fy)
    //   line 2: neutron x-direction offset  (unused here)
    //   line 3: neutron y-direction offset  (unused here)
    //   line 4: beta matching offset (added to beta_frag before Doppler correction)
    OffsetParams offsets;
    if (!offsetFile.IsNull() && offsetFile.Length() > 0)
        offsets = readOffsets(offsetFile.Data());

    // ── Open reaction tree ───────────────────────────────────────────────────
    TFile* freac = nullptr;
    TTree* treac = openTree(reactionFile, freac);
    if (!treac) return;
    std::cout << "[INFO] Reaction tree: " << treac->GetEntries() << " entries\n";

    // ── Pass 1: fit outgoing fragment PID ────────────────────────────────────
    fout_root->cd();
    Fit2DParams pidFrag = fitPID(treac, Form("h2_pid_%s", species.Data()),
                                  pw.aoqMin, pw.aoqMax,
                                  pw.zMin,   pw.zMax, 1);

    // Fall back to rectangular window if fit failed
    if (!pidFrag.valid) {
        pidFrag.muAoQ    = 0.5 * (pw.aoqMin + pw.aoqMax);
        pidFrag.sigmaAoQ = (pw.aoqMax - pw.aoqMin) / (2.0 * kPID);
        pidFrag.muZ      = 0.5 * (pw.zMin + pw.zMax);
        pidFrag.sigmaZ   = (pw.zMax - pw.zMin) / (2.0 * kPID);
    }

    // ── Optional: incoming 25F PID from unreacted file ───────────────────────
    // The eventFilter already applied a graphical incoming cut; this block
    // extracts the Gaussian parameters for diagnostic purposes.
    Fit2DParams pid25F;
    if (!unreactedFile.IsNull() && unreactedFile.Length() > 0) {
        TFile* funr = nullptr;
        TTree* tunr = openTree(unreactedFile, funr);
        if (tunr) {
            std::cout << "[INFO] Unreacted tree (25F): " << tunr->GetEntries() << " entries\n";
            fout_root->cd();
            pid25F = fitPID(tunr, "h2_pid_25F",
                             AOQ_25F_MIN, AOQ_25F_MAX,
                             Z_25F_MIN,   Z_25F_MAX, 2);
            if (pid25F.valid)
                std::cout << "[25F incoming] muAoQ=" << pid25F.muAoQ
                          << " muZ=" << pid25F.muZ << "\n";
        }
        if (funr) funr->Close();
    }

    // ── Create analysis histograms ────────────────────────────────────────────
    fout_root->cd();

    // ---- Diagnostic / QA ----
    auto* h2_pid_all = new TH2F("h2_pid_all",
        "All events: outgoing PID;AoQ;Z",
        500, pw.aoqMin-0.05, pw.aoqMax+0.05,
        500, pw.zMin-0.2,   pw.zMax+0.2);

    auto* h1_opa_after_pid = new TH1F("h1_opa_after_pid",
        "OPA after fragment PID (no OPA cut);#theta_{OPA} [rad];Counts",
        180, 0.0, TMath::Pi());

    auto* h1_beta_frag = new TH1F("h1_beta_frag",
        "Fragment #beta (after all cuts);#beta_{frag};Counts",
        200, 0.78, 0.84);

    // ---- Gamma cluster multiplicity (recomputed: # clusters passing lab-E cut) ----
    auto* h1_gamma_mult = new TH1F("h1_gamma_mult",
        "Gamma cluster multiplicity (after all cuts);M_{#gamma};Events",
        16, -0.5, 15.5);

    // ---- Raw (uncorrected) gamma energy, all multiplicities ----
    auto* h1_gamma_E_raw_all = new TH1F("h1_gamma_E_raw_all",
        "Raw #gamma energy, all mult;E_{#gamma} [MeV];Counts",
        NBINS_E, E_MIN_MEV, E_MAX_MEV);

    // ---- Doppler-corrected gamma energy: all / M==m / M<=m / M>=m ----
    auto* h1_gamma_E_corr_all = new TH1F(DataHist(SliceAll()),
        "Doppler-corrected #gamma, all mult;E_{#gamma} [MeV];Counts",
        NBINS_E, E_MIN_MEV, E_MAX_MEV);

    TH1F* h1_gamma_E_corr_mult[MULT_MAX];
    TH1F* h1_gamma_E_corr_mleq[MULT_MAX];
    for (int m = 1; m <= MULT_MAX; ++m) {
        h1_gamma_E_corr_mult[m-1] = new TH1F(
            DataHist(SliceExact(m)),
            Form("Doppler-corrected #gamma, M_{#gamma}==%d;E_{#gamma} [MeV];Counts", m),
            NBINS_E, E_MIN_MEV, E_MAX_MEV);
        h1_gamma_E_corr_mleq[m-1] = new TH1F(
            DataHist(SliceLeq(m)),
            Form("Doppler-corrected #gamma, M_{#gamma}<=%d;E_{#gamma} [MeV];Counts", m),
            NBINS_E, E_MIN_MEV, E_MAX_MEV);
    }

    // M >= m slices, m = 2..MULT_MAX (no upper cap on M)
    TH1F* h1_gamma_E_corr_mgeq[MULT_MAX + 1] = { nullptr };
    for (int m = 2; m <= MULT_MAX; ++m) {
        h1_gamma_E_corr_mgeq[m] = new TH1F(
            DataHist(SliceGeq(m)),
            Form("Doppler-corrected #gamma, M_{#gamma}>=%d;E_{#gamma} [MeV];Counts", m),
            NBINS_E, E_MIN_MEV, E_MAX_MEV);
    }

    // ---- Gamma-gamma gated spectra (companion clusters, M>=2 events) ----
    //  gate1383: spectrum of the OTHER clusters when one cluster falls in the
    //            1383 keV gate  -> cascade shows a prominent 3199 peak here.
    //  gate3199: vice versa     -> cascade shows a prominent 1383 peak here.
    auto* h1_gamma_E_corr_gate1383 = new TH1F("h1_gamma_E_corr_gate1383",
        Form("Companions of clusters in [%.2f,%.2f] MeV;E_{#gamma} [MeV];Counts",
             GATE_1383_LO, GATE_1383_HI),
        NBINS_E, E_MIN_MEV, E_MAX_MEV);

    auto* h1_gamma_E_corr_gate3199 = new TH1F("h1_gamma_E_corr_gate3199",
        Form("Companions of clusters in [%.2f,%.2f] MeV;E_{#gamma} [MeV];Counts",
             GATE_3199_LO, GATE_3199_HI),
        NBINS_E, E_MIN_MEV, E_MAX_MEV);

    // ---- M==2 diagnostics: E1 vs E2 scatter (cascade island at (1.38, 3.20))
    //      and calorimetric sum (cascade peak at 4.58 MeV) ----
    auto* h2_gamma_E1_vs_E2_mult2 = new TH2F("h2_gamma_E1_vs_E2_mult2",
        "M_{#gamma}==2;E_{low} [MeV];E_{high} [MeV]",
        NBINS_E, E_MIN_MEV, E_MAX_MEV,
        NBINS_E, E_MIN_MEV, E_MAX_MEV);

    auto* h1_gamma_Esum_mult2 = new TH1F("h1_gamma_Esum_mult2",
        "M_{#gamma}==2: E_{1}+E_{2};E_{sum} [MeV];Events",
        NBINS_E, E_MIN_MEV, E_MAX_MEV);

    // ---- 2D: Doppler-corrected energy vs theta (lab frame, degrees) ----
    auto* h2_gamma_Ecorr_vs_theta_all = new TH2F("h2_gamma_Ecorr_vs_theta_all",
        "Doppler-corrected #gamma, all mult;#theta_{lab} [deg];E_{#gamma}^{DC} [MeV]",
        NBINS_TH, TH_MIN_DEG, TH_MAX_DEG,
        NBINS_E,  E_MIN_MEV,  E_MAX_MEV);

    // ---- 2D: raw energy vs theta ----
    auto* h2_gamma_E_vs_theta_all = new TH2F("h2_gamma_E_vs_theta_all",
        "Raw #gamma energy vs #theta, all mult;#theta_{lab} [deg];E_{#gamma}^{raw} [MeV]",
        NBINS_TH, TH_MIN_DEG, TH_MAX_DEG,
        NBINS_E,  E_MIN_MEV,  E_MAX_MEV);

    // ---- 2D: Doppler-corrected energy vs theta by mult<=2 and <=3 ----
    auto* h2_gamma_Ecorr_vs_theta_mleq2 = new TH2F("h2_gamma_Ecorr_vs_theta_mleq2",
        "Doppler-corrected #gamma, M_{#gamma}<=2;#theta_{lab} [deg];E_{#gamma}^{DC} [MeV]",
        NBINS_TH, TH_MIN_DEG, TH_MAX_DEG,
        NBINS_E,  E_MIN_MEV,  E_MAX_MEV);

    auto* h2_gamma_Ecorr_vs_theta_mleq3 = new TH2F("h2_gamma_Ecorr_vs_theta_mleq3",
        "Doppler-corrected #gamma, M_{#gamma}<=3;#theta_{lab} [deg];E_{#gamma}^{DC} [MeV]",
        NBINS_TH, TH_MIN_DEG, TH_MAX_DEG,
        NBINS_E,  E_MIN_MEV,  E_MAX_MEV);

    // ── Pass 2: physics loop ─────────────────────────────────────────────────
    double Z_frag_est = 0, AoQ_frag  = 0;
    double califa_opa = 0, beta_frag  = 0;
    double px_frag    = 0, py_frag    = 0, pz_frag = 0;
    int    califa_gamma_mult = 0;
    std::vector<double>* califa_gamma_E     = nullptr;
    std::vector<double>* califa_gamma_theta = nullptr;
    std::vector<double>* califa_gamma_phi   = nullptr;

    treac->SetBranchStatus("*", 0);
    for (const char* br : {"Z_frag_est","AoQ_frag","califa_opa","beta_frag",
                           "px_frag","py_frag","pz_frag","califa_gamma_mult",
                           "califa_gamma_E","califa_gamma_theta","califa_gamma_phi"})
        treac->SetBranchStatus(br, 1);

    treac->SetBranchAddress("Z_frag_est",        &Z_frag_est);
    treac->SetBranchAddress("AoQ_frag",          &AoQ_frag);
    treac->SetBranchAddress("califa_opa",        &califa_opa);
    treac->SetBranchAddress("beta_frag",         &beta_frag);
    treac->SetBranchAddress("px_frag",           &px_frag);
    treac->SetBranchAddress("py_frag",           &py_frag);
    treac->SetBranchAddress("pz_frag",           &pz_frag);
    treac->SetBranchAddress("califa_gamma_mult", &califa_gamma_mult);
    treac->SetBranchAddress("califa_gamma_E",    &califa_gamma_E);
    treac->SetBranchAddress("califa_gamma_theta",&califa_gamma_theta);
    treac->SetBranchAddress("califa_gamma_phi",  &califa_gamma_phi);

    // ── Vertex + crystal-ID correction: detect branches, init geometry ───────
    double vtxX = 0., vtxY = 0., vtxZ = 0., vtxDca = 0.;
    std::vector<int>* gCrystalID = nullptr;
    bool vtxOK = false;      // reconstructed vertex available
    bool crysOK = false;     // per-gamma crystal ID available
    R3BCalifaGeometry* califaGeo = nullptr;

    auto hasBranch = [&](const TString& b) -> bool {
        return (!b.IsNull() && b.Length() > 0 && treac->GetBranch(b) != nullptr);
    };

    if (useVertex) {
        if (hasBranch(vtxXBranch) && hasBranch(vtxYBranch) && hasBranch(vtxZBranch)) {
            treac->SetBranchStatus(vtxXBranch, 1);
            treac->SetBranchStatus(vtxYBranch, 1);
            treac->SetBranchStatus(vtxZBranch, 1);
            treac->SetBranchAddress(vtxXBranch, &vtxX);
            treac->SetBranchAddress(vtxYBranch, &vtxY);
            treac->SetBranchAddress(vtxZBranch, &vtxZ);
            vtxOK = true;
            std::cout << "[vertex] reconstructed vertex from "
                      << vtxXBranch << "/" << vtxYBranch << "/" << vtxZBranch << "\n";
            if (hasBranch(vtxDcaBranch)) {
                treac->SetBranchStatus(vtxDcaBranch, 1);
                treac->SetBranchAddress(vtxDcaBranch, &vtxDca);
                std::cout << "[vertex] DCA quality cut: vertex_dca < " << maxVertexDCA << " cm\n";
            }
            // crystal ID for exact geometry lookup
            if (hasBranch(crystalIDBranch)) {
                treac->SetBranchStatus(crystalIDBranch, 1);
                treac->SetBranchAddress(crystalIDBranch, &gCrystalID);
                califaGeo = InitCalifaGeometry(califaGeoVersion);
                crysOK = (califaGeo != nullptr);
                std::cout << "[vertex] EXACT crystal-centre lookup via R3BCalifaGeometry("
                          << califaGeoVersion << ") — NO radius assumption.\n";
            } else {
                std::cerr << "[WARN] crystal-ID branch '" << crystalIDBranch
                          << "' not found — using (theta,phi) origin direction minus vertex.\n";
            }
        } else {
            std::cerr << "[WARN] useVertex requested but vertex branches not found — "
                         "falling back to ORIGIN-based correction.\n";
        }
    }

    // Vertex-distribution QA histograms (filled after all cuts).
    auto* h1_vertex_X = new TH1F("h1_vertex_X", "Reconstructed vertex X;X [cm];Events", 200, -5, 5);
    auto* h1_vertex_Y = new TH1F("h1_vertex_Y", "Reconstructed vertex Y;Y [cm];Events", 200, -5, 5);
    auto* h1_vertex_Z = new TH1F("h1_vertex_Z", "Reconstructed vertex Z;Z [cm];Events", 200, -10, 20);
    auto* h1_vertex_DCA = new TH1F("h1_vertex_DCA", "Vertex DCA;DCA [cm];Events", 200, 0, 2);
    auto* h2_vertex_XY = new TH2F("h2_vertex_XY", "Vertex X vs Y;X [cm];Y [cm]", 200, -5, 5, 200, -5, 5);
    auto* h2_vertex_XZ = new TH2F("h2_vertex_XZ", "Vertex Z vs X;Z [cm];X [cm]", 200, -10, 20, 200, -5, 5);
    auto* h2_vertex_YZ = new TH2F("h2_vertex_YZ", "Vertex Z vs Y;Z [cm];Y [cm]", 200, -10, 20, 200, -5, 5);

    // ── Optional NeuLAND cut: arbitrary formula on tree branches ─────────────
    // Interim 2n proxy, e.g. "neuland_mult>=1". Later replace with the real
    // 2n-selection expression without touching anything else.
    TTreeFormula* fNeu = nullptr;
    if (!neulandCut.IsNull() && neulandCut.Length() > 0) {
        fNeu = new TTreeFormula("fNeu", neulandCut, treac);
        if (fNeu->GetNdim() == 0) {
            std::cerr << "[WARN] NeuLAND cut '" << neulandCut
                      << "' is not a valid formula on this tree — cut DISABLED.\n";
            delete fNeu;
            fNeu = nullptr;
        } else {
            // (Re-)enable the branches the formula needs.
            for (Int_t c = 0; c < fNeu->GetNcodes(); ++c) {
                TLeaf* lf = fNeu->GetLeaf(c);
                if (lf && lf->GetBranch()) lf->GetBranch()->SetStatus(1);
            }
            fNeu->UpdateFormulaLeaves();
        }
    }

    const Long64_t nEntries = treac->GetEntries();
    Long64_t nPassPID = 0, nPassOPA = 0, nPassNeu = 0;
    Long64_t nMultMismatch = 0, nEventsM0 = 0, nCoinc = 0;

    // scratch arrays for accepted clusters of the current event
    std::vector<double> eCorr, eRaw, thDeg;

    for (Long64_t i = 0; i < nEntries; ++i) {
        treac->GetEntry(i);

        if (i % 100000 == 0)
            std::cout << "\r  processing event " << i << " / " << nEntries << std::flush;

        // ── Outgoing fragment PID ellipse cut ─────────────────────────────
        h2_pid_all->Fill(AoQ_frag, Z_frag_est);

        if (!insideEllipse(AoQ_frag, Z_frag_est,
                           pidFrag.muAoQ, pidFrag.sigmaAoQ,
                           pidFrag.muZ,   pidFrag.sigmaZ,  kPID))
            continue;
        ++nPassPID;

        h1_opa_after_pid->Fill(califa_opa);

        // ── OPA cut ──────────────────────────────────────────────────────
        if (califa_opa < opaMin || califa_opa > opaMax) continue;
        ++nPassOPA;

        // ── NeuLAND cut (optional, interim 2n proxy) ─────────────────────
        if (fNeu) {
            fNeu->GetNdata();
            if (fNeu->EvalInstance(0) == 0.0) continue;
        }
        ++nPassNeu;   // <-- nFragEvents: denominator, counted with NO gamma condition

        // ── Vertex QA fills (only when a reconstructed vertex is present) ────
        if (vtxOK && vtxZ > -900.0) {
            h1_vertex_DCA->Fill(vtxDca);
            if (vtxDca < maxVertexDCA) {
                h1_vertex_X->Fill(vtxX);
                h1_vertex_Y->Fill(vtxY);
                h1_vertex_Z->Fill(vtxZ);
                h2_vertex_XY->Fill(vtxX, vtxY);
                h2_vertex_XZ->Fill(vtxZ, vtxX);
                h2_vertex_YZ->Fill(vtxZ, vtxY);
            }
        }

        // ── Apply kinematic offsets (same prescription as dataAnalysis.cpp) ──
        // Beta: additive correction from beta matching (line 4 of settings file)
        const double beta_corr = beta_frag + offsets.betaMatch;

        // Fragment direction: correct the px/pz and py/pz slopes then rebuild
        // the unit vector.  Equivalent to dataAnalysis:
        //   fx_corr = px_frag/pz_frag - fragOffsetX
        const double pz_safe = (std::abs(pz_frag) > 0.0) ? pz_frag : 1.0;
        const double px_corr = px_frag - offsets.fragOffsetX * pz_safe;
        const double py_corr = py_frag - offsets.fragOffsetY * pz_safe;

        h1_beta_frag->Fill(beta_corr);

        // ── Collect ACCEPTED clusters and recompute the multiplicity M ───
        // M := number of clusters with lab energy inside GammaCfg window.
        // The vector size is the source of truth (guards against a stale
        // califa_gamma_mult branch).
        eCorr.clear(); eRaw.clear(); thDeg.clear();
        const int nCl = califa_gamma_E ? static_cast<int>(califa_gamma_E->size()) : 0;
        if (nCl != califa_gamma_mult) ++nMultMismatch;

        for (int ig = 0; ig < nCl; ++ig) {
            // energies stored in keV → convert to MeV
            const double E_raw_mev  = (*califa_gamma_E)[ig]     / 1000.0;
            if (!GammaCfg::PassCluster(E_raw_mev)) continue;

            const double theta_rad  = (*califa_gamma_theta)[ig];
            const double phi_rad    = (*califa_gamma_phi)[ig];

            // Doppler correction:
            //  * EXACT: vertex (FOOT) -> true crystal centre (from crystal ID via
            //    R3BCalifaGeometry), NO radius assumption;
            //  * else vertex minus (theta,phi) origin direction;
            //  * else legacy origin-based.
            // The DCA cut rejects events whose vertex is poorly reconstructed.
            const bool vtxUsable = vtxOK && (vtxDca < maxVertexDCA) && (vtxZ > -900.0);
            double E_corr_mev;
            if (vtxUsable) {
                int cid = -1;
                if (crysOK && gCrystalID && ig < (int)gCrystalID->size())
                    cid = (*gCrystalID)[ig];
                E_corr_mev = DopplerCorrectCrystal(
                    E_raw_mev, cid, theta_rad, phi_rad,
                    px_corr, py_corr, pz_safe, beta_corr,
                    vtxX, vtxY, vtxZ, califaGeo);
            } else {
                E_corr_mev = DopplerCorrect(E_raw_mev,
                                            theta_rad, phi_rad,
                                            px_corr, py_corr, pz_safe,
                                            beta_corr);
            }
            eRaw .push_back(E_raw_mev);
            eCorr.push_back(E_corr_mev);
            thDeg.push_back(theta_rad * TMath::RadToDeg());
        }

        const int mult = static_cast<int>(eCorr.size());
        h1_gamma_mult->Fill(static_cast<double>(mult));
        if (mult == 0) { ++nEventsM0; continue; }

        // ── Per-cluster fills ─────────────────────────────────────────────
        for (int ig = 0; ig < mult; ++ig) {
            // Raw histograms
            h1_gamma_E_raw_all->Fill(eRaw[ig]);
            h2_gamma_E_vs_theta_all->Fill(thDeg[ig], eRaw[ig]);

            // Corrected histograms (all mult)
            h1_gamma_E_corr_all->Fill(eCorr[ig]);
            h2_gamma_Ecorr_vs_theta_all->Fill(thDeg[ig], eCorr[ig]);

            // Exact multiplicity histograms (mult == m, m=1..MULT_MAX)
            if (mult >= 1 && mult <= GammaCfg::MULT_MAX)
                h1_gamma_E_corr_mult[mult-1]->Fill(eCorr[ig]);

            // Cumulative multiplicity histograms (mult <= m)
            for (int m = 1; m <= GammaCfg::MULT_MAX; ++m)
                if (mult <= m) h1_gamma_E_corr_mleq[m-1]->Fill(eCorr[ig]);

            // Reverse-cumulative histograms (mult >= m)
            for (int m = 2; m <= GammaCfg::MULT_MAX; ++m)
                if (mult >= m) h1_gamma_E_corr_mgeq[m]->Fill(eCorr[ig]);

            // 2D E vs theta for M<=2 and M<=3
            if (mult <= 2) h2_gamma_Ecorr_vs_theta_mleq2->Fill(thDeg[ig], eCorr[ig]);
            if (mult <= 3) h2_gamma_Ecorr_vs_theta_mleq3->Fill(thDeg[ig], eCorr[ig]);
        }

        // ── Gamma-gamma gated spectra + coincidence counting (M>=2) ──────
        if (mult >= 2) {
            for (int j = 0; j < mult; ++j) {
                for (int k = 0; k < mult; ++k) {
                    if (k == j) continue;
                    if (GammaCfg::In1383Gate(eCorr[j]))
                        h1_gamma_E_corr_gate1383->Fill(eCorr[k]);
                    if (GammaCfg::In3199Gate(eCorr[j]))
                        h1_gamma_E_corr_gate3199->Fill(eCorr[k]);
                }
            }
            // photopeak-photopeak coincidence pairs (each unordered pair once)
            for (int j = 0; j < mult; ++j)
                for (int k = j + 1; k < mult; ++k) {
                    const bool c1 = GammaCfg::In1383Gate(eCorr[j]) && GammaCfg::In3199Gate(eCorr[k]);
                    const bool c2 = GammaCfg::In1383Gate(eCorr[k]) && GammaCfg::In3199Gate(eCorr[j]);
                    if (c1 || c2) ++nCoinc;
                }
        }

        // ── M==2 diagnostics ──────────────────────────────────────────────
        if (mult == 2) {
            const double eLo = std::min(eCorr[0], eCorr[1]);
            const double eHi = std::max(eCorr[0], eCorr[1]);
            h2_gamma_E1_vs_E2_mult2->Fill(eLo, eHi);
            h1_gamma_Esum_mult2->Fill(eLo + eHi);
        }
    }

    std::cout << "\n\n[SUMMARY]\n";
    std::cout << "  Total events           : " << nEntries  << "\n";
    std::cout << "  After PID cut          : " << nPassPID  << "\n";
    std::cout << "  After OPA cut          : " << nPassOPA  << "\n";
    std::cout << "  After NeuLAND cut      : " << nPassNeu
              << "   <-- nFragEvents (population denominator)\n";
    std::cout << "  Events with M==0       : " << nEventsM0 << "\n";
    std::cout << "  1383(x)3199 coinc pairs: " << nCoinc    << "\n";
    if (nMultMismatch > 0)
        std::cout << "  [WARN] califa_gamma_mult != vector size in "
                  << nMultMismatch << " events (vector size used)\n";

    // ── Write all histograms + book-keeping parameters ───────────────────────
    fout_root->cd();

    h2_pid_all->Write();
    h1_opa_after_pid->Write();
    h1_beta_frag->Write();
    h1_gamma_mult->Write();

    // vertex-distribution QA
    if (vtxOK) {
        h1_vertex_X->Write();  h1_vertex_Y->Write();  h1_vertex_Z->Write();
        h1_vertex_DCA->Write();
        h2_vertex_XY->Write(); h2_vertex_XZ->Write(); h2_vertex_YZ->Write();
    }
    h1_gamma_E_raw_all->Write();
    h1_gamma_E_corr_all->Write();
    for (int m = 1; m <= GammaCfg::MULT_MAX; ++m) {
        h1_gamma_E_corr_mult[m-1]->Write();
        h1_gamma_E_corr_mleq[m-1]->Write();
    }
    for (int m = 2; m <= GammaCfg::MULT_MAX; ++m)
        h1_gamma_E_corr_mgeq[m]->Write();
    h1_gamma_E_corr_gate1383->Write();
    h1_gamma_E_corr_gate3199->Write();
    h2_gamma_E1_vs_E2_mult2->Write();
    h1_gamma_Esum_mult2->Write();
    h2_gamma_E_vs_theta_all->Write();
    h2_gamma_Ecorr_vs_theta_all->Write();
    h2_gamma_Ecorr_vs_theta_mleq2->Write();
    h2_gamma_Ecorr_vs_theta_mleq3->Write();

    // ── Doppler-correction QA: peak position vs theta ────────────────────────
    // If <base>_mean is FLAT vs theta, the correction is good and the peak width
    // is genuine detector resolution. A tilt => correction artefact broadening.
    // Corrected spectrum: expect FLAT.  Raw spectrum: expect a strong TILT
    // (that tilt is exactly the Doppler shift the correction is meant to remove).
    DopplerQA_FitLineVsTheta(h2_gamma_Ecorr_vs_theta_all,
                             "qa_corr_3199", 3.199, 2.40, 4.00, 10);
    DopplerQA_FitLineVsTheta(h2_gamma_Ecorr_vs_theta_all,
                             "qa_corr_1383", 1.383, 1.00, 1.80, 10);
    DopplerQA_FitLineVsTheta(h2_gamma_E_vs_theta_all,
                             "qa_raw_3199",  3.199, 2.00, 4.60, 10);

    TParameter<double> pFrag (GammaCfg::NFragName(),  static_cast<double>(nPassNeu));
    TParameter<double> pCoinc(GammaCfg::NCoincName(), static_cast<double>(nCoinc));
    pFrag.Write();
    pCoinc.Write();

    fout_root->Close();
    freac->Close();
    delete fNeu;

    std::cout << "[OK] Histograms saved to: " << outFile << "\n";
}

// ─────────────────────────────────────────────────────────────────────────────
//  Entry point (called by ROOT or compiled)
// ─────────────────────────────────────────────────────────────────────────────
#ifndef __CINT__
int main(int argc, char** argv)
{
    if (argc < 2) {
        std::cerr << "Usage: gammaSpectra <reactionFile.root>"
                     " [outFile.root] [unreactedFile.root] [offsetFile]"
                     " [kPID] [opaMin] [opaMax] [species] [neulandCut]\n";
        return 1;
    }

    TString reac  = argv[1];
    TString out   = (argc > 2) ? argv[2] : "gamma_spectra_22O.root";  // bare name → resolved inside gammaSpectra()
    TString unr   = (argc > 3) ? argv[3] : "";
    TString offs  = (argc > 4) ? argv[4] : "22O.txt";
    double  kPID  = (argc > 5) ? std::stod(argv[5]) : PID_K_DEFAULT;
    double  opaLo = (argc > 6) ? std::stod(argv[6]) : OPA_MIN;
    double  opaHi = (argc > 7) ? std::stod(argv[7]) : OPA_MAX;
    TString spec  = (argc > 8) ? argv[8] : "22O";
    TString neu   = (argc > 9) ? argv[9] : "";

    gammaSpectra(reac, out, unr, offs, kPID, opaLo, opaHi, spec, neu);
    return 0;
}
#endif
