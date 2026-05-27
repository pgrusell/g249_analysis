/** --------------------------------------------------------------------
 **
 **  run_sim_analysis.C
 **
 **  Post-processing of the simulation produced by runsim.C
 **  for the 25F(p,2p)24O experiment (G249).
 **
 **  Strategy (mirrors the experimental analysis chain in
 **  unpack_data_all_levels.C and R3BTrackingG249::Exec):
 **
 **    1) Read MC truth from glad.simu.root (event tree) and
 **       glad.truth.root (per-event SPO truth).
 **    2) Build a "FOOT-equivalent" outgoing track at the target by
 **       taking the MC primary vertex and initial momentum direction,
 **       and smearing it with the FOOT position resolution (75 um).
 **    3) For each X-fiber (Fi32, Fi33, Fi31) and Y-fiber (Fi30):
 **         - take the highest-energy MC point,
 **         - transform its lab (x,y,z) into the detector frame,
 **         - smear the strip coordinate by 150 um,
 **         - transform back to the lab frame.
 **    4) Pick the last X-fiber (Fi33 or Fi31) by highest fELoss.
 **    5) Build the 9-variable MDF input and evaluate PoQ.
 **    6) Charge: smeared TofD eloss-derived Z (sigma_Z = 0.2).
 **       A/Q: PoQ / (beta * gamma * AMU), using MC beta (parent).
 **    7) Recover 4-momentum in the lab frame, boost back by -beta along z
 **       and fill the rest-frame momentum distributions.
 **
 **  Usage:
 **      root -l 'run_sim_analysis.C("glad.simu.root",
 **                                  "glad.truth.root",
 **                                  "PoQ_9vars_150terms.txt",
 **                                  "sim_analysis.root")'
 **
 **/

// ----------------------------------------------------------------------------
// Constants and geometry (must match runsim.C / unpack_data_all_levels.C)
// ----------------------------------------------------------------------------
static constexpr Double_t kAMU = 0.9314940038; // GeV/c^2
static constexpr Double_t kAngleRad = 18.0 * TMath::DegToRad();

// Lab positions (cm). These are NON-const globals so they can be overridden
// from the steering call (see the optional last argument of run_sim_analysis).
//
// Defaults: the ones actually used to place the detectors in runsim.C
// (the *uncommented* block in that file).
//
// IMPORTANT CAVEAT: the experiment's tracking task R3BTrackingG249 in
// unpack_data_all_levels.C uses positions computed from Target_to_TP and
// TP_to_F* formulae, and these differ from the runsim.C placements by a few
// cm. The MDF function was trained against the experiment positions, so for
// a *closure test* you must use the same positions in both places. The
// cleanest fix is to make runsim.C match the experiment positions; the
// second-best (compatible) fix is to call run_sim_analysis() with the
// experiment positions via the kUseExperimentGeometry switch below.
static TVector3 gPosFi32(-132.754, 0.000, 696.144);
static TVector3 gPosFi30(-141.252, 0.000, 722.299);
static TVector3 gPosFi33(-138.490, 0.000, 794.696);
static TVector3 gPosFi31(-191.109, 0.000, 794.843);

// Helper: load the alternative positions used by R3BTrackingG249 in the
// experimental analysis (recomputed from the constants in
// unpack_data_all_levels.C). Use this if the MDF was trained with those.
static void UseExperimentGeometry()
{
    const Double_t A = kAngleRad;
    const Double_t C = TMath::Cos(A);
    const Double_t S = TMath::Sin(A);
    const Double_t TgF = 120.31;       // Target_to_GLAD_flange [cm]
    const Double_t T_TP = TgF + 165.4; // Target_to_TP
    const Double_t TP_F32 = 253.01834 + 175.;
    const Double_t TP_F30 = TP_F32 + 27.5;
    const Double_t TP_F31 = TP_F30 + 82.8;
    const Double_t TP_F33 = TP_F30 + 68.;
    const Double_t F33_off = 25.;
    const Double_t F31_off = -25.;

    gPosFi32.SetXYZ(0. - (TP_F32 * S) + 2.10016, 0., T_TP + TP_F32 * C + 2.56339);
    gPosFi30.SetXYZ(0. - (TP_F30 * S) - 4.35615, 0. - 1.4757,
                    T_TP + TP_F30 * C + 8.07945);
    gPosFi33.SetXYZ(F33_off * C - (TP_F33 * S) + 3.5479, 0.,
                    T_TP + TP_F33 * C + F33_off * S - 1.2248);
    gPosFi31.SetXYZ(F31_off * C - (TP_F31 * S) + 0.357029, 0.,
                    T_TP + TP_F31 * C + F31_off * S + 1.34728);

    std::cout << "[run_sim_analysis] Using EXPERIMENT geometry (matches MDF training):\n"
              << "    Fi32 = (" << gPosFi32.X() << "," << gPosFi32.Y() << "," << gPosFi32.Z() << ")\n"
              << "    Fi30 = (" << gPosFi30.X() << "," << gPosFi30.Y() << "," << gPosFi30.Z() << ")\n"
              << "    Fi33 = (" << gPosFi33.X() << "," << gPosFi33.Y() << "," << gPosFi33.Z() << ")\n"
              << "    Fi31 = (" << gPosFi31.X() << "," << gPosFi31.Y() << "," << gPosFi31.Z() << ")\n";
}

// In R3BTrackingG249 the experiment macro uses SetAnglesFib32(0, -Angle, 0),
// so only a rotation around the Y axis (by -18 deg) is applied to map local
// strip coordinates into the lab frame. We use exactly the same convention.
static const TVector3 kRotFi(0., -kAngleRad, 0.); // same for all four fibers

// Fragment identity (must match runsim.C and the rest of the chain)
static constexpr Int_t kFragZ = 8;
static constexpr Int_t kFragA = 24;
static const Double_t kFragMass = static_cast<Double_t>(kFragA) * kAMU; // GeV/c^2

// Smearing values (from the user request)
static constexpr Double_t kSigmaFootPos = 75.e-4;   // 75 um  -> cm
static constexpr Double_t kSigmaFiberPos = 150.e-4; // 150 um -> cm
static constexpr Double_t kSigmaTofdPos = 150.e-4;  // 150 um -> cm
static constexpr Double_t kSigmaTofdCharge = 0.2;   // absolute Z units

// GLAD currents (must match the experimental analysis exactly so that the
// MDF rescaling reproduces the experiment)
static constexpr Double_t kGladCurrent = 2668.0; // A (this simulation)
static constexpr Double_t kGladRefMDF = 2443.0;  // A (training current)

// ----------------------------------------------------------------------------
// Geometry helper -- same rotation convention as R3BTrackingG249::TransformPoint
// but specialized to (0, rotY, 0). Direction = +1 -> Det2Lab, -1 -> Lab2Det.
// ----------------------------------------------------------------------------
static inline TVector3 ToLab(const TVector3 &local, const TVector3 &pos)
{
    TRotation r;
    r.RotateY(kRotFi.Y());
    TVector3 v = local;
    v.Transform(r);
    v += pos;
    return v;
}

static inline TVector3 ToDet(const TVector3 &lab, const TVector3 &pos)
{
    TRotation r;
    r.RotateY(kRotFi.Y());
    TVector3 v = lab - pos;
    v.Transform(r.Inverse());
    return v;
}

// ----------------------------------------------------------------------------
// Get the highest-fELoss MC point from a (TraPoint-like) clones array for a
// given event. Returns (lab x,y,z) in cm and fills `valid`.
// We look at branches named "<prefix>fX", "<prefix>fY", "<prefix>fZ",
// "<prefix>fELoss" (split-mode reading via TTreeReaderArray).
// ----------------------------------------------------------------------------

// A small struct to carry the per-event reconstructed hit in the LAB frame.
struct FiberLabHit
{
    bool valid = false;
    TVector3 lab; // cm
    Double_t eloss = 0.;
};

// Pick the highest-eloss entry from parallel TTreeReaderArray<Double32_t>.
static FiberLabHit
PickMaxElossHit(TTreeReaderArray<Double32_t> &x,
                TTreeReaderArray<Double32_t> &y,
                TTreeReaderArray<Double32_t> &z,
                TTreeReaderArray<Double32_t> &eloss)
{
    FiberLabHit hit;
    const auto n = x.GetSize();
    if (n == 0)
        return hit;
    Int_t imax = -1;
    Double_t emax = -1.;
    for (size_t i = 0; i < n; ++i)
    {
        if (eloss[i] > emax)
        {
            emax = eloss[i];
            imax = static_cast<Int_t>(i);
        }
    }
    if (imax < 0)
        return hit;
    hit.valid = true;
    hit.lab.SetXYZ(x[imax], y[imax], z[imax]); // already in cm in MC points
    hit.eloss = emax;
    return hit;
}

// ----------------------------------------------------------------------------
// Extrapolate (X,Z) of the Fi30 hit from the X-track (Fi32 -> last X-fiber)
// intersected with the Fi30 plane. Mirrors R3BTrackingG249::FixupFiberTrack.
// f1 = Fi32 lab hit (X), f3 = last X-fiber lab hit (X).
// On entry f2.Y is the Fi30 measured Y in lab; on exit f2.{X,Z} are set.
// ----------------------------------------------------------------------------
static void FixupFib30(const TVector3 &f1, TVector3 &f2, const TVector3 &f3)
{
    // Two points defining the Fi30 plane in (X,Z), at Y=0 locally
    TVector3 e0(-1., 0., 0.);
    TVector3 e1(1., 0., 0.);
    e0 = ToLab(e0, gPosFi30);
    e1 = ToLab(e1, gPosFi30);
    const Double_t f30_slope = (e1.X() - e0.X()) / (e1.Z() - e0.Z());
    const Double_t f30_offset = e0.X() - f30_slope * e0.Z();

    const Double_t track_slope = (f3.X() - f1.X()) / (f3.Z() - f1.Z());
    const Double_t track_offset = f3.X() - track_slope * f3.Z();

    const Double_t z_new = (track_offset - f30_offset) / (f30_slope - track_slope);
    const Double_t x_new = track_slope * z_new + track_offset;
    f2.SetX(x_new);
    f2.SetZ(z_new);
}

// ----------------------------------------------------------------------------
// Main entry point
// ----------------------------------------------------------------------------
void ana(const char *simFile = "glad.simu.root",
         const char *truthFile = "glad.truth.root",
         const char *mdfFile = "PoQ_9vars_150terms.txt",
         const char *outFile = "sim_analysis.root",
         Long64_t maxEvents = -1,
         UInt_t randomSeed = 0,
         Bool_t useExperimentGeometry = kFALSE)
{
    TStopwatch timer;
    timer.Start();

    gRandom->SetSeed(randomSeed);

    if (useExperimentGeometry)
    {
        UseExperimentGeometry();
    }

    // ------------------------------------------------------------------------
    // 1) Open input files
    // ------------------------------------------------------------------------
    std::cout << "[run_sim_analysis] Opening input files\n"
              << "                   sim   : " << simFile << "\n"
              << "                   truth : " << truthFile << std::endl;

    TFile *fSim = TFile::Open(simFile, "READ");
    if (!fSim || fSim->IsZombie())
    {
        std::cerr << "ERROR: cannot open " << simFile << std::endl;
        return;
    }
    TTree *tSim = dynamic_cast<TTree *>(fSim->Get("evt"));
    if (!tSim)
    {
        std::cerr << "ERROR: tree 'evt' not found in " << simFile << std::endl;
        return;
    }

    TFile *fTruth = TFile::Open(truthFile, "READ");
    TTree *tTruth = nullptr;
    if (fTruth && !fTruth->IsZombie())
    {
        tTruth = dynamic_cast<TTree *>(fTruth->Get("truth"));
    }
    if (!tTruth)
    {
        std::cerr << "WARNING: tree 'truth' not found in " << truthFile
                  << " -- will fall back to MCTrack for the parent beta.\n";
    }

    // ------------------------------------------------------------------------
    // 2) Set up the MDF wrapper (same as R3BTrackingG249)
    // ------------------------------------------------------------------------
    std::cout << "[run_sim_analysis] Loading MDF: " << mdfFile << std::endl;
    R3BMDFWrapper *mdfPoQ = new R3BMDFWrapper(mdfFile);

    // ------------------------------------------------------------------------
    // 3) TTreeReader for the simulation tree.
    //    Branch names follow what evt->Print() showed for runsim.C output.
    //
    //    IMPORTANT: MCTrack has FOUR sub-branches with the same leaf name
    //    `fCoordinates.fX/fY/fZ` (one set for fStartVertex, one for
    //    fMomentumMass). TTreeReaderArray<Double_t> resolves these by leaf
    //    path, which is ambiguous and silently returns uninitialized memory
    //    for the wrong branch. We therefore do NOT use TTreeReaderArray for
    //    MCTrack and instead read the TClonesArray of R3BMCTrack objects
    //    directly -- which is what every R3B analysis task does.
    // ------------------------------------------------------------------------
    TTreeReader rdr(tSim);

    // Fiber MC points (these branches don't collide, TTreeReaderArray is safe)
    TTreeReaderArray<Double32_t> f30x(rdr, "Fi30Point.fX");
    TTreeReaderArray<Double32_t> f30y(rdr, "Fi30Point.fY");
    TTreeReaderArray<Double32_t> f30z(rdr, "Fi30Point.fZ");
    TTreeReaderArray<Double32_t> f30e(rdr, "Fi30Point.fELoss");

    TTreeReaderArray<Double32_t> f31x(rdr, "Fi31Point.fX");
    TTreeReaderArray<Double32_t> f31y(rdr, "Fi31Point.fY");
    TTreeReaderArray<Double32_t> f31z(rdr, "Fi31Point.fZ");
    TTreeReaderArray<Double32_t> f31e(rdr, "Fi31Point.fELoss");

    TTreeReaderArray<Double32_t> f32x(rdr, "Fi32Point.fX");
    TTreeReaderArray<Double32_t> f32y(rdr, "Fi32Point.fY");
    TTreeReaderArray<Double32_t> f32z(rdr, "Fi32Point.fZ");
    TTreeReaderArray<Double32_t> f32e(rdr, "Fi32Point.fELoss");

    TTreeReaderArray<Double32_t> f33x(rdr, "Fi33Point.fX");
    TTreeReaderArray<Double32_t> f33y(rdr, "Fi33Point.fY");
    TTreeReaderArray<Double32_t> f33z(rdr, "Fi33Point.fZ");
    TTreeReaderArray<Double32_t> f33e(rdr, "Fi33Point.fELoss");

    // TofD MC points
    TTreeReaderArray<Double32_t> tdEloss(rdr, "TofDPoint.fELoss");

    // ---- MCTrack via TClonesArray (avoids the fCoordinates.fX ambiguity) ---
    // We attach the MCTrack branch to a TClonesArray and read it through
    // R3BMCTrack accessors. We construct the TCA with a nullptr first and
    // let SetBranchAddress figure out the actual class from the file.
    TClonesArray *mcArr = nullptr;
    TBranch *mcBranch = tSim->GetBranch("MCTrack");
    if (!mcBranch)
    {
        std::cerr << "ERROR: MCTrack branch not found in " << simFile << std::endl;
        return;
    }
    tSim->SetBranchAddress("MCTrack", &mcArr);

    // ------------------------------------------------------------------------
    // 4) Truth tree readers (parallel-friend by entry index)
    // ------------------------------------------------------------------------
    Double_t truthBeta = -1.;
    Double_t truthPxRest = 0.;
    Double_t truthPyRest = 0.;
    Double_t truthPzRest = 0.;
    if (tTruth)
    {
        tTruth->SetBranchAddress("BetaParent", &truthBeta);
        tTruth->SetBranchAddress("PxRest", &truthPxRest);
        tTruth->SetBranchAddress("PyRest", &truthPyRest);
        tTruth->SetBranchAddress("PzRest", &truthPzRest);
    }

    // ------------------------------------------------------------------------
    // 5) Output histograms and tree
    // ------------------------------------------------------------------------
    TFile *fout = TFile::Open(outFile, "RECREATE");
    fout->cd();

    // ----- Reconstructed lab observables --------------------------------
    auto hPoQ = new TH1D("hPoQ", "Reconstructed P/Q;P/Q [GeV/c/e];counts", 600, -10.0, 20.0);
    auto hAoQ = new TH1D("hAoQ", "Reconstructed A/Q;A/Q;counts", 400, 1.5, 4.5);
    auto hZ = new TH1D("hZ", "TofD charge (smeared);Z;counts", 400, 5.0, 10.0);
    auto hPlab = new TH1D("hPlab", "Reconstructed |p|_{lab};|p|_{lab} [GeV/c];counts",
                          500, 0., 40.);

    // ----- Reconstructed rest-frame momentum (boost back along -z by MC beta)
    auto hPxRestRec = new TH1D("hPxRestRec", "Rec. p_{x} in parent rest frame;p_{x} [GeV/c];counts", 400, -0.5, 0.5);
    auto hPyRestRec = new TH1D("hPyRestRec", "Rec. p_{y} in parent rest frame;p_{y} [GeV/c];counts", 400, -0.5, 0.5);
    auto hPzRestRec = new TH1D("hPzRestRec", "Rec. p_{z} in parent rest frame;p_{z} [GeV/c];counts", 400, -1.0, 1.0);
    auto hPmagRestRec = new TH1D("hPmagRestRec",
                                 "Rec. |p| in parent rest frame;|p| [GeV/c];counts",
                                 400, 0., 1.0);

    // ----- Truth rest-frame momentum (from glad.truth.root) -------------
    auto hPxRestTru = new TH1D("hPxRestTru", "Truth p_{x} (rest);p_{x} [GeV/c];counts", 400, -0.5, 0.5);
    auto hPyRestTru = new TH1D("hPyRestTru", "Truth p_{y} (rest);p_{y} [GeV/c];counts", 400, -0.5, 0.5);
    auto hPzRestTru = new TH1D("hPzRestTru", "Truth p_{z} (rest);p_{z} [GeV/c];counts", 400, -1.0, 1.0);

    // ----- Residuals (rec - truth) --------------------------------------
    auto hResPx = new TH1D("hResPx", "p_{x}^{rec} - p_{x}^{truth};#Deltap_{x} [GeV/c];counts", 400, -0.3, 0.3);
    auto hResPy = new TH1D("hResPy", "p_{y}^{rec} - p_{y}^{truth};#Deltap_{y} [GeV/c];counts", 400, -0.3, 0.3);
    auto hResPz = new TH1D("hResPz", "p_{z}^{rec} - p_{z}^{truth};#Deltap_{z} [GeV/c];counts", 400, -0.5, 0.5);

    // ----- Correlation (truth vs rec) -----------------------------------
    auto h2Pz = new TH2D("h2Pz_truth_vs_rec",
                         "p_{z} rec vs truth (rest);p_{z}^{truth} [GeV/c];p_{z}^{rec} [GeV/c]",
                         200, -1.0, 1.0, 200, -1.0, 1.0);

    // Output ntuple with the per-event variables
    auto tout = new TTree("ana", "G249 simulation reconstruction");
    Double_t br_PoQ, br_AoQ, br_Z, br_beta, br_Plab;
    Double_t br_pxR_rec, br_pyR_rec, br_pzR_rec;
    Double_t br_pxR_tru, br_pyR_tru, br_pzR_tru;
    Double_t br_vx, br_vy, br_vz, br_tx, br_ty;
    tout->Branch("PoQ", &br_PoQ, "PoQ/D");
    tout->Branch("AoQ", &br_AoQ, "AoQ/D");
    tout->Branch("Z", &br_Z, "Z/D");
    tout->Branch("beta", &br_beta, "beta/D");
    tout->Branch("Plab", &br_Plab, "Plab/D");
    tout->Branch("PxRestRec", &br_pxR_rec, "PxRestRec/D");
    tout->Branch("PyRestRec", &br_pyR_rec, "PyRestRec/D");
    tout->Branch("PzRestRec", &br_pzR_rec, "PzRestRec/D");
    tout->Branch("PxRestTru", &br_pxR_tru, "PxRestTru/D");
    tout->Branch("PyRestTru", &br_pyR_tru, "PyRestTru/D");
    tout->Branch("PzRestTru", &br_pzR_tru, "PzRestTru/D");
    tout->Branch("Vx", &br_vx, "Vx/D");
    tout->Branch("Vy", &br_vy, "Vy/D");
    tout->Branch("Vz", &br_vz, "Vz/D");
    tout->Branch("TX", &br_tx, "TX/D");
    tout->Branch("TY", &br_ty, "TY/D");

    // Bookkeeping counters -- one per silent `continue` so we can diagnose
    // a zero-yield run without binary searching the macro.
    Long64_t nProcessed = 0;
    Long64_t nNoPrimary = 0;
    Long64_t nFiberMiss = 0;
    Long64_t nBadPoQ = 0;
    Long64_t nNoTofd = 0;
    Long64_t nBadBeta = 0;
    Long64_t nReconstructed = 0;
    // Track what kind of primary we ended up using, in case the PDG filter
    // didn't match anything (the fallback would otherwise be invisible).
    Long64_t nPrimaryByPdg = 0;
    Long64_t nPrimaryByFallback = 0;
    // First-event verbose dump: print sizes / values for the first few events
    // so the user can see immediately whether the branches are populated.
    constexpr Int_t kVerboseFirstN = 5;

    // ------------------------------------------------------------------------
    // 6) Event loop
    // ------------------------------------------------------------------------
    const Long64_t nEntries = (maxEvents < 0)
                                  ? tSim->GetEntries()
                                  : std::min<Long64_t>(maxEvents, tSim->GetEntries());

    std::cout << "[run_sim_analysis] Processing " << nEntries << " events\n";

    Long64_t iEv = -1;
    while (rdr.Next())
    {
        ++iEv;
        if (iEv >= nEntries)
            break;
        ++nProcessed;

        // Manually pull MCTrack for this entry (rdr.Next() doesn't touch it
        // because we attached MCTrack outside the TTreeReader).
        mcBranch->GetEntry(iEv);
        const Int_t nMC = mcArr ? mcArr->GetEntriesFast() : 0;

        if (tTruth)
        {
            tTruth->GetEntry(iEv);
        }
        else
        {
            truthBeta = -1.;
            truthPxRest = truthPyRest = truthPzRest = 0.;
        }

        // For the first few events, dump branch sizes so the user can see
        // immediately whether the input file has any content where we expect.
        if (iEv < kVerboseFirstN)
        {
            std::cout << "  [diag] ev=" << iEv
                      << "  MCTrack=" << nMC
                      << "  Fi32=" << f32x.GetSize()
                      << "  Fi30=" << f30x.GetSize()
                      << "  Fi33=" << f33x.GetSize()
                      << "  Fi31=" << f31x.GetSize()
                      << "  TofD=" << tdEloss.GetSize()
                      << "  truthBeta=" << truthBeta
                      << std::endl;
            if (nMC > 0)
            {
                auto *mc0 = dynamic_cast<R3BMCTrack *>(mcArr->At(0));
                if (mc0)
                {
                    std::cout << "         MCTrack[0] pdg=" << mc0->GetPdgCode()
                              << "  mother=" << mc0->GetMotherId()
                              << "  p=(" << mc0->GetPx() << "," << mc0->GetPy()
                              << "," << mc0->GetPz() << ")"
                              << "  v=(" << mc0->GetStartX() << "," << mc0->GetStartY()
                              << "," << mc0->GetStartZ() << ")"
                              << std::endl;
                }
            }
        }

        // -------- Find the primary fragment in MCTrack ----------------------
        // The PDG of a fully-stripped 24O is 1000000000 + 8*10000 + 24*10 = 1000080240.
        // In FairStack different versions use different "no mother" sentinels;
        // we accept anything that looks unambiguously primary (motherId<=0 OR
        // first entry) together with a PDG match.
        const Int_t expectedPdg = 1000000000 + kFragZ * 10000 + kFragA * 10;
        Int_t iPrimary = -1;
        // Pass 1: explicit primary (mother<=0) + PDG match
        for (Int_t i = 0; i < nMC; ++i)
        {
            auto *m = dynamic_cast<R3BMCTrack *>(mcArr->At(i));
            if (!m)
                continue;
            if (m->GetMotherId() <= 0 && m->GetPdgCode() == expectedPdg)
            {
                iPrimary = i;
                ++nPrimaryByPdg;
                break;
            }
        }
        // Pass 2: any PDG match
        if (iPrimary < 0)
        {
            for (Int_t i = 0; i < nMC; ++i)
            {
                auto *m = dynamic_cast<R3BMCTrack *>(mcArr->At(i));
                if (!m)
                    continue;
                if (m->GetPdgCode() == expectedPdg)
                {
                    iPrimary = i;
                    ++nPrimaryByPdg;
                    break;
                }
            }
        }
        // Pass 3: fall back to MCTrack[0] -- almost always the primary in
        // FairRoot, regardless of the motherId convention.
        if (iPrimary < 0 && nMC > 0)
        {
            iPrimary = 0;
            ++nPrimaryByFallback;
            if (iEv < kVerboseFirstN)
            {
                auto *m0 = dynamic_cast<R3BMCTrack *>(mcArr->At(0));
                std::cout << "  [diag] ev=" << iEv
                          << ": no PDG=" << expectedPdg << " in MCTrack, using index 0 (pdg="
                          << (m0 ? m0->GetPdgCode() : 0) << ")\n";
            }
        }
        if (iPrimary < 0)
        {
            ++nNoPrimary;
            continue;
        }

        // True primary vertex (cm) and momentum (GeV/c) in the lab
        auto *mcPrim = dynamic_cast<R3BMCTrack *>(mcArr->At(iPrimary));
        if (!mcPrim)
        {
            ++nNoPrimary;
            continue;
        }
        const Double_t vx_true = mcPrim->GetStartX();
        const Double_t vy_true = mcPrim->GetStartY();
        const Double_t vz_true = mcPrim->GetStartZ();
        const Double_t px_true = mcPrim->GetPx();
        const Double_t py_true = mcPrim->GetPy();
        const Double_t pz_true = mcPrim->GetPz();
        const Double_t p_true = std::sqrt(px_true * px_true + py_true * py_true + pz_true * pz_true);

        // -------- Build a "FOOT-equivalent" outgoing track ------------------
        // The experiment's MDF uses (x,y,z,tx,ty) at the target as the seed
        // (these come from the FOOT tracker, which we don't simulate here).
        // We smear the true vertex (x,y) and direction (tx,ty) with the FOOT
        // position resolution, treating it as a 4-point straight track.
        //
        // For a 4-strip FOOT layout with sigma_pos = 75 um and roughly L_FOOT
        // baseline, the angular resolution would be sigma_pos*sqrt(2)/L.
        // Without a concrete FOOT geometry in the sim we adopt the simpler,
        // physically motivated choice: position smeared by 75 um, slope
        // smeared by 75 um / L_baseline. We take L_baseline = 30 cm as a
        // representative FOOT footprint (override below if you know better).
        constexpr Double_t kFootBaseline = 30.; // cm, FOOT footprint
        const Double_t sigma_tx = kSigmaFootPos / kFootBaseline;
        const Double_t sigma_ty = kSigmaFootPos / kFootBaseline;

        Double_t vx_smear = vx_true + gRandom->Gaus(0., kSigmaFootPos);
        Double_t vy_smear = vy_true + gRandom->Gaus(0., kSigmaFootPos);
        Double_t vz_smear = vz_true; // FOOT extrapolated to target z

        Double_t tx_true = px_true / pz_true;
        Double_t ty_true = py_true / pz_true;
        Double_t tx_smear = tx_true + gRandom->Gaus(0., sigma_tx);
        Double_t ty_smear = ty_true + gRandom->Gaus(0., sigma_ty);

        // -------- Fiber hits: highest-eloss point, smear local coord -------
        FiberLabHit f32 = PickMaxElossHit(f32x, f32y, f32z, f32e);
        FiberLabHit f30 = PickMaxElossHit(f30x, f30y, f30z, f30e);
        FiberLabHit f31 = PickMaxElossHit(f31x, f31y, f31z, f31e);
        FiberLabHit f33 = PickMaxElossHit(f33x, f33y, f33z, f33e);

        // Need at least Fi32 and Fi30, AND at least one of {Fi31, Fi33}
        if (!f32.valid || !f30.valid || (!f31.valid && !f33.valid))
        {
            ++nFiberMiss;
            continue;
        }

        // -- Smear Fi32 (X-fiber): local X is the measured strip
        {
            TVector3 d = ToDet(f32.lab, gPosFi32);
            d.SetX(d.X() + gRandom->Gaus(0., kSigmaFiberPos));
            // X-fibers measure only X; the others (Y,Z in det frame) are
            // not physical, so we drop them by setting them to zero before
            // transforming back -- the lab hit is then the strip in lab.
            d.SetY(0.);
            d.SetZ(0.);
            f32.lab = ToLab(d, gPosFi32);
        }

        // -- Smear Fi30 (Y-fiber): keep only Y
        {
            TVector3 d = ToDet(f30.lab, gPosFi30);
            d.SetY(d.Y() + gRandom->Gaus(0., kSigmaFiberPos));
            d.SetX(0.);
            d.SetZ(0.);
            f30.lab = ToLab(d, gPosFi30);
        }

        // -- Smear Fi31 (X-fiber)
        if (f31.valid)
        {
            TVector3 d = ToDet(f31.lab, gPosFi31);
            d.SetX(d.X() + gRandom->Gaus(0., kSigmaFiberPos));
            d.SetY(0.);
            d.SetZ(0.);
            f31.lab = ToLab(d, gPosFi31);
        }

        // -- Smear Fi33 (X-fiber)
        if (f33.valid)
        {
            TVector3 d = ToDet(f33.lab, gPosFi33);
            d.SetX(d.X() + gRandom->Gaus(0., kSigmaFiberPos));
            d.SetY(0.);
            d.SetZ(0.);
            f33.lab = ToLab(d, gPosFi33);
        }

        // -------- Pick the last X-fiber by highest eloss (mirrors R3BTrackingG249)
        TVector3 fLast_lab;
        bool useFi33 = false;
        if (f33.valid && f31.valid)
        {
            useFi33 = (f33.eloss >= f31.eloss);
        }
        else
        {
            useFi33 = f33.valid;
        }
        fLast_lab = useFi33 ? f33.lab : f31.lab;

        // -------- Fixup Fi30: extrapolate its lab (X,Z) from the X-track ----
        TVector3 fFi30_lab = f30.lab; // Y is set; X,Z will be overwritten
        FixupFib30(f32.lab, fFi30_lab, fLast_lab);

        // -------- Build MDF input (exactly as R3BTrackingG249::Exec) --------
        // Vertex and slopes must be in cm. Smeared values.
        const Double_t startX = vx_smear;
        const Double_t startY = vy_smear;
        const Double_t startZ = vz_smear;
        const Double_t startTX = tx_smear;
        const Double_t startTY = ty_smear;

        Double_t mdf_data[9];
        int k = 0;
        mdf_data[k] = startX;
        mdf_data[++k] = startY;
        mdf_data[++k] = startZ;
        mdf_data[++k] = startTX;
        mdf_data[++k] = startTY;
        mdf_data[++k] = f32.lab.X();
        mdf_data[++k] = f32.lab.Z();
        mdf_data[++k] = (fLast_lab.X() - f32.lab.X()) / (fLast_lab.Z() - f32.lab.Z());
        mdf_data[++k] = (fFi30_lab.Y() - startY) / (fFi30_lab.Z() - startZ);

        // ====================================================================
        // Diagnostic: check that no MDF input is NaN/Inf before calling the
        // wrapper. If something is bad, report which one and skip the event.
        // This is mainly here to localize bugs; once everything is sane it
        // never fires.
        // ====================================================================
        bool inputBad = false;
        int badIdx = -1;
        for (int kk = 0; kk < 9; ++kk)
        {
            if (!std::isfinite(mdf_data[kk]))
            {
                inputBad = true;
                badIdx = kk;
                break;
            }
        }
        if (inputBad)
        {
            ++nBadPoQ;
            if (iEv < kVerboseFirstN || nBadPoQ <= 5)
            {
                static const char *names[9] = {
                    "Vx", "Vy", "Vz", "TX", "TY",
                    "f32.X", "f32.Z", "slopeX", "slopeY"};
                std::cerr << "  [diag] ev=" << iEv
                          << " : MDF input #" << badIdx
                          << " (" << names[badIdx] << ") is non-finite ("
                          << mdf_data[badIdx] << ")."
                          << "  Full inputs: [";
                for (int kk = 0; kk < 9; ++kk)
                    std::cerr << mdf_data[kk] << (kk < 8 ? "," : "]");
                std::cerr << "\n         hits: vx_true=" << vx_true
                          << " vy_true=" << vy_true
                          << " vz_true=" << vz_true
                          << " px=" << px_true
                          << " py=" << py_true
                          << " pz=" << pz_true
                          << " f32.lab=(" << f32.lab.X() << "," << f32.lab.Y() << "," << f32.lab.Z() << ")"
                          << " fLast.lab=(" << fLast_lab.X() << "," << fLast_lab.Y() << "," << fLast_lab.Z() << ")"
                          << " fFi30.lab=(" << fFi30_lab.X() << "," << fFi30_lab.Y() << "," << fFi30_lab.Z() << ")"
                          << "\n";
            }
            continue;
        }

        Double_t PoQ_raw = mdfPoQ->MDF(mdf_data);
        Double_t PoQ = PoQ_raw * (kGladCurrent / kGladRefMDF);

        // For the first few events, print MDF inputs and outputs so the user
        // can see whether the polynomial is being driven into a bad domain.
        if (iEv < kVerboseFirstN)
        {
            std::cout << "  [diag] ev=" << iEv
                      << " PoQ_raw=" << PoQ_raw << " PoQ=" << PoQ
                      << " mdf_data=[";
            for (int kk = 0; kk < 9; ++kk)
                std::cout << mdf_data[kk] << (kk < 8 ? "," : "]");
            std::cout << std::endl;
        }

        // Only veto truly non-finite results (NaN, Inf). A finite but negative
        // PoQ means the MDF is being driven outside its training domain --
        // we still want to see that in the output histograms so the user can
        // diagnose, so we keep the event and let the histograms show the tail.
        if (!std::isfinite(PoQ))
        {
            ++nBadPoQ;
            if (nBadPoQ <= 5)
            {
                std::cerr << "  [diag] ev=" << iEv
                          << " : MDF returned non-finite PoQ from finite inputs."
                          << " This usually means the MDF coefficient file failed"
                          << " to load or the polynomial overflows.\n";
            }
            continue;
        }

        // -------- TofD charge: smear and average all hits -------------------
        // In the experiment R3BTofDCal2Hit gives an effective Z per hit; here
        // we just take Z_true (=kFragZ) and smear by sigma=0.2.
        //
        // Do NOT veto the event if there are no TofD points -- the MDF has
        // already given us PoQ from the upstream tracking, and the charge
        // only enters the |p|_lab calculation. We just flag it.
        Double_t Z_smear;
        if (tdEloss.GetSize() == 0)
        {
            ++nNoTofd;
            // Use the nominal Z (still smeared) so we don't lose the event.
            Z_smear = static_cast<Double_t>(kFragZ) + gRandom->Gaus(0., kSigmaTofdCharge);
        }
        else
        {
            Z_smear = static_cast<Double_t>(kFragZ) + gRandom->Gaus(0., kSigmaTofdCharge);
        }

        // -------- Parent beta: take MC truth (BetaParent if available) -----
        // This mirrors using FRS beta in the experimental analysis.
        Double_t beta_parent = truthBeta;
        if (beta_parent <= 0. || beta_parent >= 1.)
        {
            // Fallback: derive from the MC fragment itself (since the rest-
            // frame momentum is small compared to the boost, this is close)
            const Double_t E_true = std::sqrt(p_true * p_true + kFragMass * kFragMass);
            beta_parent = pz_true / E_true;
        }
        if (beta_parent <= 0. || beta_parent >= 1.)
        {
            ++nBadBeta;
            continue;
        }
        const Double_t gamma_parent = 1. / std::sqrt(1. - beta_parent * beta_parent);

        // -------- A/Q from PoQ, beta, gamma -------------------------------
        // Use |PoQ| for derived quantities (a sign flip from the MDF being
        // driven outside its training domain shouldn't kill the boost). The
        // raw signed value is still stored in hPoQ / the ntuple.
        const Double_t PoQ_abs = std::fabs(PoQ);
        const Double_t AoQ = PoQ_abs / (beta_parent * gamma_parent * kAMU);

        // -------- Build the reconstructed lab 4-momentum --------------------
        // The MDF gives only |p|/Q. We take the direction from the smeared
        // outgoing-track at the target ( (tx, ty, 1) normalised ) -- same as
        // R3BTrackingG249, which calls fragment_momentum.SetMag(PoQ) on the
        // outgoing-FOOT momentum vector.
        TVector3 pdir(startTX, startTY, 1.);
        pdir = pdir.Unit();
        const Double_t Q_eff = Z_smear;           // assume fully stripped
        const Double_t p_recon = PoQ_abs * Q_eff; // GeV/c
        TVector3 p3_lab = pdir * p_recon;
        const Double_t E_lab = std::sqrt(p_recon * p_recon + kFragMass * kFragMass);
        TLorentzVector p4_lab(p3_lab, E_lab);

        // -------- Boost back to the parent rest frame ----------------------
        TLorentzVector p4_rest = p4_lab;
        p4_rest.Boost(0., 0., -beta_parent);

        // -------- Fill histograms ------------------------------------------
        hPoQ->Fill(PoQ);
        hAoQ->Fill(AoQ);
        hZ->Fill(Z_smear);
        hPlab->Fill(p_recon);

        hPxRestRec->Fill(p4_rest.Px());
        hPyRestRec->Fill(p4_rest.Py());
        hPzRestRec->Fill(p4_rest.Pz());
        hPmagRestRec->Fill(p4_rest.Vect().Mag());

        if (tTruth)
        {
            hPxRestTru->Fill(truthPxRest);
            hPyRestTru->Fill(truthPyRest);
            hPzRestTru->Fill(truthPzRest);
            hResPx->Fill(p4_rest.Px() - truthPxRest);
            hResPy->Fill(p4_rest.Py() - truthPyRest);
            hResPz->Fill(p4_rest.Pz() - truthPzRest);
            h2Pz->Fill(truthPzRest, p4_rest.Pz());
        }

        br_PoQ = PoQ;
        br_AoQ = AoQ;
        br_Z = Z_smear;
        br_beta = beta_parent;
        br_Plab = p_recon;
        br_pxR_rec = p4_rest.Px();
        br_pyR_rec = p4_rest.Py();
        br_pzR_rec = p4_rest.Pz();
        br_pxR_tru = truthPxRest;
        br_pyR_tru = truthPyRest;
        br_pzR_tru = truthPzRest;
        br_vx = startX;
        br_vy = startY;
        br_vz = startZ;
        br_tx = startTX;
        br_ty = startTY;
        tout->Fill();

        ++nReconstructed;

        if (nReconstructed <= 3)
        {
            std::cout << "  event " << iEv
                      << " : PoQ=" << PoQ
                      << " AoQ=" << AoQ
                      << " beta=" << beta_parent
                      << " |p|_lab=" << p_recon
                      << " p_rest=(" << p4_rest.Px() << "," << p4_rest.Py()
                      << "," << p4_rest.Pz() << ")" << std::endl;
        }
    } // event loop

    // ------------------------------------------------------------------------
    // 7) Write and close
    // ------------------------------------------------------------------------
    fout->cd();
    tout->Write();
    hPoQ->Write();
    hAoQ->Write();
    hZ->Write();
    hPlab->Write();
    hPxRestRec->Write();
    hPyRestRec->Write();
    hPzRestRec->Write();
    hPmagRestRec->Write();
    hPxRestTru->Write();
    hPyRestTru->Write();
    hPzRestTru->Write();
    hResPx->Write();
    hResPy->Write();
    hResPz->Write();
    h2Pz->Write();
    fout->Close();

    fSim->Close();
    if (fTruth)
        fTruth->Close();

    timer.Stop();
    std::cout << "\n[run_sim_analysis] Done.\n"
              << "  events processed      : " << nProcessed << "\n"
              << "  primary by PDG match  : " << nPrimaryByPdg << "\n"
              << "  primary by fallback   : " << nPrimaryByFallback << "\n"
              << "  -- gates that vetoed events --\n"
              << "  no primary fragment   : " << nNoPrimary << "\n"
              << "  fiber multiplicity    : " << nFiberMiss << "\n"
              << "  MDF returned bad PoQ  : " << nBadPoQ << "\n"
              << "  bad parent beta       : " << nBadBeta << "\n"
              << "  -- non-vetoing flags --\n"
              << "  no TofD hit (Z=8 used): " << nNoTofd << "\n"
              << "  -- result --\n"
              << "  reconstructed         : " << nReconstructed << "\n"
              << "  output                : " << outFile << "\n"
              << "  real time             : " << timer.RealTime() << " s\n"
              << "  CPU time              : " << timer.CpuTime() << " s\n";

    delete mdfPoQ;
}