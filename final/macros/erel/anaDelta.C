/** --------------------------------------------------------------------
 **
 **  anaDelta.C   (corrected & extended)
 **
 **  Post-processing of the simulation produced by runsim.C for the
 **  in-flight breakup 24O -> 23O + n (G249): detector response in Erel.
 **
 **  Changes with respect to the previous version:
 **
 **   (a) FRAGMENT MOMENTUM NOW MIRRORS THE EXPERIMENT EXACTLY.
 **       In eventFilter.C the fragment momentum is p_frag = (P/Q) * Z_FRAG
 **       with Z_FRAG = 8 *fixed* (the ToFD charge is only used for PID,
 **       never to scale the momentum). The old macro multiplied PoQ by the
 **       *smeared* ToFD charge, adding a spurious ~2.5% momentum smearing
 **       that does not exist in the experimental chain. Z_smear is kept
 **       only as a PID observable (hZ / branch Z).
 **
 **   (b) FRAGMENT BETA FROM THE MASS HYPOTHESIS (as in the experiment).
 **       beta_frag = p_frag / sqrt(p_frag^2 + M_frag^2) with
 **       M_frag = (23.015696686 - 8*m_e) * AMU, i.e. exactly the
 **       25F23O ReactionConfig of eventFilter.C and the beta used by
 **       DataAnalysis::getData for the Erel. (The old AoQ diagnostic,
 **       which uses the MC-truth beta as a stand-in for the missing TOF,
 **       is kept unchanged -- it is only a closure/PID check.)
 **
 **   (c) NEULAND RESOLUTIONS ADDED (x/y along the bar + time), matching
 **       what the digitizer provides in the experiment:
 **         - time:            Gaussian, sigma = kSigmaNeulandTime (150 ps)
 **         - along-bar coord: Gaussian, sigma = kSigmaNeulandAlong (1.5 cm)
 **           (X for horizontal paddles, Y for vertical ones)
 **         - transverse and Z: uniform +-2.5 cm, exactly as
 **           eventFilter.C::defineNeutronColumns() applies on top of the
 **           digitized hits.
 **       The 80 ns TOF cut (together with an ELoss > 6 MeV requirement) is
 **       applied to the *smeared* time and beta_n is computed from the
 **       smeared flight path and smeared time, as in the experimental
 **       chain.
 **
 **   (d) EREL COMPUTED EXACTLY AS DataAnalysis::getData():
 **       slope vectors (px/pz, py/pz, 1) for neutron and fragment (offsets
 **       and beta-match = 0 in the simulation), cos(angle) from those,
 **       Erel = sqrt(mf^2 + mn^2 + 2 g_n g_f mf mn (1 - b_n b_f cos)) -
 **              mf - mn.
 **       Also delta_beta = beta_frag - beta_neu and the lab opening angle.
 **
 **   (e) TRUTH EREL from the two primaries (23O + n) is stored per event,
 **       so the output directly gives the detector *response*:
 **       hErel (reconstructed), hErelTrue, hErelRes (rec - true).
 **
 **   (f) ScanResolutions(): nominal sigmas now really match the defaults
 **       (the old ones silently disagreed), and the scanned histogram is
 **       selectable (default now "hErel", since that is the observable of
 **       this study).
 **
 **  Usage:
 **      root -l 'anaDelta.C("delta.simu.root",
 **                           "PoQ_9vars_150terms.txt",
 **                           "sim_analysis.root")'
 **
 **      The first argument may list several input files, comma- and/or
 **      space-separated, each of which may itself use wildcards; all of
 **      them are chained together on the "evt" tree, e.g.:
 **      root -l 'anaDelta.C("range.simu.root,range.simu_1.root",
 **                           "PoQ_9vars_150terms.txt",
 **                           "sim_analysis.root")'
 **      root -l 'anaDelta.C("sim_results/range.simu*.root", ...)'
 **
 **/

// ----------------------------------------------------------------------------
// Constants and geometry (must match runsim.C / unpack_data_all_levels.C)
// ----------------------------------------------------------------------------
static constexpr Double_t kAMU = 0.9314940038; // GeV/c^2 (legacy, AoQ diag only)
static constexpr Double_t kAngleRad = 18.0 * TMath::DegToRad();

// Constants used by the experimental analysis chain (eventFilter.C /
// DataAnalysis). Use THESE for everything entering the Erel so that the
// simulated response is computed with the very same numbers as the data.
static constexpr Double_t kAMU_EF = 0.93149410242;                    // GeV/c^2
static constexpr Double_t kFragMassAMU = 23.015696686 - 8. * 0.00511; // 23O, e- stripped
static constexpr Double_t kMassNeutron = 0.939565420;                 // GeV/c^2

// Lab positions (cm). Non-const so they can be overridden (see
// UseExperimentGeometry / the useExperimentGeometry switch of anaDelta).
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

    std::cout << "[anaDelta] Using EXPERIMENT geometry (matches MDF training):\n"
              << "    Fi32 = (" << gPosFi32.X() << "," << gPosFi32.Y() << "," << gPosFi32.Z() << ")\n"
              << "    Fi30 = (" << gPosFi30.X() << "," << gPosFi30.Y() << "," << gPosFi30.Z() << ")\n"
              << "    Fi33 = (" << gPosFi33.X() << "," << gPosFi33.Y() << "," << gPosFi33.Z() << ")\n"
              << "    Fi31 = (" << gPosFi31.X() << "," << gPosFi31.Y() << "," << gPosFi31.Z() << ")\n";
}

// Rotation convention: only RotateY(-18 deg), as in R3BTrackingG249.
static const TVector3 kRotFi(0., -kAngleRad, 0.);

// Fragment identity: in delta.simu.root the generator injects the breakup
// products directly: MCTrack[0] = 23O, MCTrack[1] = neutron.
static constexpr Int_t kFragZ = 8;
static constexpr Int_t kFragA = 23;

// ---- Smearing values (mutable so ScanResolutions can vary them) -----------
static Double_t kSigmaFootPos = 75.e-4;   // 75 um  -> cm
static Double_t kSigmaFiberPos = 150.e-4; // 150 um -> cm
static Double_t kSigmaTofdCharge = 0.2;   // absolute Z units (PID only!)

// ---- NeuLAND ----------------------------------------------------------------
// Transverse/Z paddle reassignment: uniform +-2.5 cm, exactly what
// eventFilter.C adds on top of the digitized hits.
static Double_t kSigmaNeulandPaddle = 2.5; // cm (uniform half-width)
// What the digitizer itself provides in the experiment and the raw MC
// points do NOT contain -- must therefore be emulated here:
static Double_t kSigmaNeulandAlong = 1.5;  // cm, along-bar coord (from dt)
static Double_t kSigmaNeulandTime = 0.150; // ns, paddle time resolution

// static Double_t kSigmaNeulandTime = 0.; // ns, paddle time resolution

static constexpr Double_t kNeuTofMaxNs = 80.0;

static constexpr Double_t kNeuElossMinGeV = 6.e-3; // 6 MeV energy-loss threshold
static constexpr Double_t kCcmPerNs = 29.9792458;

// GLAD currents (must match the experimental analysis)
static constexpr Double_t kGladCurrent = 2668.0; // A (this simulation)
static constexpr Double_t kGladRefMDF = 2443.0;  // A (training current)

// ----------------------------------------------------------------------------
// Geometry helpers -- same convention as R3BTrackingG249::TransformPoint
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
// Fiber MC-point helpers
// ----------------------------------------------------------------------------
struct FiberLabHit
{
    bool valid = false;
    TVector3 lab; // cm
    Double_t eloss = 0.;
};

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
    hit.lab.SetXYZ(x[imax], y[imax], z[imax]); // already cm in MC points
    hit.eloss = emax;
    return hit;
}

// ----------------------------------------------------------------------------
// Fi30 X/Z fixup from the X-track -- mirrors R3BTrackingG249::FixupFiberTrack
// ----------------------------------------------------------------------------
static void FixupFib30(const TVector3 &f1, TVector3 &f2, const TVector3 &f3)
{
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
// Neutron (first hit): eventFilter.C::defineNeutronColumns() applied to the
// raw NeulandPoints, but now with a full emulation of the hit-level
// resolutions the digitizer would have provided:
//
//   * time smeared by kSigmaNeulandTime (Gaussian); the t < 80 ns cut is
//     applied to the SMEARED time (as it is on measured hits in the data);
//   * energy loss must exceed kNeuElossMinGeV (6 MeV), as in the experimental
//     neutron selection;
//   * along-bar coordinate smeared by kSigmaNeulandAlong (Gaussian):
//       horizontal paddles ((paddleId/50)%2 == 0): X is along the bar,
//       vertical paddles:                          Y is along the bar;
//   * transverse coordinate and Z reassigned uniformly within +-2.5 cm,
//     exactly the smearing eventFilter.C adds on top of R3BNeulandHit.
//
// Candidate selection ("first" hit): smallest Z among points passing the
// (smeared-)time and energy-loss cuts, as in the experimental analysis.
// ----------------------------------------------------------------------------
struct NeutronHit
{
    bool valid = false;
    Double_t tns = -1.; // smeared TOF [ns]
    TVector3 pos;       // cm, lab frame, fully smeared
};

static NeutronHit
FindFirstNeutron(TTreeReaderArray<Double32_t> &x,
                 TTreeReaderArray<Double32_t> &y,
                 TTreeReaderArray<Double32_t> &z,
                 TTreeReaderArray<Double32_t> &t,
                 TTreeReaderArray<Double32_t> &eloss,
                 TTreeReaderArray<Int_t> &paddle)
{
    NeutronHit best;
    Double_t bestZ = 1.e99;
    const auto n = x.GetSize();
    for (size_t i = 0; i < n; ++i)
    {
        const Double_t tSmear = t[i] + gRandom->Gaus(0., kSigmaNeulandTime);

        if (eloss[i] <= kNeuElossMinGeV)
            continue;
        if (z[i] >= bestZ)
            continue;

        bestZ = z[i];
        best.valid = true;
        best.tns = tSmear;

        const bool horizontal = (((paddle[i] / 50) % 2) == 0);
        if (horizontal)
        {
            // bar axis = X (timing coordinate), paddle stacking = Y/Z
            best.pos.SetXYZ(
                x[i] + gRandom->Gaus(0., kSigmaNeulandAlong),
                y[i] + gRandom->Uniform(-kSigmaNeulandPaddle, kSigmaNeulandPaddle),
                z[i] + gRandom->Uniform(-kSigmaNeulandPaddle, kSigmaNeulandPaddle));
        }
        else
        {
            // bar axis = Y, paddle stacking = X/Z
            best.pos.SetXYZ(
                x[i] + gRandom->Uniform(-kSigmaNeulandPaddle, kSigmaNeulandPaddle),
                y[i] + gRandom->Gaus(0., kSigmaNeulandAlong),
                z[i] + gRandom->Uniform(-kSigmaNeulandPaddle, kSigmaNeulandPaddle));
        }
    }
    return best;
}

// ----------------------------------------------------------------------------
// Main entry point
// ----------------------------------------------------------------------------
void anaDelta(const char *simFile = "sim_results/bg.simu.root",
              const char *mdfFile = "PoQ_9vars_150terms.txt",
              const char *outFile = "ana_results/bg_analysis.root",
              Long64_t maxEvents = 100000000,
              UInt_t randomSeed = 0,
              Bool_t useExperimentGeometry = kFALSE,
              Double_t betaFragOffset = 0.01368 - 0.00730334,
              Double_t fragPxOffset = 0.0,
              Double_t fragPyOffset = 0.0,
              Double_t neuPxOffset = 0.0,
              Double_t neuPyOffset = 0.0)
{
    TStopwatch timer;
    timer.Start();

    gRandom->SetSeed(randomSeed);

    if (useExperimentGeometry)
        UseExperimentGeometry();

    // ------------------------------------------------------------------------
    // 1) Open input file(s)
    //    'simFile' may list several files (comma- and/or space-separated),
    //    each of which may itself contain wildcards, e.g.:
    //      "range.simu.root,range.simu_1.root"
    //      "sim_results/range.simu*.root"
    //    All matching files are chained together on the "evt" tree.
    // ------------------------------------------------------------------------
    std::cout << "[anaDelta] Opening input file(s)\n";

    TChain *tSim = new TChain("evt");
    {
        TString fileList(simFile);
        std::unique_ptr<TObjArray> tokens(fileList.Tokenize(", "));
        Int_t nFilesAdded = 0;
        for (Int_t i = 0; i < tokens->GetEntriesFast(); ++i)
        {
            TString one = static_cast<TObjString *>(tokens->At(i))->GetString();
            one = one.Strip(TString::kBoth);
            if (one.IsNull())
                continue;
            const Int_t added = tSim->Add(one.Data());
            std::cout << "           sim   : " << one << " (" << added << " file(s) matched)" << std::endl;
            nFilesAdded += added;
        }
        if (nFilesAdded == 0)
        {
            std::cerr << "ERROR: no input files found matching " << simFile << std::endl;
            delete tSim;
            return;
        }
    }

    // ------------------------------------------------------------------------
    // 2) MDF wrapper (same as R3BTrackingG249)
    // ------------------------------------------------------------------------
    std::cout << "[anaDelta] Loading MDF: " << mdfFile << std::endl;
    R3BMDFWrapper *mdfPoQ = new R3BMDFWrapper(mdfFile);

    // ------------------------------------------------------------------------
    // 3) TTreeReader for the simulation tree
    //    (MCTrack read via TClonesArray to avoid the fCoordinates.fX
    //    leaf-name ambiguity -- see comment in the previous version.)
    // ------------------------------------------------------------------------
    TTreeReader rdr(tSim);

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

    TTreeReaderArray<Double32_t> tdEloss(rdr, "TofDPoint.fELoss");

    // FOOT-equivalent target trackers (TraPoint, DetId 5..8)
    TTreeReaderArray<Double32_t> trX(rdr, "TraPoint.fX");
    TTreeReaderArray<Double32_t> trY(rdr, "TraPoint.fY");
    TTreeReaderArray<Double32_t> trZ(rdr, "TraPoint.fZ");
    TTreeReaderArray<Double32_t> trE(rdr, "TraPoint.fELoss");
    TTreeReaderArray<Int_t> trDet(rdr, "TraPoint.fDetectorID");

    // Read via TTreeReaderValue (rather than a raw TBranch + SetBranchAddress)
    // so that branch reattachment across chain links (multiple input files)
    // is handled automatically.
    TTreeReaderValue<TClonesArray> mcTrackVal(rdr, "MCTrack");

    const Bool_t hasNeuland = (tSim->GetBranch("NeulandPoints") != nullptr);
    std::unique_ptr<TTreeReaderArray<Double32_t>> nlX, nlY, nlZ, nlT, nlE;
    std::unique_ptr<TTreeReaderArray<Int_t>> nlDet;
    if (hasNeuland)
    {
        nlX = std::make_unique<TTreeReaderArray<Double32_t>>(rdr, "NeulandPoints.fX");
        nlY = std::make_unique<TTreeReaderArray<Double32_t>>(rdr, "NeulandPoints.fY");
        nlZ = std::make_unique<TTreeReaderArray<Double32_t>>(rdr, "NeulandPoints.fZ");
        nlT = std::make_unique<TTreeReaderArray<Double32_t>>(rdr, "NeulandPoints.fTime");
        nlE = std::make_unique<TTreeReaderArray<Double32_t>>(rdr, "NeulandPoints.fELoss");
        nlDet = std::make_unique<TTreeReaderArray<Int_t>>(rdr, "NeulandPoints.fDetectorID");
    }
    else
    {
        std::cout << "[anaDelta] WARNING: branch 'NeulandPoints' not found in "
                  << simFile << " -- neutron variables and Erel will be unfilled.\n";
    }

    // ------------------------------------------------------------------------
    // 4) Output histograms and tree
    // ------------------------------------------------------------------------
    TFile *fout = TFile::Open(outFile, "RECREATE");
    fout->cd();

    // ----- Reconstructed lab observables --------------------------------
    auto hPoQ = new TH1D("hPoQ", "Reconstructed P/Q;P/Q [GeV/c/e];counts", 600, -10.0, 20.0);
    auto hAoQ = new TH1D("hAoQ", "Reconstructed A/Q (truth-#beta closure);A/Q;counts", 400, 1.5, 4.5);
    auto hZ = new TH1D("hZ", "TofD charge (smeared, PID only);Z;counts", 400, 5.0, 10.0);
    auto hPlab = new TH1D("hPlab", "Reconstructed |p|_{frag} = PoQ #upoint Z_{frag};|p| [GeV/c];counts",
                          500, 0., 40.);

    // ----- Fragment momentum resolutions --------------------------------
    auto hPRes = new TH1D("hPRes",
                          "Total momentum resolution;(p_{rec}-p_{tru})/p_{tru};counts",
                          400, -0.05, 0.05);
    auto hPxRes = new TH1D("hPxRes",
                           "p_{x} resolution;(p_{x}^{rec}-p_{x}^{tru})/p_{tru};counts",
                           400, -0.05, 0.05);
    auto hPyRes = new TH1D("hPyRes",
                           "p_{y} resolution;(p_{y}^{rec}-p_{y}^{tru})/p_{tru};counts",
                           400, -0.05, 0.05);
    auto hPzRes = new TH1D("hPzRes",
                           "p_{z} resolution;(p_{z}^{rec}-p_{z}^{tru})/p_{tru};counts",
                           400, -0.5, 0.5);

    // ----- Neutron (first NeuLAND hit) -----------------------------------
    auto hBetaNeu = new TH1D("hBetaNeu", "Neutron #beta (first hit);#beta_{n};counts",
                             400, 0., 1.);
    auto hPNeu = new TH1D("hPNeu", "Neutron |p| (first hit);|p|_{n} [GeV/c];counts",
                          400, 0., 2.);

    // ----- Erel / response (same binning as DataAnalysis::fErel) ---------
    auto hErel = new TH1D("hErel", "Reconstructed E_{rel};E_{rel} [MeV];counts",
                          100, -5., 10.);
    auto hErelTrue = new TH1D("hErelTrue", "MC-truth E_{rel};E_{rel} [MeV];counts",
                              100, -5., 10.);
    auto hErelFragTrue = new TH1D("hErelFragTrue",
                                  "E_{rel} with true fragment momentum + reconstructed neutron;E_{rel} [MeV];counts",
                                  100, -5., 10.);
    auto hErelRes = new TH1D("hErelRes",
                             "E_{rel} response;E_{rel}^{rec}-E_{rel}^{tru} [MeV];counts",
                             200, -5., 5.);
    auto hErelRecVsTrue = new TH2D("hErelRecVsTrue",
                                   "Response matrix;E_{rel}^{tru} [MeV];E_{rel}^{rec} [MeV]",
                                   100, 0., 10., 100, -5., 10.);
    auto hErelFragTrueRes = new TH1D("hErelFragTrueRes",
                                     "E_{rel} response (true p_{frag} + reco neutron);"
                                     "E_{rel}^{fragTrue}-E_{rel}^{tru} [MeV];counts",
                                     200, -5., 5.);
    auto hErelFragTrueRecVsTrue = new TH2D("hErelFragTrueRecVsTrue",
                                           "Response matrix (true p_{frag} + reco neutron);"
                                           "E_{rel}^{tru} [MeV];E_{rel}^{fragTrue} [MeV]",
                                           100, 0., 10., 100, -5., 10.);
    auto hDeltaBeta = new TH1D("hDeltaBeta", "#Delta#beta = #beta_{frag}-#beta_{n};#Delta#beta;counts",
                               100, -0.3, 0.3);

    // Output ntuple
    auto tout = new TTree("ana", "G249 simulation reconstruction");
    Double_t br_PoQ, br_AoQ, br_Z, br_Plab, br_betaFrag;
    Double_t br_pxRec, br_pyRec, br_pzRec;
    Double_t br_pxTrue, br_pyTrue, br_pzTrue, br_pTrue;
    Double_t br_vx, br_vy, br_vz, br_tx, br_ty;
    Double_t br_betaNeu, br_pNeu, br_pxNeu, br_pyNeu, br_pzNeu;
    Double_t br_xNeuHit, br_yNeuHit, br_zNeuHit, br_tofNeuland;
    Double_t br_Erel, br_ErelTrue, br_ErelFragTrue, br_deltaBeta, br_opaLab;
    tout->Branch("PoQ", &br_PoQ, "PoQ/D");
    tout->Branch("AoQ", &br_AoQ, "AoQ/D");
    tout->Branch("Z", &br_Z, "Z/D");
    tout->Branch("Plab", &br_Plab, "Plab/D");
    tout->Branch("BetaFrag", &br_betaFrag, "BetaFrag/D");
    tout->Branch("PxRec", &br_pxRec, "PxRec/D");
    tout->Branch("PyRec", &br_pyRec, "PyRec/D");
    tout->Branch("PzRec", &br_pzRec, "PzRec/D");
    tout->Branch("PxTrue", &br_pxTrue, "PxTrue/D");
    tout->Branch("PyTrue", &br_pyTrue, "PyTrue/D");
    tout->Branch("PzTrue", &br_pzTrue, "PzTrue/D");
    tout->Branch("PTrue", &br_pTrue, "PTrue/D");
    tout->Branch("Vx", &br_vx, "Vx/D");
    tout->Branch("Vy", &br_vy, "Vy/D");
    tout->Branch("Vz", &br_vz, "Vz/D");
    tout->Branch("TX", &br_tx, "TX/D");
    tout->Branch("TY", &br_ty, "TY/D");
    tout->Branch("BetaNeu", &br_betaNeu, "BetaNeu/D");
    tout->Branch("PNeu", &br_pNeu, "PNeu/D");
    tout->Branch("PxNeu", &br_pxNeu, "PxNeu/D");
    tout->Branch("PyNeu", &br_pyNeu, "PyNeu/D");
    tout->Branch("PzNeu", &br_pzNeu, "PzNeu/D");
    tout->Branch("XNeuHit", &br_xNeuHit, "XNeuHit/D");
    tout->Branch("YNeuHit", &br_yNeuHit, "YNeuHit/D");
    tout->Branch("ZNeuHit", &br_zNeuHit, "ZNeuHit/D");
    tout->Branch("TofNeuland", &br_tofNeuland, "TofNeuland/D");
    tout->Branch("Erel", &br_Erel, "Erel/D");                         // GeV, as DataAnalysis
    tout->Branch("ErelTrue", &br_ErelTrue, "ErelTrue/D");             // GeV
    tout->Branch("ErelFragTrue", &br_ErelFragTrue, "ErelFragTrue/D"); // GeV: true p_frag + reco neutron
    tout->Branch("DeltaBeta", &br_deltaBeta, "DeltaBeta/D");
    tout->Branch("OpaLab", &br_opaLab, "OpaLab/D");

    // Bookkeeping counters
    Long64_t nProcessed = 0;
    Long64_t nNoPrimary = 0;
    Long64_t nFootMiss = 0;
    Long64_t nFiberMiss = 0;
    Long64_t nBadPoQ = 0;
    Long64_t nNoTofd = 0;
    Long64_t nReconstructed = 0;
    Long64_t nNoNeutron = 0;
    Long64_t nBadBetaNeu = 0;
    Long64_t nErelFilled = 0;
    Long64_t nPrimaryByPdg = 0;
    Long64_t nPrimaryByFallback = 0;
    constexpr Int_t kVerboseFirstN = 5;

    // ------------------------------------------------------------------------
    // 5) Event loop
    // ------------------------------------------------------------------------
    Long64_t nEntries = (maxEvents < 0)
                            ? tSim->GetEntries()
                            : std::min<Long64_t>(maxEvents, tSim->GetEntries());

    std::cout << "[anaDelta] Processing " << nEntries << " events\n";

    const Double_t M_frag = kFragMassAMU * kAMU_EF; // GeV, mass hypothesis (= data)
    const Double_t m_n = kMassNeutron;

    Long64_t iEv = -1;
    while (rdr.Next())
    {
        ++iEv;
        if (iEv >= nEntries)
            break;
        ++nProcessed;

        TClonesArray *mcArr = mcTrackVal.Get();
        const Int_t nMC = mcArr ? mcArr->GetEntriesFast() : 0;

        if (iEv < kVerboseFirstN)
        {
            std::cout << "  [diag] ev=" << iEv
                      << "  MCTrack=" << nMC
                      << "  TraPoint=" << trX.GetSize()
                      << "  Fi32=" << f32x.GetSize()
                      << "  Fi30=" << f30x.GetSize()
                      << "  Fi33=" << f33x.GetSize()
                      << "  Fi31=" << f31x.GetSize()
                      << "  TofD=" << tdEloss.GetSize()
                      << "  Neuland=" << (hasNeuland ? (int)nlX->GetSize() : -1)
                      << std::endl;
        }

        // -------- Find the fragment primary in MCTrack ----------------------
        const Int_t expectedPdg = 1000000000 + kFragZ * 10000 + kFragA * 10;
        Int_t iPrimary = -1;
        for (Int_t i = 0; i < nMC; ++i)
        {
            auto *m = dynamic_cast<R3BMCTrack *>(mcArr->At(i));
            if (m && m->GetMotherId() <= 0 && m->GetPdgCode() == expectedPdg)
            {
                iPrimary = i;
                ++nPrimaryByPdg;
                break;
            }
        }
        if (iPrimary < 0)
        {
            for (Int_t i = 0; i < nMC; ++i)
            {
                auto *m = dynamic_cast<R3BMCTrack *>(mcArr->At(i));
                if (m && m->GetPdgCode() == expectedPdg)
                {
                    iPrimary = i;
                    ++nPrimaryByPdg;
                    break;
                }
            }
        }
        if (iPrimary < 0 && nMC > 0)
        {
            iPrimary = 0;
            ++nPrimaryByFallback;
        }
        if (iPrimary < 0)
        {
            ++nNoPrimary;
            continue;
        }

        auto *mcPrim = dynamic_cast<R3BMCTrack *>(mcArr->At(iPrimary));
        if (!mcPrim)
        {
            ++nNoPrimary;
            continue;
        }
        const Double_t px_true = mcPrim->GetPx();
        const Double_t py_true = mcPrim->GetPy();
        const Double_t pz_true = mcPrim->GetPz();
        const Double_t p_true = std::sqrt(px_true * px_true + py_true * py_true + pz_true * pz_true);

        // -------- Find the neutron primary (for truth Erel) ------------------
        Int_t iNeuTrue = -1;
        for (Int_t i = 0; i < nMC; ++i)
        {
            if (i == iPrimary)
                continue;
            auto *m = dynamic_cast<R3BMCTrack *>(mcArr->At(i));
            if (m && m->GetPdgCode() == 2112 && m->GetMotherId() <= 0)
            {
                iNeuTrue = i;
                break;
            }
        }
        if (iNeuTrue < 0)
        {
            for (Int_t i = 0; i < nMC; ++i)
            {
                if (i == iPrimary)
                    continue;
                auto *m = dynamic_cast<R3BMCTrack *>(mcArr->At(i));
                if (m && m->GetPdgCode() == 2112)
                {
                    iNeuTrue = i;
                    break;
                }
            }
        }

        // -------- FOOT-equivalent track --------------------------------------
        struct FootHit
        {
            bool ok = false;
            Double_t pos = 0., z = 0., eloss = -1.;
        };
        FootHit h5, h6, h7, h8; // 5,7 -> Y ; 6,8 -> X
        const Int_t nTra = static_cast<Int_t>(trX.GetSize());
        for (Int_t i = 0; i < nTra; ++i)
        {
            const Int_t did = trDet[i];
            if (did < 5 || did > 8)
                continue;
            FootHit *h = nullptr;
            Double_t measured = 0.;
            switch (did)
            {
            case 5:
                h = &h5;
                measured = trY[i];
                break;
            case 6:
                h = &h6;
                measured = trX[i];
                break;
            case 7:
                h = &h7;
                measured = trY[i];
                break;
            case 8:
                h = &h8;
                measured = trX[i];
                break;
            }
            if (trE[i] > h->eloss)
            {
                h->ok = true;
                h->pos = measured;
                h->z = trZ[i];
                h->eloss = trE[i];
            }
        }
        if (!(h5.ok && h6.ok && h7.ok && h8.ok))
        {
            ++nFootMiss;
            continue;
        }

        const Double_t Y5 = h5.pos + gRandom->Gaus(0., kSigmaFootPos);
        const Double_t X6 = h6.pos + gRandom->Gaus(0., kSigmaFootPos);
        const Double_t Y7 = h7.pos + gRandom->Gaus(0., kSigmaFootPos);
        const Double_t X8 = h8.pos + gRandom->Gaus(0., kSigmaFootPos);
        const Double_t Z5 = h5.z, Z6 = h6.z, Z7 = h7.z, Z8 = h8.z;

        const Double_t tx_smear = (X8 - X6) / (Z8 - Z6);
        const Double_t ty_smear = (Y7 - Y5) / (Z7 - Z5);
        const Double_t footStartX = X6;
        const Double_t footStartY = Y5;
        const Double_t footStartZ = Z5;

        // -------- Fiber hits --------------------------------------------------
        FiberLabHit f32 = PickMaxElossHit(f32x, f32y, f32z, f32e);
        FiberLabHit f30 = PickMaxElossHit(f30x, f30y, f30z, f30e);
        FiberLabHit f31 = PickMaxElossHit(f31x, f31y, f31z, f31e);
        FiberLabHit f33 = PickMaxElossHit(f33x, f33y, f33z, f33e);

        if (!f32.valid || !f30.valid || (!f31.valid && !f33.valid))
        {
            ++nFiberMiss;
            continue;
        }

        { // Fi32 (X-fiber)
            TVector3 d = ToDet(f32.lab, gPosFi32);
            d.SetX(d.X() + gRandom->Gaus(0., kSigmaFiberPos));
            d.SetY(0.);
            d.SetZ(0.);
            f32.lab = ToLab(d, gPosFi32);
        }
        { // Fi30 (Y-fiber)
            TVector3 d = ToDet(f30.lab, gPosFi30);
            d.SetY(d.Y() + gRandom->Gaus(0., kSigmaFiberPos));
            d.SetX(0.);
            d.SetZ(0.);
            f30.lab = ToLab(d, gPosFi30);
        }
        if (f31.valid)
        {
            TVector3 d = ToDet(f31.lab, gPosFi31);
            d.SetX(d.X() + gRandom->Gaus(0., kSigmaFiberPos));
            d.SetY(0.);
            d.SetZ(0.);
            f31.lab = ToLab(d, gPosFi31);
        }
        if (f33.valid)
        {
            TVector3 d = ToDet(f33.lab, gPosFi33);
            d.SetX(d.X() + gRandom->Gaus(0., kSigmaFiberPos));
            d.SetY(0.);
            d.SetZ(0.);
            f33.lab = ToLab(d, gPosFi33);
        }

        bool useFi33 = (f33.valid && f31.valid) ? (f33.eloss >= f31.eloss) : f33.valid;
        const TVector3 fLast_lab = useFi33 ? f33.lab : f31.lab;

        TVector3 fFi30_lab = f30.lab;
        FixupFib30(f32.lab, fFi30_lab, fLast_lab);

        // -------- MDF input (exactly as R3BTrackingG249::Exec) ----------------
        const Double_t startX = footStartX;
        const Double_t startY = footStartY;
        const Double_t startZ = footStartZ;
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

        bool inputBad = false;
        for (int kk = 0; kk < 9; ++kk)
            if (!std::isfinite(mdf_data[kk]))
            {
                inputBad = true;
                break;
            }
        if (inputBad)
        {
            ++nBadPoQ;
            continue;
        }

        Double_t PoQ_raw = mdfPoQ->MDF(mdf_data);
        Double_t PoQ = PoQ_raw * (kGladCurrent / kGladRefMDF);

        if (iEv < kVerboseFirstN)
        {
            std::cout << "  [diag] ev=" << iEv
                      << " PoQ_raw=" << PoQ_raw << " PoQ=" << PoQ << std::endl;
        }
        if (!std::isfinite(PoQ))
        {
            ++nBadPoQ;
            continue;
        }

        // -------- ToFD charge: PID observable ONLY ---------------------------
        // In the experimental chain (eventFilter.C) the momentum is built
        // with the FIXED nominal charge Z_FRAG = 8; the measured Z is only a
        // selection variable. We reproduce that: Z_smear feeds hZ / branch Z
        // and nothing else.
        if (tdEloss.GetSize() == 0)
            ++nNoTofd; // flagged, never vetoed
        const Double_t Z_smear =
            static_cast<Double_t>(kFragZ) + gRandom->Gaus(0., kSigmaTofdCharge);

        // -------- Fragment momentum, beta (mirror eventFilter/dataAnalysis) --
        TVector3 pdir(startTX, startTY, 1.);
        pdir = pdir.Unit();
        const Double_t PoQ_abs = std::fabs(PoQ);
        const Double_t p_frag_raw = PoQ_abs * static_cast<Double_t>(kFragZ); // Z fixed!
        const Double_t beta_frag_raw = p_frag_raw / std::sqrt(p_frag_raw * p_frag_raw + M_frag * M_frag);
        // Beta offset correction, as DataAnalysis::getData() (beta_frag += fBetaMatchValue):
        // |p_frag| is then rebuilt from the corrected beta, only the track direction
        // is kept from the raw reconstruction, so everything downstream (gamma_frag,
        // Erel, deltaBeta, momentum residuals) is recalculated from the corrected beta.
        const Double_t beta_frag = beta_frag_raw + betaFragOffset;
        const Double_t gamma_frag = 1. / std::sqrt(1. - beta_frag * beta_frag);
        const Double_t p_frag = gamma_frag * M_frag * beta_frag;
        TVector3 p3_frag = pdir * p_frag;

        // -------- A/Q closure diagnostic (uses MC-truth beta, as before) -----
        const Double_t M_true = static_cast<Double_t>(kFragA) * kAMU;
        const Double_t E_true = std::sqrt(p_true * p_true + M_true * M_true);
        const Double_t beta_true = p_true / E_true;
        const Double_t gamma_true = 1. / std::sqrt(1. - beta_true * beta_true);
        const Double_t AoQ = PoQ_abs / (beta_true * gamma_true * kAMU);

        // -------- Neutron (first hit, fully smeared) --------------------------
        Double_t betaNeu = -999., pNeu = -999.;
        Double_t pxNeu = -999., pyNeu = -999., pzNeu = -999.;
        Double_t xNeuHit = -999., yNeuHit = -999., zNeuHit = -999., tofNeuland = -999.;
        if (hasNeuland)
        {
            NeutronHit nHit = FindFirstNeutron(*nlX, *nlY, *nlZ, *nlT, *nlE, *nlDet);
            if (!nHit.valid)
            {
                ++nNoNeutron;
            }
            else
            {
                const Double_t L = nHit.pos.Mag();
                const Double_t tns = nHit.tns;
                Double_t b = (L > 0. && tns > 0.) ? L / (kCcmPerNs * tns) : -1.;
                if (b <= 0. || b >= 1.)
                {
                    ++nBadBetaNeu;
                }
                else
                {
                    betaNeu = b;
                    const Double_t gammaNeu = 1. / std::sqrt(1. - betaNeu * betaNeu);
                    pNeu = gammaNeu * m_n * betaNeu;
                    TVector3 nDir = nHit.pos.Unit();
                    TVector3 pvecNeu = nDir * pNeu;
                    pxNeu = pvecNeu.X();
                    pyNeu = pvecNeu.Y();
                    pzNeu = pvecNeu.Z();
                    xNeuHit = nHit.pos.X();
                    yNeuHit = nHit.pos.Y();
                    zNeuHit = nHit.pos.Z();
                    tofNeuland = tns;

                    hBetaNeu->Fill(betaNeu);
                    hPNeu->Fill(pNeu);
                }
            }
        }

        // -------- Erel: EXACTLY as DataAnalysis::getData() --------------------
        // Slope vectors, centered by the same px/py offsets as
        // DataAnalysis::getData() (dx_neu_off/dy_neu_off/fx_frag_off/fy_frag_off):
        //   d_neu  = (px_n/pz_n + neuPxOffset,  py_n/pz_n + neuPyOffset,  1)
        //   f_frag = (px_f/pz_f + fragPxOffset, py_f/pz_f + fragPyOffset, 1)
        //   cos    = (dx fx + dy fy + 1) / (|d| |f|)
        //   Erel   = sqrt(mf^2+mn^2+2 g_n g_f mf mn (1 - b_n b_f cos)) - mf - mn
        Double_t Erel = -999., ErelTrue = -999., ErelFragTrue = -999.,
                 deltaBeta = -999., opaLab = -999.;
        if (betaNeu > 0.)
        {
            const Double_t dx_neu = pxNeu / pzNeu + neuPxOffset;
            const Double_t dy_neu = pyNeu / pzNeu + neuPyOffset;
            const Double_t fx_frag = p3_frag.X() / p3_frag.Z() + fragPxOffset;
            const Double_t fy_frag = p3_frag.Y() / p3_frag.Z() + fragPyOffset;

            const Double_t cos_ang =
                (dx_neu * fx_frag + dy_neu * fy_frag + 1.0) /
                (std::sqrt(dx_neu * dx_neu + dy_neu * dy_neu + 1.0) *
                 std::sqrt(fx_frag * fx_frag + fy_frag * fy_frag + 1.0));

            const Double_t gamma_neu = 1. / std::sqrt(1. - betaNeu * betaNeu);

            Erel = std::sqrt(M_frag * M_frag + m_n * m_n +
                             2.0 * gamma_neu * gamma_frag * M_frag * m_n *
                                 (1.0 - betaNeu * beta_frag * cos_ang)) -
                   M_frag - m_n;

            deltaBeta = beta_frag - betaNeu;
            opaLab = TVector3(pxNeu, pyNeu, pzNeu).Angle(p3_frag);

            hErel->Fill(Erel * 1000.); // MeV, as DataAnalysis::fErel
            hDeltaBeta->Fill(deltaBeta);
            ++nErelFilled;

            // ---- Same Erel, but with the TRUE fragment momentum (MC truth)
            //      in place of the MDF-reconstructed one; the neutron stays
            //      the reconstructed one. Isolates the fragment-tracking
            //      contribution to the Erel resolution. Mass hypothesis
            //      (M_frag) kept the same as the reconstructed Erel so only
            //      the momentum vector changes.
            const Double_t fx_fragTrue = px_true / pz_true;
            const Double_t fy_fragTrue = py_true / pz_true;

            const Double_t cos_angFragTrue =
                (dx_neu * fx_fragTrue + dy_neu * fy_fragTrue + 1.0) /
                (std::sqrt(dx_neu * dx_neu + dy_neu * dy_neu + 1.0) *
                 std::sqrt(fx_fragTrue * fx_fragTrue + fy_fragTrue * fy_fragTrue + 1.0));

            const Double_t beta_fragTrue = p_true / std::sqrt(p_true * p_true + M_frag * M_frag);
            const Double_t gamma_fragTrue = 1. / std::sqrt(1. - beta_fragTrue * beta_fragTrue);

            ErelFragTrue = std::sqrt(M_frag * M_frag + m_n * m_n +
                                     2.0 * gamma_neu * gamma_fragTrue * M_frag * m_n *
                                         (1.0 - betaNeu * beta_fragTrue * cos_angFragTrue)) -
                           M_frag - m_n;

            hErelFragTrue->Fill(ErelFragTrue * 1000.); // MeV
        }

        // -------- Truth Erel from the two primaries ----------------------------
        // Same masses as the reconstruction (mass hypothesis), so the residual
        // Erel_rec - Erel_true isolates the *detector* response.
        if (iNeuTrue >= 0)
        {
            auto *mcN = dynamic_cast<R3BMCTrack *>(mcArr->At(iNeuTrue));
            if (mcN)
            {
                const Double_t pxn = mcN->GetPx();
                const Double_t pyn = mcN->GetPy();
                const Double_t pzn = mcN->GetPz();
                const Double_t pn2 = pxn * pxn + pyn * pyn + pzn * pzn;
                const Double_t E_f = std::sqrt(p_true * p_true + M_frag * M_frag);
                const Double_t E_n = std::sqrt(pn2 + m_n * m_n);
                const Double_t dot = px_true * pxn + py_true * pyn + pz_true * pzn;
                const Double_t Minv2 =
                    M_frag * M_frag + m_n * m_n + 2. * (E_f * E_n - dot);
                ErelTrue = std::sqrt(std::max(0., Minv2)) - M_frag - m_n;
                hErelTrue->Fill(ErelTrue * 1000.);
                if (Erel > -998.)
                {
                    hErelRes->Fill((Erel - ErelTrue) * 1000.);
                    hErelRecVsTrue->Fill(ErelTrue * 1000., Erel * 1000.);
                }
                if (ErelFragTrue > -998.)
                {
                    hErelFragTrueRes->Fill((ErelFragTrue - ErelTrue) * 1000.);
                    hErelFragTrueRecVsTrue->Fill(ErelTrue * 1000., ErelFragTrue * 1000.);
                }
            }
        }

        // -------- Fill histograms / tree ---------------------------------------
        hPoQ->Fill(PoQ);
        hAoQ->Fill(AoQ);
        hZ->Fill(Z_smear);
        hPlab->Fill(p_frag);
        hPRes->Fill((p_frag - p_true) / p_true);
        hPxRes->Fill((p3_frag.X() - px_true) / p_true);
        hPyRes->Fill((p3_frag.Y() - py_true) / p_true);
        hPzRes->Fill((p3_frag.Z() - pz_true) / p_true);

        br_PoQ = PoQ;
        br_AoQ = AoQ;
        br_Z = Z_smear;
        br_Plab = p_frag;
        br_betaFrag = beta_frag;
        br_pxRec = p3_frag.X();
        br_pyRec = p3_frag.Y();
        br_pzRec = p3_frag.Z();
        br_pxTrue = px_true;
        br_pyTrue = py_true;
        br_pzTrue = pz_true;
        br_pTrue = p_true;
        br_vx = startX;
        br_vy = startY;
        br_vz = startZ;
        br_tx = startTX;
        br_ty = startTY;
        br_betaNeu = betaNeu;
        br_pNeu = pNeu;
        br_pxNeu = pxNeu;
        br_pyNeu = pyNeu;
        br_pzNeu = pzNeu;
        br_xNeuHit = xNeuHit;
        br_yNeuHit = yNeuHit;
        br_zNeuHit = zNeuHit;
        br_tofNeuland = tofNeuland;
        br_Erel = Erel;
        br_ErelTrue = ErelTrue;
        br_ErelFragTrue = ErelFragTrue;
        br_deltaBeta = deltaBeta;
        br_opaLab = opaLab;
        tout->Fill();

        ++nReconstructed;

        if (nReconstructed <= 3)
        {
            std::cout << "  event " << iEv
                      << " : PoQ=" << PoQ
                      << " |p|_frag=" << p_frag
                      << " |p|_true=" << p_true
                      << " beta_frag=" << beta_frag;
            if (betaNeu > 0.)
                std::cout << " | n: beta=" << betaNeu << " |p|=" << pNeu
                          << " | Erel=" << Erel * 1000. << " MeV"
                          << " (true " << ErelTrue * 1000. << " MeV)";
            std::cout << std::endl;
        }
    } // event loop

    // ------------------------------------------------------------------------
    // 6) Write and close
    // ------------------------------------------------------------------------
    fout->cd();
    tout->Write();
    hPoQ->Write();
    hAoQ->Write();
    hZ->Write();
    hPlab->Write();
    hPRes->Write();
    hPxRes->Write();
    hPyRes->Write();
    hPzRes->Write();
    hBetaNeu->Write();
    hPNeu->Write();
    hErel->Write();
    hErelTrue->Write();
    hErelFragTrue->Write();
    hErelRes->Write();
    hErelRecVsTrue->Write();
    hErelFragTrueRes->Write();
    hErelFragTrueRecVsTrue->Write();
    hDeltaBeta->Write();
    fout->Close();

    delete tSim;

    timer.Stop();
    std::cout << "\n[anaDelta] Done.\n"
              << "  events processed      : " << nProcessed << "\n"
              << "  primary by PDG match  : " << nPrimaryByPdg << "\n"
              << "  primary by fallback   : " << nPrimaryByFallback << "\n"
              << "  -- gates that vetoed events --\n"
              << "  no primary fragment   : " << nNoPrimary << "\n"
              << "  FOOT multiplicity     : " << nFootMiss << "\n"
              << "  fiber multiplicity    : " << nFiberMiss << "\n"
              << "  MDF returned bad PoQ  : " << nBadPoQ << "\n"
              << "  -- non-vetoing flags --\n"
              << "  no TofD hit           : " << nNoTofd << "\n"
              << "  no neutron first hit  : " << nNoNeutron
              << (hasNeuland ? "" : " (no NeulandPoints branch)") << "\n"
              << "  bad neutron beta      : " << nBadBetaNeu << "\n"
              << "  -- result --\n"
              << "  fragment reconstructed: " << nReconstructed << "\n"
              << "  Erel filled           : " << nErelFilled << "\n"
              << "  output                : " << outFile << "\n"
              << "  real time             : " << timer.RealTime() << " s\n"
              << "  CPU time              : " << timer.CpuTime() << " s\n";

    delete mdfPoQ;
}

// ----------------------------------------------------------------------------
// ScanResolutions()
//
// Run the full reconstruction for 4 values of each detector's resolution
// (FOOT position, Fiber position, NeuLAND along-bar coordinate), keeping the
// others at their nominal values, and overlay the chosen histogram
// (default: hErel -- the observable of this study).
//
//   ScanResolutions()
//   ScanResolutions("delta.simu.root","PoQ_9vars_150terms.txt",5000,"hPyRes")
// ----------------------------------------------------------------------------
void ScanResolutions(
    const char *simFile = "sim_results/delta0.simu.root",
    const char *mdfFile = "PoQ_9vars_150terms.txt",
    Long64_t nEvents = 5000,
    const char *histName = "hErel")
{
    // Nominal values -- NOW really equal to the defaults above
    const Double_t nomFoot = 75.e-4;   // cm
    const Double_t nomFiber = 150.e-4; // cm
    const Double_t nomNeuAlong = 1.5;  // cm

    const Int_t nV = 4;
    Double_t footVals[nV] = {75.e-4, 100.e-4, 150.e-4, 200.e-4};
    Double_t fiberVals[nV] = {150.e-4, 200.e-4, 250.e-4, 1000.e-4};
    Double_t neuVals[nV] = {0.75, 1.5, 3.0, 5.0};

    const Int_t colors[nV] = {kBlue + 1, kGreen + 2, kRed, kViolet + 1};

    UInt_t seed = 42;
    auto doRun = [&](const char *tag) -> TH1D *
    {
        TString tmp = Form("_scanres_%s.root", tag);
        anaDelta(simFile, mdfFile, tmp.Data(), nEvents, seed);
        seed += 100;
        TFile *f = TFile::Open(tmp.Data(), "READ");
        if (!f || f->IsZombie())
        {
            std::cerr << "[ScanResolutions] ERROR: could not open " << tmp << "\n";
            return nullptr;
        }
        auto *h = dynamic_cast<TH1D *>(f->Get(histName));
        TH1D *hc = h ? static_cast<TH1D *>(h->Clone(Form("%s_%s", histName, tag))) : nullptr;
        if (hc)
            hc->SetDirectory(nullptr);
        f->Close();
        gSystem->Unlink(tmp.Data());
        return hc;
    };

    // ---- FOOT position scan ------------------------------------------------
    TH1D *hFoot[nV] = {};
    kSigmaFiberPos = nomFiber;
    kSigmaNeulandAlong = nomNeuAlong;
    for (Int_t i = 0; i < nV; ++i)
    {
        kSigmaFootPos = footVals[i];
        hFoot[i] = doRun(Form("foot%d", i));
    }
    kSigmaFootPos = nomFoot;

    // ---- Fiber position scan -------------------------------------------------
    TH1D *hFiber[nV] = {};
    for (Int_t i = 0; i < nV; ++i)
    {
        kSigmaFiberPos = fiberVals[i];
        hFiber[i] = doRun(Form("fiber%d", i));
    }
    kSigmaFiberPos = nomFiber;

    // ---- NeuLAND along-bar scan -----------------------------------------------
    TH1D *hNeu[nV] = {};
    for (Int_t i = 0; i < nV; ++i)
    {
        kSigmaNeulandAlong = neuVals[i];
        hNeu[i] = doRun(Form("neu%d", i));
    }
    kSigmaNeulandAlong = nomNeuAlong;

    auto drawScan = [&](TH1D **hArr, const char *cname,
                        Double_t *vals, const char *unit, Double_t scale)
    {
        TCanvas *c = new TCanvas(cname, cname, 900, 650);
        c->SetLeftMargin(0.13);
        TLegend *leg = new TLegend(0.58, 0.58, 0.93, 0.92);
        leg->SetBorderSize(0);
        leg->SetFillStyle(0);

        Double_t yMax = 0.;
        for (Int_t i = 0; i < nV; ++i)
            if (hArr[i])
                yMax = TMath::Max(yMax, hArr[i]->GetMaximum());

        Bool_t first = kTRUE;
        for (Int_t i = 0; i < nV; ++i)
        {
            if (!hArr[i])
                continue;
            hArr[i]->SetLineColor(colors[i]);
            hArr[i]->SetLineWidth(2);
            if (first)
            {
                hArr[i]->GetYaxis()->SetRangeUser(0., yMax * 1.15);
                hArr[i]->SetTitle(Form("%s scan (%s)", histName, cname));
                hArr[i]->Draw("HIST");
                first = kFALSE;
            }
            else
                hArr[i]->Draw("HIST SAME");
            leg->AddEntry(hArr[i],
                          Form("#sigma = %.2f %s", vals[i] * scale, unit), "l");
        }
        leg->Draw();
        c->Update();
    };

    drawScan(hFoot, "cScanFoot", footVals, "#mum", 1.e4);
    drawScan(hFiber, "cScanFiber", fiberVals, "#mum", 1.e4);
    drawScan(hNeu, "cScanNeuland", neuVals, "cm", 1.0);
}