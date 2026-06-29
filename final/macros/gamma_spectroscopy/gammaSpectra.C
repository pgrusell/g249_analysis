#include "gammaSpectra.h"

// ─────────────────────────────────────────────────────────────────────────────
//  Configuration
// ─────────────────────────────────────────────────────────────────────────────

// OPA (proton-proton opening angle) cut in radians – same as fitMomdis.C
static constexpr double OPA_MIN = 1.3;
static constexpr double OPA_MAX = 1.6;

// k-sigma for PID ellipse (same convention as CrossSections.cpp)
static constexpr double PID_K_DEFAULT = 2.5;

// 22O outgoing PID window (coarse, used to seed the 2D Gaussian fit)
static constexpr double AOQ_22O_MIN = 2.67;
static constexpr double AOQ_22O_MAX = 2.76;
static constexpr double Z_22O_MIN   = 7.6;
static constexpr double Z_22O_MAX   = 8.5;

// 25F incoming PID window (used only on the unreacted file if provided)
static constexpr double AOQ_25F_MIN = 2.71;
static constexpr double AOQ_25F_MAX = 2.77;
static constexpr double Z_25F_MIN   = 8.3;
static constexpr double Z_25F_MAX   = 9.5;

// Gamma energy histogram range — 50 keV/bin
static constexpr int    NBINS_E    = 200;   // 10 MeV / 0.05 MeV per bin
static constexpr double E_MIN_MEV  = 0.0;
static constexpr double E_MAX_MEV  = 10.0;

// Gamma theta histogram range (degrees) — 1 deg/bin, upper limit 100 deg
static constexpr int    NBINS_TH   = 100;
static constexpr double TH_MIN_DEG = 0.0;
static constexpr double TH_MAX_DEG = 100.0;

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
                  double  opaMax        = OPA_MAX)
{
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

    std::cout << "\n=== gamma spectroscopy: 23F(p,2p)22O ===\n";
    std::cout << "  Reaction file : " << reactionFile   << "\n";
    std::cout << "  Output file   : " << outFile        << "\n";
    std::cout << "  Unreacted file: " << (unreactedFile.IsNull() ? "(none)" : unreactedFile.Data()) << "\n";
    std::cout << "  Offset file   : " << offsetFile     << "\n";
    std::cout << "  kPID          : " << kPID    << "\n";
    std::cout << "  OPA cut [rad] : " << opaMin  << " - " << opaMax << "\n\n";

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

    // ── Pass 1: fit outgoing 22O PID ─────────────────────────────────────────
    fout_root->cd();
    Fit2DParams pid22O = fitPID(treac, "h2_pid_22O",
                                 AOQ_22O_MIN, AOQ_22O_MAX,
                                 Z_22O_MIN,   Z_22O_MAX, 1);

    // Fall back to rectangular window if fit failed
    if (!pid22O.valid) {
        pid22O.muAoQ    = 0.5 * (AOQ_22O_MIN + AOQ_22O_MAX);
        pid22O.sigmaAoQ = (AOQ_22O_MAX - AOQ_22O_MIN) / (2.0 * kPID);
        pid22O.muZ      = 0.5 * (Z_22O_MIN + Z_22O_MAX);
        pid22O.sigmaZ   = (Z_22O_MAX - Z_22O_MIN) / (2.0 * kPID);
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
        500, AOQ_22O_MIN-0.05, AOQ_22O_MAX+0.05,
        500, Z_22O_MIN-0.2,   Z_22O_MAX+0.2);

    auto* h1_opa_after_pid = new TH1F("h1_opa_after_pid",
        "OPA after 22O PID (no OPA cut);#theta_{OPA} [rad];Counts",
        180, 0.0, TMath::Pi());

    auto* h1_beta_frag = new TH1F("h1_beta_frag",
        "Fragment #beta (22O, after all cuts);#beta_{frag};Counts",
        200, 0.78, 0.84);

    // ---- Gamma cluster multiplicity ----
    auto* h1_gamma_mult = new TH1F("h1_gamma_mult",
        "Gamma cluster multiplicity (after all cuts);M_{#gamma};Counts",
        16, -0.5, 15.5);

    // ---- Raw (uncorrected) gamma energy, all multiplicities ----
    auto* h1_gamma_E_raw_all = new TH1F("h1_gamma_E_raw_all",
        "Raw #gamma energy, all mult;E_{#gamma} [MeV];Counts",
        NBINS_E, E_MIN_MEV, E_MAX_MEV);

    // ---- Doppler-corrected gamma energy by exact multiplicity ----
    const int NMULT = 5;
    TH1F* h1_gamma_E_corr_mult[NMULT];
    for (int m = 1; m <= NMULT; ++m) {
        h1_gamma_E_corr_mult[m-1] = new TH1F(
            Form("h1_gamma_E_corr_mult%d", m),
            Form("Doppler-corrected #gamma, M_{#gamma}==%d;E_{#gamma} [MeV];Counts", m),
            NBINS_E, E_MIN_MEV, E_MAX_MEV);
    }

    // ---- Doppler-corrected gamma energy by cumulative multiplicity (M<=m) ----
    TH1F* h1_gamma_E_corr_mleq[NMULT];
    for (int m = 1; m <= NMULT; ++m) {
        h1_gamma_E_corr_mleq[m-1] = new TH1F(
            Form("h1_gamma_E_corr_mleq%d", m),
            Form("Doppler-corrected #gamma, M_{#gamma}<=%d;E_{#gamma} [MeV];Counts", m),
            NBINS_E, E_MIN_MEV, E_MAX_MEV);
    }

    // ---- All multiplicities combined (convenience) ----
    auto* h1_gamma_E_corr_all = new TH1F("h1_gamma_E_corr_all",
        "Doppler-corrected #gamma, all mult;E_{#gamma} [MeV];Counts",
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

    const Long64_t nEntries = treac->GetEntries();
    Long64_t nPassPID = 0, nPassOPA = 0;

    for (Long64_t i = 0; i < nEntries; ++i) {
        treac->GetEntry(i);

        if (i % 100000 == 0)
            std::cout << "\r  processing event " << i << " / " << nEntries << std::flush;

        // ── Outgoing 22O PID ellipse cut ─────────────────────────────────
        h2_pid_all->Fill(AoQ_frag, Z_frag_est);

        if (!insideEllipse(AoQ_frag, Z_frag_est,
                           pid22O.muAoQ, pid22O.sigmaAoQ,
                           pid22O.muZ,   pid22O.sigmaZ,  kPID))
            continue;
        ++nPassPID;

        h1_opa_after_pid->Fill(califa_opa);

        // ── OPA cut ──────────────────────────────────────────────────────
        if (califa_opa < opaMin || califa_opa > opaMax) continue;
        ++nPassOPA;

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

        const int mult = califa_gamma_mult;
        h1_gamma_mult->Fill(static_cast<double>(mult));

        // ── Per-gamma cluster loop ────────────────────────────────────────
        for (int ig = 0; ig < mult; ++ig) {
            // energies stored in keV → convert to MeV
            const double E_raw_mev  = (*califa_gamma_E)[ig]     / 1000.0;
            const double theta_rad  = (*califa_gamma_theta)[ig];
            const double phi_rad    = (*califa_gamma_phi)[ig];
            const double theta_deg  = theta_rad * TMath::RadToDeg();

            // Doppler-correct using offset-corrected fragment direction and beta
            const double E_corr_mev = DopplerCorrect(E_raw_mev,
                                                      theta_rad, phi_rad,
                                                      px_corr, py_corr, pz_safe,
                                                      beta_corr);

            // Raw histograms
            h1_gamma_E_raw_all->Fill(E_raw_mev);
            h2_gamma_E_vs_theta_all->Fill(theta_deg, E_raw_mev);

            // Corrected histograms (all mult)
            h1_gamma_E_corr_all->Fill(E_corr_mev);
            h2_gamma_Ecorr_vs_theta_all->Fill(theta_deg, E_corr_mev);

            // Exact multiplicity histograms (mult == m, m=1..NMULT)
            if (mult >= 1 && mult <= NMULT)
                h1_gamma_E_corr_mult[mult-1]->Fill(E_corr_mev);

            // Cumulative multiplicity histograms (mult <= m, m=1..NMULT)
            for (int m = 1; m <= NMULT; ++m)
                if (mult <= m) h1_gamma_E_corr_mleq[m-1]->Fill(E_corr_mev);

            // 2D E vs theta for M<=2 and M<=3
            if (mult <= 2) h2_gamma_Ecorr_vs_theta_mleq2->Fill(theta_deg, E_corr_mev);
            if (mult <= 3) h2_gamma_Ecorr_vs_theta_mleq3->Fill(theta_deg, E_corr_mev);
        }
    }

    std::cout << "\n\n[SUMMARY]\n";
    std::cout << "  Total events      : " << nEntries  << "\n";
    std::cout << "  After PID cut     : " << nPassPID  << "\n";
    std::cout << "  After OPA cut     : " << nPassOPA  << "\n";

    // ── Write all histograms ─────────────────────────────────────────────────
    fout_root->cd();

    h2_pid_all->Write();
    h1_opa_after_pid->Write();
    h1_beta_frag->Write();
    h1_gamma_mult->Write();
    h1_gamma_E_raw_all->Write();
    h1_gamma_E_corr_all->Write();
    for (int m = 1; m <= NMULT; ++m) {
        h1_gamma_E_corr_mult[m-1]->Write();
        h1_gamma_E_corr_mleq[m-1]->Write();
    }
    h2_gamma_E_vs_theta_all->Write();
    h2_gamma_Ecorr_vs_theta_all->Write();
    h2_gamma_Ecorr_vs_theta_mleq2->Write();
    h2_gamma_Ecorr_vs_theta_mleq3->Write();

    fout_root->Close();
    freac->Close();

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
                     " [outFile.root] [unreactedFile.root]"
                     " [kPID] [opaMin] [opaMax]\n";
        return 1;
    }

    TString reac  = argv[1];
    TString out   = (argc > 2) ? argv[2] : "gamma_spectra_22O.root";  // bare name → resolved inside gammaSpectra()
    TString unr   = (argc > 3) ? argv[3] : "";
    TString offs  = (argc > 4) ? argv[4] : "22O.txt";
    double  kPID  = (argc > 5) ? std::stod(argv[5]) : PID_K_DEFAULT;
    double  opaLo = (argc > 6) ? std::stod(argv[6]) : OPA_MIN;
    double  opaHi = (argc > 7) ? std::stod(argv[7]) : OPA_MAX;

    gammaSpectra(reac, out, unr, offs, kPID, opaLo, opaHi);
    return 0;
}
#endif
