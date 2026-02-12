/*

:setting: .txt file containing the paths of the files of the setting (e.g. 25F.txt)
:reaction: name of the reaction to study (e.g. 25F23O)

*/

/////// Some constants ////////
static constexpr double AMU_GeV = 0.93149410242;
static constexpr double Z_FRAG = 8.;
static constexpr double SEL_Z_MIN = 7.6;
static constexpr double SEL_Z_MAX = 8.5;
static constexpr double MN_GeV = 0.939565420;
static constexpr double C_CM_PER_NS = 29.9792458;

/////// Helpers //////

// Trim helper
static inline void trim_inplace(std::string &s)
{
    auto notSpace = [](unsigned char c)
    { return !std::isspace(c); };
    s.erase(s.begin(), std::find_if(s.begin(), s.end(), notSpace));
    s.erase(std::find_if(s.rbegin(), s.rend(), notSpace).base(), s.end());
}

// Remove surrounding quotes "..."
static inline void strip_quotes_inplace(std::string &s)
{
    if (s.size() >= 2 && ((s.front() == '"' && s.back() == '"') ||
                          (s.front() == '\'' && s.back() == '\'')))
    {
        s = s.substr(1, s.size() - 2);
    }
}

// Reads file list from text file: one path per line, optional quotes.
// Supports comments starting with # (whole line) and inline " # comment" if you want.
static std::vector<std::string> read_file_list(const std::string &txtPath)
{
    std::ifstream in(txtPath);
    if (!in.is_open())
        throw std::runtime_error("Cannot open file list: " + txtPath);

    std::vector<std::string> files;
    files.reserve(256);

    std::string line;
    while (std::getline(in, line))
    {
        // Remove inline comments starting with #
        // (If you DON'T want inline comments, delete this block.)
        auto hashPos = line.find('#');
        if (hashPos != std::string::npos)
            line = line.substr(0, hashPos);

        trim_inplace(line);
        if (line.empty())
            continue;

        strip_quotes_inplace(line);
        trim_inplace(line);

        if (!line.empty())
            files.push_back(line);
    }

    if (files.empty())
        throw std::runtime_error("File list is empty: " + txtPath);

    return files;
}

bool isGoodIncoming(double AoQ, double Z, double minAoQ, double maxAoQ)
{
    return (Z >= 8.2 && Z <= 10.0 &&
            AoQ >= minAoQ && AoQ <= maxAoQ);
}

bool isGoodFootTrack(TClonesArray *ft)
{
    static constexpr double Zmin[8] = {1, 1, 1, 8.5, 7.5, 1, 1, 1};
    static constexpr double Zmax[8] = {10, 10, 10, 9.5, 8.5, 10, 10, 10};

    for (int i = 0; i < 8; i++)
    {
        auto *h = (R3BFootHitData *)ft->UncheckedAt(i);

        if (h->GetMulStrip() == 1)
            continue;
        if (h->GetDetId() - 1 != i)
            return false;

        double z = h->GetZCharge();
        if (z < Zmin[i] || z > Zmax[i])
            return false;
    }
    return true;
}

struct NFirst
{
    int idx;
    double tns;
    TVector3 pos;
};

////////// MAIN FUNCTION ///////////

void eventFilter(std::string setting = "", TString reaction = "", bool test = false, bool append = false)
{

    ROOT::EnableImplicitMT(25);

    ///////// Select the setting /////////
    const std::string listTxt =
        std::string(getenv("repopath")) + "/final/settings/" + setting;

    std::vector<std::string> files;
    try
    {
        files = read_file_list(listTxt);
    }
    catch (const std::exception &e)
    {
        std::cerr << "[ERROR] " << e.what() << "\n";
        return;
    }

    std::cout << "[OK] Loaded " << files.size() << " ROOT files from " << listTxt << "\n";

    ///////// Select the reaction /////////
    double SEL_AOQ_OUT_MIN = 0;
    double SEL_AOQ_OUT_MAX = 0;
    double SEL_AOQ_IN_MIN = 0;
    double SEL_AOQ_IN_MAX = 0;
    std::string outFileName = "";
    double MASS_AMU = 0;
    bool hasNeutrons = false;

    if (reaction == "25F23O")
    {
        SEL_AOQ_OUT_MIN = 2.8;
        SEL_AOQ_OUT_MAX = 2.88;
        SEL_AOQ_IN_MIN = 2.77;
        SEL_AOQ_IN_MAX = 2.79;

        outFileName = "data_23O";
        MASS_AMU = 23.015696686 - 8 * 0.00511;
        hasNeutrons = true;
    }

    if (reaction == "25F22O")
    {
        SEL_AOQ_OUT_MIN = 2.67;
        SEL_AOQ_OUT_MAX = 2.76;
        SEL_AOQ_IN_MIN = 2.77;
        SEL_AOQ_IN_MAX = 2.79;

        outFileName = "data_22O";
        MASS_AMU = 22.009965744 - 8.0 * 0.00511;
        hasNeutrons = true;
    }

    if (reaction == "25F24O")
    {
        SEL_AOQ_OUT_MIN = 2.9;
        SEL_AOQ_OUT_MAX = 3.;
        SEL_AOQ_IN_MIN = 2.77;
        SEL_AOQ_IN_MAX = 2.79;

        outFileName = "data_24O";
        MASS_AMU = 24.019861000 - 8 * 0.00511;
        hasNeutrons = false;
    }

    if (test)
        outFileName += "_test";

    outFileName += ".root";

    // Other cases
    if (outFileName.empty() || MASS_AMU <= 0.0 || (SEL_AOQ_IN_MAX <= SEL_AOQ_IN_MIN))
    {
        std::cerr << "[ERROR] Unknown/invalid reaction: " << reaction << "\n";
        return;
    }

    double M_FRAG_GeV = MASS_AMU * AMU_GeV;

    /////////// Extract data from the rootfile /////////
    std::string outFile = std::string(getenv("repopath")) +
                          "/results/final/" + outFileName;

    bool exists = (gSystem->AccessPathName(outFile.c_str()) == kFALSE);

    TFile *fout = nullptr;
    if (!append || !exists)
        fout = new TFile(outFile.c_str(), "RECREATE");
    else
        fout = new TFile(outFile.c_str(), "UPDATE");

    ROOT::RDataFrame df("evt", files);

    // Select
    ROOT::RDF::RNode base = df
                                .Filter([](TClonesArray &f)
                                        { return f.GetEntriesFast() > 0; }, {"FrsData"})
                                .Filter([](TClonesArray &m)
                                        { return m.GetEntriesFast() > 0; }, {"FragmentMDFTrack"})
                                .Filter([](TClonesArray &it)
                                        { return it.GetEntriesFast() > 0; }, {"IncomingTrackFoot"})
                                .Filter([](TClonesArray &c)
                                        { return c.GetEntriesFast() >= 2; }, {"CalifaClusterData"});

    ROOT::RDF::RNode df_cut = base;
    if (hasNeutrons)
        df_cut = df_cut.Filter([](TClonesArray &n)
                               { return n.GetEntriesFast() > 0; }, {"NeulandHits"});

    // Select the incoming in the FRS
    auto df_frs = df_cut
                      .Define("in_AoQ", [](TClonesArray &f)
                              { return ((R3BFrsData *)f.UncheckedAt(0))->GetAq(); }, {"FrsData"})
                      .Define("in_Z", [](TClonesArray &f)
                              { return ((R3BFrsData *)f.UncheckedAt(0))->GetZ(); }, {"FrsData"})
                      .Define("beta_proj", [](TClonesArray &f)
                              {
                                  auto *frs = (R3BFrsData *)f.UncheckedAt(0);
                                  return frs->GetBeta(); // ajusta el nombre si es distinto
                              },
                              {"FrsData"})
                      .Filter([&](double aoq, double z)
                              { return isGoodIncoming(aoq, z, SEL_AOQ_IN_MIN, SEL_AOQ_IN_MAX); }, {"in_AoQ", "in_Z"});

    // Kinematics of the fragment from the MDF
    auto df_frag = df_frs
                       .Define("frag", [](TClonesArray &m)
                               { return (R3BTrackingParticle *)m.UncheckedAt(0); },
                               {"FragmentMDFTrack"})

                       // PID from MDF data
                       .Define("AoQ_frag", [](R3BTrackingParticle *t)
                               { return (double)t->GetMass(); },
                               {"frag"})
                       .Define("Z_frag_est", [](R3BTrackingParticle *t)
                               { return (double)t->GetCharge(); },
                               {"frag"})
                       .Define("A_frag", [](double aoq, double z)
                               { return aoq * z; },
                               {"AoQ_frag", "Z_frag_est"})

                       // P/Q from MDF
                       .Define("P_frag_PoQ_vec", [](R3BTrackingParticle *t)
                               {
                                   return t->GetStartMomentum(); // asumido P/Q
                               },
                               {"frag"})
                       .Define("P_over_Q_frag", [](const TVector3 &pvec)
                               { return pvec.Mag(); }, {"P_frag_PoQ_vec"})
                       .Define("P_frag_dir", [](const TVector3 &pvec)
                               { return pvec.Mag() > 0 ? pvec.Unit() : TVector3(0, 0, 1); }, {"P_frag_PoQ_vec"})

                       // Momentum of the incoming proyectile from IncomingFoots
                       .Define("P_in_vec", [](TClonesArray &it)
                               {
                                   auto *trk = (R3BTrackingParticle *)it.UncheckedAt(0);
                                   return trk->GetStartMomentum(); }, {"IncomingTrackFoot"})
                       .Define("px_in", [](const TVector3 &p)
                               { return p.X(); }, {"P_in_vec"})
                       .Define("py_in", [](const TVector3 &p)
                               { return p.Y(); }, {"P_in_vec"})
                       .Define("pz_in", [](const TVector3 &p)
                               { return p.Z(); }, {"P_in_vec"})

                       // Z fixed for the fragment (not the stimated from the MDF!)
                       .Define("Z_frag", []()
                               { return Z_FRAG; }, {})

                       // Fixed mass of the fragment (not the retrived from the MDF!)
                       .Define("M_frag", [&]()
                               { return M_FRAG_GeV; }, {})

                       // Module of the momentum p = (P/Q) * Z
                       .Define("p_frag", [](double poverq, double Z)
                               {
                                   if (poverq <= 0.0 || Z <= 0.0)
                                       return 0.0;
                                   return poverq * Z; }, {"P_over_Q_frag", "Z_frag"})

                       // Beta calculated from p and M (not from FRS)
                       .Define("beta_frag", [](double p, double M)
                               {
                                   if (p <= 0.0 || M <= 0.0)
                                       return 0.0;
                                   double E = std::sqrt(p * p + M * M);
                                   return p / E; }, {"p_frag", "M_frag"})

                       // Momentum 4-vector
                       .Define("P4_frag", [](double p, double Mf, const TVector3 &dir)
                               {
                                   if (p <= 0.0)
                                       return TLorentzVector(0, 0, 0, 0);

                                   TVector3 pvec = dir * p;
                                   double E = std::sqrt(p * p + Mf * Mf);  // E = sqrt(p^2 + M^2)
                                   return TLorentzVector(pvec, E); }, {"p_frag", "M_frag", "P_frag_dir"});

    // Select only the fragment that we want
    auto df_frag_filtered = df_frag.Filter(
        [&](double aoq, double z)
        {
            return (aoq >= SEL_AOQ_OUT_MIN && aoq <= SEL_AOQ_OUT_MAX &&
                    z >= SEL_Z_MIN && z <= SEL_Z_MAX);
        },
        {"AoQ_frag", "Z_frag_est"});

    // Nuetron variables (if any)
    ROOT::RDF::RNode df_p4n = df_frag_filtered;

    if (hasNeutrons)
    {
        // Select the first neutron in Neuland
        auto df_n = df_frag_filtered
                        .Define("n_first", [](TClonesArray &nl)
                                {
                int idx = -1;
                double bestZ = 1e99;
                TVector3 bestPos(0, 0, 0);
                double bestT = -1;

                const int n = nl.GetEntriesFast();
                for (int i = 0; i < n; i++)
                {
                    auto *h = static_cast<R3BNeulandHit *>(nl.UncheckedAt(i));

                    const double z = h->GetPosition().Z();
                    const double t = h->GetT();

                    // First hit in plane inside the time window
                    if (t <= 76.0 && z < bestZ)
                    //if (t >= 0 && z < bestZ)
                    {
                        bestZ  = z;
                        bestT  = t;
                        idx    = i;

                        TVector3 pos = h->GetPosition();
                        TVector3 posSmear = pos;

                        const int paddle = h->GetPaddle();

                        // X and Y are the central value by default in the Cal2clust
                        // task for vertical and horizontal bars respectively.
                        static thread_local TRandom3 rng(0);

                        if ( ((paddle / 50) % 2) == 0 )
                        {
                            // vertical bars smear en Y y Z
                            posSmear.SetXYZ(
                                pos.X(),
                                pos.Y() + rng.Uniform(-2.5, 2.5),
                                pos.Z() + rng.Uniform(-2.5, 2.5)
                            );
                        }
                        else
                        {
                            // horizontal bars smear en X y Z
                            posSmear.SetXYZ(
                                pos.X() + rng.Uniform(-2.5, 2.5),
                                pos.Y(),
                                pos.Z() + rng.Uniform(-2.5, 2.5)
                            );
                        }

                        bestPos = posSmear;
                    }
                }

                return NFirst{idx, bestT, bestPos}; },
                                {"NeulandHits"})
                        .Filter([](const NFirst &nf)
                                { return nf.idx >= 0; },
                                {"n_first"});

        // Calculate beta of the neutron
        auto df_bn = df_n
                         .Define("beta_neu", [](const NFirst &nf)
                                 {
                             double L = nf.pos.Mag();
                             double t = nf.tns;
                             if (L <= 0 || t <= 0)
                                 return -1.0;
                             return L / (C_CM_PER_NS * t); },
                                 {"n_first"})
                         .Filter([](double b)
                                 { return b > 0 && b < 1; },
                                 {"beta_neu"});

        // Neutron direction
        auto df_dir = df_bn.Define(
            "n_dir", [](const NFirst &nf)
            { return nf.pos.Mag() > 0 ? nf.pos.Unit() : TVector3(0, 0, 1); },
            {"n_first"});

        // Calculate 4-momentum of the neutron
        df_p4n = df_dir
                     .Define("p_neu", [](double bn)
                             {
                        double g = 1.0 / std::sqrt(1.0 - bn * bn);
                        return g * MN_GeV * bn; },
                             {"beta_neu"})
                     .Define("P4_neu", [](double p, const TVector3 &dir, double bn)
                             {
                        double g = 1.0 / std::sqrt(1.0 - bn * bn);
                        TVector3 P = dir * p;
                        return TLorentzVector(P, g * MN_GeV); },
                             {"p_neu", "n_dir", "beta_neu"});
    }

    // CALIFA opening angle and cluster angles
    auto df_califa = df_p4n
                         .Define("califa_opa", [](TClonesArray &clu)
                                 {
        double opa = -999.0;

        const double EminClu = 20e3;   // 20 keV prefilter
        const double EminP   = 35e3;   // 35 keV for the selected two

        // top-2 by energy
        double e1=-1., e2=-1.;
        double th1=-999., ph1=-999.;
        double th2=-999., ph2=-999.;

        for (int i = 0; i < clu.GetEntriesFast(); ++i)
        {
            auto *hit = (R3BCalifaClusterData *)clu.UncheckedAt(i);
            const double E = hit->GetEnergy();
            if (E < EminClu) continue;

            const double th = hit->GetTheta();
            const double ph = hit->GetPhi();

            if (E > e1) { e2=e1; th2=th1; ph2=ph1; e1=E; th1=th; ph1=ph; }
            else if (E > e2) { e2=E; th2=th; ph2=ph; }
        }

        if (e1 < EminP || e2 < EminP) return opa;

        const double dphi = TVector2::Phi_mpi_pi(ph2 - ph1);
        const double dphiDeg = std::abs(dphi) * TMath::RadToDeg();
        if (std::abs(dphiDeg - 180.0) > 30.0) return opa;

        TVector3 v1; v1.SetMagThetaPhi(1.0, th1, ph1);
        TVector3 v2; v2.SetMagThetaPhi(1.0, th2, ph2);
        opa = v1.Angle(v2); // rad

        return opa; }, {"CalifaClusterData"})

                         // decide swap ONCE per event (thread-safe with IMT)
                         .Define("califa_swap", []()
                                 {
        static thread_local TRandom3 rng(0);
        return rng.Rndm() < 0.5; })

                         // helper columns: the 2 candidates (top-2) after cuts; return -999 if not valid
                         .Define("califa_th1", [](TClonesArray &clu)
                                 {
        const double EminClu = 20e3, EminP = 35e3;
        double e1=-1., e2=-1., th1=-999., ph1=-999., th2=-999., ph2=-999.;
        for (int i=0;i<clu.GetEntriesFast();++i){
            auto *hit=(R3BCalifaClusterData*)clu.UncheckedAt(i);
            double E=hit->GetEnergy(); if(E<EminClu) continue;
            double th=hit->GetTheta(), ph=hit->GetPhi();
            if(E>e1){ e2=e1; th2=th1; ph2=ph1; e1=E; th1=th; ph1=ph; }
            else if(E>e2){ e2=E; th2=th; ph2=ph; }
        }
        if (e1 < EminP || e2 < EminP) return -999.0;
        double dphiDeg = std::abs(TVector2::Phi_mpi_pi(ph2 - ph1))*TMath::RadToDeg();
        if (std::abs(dphiDeg - 180.0) > 30.0) return -999.0;
        return th1; }, {"CalifaClusterData"})
                         .Define("califa_ph1", [](TClonesArray &clu)
                                 {
        const double EminClu = 20e3, EminP = 35e3;
        double e1=-1., e2=-1., th1=-999., ph1=-999., th2=-999., ph2=-999.;
        for (int i=0;i<clu.GetEntriesFast();++i){
            auto *hit=(R3BCalifaClusterData*)clu.UncheckedAt(i);
            double E=hit->GetEnergy(); if(E<EminClu) continue;
            double th=hit->GetTheta(), ph=hit->GetPhi();
            if(E>e1){ e2=e1; th2=th1; ph2=ph1; e1=E; th1=th; ph1=ph; }
            else if(E>e2){ e2=E; th2=th; ph2=ph; }
        }
        if (e1 < EminP || e2 < EminP) return -999.0;
        double dphiDeg = std::abs(TVector2::Phi_mpi_pi(ph2 - ph1))*TMath::RadToDeg();
        if (std::abs(dphiDeg - 180.0) > 30.0) return -999.0;
        return ph1; }, {"CalifaClusterData"})
                         .Define("califa_th2", [](TClonesArray &clu)
                                 {
        const double EminClu = 20e3, EminP = 35e3;
        double e1=-1., e2=-1., th1=-999., ph1=-999., th2=-999., ph2=-999.;
        for (int i=0;i<clu.GetEntriesFast();++i){
            auto *hit=(R3BCalifaClusterData*)clu.UncheckedAt(i);
            double E=hit->GetEnergy(); if(E<EminClu) continue;
            double th=hit->GetTheta(), ph=hit->GetPhi();
            if(E>e1){ e2=e1; th2=th1; ph2=ph1; e1=E; th1=th; ph1=ph; }
            else if(E>e2){ e2=E; th2=th; ph2=ph; }
        }
        if (e1 < EminP || e2 < EminP) return -999.0;
        double dphiDeg = std::abs(TVector2::Phi_mpi_pi(ph2 - ph1))*TMath::RadToDeg();
        if (std::abs(dphiDeg - 180.0) > 30.0) return -999.0;
        return th2; }, {"CalifaClusterData"})
                         .Define("califa_ph2", [](TClonesArray &clu)
                                 {
        const double EminClu = 20e3, EminP = 35e3;
        double e1=-1., e2=-1., th1=-999., ph1=-999., th2=-999., ph2=-999.;
        for (int i=0;i<clu.GetEntriesFast();++i){
            auto *hit=(R3BCalifaClusterData*)clu.UncheckedAt(i);
            double E=hit->GetEnergy(); if(E<EminClu) continue;
            double th=hit->GetTheta(), ph=hit->GetPhi();
            if(E>e1){ e2=e1; th2=th1; ph2=ph1; e1=E; th1=th; ph1=ph; }
            else if(E>e2){ e2=E; th2=th; ph2=ph; }
        }
        if (e1 < EminP || e2 < EminP) return -999.0;
        double dphiDeg = std::abs(TVector2::Phi_mpi_pi(ph2 - ph1))*TMath::RadToDeg();
        if (std::abs(dphiDeg - 180.0) > 30.0) return -999.0;
        return ph2; }, {"CalifaClusterData"})

                         // Now map to L/R using the SAME swap for all four outputs
                         .Define("califa_theta_L", [](bool swap, double th1, double th2)
                                 {
        if (th1 < -990.0 || th2 < -990.0) return -999.0;
        return swap ? th2 : th1; }, {"califa_swap", "califa_th1", "califa_th2"})
                         .Define("califa_phi_L", [](bool swap, double ph1, double ph2)
                                 {
        if (ph1 < -990.0 || ph2 < -990.0) return -999.0;
        return swap ? ph2 : ph1; }, {"califa_swap", "califa_ph1", "califa_ph2"})
                         .Define("califa_theta_R", [](bool swap, double th1, double th2)
                                 {
        if (th1 < -990.0 || th2 < -990.0) return -999.0;
        return swap ? th1 : th2; }, {"califa_swap", "califa_th1", "califa_th2"})
                         .Define("califa_phi_R", [](bool swap, double ph1, double ph2)
                                 {
        if (ph1 < -990.0 || ph2 < -990.0) return -999.0;
        return swap ? ph1 : ph2; }, {"califa_swap", "califa_ph1", "califa_ph2"});

    auto df_califa_good = df_califa.Filter(
        [](double opa)
        {
            return opa > -990.0;
        },
        {"califa_opa"},
        "good CALIFA event");

    auto df_out = df_califa_good
                      .Define("px_frag", [](const TLorentzVector &F)
                              { return F.Px(); }, {"P4_frag"})
                      .Define("py_frag", [](const TLorentzVector &F)
                              { return F.Py(); }, {"P4_frag"})
                      .Define("pz_frag", [](const TLorentzVector &F)
                              { return F.Pz(); }, {"P4_frag"});

    if (hasNeutrons)
    {
        df_out = df_out
                     .Define("px_neu", [](const TLorentzVector &N)
                             { return N.Px(); }, {"P4_neu"})
                     .Define("py_neu", [](const TLorentzVector &N)
                             { return N.Py(); }, {"P4_neu"})
                     .Define("pz_neu", [](const TLorentzVector &N)
                             { return N.Pz(); }, {"P4_neu"})
                     .Define("x_neu_hit", [](const NFirst &nf)
                             { return nf.pos.X(); }, {"n_first"})
                     .Define("y_neu_hit", [](const NFirst &nf)
                             { return nf.pos.Y(); }, {"n_first"})
                     .Define("z_neu_hit", [](const NFirst &nf)
                             { return nf.pos.Z(); }, {"n_first"})
                     .Define("tof_neuland", [](const NFirst &nf)
                             { return nf.tns; }, {"n_first"});
    }

    // Save the data
    std::vector<std::string> cols = {
        "AoQ_frag",
        "Z_frag_est",
        "A_frag",
        "M_frag",
        "beta_frag",
        "p_frag",
        "califa_opa",
        "califa_theta_L",
        "califa_phi_L",
        "califa_theta_R",
        "califa_phi_R",
        "px_frag",
        "py_frag",
        "pz_frag",
        "beta_proj",
        "px_in",
        "py_in",
        "pz_in"};

    if (hasNeutrons)
    {
        cols = {
            "AoQ_frag",
            "Z_frag_est",
            "A_frag",
            "M_frag",
            "beta_frag",
            "beta_neu",
            "p_frag",
            "p_neu",
            "califa_opa",
            "califa_theta_L",
            "califa_phi_L",
            "califa_theta_R",
            "califa_phi_R",
            "px_frag",
            "py_frag",
            "pz_frag",
            "beta_proj",
            "px_neu",
            "py_neu",
            "pz_neu",
            "x_neu_hit",
            "y_neu_hit",
            "z_neu_hit",
            "tof_neuland",
            "px_in",
            "py_in",
            "pz_in"};
    }

    df_out.Snapshot("FilterDataTree", outFile, cols);

    std::cout << "\n[OK] TTree saved in: " << outFile << "\n";

    if (fout)
        fout->Close();
}
