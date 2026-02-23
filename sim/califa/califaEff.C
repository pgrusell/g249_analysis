
static bool DecodeIonPDG(int pdg, int &Z, int &A)
{
    if (pdg < 1000000000)
        return false;
    int code = pdg - 1000000000;
    Z = code / 10000;
    A = (code % 10000) / 10;
    return (Z > 0 && A > 0);
}

static bool HasO24(TClonesArray *mcTracks)
{
    if (!mcTracks)
        return false;
    for (int i = 0; i < mcTracks->GetEntriesFast(); ++i)
    {
        auto *trk = static_cast<R3BMCTrack *>(mcTracks->UncheckedAt(i));
        if (!trk)
            continue;
        int Z = 0, A = 0;
        if (!DecodeIonPDG(trk->GetPdgCode(), Z, A))
            continue;
        if (Z == 8 && A == 24)
            return true;
    }
    return false;
}

static bool HasAllFourFoots(TClonesArray *traPoints)
{
    if (!traPoints)
        return false;
    bool has5 = false, has6 = false, has7 = false, has8 = false;
    for (int i = 0; i < traPoints->GetEntriesFast(); ++i)
    {
        auto *p = static_cast<R3BTraPoint *>(traPoints->UncheckedAt(i));
        if (!p)
            continue;
        int id = p->GetDetectorID();
        if (id == 5)
            has5 = true;
        else if (id == 6)
            has6 = true;
        else if (id == 7)
            has7 = true;
        else if (id == 8)
            has8 = true;
    }
    return (has5 && has6 && has7 && has8);
}

struct OAResult
{
    bool ok = false;
    double opa = -999.;
    double th1 = -999., ph1 = -999.;
    double th2 = -999., ph2 = -999.;
};

static OAResult ComputeOpeningAngle(TClonesArray *clu,
                                    double EminClu = 20,
                                    double EminP = 20,
                                    double dphiWinDeg = 30.0)
{
    OAResult out;
    if (!clu)
        return out;

    double e1 = -1., e2 = -1.;
    double th1 = -999., ph1 = -999.;
    double th2 = -999., ph2 = -999.;

    for (int i = 0; i < clu->GetEntriesFast(); ++i)
    {
        auto *hit = static_cast<R3BCalifaClusterData *>(clu->UncheckedAt(i));
        if (!hit)
            continue;
        double E = hit->GetEnergy();
        if (E < EminClu)
            continue;

        double th = hit->GetTheta();
        double ph = hit->GetPhi();

        if (E > e1)
        {
            e2 = e1;
            th2 = th1;
            ph2 = ph1;
            e1 = E;
            th1 = th;
            ph1 = ph;
        }
        else if (E > e2)
        {
            e2 = E;
            th2 = th;
            ph2 = ph;
        }
    }

    if (e1 < EminP || e2 < EminP)
        return out;

    double dphi = TVector2::Phi_mpi_pi(ph2 - ph1);
    double dphiDeg = std::abs(dphi) * TMath::RadToDeg();
    if (std::abs(dphiDeg - 180.0) > dphiWinDeg)
        return out;

    TVector3 v1;
    v1.SetMagThetaPhi(1.0, th1, ph1);
    TVector3 v2;
    v2.SetMagThetaPhi(1.0, th2, ph2);

    out.ok = true;
    out.opa = v1.Angle(v2) * TMath::RadToDeg();
    out.th1 = th1;
    out.ph1 = ph1;
    out.th2 = th2;
    out.ph2 = ph2;
    return out;
}

struct OpaFitResult
{
    double signalYield = 0.0;
    double mu = 0.0;
    double sigma = 0.0;
    double chi2NDF = 0.0;
    int fitStatus = -1;
};

static OpaFitResult FitOpaDistribution(TH1F *h, bool verbose = false)
{
    OpaFitResult res;

    if (!h || h->GetEntries() < 10)
    {
        if (verbose)
            std::cout << "[FitOpaDistribution] Too few entries.\n";
        return res;
    }

    const double center = h->GetBinCenter(h->GetMaximumBin());
    const double binWidth = h->GetXaxis()->GetBinWidth(1);

    TF1 *func = new TF1("funcOpa",
                        "[0] + [1]*x + [2]*exp(-0.5*((x-[3])/[4])^2)",
                        75,
                        90);

    func->SetParameters(10.0, 0.0, h->GetMaximum(), center, 5.0, 0.0);
    func->SetParLimits(2, 0.0, 1e9);
    func->SetParLimits(3, 60.0, 100.0);
    func->SetParLimits(4, 1.0, 20.0);

    int fitSt = h->Fit(func, "RQ");

    const double A = func->GetParameter(2);
    const double mu = func->GetParameter(3);
    const double sg = std::abs(func->GetParameter(4));

    res.mu = mu;
    res.sigma = sg;
    res.fitStatus = fitSt;
    if (func->GetNDF() > 0)
        res.chi2NDF = func->GetChisquare() / func->GetNDF();

    res.signalYield = A * sg * std::sqrt(2.0 * TMath::Pi()) / binWidth;

    if (verbose)
    {
        std::cout << "\n=== OPA FIT (simulation) ===\n";
        std::cout << "Fit status : " << fitSt << "\n";
        std::cout << "chi2/NDF   : " << res.chi2NDF << "\n";
        std::cout << "mu (deg)   : " << mu << "\n";
        std::cout << "sigma (deg): " << sg << "\n";
        std::cout << "A          : " << A << "\n";
        std::cout << "signalYield: " << res.signalYield << "\n";
        std::cout << "============================\n\n";
    }

    delete func;
    return res;
}

void califaEff(const char *inFile = "sim.root",
               const char *outFile = "califaEff_out.root",
               Long64_t maxEvents = -1,
               bool verbose = false)
{

    std::cout << "vamos\n";

    TFile *fin = TFile::Open(inFile, "READ");
    if (!fin || fin->IsZombie())
    {
        std::cerr << "[ERROR] Cannot open " << inFile << "\n";
        return;
    }

    TTree *t = dynamic_cast<TTree *>(fin->Get("evt"));
    if (!t)
    {
        std::cerr << "[ERROR] Tree 'evt' not found.\n";
        fin->ls();
        return;
    }

    TClonesArray *mcTracks = nullptr;
    TClonesArray *traPoints = nullptr;
    TClonesArray *califaClusters = nullptr;

    t->SetBranchAddress("MCTrack", &mcTracks);
    t->SetBranchAddress("TraPoint", &traPoints);
    t->SetBranchAddress("CalifaClusterData", &califaClusters);

    TFile *fout = TFile::Open(outFile, "RECREATE");

    int evtId = -1;
    double opa = -999., th1 = -999., ph1 = -999., th2 = -999., ph2 = -999.;

    TTree *tout = new TTree("goodO24",
                            "24O events with 4 FOOTs + 2 CALIFA clusters + valid OA");
    tout->Branch("event", &evtId, "event/I");
    tout->Branch("opa", &opa, "opa/D");
    tout->Branch("theta1", &th1, "theta1/D");
    tout->Branch("phi1", &ph1, "phi1/D");
    tout->Branch("theta2", &th2, "theta2/D");
    tout->Branch("phi2", &ph2, "phi2/D");

    TH1F *hOpa = new TH1F("hOpa", "OPA distribution (sim);OPA (deg);Counts",
                          100, 60, 100.0);
    hOpa->SetDirectory(nullptr);

    const Long64_t nEntries = t->GetEntries();
    const Long64_t nToRead = (maxEvents > 0) ? std::min(maxEvents, nEntries) : nEntries;

    Long64_t nO24 = 0;
    Long64_t nO24_foots = 0;
    Long64_t nO24_foots_2clu = 0;
    Long64_t nO24_foots_2clu_oa = 0;

    for (Long64_t ie = 0; ie < nToRead; ++ie)
    {
        t->GetEntry(ie);

        if (!HasO24(mcTracks))
            continue;
        nO24++;

        if (!HasAllFourFoots(traPoints))
            continue;
        nO24_foots++;

        if (!califaClusters || califaClusters->GetEntriesFast() < 2)
            continue;
        nO24_foots_2clu++;

        OAResult res = ComputeOpeningAngle(califaClusters);
        if (!res.ok)
            continue;
        nO24_foots_2clu_oa++;

        hOpa->Fill(res.opa);

        evtId = (int)ie;
        opa = res.opa;
        th1 = res.th1;
        ph1 = res.ph1;
        th2 = res.th2;
        ph2 = res.ph2;
        tout->Fill();
    }

    OpaFitResult fitRes = FitOpaDistribution(hOpa, verbose);

    std::cout << "\n========== califaEff summary ==========\n";
    std::cout << "Input file:   " << inFile << "\n";
    std::cout << "Events read:  " << nToRead << " / " << nEntries << "\n";
    std::cout << "---------------------------------------\n";
    std::cout << "n24O (generated):                  " << nO24 << "\n";
    std::cout << "n24O + 4 FOOTs:                    " << nO24_foots << "\n";
    std::cout << "n24O + 4 FOOTs + >=2 clusters:     " << nO24_foots_2clu << "\n";
    std::cout << "n24O + 4 FOOTs + 2clu + valid OA:  " << nO24_foots_2clu_oa << "\n";
    std::cout << "---------------------------------------\n";
    std::cout << "--- Conteo directo (método antiguo) ---\n";
    std::cout << "Eff p2p (conteo): "
              << (double)nO24_foots_2clu_oa / (double)nO24_foots << "\n";
    std::cout << "--- Fit gaussiana + fondo (método análisis) ---\n";
    std::cout << "OPA fit mu    : " << fitRes.mu << " deg\n";
    std::cout << "OPA fit sigma : " << fitRes.sigma << " deg\n";
    std::cout << "OPA fit chi2/NDF: " << fitRes.chi2NDF << "\n";
    std::cout << "Signal yield (fit): " << fitRes.signalYield << "\n";
    std::cout << "Eff p2p (fit):      "
              << fitRes.signalYield / (double)nO24_foots << "\n";
    std::cout << "=======================================\n\n";

    fout->cd();
    hOpa->Write();
    tout->Write();
    fout->Close();
    fin->Close();
}