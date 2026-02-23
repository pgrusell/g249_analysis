struct Top2Clusters
{
    double e1 = -1, th1 = -999, ph1 = -999;
    double e2 = -1, th2 = -999, ph2 = -999;
    bool good = false;
};

Top2Clusters findTop2Clusters(TClonesArray &clu,
                              double EminClu = 20,
                              double EminP = 20)
{
    Top2Clusters c;

    for (int i = 0; i < clu.GetEntriesFast(); ++i)
    {
        auto *hit = (R3BCalifaClusterData *)clu.UncheckedAt(i);
        double E = hit->GetEnergy();
        if (E < EminClu)
            continue;
        double th = hit->GetTheta(), ph = hit->GetPhi();

        if (E > c.e1)
        {
            c.e2 = c.e1;
            c.th2 = c.th1;
            c.ph2 = c.ph1;
            c.e1 = E;
            c.th1 = th;
            c.ph1 = ph;
        }
        else if (E > c.e2)
        {
            c.e2 = E;
            c.th2 = th;
            c.ph2 = ph;
        }
    }

    if (c.e1 < EminP || c.e2 < EminP)
        return c;

    double dphiDeg = std::abs(TVector2::Phi_mpi_pi(c.ph2 - c.ph1)) * TMath::RadToDeg();
    if (std::abs(dphiDeg - 180.0) > 30.0)
        return c;

    c.good = true;

    if (gRandom->Uniform() > 0.5)
    {
        std::swap(c.e1, c.e2);
        std::swap(c.th1, c.th2);
        std::swap(c.ph1, c.ph2);
    }

    return c;
}

struct OpaCutResults
{
    double roiMin = 0.0, roiMax = 0.0;
    double mu = 0.0, sigma = 0.0, kRoi = 0.0;
    double eventsRoi = 0.0, signalInRoi = 0.0, backgroundInRoi = 0.0;
    double signalTotal = 0.0;
    double efficiency = 0.0, purity = 0.0;
    double selectedSignal = 0.0, totalSignalFromRoi = 0.0;
};

OpaCutResults fitOpaAndExtract(TH1F *h, double kRoi, bool verbose = false)
{
    OpaCutResults res;
    res.kRoi = kRoi;

    if (h->GetEntries() == 0)
    {
        if (verbose)
            std::cout << "[WARN] OPA histogram is empty.\n";
        return res;
    }

    const double center = h->GetBinCenter(h->GetMaximumBin());

    auto *func = new TF1("funcOpa",
                         "[0] + [1]*x + [2]*exp(-0.5*((x-[3])/[4])^2) + [5]*x*x",
                         0.6, 2.0);

    func->SetParameters(10.0, 0.0, 1000.0, center, 0.08, 0.0);
    func->SetParLimits(2, 0.0, 1e9);
    func->SetParLimits(3, 0.4, 2.5);
    func->SetParLimits(4, 0.01, 0.5);

    int fitStatus = h->Fit(func, "RQS");

    const double p0 = func->GetParameter(0);
    const double p1 = func->GetParameter(1);
    const double A = func->GetParameter(2);
    const double mu = func->GetParameter(3);
    const double sg = std::abs(func->GetParameter(4));
    const double p2 = func->GetParameter(5);

    res.mu = mu;
    res.sigma = sg;
    res.roiMin = mu - kRoi * sg;
    res.roiMax = mu + kRoi * sg;

    if (res.roiMin < h->GetXaxis()->GetXmin())
        res.roiMin = h->GetXaxis()->GetXmin();
    if (res.roiMax > h->GetXaxis()->GetXmax())
        res.roiMax = h->GetXaxis()->GetXmax();

    auto *sig = new TF1("sigOpa", "[0]*exp(-0.5*((x-[1])/[2])^2)", 0.0, 3.0);
    sig->SetParameters(A, mu, sg);

    auto *bkg = new TF1("bkgOpa", "[0] + [1]*x + [2]*x*x", 0.0, 3.0);
    bkg->SetParameters(p0, p1, p2);

    const double binWidth = h->GetXaxis()->GetBinWidth(1);

    res.signalInRoi = sig->Integral(res.roiMin, res.roiMax) / binWidth;
    res.backgroundInRoi = bkg->Integral(res.roiMin, res.roiMax) / binWidth;
    if (res.backgroundInRoi < 0.0)
        res.backgroundInRoi = 0.0;

    const int bin1 = h->FindBin(res.roiMin);
    const int bin2 = h->FindBin(res.roiMax);
    res.eventsRoi = h->Integral(bin1, bin2);

    res.signalTotal = A * sg * std::sqrt(2.0 * TMath::Pi()) / binWidth;

    res.efficiency = (res.signalTotal > 0.0) ? res.signalInRoi / res.signalTotal : 0.0;

    const double denom = res.signalInRoi + res.backgroundInRoi;
    res.purity = (denom > 0.0) ? res.signalInRoi / denom : 0.0;

    res.selectedSignal = res.eventsRoi * res.purity;
    res.totalSignalFromRoi = (res.efficiency > 0.0)
                                 ? res.selectedSignal / res.efficiency
                                 : 0.0;

    if (verbose)
    {
        std::cout << "\n=== OPA FIT RESULTS ===\n";
        std::cout << "Fit status            : " << fitStatus << "\n";
        std::cout << "mu                    : " << res.mu << "\n";
        std::cout << "sigma                 : " << res.sigma << "\n";
        std::cout << "ROI                   : [" << res.roiMin << ", " << res.roiMax << "]\n";
        std::cout << "Events in ROI (data)  : " << res.eventsRoi << "\n";
        std::cout << "Signal in ROI (fit)   : " << res.signalInRoi << "\n";
        std::cout << "Background in ROI     : " << res.backgroundInRoi << "\n";
        std::cout << "Signal total (Gauss)  : " << res.signalTotal << "\n";
        std::cout << "ROI efficiency        : " << res.efficiency << "\n";
        std::cout << "ROI purity            : " << res.purity << "\n";
        std::cout << "Selected signal       : " << res.selectedSignal << "\n";
        std::cout << "Total signal from ROI : " << res.totalSignalFromRoi << "\n";
        std::cout << "=======================\n\n";
    }

    delete func;
    delete sig;
    delete bkg;

    return res;
}

static constexpr int PDG_24O = 1000080240;
static constexpr int PDG_PROTON = 2212;

double analyzeOneFile(const char *filename, double kRoi = 3.5, bool verbose = true)
{
    TFile *f = TFile::Open(filename, "READ");
    if (!f || f->IsZombie())
    {
        std::cerr << "[ERROR] Cannot open " << filename << "\n";
        return -1.0;
    }

    TTree *tree = (TTree *)f->Get("evt");
    if (!tree)
    {
        std::cerr << "[ERROR] TTree 'evt' not found in " << filename << "\n";
        f->Close();
        return -1.0;
    }

    const Long64_t nEntries = tree->GetEntries();
    if (verbose)
        std::cout << "[INFO] " << filename << " — entries: " << nEntries << "\n";

    TClonesArray *mcTracks = nullptr;
    TClonesArray *califaClusters = nullptr;
    tree->SetBranchAddress("MCTrack", &mcTracks);
    tree->SetBranchAddress("CalifaClusterData", &califaClusters);

    Long64_t nSelected24O = 0, nGenP2P = 0, nGoodCalifa = 0;

    TString hName = TString::Format("hOpa_%s", gSystem->BaseName(filename));
    auto *hOpa = new TH1F(hName, "Opening angle;#theta_{opa} [rad];Counts",
                          100, 0, 3);

    for (Long64_t iEntry = 0; iEntry < nEntries; ++iEntry)
    {
        tree->GetEntry(iEntry);

        if (!mcTracks || mcTracks->GetEntriesFast() == 0)
            continue;

        // 1) Check for 24O
        bool has24O = false;
        for (int i = 0; i < mcTracks->GetEntriesFast(); ++i)
        {
            if (((R3BMCTrack *)mcTracks->UncheckedAt(i))->GetPdgCode() == PDG_24O)
            {
                has24O = true;
                break;
            }
        }
        if (!has24O)
            continue;
        nSelected24O++;

        // 2) Count primary protons => p2p tag
        int nPrimaryProtons = 0;
        for (int i = 0; i < mcTracks->GetEntriesFast(); ++i)
        {
            auto *trk = (R3BMCTrack *)mcTracks->UncheckedAt(i);
            if (trk->GetPdgCode() == PDG_PROTON && trk->GetMotherId() == -1)
                nPrimaryProtons++;
        }
        bool isP2P = (nPrimaryProtons >= 2);
        if (isP2P)
            nGenP2P++;

        // 3) CALIFA reconstruction
        if (!califaClusters || califaClusters->GetEntriesFast() < 2)
            continue;

        Top2Clusters top2 = findTop2Clusters(*califaClusters);
        if (!top2.good)
            continue;

        TVector3 v1, v2;
        v1.SetMagThetaPhi(1.0, top2.th1, top2.ph1);
        v2.SetMagThetaPhi(1.0, top2.th2, top2.ph2);
        double opa = v1.Angle(v2);

        hOpa->Fill(opa);
        if (isP2P)
            nGoodCalifa++;
    }

    // 4) Fit OPA
    OpaCutResults opaRes = fitOpaAndExtract(hOpa, kRoi, verbose);
    double measP2P = opaRes.totalSignalFromRoi;

    double effP2P = (nGenP2P > 0) ? measP2P / (double)nGenP2P : 0.0;
    double effRaw = (nGenP2P > 0) ? (double)nGoodCalifa / (double)nGenP2P : 0.0;

    if (verbose)
    {
        std::cout << "  Events with 24O              : " << nSelected24O << "\n";
        std::cout << "  Generated p2p                : " << nGenP2P << "\n";
        std::cout << "  p2p passing CALIFA (truth)   : " << nGoodCalifa << "\n";
        std::cout << "  Measured p2p (OPA corrected)  : " << measP2P << "\n";
        std::cout << "  eff(p2p) raw                   : " << effRaw << "\n";
        std::cout << "  eff(p2p) fit                   : " << effP2P << "\n";
    }

    delete hOpa;
    f->Close();
    delete f;

    return effP2P;
}

void loopAll(const char *basedir = "/nucl_lustre/pablogrusell/g249/g249_analysis/sim/califa/",
             double kRoi = 3.5,
             bool verbose = true)
{
    const int nFiles = 45;

    auto *hEff = new TH1F("hEff",
                          "CALIFA p2p efficiency per file;"
                          "#varepsilon(p2p);Number of files",
                          50, 0.0, 1.0);
    hEff->SetLineColor(kBlue);
    hEff->SetLineWidth(2);
    hEff->SetFillColor(kBlue - 9);
    hEff->SetFillStyle(1001);

    auto *grEff = new TGraphErrors();
    grEff->SetName("grEff");
    grEff->SetTitle("CALIFA p2p efficiency vs file index;"
                    "File index i;#varepsilon(p2p)");
    grEff->SetMarkerStyle(20);
    grEff->SetMarkerSize(0.9);
    grEff->SetMarkerColor(kBlue);
    grEff->SetLineColor(kBlue);

    std::vector<double> effValues;
    int nGood = 0;

    for (int i = 0; i < nFiles; ++i)
    {
        TString fname = TString::Format("%s/sim_%d.root", basedir, i);

        std::cout << "  Processing file " << i << " / " << nFiles - 1
                  << " : " << fname << "\n";

        double eff = analyzeOneFile(fname.Data(), kRoi, verbose);

        if (eff < 0.0)
        {
            std::cerr << "[WARN] Skipping file " << i << "\n";
            continue;
        }

        hEff->Fill(eff);
        grEff->SetPoint(nGood, (double)i, eff);
        effValues.push_back(eff);
        nGood++;
    }

    double meanEff = 0.0, rmsEff = 0.0;
    if (!effValues.empty())
    {
        for (auto v : effValues)
            meanEff += v;
        meanEff /= effValues.size();
        for (auto v : effValues)
            rmsEff += (v - meanEff) * (v - meanEff);
        rmsEff = std::sqrt(rmsEff / effValues.size());
    }

    std::cout << "  SUMMARY: CALIFA p2p EFFICIENCY OVER " << nGood << " FILES\n";
    std::cout << "  Mean eff(p2p)   : " << meanEff << "\n";
    std::cout << "  RMS           : " << rmsEff << "\n";

    if (!effValues.empty())
    {
        std::cout << "  Min           : " << *std::min_element(effValues.begin(), effValues.end()) << "\n";
        std::cout << "  Max           : " << *std::max_element(effValues.begin(), effValues.end()) << "\n";
    }
    std::cout << "═══════════════════════════════════════════════════════════════\n\n";

    auto *c1 = new TCanvas("cEffHist", "p2p Efficiency Distribution", 800, 600);
    hEff->Draw();

    auto *lineMean = new TLine(meanEff, 0, meanEff, hEff->GetMaximum() * 1.05);
    lineMean->SetLineColor(kRed);
    lineMean->SetLineStyle(2);
    lineMean->SetLineWidth(2);
    lineMean->Draw("same");

    auto *leg1 = new TLegend(0.55, 0.75, 0.88, 0.88);
    leg1->AddEntry(hEff, "Efficiency per file", "f");
    leg1->AddEntry(lineMean, TString::Format("Mean = %.4f", meanEff), "l");
    leg1->SetBorderSize(0);
    leg1->Draw();
    c1->SaveAs("califaEffp2p_distribution.png");

    auto *c2 = new TCanvas("cEffVsFile", "p2p Efficiency vs File Index", 900, 600);
    grEff->Draw("AP");

    auto *band = new TBox(0, meanEff - rmsEff, nFiles - 1, meanEff + rmsEff);
    band->SetFillColorAlpha(kRed, 0.15);
    band->SetLineColor(kRed);
    band->SetLineStyle(2);
    band->Draw("same");

    auto *lineMean2 = new TLine(0, meanEff, nFiles - 1, meanEff);
    lineMean2->SetLineColor(kRed);
    lineMean2->SetLineStyle(2);
    lineMean2->SetLineWidth(2);
    lineMean2->Draw("same");

    grEff->Draw("P same");

    auto *leg2 = new TLegend(0.55, 0.75, 0.88, 0.88);
    leg2->AddEntry(grEff, "#varepsilon(p2p) per file", "p");
    leg2->AddEntry(lineMean2, TString::Format("Mean = %.4f", meanEff), "l");
    leg2->AddEntry(band, TString::Format("RMS = %.4f", rmsEff), "f");
    leg2->SetBorderSize(0);
    leg2->Draw();
    c2->SaveAs("califaEffp2p_vs_fileindex.png");

    std::cout << "[INFO] Plots saved.\n";
}

void califaEffp2p(const char *filename = "/nucl_lustre/pablogrusell/g249/g249_analysis/sim/califa/sim1234.root",
                  double kRoi = 3.5,
                  bool verbose = true)
{
    double eff = analyzeOneFile(filename, kRoi, verbose);
}