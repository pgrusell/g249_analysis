struct MomentaDist
{
    TH1F *Qz = nullptr;
    TH1F *Qt = nullptr;
    TH1F *Qy = nullptr;
    TH1F *Q = nullptr;
};

// This function will assume a .txt without any text and with a structure like Qi Qz Qt Qy Q
MomentaDist getMomentaDistFromtxt(std::string txtFile, std::string histName)
{

    MomentaDist momdis;

    std::vector<double> Qi;
    std::vector<double> Qz;
    std::vector<double> Qt;
    std::vector<double> Qy;
    std::vector<double> Q;

    ifstream f(txtFile.c_str());

    std::string s;

    while (getline(f, s))
    {
        std::istringstream ss(s);
        double qi, qz, qt, qy, q;

        if (!(ss >> qi >> qz >> qt >> qy >> q))
            continue;

        Qi.push_back(qi);
        Qz.push_back(qz);
        Qt.push_back(qt);
        Qy.push_back(qy);
        Q.push_back(q);
    }

    const int nBins = Qi.size();
    const double maxBin = Qi[nBins - 1] + (Qi[1] - Qi[0]) / 2.;

    momdis.Qz = new TH1F((std::string("Qz") + histName).c_str(), "Qz", nBins, -maxBin, maxBin);
    momdis.Qt = new TH1F((std::string("Qt") + histName).c_str(), "Qt", nBins, -maxBin, maxBin);
    momdis.Qy = new TH1F((std::string("Qy") + histName).c_str(), "Qy", nBins, -maxBin, maxBin);
    momdis.Q = new TH1F((std::string("Q") + histName).c_str(), "Q", nBins, -maxBin, maxBin);

    for (int i = 0; i < nBins; i++)
    {
        momdis.Qz->SetBinContent(i + 1, Qz[i]);
        momdis.Qy->SetBinContent(i + 1, Qy[i]);

        if (Qi[i] > 0)
        {
            momdis.Qt->SetBinContent(i + 1, Qt[i]);
            momdis.Q->SetBinContent(i + 1, Q[i]);
        }
    }

    return momdis;
}

TH1F *getMomentaDistFromRoot(std::string rootFile, int nBins = 50, double maxBin = 300)
{

    auto *h = new TH1F("hExp", "hExp", nBins, -maxBin, maxBin);

    auto *f = new TFile(rootFile.c_str());
    auto *tr = static_cast<TTree *>(f->Get("KinTree"));

    double py;

    tr->SetBranchAddress("py_rf_rot", &py);

    double sum = 0.0;
    Long64_t n = 0;

    const int nEntries = tr->GetEntries();

    for (Long64_t i = 0; i < nEntries; ++i)
    {
        tr->GetEntry(i);
        sum += py * 1000.0;
        ++n;
    }

    const double offset = sum / n;

    for (int i = 0; i < nEntries; ++i)
    {
        tr->GetEntry(i);

        h->Fill(py * 1000 - offset);
    }

    return h;
}

MomentaDist reampleMomDis(MomentaDist momdisTheo, TH1F *momdisExp, std::string histName, int nSamples = 10000000)
{
    MomentaDist resampled;

    const int nBins = momdisExp->GetNbinsX();
    const int maxVal = momdisExp->GetXaxis()->GetXmax();

    resampled.Qz = new TH1F((std::string("Qz_resampled") + histName).c_str(), "Qz_resampled", nBins, -maxVal, maxVal);
    resampled.Qy = new TH1F((std::string("Qy_resampled") + histName).c_str(), "Qy_resampled", nBins, -maxVal, maxVal);
    resampled.Qt = new TH1F((std::string("Qt_resampled") + histName).c_str(), "Qt_resampled", nBins, -maxVal, maxVal);
    resampled.Q = new TH1F((std::string("Q_resampled") + histName).c_str(), "Q_resampled", nBins, -maxVal, maxVal);

    for (int i = 0; i < nSamples; i++)
    {
        resampled.Qz->Fill(momdisTheo.Qz->GetRandom());
        resampled.Qy->Fill(momdisTheo.Qy->GetRandom());
        resampled.Qt->Fill(momdisTheo.Qt->GetRandom());
        resampled.Q->Fill(momdisTheo.Q->GetRandom());
    }

    return resampled;
}

void fitMomentaDist()
{
    // In files
    std::vector<std::string> inFilesTheo = {
        "/nucl_lustre/pablogrusell/g249/g249_analysis/theory/25F(p,2p)-2s12-Sp-35.392.txt",
        "/nucl_lustre/pablogrusell/g249/g249_analysis/theory/25F(p,2p)-2s12-Sp-33.000.txt",
        "/nucl_lustre/pablogrusell/g249/g249_analysis/theory/25F(p,2p)-2s12-Sp-30.930.txt",

        "/nucl_lustre/pablogrusell/g249/g249_analysis/theory/25F(p,2p)-2s12-Sp-29.850.txt",
        "/nucl_lustre/pablogrusell/g249/g249_analysis/theory/25F(p,2p)-1p32-Sp-35.392.txt",
        "/nucl_lustre/pablogrusell/g249/g249_analysis/theory/25F(p,2p)-1p32-Sp-33.000.txt",

        "/nucl_lustre/pablogrusell/g249/g249_analysis/theory/25F(p,2p)-1p12-Sp-35.392.txt",
        "/nucl_lustre/pablogrusell/g249/g249_analysis/theory/25F(p,2p)-1p12-Sp-33.000.txt",
        "/nucl_lustre/pablogrusell/g249/g249_analysis/theory/25F(p,2p)-1p12-Sp-30.930.txt",

        "/nucl_lustre/pablogrusell/g249/g249_analysis/theory/25F(p,2p)-1d52-Sp-35.392.txt",
        "/nucl_lustre/pablogrusell/g249/g249_analysis/theory/25F(p,2p)-1d52-Sp-33.000.txt",
        "/nucl_lustre/pablogrusell/g249/g249_analysis/theory/25F(p,2p)-1d52-Sp-30.930.txt",

        "/nucl_lustre/pablogrusell/g249/g249_analysis/theory/25F(p,2p)-1d52-Sp-29.850.txt",
        "/nucl_lustre/pablogrusell/g249/g249_analysis/theory/25F(p,2p)-1d52-Sp-25.509.txt",
    };

    std::string inFileExp = "/nucl_lustre/pablogrusell/g249/g249_analysis/results/final/24O_analyzed_test.root";

    // Read the histograms from the .txt file
    auto momdis_2s12_35 = getMomentaDistFromtxt(inFilesTheo[0], "2s12_35");
    auto momdis_2s12_33 = getMomentaDistFromtxt(inFilesTheo[1], "2s12_33");
    auto momdis_2s12_25 = getMomentaDistFromtxt(inFilesTheo[2], "2s12_25");

    auto momdis_2s12_29 = getMomentaDistFromtxt(inFilesTheo[3], "2s12_29");

    auto momdis_1p32_35 = getMomentaDistFromtxt(inFilesTheo[4], "1p32_35");
    auto momdis_1p32_33 = getMomentaDistFromtxt(inFilesTheo[5], "1p32_33");

    auto momdis_1p12_35 = getMomentaDistFromtxt(inFilesTheo[6], "1p12_35");
    auto momdis_1p12_33 = getMomentaDistFromtxt(inFilesTheo[7], "1p12_33");
    auto momdis_1p12_30 = getMomentaDistFromtxt(inFilesTheo[8], "1p12_30");

    auto momdis_1d52_35 = getMomentaDistFromtxt(inFilesTheo[9], "1d52_35");
    auto momdis_1d52_33 = getMomentaDistFromtxt(inFilesTheo[10], "1d52_33");
    auto momdis_1d52_30 = getMomentaDistFromtxt(inFilesTheo[11], "1d52_30");

    auto momdis_1d52_29 = getMomentaDistFromtxt(inFilesTheo[12], "1d52_29");
    auto momdis_1d52_25 = getMomentaDistFromtxt(inFilesTheo[13], "1d52_25");

    // Read the experimental histogram from the .root file
    auto momdis_exp = getMomentaDistFromRoot(inFileExp);

    // Resample the theoretical histogram to the bins of the experimental one
    auto momdis_2s12_35_resampled = reampleMomDis(momdis_2s12_35, momdis_exp, "2s12_35");
    auto momdis_2s12_33_resampled = reampleMomDis(momdis_2s12_33, momdis_exp, "2s12_33");
    auto momdis_2s12_29_resampled = reampleMomDis(momdis_2s12_29, momdis_exp, "2s12_29");
    auto momdis_2s12_25_resampled = reampleMomDis(momdis_2s12_25, momdis_exp, "2s12_25");

    auto momdis_1p32_35_resampled = reampleMomDis(momdis_1p32_35, momdis_exp, "1p32_35");
    auto momdis_1p32_33_resampled = reampleMomDis(momdis_1p32_33, momdis_exp, "1p32_33");

    auto momdis_1p12_35_resampled = reampleMomDis(momdis_1p12_35, momdis_exp, "1p12_35");
    auto momdis_1p12_33_resampled = reampleMomDis(momdis_1p12_33, momdis_exp, "1p12_33");
    auto momdis_1p12_30_resampled = reampleMomDis(momdis_1p12_30, momdis_exp, "1p12_30");

    auto momdis_1d52_35_resampled = reampleMomDis(momdis_1d52_35, momdis_exp, "1d52_35");
    auto momdis_1d52_33_resampled = reampleMomDis(momdis_1d52_33, momdis_exp, "1d52_33");
    auto momdis_1d52_30_resampled = reampleMomDis(momdis_1d52_30, momdis_exp, "1d52_30");
    auto momdis_1d52_29_resampled = reampleMomDis(momdis_1d52_29, momdis_exp, "1d52_29");
    auto momdis_1d52_25_resampled = reampleMomDis(momdis_1d52_25, momdis_exp, "1d52_25");

    // Fit
    TObjArray *mc = new TObjArray(2);
    mc->Add(momdis_1d52_29_resampled.Qy);
    mc->Add(momdis_1d52_25_resampled.Qy);

    // mc->Add(momdis_1p12_30_resampled.Qy);
    // mc->Add(momdis_1p12_30_resampled.Qy);

    TFractionFitter *fit = new TFractionFitter(momdis_exp, mc);

    fit->Constrain(0, 0.0, 1.0);
    fit->Constrain(1, 0.0, 1.0);
    // fit->Constrain(2, 0.0, 1.0);
    // fit->Constrain(3, 0.0, 1.0);

    Int_t status = fit->Fit();
    std::cout << "fit status: " << status << std::endl;

    if (status == 0)
    { // check on fit status
        TH1F *result = (TH1F *)fit->GetPlot();
        momdis_exp->Draw("Ep");
        result->Draw("same");
    }
}