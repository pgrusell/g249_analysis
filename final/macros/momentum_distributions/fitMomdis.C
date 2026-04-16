struct MomentaDist
{
    TH1F *Qz = nullptr;
    TH1F *Qt = nullptr;
    TH1F *Qy = nullptr;
    TH1F *Q = nullptr;
};

// This function will assume a .txt without any text and with a structure like Qi Qt
MomentaDist getMomentaDistFromtxt(std::string txtFile, std::string histName)
{

    MomentaDist momdis;

    std::vector<double> Qi;
    std::vector<double> Qt;

    ifstream f(txtFile.c_str());

    std::string s;

    while (getline(f, s))
    {
        std::istringstream ss(s);
        double qi, qt;

        if (!(ss >> qi >> qt))
            continue;

        Qi.push_back(qi);
        Qt.push_back(qt);
    }

    const int nBins = Qi.size();
    const double maxBin = Qi[nBins - 1] + (Qi[1] - Qi[0]) / 2.;

    momdis.Qt = new TH1F((std::string("Qt") + histName).c_str(), "Qt", nBins, -maxBin, maxBin);

    for (int i = 0; i < nBins; i++)
    {
        momdis.Qt->SetBinContent(i + 1, Qt[i]);
    }

    return momdis;
}

TH1F *getMomentaDistFromRoot(std::string rootFile, int nBins = 50, double maxBin = 300)
{

    auto *h = new TH1F("hExp", "hExp", nBins, -maxBin, maxBin);

    auto *f = new TFile(rootFile.c_str());
    auto *tr = static_cast<TTree *>(f->Get("KinTree"));

    double py;

    tr->SetBranchAddress("px_rf_rot", &py);

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

    // Normalize to unit area
    double integral = h->Integral();
    if (integral > 0)
        h->Scale(1.0 / integral);

    return h;
}

MomentaDist reampleMomDis(MomentaDist momdisTheo, TH1F *momdisExp, std::string histName, int nSamples = 10000000)
{
    MomentaDist resampled;

    const int nBins = momdisExp->GetNbinsX();
    const int maxVal = momdisExp->GetXaxis()->GetXmax();

    resampled.Qt = new TH1F((std::string("Qt_resampled") + histName).c_str(), "Qt_resampled", nBins, -maxVal, maxVal);

    for (int i = 0; i < nSamples; i++)
        resampled.Qt->Fill(momdisTheo.Qt->GetRandom());

    return resampled;
}

void fitMomdis()
{
    // In files
    std::vector<std::string> inFilesTheo = {
        "/nucl_lustre/pablogrusell/g249/g249_analysis/theory/JT/25F/sigt_s12-gs.txt",
        "/nucl_lustre/pablogrusell/g249/g249_analysis/theory/JT/25F/sigt_d52-gs.txt",
    };

    std::string inFileExp = "/nucl_lustre/pablogrusell/g249/g249_analysis/results/final/24O_analyzed_test.root";

    // Read the histograms from the .txt file
    auto momdis_2s12_gs = getMomentaDistFromtxt(inFilesTheo[0], "2s12_gs");
    auto momdis_1d52_gs = getMomentaDistFromtxt(inFilesTheo[1], "1d52_gs");

    // Read the experimental histogram from the .root file
    auto momdis_exp = getMomentaDistFromRoot(inFileExp);

    // Resample the theoretical histogram to the bins of the experimental one
    auto momdis_2s12_12_resampled = reampleMomDis(momdis_2s12_gs, momdis_exp, "2s12_12", 10000000);
    auto momdis_1d52_15_resampled = reampleMomDis(momdis_1d52_gs, momdis_exp, "1d52_15", 10000000);

    // Fit
    TObjArray *mc = new TObjArray(2);
    mc->Add(momdis_2s12_12_resampled.Qt);
    mc->Add(momdis_1d52_15_resampled.Qt);

    momdis_2s12_12_resampled.Qt->Draw();

    TFractionFitter *fit = new TFractionFitter(momdis_exp, mc);

    fit->Constrain(0, 0.0, 1.0);
    fit->Constrain(1, 0.0, 1.0);

    Int_t status = fit->Fit();

    if (status == 0)
    { // check on fit status
        Double_t chi2 = fit->GetChisquare();
        Int_t ndf = fit->GetNDF();
        std::cout << "chi2 = " << chi2 << "  ndf = " << ndf
                  << "  chi2/ndf = " << chi2 / ndf << std::endl;

        Double_t frac0, err0, frac1, err1;
        fit->GetResult(0, frac0, err0);
        fit->GetResult(1, frac1, err1);
        std::cout << "fraction[0] (2s1/2) = " << frac0 << " +/- " << err0 << std::endl;
        std::cout << "fraction[1] (1d5/2) = " << frac1 << " +/- " << err1 << std::endl;

        TH1F *result = (TH1F *)fit->GetPlot();

        // Build individual components scaled to the fit result
        TH1F *h_s12 = (TH1F *)fit->GetMCPrediction(0);
        TH1F *h_d52 = (TH1F *)fit->GetMCPrediction(1);

        TH1F *comp_s12 = (TH1F *)h_s12->Clone("comp_s12");
        TH1F *comp_d52 = (TH1F *)h_d52->Clone("comp_d52");
        comp_s12->Scale(frac0 * momdis_exp->Integral() / comp_s12->Integral());
        comp_d52->Scale(frac1 * momdis_exp->Integral() / comp_d52->Integral());

        comp_s12->SetLineColor(kBlue);
        comp_s12->SetLineStyle(2);
        comp_s12->SetLineWidth(2);

        comp_d52->SetLineColor(kGreen + 2);
        comp_d52->SetLineStyle(2);
        comp_d52->SetLineWidth(2);

        result->SetLineColor(kRed);
        result->SetLineWidth(2);

        momdis_exp->Draw("Ep");
        result->Draw("same");
        comp_s12->Draw("same hist");
        comp_d52->Draw("same hist");

        auto *leg = new TLegend(0.65, 0.7, 0.89, 0.89);
        leg->AddEntry(momdis_exp, "Exp", "ep");
        leg->AddEntry(result, "Fit total", "l");
        leg->AddEntry(comp_s12, Form("2s_{1/2} (%.1f%%)", frac0 * 100), "l");
        leg->AddEntry(comp_d52, Form("1d_{5/2} (%.1f%%)", frac1 * 100), "l");
        leg->Draw();
    }
}