struct MomentaDist
{
    TH1F *Qz = nullptr;
    TH1F *Qt = nullptr;
    TH1F *Qy = nullptr;
    TH1F *Q = nullptr;
};

MomentaDist getMomentaDistFromtxt(std::string txtFile, std::string histName)
{
    MomentaDist momdis;
    std::vector<double> Qi, Qt;
    std::ifstream f(txtFile.c_str());
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
    if (nBins < 2)
    {
        std::cerr << "[getMomentaDistFromtxt] Not enough points in " << txtFile << "\n";
        return momdis;
    }
    const double maxBin = Qi[nBins - 1] + (Qi[1] - Qi[0]) / 2.;
    momdis.Qt = new TH1F((std::string("Qt") + histName).c_str(), "Qt",
                         nBins, -maxBin, maxBin);
    for (int i = 0; i < nBins; i++)
        momdis.Qt->SetBinContent(i + 1, Qt[i]);

    return momdis;
}

// ---------------------------------------------------------------------------
// column == "y" -> dS/dQ_y : transverse component (px, or py by symmetry)
// column == "z" -> dS/dQ_z : longitudinal component (pz)
// column == "t" -> dS/dQ_t : transverse magnitude (pT); only defined/tabulated
//                  for Q_i > 0, so the output histogram covers [0, maxBin].
// ---------------------------------------------------------------------------
MomentaDist getMomentaDistFromtxtMultiCol(std::string txtFile, std::string histName,
                                          std::string column = "y")
{
    MomentaDist momdis;
    std::vector<double> Qi, Qsel;
    std::ifstream f(txtFile.c_str());
    if (!f.is_open())
    {
        std::cerr << "[getMomentaDistFromtxtMultiCol] Cannot open " << txtFile << "\n";
        return momdis;
    }
    std::string s;

    while (getline(f, s))
    {
        std::istringstream ss(s);
        double qi, qz, qt, qy, q;
        // A valid data line must parse as exactly 5 doubles.
        // Header / label / units / blank lines fail this and are skipped.
        if (!(ss >> qi >> qz >> qt >> qy >> q))
            continue;
        if (column == "t" && qi <= 0.0)
            continue; // dS/dQ_t (pT) is only tabulated for Q_i > 0

        Qi.push_back(qi);
        Qsel.push_back(column == "z" ? qz : (column == "t" ? qt : qy));
    }

    const int nBins = Qi.size();
    if (nBins < 2)
    {
        std::cerr << "[getMomentaDistFromtxtMultiCol] Not enough points in "
                  << txtFile << "\n";
        return momdis;
    }
    const double halfBin = (Qi[1] - Qi[0]) / 2.;
    const double xlo = (column == "t") ? 0.0 : -(Qi[nBins - 1] + halfBin);
    const double xhi = Qi[nBins - 1] + halfBin;
    momdis.Qt = new TH1F((std::string("Q") + column + histName).c_str(), "Q",
                         nBins, xlo, xhi);
    for (int i = 0; i < nBins; i++)
        momdis.Qt->SetBinContent(i + 1, Qsel[i]);

    return momdis;
}

// ---------------------------------------------------------------------------
//   format == 0  ->  two-column  (legacy: Qi Qt)             -> extracts Qt
//   format == 1  ->  multi-column (JT/CB 5-col)               -> extracts `column`
// ---------------------------------------------------------------------------
MomentaDist getMomentaDist(std::string txtFile, std::string histName, int format,
                           std::string column = "y")
{
    if (format == 1)
        return getMomentaDistFromtxtMultiCol(txtFile, histName, column);
    return getMomentaDistFromtxt(txtFile, histName); // default: legacy
}

// ---------------------------------------------------------------------------
// Read px_rf_rot or py_rf_rot from the experimental tree.
// ---------------------------------------------------------------------------
TH1F *getMomentaDistFromRoot(std::string rootFile, std::string branch,
                             bool useErelCut,
                             double erelMin, double erelMax,
                             int nBins, double maxBin,
                             std::string histName)
{
    auto *h = new TH1F(histName.c_str(), histName.c_str(), nBins, -maxBin, maxBin);
    h->SetDirectory(nullptr);

    auto *f = new TFile(rootFile.c_str());
    auto *tr = static_cast<TTree *>(f->Get("KinTree"));

    double p;
    tr->SetBranchAddress(branch.c_str(), &p);

    double erel = 0.0;
    if (useErelCut)
        tr->SetBranchAddress("Erel", &erel);

    double opa;
    tr->SetBranchAddress("califa_opa", &opa);

    const double opamin = 1.3;
    const double opamax = 1.62;

    // const double opamin = 0.;
    // const double opamax = 3;

    // const double opamin = 1.38;
    // const double opamax = 1.54;

    double sum = 0.0;
    Long64_t n = 0;
    const Long64_t nEntries = tr->GetEntries();

    for (Long64_t i = 0; i < nEntries; ++i)
    {
        tr->GetEntry(i);
        if (useErelCut && (erel * 1000.0 < erelMin || erel * 1000.0 > erelMax))
            continue;
        if (opa > opamin && opa < opamax)
        {
            sum += p * 1000.0;
            ++n;
        }
    }

    const double offset = (n > 0) ? sum / n : 0.0;
    if (n == 0)
        std::cerr << "[getMomentaDistFromRoot] WARNING: no events passed the cuts ("
                  << branch << ").\n";

    for (Long64_t i = 0; i < nEntries; ++i)
    {
        tr->GetEntry(i);
        if (useErelCut && (erel * 1000.0 < erelMin || erel * 1000.0 > erelMax))
            continue;
        if (opa > opamin && opa < opamax)
            h->Fill(p * 1000.0 - offset);
    }

    f->Close();
    delete f;
    return h;
}

// ---------------------------------------------------------------------------
// Rebin a theoretical histogram onto the binning of a reference histogram,
// preserving integrated content via overlap integration. Output normalized.
// ---------------------------------------------------------------------------
TH1F *buildTemplate(TH1F *hTheo, TH1F *hRef, std::string name)
{
    const int nBins = hRef->GetNbinsX();
    const double xmin = hRef->GetXaxis()->GetXmin();
    const double xmax = hRef->GetXaxis()->GetXmax();

    TH1F *hOut = new TH1F(name.c_str(), name.c_str(), nBins, xmin, xmax);
    hOut->SetDirectory(nullptr);

    for (int i = 1; i <= nBins; i++)
    {
        const double lo = hOut->GetXaxis()->GetBinLowEdge(i);
        const double hi = hOut->GetXaxis()->GetBinUpEdge(i);

        double integral = 0.0;
        const int nTheoBins = hTheo->GetNbinsX();
        for (int j = 1; j <= nTheoBins; j++)
        {
            const double tLo = hTheo->GetXaxis()->GetBinLowEdge(j);
            const double tHi = hTheo->GetXaxis()->GetBinUpEdge(j);
            const double overlap = std::max(0.0, std::min(hi, tHi) - std::max(lo, tLo));
            if (overlap <= 0)
                continue;
            const double tWidth = tHi - tLo;
            integral += hTheo->GetBinContent(j) * (overlap / tWidth);
        }
        hOut->SetBinContent(i, integral);
    }

    const double tot = hOut->Integral();
    if (tot > 0)
        hOut->Scale(1.0 / tot);
    return hOut;
}

// ---------------------------------------------------------------------------
// Wrap a histogram into a RooHistPdf (interpolation order 0).
// ---------------------------------------------------------------------------
RooHistPdf *create_pdf_from_histogram(RooRealVar &x, TH1F *h)
{
    std::string name = h->GetName();
    auto *dh = new RooDataHist((name + "_dh").c_str(), (name + "_dh").c_str(), x, h);
    auto *pdf = new RooHistPdf((name + "_pdf").c_str(), (name + "_pdf").c_str(),
                               x, *dh, 0);
    return pdf;
}

std::vector<double> make_fit(RooDataHist *data,
                             std::vector<RooHistPdf *> theory,
                             RooRealVar &x,
                             bool doPlot = false)
{
    const int nTerms = theory.size() - 1;

    std::vector<RooRealVar *> coefVars;
    RooArgList coeffs("coeffs");
    for (int i = 0; i < nTerms; i++)
    {
        std::string name = "val" + std::to_string(i);
        const double init = 1.0 / (nTerms + 1 - i);

        auto *v = new RooRealVar(name.c_str(), name.c_str(), init, 0., 1.);
        coefVars.push_back(v);
        coeffs.add(*v);
    }

    RooArgList pdfs("pdfs");
    for (int i = 0; i < nTerms + 1; i++)
        pdfs.add(*theory[i]);

    RooAddPdf model("model_fit", "model_fit", pdfs, coeffs, kTRUE);
    model.fitTo(*data, RooFit::Minos(kTRUE), RooFit::PrintLevel(-1));

    std::vector<double> falseCoefs;
    for (auto &coef : coeffs)
        falseCoefs.push_back(static_cast<RooRealVar *>(coef)->getVal());

    std::vector<double> trueCoefs(nTerms + 1);
    for (int i = 0; i < nTerms + 1; i++)
    {
        double ai = (i == nTerms) ? 1.0 : falseCoefs[i];
        for (int j = 0; j < i; j++)
            ai *= (1.0 - falseCoefs[j]);
        trueCoefs[i] = ai;
    }

    if (doPlot)
    {
        RooPlot *xframe = x.frame();
        data->plotOn(xframe);
        model.plotOn(xframe);
        xframe->Draw();
    }

    for (auto *v : coefVars)
        delete v;
    return trueCoefs;
}

TH1F *boostrapPoisson(TH1F *h)
{
    TH1F *h2 = static_cast<TH1F *>(h->Clone());
    h2->SetDirectory(nullptr);
    for (int i = 1; i <= h2->GetNbinsX(); i++)
        h2->SetBinContent(i, gRandom->Poisson(h2->GetBinContent(i)));
    return h2;
}

// ---------------------------------------------------------------------------
// 16% / 84% quantiles of a histogram.
// ---------------------------------------------------------------------------
std::vector<double> getQuantiles(TH1F *h)
{
    constexpr Int_t n = 2;
    Double_t p[n] = {0.16, 0.84};
    Double_t xp[n];
    h->GetQuantiles(n, xp, p);
    return {xp[0], xp[1]};
}

inline std::vector<double> getUncertainties(TH1F *h) { return getQuantiles(h); }

std::vector<double> calculateUncertainties(RooDataHist *data,
                                           std::vector<RooHistPdf *> theory,
                                           RooRealVar &x,
                                           std::vector<TH1F *> &parDistributions,
                                           int nToys = 10000,
                                           std::string tag = "boot")
{
    parDistributions.clear();
    for (int i = 0; i < (int)theory.size(); i++)
    {
        std::string name = tag + "_var_" + std::to_string(i);
        auto *h = new TH1F(name.c_str(), name.c_str(), 200, 0, 1);
        h->SetDirectory(nullptr);
        parDistributions.push_back(h);
    }

    auto *hData = static_cast<TH1F *>(data->createHistogram(
        ("TData_" + tag).c_str(), x));
    hData->SetDirectory(nullptr);

    for (int i = 0; i < nToys; i++)
    {
        auto *hToy = boostrapPoisson(hData);
        std::string dhName = "dh_" + tag + "_" + std::to_string(i);
        auto *dhToy = new RooDataHist(dhName.c_str(), dhName.c_str(),
                                      RooArgList(x), RooFit::Import(*hToy));
        std::vector<double> results = make_fit(dhToy, theory, x);
        for (int j = 0; j < (int)results.size(); j++)
            parDistributions[j]->Fill(results[j]);
        delete dhToy;
        delete hToy;
    }
    delete hData;

    std::vector<double> unc;
    for (int i = 0; i < (int)theory.size(); i++)
    {
        std::vector<double> res = getQuantiles(parDistributions[i]);
        unc.push_back(res[0]);
        unc.push_back(res[1]);
    }
    return unc;
}

std::vector<TH1F *> sampleFromTruth(RooAbsPdf &trueModel,
                                    std::vector<RooHistPdf *> theory,
                                    RooRealVar &x,
                                    int nEventsPerToy,
                                    int nToys = 10000)
{
    std::vector<TH1F *> dists;
    for (int i = 0; i < (int)theory.size(); i++)
    {
        std::string name = "truth_var_" + std::to_string(i);
        auto *h = new TH1F(name.c_str(), name.c_str(), 200, 0, 1);
        h->SetDirectory(nullptr);
        dists.push_back(h);
    }

    for (int i = 0; i < nToys; i++)
    {
        RooDataSet *toyData = trueModel.generate(x, nEventsPerToy);
        RooDataHist *toyBinned = toyData->binnedClone();
        std::vector<double> results = make_fit(toyBinned, theory, x);
        for (int j = 0; j < (int)results.size(); j++)
            dists[j]->Fill(results[j]);
        delete toyBinned;
        delete toyData;
    }
    return dists;
}

RooAddPdf *buildFixedModel(std::vector<RooHistPdf *> theory,
                           const std::vector<double> &phi,
                           std::vector<RooRealVar *> &storedCoefs,
                           std::string tag)
{
    const int N = (int)theory.size();
    // phi -> recursive a_k
    std::vector<double> aRec(N - 1);
    double remain = 1.0;
    for (int k = 0; k < N - 1; k++)
    {
        aRec[k] = (remain > 0) ? phi[k] / remain : 0.0;
        if (aRec[k] < 0.0)
            aRec[k] = 0.0;
        if (aRec[k] > 1.0)
            aRec[k] = 1.0;
        remain -= phi[k];
    }

    RooArgList coeffs(("coeffs_fixed_" + tag).c_str());
    storedCoefs.clear();
    for (int k = 0; k < N - 1; k++)
    {
        auto *v = new RooRealVar(Form("aFixed_%s_%d", tag.c_str(), k),
                                 Form("aFixed_%s_%d", tag.c_str(), k),
                                 aRec[k], 0., 1.);
        v->setConstant(true);
        storedCoefs.push_back(v);
        coeffs.add(*v);
    }

    RooArgList pdfs(("pdfs_fixed_" + tag).c_str());
    for (int k = 0; k < N; k++)
        pdfs.add(*theory[k]);

    return new RooAddPdf(("model_fixed_" + tag).c_str(),
                         ("model_fixed_" + tag).c_str(),
                         pdfs, coeffs, kTRUE);
}

void computeChi2(TH1F *h_data, RooAbsPdf *fitModel, RooRealVar &x,
                 int nC, std::string rangeTag,
                 double &chi2, int &ndf, int &nBinsUsed)
{
    const double Ntot = h_data->Integral();
    chi2 = 0.0;
    nBinsUsed = 0;

    for (int i = 1; i <= h_data->GetNbinsX(); i++)
    {
        const double lo = h_data->GetXaxis()->GetBinLowEdge(i);
        const double hi = h_data->GetXaxis()->GetBinUpEdge(i);
        std::string rname = "binChi2_" + rangeTag + "_" + std::to_string(i);
        x.setRange(rname.c_str(), lo, hi);
        RooAbsReal *binInt = fitModel->createIntegral(
            RooArgSet(x), RooFit::NormSet(RooArgSet(x)),
            RooFit::Range(rname.c_str()));
        const double Nexp = Ntot * binInt->getVal();
        delete binInt;

        if (Nexp <= 0)
            continue;
        const double Nobs = h_data->GetBinContent(i);
        const double d = Nobs - Nexp;
        chi2 += (d * d) / Nexp;
        nBinsUsed++;
    }
    ndf = nBinsUsed - nC;
}

void drawFitCanvas(const std::string &cName, const std::string &cTitle,
                   const std::string &axisTitle,
                   TH1F *h_data, RooDataHist *data,
                   std::vector<RooHistPdf *> theory,
                   RooRealVar &x,
                   const std::vector<double> &phi,
                   const std::vector<double> &lo16,
                   const std::vector<double> &hi84,
                   const std::vector<std::string> &labels,
                   const std::vector<int> &colors,
                   double chi2, int ndf,
                   std::string tag)
{
    using namespace RooFit;
    const int nC = (int)theory.size();

    std::vector<RooRealVar *> fixedCoefs;
    RooAddPdf *modelFixed = buildFixedModel(theory, phi, fixedCoefs, tag);

    TCanvas *c = new TCanvas(cName.c_str(), cTitle.c_str(), 800, 800);

    TPad *padTop = new TPad(("padTop_" + tag).c_str(), "padTop", 0.0, 0.30, 1.0, 1.0);
    padTop->SetBottomMargin(0.02);
    padTop->Draw();
    TPad *padBot = new TPad(("padBot_" + tag).c_str(), "padBot", 0.0, 0.0, 1.0, 0.30);
    padBot->SetTopMargin(0.02);
    padBot->SetBottomMargin(0.32);
    padBot->Draw();

    padTop->cd();
    RooPlot *xframe = x.frame(Title(cTitle.c_str()));
    data->plotOn(xframe, Name(("data_" + tag).c_str()));
    modelFixed->plotOn(xframe, LineColor(kRed), LineWidth(2),
                       Name(("total_" + tag).c_str()));

    if (nC > 1)
        for (int k = 0; k < nC; k++)
            modelFixed->plotOn(xframe,
                               Components(*theory[k]),
                               LineStyle(kDashed),
                               LineColor(colors[k % colors.size()]),
                               LineWidth(2),
                               Name(Form("comp_%s_%d", tag.c_str(), k)));

    xframe->GetXaxis()->SetLabelSize(0.0);
    xframe->GetXaxis()->SetTitleSize(0.0);
    xframe->Draw();

    auto *leg = new TLegend(0.55, 0.55, 0.89, 0.89);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->AddEntry(xframe->findObject(("data_" + tag).c_str()), "Exp", "ep");
    leg->AddEntry(xframe->findObject(("total_" + tag).c_str()), "Fit total", "l");
    for (int k = 0; k < nC; k++)
    {
        const char *entryName = (nC > 1)
                                    ? Form("comp_%s_%d", tag.c_str(), k)
                                    : ("total_" + tag).c_str();
        const double minus = phi[k] - lo16[k];
        const double plus = hi84[k] - phi[k];
        leg->AddEntry(xframe->findObject(entryName),
                      Form("%s (%.1f_{-%.1f}^{+%.1f})%%",
                           labels[k].c_str(),
                           phi[k] * 100, minus * 100, plus * 100),
                      "l");
    }
    leg->AddEntry((TObject *)nullptr,
                  Form("#chi^{2}/ndf = %.2f/%d = %.2f",
                       chi2, ndf, (ndf > 0 ? chi2 / ndf : -1.0)),
                  "");
    leg->Draw();

    // Residual / pull plot
    padBot->cd();
    TH1F *h_res = new TH1F(("h_res_" + tag).c_str(),
                           Form(";%s;(N_{obs} - N_{exp})/#sqrt{N_{exp}}",
                                axisTitle.c_str()),
                           h_data->GetNbinsX(),
                           h_data->GetXaxis()->GetXmin(),
                           h_data->GetXaxis()->GetXmax());
    h_res->SetDirectory(nullptr);

    const double Ntot = h_data->Integral();
    double maxAbsR = 0.0;
    for (int i = 1; i <= h_data->GetNbinsX(); i++)
    {
        const double lo = h_data->GetXaxis()->GetBinLowEdge(i);
        const double hi = h_data->GetXaxis()->GetBinUpEdge(i);
        std::string rname = "binRes_" + tag + "_" + std::to_string(i);
        x.setRange(rname.c_str(), lo, hi);
        RooAbsReal *binInt = modelFixed->createIntegral(
            RooArgSet(x), NormSet(RooArgSet(x)),
            Range(rname.c_str()));
        const double Nexp = Ntot * binInt->getVal();
        delete binInt;

        if (Nexp <= 0)
            continue;
        const double Nobs = h_data->GetBinContent(i);
        const double r = (Nobs - Nexp) / std::sqrt(Nexp);
        h_res->SetBinContent(i, r);
        h_res->SetBinError(i, 1.0);
        if (std::abs(r) > maxAbsR)
            maxAbsR = std::abs(r);
    }

    const double yRange = std::max(3.5, std::ceil(maxAbsR) + 0.5);
    h_res->GetYaxis()->SetRangeUser(-yRange, yRange);
    h_res->GetYaxis()->SetNdivisions(505);
    h_res->GetYaxis()->SetTitleSize(0.10);
    h_res->GetYaxis()->SetTitleOffset(0.45);
    h_res->GetYaxis()->SetLabelSize(0.09);
    h_res->GetXaxis()->SetTitleSize(0.12);
    h_res->GetXaxis()->SetTitleOffset(1.10);
    h_res->GetXaxis()->SetLabelSize(0.10);
    h_res->SetMarkerStyle(20);
    h_res->SetMarkerSize(0.7);
    h_res->SetLineColor(kBlack);
    h_res->Draw("E1");

    TLine *l0 = new TLine(h_data->GetXaxis()->GetXmin(), 0,
                          h_data->GetXaxis()->GetXmax(), 0);
    l0->SetLineColor(kRed);
    l0->SetLineWidth(2);
    l0->Draw();
    for (int s : {-1, 1})
    {
        TLine *ls = new TLine(h_data->GetXaxis()->GetXmin(), s,
                              h_data->GetXaxis()->GetXmax(), s);
        ls->SetLineColor(kGray + 1);
        ls->SetLineStyle(2);
        ls->Draw();
    }

    c->cd();
    c->Update();
}

void drawBootstrapCanvas(const std::string &cName, const std::string &cTitle,
                         std::vector<TH1F *> bootDists,
                         const std::vector<double> &phi,
                         const std::vector<double> &lo16,
                         const std::vector<double> &hi84,
                         const std::vector<std::string> &labels)
{
    const int nPars = (int)bootDists.size();
    auto *c = new TCanvas(cName.c_str(), cTitle.c_str(), 1400, 450);
    c->Divide(nPars, 1);

    for (int i = 0; i < nPars; i++)
    {
        c->cd(i + 1);
        TH1F *hb = bootDists[i];
        if (hb->Integral() > 0)
            hb->Scale(1.0 / hb->Integral());

        hb->SetLineColor(kBlue + 1);
        hb->SetLineWidth(2);
        std::string title = labels[i] + " (#phi_{" + std::to_string(i) +
                            "});fitted value;normalized";
        hb->SetTitle(title.c_str());

        double ymax = hb->GetMaximum() * 1.25;
        hb->SetMaximum(ymax);
        hb->SetMinimum(0);
        hb->Draw("HIST");

        auto *lNom = new TLine(phi[i], 0, phi[i], ymax);
        lNom->SetLineColor(kBlack);
        lNom->SetLineStyle(2);
        lNom->SetLineWidth(2);
        lNom->Draw();

        auto *lLo = new TLine(lo16[i], 0, lo16[i], ymax);
        lLo->SetLineColor(kGreen + 2);
        lLo->SetLineStyle(3);
        lLo->SetLineWidth(2);
        lLo->Draw();

        auto *lHi = new TLine(hi84[i], 0, hi84[i], ymax);
        lHi->SetLineColor(kGreen + 2);
        lHi->SetLineStyle(3);
        lHi->SetLineWidth(2);
        lHi->Draw();

        auto *leg = new TLegend(0.55, 0.68, 0.88, 0.88);
        leg->SetBorderSize(0);
        leg->SetFillStyle(0);
        leg->AddEntry(hb, "Bootstrap", "l");
        leg->AddEntry(lNom, "Nominal fit", "l");
        leg->AddEntry(lLo, "16% / 84% q.", "l");
        leg->Draw();
    }
    c->Update();
}

void test(int nEvents = 6000, int nToys = 10000, double noiseSigma = 5.0)
{
    RooMsgService::instance().setGlobalKillBelow(RooFit::ERROR);

    // -----------------------------------------------------------------------
    // 1) Bertulani (CB) calculation: same three states, read as px (dS/dQ_y)
    //    and pz (dS/dQ_z) for the generator, and as pT (dS/dQ_t) for the fit
    //    templates. Using the same basis for generation and fit makes this a
    //    closure test of the pT fitting procedure.
    // -----------------------------------------------------------------------
    std::vector<std::string> inFiles = {
        "/nucl_lustre/pablogrusell/g249/g249_analysis/theory/CB/momdis_1p12_19_215.txt",
        "/nucl_lustre/pablogrusell/g249/g249_analysis/theory/CB/momdis_2s12_14_455.txt",
        "/nucl_lustre/pablogrusell/g249/g249_analysis/theory/CB/momdis_1d52_14_455.txt"};
    std::vector<std::string> labels = {"1p_{1/2} [CB]", "2s_{1/2} [CB]", "1d_{5/2} [CB]"};
    const int nC = (int)inFiles.size();

    const double trueF_p12 = 0.1;                               // p_{1/2}
    const double trueF_d52 = 0.9;                               // d_{5/2}
    std::vector<double> trueVals = {trueF_p12, 0.0, trueF_d52}; // p12, s12, d52

    const int nBinsRef = 40;
    const double maxBin = 300.0;

    RooRealVar px("px", "p_{x} [MeV/c]", -maxBin, maxBin);
    RooRealVar pz("pz", "p_{z} [MeV/c]", -maxBin, maxBin);
    RooRealVar pT("pT", "p_{T} [MeV/c]", 0.0, maxBin);

    auto *hRefPx = new TH1F("hRefPx", "binning reference px", nBinsRef, -maxBin, maxBin);
    hRefPx->SetDirectory(nullptr);
    auto *hRefPz = new TH1F("hRefPz", "binning reference pz", nBinsRef, -maxBin, maxBin);
    hRefPz->SetDirectory(nullptr);
    auto *hRefPt = new TH1F("hRefPt", "binning reference pT", nBinsRef / 2, 0.0, maxBin);
    hRefPt->SetDirectory(nullptr);

    std::vector<TH1F *> tmplPx(nC), tmplPz(nC), tmplPt(nC);
    std::vector<RooHistPdf *> theoryPx(nC), theoryPz(nC), theoryPt(nC);
    for (int k = 0; k < nC; k++)
    {
        auto mdPx = getMomentaDist(inFiles[k], Form("px_%d", k), 1, "y");
        auto mdPz = getMomentaDist(inFiles[k], Form("pz_%d", k), 1, "z");
        auto mdPt = getMomentaDist(inFiles[k], Form("pt_%d", k), 1, "t");

        tmplPx[k] = buildTemplate(mdPx.Qt, hRefPx, Form("tmpl_px_%d", k));
        tmplPz[k] = buildTemplate(mdPz.Qt, hRefPz, Form("tmpl_pz_%d", k));
        tmplPt[k] = buildTemplate(mdPt.Qt, hRefPt, Form("tmpl_pt_%d", k));

        theoryPx[k] = create_pdf_from_histogram(px, tmplPx[k]);
        theoryPz[k] = create_pdf_from_histogram(pz, tmplPz[k]);
        theoryPt[k] = create_pdf_from_histogram(pT, tmplPt[k]);
    }

    // -----------------------------------------------------------------------
    // Generator model: 0.1 p12 + 0.9 d52, in the px/pz (CB) basis.
    // -----------------------------------------------------------------------
    RooRealVar frac_p12("frac_p12", "f(p12)", trueF_p12, 0., 1.);
    RooRealVar frac_d52("frac_d52", "f(d52)", trueF_d52, 0., 1.);

    RooAddPdf true_model_px("true_model_px", "0.1 p12 + 0.9 d52 [CB, px basis]",
                            RooArgList(*theoryPx[0], *theoryPx[2]),
                            RooArgList(frac_p12, frac_d52));
    RooAddPdf true_model_pz("true_model_pz", "0.1 p12 + 0.9 d52 [CB, pz basis]",
                            RooArgList(*theoryPz[0], *theoryPz[2]),
                            RooArgList(frac_p12, frac_d52));

    // -----------------------------------------------------------------------
    // 2) Generate toy px and pz from the Bertulani model. py is obtained by
    //    conservation: azimuthal symmetry around the beam (z) axis means px
    //    and py follow the same 1D distribution, so py is drawn independently
    //    from the same px model rather than from a separate template.
    // -----------------------------------------------------------------------
    RooDataSet *dataPx = true_model_px.generate(px, nEvents, RooFit::AutoBinned(kFALSE));
    RooDataSet *dataPy = true_model_px.generate(px, nEvents, RooFit::AutoBinned(kFALSE)); // independent draw -> py
    RooDataSet *dataPz = true_model_pz.generate(pz, nEvents, RooFit::AutoBinned(kFALSE));

    // -----------------------------------------------------------------------
    // 3) Add Gaussian resolution noise to px, py; 4) compute pT.
    // -----------------------------------------------------------------------
    auto *hPt = new TH1F("hPt_toy", "toy p_{T}", nBinsRef / 2, 0.0, maxBin);
    hPt->SetDirectory(nullptr);

    for (int i = 0; i < nEvents; i++)
    {
        const double pxVal = dataPx->get(i)->getRealValue("px") + gRandom->Gaus(0.0, noiseSigma);
        const double pyVal = dataPy->get(i)->getRealValue("px") + gRandom->Gaus(0.0, noiseSigma);
        const double ptVal = std::sqrt(pxVal * pxVal + pyVal * pyVal) + 10;
        if (ptVal < maxBin)
            hPt->Fill(ptVal);
    }

    // -----------------------------------------------------------------------
    // 5) Bertulani pT templates (dS/dQ_t) + 6) fit.
    // -----------------------------------------------------------------------
    RooDataHist dataPt("dataPt_toy", "toy data pT", RooArgList(pT), hPt);

    std::vector<double> params = make_fit(&dataPt, theoryPt, pT, /*doPlot=*/true);

    std::vector<TH1F *> bootDists;
    std::vector<double> uncs = calculateUncertainties(&dataPt, theoryPt, pT,
                                                      bootDists, nToys, "boot_test");

    std::cout << "\n=== test(): px,py generated from Bertulani (CB) px basis "
                 "(0.1 p12 + 0.9 d52), noise sigma = "
              << noiseSigma
              << " MeV/c, fit pT = sqrt(px^2+py^2) with Bertulani pT "
                 "(dS/dQ_t) templates (p12+s12+d52) ===\n";
    for (int i = 0; i < (int)params.size(); i++)
        std::cout << labels[i] << " (phi_" << i << "): "
                  << params[i]
                  << "  -" << params[i] - uncs[2 * i]
                  << "  +" << uncs[2 * i + 1] - params[i]
                  << "   (true ~ " << trueVals[i] << ")\n";

    std::vector<int> colors = {kBlue, kMagenta, kGreen + 2};
    std::vector<double> lo16(nC), hi84(nC);
    for (int k = 0; k < nC; k++)
    {
        lo16[k] = uncs[2 * k];
        hi84[k] = uncs[2 * k + 1];
    }

    std::vector<RooRealVar *> dummyCoefs;
    std::unique_ptr<RooAddPdf> modelFixed(
        buildFixedModel(theoryPt, params, dummyCoefs, "chi2_test"));

    double chi2;
    int ndf, nbu;
    computeChi2(hPt, modelFixed.get(), pT, nC, "test", chi2, ndf, nbu);

    drawFitCanvas("c_fit_test_pT",
                  "Toy p_{T} (from noisy px, py) fitted with Bertulani p_{T} templates",
                  "p_{T} [MeV/c]",
                  hPt, &dataPt, theoryPt, pT,
                  params, lo16, hi84,
                  labels, colors, chi2, ndf, "test");

    drawBootstrapCanvas("c_boot_test_pT", "Bootstrap distributions p_{T}",
                        bootDists, params, lo16, hi84, labels);

    std::cout << "\nchi2/ndf = " << chi2 << " / " << ndf
              << " = " << (ndf > 0 ? chi2 / ndf : -1.) << "\n";

    delete dataPx;
    delete dataPy;
    delete dataPz;
}

void fitMomdis(double erelMin = 1.5, double erelMax = 4,
               int nToys = 1000, int theoFormat = 0)
{
    using namespace RooFit;
    RooMsgService::instance().setGlobalKillBelow(RooFit::ERROR);

    // ---------------------------------------------------------------------
    // Inputs
    // ---------------------------------------------------------------------
    std::vector<std::string> inFilesTheo = {
        "/nucl_lustre/pablogrusell/g249/g249_analysis/theory/JT/25F/sigt_d52-4.txt",
        "/nucl_lustre/pablogrusell/g249/g249_analysis/theory/JT/25F/sigt_s12-4.txt",
        "/nucl_lustre/pablogrusell/g249/g249_analysis/theory/JT/25F/sigt_p12-4.txt"};

    // std::vector<std::string> inFilesTheo = {
    //     "/nucl_lustre/pablogrusell/g249/g249_analysis/theory/CB/momdis_1d52_14_455.txt",
    //     "/nucl_lustre/pablogrusell/g249/g249_analysis/theory/CB/momdis_2s12_14_455.txt"};

    std::vector<std::string> labels = {"1d_{5/2}", "2s_{1/2}", "1p_{1/2}"};
    std::vector<int> colors = {kBlue, kGreen + 2, kMagenta, kCyan + 1};

    // std::string inFileExp =
    //     "/nucl_lustre/pablogrusell/g249/g249_analysis/results/final/24O_analyzed_test.root";

    std::string inFileExp =
        "/nucl_lustre/pablogrusell/g249/g249_analysis/results/final/23O_analyzed_test.root";

    const int nC = (int)inFilesTheo.size();
    const int nBinsRef = 5;
    const double maxBin = 100.0;
    const bool useErelCut = !(erelMin < 0 && erelMax < 0);

    // ---------------------------------------------------------------------
    // Experimental histograms for px and py
    // ---------------------------------------------------------------------
    TH1F *h_px = getMomentaDistFromRoot(inFileExp, "px_rf_rot", useErelCut,
                                        erelMin, erelMax, nBinsRef, maxBin,
                                        "h_exp_px");
    TH1F *h_py = getMomentaDistFromRoot(inFileExp, "py_rf_rot", useErelCut,
                                        erelMin, erelMax, nBinsRef, maxBin,
                                        "h_exp_py");

    // ---------------------------------------------------------------------
    // Templates on the experimental binning.
    // For theoFormat==1 the reader returns the dS/dQ_y column (in .Qt),
    // which is used as the template for BOTH the px and py fits.
    // ---------------------------------------------------------------------
    RooRealVar x("x", "p [MeV/c]", -maxBin, maxBin);

    std::vector<TH1F *> templates(nC);
    std::vector<RooHistPdf *> theory(nC);
    for (int k = 0; k < nC; k++)
    {
        auto md = getMomentaDist(inFilesTheo[k], Form("theo_%d", k), theoFormat);
        templates[k] = buildTemplate(md.Qt, h_px, Form("tmpl_%d", k));
        theory[k] = create_pdf_from_histogram(x, templates[k]);
    }

    // =====================================================================
    // PX FIT
    // =====================================================================
    std::cout << "\n############### FIT px ###############\n";
    RooDataHist data_px("data_px", "experimental data px",
                        RooArgList(x), h_px);
    std::vector<double> phi_px = make_fit(&data_px, theory, x);

    std::vector<TH1F *> boot_px;
    std::vector<double> q_px = calculateUncertainties(&data_px, theory, x,
                                                      boot_px, nToys, "boot_px");
    std::vector<double> lo16_px(nC), hi84_px(nC);
    for (int k = 0; k < nC; k++)
    {
        lo16_px[k] = q_px[2 * k];
        hi84_px[k] = q_px[2 * k + 1];
    }

    std::vector<RooRealVar *> dummy_px;
    std::unique_ptr<RooAddPdf> model_px(buildFixedModel(theory, phi_px, dummy_px, "chi2_px"));
    double chi2_px;
    int ndf_px, nbu_px;
    computeChi2(h_px, model_px.get(), x, nC, "px", chi2_px, ndf_px, nbu_px);

    std::cout << "N(px)    = " << h_px->Integral() << "\n";
    std::cout << "chi2/ndf = " << chi2_px << " / " << ndf_px
              << " = " << (ndf_px > 0 ? chi2_px / ndf_px : -1.) << "\n";
    for (int k = 0; k < nC; k++)
        std::cout << "phi(" << labels[k] << ") = " << phi_px[k]
                  << "  - " << phi_px[k] - lo16_px[k]
                  << "  + " << hi84_px[k] - phi_px[k] << "\n";

    // =====================================================================
    // PY FIT
    // =====================================================================
    std::cout << "\n############### FIT py ###############\n";
    RooDataHist data_py("data_py", "experimental data py",
                        RooArgList(x), h_py);
    std::vector<double> phi_py = make_fit(&data_py, theory, x);

    std::vector<TH1F *> boot_py;
    std::vector<double> q_py = calculateUncertainties(&data_py, theory, x,
                                                      boot_py, nToys, "boot_py");
    std::vector<double> lo16_py(nC), hi84_py(nC);
    for (int k = 0; k < nC; k++)
    {
        lo16_py[k] = q_py[2 * k];
        hi84_py[k] = q_py[2 * k + 1];
    }

    std::vector<RooRealVar *> dummy_py;
    std::unique_ptr<RooAddPdf> model_py(buildFixedModel(theory, phi_py, dummy_py, "chi2_py"));
    double chi2_py;
    int ndf_py, nbu_py;
    computeChi2(h_py, model_py.get(), x, nC, "py", chi2_py, ndf_py, nbu_py);

    std::cout << "N(py)    = " << h_py->Integral() << "\n";
    std::cout << "chi2/ndf = " << chi2_py << " / " << ndf_py
              << " = " << (ndf_py > 0 ? chi2_py / ndf_py : -1.) << "\n";
    for (int k = 0; k < nC; k++)
        std::cout << "phi(" << labels[k] << ") = " << phi_py[k]
                  << "  - " << phi_py[k] - lo16_py[k]
                  << "  + " << hi84_py[k] - phi_py[k] << "\n";

    // =====================================================================
    // Plots
    // =====================================================================
    drawFitCanvas("c_px", "Template fit p_{x}", "p_{x} [MeV/c]",
                  h_px, &data_px, theory, x,
                  phi_px, lo16_px, hi84_px,
                  labels, colors, chi2_px, ndf_px, "px");

    drawFitCanvas("c_py", "Template fit p_{y}", "p_{y} [MeV/c]",
                  h_py, &data_py, theory, x,
                  phi_py, lo16_py, hi84_py,
                  labels, colors, chi2_py, ndf_py, "py");

    drawBootstrapCanvas("c_boot_px", "Bootstrap distributions p_{x}",
                        boot_px, phi_px, lo16_px, hi84_px, labels);
    drawBootstrapCanvas("c_boot_py", "Bootstrap distributions p_{y}",
                        boot_py, phi_py, lo16_py, hi84_py, labels);
}