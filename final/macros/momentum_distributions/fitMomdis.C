// ===========================================================================
// COMMON HELPERS (used by both test() and fitMomdisBootstrap())
// ===========================================================================

// ---------------------------------------------------------------------------
// Container for theoretical momentum distributions.
// ---------------------------------------------------------------------------
struct MomentaDist
{
    TH1F *Qz = nullptr;
    TH1F *Qt = nullptr;
    TH1F *Qy = nullptr;
    TH1F *Q = nullptr;
};

// ---------------------------------------------------------------------------
// Read a theoretical momentum distribution from a two-column text file.
// (LEGACY format: Qi  Qt). Extracts the transverse component into .Qt
// ---------------------------------------------------------------------------
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
// Read a theoretical momentum distribution from a MULTI-COLUMN text file
// (the "JT 5-column" format):
//
//   <total cross section>
//   <blank line(s)>
//   Momentum distributions
//     Q_i      dS/dQ_z   dS/dQ_t   dS/dQ_y   dS/dQ
//   (MeV/C)  (mb/MeV/c) (mb/MeV/c) (mb/MeV/c) (mb/MeV/c)
//   -986.635   0.00000  0.95E-05  0.17E-06  0.24E-05
//   ...
//
// Column layout:  Q_i  dS/dQ_z  dS/dQ_t  dS/dQ_y  dS/dQ
//
// For this format the fit (both px and py) is performed against the
// dS/dQ_y column (4th column). We store that distribution in the .Qt
// slot of MomentaDist, because that is the field the downstream code
// (buildTemplate / create_pdf_from_histogram) reads from. The slot name
// is historical; the *content* here is dS/dQ_y.
// ---------------------------------------------------------------------------
MomentaDist getMomentaDistFromtxtMultiCol(std::string txtFile, std::string histName)
{
    MomentaDist momdis;
    std::vector<double> Qi, Qy;
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
        Qi.push_back(qi);
        Qy.push_back(qy); // dS/dQ_y -> used for BOTH px and py fits
    }

    const int nBins = Qi.size();
    if (nBins < 2)
    {
        std::cerr << "[getMomentaDistFromtxtMultiCol] Not enough points in "
                  << txtFile << "\n";
        return momdis;
    }
    const double maxBin = Qi[nBins - 1] + (Qi[1] - Qi[0]) / 2.;
    // Stored in .Qt so downstream code is unchanged, but content is dS/dQ_y.
    momdis.Qt = new TH1F((std::string("Qy") + histName).c_str(), "Qy",
                         nBins, -maxBin, maxBin);
    for (int i = 0; i < nBins; i++)
        momdis.Qt->SetBinContent(i + 1, Qy[i]);

    return momdis;
}

// ---------------------------------------------------------------------------
// Dispatcher: pick the reader based on a format flag.
//   format == 0  ->  two-column  (legacy: Qi Qt)             -> extracts Qt
//   format == 1  ->  multi-column (JT 5-col)                 -> extracts Qy
// In both cases the chosen distribution is returned in momdis.Qt.
// ---------------------------------------------------------------------------
MomentaDist getMomentaDist(std::string txtFile, std::string histName, int format)
{
    if (format == 1)
        return getMomentaDistFromtxtMultiCol(txtFile, histName);
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

// ---------------------------------------------------------------------------
// Fit a binned dataset to a sum of templates with the recursive-fraction
// RooAddPdf convention. Returns physical fractions (phi_0, ..., phi_{N-1}).
// ---------------------------------------------------------------------------
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

// ---------------------------------------------------------------------------
// Poissonian bootstrap of an observed histogram.
// ---------------------------------------------------------------------------
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

// Keep the original API name (used by test()) as an alias.
inline std::vector<double> getUncertainties(TH1F *h) { return getQuantiles(h); }

// ---------------------------------------------------------------------------
// Bootstrap: Poisson fluctuations on the observed histogram, refit each toy.
// Fills parDistributions and returns interleaved {q16,q84} per parameter.
// 'tag' makes RooFit object names unique across calls.
// ---------------------------------------------------------------------------
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

// ---------------------------------------------------------------------------
// Truth resampling: independent datasets from the true model, refit each.
// Returns one fitted-value distribution per parameter.
// ---------------------------------------------------------------------------
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

// ---------------------------------------------------------------------------
// Build an RooAddPdf with FIXED recursive coefficients matching the given
// physical fractions. Used to plot the best-fit model and compute chi2.
// ---------------------------------------------------------------------------
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

// ---------------------------------------------------------------------------
// Analytic chi2: sum_i (N_obs - N_exp)^2 / N_exp ; ndf = nBinsUsed - nC.
// ---------------------------------------------------------------------------
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

// ---------------------------------------------------------------------------
// Fit + residual canvas for one momentum component (used by real-data fit).
// ---------------------------------------------------------------------------
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

// ---------------------------------------------------------------------------
// Bootstrap distributions on one canvas (real-data version: bootstrap only).
// ---------------------------------------------------------------------------
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

// ===========================================================================
// 1) test(): toy validation
//
// Generates pseudo-data from a known mixture of the three theoretical
// templates (d_{5/2}, s_{1/2}, p_{1/2}) and compares the bootstrap-based
// uncertainty with the truth-resampling distribution. Useful to verify
// that the bootstrap correctly reproduces the spread expected from
// independent realizations of the true model.
//
// 'theoFormat' selects the theory-file reader (0 = legacy 2-col / Qt,
//              1 = JT multi-col / Qy).
// ===========================================================================
void test(int theoFormat_gen = 0, // 1 = JT multi-col  (generador)
          int theoFormat_fit = 1) // 0 = legacy 2-col  (ajuste CB)
{
    RooMsgService::instance().setGlobalKillBelow(RooFit::ERROR);

    // -----------------------------------------------------------------------
    // 1) Templates GENERADORES: Tostevin (JT, formato multi-columna)
    //    Mezcla verdadera: 0.7 * p_{1/2} + 0.3 * d_{5/2}
    // -----------------------------------------------------------------------
    std::vector<std::string> inFilesGen = {
        "/nucl_lustre/pablogrusell/g249/g249_analysis/theory/JT/25F/sigt_p12-gs.txt",
        "/nucl_lustre/pablogrusell/g249/g249_analysis/theory/JT/25F/sigt_d52-gs.txt"};
    std::vector<std::string> labelsGen = {"1p_{1/2} [JT]", "1d_{5/2} [JT]"};

    const double trueF_p12 = 0.1; // p_{1/2}
    const double trueF_d52 = 0.9; // d_{5/2}

    // -----------------------------------------------------------------------
    // 2) Templates de AJUSTE: Bertulani (CB, formato legacy 2-col)
    //    p_{1/2} + s_{1/2} + d_{5/2}
    // -----------------------------------------------------------------------
    std::vector<std::string> inFilesFit = {
        "/nucl_lustre/pablogrusell/g249/g249_analysis/theory/CB/momdis_1p12_19_215.txt",
        "/nucl_lustre/pablogrusell/g249/g249_analysis/theory/CB/momdis_2s12_14_455.txt",
        "/nucl_lustre/pablogrusell/g249/g249_analysis/theory/CB/momdis_1d52_14_455.txt"};
    std::vector<std::string> labelsFit = {"1p_{1/2} [CB]", "2s_{1/2} [CB]", "1d_{5/2} [CB]"};

    // -----------------------------------------------------------------------
    // Binning de referencia (igual que en fitMomdis)
    // -----------------------------------------------------------------------
    const int nBinsRef = 40;
    const double maxBin = 300.0;
    RooRealVar x("x", "p_{x} [MeV/c]", -maxBin, maxBin);

    auto *hRef = new TH1F("hRef", "binning reference", nBinsRef, -maxBin, maxBin);
    hRef->SetDirectory(nullptr);

    // -----------------------------------------------------------------------
    // Cargar templates generadores (JT) y construir el modelo verdadero
    // -----------------------------------------------------------------------
    const int nGen = (int)inFilesGen.size();
    std::vector<TH1F *> tmplGen(nGen);
    std::vector<RooHistPdf *> theoryGen(nGen);
    for (int k = 0; k < nGen; k++)
    {
        auto md = getMomentaDist(inFilesGen[k], Form("gen_%d", k), theoFormat_gen);
        tmplGen[k] = buildTemplate(md.Qt, hRef, Form("tmpl_gen_%d", k));
        theoryGen[k] = create_pdf_from_histogram(x, tmplGen[k]);
    }

    // Modelo generador: RooAddPdf NO recursivo -> fracciones físicas directas
    // Índice 0 = p_{1/2}, índice 1 = d_{5/2}
    RooRealVar frac_p12("frac_p12", "f(p12)", trueF_p12, 0., 1.);
    RooRealVar frac_d52("frac_d52", "f(d52)", trueF_d52, 0., 1.);

    RooAddPdf true_model("true_model", "0.7 p12 + 0.3 d52 [JT]",
                         RooArgList(*theoryGen[0], *theoryGen[1]),
                         RooArgList(frac_p12, frac_d52));

    const int nEvents = 6000;

    // -----------------------------------------------------------------------
    // Cargar templates de ajuste (CB)
    // -----------------------------------------------------------------------
    const int nFit = (int)inFilesFit.size();
    std::vector<TH1F *> tmplFit(nFit);
    std::vector<RooHistPdf *> theoryFit(nFit);
    for (int k = 0; k < nFit; k++)
    {
        auto md = getMomentaDist(inFilesFit[k], Form("fit_%d", k), theoFormat_fit);
        tmplFit[k] = buildTemplate(md.Qt, hRef, Form("tmpl_fit_%d", k));
        theoryFit[k] = create_pdf_from_histogram(x, tmplFit[k]);
    }

    // -----------------------------------------------------------------------
    // Generar dataset observado y ajustar con CB
    // -----------------------------------------------------------------------
    RooDataSet *data = true_model.generate(x, nEvents);
    RooDataHist *binnedData = data->binnedClone();

    // Ajuste nominal
    std::vector<double> params = make_fit(binnedData, theoryFit, x, /*doPlot=*/true);

    // Bootstrap (incertidumbres asimétricas)
    std::vector<TH1F *> bootDists;
    std::vector<double> uncs = calculateUncertainties(binnedData, theoryFit, x,
                                                      bootDists, 10000, "boot_test");

    // Truth resampling (con el modelo generador JT, ajustado con CB)
    std::vector<TH1F *> truthDists = sampleFromTruth(true_model, theoryFit, x,
                                                     nEvents, 10000);

    // -----------------------------------------------------------------------
    // Resumen por consola
    // -----------------------------------------------------------------------
    // Valores "verdaderos" en el espacio CB: la señal p12 debería
    // absorberse en p12-CB, la de d52 en d52-CB, y s12 ~ 0.
    std::vector<double> trueVals = {trueF_p12, 0.0, trueF_d52}; // p12, s12, d52

    std::cout << "\n=== test() JT->CB: 0.7 p12 + 0.3 d52 generado con JT, "
                 "ajustado con CB (p12+s12+d52) ===\n";
    for (int i = 0; i < (int)params.size(); i++)
        std::cout << labelsFit[i] << " (phi_" << i << "): "
                  << params[i]
                  << "  -" << params[i] - uncs[2 * i]
                  << "  +" << uncs[2 * i + 1] - params[i]
                  << "   (true ~ " << trueVals[i] << ")\n";

    // -----------------------------------------------------------------------
    // Canvas de comparación: bootstrap vs truth resampling
    // -----------------------------------------------------------------------
    const int nPars = (int)params.size();
    auto *cmp = new TCanvas("cmp_JT_CB",
                            "JT generated, CB fitted: bootstrap vs truth",
                            1400, 450);
    cmp->Divide(nPars, 1);

    for (int i = 0; i < nPars; i++)
    {
        cmp->cd(i + 1);
        TH1F *hb = bootDists[i];
        TH1F *ht = truthDists[i];

        if (hb->Integral() > 0)
            hb->Scale(1.0 / hb->Integral());
        if (ht->Integral() > 0)
            ht->Scale(1.0 / ht->Integral());

        hb->SetLineColor(kBlue + 1);
        hb->SetLineWidth(2);
        ht->SetLineColor(kRed + 1);
        ht->SetLineWidth(2);

        std::string title = labelsFit[i] + " (#phi_{" + std::to_string(i) + "}" + ");fitted value;normalized";
        hb->SetTitle(title.c_str());

        double ymax = std::max(hb->GetMaximum(), ht->GetMaximum()) * 1.25;
        hb->SetMaximum(ymax);
        hb->SetMinimum(0);
        hb->Draw("HIST");
        ht->Draw("HIST SAME");

        // Valor nominal del ajuste
        auto *lNom = new TLine(params[i], 0, params[i], ymax);
        lNom->SetLineColor(kBlack);
        lNom->SetLineStyle(2);
        lNom->SetLineWidth(2);
        lNom->Draw();

        // Valor verdadero esperado
        auto *lTrue = new TLine(trueVals[i], 0, trueVals[i], ymax);
        lTrue->SetLineColor(kGreen + 2);
        lTrue->SetLineStyle(1);
        lTrue->SetLineWidth(2);
        lTrue->Draw();

        // Cuantiles 16 / 84 del bootstrap
        auto *lLo = new TLine(uncs[2 * i], 0, uncs[2 * i], ymax);
        auto *lHi = new TLine(uncs[2 * i + 1], 0, uncs[2 * i + 1], ymax);
        lLo->SetLineColor(kOrange + 1);
        lLo->SetLineStyle(3);
        lLo->SetLineWidth(2);
        lHi->SetLineColor(kOrange + 1);
        lHi->SetLineStyle(3);
        lHi->SetLineWidth(2);
        lLo->Draw();
        lHi->Draw();

        auto *leg = new TLegend(0.52, 0.62, 0.88, 0.88);
        leg->SetBorderSize(0);
        leg->SetFillStyle(0);
        leg->AddEntry(hb, "Bootstrap", "l");
        leg->AddEntry(ht, "Truth resample", "l");
        leg->AddEntry(lNom, "Nominal fit", "l");
        leg->AddEntry(lTrue, "True value", "l");
        leg->AddEntry(lLo, "16% / 84% q.", "l");
        leg->Draw();
    }
    cmp->Update();

    // -----------------------------------------------------------------------
    // Canvas con el ajuste sobre los pseudo-datos (RooPlot)
    // -----------------------------------------------------------------------
    std::vector<int> colors = {kBlue, kMagenta, kGreen + 2};
    std::vector<double> lo16(nFit), hi84(nFit);
    for (int k = 0; k < nFit; k++)
    {
        lo16[k] = uncs[2 * k];
        hi84[k] = uncs[2 * k + 1];
    }

    // Chi2 del ajuste nominal
    std::vector<RooRealVar *> dummyCoefs;
    std::unique_ptr<RooAddPdf> modelFixed(
        buildFixedModel(theoryFit, params, dummyCoefs, "chi2_test"));

    TH1F *hData = static_cast<TH1F *>(binnedData->createHistogram("hData_test", x));
    hData->SetDirectory(nullptr);

    double chi2;
    int ndf, nbu;
    computeChi2(hData, modelFixed.get(), x, nFit, "test", chi2, ndf, nbu);

    drawFitCanvas("c_fit_JT_CB",
                  "JT pseudo-data fitted with CB templates",
                  "p_{x} [MeV/c]",
                  hData, binnedData, theoryFit, x,
                  params, lo16, hi84,
                  labelsFit, colors, chi2, ndf, "test");

    std::cout << "\nchi2/ndf = " << chi2 << " / " << ndf
              << " = " << (ndf > 0 ? chi2 / ndf : -1.) << "\n";

    delete hData;
    delete data;
    delete binnedData;
}

// ===========================================================================
// 2) fitMomdis(): real-data analysis
//
// Reads px_rf_rot and py_rf_rot from the experimental ROOT file, fits each
// to the theoretical templates, and obtains asymmetric uncertainties
// from a Poisson bootstrap of the observed histogram. Also computes the
// analytic chi2 and draws residuals.
//
// 'theoFormat' selects the theory-file reader:
//   0 = legacy two-column  (Qi Qt)        -> uses the Qt column
//   1 = JT multi-column    (Qi Qz Qt Qy Q)-> uses the Qy column for
//                                            BOTH the px and py fits
// ===========================================================================
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
    const int nBinsRef = 30;
    const double maxBin = 300.0;
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