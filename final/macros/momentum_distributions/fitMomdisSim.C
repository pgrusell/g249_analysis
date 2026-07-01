// ===========================================================================
//  fitMomdisSim.C
//
//  Template fit of the experimental momentum distributions (px, py) using
//  SIMULATED templates read from the reconstructed histogram (hPxRestRec)
//  of the sim_analysis ROOT files:
//
//      sim_analysis_d52.root  ->  1d_{5/2} template
//      sim_analysis_s12.root  ->  2s_{1/2} template
//
//  The reconstructed histogram is used (not the truth) so the templates carry
//  the same detector resolution as the experimental data they are fitted to.
//
//  Pipeline (per component):
//      hPxRestRec [GeV/c]  --rescale axis x1000-->  [MeV/c]
//                          --buildTemplate--------> experimental binning
//                          --create_pdf_from_histogram--> RooHistPdf
//  then a recursive-fraction RooAddPdf fit, Poisson bootstrap for asymmetric
//  errors, analytic chi2, and fit + residual + bootstrap canvases.
//
//  Usage:
//      root -l
//      .L fitMomdisSim.C
//      fitMomdisSim(1.5, 4, 1000);                 // default: hPxRestRec
//      fitMomdisSim(1.5, 4, 1000, "hPxRestTru");   // truth, for comparison
// ===========================================================================

// ---------------------------------------------------------------------------
// Read px_rf_rot or py_rf_rot from the experimental tree (in MeV/c, centred).
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
// Read a template distribution directly from a histogram inside a ROOT file
// (e.g. hPxRestRec from sim_analysis_d52.root).
//
// IMPORTANT (units): the sim histograms are in GeV/c (axis -0.5..0.5), while
// the rest of the macro works in MeV/c. We rescale the axis by `unitScale`
// (default 1000: GeV/c -> MeV/c). Content is copied bin-by-bin. buildTemplate()
// later re-bins this onto the experimental binning and normalises, so the
// absolute scale here is irrelevant -- only the *shape* and *axis units* matter.
//
// `recenter` optionally shifts the axis to subtract the histogram mean, the
// same way getMomentaDistFromRoot centres the experimental data. For a px/py
// in the rest frame the template is ~symmetric, so it defaults to false.
// ---------------------------------------------------------------------------
TH1F *getTemplateFromRootHist(const std::string &rootFile,
                              const std::string &histName,
                              const std::string &outName,
                              double unitScale = 1000.0,
                              bool recenter = false)
{
    TFile *f = TFile::Open(rootFile.c_str(), "READ");
    if (!f || f->IsZombie())
    {
        std::cerr << "[getTemplateFromRootHist] cannot open " << rootFile << "\n";
        return nullptr;
    }

    TH1 *hSrc = dynamic_cast<TH1 *>(f->Get(histName.c_str()));
    if (!hSrc)
    {
        std::cerr << "[getTemplateFromRootHist] '" << histName
                  << "' not found in " << rootFile << "\n";
        f->Close();
        delete f;
        return nullptr;
    }

    const int nb = hSrc->GetNbinsX();
    const double xlo = hSrc->GetXaxis()->GetXmin() * unitScale;
    const double xhi = hSrc->GetXaxis()->GetXmax() * unitScale;
    const double shift = recenter ? hSrc->GetMean() * unitScale : 0.0;

    // Source is TH1D in the sim output; copy into TH1F because the downstream
    // code (buildTemplate / RooDataHist) is written for TH1F.
    TH1F *hOut = new TH1F(outName.c_str(), outName.c_str(),
                          nb, xlo - shift, xhi - shift);
    hOut->SetDirectory(nullptr);
    for (int i = 1; i <= nb; ++i)
        hOut->SetBinContent(i, hSrc->GetBinContent(i));

    if (hOut->Integral() <= 0)
        std::cerr << "[getTemplateFromRootHist] WARNING: '" << histName
                  << "' in " << rootFile << " is empty.\n";

    f->Close();
    delete f;
    return hOut;
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

inline std::vector<double> getUncertainties(TH1F *h) { return getQuantiles(h); }

// ---------------------------------------------------------------------------
// Bootstrap: Poisson fluctuations on the observed histogram, refit each toy.
// Fills parDistributions and returns interleaved {q16,q84} per parameter.
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
// Build an RooAddPdf with FIXED recursive coefficients matching the given
// physical fractions. Used to plot the best-fit model and compute chi2.
// ---------------------------------------------------------------------------
RooAddPdf *buildFixedModel(std::vector<RooHistPdf *> theory,
                           const std::vector<double> &phi,
                           std::vector<RooRealVar *> &storedCoefs,
                           std::string tag)
{
    const int N = (int)theory.size();
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
// Fit + residual canvas for one momentum component.
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
// Bootstrap distributions on one canvas.
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
//  Main: fit experimental px / py to simulated d52 / s12 templates.
//
//   simHistName : histogram to use as template inside each sim file.
//                 Default "hPxRestRec" (reconstructed -> carries the detector
//                 resolution, matching the experimental data). Pass
//                 "hPxRestTru" to compare against the unsmeared truth.
//   recenter    : subtract each template's mean before fitting (default off).
// ===========================================================================
void fitMomdisSim(double erelMin = -1, double erelMax = -1,
                  int nToys = 1000,
                  const std::string &simHistName = "hPxRestRec",
                  bool recenter = false,
                  const std::string &simFileD52 = "/nucl_lustre/pablogrusell/g249/g249_analysis/sim/full/sim_analysis_d52.root",
                  const std::string &simFileS12 = "/nucl_lustre/pablogrusell/g249/g249_analysis/sim/full/sim_analysis_s12.root")
{
    using namespace RooFit;
    RooMsgService::instance().setGlobalKillBelow(RooFit::ERROR);

    // ---- Template sources: the two simulated ROOT files ------------------
    std::vector<std::string> simFiles = {simFileD52, simFileS12};
    std::vector<std::string> labels = {"1d_{5/2}", "2s_{1/2}"};
    std::vector<int> colors = {kBlue, kGreen + 2, kMagenta, kCyan + 1};

    // ---- Experimental data file ------------------------------------------
    std::string inFileExp =
        "/nucl_lustre/pablogrusell/g249/g249_analysis/results/final/24O_analyzed_test.root";

    const int nC = (int)simFiles.size();
    const int nBinsRef = 30;
    const double maxBin = 300.0; // MeV/c
    const bool useErelCut = !(erelMin < 0 && erelMax < 0);

    // ---- Experimental histograms for px and py (in MeV/c, centred) -------
    TH1F *h_px = getMomentaDistFromRoot(inFileExp, "px_rf_rot", useErelCut,
                                        erelMin, erelMax, nBinsRef, maxBin,
                                        "h_exp_px");
    TH1F *h_py = getMomentaDistFromRoot(inFileExp, "py_rf_rot", useErelCut,
                                        erelMin, erelMax, nBinsRef, maxBin,
                                        "h_exp_py");

    RooRealVar x("x", "p [MeV/c]", -maxBin, maxBin);

    // ---- Build templates from the simulated histograms -------------------
    // read hPxRestRec (GeV/c -> MeV/c), re-bin onto experimental binning,
    // wrap into a RooHistPdf.
    std::vector<TH1F *> templates(nC);
    std::vector<RooHistPdf *> theory(nC);
    for (int k = 0; k < nC; k++)
    {
        TH1F *raw = getTemplateFromRootHist(simFiles[k], simHistName,
                                            Form("raw_sim_%d", k),
                                            1000.0, recenter);
        if (!raw || raw->Integral() <= 0)
        {
            std::cerr << "Aborting: template " << k << " ("
                      << simFiles[k] << ":" << simHistName
                      << ") is missing or empty.\n";
            return;
        }
        templates[k] = buildTemplate(raw, h_px, Form("tmpl_sim_%d", k));
        theory[k] = create_pdf_from_histogram(x, templates[k]);
        delete raw;
    }

    // =====================================================================
    // PX FIT
    // =====================================================================
    std::cout << "\n############### FIT px (sim templates: " << simHistName
              << ") ###############\n";
    RooDataHist data_px("data_px", "experimental data px",
                        RooArgList(x), h_px);
    std::vector<double> phi_px = make_fit(&data_px, theory, x);

    std::vector<TH1F *> boot_px;
    std::vector<double> q_px = calculateUncertainties(&data_px, theory, x,
                                                      boot_px, nToys, "boot_px_sim");
    std::vector<double> lo16_px(nC), hi84_px(nC);
    for (int k = 0; k < nC; k++)
    {
        lo16_px[k] = q_px[2 * k];
        hi84_px[k] = q_px[2 * k + 1];
    }

    std::vector<RooRealVar *> dummy_px;
    std::unique_ptr<RooAddPdf> model_px(
        buildFixedModel(theory, phi_px, dummy_px, "chi2_px_sim"));
    double chi2_px;
    int ndf_px, nbu_px;
    computeChi2(h_px, model_px.get(), x, nC, "px_sim", chi2_px, ndf_px, nbu_px);

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
    std::cout << "\n############### FIT py (sim templates: " << simHistName
              << ") ###############\n";
    RooDataHist data_py("data_py", "experimental data py",
                        RooArgList(x), h_py);
    std::vector<double> phi_py = make_fit(&data_py, theory, x);

    std::vector<TH1F *> boot_py;
    std::vector<double> q_py = calculateUncertainties(&data_py, theory, x,
                                                      boot_py, nToys, "boot_py_sim");
    std::vector<double> lo16_py(nC), hi84_py(nC);
    for (int k = 0; k < nC; k++)
    {
        lo16_py[k] = q_py[2 * k];
        hi84_py[k] = q_py[2 * k + 1];
    }

    std::vector<RooRealVar *> dummy_py;
    std::unique_ptr<RooAddPdf> model_py(
        buildFixedModel(theory, phi_py, dummy_py, "chi2_py_sim"));
    double chi2_py;
    int ndf_py, nbu_py;
    computeChi2(h_py, model_py.get(), x, nC, "py_sim", chi2_py, ndf_py, nbu_py);

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
    drawFitCanvas("c_px_sim", "Template fit p_{x} (sim templates)", "p_{x} [MeV/c]",
                  h_px, &data_px, theory, x,
                  phi_px, lo16_px, hi84_px,
                  labels, colors, chi2_px, ndf_px, "px_sim");

    drawFitCanvas("c_py_sim", "Template fit p_{y} (sim templates)", "p_{y} [MeV/c]",
                  h_py, &data_py, theory, x,
                  phi_py, lo16_py, hi84_py,
                  labels, colors, chi2_py, ndf_py, "py_sim");

    drawBootstrapCanvas("c_boot_px_sim", "Bootstrap distributions p_{x}",
                        boot_px, phi_px, lo16_px, hi84_px, labels);
    drawBootstrapCanvas("c_boot_py_sim", "Bootstrap distributions p_{y}",
                        boot_py, phi_py, lo16_py, hi84_py, labels);
}