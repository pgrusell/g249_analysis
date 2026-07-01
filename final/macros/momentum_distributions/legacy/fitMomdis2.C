#include <RooRealVar.h>
#include <RooDataHist.h>
#include <RooHistPdf.h>
#include <RooAddPdf.h>
#include <RooFormulaVar.h>
#include <RooArgList.h>
#include <RooArgSet.h>
#include <RooPlot.h>
#include <RooFitResult.h>
#include <RooMsgService.h>
#include <RooChi2Var.h>
#include <TRandom3.h>

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
        momdis.Qt->SetBinContent(i + 1, Qt[i]);

    return momdis;
}

// ---------------------------------------------------------------------------
// Read px_rf_rot from the experimental tree (identical to recursive version)
// ---------------------------------------------------------------------------
TH1F *getMomentaDistFromRoot(std::string rootFile, bool useErelCut,
                             double erelMin = 0.2, double erelMax = 1.,
                             int nBins = 40, double maxBin = 300)
{
    auto *h = new TH1F("hExp", "hExp", nBins, -maxBin, maxBin);

    auto *f = new TFile(rootFile.c_str());
    auto *tr = static_cast<TTree *>(f->Get("KinTree"));

    double px;
    tr->SetBranchAddress("px_rf_rot", &px);

    double erel = 0.0;
    if (useErelCut)
        tr->SetBranchAddress("Erel", &erel);

    double opa;
    tr->SetBranchAddress("califa_opa", &opa);

    double sum = 0.0;
    Long64_t n = 0;
    const Long64_t nEntries = tr->GetEntries();

    double opamin = 1.3;
    double opamax = 1.6;

    for (Long64_t i = 0; i < nEntries; ++i)
    {
        tr->GetEntry(i);

        if (useErelCut && (erel * 1000.0 < erelMin || erel * 1000.0 > erelMax))
            continue;

        if (opa > opamin && opa < opamax)
        {
            sum += px * 1000.0;
            ++n;
        }
    }

    const double offset = (n > 0) ? sum / n : 0.0;
    if (n == 0)
        std::cerr << "[getMomentaDistFromRoot] WARNING: no events passed the cuts.\n";

    for (Long64_t i = 0; i < nEntries; ++i)
    {
        tr->GetEntry(i);

        if (useErelCut && (erel * 1000.0 < erelMin || erel * 1000.0 > erelMax))
            continue;

        if (opa > opamin && opa < opamax)
            h->Fill(px * 1000.0 - offset);
    }

    return h;
}

TH1F *buildTemplate(TH1F *hTheo, TH1F *hExp, std::string name)
{
    const int nBins = hExp->GetNbinsX();
    const double xmin = hExp->GetXaxis()->GetXmin();
    const double xmax = hExp->GetXaxis()->GetXmax();

    TH1F *hOut = new TH1F(name.c_str(), name.c_str(), nBins, xmin, xmax);

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

// =============================================================================
// MAIN
// =============================================================================
void fitMomdis2(double erelMin = 4, double erelMax = 7.5)
{
    using namespace RooFit;
    RooMsgService::instance().setGlobalKillBelow(RooFit::WARNING);

    // =====================================================================
    // Inputs
    // =====================================================================
    std::vector<std::string> inFilesTheo = {
        "/nucl_lustre/pablogrusell/g249/g249_analysis/theory/JT/25F/sigt_d52-gs.txt",
        "/nucl_lustre/pablogrusell/g249/g249_analysis/theory/JT/25F/sigt_s12-gs.txt",
        "/nucl_lustre/pablogrusell/g249/g249_analysis/theory/JT/25F/sigt_p12-gs.txt"};

    std::vector<std::string> labels = {"1d_{5/2}", "1s_{1/2}", "1p_{1/2}"};
    std::vector<int> colors = {kBlue, kGreen + 2, kMagenta, kCyan + 1, kOrange + 7, kViolet};

    std::string inFileExp = "/nucl_lustre/pablogrusell/g249/g249_analysis/results/final/23O_analyzed_test.root";

    const int nC = (int)inFilesTheo.size();
    if ((int)labels.size() != nC)
    {
        std::cerr << "ERROR: labels.size() != inFilesTheo.size() ("
                  << labels.size() << " vs " << nC << ")\n";
        return;
    }
    const bool useErelCut = !(erelMin < 0 && erelMax < 0);

    auto *h_exp = getMomentaDistFromRoot(inFileExp, useErelCut, erelMin, erelMax);

    std::vector<TH1F *> templates(nC);
    for (int k = 0; k < nC; k++)
    {
        auto momdis = getMomentaDistFromtxt(inFilesTheo[k], Form("theo_%d", k));
        templates[k] = buildTemplate(momdis.Qt, h_exp, Form("tmpl_%d", k));
    }

    // =====================================================================
    // Build the EXTENDED RooFit model
    //
    // c_k : free yield for component k, range [-0.05*Ntot, 10*Ntot].
    //       Slight negative lower bound so MINOS can find the lower contour
    //       for components at the physical boundary (truth ~ 0).
    // model = RooAddPdf(pdfList, yieldList) — extended automatically when
    //         #coeffs == #pdfs.
    // Ntot  : RooFormulaVar = sum(c_k).
    // phi_k : RooFormulaVar = c_k / Ntot.
    // =====================================================================
    const double xmin = h_exp->GetXaxis()->GetXmin();
    const double xmax = h_exp->GetXaxis()->GetXmax();
    const double Nseed = std::max(1.0, h_exp->Integral());

    RooRealVar x("x", "p_{x} [MeV/c]", xmin, xmax);

    std::vector<RooDataHist *> dh_tmpl(nC);
    std::vector<RooHistPdf *> pdf_tmpl(nC);
    for (int k = 0; k < nC; k++)
    {
        dh_tmpl[k] = new RooDataHist(Form("dh_tmpl_%d", k), Form("dh_tmpl_%d", k),
                                     RooArgList(x), templates[k]);
        pdf_tmpl[k] = new RooHistPdf(Form("pdf_tmpl_%d", k), Form("pdf %s", labels[k].c_str()),
                                     RooArgSet(x), *dh_tmpl[k], 0);
    }

    RooArgList pdfList;
    for (int k = 0; k < nC; k++)
        pdfList.add(*pdf_tmpl[k]);

    // Free yields
    std::vector<RooRealVar *> c(nC, nullptr);
    RooArgList yieldList;
    for (int k = 0; k < nC; k++)
    {
        const double init = Nseed / nC;
        c[k] = new RooRealVar(Form("c%d", k),
                              Form("yield %s", labels[k].c_str()),
                              init, 0., 10.0 * Nseed);
        yieldList.add(*c[k]);
    }

    RooAddPdf *model = new RooAddPdf("model", "extended sum of templates",
                                     pdfList, yieldList);
    RooAbsPdf *fitModel = model;

    // Total yield Ntot = sum(c_k)
    RooFormulaVar *Ntot_fv = nullptr;
    {
        std::string sumExpr;
        RooArgList sumDeps;
        for (int k = 0; k < nC; k++)
        {
            if (k > 0)
                sumExpr += "+";
            sumExpr += Form("@%d", k);
            sumDeps.add(*c[k]);
        }
        Ntot_fv = new RooFormulaVar("Ntot", "total yield",
                                    sumExpr.c_str(), sumDeps);
    }

    // Physical fractions phi_k = c_k / sum(c_j)
    std::vector<RooFormulaVar *> phi(nC, nullptr);
    for (int k = 0; k < nC; k++)
    {
        std::string num = Form("@%d", k);
        std::string den;
        RooArgList deps;
        for (int j = 0; j < nC; j++)
        {
            if (j > 0)
                den += "+";
            den += Form("@%d", j);
            deps.add(*c[j]);
        }
        const std::string expr = num + "/(" + den + ")";
        phi[k] = new RooFormulaVar(Form("phi_%d", k),
                                   Form("physical fraction %s", labels[k].c_str()),
                                   expr.c_str(), deps);
    }

    // =====================================================================
    // Nominal fit on real data (px)
    // =====================================================================
    RooDataHist data("data", "experimental data", RooArgList(x), h_exp);

    std::unique_ptr<RooFitResult> fitRes(
        fitModel->fitTo(data,
                        Extended(true),
                        Save(true),
                        PrintLevel(-1),
                        Minos(true)));

    std::vector<double> frac_nom(nC, 0.0);
    std::vector<double> frac_err(nC, 0.0);
    for (int k = 0; k < nC; k++)
    {
        frac_nom[k] = phi[k]->getVal();
        frac_err[k] = (fitRes ? phi[k]->getPropagatedError(*fitRes) : 0.0);
    }

    const double Nfit = Ntot_fv->getVal();
    const double Nfit_err = (fitRes ? Ntot_fv->getPropagatedError(*fitRes) : 0.0);
    const double Ndata = h_exp->Integral();

    std::cout << "\n=== Nominal RooFit result (extended ML, px) ===\n";
    std::cout << "N (data)  = " << Ndata << "\n";
    std::cout << "N (fit)   = " << Nfit << "  +/- " << Nfit_err << "\n";
    std::cout << "--- yields ---\n";
    for (int k = 0; k < nC; k++)
        std::cout << "c_" << k << " (" << labels[k] << ") = "
                  << c[k]->getVal()
                  << "  +/- " << c[k]->getError()
                  << "  (MINOS +" << c[k]->getAsymErrorHi()
                  << " / " << c[k]->getAsymErrorLo() << ")\n";
    std::cout << "--- physical fractions ---\n";
    for (int k = 0; k < nC; k++)
        std::cout << "phi(" << labels[k] << ") = " << frac_nom[k]
                  << " +/- " << frac_err[k] << " (propagated)\n";
    if (fitRes)
        std::cout << "MIGRAD status = " << fitRes->status()
                  << "  covQual = " << fitRes->covQual()
                  << "  minNll = " << fitRes->minNll() << "\n";

    // ---------------------------------------------------------------------
    // chi2 manual:  sum_i (N_obs_i - N_exp_i)^2 / N_exp_i
    // For an extended model, N_exp_i = N_fit * pdf_int(bin), where N_fit is
    // the total yield from the fit (NOT the data integral). With Extended,
    // sum_i N_exp_i is itself a fit result, not a constraint, so:
    //   ndf = N_bins_used - nC          (nC free yields, no extra constraint)
    // ---------------------------------------------------------------------
    double chi2 = 0.0;
    int nBinsUsed = 0;
    for (int i = 1; i <= h_exp->GetNbinsX(); i++)
    {
        const double lo = h_exp->GetXaxis()->GetBinLowEdge(i);
        const double hi = h_exp->GetXaxis()->GetBinUpEdge(i);

        x.setRange("binChi2", lo, hi);
        RooAbsReal *binInt = fitModel->createIntegral(RooArgSet(x), NormSet(RooArgSet(x)),
                                                      Range("binChi2"));
        const double Nexp = Nfit * binInt->getVal();
        delete binInt;

        if (Nexp <= 0)
            continue;

        const double Nobs = h_exp->GetBinContent(i);
        const double d = Nobs - Nexp;
        chi2 += (d * d) / Nexp;
        nBinsUsed++;
    }

    const int ndf = nBinsUsed - nC;
    std::cout << "chi2 = " << chi2 << "  ndf = " << ndf
              << "  (= " << nBinsUsed << " bins - " << nC
              << " free yields)\n"
              << "  chi2/ndf = " << (ndf > 0 ? chi2 / ndf : -1.0) << "\n";
    std::cout << "================================================\n\n";

    // =====================================================================
    // Plots (px) — same layout as the recursive version
    // =====================================================================
    TCanvas *c1 = new TCanvas("c1", "RooFit template fit (extended)", 800, 800);

    TPad *padTop = new TPad("padTop", "padTop", 0.0, 0.30, 1.0, 1.0);
    padTop->SetBottomMargin(0.02);
    padTop->Draw();
    TPad *padBot = new TPad("padBot", "padBot", 0.0, 0.0, 1.0, 0.30);
    padBot->SetTopMargin(0.02);
    padBot->SetBottomMargin(0.32);
    padBot->Draw();

    padTop->cd();

    RooPlot *xframe = x.frame(Title("Template fit (extended ML)"));
    data.plotOn(xframe, Name("data"));
    fitModel->plotOn(xframe, LineColor(kRed), LineWidth(2), Name("total"));

    for (int k = 0; k < nC; k++)
        model->plotOn(xframe,
                      Components(*pdf_tmpl[k]),
                      LineStyle(kDashed),
                      LineColor(colors[k % colors.size()]),
                      LineWidth(2),
                      Name(Form("comp_%d", k)));

    xframe->GetXaxis()->SetLabelSize(0.0);
    xframe->GetXaxis()->SetTitleSize(0.0);
    xframe->Draw();

    auto *leg = new TLegend(0.55, 0.60, 0.89, 0.89);
    leg->AddEntry(xframe->findObject("data"), "Exp", "ep");
    leg->AddEntry(xframe->findObject("total"),
                  Form("Fit total (N = %.0f #pm %.0f)", Nfit, Nfit_err), "l");
    for (int k = 0; k < nC; k++)
        leg->AddEntry(xframe->findObject(Form("comp_%d", k)),
                      Form("%s (%.1f #pm %.1f)%%",
                           labels[k].c_str(), frac_nom[k] * 100, frac_err[k] * 100),
                      "l");
    leg->AddEntry((TObject *)nullptr,
                  Form("#chi^{2}/ndf = %.2f/%d = %.2f", chi2, ndf,
                       (ndf > 0 ? chi2 / ndf : -1.0)),
                  "");
    leg->Draw();

    // Residual plot
    padBot->cd();

    TH1F *h_res = new TH1F("h_res", ";p_{x} [MeV/c];(N_{obs} - N_{exp})/#sqrt{N_{exp}}",
                           h_exp->GetNbinsX(),
                           h_exp->GetXaxis()->GetXmin(),
                           h_exp->GetXaxis()->GetXmax());

    double maxAbsR = 0.0;
    for (int i = 1; i <= h_exp->GetNbinsX(); i++)
    {
        const double lo = h_exp->GetXaxis()->GetBinLowEdge(i);
        const double hi = h_exp->GetXaxis()->GetBinUpEdge(i);

        x.setRange("binRange", lo, hi);
        RooAbsReal *binInt = fitModel->createIntegral(RooArgSet(x), NormSet(RooArgSet(x)),
                                                      Range("binRange"));
        const double Nexp = Nfit * binInt->getVal();
        delete binInt;

        const double Nobs = h_exp->GetBinContent(i);
        if (Nexp <= 0)
            continue;

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

    TLine *l0 = new TLine(h_exp->GetXaxis()->GetXmin(), 0,
                          h_exp->GetXaxis()->GetXmax(), 0);
    l0->SetLineColor(kRed);
    l0->SetLineWidth(2);
    l0->Draw();
    for (int s : {-1, 1})
    {
        TLine *ls = new TLine(h_exp->GetXaxis()->GetXmin(), s,
                              h_exp->GetXaxis()->GetXmax(), s);
        ls->SetLineColor(kGray + 1);
        ls->SetLineStyle(2);
        ls->Draw();
    }

    c1->cd();

    // =====================================================================
    // py — same procedure
    // =====================================================================
    const double opamin = 1.3;
    const double opamax = 1.6;
    const int nBins = h_exp->GetNbinsX();
    const double maxBin = h_exp->GetXaxis()->GetXmax();

    auto *f_py = new TFile(inFileExp.c_str());
    auto *tr_py = static_cast<TTree *>(f_py->Get("KinTree"));

    double py;
    tr_py->SetBranchAddress("py_rf_rot", &py);

    double erel_py = 0.0;
    if (useErelCut)
        tr_py->SetBranchAddress("Erel", &erel_py);

    double opa_py;
    tr_py->SetBranchAddress("califa_opa", &opa_py);

    double sum_py = 0.0;
    Long64_t n_py = 0;
    const Long64_t nEntries_py = tr_py->GetEntries();

    for (Long64_t i = 0; i < nEntries_py; ++i)
    {
        tr_py->GetEntry(i);

        if (useErelCut && (erel_py * 1000.0 < erelMin || erel_py * 1000.0 > erelMax))
            continue;

        if (opa_py > opamin && opa_py < opamax)
        {
            sum_py += py * 1000.0;
            ++n_py;
        }
    }

    const double offset_py = (n_py > 0) ? sum_py / n_py : 0.0;
    auto *h_exp_py = new TH1F("hExp_py", "hExp_py", nBins, -maxBin, maxBin);

    for (Long64_t i = 0; i < nEntries_py; ++i)
    {
        tr_py->GetEntry(i);

        if (useErelCut && (erel_py * 1000.0 < erelMin || erel_py * 1000.0 > erelMax))
            continue;

        if (opa_py > opamin && opa_py < opamax)
            h_exp_py->Fill(py * 1000.0 - offset_py);
    }

    // Refit on py — keep current c_k values from the px fit as seed
    RooDataHist data_py("data_py", "experimental data py", RooArgList(x), h_exp_py);

    std::unique_ptr<RooFitResult> fitRes_py(
        fitModel->fitTo(data_py,
                        Extended(true),
                        Save(true),
                        PrintLevel(-1),
                        Minos(true)));

    std::vector<double> frac_nom_py(nC, 0.0);
    std::vector<double> frac_err_py(nC, 0.0);
    for (int k = 0; k < nC; k++)
    {
        frac_nom_py[k] = phi[k]->getVal();
        frac_err_py[k] = (fitRes_py ? phi[k]->getPropagatedError(*fitRes_py) : 0.0);
    }

    const double Nfit_py = Ntot_fv->getVal();
    const double Nfit_py_err = (fitRes_py ? Ntot_fv->getPropagatedError(*fitRes_py) : 0.0);
    const double Ndata_py = h_exp_py->Integral();

    std::cout << "\n=== Nominal RooFit result (extended ML, py) ===\n";
    std::cout << "N (data)  = " << Ndata_py << "\n";
    std::cout << "N (fit)   = " << Nfit_py << "  +/- " << Nfit_py_err << "\n";
    std::cout << "--- yields ---\n";
    for (int k = 0; k < nC; k++)
        std::cout << "c_" << k << " (" << labels[k] << ") = "
                  << c[k]->getVal()
                  << "  +/- " << c[k]->getError()
                  << "  (MINOS +" << c[k]->getAsymErrorHi()
                  << " / " << c[k]->getAsymErrorLo() << ")\n";
    std::cout << "--- physical fractions ---\n";
    for (int k = 0; k < nC; k++)
        std::cout << "phi(" << labels[k] << ") = " << frac_nom_py[k]
                  << " +/- " << frac_err_py[k] << " (propagated)\n";
    if (fitRes_py)
        std::cout << "MIGRAD status = " << fitRes_py->status()
                  << "  covQual = " << fitRes_py->covQual()
                  << "  minNll = " << fitRes_py->minNll() << "\n";

    // chi2 (py)
    double chi2_py = 0.0;
    int nBinsUsed_py = 0;
    for (int i = 1; i <= h_exp_py->GetNbinsX(); i++)
    {
        const double lo = h_exp_py->GetXaxis()->GetBinLowEdge(i);
        const double hi = h_exp_py->GetXaxis()->GetBinUpEdge(i);
        x.setRange("binChi2_py", lo, hi);
        RooAbsReal *binInt = fitModel->createIntegral(RooArgSet(x), NormSet(RooArgSet(x)),
                                                      Range("binChi2_py"));
        const double Nexp = Nfit_py * binInt->getVal();
        delete binInt;
        if (Nexp <= 0)
            continue;
        const double Nobs = h_exp_py->GetBinContent(i);
        const double d = Nobs - Nexp;
        chi2_py += (d * d) / Nexp;
        nBinsUsed_py++;
    }
    const int ndf_py = nBinsUsed_py - nC;
    std::cout << "chi2 = " << chi2_py << "  ndf = " << ndf_py
              << "  (= " << nBinsUsed_py << " bins - " << nC << " free yields)\n"
              << "  chi2/ndf = " << (ndf_py > 0 ? chi2_py / ndf_py : -1.0) << "\n";
    std::cout << "================================================\n\n";

    // ---- Canvas for py ----
    TCanvas *c_py = new TCanvas("c_py", "RooFit template fit (extended, py)", 800, 800);

    TPad *padTop_py = new TPad("padTop_py", "padTop_py", 0.0, 0.30, 1.0, 1.0);
    padTop_py->SetBottomMargin(0.02);
    padTop_py->Draw();

    TPad *padBot_py = new TPad("padBot_py", "padBot_py", 0.0, 0.0, 1.0, 0.30);
    padBot_py->SetTopMargin(0.02);
    padBot_py->SetBottomMargin(0.32);
    padBot_py->Draw();

    padTop_py->cd();

    RooPlot *xframe_py = x.frame(Title("Template fit (extended ML, py)"));
    data_py.plotOn(xframe_py, Name("data_py"));
    fitModel->plotOn(xframe_py, LineColor(kRed), LineWidth(2), Name("total_py"));

    for (int k = 0; k < nC; k++)
        model->plotOn(xframe_py,
                      Components(*pdf_tmpl[k]),
                      LineStyle(kDashed),
                      LineColor(colors[k % colors.size()]),
                      LineWidth(2),
                      Name(Form("comp_py_%d", k)));

    xframe_py->GetXaxis()->SetLabelSize(0.0);
    xframe_py->GetXaxis()->SetTitleSize(0.0);
    xframe_py->Draw();

    auto *leg_py = new TLegend(0.55, 0.60, 0.89, 0.89);
    leg_py->AddEntry(xframe_py->findObject("data_py"), "Exp", "ep");
    leg_py->AddEntry(xframe_py->findObject("total_py"),
                     Form("Fit total (N = %.0f #pm %.0f)", Nfit_py, Nfit_py_err), "l");
    for (int k = 0; k < nC; k++)
        leg_py->AddEntry(xframe_py->findObject(Form("comp_py_%d", k)),
                         Form("%s (%.1f #pm %.1f)%%",
                              labels[k].c_str(), frac_nom_py[k] * 100, frac_err_py[k] * 100),
                         "l");
    leg_py->AddEntry((TObject *)nullptr,
                     Form("#chi^{2}/ndf = %.2f/%d = %.2f", chi2_py, ndf_py,
                          (ndf_py > 0 ? chi2_py / ndf_py : -1.0)),
                     "");
    leg_py->Draw();

    // Residual plot for py
    padBot_py->cd();

    TH1F *h_res_py = new TH1F("h_res_py", ";p_{y} [MeV/c];(N_{obs} - N_{exp})/#sqrt{N_{exp}}",
                              h_exp_py->GetNbinsX(),
                              h_exp_py->GetXaxis()->GetXmin(),
                              h_exp_py->GetXaxis()->GetXmax());

    maxAbsR = 0.0;
    for (int i = 1; i <= h_exp_py->GetNbinsX(); i++)
    {
        const double lo = h_exp_py->GetXaxis()->GetBinLowEdge(i);
        const double hi = h_exp_py->GetXaxis()->GetBinUpEdge(i);

        x.setRange("binRange_py", lo, hi);
        RooAbsReal *binInt = fitModel->createIntegral(RooArgSet(x), NormSet(RooArgSet(x)),
                                                      Range("binRange_py"));
        const double Nexp = Nfit_py * binInt->getVal();
        delete binInt;

        const double Nobs = h_exp_py->GetBinContent(i);
        if (Nexp <= 0)
            continue;

        const double r = (Nobs - Nexp) / std::sqrt(Nexp);
        h_res_py->SetBinContent(i, r);
        h_res_py->SetBinError(i, 1.0);
        if (std::abs(r) > maxAbsR)
            maxAbsR = std::abs(r);
    }

    const double yRange_py = std::max(3.5, std::ceil(maxAbsR) + 0.5);
    h_res_py->GetYaxis()->SetRangeUser(-yRange_py, yRange_py);
    h_res_py->GetYaxis()->SetNdivisions(505);
    h_res_py->GetYaxis()->SetTitleSize(0.10);
    h_res_py->GetYaxis()->SetTitleOffset(0.45);
    h_res_py->GetYaxis()->SetLabelSize(0.09);
    h_res_py->GetXaxis()->SetTitleSize(0.12);
    h_res_py->GetXaxis()->SetTitleOffset(1.10);
    h_res_py->GetXaxis()->SetLabelSize(0.10);
    h_res_py->SetMarkerStyle(20);
    h_res_py->SetMarkerSize(0.7);
    h_res_py->SetLineColor(kBlack);
    h_res_py->Draw("E1");

    TLine *l0_py = new TLine(h_exp_py->GetXaxis()->GetXmin(), 0,
                             h_exp_py->GetXaxis()->GetXmax(), 0);
    l0_py->SetLineColor(kRed);
    l0_py->SetLineWidth(2);
    l0_py->Draw();
    for (int s : {-1, 1})
    {
        TLine *ls_py = new TLine(h_exp_py->GetXaxis()->GetXmin(), s,
                                 h_exp_py->GetXaxis()->GetXmax(), s);
        ls_py->SetLineColor(kGray + 1);
        ls_py->SetLineStyle(2);
        ls_py->Draw();
    }
}