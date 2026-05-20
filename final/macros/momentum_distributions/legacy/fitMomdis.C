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
// Read px_rf_rot from the experimental tree.
// erelMin, erelMax: cuts on Erel (in MeV, applied as Erel*1000 from the tree
//                   which stores it in GeV — keep the original convention).
//                   Set BOTH to -1 to disable the Erel cut entirely.
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

    double opamin = 1.3; // 1.25
    double opamax = 1.6; // 1.65

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

// Small helper to compute summary stats of a bootstrap sample
std::tuple<double, double, double, double, double> toyStats(std::vector<double> v)
{
    std::sort(v.begin(), v.end());
    const int n = v.size();
    double mean = 0.0;
    for (double x : v)
        mean += x;
    mean /= n;
    double var = 0.0;
    for (double x : v)
        var += (x - mean) * (x - mean);
    const double sigma = std::sqrt(var / (n - 1));
    const double p16 = v[(int)(0.16 * n)];
    const double p50 = v[(int)(0.50 * n)];
    const double p84 = v[(int)(0.84 * n)];
    return std::make_tuple(mean, sigma, p16, p50, p84);
}

// ---------------------------------------------------------------------------
// Convert recursive fractions (a_0, a_1, ..., a_{N-2}) into physical fractions
// (phi_0, ..., phi_{N-1}) following the RooAddPdf recursive convention:
//
//   phi_0     = a_0
//   phi_1     = (1 - a_0) * a_1
//   phi_2     = (1 - a_0) * (1 - a_1) * a_2
//   ...
//   phi_{N-1} = (1 - a_0) * (1 - a_1) * ... * (1 - a_{N-2})
// ---------------------------------------------------------------------------
std::vector<double> recursiveToPhysical(const std::vector<double> &a)
{
    const int N = (int)a.size() + 1;
    std::vector<double> phi(N);
    double remain = 1.0;
    for (int k = 0; k < N - 1; k++)
    {
        phi[k] = remain * a[k];
        remain *= (1.0 - a[k]);
    }
    phi[N - 1] = remain;
    return phi;
}

// void fitMomdis(double erelMin = 2.2, double erelMax = 4)
void fitMomdis(double erelMin = 4, double erelMax = 7.5)
{
    using namespace RooFit;

    RooMsgService::instance().setGlobalKillBelow(RooFit::WARNING);

    // =====================================================================
    // Inputs — reduce inFilesTheo to a single entry for a single-component fit
    // =====================================================================
    std::vector<std::string> inFilesTheo = {
        "/nucl_lustre/pablogrusell/g249/g249_analysis/theory/JT/25F/sigt_d52-gs.txt",
        "/nucl_lustre/pablogrusell/g249/g249_analysis/theory/JT/25F/sigt_s12-gs.txt",
        "/nucl_lustre/pablogrusell/g249/g249_analysis/theory/JT/25F/sigt_p12-gs.txt"};

    std::vector<std::string> labels = {"1d_{5/2}", "1s_{1/2}", "1p_{1/2}"};
    std::vector<int> colors = {kBlue, kGreen + 2, kMagenta, kCyan + 1, kOrange + 7, kViolet};

    std::string inFileExp = "/nucl_lustre/pablogrusell/g249/g249_analysis/results/final/23O_analyzed_test.root";

    const int nC = (int)inFilesTheo.size();
    const bool useErelCut = !(erelMin < 0 && erelMax < 0);

    // Read experimental histogram and build templates in the same binning
    auto *h_exp = getMomentaDistFromRoot(inFileExp, useErelCut, erelMin, erelMax);

    std::vector<TH1F *> templates(nC);
    for (int k = 0; k < nC; k++)
    {
        auto momdis = getMomentaDistFromtxt(inFilesTheo[k], Form("theo_%d", k));
        templates[k] = buildTemplate(momdis.Qt, h_exp, Form("tmpl_%d", k));
    }

    // =====================================================================
    // Build the RooFit model
    //
    // nC == 1: fitModel = single RooHistPdf, no free fraction parameters.
    // nC  > 1: fitModel = RooAddPdf with recursive fractions (original logic).
    // =====================================================================
    const double xmin = h_exp->GetXaxis()->GetXmin();
    const double xmax = h_exp->GetXaxis()->GetXmax();

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

    // Recursive fractions and composite model (only built when nC > 1)
    std::vector<RooRealVar *> a(nC > 1 ? nC - 1 : 0, nullptr);
    RooArgList fracList, pdfList;
    for (int k = 0; k < nC; k++)
        pdfList.add(*pdf_tmpl[k]);

    RooAddPdf *model = nullptr;    // non-null only when nC > 1
    RooAbsPdf *fitModel = nullptr; // always set below

    if (nC > 1)
    {
        for (int k = 0; k < nC - 1; k++)
        {
            const double init = 1.0 / (nC - k);
            a[k] = new RooRealVar(Form("a%d", k), Form("recursive frac a_%d", k), init, 0.0, 1.0);
            fracList.add(*a[k]);
        }
        model = new RooAddPdf("model", "sum of templates", pdfList, fracList, /*recursiveFractions=*/true);
        fitModel = model;
    }
    else
    {
        fitModel = pdf_tmpl[0];
    }

    // Physical fraction formulas — propagate covariance (only needed when nC > 1)
    std::vector<RooFormulaVar *> phi(nC > 1 ? nC : 0, nullptr);
    if (nC > 1)
    {
        for (int k = 0; k < nC; k++)
        {
            std::string expr;
            RooArgList deps;
            if (k < nC - 1)
            {
                // phi_k = (1-a_0)*...*(1-a_{k-1}) * a_k
                for (int j = 0; j < k; j++)
                {
                    expr += Form("(1-@%d)*", j);
                    deps.add(*a[j]);
                }
                expr += Form("@%d", k);
                deps.add(*a[k]);
            }
            else
            {
                // phi_{N-1} = prod_{j=0}^{N-2} (1-a_j)
                for (int j = 0; j < nC - 1; j++)
                {
                    if (j > 0)
                        expr += "*";
                    expr += Form("(1-@%d)", j);
                    deps.add(*a[j]);
                }
            }
            phi[k] = new RooFormulaVar(Form("phi_%d", k),
                                       Form("physical fraction %s", labels[k].c_str()),
                                       expr.c_str(), deps);
        }
    }

    // =====================================================================
    // Nominal fit on real data
    // =====================================================================
    RooDataHist data("data", "experimental data", RooArgList(x), h_exp);

    // nC == 1: no free parameters, nothing to minimize — skip fitTo
    std::unique_ptr<RooFitResult> fitRes;
    if (nC > 1)
        fitRes.reset(fitModel->fitTo(data, Save(true), PrintLevel(-1), Minos(true)));

    // Physical fractions: trivially 1.0 when nC == 1
    std::vector<double> frac_nom(nC, 1.0);
    std::vector<double> frac_err(nC, 0.0);
    if (nC > 1)
        for (int k = 0; k < nC; k++)
        {
            frac_nom[k] = phi[k]->getVal();
            frac_err[k] = phi[k]->getPropagatedError(*fitRes);
        }

    const double Ntot = h_exp->Integral();

    std::cout << "\n=== Nominal RooFit result ===\n";
    std::cout << "N (data) = " << Ntot << "\n";
    if (nC > 1)
    {
        for (int k = 0; k < nC - 1; k++)
            std::cout << "a_" << k << " = " << a[k]->getVal() << " +/- " << a[k]->getError() << "\n";
    }
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
    // Bins with N_exp <= 0 are skipped.
    //
    // Degrees of freedom:
    //   ndf = N_bins_used - nC
    // The fit has (nC - 1) free fraction parameters, plus 1 implicit
    // constraint from fixing the total normalization to the data
    // (N_exp_i is built as Ntot * pdf_int, so sum_i N_exp_i = sum_i N_obs_i
    // is enforced by construction, consuming one extra degree of freedom).
    // For nC == 1 this gives ndf = N_bins_used - 1, consistent with the
    // single-component case where only the global normalization is fixed.
    // ---------------------------------------------------------------------
    const double Ntot_for_chi2 = h_exp->Integral();
    double chi2 = 0.0;
    int nBinsUsed = 0;
    for (int i = 1; i <= h_exp->GetNbinsX(); i++)
    {
        const double lo = h_exp->GetXaxis()->GetBinLowEdge(i);
        const double hi = h_exp->GetXaxis()->GetBinUpEdge(i);

        x.setRange("binChi2", lo, hi);
        RooAbsReal *binInt = fitModel->createIntegral(RooArgSet(x), NormSet(RooArgSet(x)),
                                                      Range("binChi2"));
        const double Nexp = Ntot_for_chi2 * binInt->getVal();
        delete binInt;

        if (Nexp <= 0)
            continue;

        const double Nobs = h_exp->GetBinContent(i);
        const double d = Nobs - Nexp;
        chi2 += (d * d) / Nexp;
        nBinsUsed++;
    }

    const int nPar = (nC > 1) ? nC - 1 : 0; // free parameters in the fit
    const int ndf = nBinsUsed - nC;         // = nBinsUsed - nPar - 1
    std::cout << "chi2 = " << chi2 << "  ndf = " << ndf
              << "  (= " << nBinsUsed << " bins - " << nPar
              << " free pars - 1 norm constraint)\n"
              << "  chi2/ndf = " << (ndf > 0 ? chi2 / ndf : -1.0) << "\n";
    std::cout << "=============================\n\n";

    // =====================================================================
    // Plots
    // =====================================================================
    TCanvas *c1 = new TCanvas("c1", "RooFit template fit", 800, 800);

    // Top pad: fit
    TPad *padTop = new TPad("padTop", "padTop", 0.0, 0.30, 1.0, 1.0);
    padTop->SetBottomMargin(0.02);
    padTop->Draw();
    // Bottom pad: residuals (pulls)
    TPad *padBot = new TPad("padBot", "padBot", 0.0, 0.0, 1.0, 0.30);
    padBot->SetTopMargin(0.02);
    padBot->SetBottomMargin(0.32);
    padBot->Draw();

    padTop->cd();

    RooPlot *xframe = x.frame(Title("Template fit"));
    data.plotOn(xframe, Name("data"));
    fitModel->plotOn(xframe, LineColor(kRed), LineWidth(2), Name("total"));

    // Individual components drawn dashed (only meaningful when nC > 1)
    if (nC > 1)
        for (int k = 0; k < nC; k++)
            model->plotOn(xframe,
                          Components(*pdf_tmpl[k]),
                          LineStyle(kDashed),
                          LineColor(colors[k % colors.size()]),
                          LineWidth(2),
                          Name(Form("comp_%d", k)));

    xframe->GetXaxis()->SetLabelSize(0.0); // hide x labels on top pad
    xframe->GetXaxis()->SetTitleSize(0.0);
    xframe->Draw();

    auto *leg = new TLegend(0.55, 0.60, 0.89, 0.89);
    leg->AddEntry(xframe->findObject("data"), "Exp", "ep");
    leg->AddEntry(xframe->findObject("total"), "Fit total", "l");
    for (int k = 0; k < nC; k++)
    {
        // When nC == 1 there is no separate "comp_0" curve; point the legend at "total"
        const char *entryName = (nC > 1) ? Form("comp_%d", k) : "total";
        leg->AddEntry(xframe->findObject(entryName),
                      Form("%s (%.1f #pm %.1f)%%",
                           labels[k].c_str(), frac_nom[k] * 100, frac_err[k] * 100),
                      "l");
    }
    leg->AddEntry((TObject *)nullptr,
                  Form("#chi^{2}/ndf = %.2f/%d = %.2f", chi2, ndf,
                       (ndf > 0 ? chi2 / ndf : -1.0)),
                  "");
    leg->Draw();

    // ---------------------------------------------------------------------
    // Residual / pull plot:  r_i = (N_obs - N_exp) / sqrt(N_exp)
    // N_exp is the model PDF integrated over each bin, normalized to Ntot.
    // Bins with N_exp == 0 are skipped (would divide by zero).
    // ---------------------------------------------------------------------
    padBot->cd();

    TH1F *h_res = new TH1F("h_res", ";p_{x} [MeV/c];(N_{obs} - N_{exp})/#sqrt{N_{exp}}",
                           h_exp->GetNbinsX(),
                           h_exp->GetXaxis()->GetXmin(),
                           h_exp->GetXaxis()->GetXmax());

    // Get N_exp per bin by integrating the pdf bin-by-bin
    RooAbsReal *modelInt = fitModel->createIntegral(RooArgSet(x), NormSet(RooArgSet(x)));
    const double normFactor = Ntot; // pdf is normalized to 1 over [xmin, xmax]
    (void)modelInt;                 // total integral = 1, kept for clarity

    double maxAbsR = 0.0;
    for (int i = 1; i <= h_exp->GetNbinsX(); i++)
    {
        const double lo = h_exp->GetXaxis()->GetBinLowEdge(i);
        const double hi = h_exp->GetXaxis()->GetBinUpEdge(i);

        x.setRange("binRange", lo, hi);
        RooAbsReal *binInt = fitModel->createIntegral(RooArgSet(x), NormSet(RooArgSet(x)),
                                                      Range("binRange"));
        const double Nexp = normFactor * binInt->getVal();
        delete binInt;

        const double Nobs = h_exp->GetBinContent(i);
        if (Nexp <= 0)
            continue;

        const double r = (Nobs - Nexp) / std::sqrt(Nexp);
        h_res->SetBinContent(i, r);
        h_res->SetBinError(i, 1.0); // by construction
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

    // Reference lines at 0 and +/- 1 sigma
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
    // Now make the same plots for py_rf_rot
    // =====================================================================
    const double opamin = 1.3;
    const double opamax = 1.6;
    const int nBins = h_exp->GetNbinsX();
    const double maxBin = h_exp->GetXaxis()->GetXmax();

    // Manually read py instead of px
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

    // Rebuild model for py with same templates
    RooDataHist data_py("data_py", "experimental data py", RooArgList(x), h_exp_py);

    // Refit with same initial values as nominal fit
    std::unique_ptr<RooFitResult> fitRes_py;
    if (nC > 1)
    {
        for (int k = 0; k < nC - 1; k++)
            a[k]->setVal(a[k]->getVal()); // Use current value from nominal fit as seed
        fitRes_py.reset(model->fitTo(data_py, Save(true), PrintLevel(-1), Minos(true)));
    }

    // Physical fractions for py
    std::vector<double> frac_nom_py(nC, 1.0);
    std::vector<double> frac_err_py(nC, 0.0);
    if (nC > 1)
        for (int k = 0; k < nC; k++)
        {
            frac_nom_py[k] = phi[k]->getVal();
            frac_err_py[k] = phi[k]->getPropagatedError(*fitRes_py);
        }

    const double Ntot_py = h_exp_py->Integral();

    std::cout << "\n=== Nominal RooFit result (py) ===\n";
    std::cout << "N (data) = " << Ntot_py << "\n";
    if (nC > 1)
    {
        for (int k = 0; k < nC - 1; k++)
            std::cout << "a_" << k << " = " << a[k]->getVal() << " +/- " << a[k]->getError() << "\n";
    }
    std::cout << "--- physical fractions ---\n";
    for (int k = 0; k < nC; k++)
        std::cout << "phi(" << labels[k] << ") = " << frac_nom_py[k]
                  << " +/- " << frac_err_py[k] << " (propagated)\n";
    if (fitRes_py)
        std::cout << "MIGRAD status = " << fitRes_py->status()
                  << "  covQual = " << fitRes_py->covQual()
                  << "  minNll = " << fitRes_py->minNll() << "\n";

    // chi2 for py
    const double Ntot_py_for_chi2 = h_exp_py->Integral();
    double chi2_py = 0.0;
    int nBinsUsed_py = 0;
    for (int i = 1; i <= h_exp_py->GetNbinsX(); i++)
    {
        const double lo = h_exp_py->GetXaxis()->GetBinLowEdge(i);
        const double hi = h_exp_py->GetXaxis()->GetBinUpEdge(i);
        x.setRange("binChi2_py", lo, hi);
        RooAbsReal *binInt_py = fitModel->createIntegral(RooArgSet(x), NormSet(RooArgSet(x)),
                                                         Range("binChi2_py"));
        const double Nexp_py = Ntot_py_for_chi2 * binInt_py->getVal();
        delete binInt_py;
        if (Nexp_py <= 0)
            continue;
        const double d_py = h_exp_py->GetBinContent(i) - Nexp_py;
        chi2_py += (d_py * d_py) / Nexp_py;
        nBinsUsed_py++;
    }
    const int ndf_py = nBinsUsed_py - nC;
    std::cout << "chi2 (py) = " << chi2_py << "  ndf = " << ndf_py
              << "  chi2/ndf = " << (ndf_py > 0 ? chi2_py / ndf_py : -1.0) << "\n";

    // ---- Canvas for py ----
    TCanvas *c_py = new TCanvas("c_py", "RooFit template fit (py)", 800, 800);

    TPad *padTop_py = new TPad("padTop_py", "padTop_py", 0.0, 0.30, 1.0, 1.0);
    padTop_py->SetBottomMargin(0.02);
    padTop_py->Draw();

    TPad *padBot_py = new TPad("padBot_py", "padBot_py", 0.0, 0.0, 1.0, 0.30);
    padBot_py->SetTopMargin(0.02);
    padBot_py->SetBottomMargin(0.32);
    padBot_py->Draw();

    padTop_py->cd();

    RooPlot *xframe_py = x.frame(Title("Template fit (py)"));
    data_py.plotOn(xframe_py, Name("data_py"));
    fitModel->plotOn(xframe_py, LineColor(kRed), LineWidth(2), Name("total_py"));

    if (nC > 1)
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
    leg_py->AddEntry(xframe_py->findObject("total_py"), "Fit total", "l");
    for (int k = 0; k < nC; k++)
    {
        const char *entryName_py = (nC > 1) ? Form("comp_py_%d", k) : "total_py";
        leg_py->AddEntry(xframe_py->findObject(entryName_py),
                         Form("%s (%.1f #pm %.1f)%%",
                              labels[k].c_str(), frac_nom_py[k] * 100, frac_err_py[k] * 100),
                         "l");
    }
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
        const double Nexp = Ntot_py * binInt->getVal();
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