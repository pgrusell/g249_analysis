#include <RooRealVar.h>
#include <RooDataHist.h>
#include <RooHistPdf.h>
#include <RooAddPdf.h>
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

TH1F *getMomentaDistFromRoot(std::string rootFile, int nBins = 50, double maxBin = 300)
{
    auto *h = new TH1F("hExp", "hExp", nBins, -maxBin, maxBin);

    auto *f = new TFile(rootFile.c_str());
    auto *tr = static_cast<TTree *>(f->Get("KinTree"));

    double py;
    tr->SetBranchAddress("py_rf_rot", &py);

    double erel;
    tr->SetBranchAddress("Erel", &erel);

    double sum = 0.0;
    Long64_t n = 0;
    const Long64_t nEntries = tr->GetEntries();

    for (Long64_t i = 0; i < nEntries; ++i)
    {
        tr->GetEntry(i);

        if (erel * 1000 < 7 && erel * 1000 > 5)
        {
            sum += py * 1000.0;
            ++n;
        }
    }

    const double offset = (n > 0) ? sum / n : 0.0;

    for (Long64_t i = 0; i < nEntries; ++i)
    {
        tr->GetEntry(i);

        if (erel * 1000 < 7 && erel * 1000 > 5)
            h->Fill(py * 1000 - offset);
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

void fitMomdis()
{
    using namespace RooFit;

    // Silence RooFit chatter in the bootstrap loop
    RooMsgService::instance().setGlobalKillBelow(RooFit::WARNING);

    // =====================================================================
    // Inputs
    // =====================================================================
    std::vector<std::string> inFilesTheo = {
        "/nucl_lustre/pablogrusell/g249/g249_analysis/theory/JT/25F/sigt_d52-3.txt",
        "/nucl_lustre/pablogrusell/g249/g249_analysis/theory/JT/25F/sigt_p12-3.txt",
        "/nucl_lustre/pablogrusell/g249/g249_analysis/theory/JT/25F/sigt_s12-3.txt"};

    std::vector<std::string> labels = {"1d_{5/2}", "p_{1/2}", "s_{1/2}"};
    std::vector<int> colors = {kBlue, kGreen + 2, kMagenta, kCyan + 1, kOrange + 7, kViolet};

    std::string inFileExp = "/nucl_lustre/pablogrusell/g249/g249_analysis/results/final/23O_analyzed_test.root";

    const int nC = (int)inFilesTheo.size();

    // Read experimental histogram and build templates in the same binning
    auto *h_exp = getMomentaDistFromRoot(inFileExp);

    std::vector<TH1F *> templates(nC);
    for (int k = 0; k < nC; k++)
    {
        auto momdis = getMomentaDistFromtxt(inFilesTheo[k], Form("theo_%d", k));
        templates[k] = buildTemplate(momdis.Qt, h_exp, Form("tmpl_%d", k));
    }

    // =====================================================================
    // Build the RooFit model
    // =====================================================================
    const double xmin = h_exp->GetXaxis()->GetXmin();
    const double xmax = h_exp->GetXaxis()->GetXmax();

    // 1. Observable
    RooRealVar x("x", "p_{x} [MeV/c]", xmin, xmax);

    // 2. Convert each template TH1F into a RooHistPdf
    //    We keep the RooDataHist objects alive in a vector to avoid them going
    //    out of scope while the PDFs still reference them.
    std::vector<RooDataHist *> dh_tmpl(nC);
    std::vector<RooHistPdf *> pdf_tmpl(nC);
    for (int k = 0; k < nC; k++)
    {
        dh_tmpl[k] = new RooDataHist(Form("dh_tmpl_%d", k),
                                     Form("dh_tmpl_%d", k),
                                     RooArgList(x), templates[k]);
        // order=0 -> piecewise-constant interpolation (matches the binned fit)
        pdf_tmpl[k] = new RooHistPdf(Form("pdf_tmpl_%d", k),
                                     Form("pdf %s", labels[k].c_str()),
                                     RooArgSet(x), *dh_tmpl[k], 0);
    }

    // 3. Fractions: nC-1 free parameters, the last is 1 - sum automatically
    std::vector<RooRealVar *> frac(nC - 1);
    RooArgList fracList;
    for (int k = 0; k < nC - 1; k++)
    {
        frac[k] = new RooRealVar(Form("f%d", k + 1),
                                 Form("fraction %s", labels[k].c_str()),
                                 1.0 / nC, 0.0, 1.0);
        fracList.add(*frac[k]);
    }

    // 4. Composite model (recursive=false: last fraction = 1 - sum)
    RooArgList pdfList;
    for (int k = 0; k < nC; k++)
        pdfList.add(*pdf_tmpl[k]);

    RooAddPdf model("model", "sum of templates", pdfList, fracList);

    // =====================================================================
    // Nominal fit on real data
    // =====================================================================
    RooDataHist data("data", "experimental data", RooArgList(x), h_exp);

    // Extended=false because we are fitting fractions, not absolute yields.
    // SumW2Error(true) is appropriate when the data histogram is a count
    // histogram with Poisson statistics.
    std::unique_ptr<RooFitResult> fitRes(
        model.fitTo(data, Save(true), PrintLevel(-1), Minos(true)));

    // Collect nominal fractions (including the derived last one)
    std::vector<double> frac_nom(nC);
    std::vector<double> frac_err(nC);
    double sumF = 0.0;
    for (int k = 0; k < nC - 1; k++)
    {
        frac_nom[k] = frac[k]->getVal();
        frac_err[k] = frac[k]->getError();
        sumF += frac_nom[k];
    }
    frac_nom[nC - 1] = 1.0 - sumF;

    // Error on the last (derived) fraction via error propagation with cov matrix
    {
        const TMatrixDSym &cov = fitRes->covarianceMatrix();
        // variance of sum of free fractions
        double var = 0.0;
        for (int i = 0; i < nC - 1; i++)
            for (int j = 0; j < nC - 1; j++)
                var += cov(i, j);
        frac_err[nC - 1] = std::sqrt(std::max(0.0, var));
    }

    const double Ntot = h_exp->Integral();

    std::cout << "\n=== Nominal RooFit result ===\n";
    std::cout << "N (data) = " << Ntot << "\n";
    for (int k = 0; k < nC; k++)
        std::cout << "f(" << labels[k] << ") = " << frac_nom[k]
                  << " +/- " << frac_err[k] << " (Minos/cov)\n";
    std::cout << "MIGRAD status = " << fitRes->status()
              << "  covQual = " << fitRes->covQual()
              << "  minNll = " << fitRes->minNll() << "\n";

    // Chi2 from the binned model vs binned data
    RooChi2Var chi2Var("chi2Var", "chi2Var", model, data, DataError(RooAbsData::Poisson));
    const double chi2 = chi2Var.getVal();
    const int nPar = nC - 1; // free fractions in RooAddPdf
    const int ndf = h_exp->GetNbinsX() - nPar;
    const double chi2ndf = (ndf > 0) ? (chi2 / ndf) : -1.0;

    std::cout << "chi2 = " << chi2
              << "  ndf = " << ndf
              << "  chi2/ndf = " << chi2ndf << "\n";
    std::cout << "============================\n\n";

    // =====================================================================
    // Bootstrap: Poisson toys over the experimental histogram
    // =====================================================================
    const int nToys = 1000;
    TRandom3 rng(0);

    const int nBins = h_exp->GetNbinsX();
    TH1F *h_toy = (TH1F *)h_exp->Clone("h_toy");

    // Histograms to visualise the bootstrap distributions
    std::vector<TH1F *> hF(nC);
    for (int k = 0; k < nC; k++)
        hF[k] = new TH1F(Form("hF%d", k + 1),
                         Form("Bootstrap f(%s);f_{%d};toys", labels[k].c_str(), k + 1),
                         100, -0.1, 1.1);

    std::vector<std::vector<double>> v_frac(nC);
    for (int k = 0; k < nC; k++)
        v_frac[k].reserve(nToys);

    int nFailed = 0;
    std::cout << "Running " << nToys << " bootstrap toys...\n";

    for (int t = 0; t < nToys; t++)
    {
        // Poisson fluctuation of each bin
        for (int i = 1; i <= nBins; i++)
        {
            const double mu = h_exp->GetBinContent(i);
            const double k = rng.Poisson(mu);
            h_toy->SetBinContent(i, k);
        }

        // Re-create the RooDataHist for the toy
        RooDataHist d_toy("d_toy", "toy", RooArgList(x), h_toy);

        // Reset fit parameters to the nominal values
        for (int k = 0; k < nC - 1; k++)
            frac[k]->setVal(frac_nom[k]);

        std::unique_ptr<RooFitResult> r(
            model.fitTo(d_toy, Save(true), PrintLevel(-1), Warnings(false)));

        if (!r || r->status() != 0)
        {
            nFailed++;
            continue;
        }

        double sum_t = 0.0;
        for (int k = 0; k < nC - 1; k++)
        {
            const double fv = frac[k]->getVal();
            v_frac[k].push_back(fv);
            hF[k]->Fill(fv);
            sum_t += fv;
        }
        const double flast = 1.0 - sum_t;
        v_frac[nC - 1].push_back(flast);
        hF[nC - 1]->Fill(flast);

        if ((t + 1) % 100 == 0)
            std::cout << "  toy " << t + 1 << "/" << nToys << "\n";
    }
    std::cout << "Failed toys: " << nFailed << "/" << nToys << "\n\n";

    // =====================================================================
    // Bootstrap statistics
    // =====================================================================
    std::vector<double> sigF_bs(nC);
    std::cout << "=== Bootstrap results (" << v_frac[0].size() << " toys) ===\n";
    std::cout << std::fixed;
    std::cout.precision(4);
    for (int k = 0; k < nC; k++)
    {
        auto [mean, sig, p16, p50, p84] = toyStats(v_frac[k]);
        sigF_bs[k] = sig;
        std::cout << "f(" << labels[k] << ") : mean = " << mean
                  << "  sigma = " << sig
                  << "   median = " << p50
                  << "  [68% CI: " << p16 << ", " << p84 << "]\n";
    }
    std::cout << "=========================================\n\n";

    std::cout << "=== Comparison: Minos vs bootstrap sigma ===\n";
    for (int k = 0; k < nC; k++)
        std::cout << "f(" << labels[k] << ") : nominal = " << frac_nom[k]
                  << "   Minos = " << frac_err[k]
                  << "   bootstrap = " << sigF_bs[k] << "\n";
    std::cout << "===========================================\n\n";

    // =====================================================================
    // Plots
    // =====================================================================

    // ---- Canvas 1: fit result with RooFit's native plotting ----
    TCanvas *c1 = new TCanvas("c1", "RooFit template fit", 800, 600);

    RooPlot *xframe = x.frame(Title("Template fit (RooFit)"));
    data.plotOn(xframe, Name("data"));
    model.plotOn(xframe, LineColor(kRed), LineWidth(2), Name("total"));

    // Individual components drawn dashed
    for (int k = 0; k < nC; k++)
    {
        model.plotOn(xframe,
                     Components(*pdf_tmpl[k]),
                     LineStyle(kDashed),
                     LineColor(colors[k % colors.size()]),
                     LineWidth(2),
                     Name(Form("comp_%d", k)));
    }
    xframe->Draw();

    auto *leg = new TLegend(0.55, 0.60, 0.89, 0.89);
    leg->AddEntry(xframe->findObject("data"), "Exp", "ep");
    leg->AddEntry(xframe->findObject("total"), "Fit total", "l");
    for (int k = 0; k < nC; k++)
        leg->AddEntry(xframe->findObject(Form("comp_%d", k)),
                      Form("%s (%.1f #pm %.1f)%%",
                           labels[k].c_str(),
                           frac_nom[k] * 100, sigF_bs[k] * 100),
                      "l");
    leg->Draw();

    // ---- Canvas 2: bootstrap distributions of each fraction ----
    TCanvas *c2 = new TCanvas("c2", "Bootstrap distributions", 400 * nC, 400);
    c2->Divide(nC, 1);
    for (int k = 0; k < nC; k++)
    {
        c2->cd(k + 1);
        hF[k]->SetLineColor(colors[k % colors.size()]);
        hF[k]->SetFillColorAlpha(colors[k % colors.size()], 0.3);
        hF[k]->Draw();

        TLine *lnom = new TLine(frac_nom[k], 0, frac_nom[k], hF[k]->GetMaximum());
        lnom->SetLineColor(kRed);
        lnom->SetLineStyle(2);
        lnom->SetLineWidth(2);
        lnom->Draw();
    }
}