#include <TMinuit.h>
#include <TRandom3.h>

struct MomentaDist
{
    TH1F *Qz = nullptr;
    TH1F *Qt = nullptr;
    TH1F *Qy = nullptr;
    TH1F *Q = nullptr;
};

namespace FitGlobals
{
    TH1F *gExp = nullptr;
    TH1F *gT1 = nullptr;
    TH1F *gT2 = nullptr;
}

MomentaDist getMomentaDistFromtxt(std::string txtFile, std::string histName)
{
    MomentaDist momdis;
    std::vector<double> Qi, Qt;
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

TH1F *buildTemplate(TH1F *hTheo, int nBins, double xmin, double xmax, std::string name)
{
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

// Generate a fake experimental histogram from the templates with
// known a, b fractions (a + b = 1) and total area N, with Poisson noise.
TH1F *generateFakeData(TH1F *tmpl1, TH1F *tmpl2, double a, double b, double N,
                       std::string name, TRandom3 *rng)
{
    const int nBins = tmpl1->GetNbinsX();
    const double xmin = tmpl1->GetXaxis()->GetXmin();
    const double xmax = tmpl1->GetXaxis()->GetXmax();

    TH1F *hFake = new TH1F(name.c_str(), name.c_str(), nBins, xmin, xmax);

    for (int i = 1; i <= nBins; i++)
    {
        // Expected mean in bin i
        const double mu = N * (a * tmpl1->GetBinContent(i) + b * tmpl2->GetBinContent(i));
        // Poisson draw
        const double k = rng->Poisson(mu);
        hFake->SetBinContent(i, k);
        hFake->SetBinError(i, (k > 0) ? std::sqrt(k) : 1.0);
    }

    return hFake;
}

// Event-by-event toy MC:
// Sample events from the TRUE mixture (a * hTheo1 + b * hTheo2) built from the
// ORIGINAL fine-binned theoretical distributions, and fill a histogram in the
// fit binning. The total number of events fluctuates Poisson around N.
TH1F *generateToyEventByEvent(TH1F *hTheo1, TH1F *hTheo2,
                              double a, double N,
                              int nBinsFit, double xmin, double xmax,
                              std::string name, TRandom3 *rng)
{
    TH1F *hToy = new TH1F(name.c_str(), name.c_str(), nBinsFit, xmin, xmax);

    // Total number of events fluctuates Poisson around N
    const int nEvents = rng->Poisson(N);

    for (int ev = 0; ev < nEvents; ev++)
    {
        // Choose component with probability a (s1/2) or 1-a (d5/2)
        TH1F *src = (rng->Uniform() < a) ? hTheo1 : hTheo2;
        // Sample x from the original (fine-binned) theoretical distribution
        const double x = src->GetRandom();
        hToy->Fill(x);
    }

    // Assign Poisson errors
    for (int i = 1; i <= nBinsFit; i++)
    {
        const double k = hToy->GetBinContent(i);
        hToy->SetBinError(i, (k > 0) ? std::sqrt(k) : 1.0);
    }

    return hToy;
}

void chi2Function(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag)
{
    const double N = par[0];
    const double f1 = par[1];
    const double f2 = 1.0 - f1;

    double chi2 = 0.0;
    const int nBins = FitGlobals::gExp->GetNbinsX();

    for (int i = 1; i <= nBins; i++)
    {
        const double data = FitGlobals::gExp->GetBinContent(i);
        const double err = FitGlobals::gExp->GetBinError(i);
        if (err <= 0)
            continue;

        const double model = N * (f1 * FitGlobals::gT1->GetBinContent(i) + f2 * FitGlobals::gT2->GetBinContent(i));
        const double diff = data - model;
        chi2 += (diff * diff) / (err * err);
    }
    f = chi2;
}

struct FitResult
{
    double N;
    double f1;
    double chi2;
    bool converged;
};

FitResult doFit(double N_init, double f1_init, bool verbose = false)
{
    TMinuit *minuit = new TMinuit(2);
    minuit->SetFCN(chi2Function);
    minuit->SetPrintLevel(verbose ? 1 : -1);

    double arglist[2];
    int ierflg = 0;
    arglist[0] = 1.0;
    minuit->mnexcm("SET ERR", arglist, 1, ierflg);

    minuit->mnparm(0, "N", N_init, N_init * 0.01, 0.0, 10.0 * N_init, ierflg);
    minuit->mnparm(1, "f1", f1_init, 0.05, 0.0, 1.0, ierflg);

    arglist[0] = 5000;
    arglist[1] = 0.01;
    minuit->mnexcm("MIGRAD", arglist, 2, ierflg);

    FitResult res;
    double dummy;
    minuit->GetParameter(0, res.N, dummy);
    minuit->GetParameter(1, res.f1, dummy);

    double edm, errdef;
    int nvpar, nparx, icstat;
    minuit->mnstat(res.chi2, edm, errdef, nvpar, nparx, icstat);
    res.converged = (ierflg == 0);

    delete minuit;
    return res;
}

void testFitMomdis()
{
    // =====================================================================
    // CONFIGURATION: true values used to generate fake data
    // =====================================================================
    const double a_true = 0.07; // fraction s1/2
    const double b_true = 0.93; // fraction d5/2  (must be 1 - a_true)
    const double N_true = 6000; // total counts

    const int fitNBins = 50;
    const double fitMax = 300.0;

    const int nToys = 1000;

    std::cout << "\n========== TEST CONFIG ==========\n";
    std::cout << "a_true (s1/2) = " << a_true << "\n";
    std::cout << "b_true (d5/2) = " << b_true << "\n";
    std::cout << "N_true        = " << N_true << "\n";
    std::cout << "=================================\n\n";

    // =====================================================================
    // Load theoretical templates (from txt) and build them in exp binning
    // =====================================================================
    std::vector<std::string> inFilesTheo = {
        "/nucl_lustre/pablogrusell/g249/g249_analysis/theory/JT/25F/sigt_s12-gs.txt",
        "/nucl_lustre/pablogrusell/g249/g249_analysis/theory/JT/25F/sigt_d52-gs.txt",
    };

    auto momdis_2s12_gs = getMomentaDistFromtxt(inFilesTheo[0], "2s12_gs");
    auto momdis_1d52_gs = getMomentaDistFromtxt(inFilesTheo[1], "1d52_gs");

    TH1F *tmpl_s12 = buildTemplate(momdis_2s12_gs.Qt, fitNBins, -fitMax, fitMax, "tmpl_s12");
    TH1F *tmpl_d52 = buildTemplate(momdis_1d52_gs.Qt, fitNBins, -fitMax, fitMax, "tmpl_d52");

    FitGlobals::gT1 = tmpl_s12;
    FitGlobals::gT2 = tmpl_d52;

    // =====================================================================
    // Generate fake experimental data
    // =====================================================================
    TRandom3 rng(0);
    TH1F *hFake = generateFakeData(tmpl_s12, tmpl_d52, a_true, b_true, N_true,
                                   "hFake", &rng);

    std::cout << "Generated fake data: entries = " << hFake->GetEntries()
              << "  integral = " << hFake->Integral() << "\n\n";

    // =====================================================================
    // Nominal fit on fake data
    // =====================================================================
    FitGlobals::gExp = hFake;

    const double N0 = hFake->Integral();
    FitResult nominal = doFit(N0, 0.5, true);

    std::cout << "\n=== Nominal fit on fake data ===\n";
    std::cout << "N     = " << nominal.N << "   (truth: " << N_true << ")\n";
    std::cout << "f1    = " << nominal.f1 << "   (truth: " << a_true << ")\n";
    std::cout << "f2    = " << 1.0 - nominal.f1 << "   (truth: " << b_true << ")\n";
    std::cout << "chi2  = " << nominal.chi2 << "  (ndf = " << fitNBins - 2 << ")\n";
    std::cout << "================================\n\n";

    // =====================================================================
    // METHOD 1 — Bootstrap: Poisson variations over the binned fake data
    // =====================================================================
    TH1F *hFakeToy = (TH1F *)hFake->Clone("hFakeToy");
    FitGlobals::gExp = hFakeToy;

    TH1F *hN = new TH1F("hN", "Bootstrap N;N;toys", 80,
                        nominal.N * 0.9, nominal.N * 1.1);
    TH1F *hF1 = new TH1F("hF1", "Bootstrap f(2s_{1/2});f_{1};toys", 80, 0.0, 1.0);
    TH1F *hF2 = new TH1F("hF2", "Bootstrap f(1d_{5/2});f_{2};toys", 80, 0.0, 1.0);

    std::vector<double> v_N, v_f1, v_f2;
    v_N.reserve(nToys);
    v_f1.reserve(nToys);
    v_f2.reserve(nToys);

    int nFailed = 0;
    const int nBins = hFake->GetNbinsX();

    std::cout << "Running " << nToys << " BOOTSTRAP toys (Poisson per bin)...\n";

    for (int t = 0; t < nToys; t++)
    {
        for (int i = 1; i <= nBins; i++)
        {
            const double mu = hFake->GetBinContent(i);
            const double k = rng.Poisson(mu);
            hFakeToy->SetBinContent(i, k);
            hFakeToy->SetBinError(i, (k > 0) ? std::sqrt(k) : 1.0);
        }

        FitResult r = doFit(nominal.N, nominal.f1, false);
        if (!r.converged)
        {
            nFailed++;
            continue;
        }

        const double f2 = 1.0 - r.f1;
        v_N.push_back(r.N);
        v_f1.push_back(r.f1);
        v_f2.push_back(f2);

        hN->Fill(r.N);
        hF1->Fill(r.f1);
        hF2->Fill(f2);

        if ((t + 1) % 100 == 0)
            std::cout << "  bootstrap toy " << t + 1 << "/" << nToys << "\n";
    }

    std::cout << "Failed bootstrap toys: " << nFailed << "/" << nToys << "\n\n";

    // =====================================================================
    // METHOD 2 — Event-by-event toy MC:
    //   Sample events from the TRUE mixture (a * s1/2 + b * d5/2), built from
    //   the ORIGINAL theoretical distributions (before rebinning), rebin into
    //   the fit histogram, and fit. Measures the pure statistical variance
    //   of the sampling process, independent of the bootstrap rebinning.
    // =====================================================================
    TH1F *hN_mc = new TH1F("hN_mc", "ToyMC N;N;toys", 80,
                           nominal.N * 0.9, nominal.N * 1.1);
    TH1F *hF1_mc = new TH1F("hF1_mc", "ToyMC f(2s_{1/2});f_{1};toys", 80, 0.0, 1.0);
    TH1F *hF2_mc = new TH1F("hF2_mc", "ToyMC f(1d_{5/2});f_{2};toys", 80, 0.0, 1.0);

    std::vector<double> v_N_mc, v_f1_mc, v_f2_mc;
    v_N_mc.reserve(nToys);
    v_f1_mc.reserve(nToys);
    v_f2_mc.reserve(nToys);

    int nFailed_mc = 0;

    std::cout << "Running " << nToys << " TOY-MC toys (event-by-event from true mixture)...\n";

    for (int t = 0; t < nToys; t++)
    {
        TH1F *hToy = generateToyEventByEvent(momdis_2s12_gs.Qt, momdis_1d52_gs.Qt,
                                             a_true, N_true,
                                             fitNBins, -fitMax, fitMax,
                                             Form("hToyMC_%d", t), &rng);
        FitGlobals::gExp = hToy;

        const double Ninit = (hToy->Integral() > 0) ? hToy->Integral() : N_true;
        FitResult r = doFit(Ninit, nominal.f1, false);

        if (!r.converged)
        {
            nFailed_mc++;
            delete hToy;
            continue;
        }

        const double f2 = 1.0 - r.f1;
        v_N_mc.push_back(r.N);
        v_f1_mc.push_back(r.f1);
        v_f2_mc.push_back(f2);

        hN_mc->Fill(r.N);
        hF1_mc->Fill(r.f1);
        hF2_mc->Fill(f2);

        delete hToy;

        if ((t + 1) % 100 == 0)
            std::cout << "  toy-MC " << t + 1 << "/" << nToys << "\n";
    }
    std::cout << "Failed toy-MC toys: " << nFailed_mc << "/" << nToys << "\n\n";

    // Restore nominal pointer
    FitGlobals::gExp = hFake;

    // =====================================================================
    // Statistics (both methods)
    // =====================================================================
    auto stats = [](std::vector<double> v)
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
    };

    auto [meanN, sigN, p16N, p50N, p84N] = stats(v_N);
    auto [meanF1, sigF1, p16F1, p50F1, p84F1] = stats(v_f1);
    auto [meanF2, sigF2, p16F2, p50F2, p84F2] = stats(v_f2);

    auto [meanN_mc, sigN_mc, p16N_mc, p50N_mc, p84N_mc] = stats(v_N_mc);
    auto [meanF1_mc, sigF1_mc, p16F1_mc, p50F1_mc, p84F1_mc] = stats(v_f1_mc);
    auto [meanF2_mc, sigF2_mc, p16F2_mc, p50F2_mc, p84F2_mc] = stats(v_f2_mc);

    std::cout << std::fixed;
    std::cout.precision(4);

    std::cout << "=== BOOTSTRAP results (" << v_N.size() << " toys) ===\n";
    std::cout << "N      : mean = " << meanN << "  sigma = " << sigN
              << "   median = " << p50N << "  [68% CI: " << p16N << ", " << p84N << "]   truth = " << N_true << "\n";
    std::cout << "f(s12) : mean = " << meanF1 << "  sigma = " << sigF1
              << "   median = " << p50F1 << "  [68% CI: " << p16F1 << ", " << p84F1 << "]   truth = " << a_true << "\n";
    std::cout << "f(d52) : mean = " << meanF2 << "  sigma = " << sigF2
              << "   median = " << p50F2 << "  [68% CI: " << p16F2 << ", " << p84F2 << "]   truth = " << b_true << "\n";
    std::cout << "=========================================\n\n";

    std::cout << "=== TOY-MC results (" << v_N_mc.size() << " toys) ===\n";
    std::cout << "N      : mean = " << meanN_mc << "  sigma = " << sigN_mc
              << "   median = " << p50N_mc << "  [68% CI: " << p16N_mc << ", " << p84N_mc << "]   truth = " << N_true << "\n";
    std::cout << "f(s12) : mean = " << meanF1_mc << "  sigma = " << sigF1_mc
              << "   median = " << p50F1_mc << "  [68% CI: " << p16F1_mc << ", " << p84F1_mc << "]   truth = " << a_true << "\n";
    std::cout << "f(d52) : mean = " << meanF2_mc << "  sigma = " << sigF2_mc
              << "   median = " << p50F2_mc << "  [68% CI: " << p16F2_mc << ", " << p84F2_mc << "]   truth = " << b_true << "\n";
    std::cout << "=========================================\n\n";

    std::cout << "=== SIGMA COMPARISON (bootstrap vs toy-MC) ===\n";
    std::cout << "            bootstrap        toy-MC        ratio(MC/BS)\n";
    std::cout << "N       :   " << sigN << "        " << sigN_mc << "        " << sigN_mc / sigN << "\n";
    std::cout << "f(s12)  :   " << sigF1 << "        " << sigF1_mc << "        " << sigF1_mc / sigF1 << "\n";
    std::cout << "f(d52)  :   " << sigF2 << "        " << sigF2_mc << "        " << sigF2_mc / sigF2 << "\n";
    std::cout << "===============================================\n\n";

    // Pulls: (nominal - truth) / sigma
    std::cout << "=== Pulls (nominal - truth) / sigma ===\n";
    std::cout << "            bootstrap         toy-MC\n";
    std::cout << "N       :   " << (nominal.N - N_true) / sigN
              << "          " << (nominal.N - N_true) / sigN_mc << "\n";
    std::cout << "f(s12)  :   " << (nominal.f1 - a_true) / sigF1
              << "          " << (nominal.f1 - a_true) / sigF1_mc << "\n";
    std::cout << "f(d52)  :   " << ((1.0 - nominal.f1) - b_true) / sigF2
              << "          " << ((1.0 - nominal.f1) - b_true) / sigF2_mc << "\n";
    std::cout << "=======================================\n\n";

    std::cout << "Final result:\n";
    std::cout << "  f(2s1/2) = " << nominal.f1 << " +/- " << sigF1 << " (BS) / " << sigF1_mc << " (MC)   truth: " << a_true << "\n";
    std::cout << "  f(1d5/2) = " << 1.0 - nominal.f1 << " +/- " << sigF2 << " (BS) / " << sigF2_mc << " (MC)   truth: " << b_true << "\n\n";

    // =====================================================================
    // Plots
    // =====================================================================
    TCanvas *c1 = new TCanvas("c1", "Test fit", 800, 600);

    TH1F *comp_s12 = (TH1F *)tmpl_s12->Clone("comp_s12");
    TH1F *comp_d52 = (TH1F *)tmpl_d52->Clone("comp_d52");
    TH1F *total = (TH1F *)tmpl_s12->Clone("fit_total");

    comp_s12->Scale(nominal.N * nominal.f1);
    comp_d52->Scale(nominal.N * (1.0 - nominal.f1));
    for (int i = 1; i <= total->GetNbinsX(); i++)
    {
        total->SetBinContent(i, comp_s12->GetBinContent(i) + comp_d52->GetBinContent(i));
        total->SetBinError(i, 0);
    }

    comp_s12->SetLineColor(kBlue);
    comp_s12->SetLineStyle(2);
    comp_s12->SetLineWidth(2);
    comp_d52->SetLineColor(kGreen + 2);
    comp_d52->SetLineStyle(2);
    comp_d52->SetLineWidth(2);
    total->SetLineColor(kRed);
    total->SetLineWidth(2);

    hFake->SetTitle(Form("Fake data (a=%.2f, b=%.2f, N=%.0f)", a_true, b_true, N_true));
    hFake->Draw("Ep");
    total->Draw("same hist");
    comp_s12->Draw("same hist");
    comp_d52->Draw("same hist");

    auto *leg = new TLegend(0.60, 0.65, 0.89, 0.89);
    leg->AddEntry(hFake, "Fake data", "ep");
    leg->AddEntry(total, "Fit total", "l");
    leg->AddEntry(comp_s12, Form("2s_{1/2} (%.1f #pm %.1f)%%", nominal.f1 * 100, sigF1 * 100), "l");
    leg->AddEntry(comp_d52, Form("1d_{5/2} (%.1f #pm %.1f)%%", (1.0 - nominal.f1) * 100, sigF2 * 100), "l");
    leg->Draw();

    // ---------- Canvas 2: bootstrap distributions ----------
    TCanvas *c2 = new TCanvas("c2", "Bootstrap distributions", 1200, 400);
    c2->Divide(3, 1);

    c2->cd(1);
    hN->SetLineColor(kBlack);
    hN->SetFillColorAlpha(kGray, 0.5);
    hN->Draw();
    TLine *lN = new TLine(N_true, 0, N_true, hN->GetMaximum());
    lN->SetLineColor(kRed);
    lN->SetLineStyle(2);
    lN->SetLineWidth(2);
    lN->Draw();

    c2->cd(2);
    hF1->SetLineColor(kBlue);
    hF1->SetFillColorAlpha(kBlue, 0.3);
    hF1->Draw();
    TLine *lF1 = new TLine(a_true, 0, a_true, hF1->GetMaximum());
    lF1->SetLineColor(kRed);
    lF1->SetLineStyle(2);
    lF1->SetLineWidth(2);
    lF1->Draw();

    c2->cd(3);
    hF2->SetLineColor(kGreen + 2);
    hF2->SetFillColorAlpha(kGreen + 2, 0.3);
    hF2->Draw();
    TLine *lF2 = new TLine(b_true, 0, b_true, hF2->GetMaximum());
    lF2->SetLineColor(kRed);
    lF2->SetLineStyle(2);
    lF2->SetLineWidth(2);
    lF2->Draw();

    // ---------- Canvas 3: toy-MC distributions ----------
    TCanvas *c3 = new TCanvas("c3", "Toy-MC distributions", 1200, 400);
    c3->Divide(3, 1);

    c3->cd(1);
    hN_mc->SetLineColor(kBlack);
    hN_mc->SetFillColorAlpha(kGray, 0.5);
    hN_mc->Draw();
    TLine *lN_mc = new TLine(N_true, 0, N_true, hN_mc->GetMaximum());
    lN_mc->SetLineColor(kRed);
    lN_mc->SetLineStyle(2);
    lN_mc->SetLineWidth(2);
    lN_mc->Draw();

    c3->cd(2);
    hF1_mc->SetLineColor(kBlue);
    hF1_mc->SetFillColorAlpha(kBlue, 0.3);
    hF1_mc->Draw();
    TLine *lF1_mc = new TLine(a_true, 0, a_true, hF1_mc->GetMaximum());
    lF1_mc->SetLineColor(kRed);
    lF1_mc->SetLineStyle(2);
    lF1_mc->SetLineWidth(2);
    lF1_mc->Draw();

    c3->cd(3);
    hF2_mc->SetLineColor(kGreen + 2);
    hF2_mc->SetFillColorAlpha(kGreen + 2, 0.3);
    hF2_mc->Draw();
    TLine *lF2_mc = new TLine(b_true, 0, b_true, hF2_mc->GetMaximum());
    lF2_mc->SetLineColor(kRed);
    lF2_mc->SetLineStyle(2);
    lF2_mc->SetLineWidth(2);
    lF2_mc->Draw();

    // ---------- Canvas 4: overlay bootstrap vs toy-MC ----------
    TCanvas *c4 = new TCanvas("c4", "Bootstrap vs Toy-MC overlay", 1200, 400);
    c4->Divide(3, 1);

    auto drawOverlay = [&](TH1F *hBS, TH1F *hMC, double truth, const char *xtitle)
    {
        // Normalise to unit area so shapes are directly comparable
        TH1F *hBSn = (TH1F *)hBS->Clone(Form("%s_norm", hBS->GetName()));
        TH1F *hMCn = (TH1F *)hMC->Clone(Form("%s_norm", hMC->GetName()));
        if (hBSn->Integral() > 0)
            hBSn->Scale(1.0 / hBSn->Integral());
        if (hMCn->Integral() > 0)
            hMCn->Scale(1.0 / hMCn->Integral());

        hBSn->SetLineColor(kBlack);
        hBSn->SetFillColorAlpha(kGray, 0.4);
        hMCn->SetLineColor(kAzure + 1);
        hMCn->SetFillColorAlpha(kAzure + 1, 0.3);
        hBSn->SetLineWidth(2);
        hMCn->SetLineWidth(2);

        const double ymax = 1.2 * std::max(hBSn->GetMaximum(), hMCn->GetMaximum());
        hBSn->SetMaximum(ymax);
        hBSn->GetXaxis()->SetTitle(xtitle);
        hBSn->GetYaxis()->SetTitle("normalised toys");
        hBSn->SetTitle("");
        hBSn->Draw("hist");
        hMCn->Draw("hist same");

        TLine *lt = new TLine(truth, 0, truth, ymax);
        lt->SetLineColor(kRed);
        lt->SetLineStyle(2);
        lt->SetLineWidth(2);
        lt->Draw();

        auto *lg = new TLegend(0.60, 0.70, 0.89, 0.89);
        lg->AddEntry(hBSn, "Bootstrap", "f");
        lg->AddEntry(hMCn, "Toy-MC", "f");
        lg->AddEntry(lt, "truth", "l");
        lg->Draw();
    };

    c4->cd(1);
    drawOverlay(hN, hN_mc, N_true, "N");
    c4->cd(2);
    drawOverlay(hF1, hF1_mc, a_true, "f(2s_{1/2})");
    c4->cd(3);
    drawOverlay(hF2, hF2_mc, b_true, "f(1d_{5/2})");
}