#include <TMinuit.h>
#include <TRandom3.h>
#include <TH1F.h>
#include <TCanvas.h>
#include <TLine.h>
#include <TLegend.h>
#include <TFile.h>
#include <TTree.h>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <iostream>
#include <iomanip>
#include <algorithm>
#include <tuple>
#include <cmath>

// =============================================================================
// 4-component extension of testFitMomdis.C
//
// Goal: validate whether bin-wise Poisson bootstrap is a faithful estimator
// of the statistical uncertainty by comparing it against the gold standard:
// event-by-event Monte-Carlo regeneration from the TRUE underlying mixture.
//
// Components: 1d_{5/2} , 2s_{1/2} , 1p_{1/2} , 1p_{3/2}
// =============================================================================

struct MomentaDist
{
    TH1F *Qz = nullptr;
    TH1F *Qt = nullptr;
    TH1F *Qy = nullptr;
    TH1F *Q = nullptr;
};

// ---------------------------------------------------------------------------
// Globals consumed by the TMinuit FCN. We keep the same pattern as the
// original (a single TH1F* per template) but generalised to N components.
// ---------------------------------------------------------------------------
namespace FitGlobals
{
    TH1F *gExp = nullptr;
    std::vector<TH1F *> gT; // ordered template histograms (size = N components)
}

// ---------------------------------------------------------------------------
// I/O: read a 2-column txt (q, dsigma/dq) into a TH1F
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

    const int nBins = (int)Qi.size();
    if (nBins < 2)
    {
        std::cerr << "ERROR: " << txtFile << " has fewer than 2 points\n";
        return momdis;
    }
    const double maxBin = Qi[nBins - 1] + (Qi[1] - Qi[0]) / 2.;

    momdis.Qt = new TH1F((std::string("Qt") + histName).c_str(),
                         "Qt", nBins, -maxBin, maxBin);
    for (int i = 0; i < nBins; i++)
        momdis.Qt->SetBinContent(i + 1, Qt[i]);
    return momdis;
}

// ---------------------------------------------------------------------------
// Rebin a fine theory histogram into the fit binning, normalising to unit area
// ---------------------------------------------------------------------------
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

// ---------------------------------------------------------------------------
// Recursive <-> physical fraction conversion
//   phi_0     = a_0
//   phi_1     = (1-a_0)*a_1
//   phi_2     = (1-a_0)*(1-a_1)*a_2
//   phi_{N-1} = (1-a_0)*(1-a_1)*...*(1-a_{N-2})
// All phi's automatically lie in [0,1] and sum to 1 if all a_k in [0,1].
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

std::vector<double> physicalToRecursive(const std::vector<double> &phi)
{
    const int N = (int)phi.size();
    std::vector<double> a(N - 1);
    double remain = 1.0;
    for (int k = 0; k < N - 1; k++)
    {
        if (remain <= 0)
        {
            a[k] = 0.0;
            continue;
        }
        a[k] = phi[k] / remain;
        if (a[k] < 0.0)
            a[k] = 0.0;
        if (a[k] > 1.0)
            a[k] = 1.0;
        remain *= (1.0 - a[k]);
    }
    return a;
}

// ---------------------------------------------------------------------------
// Generate fake data from the binned templates (true mixture + Poisson noise)
//   Used for: nominal "fake experimental" histogram and bootstrap source.
// ---------------------------------------------------------------------------
TH1F *generateFakeDataN(const std::vector<TH1F *> &tmpls,
                        const std::vector<double> &phi_true,
                        double N, std::string name, TRandom3 *rng)
{
    const int nBins = tmpls[0]->GetNbinsX();
    const double xmin = tmpls[0]->GetXaxis()->GetXmin();
    const double xmax = tmpls[0]->GetXaxis()->GetXmax();

    TH1F *hFake = new TH1F(name.c_str(), name.c_str(), nBins, xmin, xmax);

    for (int i = 1; i <= nBins; i++)
    {
        double mu = 0.0;
        for (size_t k = 0; k < tmpls.size(); k++)
            mu += phi_true[k] * tmpls[k]->GetBinContent(i);
        mu *= N;
        const double kcount = rng->Poisson(mu);
        hFake->SetBinContent(i, kcount);
        hFake->SetBinError(i, (kcount > 0) ? std::sqrt(kcount) : 1.0);
    }

    return hFake;
}

// ---------------------------------------------------------------------------
// Event-by-event toy MC: sample events from the true mixture using the
// FINE-binned theoretical histograms (before rebinning into fit binning).
// This is the gold-standard reference for the statistical fluctuation.
// ---------------------------------------------------------------------------
TH1F *generateToyEventByEventN(const std::vector<TH1F *> &theos,
                               const std::vector<double> &phi_true,
                               double N,
                               int nBinsFit, double xmin, double xmax,
                               std::string name, TRandom3 *rng)
{
    TH1F *hToy = new TH1F(name.c_str(), name.c_str(), nBinsFit, xmin, xmax);

    // Cumulative for component selection
    const int nC = (int)theos.size();
    std::vector<double> cum(nC);
    cum[0] = phi_true[0];
    for (int k = 1; k < nC; k++)
        cum[k] = cum[k - 1] + phi_true[k];

    const int nEvents = rng->Poisson(N);

    for (int ev = 0; ev < nEvents; ev++)
    {
        const double u = rng->Uniform();
        int kSel = nC - 1;
        for (int k = 0; k < nC; k++)
        {
            if (u < cum[k])
            {
                kSel = k;
                break;
            }
        }
        const double x = theos[kSel]->GetRandom();
        hToy->Fill(x);
    }

    for (int i = 1; i <= nBinsFit; i++)
    {
        const double k = hToy->GetBinContent(i);
        hToy->SetBinError(i, (k > 0) ? std::sqrt(k) : 1.0);
    }
    return hToy;
}

// ---------------------------------------------------------------------------
// Chi2 FCN for TMinuit. Parameters (recursive fractions):
//   par[0]      = N        (total counts)
//   par[1..N-1] = a_0..a_{N-2}    (recursive fractions)
// ---------------------------------------------------------------------------
void chi2FunctionN(Int_t & /*npar*/, Double_t * /*gin*/, Double_t &f,
                   Double_t *par, Int_t /*iflag*/)
{
    const int nC = (int)FitGlobals::gT.size();
    const double N = par[0];

    std::vector<double> a(nC - 1);
    for (int k = 0; k < nC - 1; k++)
        a[k] = par[1 + k];
    std::vector<double> phi = recursiveToPhysical(a);

    double chi2 = 0.0;
    const int nBins = FitGlobals::gExp->GetNbinsX();
    for (int i = 1; i <= nBins; i++)
    {
        const double data = FitGlobals::gExp->GetBinContent(i);
        const double err = FitGlobals::gExp->GetBinError(i);
        if (err <= 0)
            continue;

        double model = 0.0;
        for (int k = 0; k < nC; k++)
            model += phi[k] * FitGlobals::gT[k]->GetBinContent(i);
        model *= N;

        const double diff = data - model;
        chi2 += (diff * diff) / (err * err);
    }
    f = chi2;
}

struct FitResultN
{
    double N;
    std::vector<double> phi; // physical fractions
    double chi2;
    bool converged;
};

FitResultN doFitN(double N_init, const std::vector<double> &phi_init, bool verbose = false)
{
    const int nC = (int)FitGlobals::gT.size();
    const int nPar = 1 + (nC - 1); // N + (nC-1) recursive fractions

    TMinuit *minuit = new TMinuit(nPar);
    minuit->SetFCN(chi2FunctionN);
    minuit->SetPrintLevel(verbose ? 1 : -1);

    double arglist[2];
    int ierflg = 0;
    arglist[0] = 1.0;
    minuit->mnexcm("SET ERR", arglist, 1, ierflg);

    // Convert seed physical fractions to recursive
    std::vector<double> a_init = physicalToRecursive(phi_init);

    minuit->mnparm(0, "N", N_init, std::max(N_init * 0.01, 1.0),
                   0.0, 10.0 * N_init, ierflg);
    for (int k = 0; k < nC - 1; k++)
        minuit->mnparm(1 + k, Form("a%d", k), a_init[k], 0.05, 0.0, 1.0, ierflg);

    arglist[0] = 5000;
    arglist[1] = 0.01;
    minuit->mnexcm("MIGRAD", arglist, 2, ierflg);

    FitResultN res;
    res.phi.assign(nC, 0.0);
    double dummy;

    minuit->GetParameter(0, res.N, dummy);
    std::vector<double> a(nC - 1);
    for (int k = 0; k < nC - 1; k++)
        minuit->GetParameter(1 + k, a[k], dummy);
    res.phi = recursiveToPhysical(a);

    double edm, errdef;
    int nvpar, nparx, icstat;
    minuit->mnstat(res.chi2, edm, errdef, nvpar, nparx, icstat);
    res.converged = (ierflg == 0);

    delete minuit;
    return res;
}

// ---------------------------------------------------------------------------
// Stats helper
// ---------------------------------------------------------------------------
std::tuple<double, double, double, double, double> toyStats(std::vector<double> v)
{
    if (v.empty())
        return {0, 0, 0, 0, 0};
    std::sort(v.begin(), v.end());
    const int n = (int)v.size();
    double mean = 0.0;
    for (double x : v)
        mean += x;
    mean /= n;
    double var = 0.0;
    for (double x : v)
        var += (x - mean) * (x - mean);
    const double sigma = (n > 1) ? std::sqrt(var / (n - 1)) : 0.0;
    const double p16 = v[(int)(0.16 * n)];
    const double p50 = v[(int)(0.50 * n)];
    const double p84 = v[(int)(0.84 * n)];
    return std::make_tuple(mean, sigma, p16, p50, p84);
}

// =============================================================================
// MAIN
// =============================================================================
void testFitMomdis()
{
    // -------------------------------------------------------------------------
    // CONFIGURATION
    // -------------------------------------------------------------------------
    // Truth used to generate the fake data. Make sure they sum to 1.
    // Order: 1d_{5/2} , 2s_{1/2} , 1p_{1/2} , 1p_{3/2}
    std::vector<std::string> labels = {"1p_{3/2}", "1p_{1/2}", "1s_{1/2}", "1d_{5/2}"};
    std::vector<double> phi_true = {0.5, 0.5, 0.0, 0.0};

    // sanity check
    {
        double s = 0.0;
        for (double v : phi_true)
            s += v;
        if (std::abs(s - 1.0) > 1e-6)
        {
            std::cerr << "ERROR: phi_true does not sum to 1, sum = " << s << "\n";
            return;
        }
    }

    const double N_true = 6000;
    const int nToys = 1000;
    const int fitNBins = 50;
    const double fitMax = 300.0;

    const int nC = (int)labels.size();

    std::cout << "\n========== TEST CONFIG (" << nC << " components) ==========\n";
    for (int k = 0; k < nC; k++)
        std::cout << "phi_true(" << labels[k] << ") = " << phi_true[k] << "\n";
    std::cout << "N_true = " << N_true << "\n";
    std::cout << "nToys  = " << nToys << "\n";
    std::cout << "================================================\n\n";

    // -------------------------------------------------------------------------
    // Load theoretical templates
    // -------------------------------------------------------------------------
    std::vector<std::string> inFilesTheo = {
        "/nucl_lustre/pablogrusell/g249/g249_analysis/theory/JT/25F/sigt_p32-gs.txt",
        "/nucl_lustre/pablogrusell/g249/g249_analysis/theory/JT/25F/sigt_p12-gs.txt",
        "/nucl_lustre/pablogrusell/g249/g249_analysis/theory/JT/25F/sigt_s12-gs.txt",
        "/nucl_lustre/pablogrusell/g249/g249_analysis/theory/JT/25F/sigt_d52-gs.txt",
    };

    std::vector<MomentaDist> theos(nC);
    std::vector<TH1F *> tmpls(nC);
    for (int k = 0; k < nC; k++)
    {
        theos[k] = getMomentaDistFromtxt(inFilesTheo[k], Form("theo_%d", k));
        tmpls[k] = buildTemplate(theos[k].Qt, fitNBins, -fitMax, fitMax,
                                 Form("tmpl_%d", k));
    }

    FitGlobals::gT = tmpls;

    // -------------------------------------------------------------------------
    // Generate one fake "experimental" histogram from the true mixture
    // -------------------------------------------------------------------------
    TRandom3 rng(0);
    TH1F *hFake = generateFakeDataN(tmpls, phi_true, N_true, "hFake", &rng);

    std::cout << "Generated fake data: entries = " << hFake->GetEntries()
              << "  integral = " << hFake->Integral() << "\n\n";

    // -------------------------------------------------------------------------
    // Nominal fit on fake data (seeded at flat physical fractions)
    // -------------------------------------------------------------------------
    FitGlobals::gExp = hFake;

    std::vector<double> phi_seed(nC, 1.0 / nC);
    FitResultN nominal = doFitN(hFake->Integral(), phi_seed, true);

    std::cout << "\n=== Nominal fit on fake data ===\n";
    std::cout << "N    = " << nominal.N << "   (truth: " << N_true << ")\n";
    for (int k = 0; k < nC; k++)
        std::cout << "phi(" << labels[k] << ") = " << nominal.phi[k]
                  << "   (truth: " << phi_true[k] << ")\n";
    std::cout << "chi2 = " << nominal.chi2 << "  (ndf = "
              << fitNBins - (1 + (nC - 1)) << ")\n";
    std::cout << "================================\n\n";

    // -------------------------------------------------------------------------
    // METHOD 1 - BOOTSTRAP: Poisson fluctuation of the nominal fake histogram
    // -------------------------------------------------------------------------
    TH1F *hFakeToy = (TH1F *)hFake->Clone("hFakeToy");

    std::vector<std::vector<double>> v_phi_bs(nC);
    std::vector<double> v_N_bs;
    v_N_bs.reserve(nToys);
    for (int k = 0; k < nC; k++)
        v_phi_bs[k].reserve(nToys);

    int nFailedBS = 0;
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
        FitGlobals::gExp = hFakeToy;

        FitResultN r = doFitN(nominal.N, nominal.phi, false);
        if (!r.converged)
        {
            nFailedBS++;
            continue;
        }
        v_N_bs.push_back(r.N);
        for (int k = 0; k < nC; k++)
            v_phi_bs[k].push_back(r.phi[k]);

        if ((t + 1) % 100 == 0)
            std::cout << "  bootstrap toy " << t + 1 << "/" << nToys << "\n";
    }
    std::cout << "Failed bootstrap toys: " << nFailedBS << "/" << nToys << "\n\n";

    // -------------------------------------------------------------------------
    // METHOD 2 - TOY-MC: regenerate events from the TRUE mixture
    // -------------------------------------------------------------------------
    std::vector<TH1F *> theosFine(nC);
    for (int k = 0; k < nC; k++)
        theosFine[k] = theos[k].Qt;

    std::vector<std::vector<double>> v_phi_mc(nC);
    std::vector<double> v_N_mc;
    v_N_mc.reserve(nToys);
    for (int k = 0; k < nC; k++)
        v_phi_mc[k].reserve(nToys);

    int nFailedMC = 0;

    std::cout << "Running " << nToys << " TOY-MC toys (event-by-event from truth)...\n";

    for (int t = 0; t < nToys; t++)
    {
        TH1F *hToy = generateToyEventByEventN(theosFine, phi_true, N_true,
                                              fitNBins, -fitMax, fitMax,
                                              Form("hToyMC_%d", t), &rng);
        FitGlobals::gExp = hToy;

        const double Ninit = (hToy->Integral() > 0) ? hToy->Integral() : N_true;
        FitResultN r = doFitN(Ninit, nominal.phi, false);

        if (!r.converged)
        {
            nFailedMC++;
            delete hToy;
            continue;
        }

        v_N_mc.push_back(r.N);
        for (int k = 0; k < nC; k++)
            v_phi_mc[k].push_back(r.phi[k]);

        delete hToy;

        if ((t + 1) % 100 == 0)
            std::cout << "  toy-MC " << t + 1 << "/" << nToys << "\n";
    }
    std::cout << "Failed toy-MC toys: " << nFailedMC << "/" << nToys << "\n\n";

    FitGlobals::gExp = hFake;

    // -------------------------------------------------------------------------
    // Summary statistics
    // -------------------------------------------------------------------------
    std::cout << std::fixed << std::setprecision(4);

    auto [meanN_bs, sigN_bs, p16N_bs, p50N_bs, p84N_bs] = toyStats(v_N_bs);
    auto [meanN_mc, sigN_mc, p16N_mc, p50N_mc, p84N_mc] = toyStats(v_N_mc);

    std::cout << "=== BOOTSTRAP results (" << v_N_bs.size() << " toys) ===\n";
    std::cout << "N      : mean = " << meanN_bs << "  sigma = " << sigN_bs
              << "  median = " << p50N_bs
              << "  [68% CI: " << p16N_bs << ", " << p84N_bs << "]"
              << "  truth = " << N_true << "\n";
    for (int k = 0; k < nC; k++)
    {
        auto [m, s, p16, p50, p84] = toyStats(v_phi_bs[k]);
        std::cout << "phi(" << labels[k] << ") : mean = " << m
                  << "  sigma = " << s
                  << "  median = " << p50
                  << "  [68% CI: " << p16 << ", " << p84 << "]"
                  << "  truth = " << phi_true[k] << "\n";
    }
    std::cout << "==========================================\n\n";

    std::cout << "=== TOY-MC results (" << v_N_mc.size() << " toys) ===\n";
    std::cout << "N      : mean = " << meanN_mc << "  sigma = " << sigN_mc
              << "  median = " << p50N_mc
              << "  [68% CI: " << p16N_mc << ", " << p84N_mc << "]"
              << "  truth = " << N_true << "\n";
    for (int k = 0; k < nC; k++)
    {
        auto [m, s, p16, p50, p84] = toyStats(v_phi_mc[k]);
        std::cout << "phi(" << labels[k] << ") : mean = " << m
                  << "  sigma = " << s
                  << "  median = " << p50
                  << "  [68% CI: " << p16 << ", " << p84 << "]"
                  << "  truth = " << phi_true[k] << "\n";
    }
    std::cout << "==========================================\n\n";

    // -------------------------------------------------------------------------
    // Sigma comparison: bootstrap vs toy-MC. THIS is the headline test.
    // ratio MC/BS > 1  -> bootstrap UNDER-estimates the true uncertainty.
    // ratio MC/BS = 1  -> bootstrap is faithful.
    // ratio MC/BS < 1  -> bootstrap OVER-estimates the uncertainty.
    // -------------------------------------------------------------------------
    std::cout << "=== SIGMA COMPARISON (faithfulness test) ===\n";
    std::cout << "                  bootstrap     toy-MC      ratio MC/BS\n";

    auto fmt = [](double x)
    { std::ostringstream o;
                              o << std::fixed << std::setprecision(4) << x;
                              return o.str(); };

    std::cout << "N          :   " << std::setw(10) << fmt(sigN_bs)
              << "   " << std::setw(10) << fmt(sigN_mc)
              << "   " << std::setw(10)
              << fmt((sigN_bs > 0) ? sigN_mc / sigN_bs : 0.0) << "\n";
    for (int k = 0; k < nC; k++)
    {
        auto [m1, s1, a1, b1, c1] = toyStats(v_phi_bs[k]);
        auto [m2, s2, a2, b2, c2] = toyStats(v_phi_mc[k]);
        std::cout << "phi(" << labels[k] << ")  :   "
                  << std::setw(10) << fmt(s1)
                  << "   " << std::setw(10) << fmt(s2)
                  << "   " << std::setw(10)
                  << fmt((s1 > 0) ? s2 / s1 : 0.0) << "\n";
    }
    std::cout << "===========================================\n\n";

    // Pulls of the nominal fit against truth
    std::cout << "=== Pulls (nominal - truth) / sigma ===\n";
    std::cout << "                  bootstrap     toy-MC\n";
    std::cout << "N          :   " << std::setw(10)
              << fmt((sigN_bs > 0) ? (nominal.N - N_true) / sigN_bs : 0.0)
              << "   " << std::setw(10)
              << fmt((sigN_mc > 0) ? (nominal.N - N_true) / sigN_mc : 0.0) << "\n";
    for (int k = 0; k < nC; k++)
    {
        auto [m1, s1, a1, b1, c1] = toyStats(v_phi_bs[k]);
        auto [m2, s2, a2, b2, c2] = toyStats(v_phi_mc[k]);
        std::cout << "phi(" << labels[k] << ")  :   "
                  << std::setw(10)
                  << fmt((s1 > 0) ? (nominal.phi[k] - phi_true[k]) / s1 : 0.0)
                  << "   " << std::setw(10)
                  << fmt((s2 > 0) ? (nominal.phi[k] - phi_true[k]) / s2 : 0.0)
                  << "\n";
    }
    std::cout << "=======================================\n\n";

    // -------------------------------------------------------------------------
    // PLOTS
    // -------------------------------------------------------------------------
    std::vector<int> colors = {kBlue, kGreen + 2, kMagenta, kCyan + 1};

    // ---- Canvas 1: nominal fit on fake data
    TCanvas *c1 = new TCanvas("c1", "Nominal fit", 900, 650);

    TH1F *total = (TH1F *)tmpls[0]->Clone("fit_total");
    for (int i = 1; i <= total->GetNbinsX(); i++)
        total->SetBinContent(i, 0);
    total->SetBinError(0, 0);

    std::vector<TH1F *> comps(nC);
    for (int k = 0; k < nC; k++)
    {
        comps[k] = (TH1F *)tmpls[k]->Clone(Form("comp_%d", k));
        comps[k]->Scale(nominal.N * nominal.phi[k]);
        comps[k]->SetLineColor(colors[k]);
        comps[k]->SetLineStyle(2);
        comps[k]->SetLineWidth(2);
        for (int i = 1; i <= total->GetNbinsX(); i++)
            total->SetBinContent(i,
                                 total->GetBinContent(i) + comps[k]->GetBinContent(i));
    }
    total->SetLineColor(kRed);
    total->SetLineWidth(2);

    hFake->SetTitle("Fake data vs fit (4 components)");
    hFake->Draw("Ep");
    total->Draw("same hist");
    for (int k = 0; k < nC; k++)
        comps[k]->Draw("same hist");

    auto *leg = new TLegend(0.55, 0.60, 0.89, 0.89);
    leg->AddEntry(hFake, "Fake data", "ep");
    leg->AddEntry(total, "Fit total", "l");
    for (int k = 0; k < nC; k++)
    {
        auto [m, s, p16, p50, p84] = toyStats(v_phi_bs[k]);
        leg->AddEntry(comps[k],
                      Form("%s (%.1f #pm %.1f)%%   truth %.1f%%",
                           labels[k].c_str(),
                           nominal.phi[k] * 100, s * 100, phi_true[k] * 100),
                      "l");
    }
    leg->Draw();

    // ---- Canvas 2: bootstrap distributions of each fraction
    TCanvas *c2 = new TCanvas("c2", "Bootstrap distributions", 300 * nC, 350);
    c2->Divide(nC, 1);

    std::vector<TH1F *> hPhiBS(nC);
    for (int k = 0; k < nC; k++)
    {
        hPhiBS[k] = new TH1F(Form("hPhiBS_%d", k),
                             Form("Bootstrap #phi(%s);#phi;toys", labels[k].c_str()),
                             80, 0.0, 1.0);
        for (double v : v_phi_bs[k])
            hPhiBS[k]->Fill(v);

        c2->cd(k + 1);
        hPhiBS[k]->SetLineColor(colors[k]);
        hPhiBS[k]->SetFillColorAlpha(colors[k], 0.3);
        hPhiBS[k]->Draw();

        TLine *lt = new TLine(phi_true[k], 0, phi_true[k], hPhiBS[k]->GetMaximum());
        lt->SetLineColor(kRed);
        lt->SetLineStyle(2);
        lt->SetLineWidth(2);
        lt->Draw();
    }

    // ---- Canvas 3: toy-MC distributions
    TCanvas *c3 = new TCanvas("c3", "Toy-MC distributions", 300 * nC, 350);
    c3->Divide(nC, 1);

    std::vector<TH1F *> hPhiMC(nC);
    for (int k = 0; k < nC; k++)
    {
        hPhiMC[k] = new TH1F(Form("hPhiMC_%d", k),
                             Form("Toy-MC #phi(%s);#phi;toys", labels[k].c_str()),
                             80, 0.0, 1.0);
        for (double v : v_phi_mc[k])
            hPhiMC[k]->Fill(v);

        c3->cd(k + 1);
        hPhiMC[k]->SetLineColor(colors[k]);
        hPhiMC[k]->SetFillColorAlpha(colors[k], 0.3);
        hPhiMC[k]->Draw();

        TLine *lt = new TLine(phi_true[k], 0, phi_true[k], hPhiMC[k]->GetMaximum());
        lt->SetLineColor(kRed);
        lt->SetLineStyle(2);
        lt->SetLineWidth(2);
        lt->Draw();
    }

    // ---- Canvas 4: overlay of bootstrap vs toy-MC for each fraction
    TCanvas *c4 = new TCanvas("c4", "Overlay BS vs MC", 300 * nC, 350);
    c4->Divide(nC, 1);

    for (int k = 0; k < nC; k++)
    {
        c4->cd(k + 1);

        TH1F *hBS = (TH1F *)hPhiBS[k]->Clone(Form("BSn_%d", k));
        TH1F *hMC = (TH1F *)hPhiMC[k]->Clone(Form("MCn_%d", k));
        if (hBS->Integral() > 0)
            hBS->Scale(1.0 / hBS->Integral());
        if (hMC->Integral() > 0)
            hMC->Scale(1.0 / hMC->Integral());

        hBS->SetLineColor(kBlack);
        hBS->SetFillColorAlpha(kGray, 0.4);
        hBS->SetLineWidth(2);
        hMC->SetLineColor(colors[k]);
        hMC->SetFillColorAlpha(colors[k], 0.3);
        hMC->SetLineWidth(2);

        const double ymax = 1.2 * std::max(hBS->GetMaximum(), hMC->GetMaximum());
        hBS->SetMaximum(ymax);
        hBS->SetTitle(Form("%s : BS vs MC", labels[k].c_str()));
        hBS->GetXaxis()->SetTitle(Form("#phi(%s)", labels[k].c_str()));
        hBS->GetYaxis()->SetTitle("normalised toys");

        hBS->Draw("hist");
        hMC->Draw("hist same");

        TLine *lt = new TLine(phi_true[k], 0, phi_true[k], ymax);
        lt->SetLineColor(kRed);
        lt->SetLineStyle(2);
        lt->SetLineWidth(2);
        lt->Draw();

        auto *lg = new TLegend(0.55, 0.65, 0.89, 0.89);
        lg->AddEntry(hBS, "Bootstrap", "f");
        lg->AddEntry(hMC, "Toy-MC", "f");
        lg->AddEntry(lt, "truth", "l");
        lg->Draw();
    }
}