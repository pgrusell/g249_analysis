#include <TMinuit.h>
#include <TRandom3.h>
#include <TH1F.h>
#include <TCanvas.h>
#include <TLine.h>
#include <TLegend.h>
#include <TFile.h>
#include <TTree.h>
#include <TF1.h>
#include <TMath.h>
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

    // Symmetric (parabolic / HESSE-like) errors returned by Minuit
    double N_err;
    std::vector<double> phi_err;

    // Asymmetric MINOS errors (positive and negative branches)
    // For the recursive parameters first, then propagated to phi.
    double N_eplus;
    double N_eminus;
    std::vector<double> phi_eplus;
    std::vector<double> phi_eminus;

    double chi2;
    int ndf;
    bool converged;
};

// ---------------------------------------------------------------------------
// Propagate MINOS errors from the recursive parameters a_k to the physical
// fractions phi_k via numerical differentiation (Jacobian).
// We treat eplus and eminus as 1-sigma intervals and propagate them as if
// they were symmetric around the recursive parameter -- this is the standard
// linear propagation and is fine when the asymmetry is moderate.
// ---------------------------------------------------------------------------
static void propagateMinosToPhi(const std::vector<double> &a_central,
                                const std::vector<double> &a_eplus,
                                const std::vector<double> &a_eminus,
                                std::vector<double> &phi_eplus,
                                std::vector<double> &phi_eminus)
{
    const int nA = (int)a_central.size();
    const int nC = nA + 1;

    // Jacobian J[k][j] = d phi_k / d a_j   evaluated at a_central
    // Computed by forward finite differences.
    std::vector<std::vector<double>> J(nC, std::vector<double>(nA, 0.0));
    const std::vector<double> phi0 = recursiveToPhysical(a_central);

    for (int j = 0; j < nA; j++)
    {
        const double h = std::max(1e-6, 1e-4 * std::max(1.0, std::fabs(a_central[j])));
        std::vector<double> a_pert = a_central;
        a_pert[j] += h;
        if (a_pert[j] > 1.0)
            a_pert[j] = 1.0;
        if (a_pert[j] < 0.0)
            a_pert[j] = 0.0;
        const double dh = a_pert[j] - a_central[j];
        if (std::fabs(dh) < 1e-15)
            continue;
        const std::vector<double> phi_pert = recursiveToPhysical(a_pert);
        for (int k = 0; k < nC; k++)
            J[k][j] = (phi_pert[k] - phi0[k]) / dh;
    }

    // Linear (uncorrelated-in-quadrature) propagation. This is an
    // approximation: it ignores correlations between the a_k. For a more
    // rigorous treatment one would need Minuit's full covariance matrix.
    phi_eplus.assign(nC, 0.0);
    phi_eminus.assign(nC, 0.0);
    for (int k = 0; k < nC; k++)
    {
        double sp = 0.0, sm = 0.0;
        for (int j = 0; j < nA; j++)
        {
            // For each component of phi, the sign of dphi/da decides
            // which branch of the asymmetric error contributes to + and -.
            const double Jkj = J[k][j];
            const double ep = a_eplus[j];
            const double em = std::fabs(a_eminus[j]);
            if (Jkj >= 0)
            {
                sp += (Jkj * ep) * (Jkj * ep);
                sm += (Jkj * em) * (Jkj * em);
            }
            else
            {
                sp += (Jkj * em) * (Jkj * em);
                sm += (Jkj * ep) * (Jkj * ep);
            }
        }
        phi_eplus[k] = std::sqrt(sp);
        phi_eminus[k] = std::sqrt(sm);
    }
}

FitResultN doFitN(double N_init, const std::vector<double> &phi_init,
                  bool verbose = false, bool runMinos = false)
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
    res.phi_err.assign(nC, 0.0);
    res.phi_eplus.assign(nC, 0.0);
    res.phi_eminus.assign(nC, 0.0);
    res.N_eplus = 0.0;
    res.N_eminus = 0.0;

    double dummy;
    double parErr;

    // Central values + parabolic (HESSE) errors
    minuit->GetParameter(0, res.N, res.N_err);
    std::vector<double> a(nC - 1);
    std::vector<double> a_err(nC - 1);
    for (int k = 0; k < nC - 1; k++)
        minuit->GetParameter(1 + k, a[k], a_err[k]);
    res.phi = recursiveToPhysical(a);

    // Propagate parabolic errors on a_k to phi_k via the Jacobian (treating
    // a_err as symmetric eplus = eminus = a_err).
    {
        std::vector<double> a_eplus = a_err;
        std::vector<double> a_eminus = a_err;
        std::vector<double> phi_ep, phi_em;
        propagateMinosToPhi(a, a_eplus, a_eminus, phi_ep, phi_em);
        // For symmetric input the two propagated branches are equal -> use either.
        for (int k = 0; k < nC; k++)
            res.phi_err[k] = phi_ep[k];
    }

    // -----------------------------------------------------------------------
    // MINOS: asymmetric errors. We run it once, then read eplus/eminus per
    // parameter. Note that for fractions on the boundary [0,1] MINOS may
    // hit the limit and return the parabolic error instead.
    // -----------------------------------------------------------------------
    if (runMinos)
    {
        arglist[0] = 5000;
        // call MINOS for all parameters
        minuit->mnexcm("MINOS", arglist, 1, ierflg);

        double eplus, eminus, eparab, gcc;
        // Parameter 0 = N
        minuit->mnerrs(0, eplus, eminus, eparab, gcc);
        res.N_eplus = eplus;
        res.N_eminus = std::fabs(eminus);

        // Parameters 1..nC-1 = a_k
        std::vector<double> a_eplus(nC - 1, 0.0);
        std::vector<double> a_eminus(nC - 1, 0.0);
        for (int k = 0; k < nC - 1; k++)
        {
            minuit->mnerrs(1 + k, eplus, eminus, eparab, gcc);
            // Fall back to parabolic if MINOS returned 0 (limits hit etc.)
            a_eplus[k] = (eplus != 0.0) ? eplus : eparab;
            a_eminus[k] = (eminus != 0.0) ? std::fabs(eminus) : eparab;
        }

        // Propagate to physical fractions
        propagateMinosToPhi(a, a_eplus, a_eminus,
                            res.phi_eplus, res.phi_eminus);
    }

    double edm, errdef;
    int nvpar, nparx, icstat;
    minuit->mnstat(res.chi2, edm, errdef, nvpar, nparx, icstat);
    res.converged = (ierflg == 0);

    // ndf = nBins (with err>0) - nFreePars
    int nBinsUsed = 0;
    for (int i = 1; i <= FitGlobals::gExp->GetNbinsX(); i++)
        if (FitGlobals::gExp->GetBinError(i) > 0)
            nBinsUsed++;
    res.ndf = nBinsUsed - nPar;

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
    std::vector<double> phi_true = {0.5, 0.0, 0.5, 0.0};

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

    const double N_true = 600000;
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
    // Nominal fit on fake data (seeded at flat physical fractions) + MINOS
    // -------------------------------------------------------------------------
    FitGlobals::gExp = hFake;

    std::vector<double> phi_seed(nC, 1.0 / nC);
    FitResultN nominal = doFitN(hFake->Integral(), phi_seed,
                                /*verbose*/ true, /*runMinos*/ true);

    std::cout << "\n=== Nominal fit on fake data ===\n";
    std::cout << "N    = " << nominal.N
              << "  +/- " << nominal.N_err
              << "  (MINOS +" << nominal.N_eplus
              << " / -" << nominal.N_eminus << ")"
              << "   (truth: " << N_true << ")\n";
    for (int k = 0; k < nC; k++)
        std::cout << "phi(" << labels[k] << ") = " << nominal.phi[k]
                  << "  +/- " << nominal.phi_err[k]
                  << "  (MINOS +" << nominal.phi_eplus[k]
                  << " / -" << nominal.phi_eminus[k] << ")"
                  << "   (truth: " << phi_true[k] << ")\n";
    std::cout << "chi2 = " << nominal.chi2 << "  (ndf = " << nominal.ndf << ")\n";
    std::cout << "================================\n\n";

    // -------------------------------------------------------------------------
    // METHOD 1 - BOOTSTRAP: Poisson fluctuation of the nominal fake histogram
    // -------------------------------------------------------------------------
    TH1F *hFakeToy = (TH1F *)hFake->Clone("hFakeToy");

    std::vector<std::vector<double>> v_phi_bs(nC);
    std::vector<double> v_N_bs;
    std::vector<double> v_chi2_bs;
    std::vector<int> v_ndf_bs;
    // Also collect the MINOS errors per toy for diagnostic plots
    std::vector<std::vector<double>> v_phi_eplus_bs(nC);
    std::vector<std::vector<double>> v_phi_eminus_bs(nC);
    std::vector<double> v_N_eplus_bs;
    std::vector<double> v_N_eminus_bs;

    v_N_bs.reserve(nToys);
    v_chi2_bs.reserve(nToys);
    v_ndf_bs.reserve(nToys);
    v_N_eplus_bs.reserve(nToys);
    v_N_eminus_bs.reserve(nToys);
    for (int k = 0; k < nC; k++)
    {
        v_phi_bs[k].reserve(nToys);
        v_phi_eplus_bs[k].reserve(nToys);
        v_phi_eminus_bs[k].reserve(nToys);
    }

    int nFailedBS = 0;
    const int nBins = hFake->GetNbinsX();

    std::cout << "Running " << nToys << " BOOTSTRAP toys (Poisson per bin) with MINOS...\n";

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

        FitResultN r = doFitN(nominal.N, nominal.phi,
                              /*verbose*/ false, /*runMinos*/ true);
        if (!r.converged)
        {
            nFailedBS++;
            continue;
        }
        v_N_bs.push_back(r.N);
        v_chi2_bs.push_back(r.chi2);
        v_ndf_bs.push_back(r.ndf);
        v_N_eplus_bs.push_back(r.N_eplus);
        v_N_eminus_bs.push_back(r.N_eminus);
        for (int k = 0; k < nC; k++)
        {
            v_phi_bs[k].push_back(r.phi[k]);
            v_phi_eplus_bs[k].push_back(r.phi_eplus[k]);
            v_phi_eminus_bs[k].push_back(r.phi_eminus[k]);
        }

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
    std::vector<double> v_chi2_mc;
    std::vector<int> v_ndf_mc;
    std::vector<std::vector<double>> v_phi_eplus_mc(nC);
    std::vector<std::vector<double>> v_phi_eminus_mc(nC);
    std::vector<double> v_N_eplus_mc;
    std::vector<double> v_N_eminus_mc;

    v_N_mc.reserve(nToys);
    v_chi2_mc.reserve(nToys);
    v_ndf_mc.reserve(nToys);
    v_N_eplus_mc.reserve(nToys);
    v_N_eminus_mc.reserve(nToys);
    for (int k = 0; k < nC; k++)
    {
        v_phi_mc[k].reserve(nToys);
        v_phi_eplus_mc[k].reserve(nToys);
        v_phi_eminus_mc[k].reserve(nToys);
    }

    int nFailedMC = 0;

    std::cout << "Running " << nToys << " TOY-MC toys (event-by-event from truth) with MINOS...\n";

    for (int t = 0; t < nToys; t++)
    {
        TH1F *hToy = generateToyEventByEventN(theosFine, phi_true, N_true,
                                              fitNBins, -fitMax, fitMax,
                                              Form("hToyMC_%d", t), &rng);
        FitGlobals::gExp = hToy;

        const double Ninit = (hToy->Integral() > 0) ? hToy->Integral() : N_true;
        FitResultN r = doFitN(Ninit, nominal.phi,
                              /*verbose*/ false, /*runMinos*/ true);

        if (!r.converged)
        {
            nFailedMC++;
            delete hToy;
            continue;
        }

        v_N_mc.push_back(r.N);
        v_chi2_mc.push_back(r.chi2);
        v_ndf_mc.push_back(r.ndf);
        v_N_eplus_mc.push_back(r.N_eplus);
        v_N_eminus_mc.push_back(r.N_eminus);
        for (int k = 0; k < nC; k++)
        {
            v_phi_mc[k].push_back(r.phi[k]);
            v_phi_eplus_mc[k].push_back(r.phi_eplus[k]);
            v_phi_eminus_mc[k].push_back(r.phi_eminus[k]);
        }

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

    // -------------------------------------------------------------------------
    // Comparison: spread of the toy distribution vs MEAN MINOS error
    // (faithfulness of MINOS as a per-fit error estimate).
    // ratio sigma_toys / <minos> ~ 1  -> MINOS errors are correctly sized.
    // -------------------------------------------------------------------------
    auto meanOf = [](const std::vector<double> &v)
    {
        if (v.empty())
            return 0.0;
        double s = 0.0;
        for (double x : v)
            s += x;
        return s / v.size();
    };

    std::cout << "=== MINOS errors vs toy spread ===\n";
    std::cout << "                <e+ MINOS>   <e- MINOS>   sigma(toys)   ratio sig/<e_avg>\n";
    {
        const double eP = meanOf(v_N_eplus_bs);
        const double eM = meanOf(v_N_eminus_bs);
        const double eAvg = 0.5 * (eP + eM);
        std::cout << "[BS] N      :  " << std::setw(10) << fmt(eP)
                  << "   " << std::setw(10) << fmt(eM)
                  << "   " << std::setw(10) << fmt(sigN_bs)
                  << "   " << std::setw(10)
                  << fmt(eAvg > 0 ? sigN_bs / eAvg : 0.0) << "\n";
    }
    for (int k = 0; k < nC; k++)
    {
        auto [m, s, p16, p50, p84] = toyStats(v_phi_bs[k]);
        const double eP = meanOf(v_phi_eplus_bs[k]);
        const double eM = meanOf(v_phi_eminus_bs[k]);
        const double eAvg = 0.5 * (eP + eM);
        std::cout << "[BS] phi(" << labels[k] << "):  "
                  << std::setw(10) << fmt(eP)
                  << "   " << std::setw(10) << fmt(eM)
                  << "   " << std::setw(10) << fmt(s)
                  << "   " << std::setw(10)
                  << fmt(eAvg > 0 ? s / eAvg : 0.0) << "\n";
    }
    {
        const double eP = meanOf(v_N_eplus_mc);
        const double eM = meanOf(v_N_eminus_mc);
        const double eAvg = 0.5 * (eP + eM);
        std::cout << "[MC] N      :  " << std::setw(10) << fmt(eP)
                  << "   " << std::setw(10) << fmt(eM)
                  << "   " << std::setw(10) << fmt(sigN_mc)
                  << "   " << std::setw(10)
                  << fmt(eAvg > 0 ? sigN_mc / eAvg : 0.0) << "\n";
    }
    for (int k = 0; k < nC; k++)
    {
        auto [m, s, p16, p50, p84] = toyStats(v_phi_mc[k]);
        const double eP = meanOf(v_phi_eplus_mc[k]);
        const double eM = meanOf(v_phi_eminus_mc[k]);
        const double eAvg = 0.5 * (eP + eM);
        std::cout << "[MC] phi(" << labels[k] << "):  "
                  << std::setw(10) << fmt(eP)
                  << "   " << std::setw(10) << fmt(eM)
                  << "   " << std::setw(10) << fmt(s)
                  << "   " << std::setw(10)
                  << fmt(eAvg > 0 ? s / eAvg : 0.0) << "\n";
    }
    std::cout << "===================================\n\n";

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

    // -------------------------------------------------------------------------
    // ---- Canvas 5: chi2 distribution of the toy fits (BS and MC)
    // -------------------------------------------------------------------------
    // Use the median ndf of each population for the reference chi2 PDF.
    auto medianInt = [](std::vector<int> v) -> int
    {
        if (v.empty())
            return 0;
        std::sort(v.begin(), v.end());
        return v[v.size() / 2];
    };
    const int ndf_bs = medianInt(v_ndf_bs);
    const int ndf_mc = medianInt(v_ndf_mc);

    // Histogram range chosen with a sensible upper edge
    double chi2_max = 0.0;
    for (double x : v_chi2_bs)
        if (x > chi2_max)
            chi2_max = x;
    for (double x : v_chi2_mc)
        if (x > chi2_max)
            chi2_max = x;
    if (chi2_max <= 0)
        chi2_max = 2.0 * std::max(ndf_bs, ndf_mc);
    chi2_max *= 1.05;

    const int nBinsChi2 = 60;

    TH1F *hChi2BS = new TH1F("hChi2BS",
                             "Bootstrap #chi^{2};#chi^{2};toys",
                             nBinsChi2, 0.0, chi2_max);
    TH1F *hChi2MC = new TH1F("hChi2MC",
                             "Toy-MC #chi^{2};#chi^{2};toys",
                             nBinsChi2, 0.0, chi2_max);
    for (double x : v_chi2_bs)
        hChi2BS->Fill(x);
    for (double x : v_chi2_mc)
        hChi2MC->Fill(x);

    TCanvas *c5 = new TCanvas("c5", "Chi2 distributions", 1000, 450);
    c5->Divide(2, 1);

    // ---- left pad: chi2
    c5->cd(1);
    TH1F *hChi2BSn = (TH1F *)hChi2BS->Clone("hChi2BSn");
    TH1F *hChi2MCn = (TH1F *)hChi2MC->Clone("hChi2MCn");
    if (hChi2BSn->Integral() > 0)
        hChi2BSn->Scale(1.0 / hChi2BSn->Integral("width"));
    if (hChi2MCn->Integral() > 0)
        hChi2MCn->Scale(1.0 / hChi2MCn->Integral("width"));

    hChi2BSn->SetLineColor(kBlack);
    hChi2BSn->SetFillColorAlpha(kGray, 0.4);
    hChi2BSn->SetLineWidth(2);
    hChi2MCn->SetLineColor(kBlue + 1);
    hChi2MCn->SetFillColorAlpha(kBlue, 0.25);
    hChi2MCn->SetLineWidth(2);

    const double ymax5 = 1.3 * std::max(hChi2BSn->GetMaximum(), hChi2MCn->GetMaximum());
    hChi2BSn->SetMaximum(ymax5);
    hChi2BSn->SetTitle("#chi^{2} of toy fits");
    hChi2BSn->GetXaxis()->SetTitle("#chi^{2}");
    hChi2BSn->GetYaxis()->SetTitle("normalised toys / unit #chi^{2}");
    hChi2BSn->Draw("hist");
    hChi2MCn->Draw("hist same");

    // Reference chi2 pdf for the median ndf (use BS ndf)
    TF1 *fChi2 = new TF1("fChi2",
                         "ROOT::Math::chisquared_pdf(x,[0])",
                         0.0, chi2_max);
    fChi2->SetParameter(0, ndf_bs);
    fChi2->SetLineColor(kRed);
    fChi2->SetLineWidth(2);
    fChi2->Draw("same");

    auto *leg5a = new TLegend(0.55, 0.65, 0.89, 0.89);
    leg5a->AddEntry(hChi2BSn, "Bootstrap", "f");
    leg5a->AddEntry(hChi2MCn, "Toy-MC", "f");
    leg5a->AddEntry(fChi2, Form("#chi^{2} pdf, ndf=%d", ndf_bs), "l");
    leg5a->Draw();

    // ---- right pad: chi2 / ndf  (better for visual comparison)
    c5->cd(2);
    const double rmax = chi2_max / std::max(1, std::min(ndf_bs, ndf_mc));
    TH1F *hRedBS = new TH1F("hRedBS",
                            "Reduced #chi^{2};#chi^{2}/ndf;toys",
                            nBinsChi2, 0.0, rmax);
    TH1F *hRedMC = new TH1F("hRedMC",
                            "Reduced #chi^{2};#chi^{2}/ndf;toys",
                            nBinsChi2, 0.0, rmax);
    for (size_t i = 0; i < v_chi2_bs.size(); i++)
        if (v_ndf_bs[i] > 0)
            hRedBS->Fill(v_chi2_bs[i] / v_ndf_bs[i]);
    for (size_t i = 0; i < v_chi2_mc.size(); i++)
        if (v_ndf_mc[i] > 0)
            hRedMC->Fill(v_chi2_mc[i] / v_ndf_mc[i]);

    if (hRedBS->Integral() > 0)
        hRedBS->Scale(1.0 / hRedBS->Integral("width"));
    if (hRedMC->Integral() > 0)
        hRedMC->Scale(1.0 / hRedMC->Integral("width"));

    hRedBS->SetLineColor(kBlack);
    hRedBS->SetFillColorAlpha(kGray, 0.4);
    hRedBS->SetLineWidth(2);
    hRedMC->SetLineColor(kBlue + 1);
    hRedMC->SetFillColorAlpha(kBlue, 0.25);
    hRedMC->SetLineWidth(2);

    const double ymax5b = 1.3 * std::max(hRedBS->GetMaximum(), hRedMC->GetMaximum());
    hRedBS->SetMaximum(ymax5b);
    hRedBS->SetTitle("Reduced #chi^{2} of toy fits");
    hRedBS->GetXaxis()->SetTitle("#chi^{2} / ndf");
    hRedBS->GetYaxis()->SetTitle("normalised toys");
    hRedBS->Draw("hist");
    hRedMC->Draw("hist same");

    TLine *l1 = new TLine(1.0, 0.0, 1.0, ymax5b);
    l1->SetLineColor(kRed);
    l1->SetLineStyle(2);
    l1->SetLineWidth(2);
    l1->Draw();

    auto *leg5b = new TLegend(0.55, 0.70, 0.89, 0.89);
    leg5b->AddEntry(hRedBS, "Bootstrap", "f");
    leg5b->AddEntry(hRedMC, "Toy-MC", "f");
    leg5b->AddEntry(l1, "expected = 1", "l");
    leg5b->Draw();

    // Print mean reduced chi2 to stdout
    double meanRedBS = 0.0, meanRedMC = 0.0;
    int nRedBS = 0, nRedMC = 0;
    for (size_t i = 0; i < v_chi2_bs.size(); i++)
        if (v_ndf_bs[i] > 0)
        {
            meanRedBS += v_chi2_bs[i] / v_ndf_bs[i];
            nRedBS++;
        }
    for (size_t i = 0; i < v_chi2_mc.size(); i++)
        if (v_ndf_mc[i] > 0)
        {
            meanRedMC += v_chi2_mc[i] / v_ndf_mc[i];
            nRedMC++;
        }
    if (nRedBS)
        meanRedBS /= nRedBS;
    if (nRedMC)
        meanRedMC /= nRedMC;
    std::cout << "=== chi2 summary ===\n";
    std::cout << "Bootstrap: median ndf = " << ndf_bs
              << "   <chi2/ndf> = " << meanRedBS << "\n";
    std::cout << "Toy-MC   : median ndf = " << ndf_mc
              << "   <chi2/ndf> = " << meanRedMC << "\n";
    std::cout << "====================\n\n";

    // -------------------------------------------------------------------------
    // ---- Canvas 6 & 7: MINOS error distributions for each fraction (and N)
    //   Vertical reference lines:
    //     - red dashed   = sigma of the toy distribution (gold-standard)
    //     - dark green   = mean MINOS error
    // -------------------------------------------------------------------------
    auto buildErrorHists = [&](const std::vector<std::vector<double>> &ePlus,
                               const std::vector<std::vector<double>> &eMinus,
                               const std::vector<double> &eNplus,
                               const std::vector<double> &eNminus,
                               const std::vector<std::vector<double>> &v_phi,
                               const std::vector<double> &v_N,
                               const std::string &tag,
                               int colorTag) -> TCanvas *
    {
        TCanvas *cc = new TCanvas(Form("cMinos_%s", tag.c_str()),
                                  Form("MINOS errors %s", tag.c_str()),
                                  300 * (nC + 1), 350);
        cc->Divide(nC + 1, 1);

        // ---- N
        cc->cd(1);
        double xmaxN = 0.0;
        for (double x : eNplus)
            if (x > xmaxN)
                xmaxN = x;
        for (double x : eNminus)
            if (x > xmaxN)
                xmaxN = x;
        if (xmaxN <= 0)
            xmaxN = 1.0;
        xmaxN *= 1.2;

        TH1F *hEpN = new TH1F(Form("hEpN_%s", tag.c_str()),
                              Form("MINOS error N (%s);MINOS error;toys", tag.c_str()),
                              60, 0.0, xmaxN);
        TH1F *hEmN = new TH1F(Form("hEmN_%s", tag.c_str()),
                              "", 60, 0.0, xmaxN);
        for (double x : eNplus)
            hEpN->Fill(x);
        for (double x : eNminus)
            hEmN->Fill(x);

        hEpN->SetLineColor(colorTag);
        hEpN->SetFillColorAlpha(colorTag, 0.30);
        hEpN->SetLineWidth(2);
        hEmN->SetLineColor(colorTag + 2);
        hEmN->SetLineStyle(2);
        hEmN->SetLineWidth(2);

        const double ymN = 1.25 * std::max(hEpN->GetMaximum(), hEmN->GetMaximum());
        hEpN->SetMaximum(ymN);
        hEpN->Draw("hist");
        hEmN->Draw("hist same");

        // toy spread reference
        auto [mN, sN, p16N, p50N, p84N] = toyStats(v_N);
        TLine *lToyN = new TLine(sN, 0, sN, ymN);
        lToyN->SetLineColor(kRed);
        lToyN->SetLineStyle(2);
        lToyN->SetLineWidth(2);
        lToyN->Draw();

        const double meanEpN = meanOf(eNplus);
        TLine *lMeanN = new TLine(meanEpN, 0, meanEpN, ymN);
        lMeanN->SetLineColor(kGreen + 3);
        lMeanN->SetLineStyle(1);
        lMeanN->SetLineWidth(2);
        lMeanN->Draw();

        auto *legN = new TLegend(0.45, 0.65, 0.89, 0.89);
        legN->AddEntry(hEpN, "MINOS e+", "f");
        legN->AddEntry(hEmN, "MINOS e-", "l");
        legN->AddEntry(lToyN, Form("#sigma(toys) = %.2f", sN), "l");
        legN->AddEntry(lMeanN, Form("<e+> = %.2f", meanEpN), "l");
        legN->Draw();

        // ---- per fraction
        for (int k = 0; k < nC; k++)
        {
            cc->cd(k + 2);

            double xmax = 0.0;
            for (double x : ePlus[k])
                if (x > xmax)
                    xmax = x;
            for (double x : eMinus[k])
                if (x > xmax)
                    xmax = x;
            if (xmax <= 0)
                xmax = 0.05;
            xmax *= 1.2;

            TH1F *hEp = new TH1F(Form("hEp_%s_%d", tag.c_str(), k),
                                 Form("MINOS error #phi(%s) [%s];MINOS error;toys",
                                      labels[k].c_str(), tag.c_str()),
                                 60, 0.0, xmax);
            TH1F *hEm = new TH1F(Form("hEm_%s_%d", tag.c_str(), k),
                                 "", 60, 0.0, xmax);
            for (double x : ePlus[k])
                hEp->Fill(x);
            for (double x : eMinus[k])
                hEm->Fill(x);

            hEp->SetLineColor(colors[k]);
            hEp->SetFillColorAlpha(colors[k], 0.30);
            hEp->SetLineWidth(2);
            hEm->SetLineColor(colors[k] + 2);
            hEm->SetLineStyle(2);
            hEm->SetLineWidth(2);

            const double ym = 1.25 * std::max(hEp->GetMaximum(), hEm->GetMaximum());
            hEp->SetMaximum(ym);
            hEp->Draw("hist");
            hEm->Draw("hist same");

            auto [mP, sP, p16P, p50P, p84P] = toyStats(v_phi[k]);
            TLine *lToy = new TLine(sP, 0, sP, ym);
            lToy->SetLineColor(kRed);
            lToy->SetLineStyle(2);
            lToy->SetLineWidth(2);
            lToy->Draw();

            const double meanEp = meanOf(ePlus[k]);
            TLine *lMean = new TLine(meanEp, 0, meanEp, ym);
            lMean->SetLineColor(kGreen + 3);
            lMean->SetLineStyle(1);
            lMean->SetLineWidth(2);
            lMean->Draw();

            auto *lg = new TLegend(0.40, 0.65, 0.89, 0.89);
            lg->AddEntry(hEp, "MINOS e+", "f");
            lg->AddEntry(hEm, "MINOS e-", "l");
            lg->AddEntry(lToy, Form("#sigma(toys) = %.4f", sP), "l");
            lg->AddEntry(lMean, Form("<e+> = %.4f", meanEp), "l");
            lg->Draw();
        }
        return cc;
    };

    TCanvas *c6 = buildErrorHists(v_phi_eplus_bs, v_phi_eminus_bs,
                                  v_N_eplus_bs, v_N_eminus_bs,
                                  v_phi_bs, v_N_bs,
                                  "BS", kBlack);
    TCanvas *c7 = buildErrorHists(v_phi_eplus_mc, v_phi_eminus_mc,
                                  v_N_eplus_mc, v_N_eminus_mc,
                                  v_phi_mc, v_N_mc,
                                  "MC", kBlue);

    (void)c6;
    (void)c7; // silence unused warnings if compiled standalone
}