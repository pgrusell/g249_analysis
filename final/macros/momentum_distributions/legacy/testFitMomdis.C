#include <RooRealVar.h>
#include <RooDataHist.h>
#include <RooHistPdf.h>
#include <RooAddPdf.h>
#include <RooFormulaVar.h>
#include <RooArgList.h>
#include <RooArgSet.h>
#include <RooFitResult.h>
#include <RooMsgService.h>
#include <RooChi2Var.h>
#include <TRandom3.h>
#include <TH1F.h>
#include <TCanvas.h>
#include <TLine.h>
#include <TLegend.h>
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
#include <memory>

// =============================================================================
// VERSION 1 (recursive, N fixed): RooAddPdf NON-extended, N_total = data integral.
// Fractions parametrised recursively a_k in [0,1]; physical phi_k as
// RooFormulaVar.  MINOS is run on a_k; asymmetric errors on phi_k come from
// a numerical Jacobian (a_k -> phi_k), N error is taken from the data integral
// (N is not a fit parameter here, so its statistical uncertainty is sqrt(N)).
//
// Bootstrap (Poisson-fluctuation of fake data) and toy-MC (event-by-event from
// truth) are both implemented. Output format matches the "extended" version
// so the two methods can be compared line by line.
//
// Components: configured below in main()
// =============================================================================

struct MomentaDist
{
    TH1F *Qz = nullptr;
    TH1F *Qt = nullptr;
    TH1F *Qy = nullptr;
    TH1F *Q = nullptr;
};

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
// Recursive <-> physical fraction conversions
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
        if (remain <= 1e-15)
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
// Generate fake data from binned templates (Poisson noise on bin contents)
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
// Event-by-event toy MC from the truth
// ---------------------------------------------------------------------------
TH1F *generateToyEventByEventN(const std::vector<TH1F *> &theos,
                               const std::vector<double> &phi_true,
                               double N,
                               int nBinsFit, double xmin, double xmax,
                               std::string name, TRandom3 *rng)
{
    TH1F *hToy = new TH1F(name.c_str(), name.c_str(), nBinsFit, xmin, xmax);
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
            if (u < cum[k])
            {
                kSel = k;
                break;
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
// FitResultN: SAME schema as the extended version, so prints/plots are common
// ---------------------------------------------------------------------------
struct FitResultN
{
    double N;
    std::vector<double> phi;
    double N_err;
    std::vector<double> phi_err;
    double N_eplus;
    double N_eminus;
    std::vector<double> phi_eplus;
    std::vector<double> phi_eminus;
    double chi2;
    int ndf;
    bool converged;
};

// ---------------------------------------------------------------------------
// RooFit model container: recursive non-extended
// ---------------------------------------------------------------------------
struct RooFitModelN
{
    RooRealVar *x = nullptr;
    std::vector<RooHistPdf *> pdf_tmpls;
    std::vector<RooRealVar *> a;      // (nC-1) recursive fractions
    std::vector<RooFormulaVar *> phi; // nC physical fractions
    RooAddPdf *model = nullptr;       // NON-extended
    int nC = 0;
};

// ---------------------------------------------------------------------------
// Numerical Jacobian d phi_k / d a_j (analytic version below would also work,
// but we keep it numerical to mirror the extended-version code)
// ---------------------------------------------------------------------------
static std::vector<double> aToPhi(const std::vector<double> &a)
{
    return recursiveToPhysical(a);
}

static void jacobianAtoPhi(const std::vector<double> &a_central,
                           std::vector<std::vector<double>> &J)
{
    const int Na = (int)a_central.size();
    const int Np = Na + 1;
    J.assign(Np, std::vector<double>(Na, 0.0));
    const std::vector<double> phi0 = aToPhi(a_central);
    for (int j = 0; j < Na; j++)
    {
        const double h = std::max(1e-6,
                                  1e-4 * std::max(1.0, std::fabs(a_central[j])));
        std::vector<double> a_pert = a_central;
        a_pert[j] += h;
        if (a_pert[j] > 1.0)
            a_pert[j] = 1.0;
        const double dh = a_pert[j] - a_central[j];
        if (std::fabs(dh) < 1e-15)
            continue;
        const std::vector<double> phi_pert = aToPhi(a_pert);
        for (int k = 0; k < Np; k++)
            J[k][j] = (phi_pert[k] - phi0[k]) / dh;
    }
}

// ---------------------------------------------------------------------------
// Sign-aware propagation of MINOS asymmetric errors a -> phi
// ---------------------------------------------------------------------------
static void propagateMinosAtoPhi(const std::vector<double> &a_central,
                                 const std::vector<double> &a_eplus,
                                 const std::vector<double> &a_eminus,
                                 std::vector<double> &phi_eplus,
                                 std::vector<double> &phi_eminus)
{
    const int Na = (int)a_central.size();
    const int Np = Na + 1;
    std::vector<std::vector<double>> J;
    jacobianAtoPhi(a_central, J);
    phi_eplus.assign(Np, 0.0);
    phi_eminus.assign(Np, 0.0);
    for (int k = 0; k < Np; k++)
    {
        double sp = 0.0, sm = 0.0;
        for (int j = 0; j < Na; j++)
        {
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

// ---------------------------------------------------------------------------
// doFitN: RooFit fit + chi2 + MINOS, with the same FitResultN interface.
// Here N is FIXED to hExp->Integral() (model is non-extended).
// ---------------------------------------------------------------------------
FitResultN doFitN(TH1F *hExp, RooFitModelN &M,
                  double /*N_init*/, const std::vector<double> &phi_seed,
                  bool verbose = false, bool runMinos = false)
{
    using namespace RooFit;
    const int nC = M.nC;
    const std::vector<double> a_seed = physicalToRecursive(phi_seed);
    for (int k = 0; k < nC - 1; k++)
    {
        const double v = std::min(1.0, std::max(0.0, a_seed[k]));
        M.a[k]->setVal(v);
    }
    RooDataHist data("data_fit", "data", RooArgList(*M.x), hExp);
    std::unique_ptr<RooFitResult> fitRes(
        M.model->fitTo(data,
                       Save(true),
                       PrintLevel(verbose ? 1 : -1),
                       Warnings(false),
                       Minos(runMinos)));

    FitResultN res;
    res.phi.assign(nC, 0.0);
    res.phi_err.assign(nC, 0.0);
    res.phi_eplus.assign(nC, 0.0);
    res.phi_eminus.assign(nC, 0.0);
    res.N_eplus = 0.0;
    res.N_eminus = 0.0;

    // N is FIXED (non-extended): take it from the data integral
    res.N = hExp->Integral();
    // The "statistical" uncertainty associated with N here is sqrt(N) (Poisson)
    res.N_err = (res.N > 0) ? std::sqrt(res.N) : 0.0;
    res.N_eplus = res.N_err;
    res.N_eminus = res.N_err;

    for (int k = 0; k < nC; k++)
    {
        res.phi[k] = M.phi[k]->getVal();
        if (fitRes)
            res.phi_err[k] = M.phi[k]->getPropagatedError(*fitRes);
    }

    if (runMinos && fitRes)
    {
        std::vector<double> a_c(nC - 1), a_ep(nC - 1), a_em(nC - 1);
        for (int k = 0; k < nC - 1; k++)
        {
            a_c[k] = M.a[k]->getVal();
            a_ep[k] = M.a[k]->getAsymErrorHi();
            a_em[k] = std::fabs(M.a[k]->getAsymErrorLo());
            if (a_ep[k] == 0.0)
                a_ep[k] = M.a[k]->getError();
            if (a_em[k] == 0.0)
                a_em[k] = M.a[k]->getError();
        }
        propagateMinosAtoPhi(a_c, a_ep, a_em, res.phi_eplus, res.phi_eminus);
    }

    // chi2 (NON-extended): compares N_data * P(bin) to data, with Poisson errors
    RooChi2Var chi2Var("chi2Var_fit", "chi2", *M.model, data,
                       DataError(RooAbsData::Poisson));
    res.chi2 = chi2Var.getVal();

    int nBinsUsed = 0;
    for (int i = 1; i <= hExp->GetNbinsX(); i++)
        if (hExp->GetBinContent(i) > 0)
            nBinsUsed++;
    res.ndf = nBinsUsed - (nC - 1); // (nC-1) free params: a_0..a_{nC-2}

    res.converged = (fitRes && fitRes->status() == 0);
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
    using namespace RooFit;
    RooMsgService::instance().setGlobalKillBelow(RooFit::WARNING);

    // -------------------------------------------------------------------------
    // CONFIGURATION (must match the extended version for a fair comparison)
    // -------------------------------------------------------------------------
    std::vector<std::string> labels = {"1p_{3/2}", "1s_{1/2}", "1d_{5/2}"};
    std::vector<double> phi_true = {0.6, 0.1, 0.3};

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
    const int nToys = 5000;
    const int fitNBins = 50;
    const double fitMax = 300.0;

    const int nC = (int)labels.size();

    std::cout << "\n========== TEST CONFIG (recursive, N fixed, " << nC << " components) ==========\n";
    for (int k = 0; k < nC; k++)
        std::cout << "phi_true(" << labels[k] << ") = " << phi_true[k] << "\n";
    std::cout << "N_true = " << N_true << "\n";
    std::cout << "nToys  = " << nToys << "\n";
    std::cout << "================================================\n\n";

    // -------------------------------------------------------------------------
    // Load templates (paths must match the extended version)
    // -------------------------------------------------------------------------
    std::vector<std::string> inFilesTheo = {
        "/nucl_lustre/pablogrusell/g249/g249_analysis/theory/JT/25F/sigt_p32-gs.txt",
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

    // -------------------------------------------------------------------------
    // Build the recursive non-extended model
    // -------------------------------------------------------------------------
    RooFitModelN M;
    M.nC = nC;
    M.x = new RooRealVar("x", "p_{x} [MeV/c]", -fitMax, fitMax);

    M.pdf_tmpls.resize(nC);
    RooArgList pdfList;
    for (int k = 0; k < nC; k++)
    {
        auto *dh = new RooDataHist(Form("dh_tmpl_%d", k), Form("dh_tmpl_%d", k),
                                   RooArgList(*M.x), tmpls[k]);
        M.pdf_tmpls[k] = new RooHistPdf(Form("pdf_tmpl_%d", k),
                                        Form("pdf %s", labels[k].c_str()),
                                        RooArgSet(*M.x), *dh, 0);
        pdfList.add(*M.pdf_tmpls[k]);
    }

    M.a.resize(nC - 1);
    RooArgList fracList;
    {
        std::vector<double> a_seed = physicalToRecursive(std::vector<double>(nC, 1.0 / nC));
        for (int k = 0; k < nC - 1; k++)
        {
            M.a[k] = new RooRealVar(Form("a%d", k),
                                    Form("recursive frac a_%d", k),
                                    a_seed[k], 0.0, 1.0);
            fracList.add(*M.a[k]);
        }
    }
    M.model = new RooAddPdf("model", "recursive sum of templates",
                            pdfList, fracList, /*recursiveFractions=*/true);

    // Physical fractions phi_k as RooFormulaVar
    M.phi.resize(nC);
    for (int k = 0; k < nC; k++)
    {
        std::string expr;
        RooArgList deps;
        if (k < nC - 1)
        {
            for (int j = 0; j < k; j++)
            {
                expr += Form("(1-@%d)*", j);
                deps.add(*M.a[j]);
            }
            expr += Form("@%d", k);
            deps.add(*M.a[k]);
        }
        else
        {
            for (int j = 0; j < nC - 1; j++)
            {
                if (j > 0)
                    expr += "*";
                expr += Form("(1-@%d)", j);
                deps.add(*M.a[j]);
            }
        }
        M.phi[k] = new RooFormulaVar(Form("phi_%d", k),
                                     Form("physical fraction %s", labels[k].c_str()),
                                     expr.c_str(), deps);
    }

    // -------------------------------------------------------------------------
    // Generate fake data + nominal fit
    // -------------------------------------------------------------------------
    TRandom3 rng(0);
    TH1F *hFake = generateFakeDataN(tmpls, phi_true, N_true, "hFake", &rng);
    std::cout << "Generated fake data: entries = " << hFake->GetEntries()
              << "  integral = " << hFake->Integral() << "\n\n";

    std::vector<double> phi_seed(nC, 1.0 / nC);
    FitResultN nominal = doFitN(hFake, M, hFake->Integral(), phi_seed,
                                /*verbose*/ true, /*runMinos*/ true);

    std::cout << "\n=== Nominal fit on fake data (recursive, N fixed) ===\n";
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
    // BOOTSTRAP toys
    // -------------------------------------------------------------------------
    TH1F *hFakeToy = (TH1F *)hFake->Clone("hFakeToy");
    std::vector<std::vector<double>> v_phi_bs(nC);
    std::vector<double> v_N_bs;
    std::vector<double> v_chi2_bs;
    std::vector<int> v_ndf_bs;
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
    std::cout << "Running " << nToys << " BOOTSTRAP toys (recursive) with MINOS...\n";

    for (int t = 0; t < nToys; t++)
    {
        for (int i = 1; i <= nBins; i++)
        {
            const double mu = hFake->GetBinContent(i);
            const double k = rng.Poisson(mu);
            hFakeToy->SetBinContent(i, k);
            hFakeToy->SetBinError(i, (k > 0) ? std::sqrt(k) : 1.0);
        }
        FitResultN r = doFitN(hFakeToy, M, nominal.N, nominal.phi,
                              false, /*runMinos*/ true);
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
    // TOY-MC: regenerate from truth
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
    std::cout << "Running " << nToys << " TOY-MC toys (recursive) with MINOS...\n";

    for (int t = 0; t < nToys; t++)
    {
        TH1F *hToy = generateToyEventByEventN(theosFine, phi_true, N_true,
                                              fitNBins, -fitMax, fitMax,
                                              Form("hToyMC_%d", t), &rng);
        const double Ninit = (hToy->Integral() > 0) ? hToy->Integral() : N_true;
        FitResultN r = doFitN(hToy, M, Ninit, nominal.phi,
                              false, /*runMinos*/ true);
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

    // -------------------------------------------------------------------------
    // SUMMARY (identical layout to the extended version)
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

    auto fmt = [](double x)
    { std::ostringstream o;
      o << std::fixed << std::setprecision(4) << x;
      return o.str(); };

    std::cout << "=== SIGMA COMPARISON (faithfulness test) ===\n";
    std::cout << "                  bootstrap     toy-MC      ratio MC/BS\n";
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
    // PLOTS (identical to the extended-version plots)
    // -------------------------------------------------------------------------
    std::vector<int> colors = {kBlue, kGreen + 2, kMagenta, kCyan + 1};

    // Canvas 1: nominal fit
    TCanvas *c1 = new TCanvas("c1_rec", "Nominal fit (recursive)", 900, 650);
    TH1F *total = (TH1F *)tmpls[0]->Clone("fit_total_rec");
    for (int i = 1; i <= total->GetNbinsX(); i++)
        total->SetBinContent(i, 0);
    std::vector<TH1F *> comps(nC);
    for (int k = 0; k < nC; k++)
    {
        comps[k] = (TH1F *)tmpls[k]->Clone(Form("comp_rec_%d", k));
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
    hFake->SetTitle("Fake data vs fit (recursive, N fixed)");
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
    leg->AddEntry((TObject *)nullptr,
                  Form("#chi^{2}/ndf = %.2f/%d = %.2f", nominal.chi2, nominal.ndf,
                       (nominal.ndf > 0 ? nominal.chi2 / nominal.ndf : -1.0)),
                  "");
    leg->Draw();

    // Canvas 2: bootstrap distributions
    TCanvas *c2 = new TCanvas("c2_rec", "Bootstrap distributions (recursive)",
                              300 * nC, 350);
    c2->Divide(nC, 1);
    std::vector<TH1F *> hPhiBS(nC);
    for (int k = 0; k < nC; k++)
    {
        hPhiBS[k] = new TH1F(Form("hPhiBS_rec_%d", k),
                             Form("Bootstrap a(%s);a;toys", labels[k].c_str()),
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

    // ---- Canvas 3: toy-MC distributions + nominal MINOS interval ----
    TCanvas *c3 = new TCanvas("c3_rec", "Toy-MC distributions (recursive)",
                              300 * nC, 350);
    c3->Divide(nC, 1);
    std::vector<TH1F *> hPhiMC(nC);
    for (int k = 0; k < nC; k++)
    {
        hPhiMC[k] = new TH1F(Form("hPhiMC_rec_%d", k),
                             Form("Toy-MC a(%s);a;toys", labels[k].c_str()),
                             80, 0.0, 1.0);
        for (double v : v_phi_mc[k])
            hPhiMC[k]->Fill(v);

        c3->cd(k + 1);
        hPhiMC[k]->SetStats(0);
        hPhiMC[k]->SetLineColor(colors[k]);
        hPhiMC[k]->SetFillColorAlpha(colors[k], 0.3);
        hPhiMC[k]->Draw();

        const double ymax = hPhiMC[k]->GetMaximum();

        // Compute the two dispersions to compare
        auto [m_toy, s_toy, p16, p50, p84] = toyStats(v_phi_mc[k]);
        const double sMinos = 0.5 * (nominal.phi_eplus[k] + nominal.phi_eminus[k]);
        const double ratio = (s_toy > 0) ? sMinos / s_toy : 0.0;

        // truth (red dashed)
        TLine *lt = new TLine(phi_true[k], 0, phi_true[k], ymax);
        lt->SetLineColor(kRed);
        lt->SetLineStyle(2);
        lt->SetLineWidth(2);
        lt->Draw();

        // // MINOS interval from the nominal fit (blue dashed)
        // const double xLo = nominal.phi[k] - nominal.phi_eminus[k];
        // const double xHi = nominal.phi[k] + nominal.phi_eplus[k];
        // TLine *lLo = new TLine(xLo, 0, xLo, ymax);
        // TLine *lHi = new TLine(xHi, 0, xHi, ymax);
        // for (TLine *l : {lLo, lHi})
        // {
        //     l->SetLineColor(kBlue + 1);
        //     l->SetLineStyle(2);
        //     l->SetLineWidth(2);
        //     l->Draw();
        // }

        // // nominal central value (solid blue)
        // TLine *lNom = new TLine(nominal.phi[k], 0, nominal.phi[k], ymax);
        // lNom->SetLineColor(kBlue + 1);
        // lNom->SetLineStyle(1);
        // lNom->SetLineWidth(2);
        // lNom->Draw();

        auto *lg = new TLegend(0.55, 0.62, 0.92, 0.92);
        lg->SetTextSize(0.030);
        lg->AddEntry(hPhiMC[k], "Toy-MC", "f");
        lg->AddEntry(lt, Form("truth = %.3f", phi_true[k]), "l");
        // lg->AddEntry(lNom, Form("nominal #hat{a} = %.3f", nominal.phi[k]), "l");
        // lg->AddEntry(lLo,
        //              Form("MINOS [-%.3f, +%.3f]",
        //                   nominal.phi_eminus[k], nominal.phi_eplus[k]),
        //              "l");
        lg->AddEntry((TObject *)nullptr,
                     Form("#sigma_{MINOS} = %.4f", sMinos), "");
        lg->AddEntry((TObject *)nullptr,
                     Form("#sigma_{toys}  = %.4f", s_toy), "");
        lg->AddEntry((TObject *)nullptr,
                     Form("ratio = %.2f", ratio), "");
        lg->Draw();
    }

    // Canvas 4: overlay BS vs MC
    TCanvas *c4 = new TCanvas("c4_rec", "Overlay BS vs MC (recursive)",
                              300 * nC, 350);
    c4->Divide(nC, 1);
    for (int k = 0; k < nC; k++)
    {
        c4->cd(k + 1);
        TH1F *hBS = (TH1F *)hPhiBS[k]->Clone(Form("BSn_rec_%d", k));
        TH1F *hMC = (TH1F *)hPhiMC[k]->Clone(Form("MCn_rec_%d", k));
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
        hBS->SetTitle(Form("%s : BS vs MC (recursive)", labels[k].c_str()));
        hBS->GetXaxis()->SetTitle(Form("a(%s)", labels[k].c_str()));
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

    // Canvas 5: chi2 distributions
    auto medianInt = [](std::vector<int> v) -> int
    {
        if (v.empty())
            return 0;
        std::sort(v.begin(), v.end());
        return v[v.size() / 2];
    };
    const int ndf_bs = medianInt(v_ndf_bs);
    const int ndf_mc = medianInt(v_ndf_mc);

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

    TH1F *hChi2BS = new TH1F("hChi2BS_rec",
                             "Bootstrap #chi^{2};#chi^{2};toys",
                             nBinsChi2, 0.0, chi2_max);
    TH1F *hChi2MC = new TH1F("hChi2MC_rec",
                             "Toy-MC #chi^{2};#chi^{2};toys",
                             nBinsChi2, 0.0, chi2_max);
    for (double x : v_chi2_bs)
        hChi2BS->Fill(x);
    for (double x : v_chi2_mc)
        hChi2MC->Fill(x);

    TCanvas *c5 = new TCanvas("c5_rec", "Chi2 distributions (recursive)", 1000, 450);
    c5->Divide(2, 1);
    c5->cd(1);
    TH1F *hChi2BSn = (TH1F *)hChi2BS->Clone("hChi2BSn_rec");
    TH1F *hChi2MCn = (TH1F *)hChi2MC->Clone("hChi2MCn_rec");
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
    hChi2BSn->SetTitle("#chi^{2} of toy fits (recursive)");
    hChi2BSn->GetXaxis()->SetTitle("#chi^{2}");
    hChi2BSn->GetYaxis()->SetTitle("normalised toys / unit #chi^{2}");
    hChi2BSn->Draw("hist");
    hChi2MCn->Draw("hist same");
    TF1 *fChi2 = new TF1("fChi2_rec",
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

    c5->cd(2);
    const double rmax = chi2_max / std::max(1, std::min(ndf_bs, ndf_mc));
    TH1F *hRedBS = new TH1F("hRedBS_rec",
                            "Reduced #chi^{2};#chi^{2}/ndf;toys",
                            nBinsChi2, 0.0, rmax);
    TH1F *hRedMC = new TH1F("hRedMC_rec",
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
    hRedBS->SetTitle("Reduced #chi^{2} of toy fits (recursive)");
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

    // Canvas 6 & 7: MINOS error distributions
    auto buildErrorHists = [&](const std::vector<std::vector<double>> &ePlus,
                               const std::vector<std::vector<double>> &eMinus,
                               const std::vector<double> &eNplus,
                               const std::vector<double> &eNminus,
                               const std::vector<std::vector<double>> &v_phi,
                               const std::vector<double> &v_N,
                               const std::string &tag,
                               int colorTag) -> TCanvas *
    {
        TCanvas *cc = new TCanvas(Form("cMinos_rec_%s", tag.c_str()),
                                  Form("MINOS errors %s (recursive)", tag.c_str()),
                                  300 * (nC + 1), 350);
        cc->Divide(nC + 1, 1);
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
        TH1F *hEpN = new TH1F(Form("hEpN_rec_%s", tag.c_str()),
                              Form("MINOS error N (%s);MINOS error;toys", tag.c_str()),
                              60, 0.0, xmaxN);
        TH1F *hEmN = new TH1F(Form("hEmN_rec_%s", tag.c_str()),
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
            TH1F *hEp = new TH1F(Form("hEp_rec_%s_%d", tag.c_str(), k),
                                 Form("MINOS error a(%s) [%s];MINOS error;toys",
                                      labels[k].c_str(), tag.c_str()),
                                 60, 0.0, xmax);
            TH1F *hEm = new TH1F(Form("hEm_rec_%s_%d", tag.c_str(), k),
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
                                  v_phi_bs, v_N_bs, "BS", kBlack);
    TCanvas *c7 = buildErrorHists(v_phi_eplus_mc, v_phi_eminus_mc,
                                  v_N_eplus_mc, v_N_eminus_mc,
                                  v_phi_mc, v_N_mc, "MC", kBlue);
    (void)c6;
    (void)c7;
}