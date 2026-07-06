// ─────────────────────────────────────────────────────────────────────────────
//  subtract_background.C  —  22O population extraction from CALIFA gamma data
//
//  Multiplicity-resolved SIMULTANEOUS Poisson-likelihood fit of the
//  Doppler-corrected gamma spectra in the disjoint slices
//
//        M == 1 ,  M == 2 ,  M >= 3
//
//  with TWO signal components sharing ONE amplitude each across all slices:
//
//        N_direct   : decays of the 2+ (3199 keV) populated directly
//                     (template hResp_direct_3199_<slice>)
//        N_cascade  : decays of the 3+ (4582 keV), i.e. 1383+3199 cascade
//                     (template hResp_casc_1383_3199_<slice>)
//
//  Templates are normalised PER SIMULATED DECAY (histogram / nDecays_<state>),
//  so the fitted amplitudes are directly "number of decays of that state in
//  the data sample" — no efficiency correction is applied afterwards.
//  The relative slice populations of each template (how much of it sits at
//  M==1 vs M==2) is exactly the information that breaks the direct-vs-cascade
//  degeneracy of the old 1D fit; therefore the slices are NEVER individually
//  re-normalised.
//
//  BACKGROUND MODES (bgMode):
//    "23O"  (default) : per-slice background SHAPES taken from the 23O-gated
//                       gamma_spectra file (same reaction topology, ZERO bound
//                       excited states -> pure background), one free
//                       normalisation per slice.
//    "24O"            : same, using the 24O-gated file (cross-check /
//                       systematic).
//    "expo"           : free double exponential  A1*exp(k1*E) + A2*exp(k2*E),
//                       slopes k1,k2 shared across slices, amplitudes free
//                       per slice.  (The old FIXED-shape background is gone.)
//
//  TEMPLATE SHAPE FIT (doShapeFit, ON by default):
//    The simulated response peaks are usually slightly OFFSET from the data
//    (small beta/gain mismatch in the Doppler correction). By default the fit
//    applies ONLY a shared rigid energy SHIFT to both signal templates,
//         E_data = E_sim + shift ,
//    fitted and shared across all slices. This registers the 1383 and 3199
//    peaks onto the data without regenerating templates and without distorting
//    the amplitude split (a rigid shift preserves each template's integral, so
//    the fitted amplitudes remain "number of decays").
//
//    Two extra transforms are available but OFF by default, because each tends
//    to rail to its bound and drag the direct-vs-cascade split to a wrong
//    minimum:
//       fitStretch=kTRUE : also fit a multiplicative gain (E_data=shift+stretch*E_sim)
//       fitBroaden=kTRUE : also fit an extra Gaussian smearing of width sigExtra
//    Only enable these if a pure shift visibly fails to register BOTH peaks; if
//    it does, the right fix is usually the simulation beta, not more fit knobs.
//    Set doShapeFit=kFALSE for the legacy fixed-template fit.
//
//  GOODNESS OF FIT: the minimiser minimises the POISSON DEVIANCE (-2 lnL), the
//  correct statistic for counts. The legend/printout report it as an effective
//  chi2/ndf (global, and per slice) so a bad fit is obvious at a glance.
//
//  Populations: with nFragEvents (written by gammaSpectra.C = events after
//  PID+OPA(+NeuLAND) with NO gamma condition):
//        f(2+ direct) = N_direct  / nFragEvents
//        f(3+)        = N_cascade / nFragEvents
//        f(g.s.)      = 1 - f(2+ direct) - f(3+)
//
//  Higher-lying 22O states are intentionally NOT fitted: their feeding
//  thresholds in 24O* sit at ~11.8-13.9 MeV (at/near S3n), so they are
//  expected negligible; any unresolved feeder of the 2+ is absorbed in
//  N_direct (state this in the write-up).
//
//  USAGE (after gammaSpectra.C on 22O/23O/24O and buildAll22OResponseSlices):
//
//    root -l -b -q 'subtract_background.C'                       // all defaults
//    root -l -b -q 'subtract_background.C("gamma_spectra_22O.root",
//                     "gamma_response_functions.root","expo")'   // free 2-expo bg
//
//    // ALWAYS run the closure test once before trusting data numbers:
//    root -l -b
//      .L subtract_background.C+
//      closureTest();                       // matched-model toy validation
//
//  Requires gammaSpectra.h (GammaCfg) in the include path.
// ─────────────────────────────────────────────────────────────────────────────

#include "gammaSpectra.h"

#include <TFile.h>
#include <TH1F.h>
#include <TH1D.h>
#include <TString.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TPaveText.h>
#include <TLegend.h>
#include <TMath.h>
#include <TAxis.h>
#include <TROOT.h>
#include <TRandom3.h>
#include <TParameter.h>

#include <Math/Minimizer.h>
#include <Math/Factory.h>
#include <Math/Functor.h>

#include <iostream>
#include <iomanip>
#include <vector>
#include <functional>
#include <cmath>
#include <algorithm>

// ─────────────────────────────────────────────────────────────────────────────
//  Old fixed 24O double-exponential shape — kept ONLY as the default truth
//  shape of the closure test and as seed values for the "expo" mode.
// ─────────────────────────────────────────────────────────────────────────────
static const Double_t kBgP0 = 6.675;
static const Double_t kBgP1 = -2.45;
static const Double_t kBgP2 = 5.204;
static const Double_t kBgP3 = -0.4046;

static Double_t EvalOldFixedBackgroundShape(Double_t x)
{
    return TMath::Exp(kBgP0 + kBgP1 * x) + TMath::Exp(kBgP2 + kBgP3 * x);
}

// ─────────────────────────────────────────────────────────────────────────────
//  Small utilities
// ─────────────────────────────────────────────────────────────────────────────
namespace
{

// Fit-range bins (inclusive); identical binning enforced on all inputs.
int gBinLo = 1;
int gBinHi = 1;

// ─────────────────────────────────────────────────────────────────────────────
//  Template energy warp + Gaussian broadening.
//
//  The simulated response peaks are systematically offset from the data peaks
//  (mainly a small beta / gain mismatch between sim and experiment) and can be
//  a touch too narrow. Rather than regenerating templates by trial and error,
//  the fit applies a shared energy transform to EVERY template:
//
//      E_sim  ->  E_data = shift + stretch * E_sim              (registration)
//      then convolve with a Gaussian of width extraSigma        (broadening)
//
//  Equivalently, the model value at a DATA energy x is obtained by sampling
//  the template at the pre-image  (x - shift)/stretch , smeared by extraSigma.
//  shift, stretch, extraSigma are FITTED and SHARED across all slices and both
//  signal templates, so they cannot soak up the direct-vs-cascade information
//  (that lives in the multiplicity dimension).
//
//  gDoShapeFit master-enables the shape transform. The three components are
//  then individually switchable:
//    gFitShift   : rigid energy shift        (recommended; physical for a small
//                                              Doppler/beta offset)
//    gFitStretch : multiplicative gain       (OFF by default — a stretch scales
//                                              the two peaks by DIFFERENT amounts
//                                              and tends to fight the shift)
//    gFitBroaden : extra Gaussian smearing   (OFF by default — it easily rails
//                                              to its bound and props up a wrong
//                                              amplitude split)
//  Default = shift only.
// ─────────────────────────────────────────────────────────────────────────────
bool gDoShapeFit = true;
bool gFitShift   = true;
bool gFitStretch = false;
bool gFitBroaden = false;

// Linear-interpolated template lookup at an arbitrary energy (0 outside range).
inline double TemplateAt(const TH1D* h, double e)
{
    if (!h) return 0.0;
    const TAxis* ax = h->GetXaxis();
    if (e <= ax->GetXmin() || e >= ax->GetXmax()) return 0.0;
    const int    b  = ax->FindBin(e);
    const double xc = ax->GetBinCenter(b);
    const int    bN = (e >= xc) ? b + 1 : b - 1;
    const double y0 = h->GetBinContent(b);
    if (bN < 1 || bN > h->GetNbinsX()) return y0;
    const double xN = ax->GetBinCenter(bN);
    const double yN = h->GetBinContent(bN);
    const double t  = (e - xc) / (xN - xc);
    return y0 + t * (yN - y0);
}

// Warped + Gaussian-broadened template value at DATA energy x.
//   shift, stretch : E_data = shift + stretch * E_sim
//   sigmaExtra     : additional Gaussian smearing [MeV]; 0 => pure warp
// The Gaussian convolution is done by a short fixed-node quadrature in data
// energy; template bin width sets a sensible node spacing.
inline double WarpBroadenTemplate(const TH1D* h, double x,
                                  double shift, double stretch, double sigmaExtra)
{
    if (!h) return 0.0;
    if (stretch <= 1e-6) stretch = 1e-6;

    // Pure warp (no extra broadening): sample template at the pre-image of x.
    if (sigmaExtra < 1e-6) {
        const double ePre = (x - shift) / stretch;
        return TemplateAt(h, ePre) / stretch; // /stretch keeps the integral
    }

    // Convolve with Gaussian of width sigmaExtra in DATA energy.
    const double dx  = h->GetXaxis()->GetBinWidth(1);
    const double stepMin = 0.35 * dx;               // don't undersample the template
    const int    nHalf   = 6;                        // +-4 sigma-ish coverage
    double step = (2.0 * sigmaExtra) / nHalf;
    if (step < stepMin) step = stepMin;

    double acc = 0.0, wsum = 0.0;
    for (int k = -nHalf; k <= nHalf; ++k) {
        const double xk = x + k * step;
        const double w  = std::exp(-0.5 * (k * step / sigmaExtra) * (k * step / sigmaExtra));
        const double ePre = (xk - shift) / stretch;
        acc  += w * TemplateAt(h, ePre) / stretch;
        wsum += w;
    }
    return (wsum > 0.) ? acc / wsum : 0.0;
}

struct FitSlice
{
    TString tag;              // "mult1", "mult2", "mgeq3", or "all"
    TH1D*   data    = nullptr;
    TH1D*   rDir    = nullptr; // per-decay direct-3199 template
    TH1D*   rCas    = nullptr; // per-decay cascade template
    TH1D*   bgShape = nullptr; // unit-normalised (in fit range); nullptr in expo mode
};

struct FitOutcome
{
    bool   ok      = false;
    int    status  = -1;
    double Ndir = 0., eNdir = 0.;
    double Ncas = 0., eNcas = 0.;
    double corr = 0.;              // correlation(N_direct, N_cascade)
    double nllMin = 0.;            // -2 lnL at minimum (total, all slices)
    double shift = 0., stretch = 1., sigExtra = 0.; // fitted shape params
    double eShift = 0., eStretch = 0., eSigExtra = 0.;
    int    nFreePar = 0;          // number of free parameters
    std::vector<double> par, err;  // full parameter vector
    std::vector<double> dev;       // Poisson deviance per slice
    std::vector<double> nbins;     // bins per slice in range
    bool   hasMinos = false;
    double NdirLo = 0., NdirUp = 0., NcasLo = 0., NcasUp = 0.;

    double totalDev() const { double s=0; for(double d:dev) s+=d; return s; }
    double totalBins() const { double s=0; for(double n:nbins) s+=n; return s; }
    double ndf() const { return totalBins() - nFreePar; }
    double chi2ndf() const { const double d=ndf(); return d>0 ? totalDev()/d : -1.; }
};

// Parameter layout:
//   [0]=N_direct  [1]=N_cascade  [2]=shift  [3]=stretch  [4]=sigExtra
//   hist bg : [5+is]                         (one norm per slice)
//   expo bg : [5]=k1 [6]=k2  [7+2is]=A1_is [8+2is]=A2_is
constexpr int kNShape = 3;              // shift, stretch, sigExtra
constexpr int kBgStart = 2 + kNShape;   // = 5
int NPar(bool expo, int nS) { return expo ? (kBgStart + 2 + 2 * nS)
                                          : (kBgStart + nS); }

// Deep copy any TH1 into a detached TH1D (uniform arithmetic type everywhere).
TH1D* CloneAsTH1D(TH1* h, const TString& newName)
{
    if (!h) return nullptr;
    auto* out = new TH1D(newName, h->GetTitle(),
                         h->GetNbinsX(),
                         h->GetXaxis()->GetXmin(),
                         h->GetXaxis()->GetXmax());
    out->Sumw2();
    out->SetDirectory(nullptr);
    for (int b = 0; b <= h->GetNbinsX() + 1; ++b) {
        out->SetBinContent(b, h->GetBinContent(b));
        out->SetBinError(b, h->GetBinError(b));
    }
    return out;
}

TH1D* GetTH1(TFile* f, const TString& name, bool required = true)
{
    if (!f) return nullptr;
    TH1* h = dynamic_cast<TH1*>(f->Get(name));
    if (!h) {
        if (required)
            std::cerr << "[ERROR] histogram '" << name << "' not found in "
                      << f->GetName() << "\n"
                      << "        -> regenerate the file with the UPDATED "
                         "gammaSpectra.C / buildGammaResponseSlices\n";
        return nullptr;
    }
    return CloneAsTH1D(h, name + "_loc");
}

double GetParamD(TFile* f, const TString& name, double fallback = -1.)
{
    if (!f) return fallback;
    auto* p = dynamic_cast<TParameter<double>*>(f->Get(name));
    return p ? p->GetVal() : fallback;
}

bool SameBinning(const TH1* a, const TH1* b)
{
    if (!a || !b) return false;
    return a->GetNbinsX() == b->GetNbinsX()
        && std::abs(a->GetXaxis()->GetXmin() - b->GetXaxis()->GetXmin()) < 1e-9
        && std::abs(a->GetXaxis()->GetXmax() - b->GetXaxis()->GetXmax()) < 1e-9;
}

// Bare filenames -> $repopath/results/final/ (same convention as gammaSpectra.C)
TString ResolveDataPath(TString file)
{
    if (!file.Contains("/")) {
        const char* repopath = getenv("repopath");
        if (repopath) file = TString(repopath) + "/results/final/" + file;
    }
    return file;
}

// Model expectation in one bin of one slice.
//   p[2]=shift, p[3]=stretch, p[4]=sigExtra warp+broaden BOTH signal templates.
double ModelMu(const FitSlice& sl, int isl, const double* p, bool expo, int b)
{
    const double x = sl.data->GetXaxis()->GetBinCenter(b);

    double dir, cas;
    if (gDoShapeFit) {
        const double shift = p[2], stretch = p[3], sig = p[4];
        dir = WarpBroadenTemplate(sl.rDir, x, shift, stretch, sig);
        cas = WarpBroadenTemplate(sl.rCas, x, shift, stretch, sig);
    } else {
        dir = sl.rDir->GetBinContent(b);
        cas = sl.rCas->GetBinContent(b);
    }
    double mu = p[0] * dir + p[1] * cas;

    if (!expo) {
        mu += p[kBgStart + isl] * sl.bgShape->GetBinContent(b);
    } else {
        mu += p[kBgStart + 1 + 2 * isl] * std::exp(p[kBgStart]     * x)
            + p[kBgStart + 2 + 2 * isl] * std::exp(p[kBgStart + 1] * x);
    }
    return mu;
}

// ─────────────────────────────────────────────────────────────────────────────
//  Core: simultaneous Poisson-likelihood fit over an arbitrary set of slices.
//  Parameter layout (kBgStart = 5):
//    p[0]=N_direct  p[1]=N_cascade                      (shared amplitudes)
//    p[2]=shift  p[3]=stretch  p[4]=sigExtra            (shared template shape)
//    hist bg : p[5+is]                                  (one norm per slice)
//    expo bg : p[5]=k1  p[6]=k2  p[7+2is]=A1_is  p[8+2is]=A2_is
// ─────────────────────────────────────────────────────────────────────────────
FitOutcome RunCombinedFit(std::vector<FitSlice>& S,
                          bool expoMode,
                          bool doMinos     = false,
                          int  printLevel  = 0)
{
    FitOutcome R;
    const int nS   = static_cast<int>(S.size());
    const int npar = NPar(expoMode, nS);
    if (nS == 0) return R;

    // ---- -2 ln L (Poisson) ---------------------------------------------------
    std::function<double(const double*)> nll = [&S, expoMode, nS](const double* p) -> double
    {
        double sum = 0.0;
        for (int is = 0; is < nS; ++is) {
            const FitSlice& sl = S[is];
            for (int b = gBinLo; b <= gBinHi; ++b) {
                double mu = ModelMu(sl, is, p, expoMode, b);
                if (mu < 1e-9) mu = 1e-9;
                const double n = sl.data->GetBinContent(b);
                sum += 2.0 * (mu - n);
                if (n > 0.) sum += 2.0 * n * std::log(n / mu);
            }
        }
        return sum;
    };

    // ---- Seeds -----------------------------------------------------------------
    double T = 0.;
    std::vector<double> Ti(nS, 0.);
    for (int is = 0; is < nS; ++is) {
        Ti[is] = S[is].data->Integral(gBinLo, gBinHi);
        T += Ti[is];
    }
    double Idir = 0., Icas = 0.; // per-decay in-range yields, summed over slices
    for (int is = 0; is < nS; ++is) {
        Idir += S[is].rDir->Integral(gBinLo, gBinHi);
        Icas += S[is].rCas->Integral(gBinLo, gBinHi);
    }
    const double Ndir0 = (Idir > 0.) ? 0.35 * T / Idir : 100.;
    const double Ncas0 = (Icas > 0.) ? 0.15 * T / Icas : 50.;

    // ---- Minimizer ---------------------------------------------------------------
    ROOT::Math::Minimizer* min = ROOT::Math::Factory::CreateMinimizer("Minuit2", "Migrad");
    if (!min) min = ROOT::Math::Factory::CreateMinimizer("Minuit", "Migrad");
    if (!min) {
        std::cerr << "[ERROR] no minimizer available (Minuit2/Minuit)\n";
        return R;
    }
    min->SetMaxFunctionCalls(500000);
    min->SetTolerance(0.1);
    min->SetPrintLevel(printLevel);
    min->SetErrorDef(1.0);           // FCN = -2 lnL  ->  1-sigma at Delta = 1

    ROOT::Math::Functor fcn(nll, npar);
    min->SetFunction(fcn);

    auto stepOf = [](double v) { return 0.1 * std::abs(v) + 1.0; };

    min->SetLimitedVariable(0, "N_direct",  Ndir0, stepOf(Ndir0), 0., 1.e12);
    min->SetLimitedVariable(1, "N_cascade", Ncas0, stepOf(Ncas0), 0., 1.e12);

    // ---- Shared template shape parameters -----------------------------------
    //   shift   : additive energy offset [MeV]  (data = shift + stretch*sim)
    //   stretch : multiplicative gain  (~1)
    //   sigExtra: extra Gaussian broadening [MeV] (>=0)
    // By default ONLY the rigid shift is free; stretch and broadening are fixed
    // (they tend to rail to their bounds and drag the amplitude split off).
    if (gDoShapeFit && gFitShift)
        min->SetLimitedVariable(2, "shift", 0.0, 0.01, -0.80, 0.80);
    else
        min->SetFixedVariable(2, "shift", 0.0);

    if (gDoShapeFit && gFitStretch)
        min->SetLimitedVariable(3, "stretch", 1.0, 0.005, 0.90, 1.10);
    else
        min->SetFixedVariable(3, "stretch", 1.0);

    if (gDoShapeFit && gFitBroaden)
        min->SetLimitedVariable(4, "sigExtra", 0.05, 0.01, 0.0, 0.40);
    else
        min->SetFixedVariable(4, "sigExtra", 0.0);

    if (!expoMode) {
        for (int is = 0; is < nS; ++is) {
            const double b0 = std::max(0.5 * Ti[is], 1.0);
            min->SetLimitedVariable(kBgStart + is, Form("bg_%s", S[is].tag.Data()),
                                    b0, stepOf(b0), 0., 1.e12);
        }
    } else {
        const double k1seed = kBgP1;   // -2.45
        const double k2seed = kBgP3;   // -0.4046
        min->SetLimitedVariable(kBgStart,     "k1", k1seed, 0.1,  -15.0, -0.30);
        min->SetLimitedVariable(kBgStart + 1, "k2", k2seed, 0.05,  -3.0, -0.01);
        for (int is = 0; is < nS; ++is) {
            const double xLo  = S[is].data->GetXaxis()->GetBinCenter(gBinLo);
            const double cLo  = std::max(S[is].data->GetBinContent(gBinLo), 1.0);
            const int    bMid = S[is].data->GetXaxis()->FindBin(
                                    0.5 * (S[is].data->GetXaxis()->GetBinCenter(gBinLo)
                                         + S[is].data->GetXaxis()->GetBinCenter(gBinHi)));
            const double cMid = std::max(S[is].data->GetBinContent(bMid), 0.5);
            const double xMid = S[is].data->GetXaxis()->GetBinCenter(bMid);
            const double A1s  = 0.6 * cLo  * std::exp(-k1seed * xLo);
            const double A2s  = 0.8 * cMid * std::exp(-k2seed * xMid);
            min->SetLimitedVariable(kBgStart + 1 + 2 * is, Form("A1_%s", S[is].tag.Data()),
                                    A1s, stepOf(A1s), 0., 1.e15);
            min->SetLimitedVariable(kBgStart + 2 + 2 * is, Form("A2_%s", S[is].tag.Data()),
                                    A2s, stepOf(A2s), 0., 1.e12);
        }
    }

    const bool okMin = min->Minimize();
    min->Hesse();

    R.status = min->Status();
    R.ok     = okMin && (R.status == 0 || R.status == 1);
    R.nllMin = min->MinValue();

    const double* x = min->X();
    const double* e = min->Errors();
    R.par.assign(x, x + npar);
    R.err.assign(e, e + npar);
    R.Ndir  = x[0];  R.eNdir = e[0];
    R.Ncas  = x[1];  R.eNcas = e[1];
    R.corr  = min->Correlation(0, 1);
    R.shift   = x[2]; R.eShift   = e[2];
    R.stretch = x[3]; R.eStretch = e[3];
    R.sigExtra= x[4]; R.eSigExtra= e[4];
    // number of ACTUALLY free shape params (for ndf): shift/stretch/broaden
    const int nShapeFree = gDoShapeFit
        ? ((gFitShift ? 1 : 0) + (gFitStretch ? 1 : 0) + (gFitBroaden ? 1 : 0)) : 0;
    R.nFreePar = (npar - kNShape) + nShapeFree; // amplitudes+bg + free shape params

    if (doMinos) {
        double lo = 0., up = 0.;
        if (min->GetMinosError(0, lo, up)) { R.NdirLo = lo; R.NdirUp = up; R.hasMinos = true; }
        if (min->GetMinosError(1, lo, up)) { R.NcasLo = lo; R.NcasUp = up; }
    }

    // ---- Poisson deviance per slice ------------------------------------------------
    R.dev.assign(nS, 0.);
    R.nbins.assign(nS, 0.);
    for (int is = 0; is < nS; ++is) {
        double d = 0.;
        for (int b = gBinLo; b <= gBinHi; ++b) {
            double mu = ModelMu(S[is], is, x, expoMode, b);
            if (mu < 1e-9) mu = 1e-9;
            const double n = S[is].data->GetBinContent(b);
            d += 2.0 * (mu - n);
            if (n > 0.) d += 2.0 * n * std::log(n / mu);
        }
        R.dev[is]   = d;
        R.nbins[is] = gBinHi - gBinLo + 1;
    }

    delete min;
    return R;
}

// ─────────────────────────────────────────────────────────────────────────────
//  Drawing: one canvas per slice with data + components + total + numbers.
// ─────────────────────────────────────────────────────────────────────────────
void DrawSlice(TFile* fout,
               FitSlice& sl, int isl,
               const FitOutcome& R, bool expoMode,
               double xminFit, double xmaxFit,
               const TString& canvasName, const TString& label)
{
    if (!fout) return;
    fout->cd();

    const double* p = R.par.data();
    const double shift = R.shift, stretch = R.stretch, sig = R.sigExtra;

    // Build the drawn components from the SAME warped+broadened evaluation used
    // by the fit, so the red/blue curves sit exactly where the model put them.
    auto* hDir = static_cast<TH1D*>(sl.data->Clone(Form("hFit_dir_%s",  sl.tag.Data())));
    auto* hCas = static_cast<TH1D*>(sl.data->Clone(Form("hFit_cas_%s",  sl.tag.Data())));
    auto* hBg  = static_cast<TH1D*>(sl.data->Clone(Form("hFit_bg_%s",   sl.tag.Data())));
    hDir->SetDirectory(nullptr); hDir->Reset();
    hCas->SetDirectory(nullptr); hCas->Reset();
    hBg ->SetDirectory(nullptr); hBg ->Reset();

    for (int b = 1; b <= hDir->GetNbinsX(); ++b) {
        const double xc = hDir->GetXaxis()->GetBinCenter(b);
        double dir, cas;
        if (gDoShapeFit) {
            dir = R.Ndir * WarpBroadenTemplate(sl.rDir, xc, shift, stretch, sig);
            cas = R.Ncas * WarpBroadenTemplate(sl.rCas, xc, shift, stretch, sig);
        } else {
            dir = R.Ndir * sl.rDir->GetBinContent(b);
            cas = R.Ncas * sl.rCas->GetBinContent(b);
        }
        double bg = 0.;
        if (!expoMode) bg = p[kBgStart + isl] * sl.bgShape->GetBinContent(b);
        else bg = p[kBgStart + 1 + 2 * isl] * std::exp(p[kBgStart]     * xc)
                + p[kBgStart + 2 + 2 * isl] * std::exp(p[kBgStart + 1] * xc);
        hDir->SetBinContent(b, dir);
        hCas->SetBinContent(b, cas);
        hBg ->SetBinContent(b, bg);
    }

    auto* hTot = static_cast<TH1D*>(hDir->Clone(Form("hFit_total_%s", sl.tag.Data())));
    hTot->SetDirectory(nullptr);
    hTot->Add(hCas);
    hTot->Add(hBg);

    TCanvas* c = new TCanvas(canvasName, canvasName, 1100, 700);
    sl.data->SetMarkerStyle(24);
    sl.data->SetMarkerSize(0.65);
    sl.data->SetLineColor(kBlack);
    sl.data->GetXaxis()->SetRangeUser(xminFit, xmaxFit);
    sl.data->Draw("E");

    hBg->SetLineColor(kGray + 2);  hBg->SetLineStyle(2);  hBg->SetLineWidth(2);  hBg->Draw("hist same");
    hDir->SetLineColor(kRed);      hDir->SetLineWidth(2);                        hDir->Draw("hist same");
    hCas->SetLineColor(kBlue);     hCas->SetLineWidth(2);                        hCas->Draw("hist same");
    hTot->SetLineColor(kBlack);    hTot->SetLineWidth(2);                        hTot->Draw("hist same");

    const double yDir = hDir->Integral(gBinLo, gBinHi);
    const double yCas = hCas->Integral(gBinLo, gBinHi);
    const double yBg  = hBg->Integral(gBinLo, gBinHi);

    auto* pt = new TPaveText(0.52, 0.45, 0.99, 0.92, "NDC");
    pt->SetFillColor(0);
    pt->SetBorderSize(1);
    pt->SetTextSize(0.022);
    pt->SetTextAlign(12);
    pt->AddText(label);
    pt->AddText(Form("slice Poisson dev./bins = %.1f / %.0f = %.2f",
                     R.dev[isl], R.nbins[isl],
                     R.nbins[isl] > 0 ? R.dev[isl] / R.nbins[isl] : -1.));
    pt->AddText(Form("global #chi^{2}_{Poisson}/ndf = %.1f / %.0f = %.2f",
                     R.totalDev(), R.ndf(), R.chi2ndf()));
    pt->AddText(Form("N_{direct 2^{+}}  = %.0f #pm %.0f", R.Ndir, R.eNdir));
    pt->AddText(Form("N_{cascade 3^{+}} = %.0f #pm %.0f", R.Ncas, R.eNcas));
    if (gDoShapeFit && gFitShift)
        pt->AddText(Form("shift = %.3f #pm %.3f MeV", R.shift, R.eShift));
    if (gDoShapeFit && gFitStretch)
        pt->AddText(Form("stretch = %.4f #pm %.4f", R.stretch, R.eStretch));
    if (gDoShapeFit && gFitBroaden)
        pt->AddText(Form("extra #sigma = %.3f #pm %.3f MeV", R.sigExtra, R.eSigExtra));
    pt->AddText(Form("in-range counts: dir %.0f, cas %.0f, bg %.0f", yDir, yCas, yBg));
    pt->Draw("same");

    c->Write();

    hDir->Write(); hCas->Write(); hBg->Write(); hTot->Write();
}

void PrintOutcome(const char* head, const std::vector<FitSlice>& S,
                  const FitOutcome& R, bool expoMode)
{
    std::cout << "\n──── " << head << " ────────────────────────────────────────\n";
    std::cout << "  status " << R.status << (R.ok ? "  (OK)" : "  (CHECK!)") << "\n";
    std::cout << "  N_direct(2+)  = " << std::fixed << std::setprecision(1)
              << R.Ndir << " +- " << R.eNdir;
    if (R.hasMinos) std::cout << "   MINOS [" << R.NdirLo << ", +" << R.NdirUp << "]";
    std::cout << "\n  N_cascade(3+) = " << R.Ncas << " +- " << R.eNcas;
    if (R.hasMinos) std::cout << "   MINOS [" << R.NcasLo << ", +" << R.NcasUp << "]";
    std::cout << "\n  corr(N_dir, N_cas) = " << std::setprecision(3) << R.corr << "\n";
    if (gDoShapeFit && (gFitShift || gFitStretch || gFitBroaden)) {
        std::cout << std::setprecision(4) << "  template shape: ";
        if (gFitShift)   std::cout << "shift = "   << R.shift   << " +- " << R.eShift   << " MeV  ";
        if (gFitStretch) std::cout << "stretch = " << R.stretch << " +- " << R.eStretch << "  ";
        if (gFitBroaden) std::cout << "extra sigma = " << R.sigExtra << " +- " << R.eSigExtra << " MeV";
        std::cout << "\n";
    }
    std::cout << std::setprecision(1)
              << "  global Poisson chi2/ndf = " << R.totalDev() << " / " << R.ndf()
              << " = " << std::setprecision(2) << R.chi2ndf() << "\n";
    for (size_t is = 0; is < S.size(); ++is)
        std::cout << "  slice " << std::left << std::setw(6) << S[is].tag
                  << " deviance/bins = " << std::setprecision(1)
                  << R.dev[is] << " / " << R.nbins[is] << "\n";
    if (!expoMode) {
        for (size_t is = 0; is < S.size(); ++is)
            std::cout << "  bg counts (" << S[is].tag << ") = "
                      << R.par[kBgStart + is] << " +- " << R.err[kBgStart + is] << "\n";
    } else {
        std::cout << "  k1 = " << R.par[kBgStart] << " +- " << R.err[kBgStart]
                  << "   k2 = " << R.par[kBgStart + 1] << " +- " << R.err[kBgStart + 1] << "\n";
    }
    std::cout << std::setprecision(6);
}

} // anonymous namespace

// ═════════════════════════════════════════════════════════════════════════════
//  MAIN MACRO
// ═════════════════════════════════════════════════════════════════════════════
void subtract_background(TString dataFile22O = "gamma_spectra_22O.root",
                         TString respFile    = "gamma_response_functions.root",
                         TString bgMode      = "23O",   // "23O" | "24O" | "expo"
                         TString outFile     = "fit_populations_22O.root",
                         Double_t xminFit    = 0.5,
                         Double_t xmaxFit    = 10.0,
                         Bool_t doPerSliceCrossCheck = kTRUE,
                         Bool_t doAllMultFit = kTRUE,
                         Bool_t doMinos      = kFALSE,
                         Int_t  nSmoothBg    = 0,       // optional smoothing of bg SHAPES
                         TString bgFileOverride = "",   // explicit 23O/24O gamma file
                         Bool_t doShapeFit   = kTRUE,   // fit template shift (see below)
                         Bool_t fitStretch   = kFALSE,  // also fit a multiplicative gain
                         Bool_t fitBroaden   = kFALSE)  // also fit extra Gaussian sigma
{
    using namespace GammaCfg;
    gROOT->SetBatch(kTRUE);
    gDoShapeFit = doShapeFit;   // master enable
    gFitShift   = doShapeFit;   // rigid shift is the default transform
    gFitStretch = doShapeFit && fitStretch;
    gFitBroaden = doShapeFit && fitBroaden;

    const bool expoMode = (bgMode == "expo");
    if (!expoMode && bgMode != "23O" && bgMode != "24O") {
        std::cerr << "[ERROR] bgMode must be \"23O\", \"24O\" or \"expo\" (got "
                  << bgMode << ")\n";
        return;
    }

    dataFile22O = ResolveDataPath(dataFile22O);

    std::cout << "\n=== 22O population fit ===\n"
              << "  data (22O) : " << dataFile22O << "\n"
              << "  responses  : " << respFile    << "\n"
              << "  background : " << bgMode      << "\n"
              << "  shape fit  : "
              << (doShapeFit ? (TString("ON (shift") + (fitStretch ? "+stretch" : "")
                                + (fitBroaden ? "+broaden" : "") + ")")
                             : TString("OFF"))
              << "\n"
              << "  fit range  : [" << xminFit << ", " << xmaxFit << "] MeV\n\n";

    // ── Open inputs ───────────────────────────────────────────────────────────
    TFile* fData = TFile::Open(dataFile22O, "READ");
    if (!fData || fData->IsZombie()) {
        std::cerr << "[ERROR] cannot open " << dataFile22O << "\n";
        return;
    }
    TFile* fResp = TFile::Open(respFile, "READ");
    if (!fResp || fResp->IsZombie()) {
        std::cerr << "[ERROR] cannot open " << respFile << "\n";
        fData->Close();
        return;
    }

    TFile* fBg = nullptr;
    TString bgFile = bgFileOverride;
    if (!expoMode) {
        if (bgFile.IsNull()) {
            bgFile = dataFile22O;
            bgFile.ReplaceAll("22O", bgMode);   // gamma_spectra_23O.root / _24O.root
        } else {
            bgFile = ResolveDataPath(bgFile);
        }
        fBg = TFile::Open(bgFile, "READ");
        if (!fBg || fBg->IsZombie()) {
            std::cerr << "[ERROR] cannot open background file " << bgFile << "\n"
                      << "        run gammaSpectra.C on the " << bgMode
                      << "-gated data first, or use bgMode=\"expo\".\n";
            fData->Close(); fResp->Close();
            return;
        }
        std::cout << "  bg file    : " << bgFile << "\n";
    }

    // ── Book-keeping parameters ────────────────────────────────────────────────
    const double nFrag  = GetParamD(fData, NFragName());
    const double nCoinc = GetParamD(fData, NCoincName());
    const double nDecDir = GetParamD(fResp, NDecaysName(StateDirect()));
    const double nDecCas = GetParamD(fResp, NDecaysName(StateCascade()));
    if (nDecDir <= 0. || nDecCas <= 0.) {
        std::cerr << "[ERROR] nDecays parameters missing in " << respFile
                  << " — rebuild templates with buildGammaResponseSlices.\n";
        return;
    }

    // ── Load slices ────────────────────────────────────────────────────────────
    const std::vector<TString> tags = { SliceExact(1), SliceExact(2), SliceGeq(3) };
    std::vector<FitSlice> S;

    for (size_t is = 0; is < tags.size(); ++is) {
        FitSlice sl;
        sl.tag  = tags[is];
        sl.data = GetTH1(fData, DataHist(tags[is]));
        sl.rDir = GetTH1(fResp, RespHist(StateDirect(),  tags[is]));
        sl.rCas = GetTH1(fResp, RespHist(StateCascade(), tags[is]));
        if (!expoMode) sl.bgShape = GetTH1(fBg, DataHist(tags[is]));
        if (!sl.data || !sl.rDir || !sl.rCas || (!expoMode && !sl.bgShape)) return;
        S.push_back(sl);
    }

    // "all" slice for display + optional degeneracy demonstration
    FitSlice sAll;
    sAll.tag  = SliceAll();
    sAll.data = GetTH1(fData, DataHist(SliceAll()));
    sAll.rDir = GetTH1(fResp, RespHist(StateDirect(),  SliceAll()));
    sAll.rCas = GetTH1(fResp, RespHist(StateCascade(), SliceAll()));
    if (!expoMode) sAll.bgShape = GetTH1(fBg, DataHist(SliceAll()));
    const bool haveAll = sAll.data && sAll.rDir && sAll.rCas && (expoMode || sAll.bgShape);

    // ── Binning consistency ──────────────────────────────────────────────────
    for (auto& sl : S) {
        if (!SameBinning(sl.data, S[0].data) || !SameBinning(sl.rDir, S[0].data)
            || !SameBinning(sl.rCas, S[0].data)
            || (!expoMode && !SameBinning(sl.bgShape, S[0].data))) {
            std::cerr << "[ERROR] binning mismatch in slice " << sl.tag
                      << " — regenerate data/templates/background with the "
                         "shared GammaCfg binning.\n";
            return;
        }
    }
    if (S[0].data->GetNbinsX() != NBINS_E)
        std::cerr << "[WARN] binning differs from GammaCfg::NBINS_E — proceeding "
                     "since all inputs agree with each other.\n";

    gBinLo = S[0].data->GetXaxis()->FindBin(xminFit + 1e-6);
    gBinHi = S[0].data->GetXaxis()->FindBin(xmaxFit - 1e-6);

    // ── Normalise templates per decay; normalise bg shapes in fit range ───────
    for (auto& sl : S) {
        sl.rDir->Scale(1.0 / nDecDir);
        sl.rCas->Scale(1.0 / nDecCas);
        if (sl.bgShape) {
            if (nSmoothBg > 0) sl.bgShape->Smooth(nSmoothBg);
            const double I = sl.bgShape->Integral(gBinLo, gBinHi);
            if (I <= 0.) {
                std::cerr << "[ERROR] empty background shape for slice " << sl.tag
                          << " — use bgMode=\"expo\".\n";
                return;
            }
            sl.bgShape->Scale(1.0 / I);
        }
    }
    if (haveAll) {
        sAll.rDir->Scale(1.0 / nDecDir);
        sAll.rCas->Scale(1.0 / nDecCas);
        if (sAll.bgShape) {
            if (nSmoothBg > 0) sAll.bgShape->Smooth(nSmoothBg);
            const double I = sAll.bgShape->Integral(gBinLo, gBinHi);
            if (I > 0.) sAll.bgShape->Scale(1.0 / I);
        }
    }

    std::cout << "  <clusters/decay in range>  direct: "
              << S[0].rDir->Integral(gBinLo, gBinHi) + S[1].rDir->Integral(gBinLo, gBinHi)
                 + S[2].rDir->Integral(gBinLo, gBinHi)
              << "   cascade: "
              << S[0].rCas->Integral(gBinLo, gBinHi) + S[1].rCas->Integral(gBinLo, gBinHi)
                 + S[2].rCas->Integral(gBinLo, gBinHi) << "\n";

    // ── THE simultaneous fit ─────────────────────────────────────────────────
    FitOutcome R = RunCombinedFit(S, expoMode, doMinos, 1);
    PrintOutcome("SIMULTANEOUS FIT  (mult1 + mult2 + mgeq3, shared amplitudes)",
                 S, R, expoMode);

    // ── Populations ───────────────────────────────────────────────────────────
    double f2 = -1., f3 = -1., fgs = -1., ef2 = 0., ef3 = 0.;
    if (nFrag > 0.) {
        f2  = R.Ndir / nFrag;  ef2 = R.eNdir / nFrag;
        f3  = R.Ncas / nFrag;  ef3 = R.eNcas / nFrag;
        fgs = 1.0 - f2 - f3;
        std::cout << "\n──── POPULATIONS (nFragEvents = " << (long long)nFrag << ") ────\n"
                  << "  f(2+ direct feeding) = " << std::fixed << std::setprecision(4)
                  << f2 << " +- " << ef2 << "\n"
                  << "  f(3+  4582 keV)      = " << f3 << " +- " << ef3 << "\n"
                  << "  f(ground state)      = " << fgs
                  << "   (1 - f2 - f3; error ~ quad. sum)\n"
                  << "  [N.B. f(2+ direct) includes any unresolved feeders that "
                     "bypass the 3+; total 2+ decays = N_dir + N_cas]\n"
                  << std::setprecision(6);
    } else {
        std::cerr << "[WARN] nFragEvents not found in data file — populations "
                     "not computed (re-run updated gammaSpectra.C).\n";
    }
    if (nCoinc >= 0.)
        std::cout << "  1383(x)3199 coincidence pairs in data: " << (long long)nCoinc
                  << "   (quick check: ~ N_cas * eff1383 * eff3199)\n";

    // ── Per-slice cross-check fits (estimator consistency) ────────────────────
    std::vector<FitOutcome> Rslice;
    if (doPerSliceCrossCheck) {
        std::cout << "\n──── PER-SLICE CROSS-CHECK (each slice fitted alone) ────\n";
        std::cout << std::left << std::setw(10) << "slice"
                  << std::setw(24) << "N_direct" << std::setw(24) << "N_cascade"
                  << "status\n";
        for (size_t is = 0; is < S.size(); ++is) {
            std::vector<FitSlice> one = { S[is] };
            FitOutcome r1 = RunCombinedFit(one, expoMode, false, 0);
            Rslice.push_back(r1);
            const TString colD = Form("%.0f +- %.0f", r1.Ndir, r1.eNdir);
            const TString colC = Form("%.0f +- %.0f", r1.Ncas, r1.eNcas);
            std::cout << std::setw(10) << S[is].tag.Data()
                      << std::setw(24) << colD.Data()
                      << std::setw(24) << colC.Data()
                      << r1.status << "\n";
        }
        const TString colDc = Form("%.0f +- %.0f", R.Ndir, R.eNdir);
        const TString colCc = Form("%.0f +- %.0f", R.Ncas, R.eNcas);
        std::cout << std::setw(10) << "COMBINED"
                  << std::setw(24) << colDc.Data()
                  << std::setw(24) << colCc.Data()
                  << R.status << "\n"
                  << "  (mutual agreement validates the decay-scheme/response model;\n"
                  << "   the M==1-only fit is expected to be near-degenerate)\n";
    }

    // ── Optional: all-multiplicity single fit — the old, DEGENERATE config ────
    FitOutcome Rall;
    bool didAll = false;
    if (doAllMultFit && haveAll) {
        std::vector<FitSlice> vAll = { sAll };
        Rall = RunCombinedFit(vAll, expoMode, false, 0);
        didAll = true;
        PrintOutcome("ALL-MULT 1D FIT (demonstration — expect large corr ~ -1)",
                     vAll, Rall, expoMode);
    }

    // ── Output file: canvases, components, result parameters ──────────────────
    TFile* fOut = TFile::Open(outFile, "RECREATE");
    for (size_t is = 0; is < S.size(); ++is)
        DrawSlice(fOut, S[is], (int)is, R, expoMode, xminFit, xmaxFit,
                  Form("c_fit_22O_%s", S[is].tag.Data()),
                  Form("Simultaneous fit — slice %s", S[is].tag.Data()));
    if (didAll) {
        // draw all-mult with its own outcome (single-slice layout: isl = 0)
        DrawSlice(fOut, sAll, 0, Rall, expoMode, xminFit, xmaxFit,
                  "c_fit_22O_all", "All multiplicities (degenerate 1D demo)");
    }

    fOut->cd();
    TParameter<double>("N_direct",      R.Ndir ).Write();
    TParameter<double>("N_direct_err",  R.eNdir).Write();
    TParameter<double>("N_cascade",     R.Ncas ).Write();
    TParameter<double>("N_cascade_err", R.eNcas).Write();
    TParameter<double>("corr_dir_cas",  R.corr ).Write();
    TParameter<double>("nFragEvents",   nFrag  ).Write();
    if (f2 >= 0.) {
        TParameter<double>("f_2plus_direct",     f2 ).Write();
        TParameter<double>("f_2plus_direct_err", ef2).Write();
        TParameter<double>("f_3plus",            f3 ).Write();
        TParameter<double>("f_3plus_err",        ef3).Write();
        TParameter<double>("f_gs",               fgs).Write();
    }
    fOut->Close();

    fData->Close();
    fResp->Close();
    if (fBg) fBg->Close();

    std::cout << "\n[OK] results written to " << outFile << "\n";
}

// ═════════════════════════════════════════════════════════════════════════════
//  CLOSURE TEST — run BEFORE trusting data numbers.
//
//  Builds pseudo-data (Poisson toys) from the response templates plus a known
//  background, fits every toy with the SAME machinery, and reports the pull
//  distributions   (N_fit - N_true) / sigma_fit   for both amplitudes.
//  Unbiased means ~0 +- RMS ~1.
//
//  bgShapeFile: a gamma_spectra_23O.root to use as the truth background shape
//               (recommended); "" -> the old fixed double-exponential shape.
// ═════════════════════════════════════════════════════════════════════════════
void closureTest(TString respFile    = "gamma_response_functions.root",
                 TString bgShapeFile = "",
                 Double_t NdirTrue   = 3000.,
                 Double_t NcasTrue   = 1200.,
                 Double_t bgMult1    = 4000.,
                 Double_t bgMult2    = 1200.,
                 Double_t bgMgeq3    = 250.,
                 Int_t    nToys      = 200,
                 UInt_t   seed       = 12345,
                 Double_t xminFit    = 0.5,
                 Double_t xmaxFit    = 10.0,
                 Bool_t   shapeFitInToys = kFALSE) // toys are generated from raw
{                                                  // templates => shape is exactly
    using namespace GammaCfg;                      // (0,1,0); fix it by default so
    gROOT->SetBatch(kTRUE);                        // the test isolates the amplitudes
    gDoShapeFit = shapeFitInToys;
    gFitShift   = shapeFitInToys;
    gFitStretch = false;
    gFitBroaden = false;

    TFile* fResp = TFile::Open(respFile, "READ");
    if (!fResp || fResp->IsZombie()) {
        std::cerr << "[ERROR] cannot open " << respFile << "\n";
        return;
    }
    const double nDecDir = GetParamD(fResp, NDecaysName(StateDirect()));
    const double nDecCas = GetParamD(fResp, NDecaysName(StateCascade()));
    if (nDecDir <= 0. || nDecCas <= 0.) {
        std::cerr << "[ERROR] nDecays parameters missing — rebuild templates.\n";
        return;
    }

    TFile* fBg = nullptr;
    if (!bgShapeFile.IsNull()) {
        fBg = TFile::Open(ResolveDataPath(bgShapeFile), "READ");
        if (!fBg || fBg->IsZombie()) {
            std::cerr << "[WARN] cannot open " << bgShapeFile
                      << " — falling back to fixed double-expo truth shape.\n";
            fBg = nullptr;
        }
    }

    const std::vector<TString>  tags   = { SliceExact(1), SliceExact(2), SliceGeq(3) };
    const std::vector<Double_t> bgCnts = { bgMult1, bgMult2, bgMgeq3 };

    // Master slices (templates + truth bg shape); toy data hists swapped per toy.
    std::vector<FitSlice> S;
    for (size_t is = 0; is < tags.size(); ++is) {
        FitSlice sl;
        sl.tag  = tags[is];
        sl.rDir = GetTH1(fResp, RespHist(StateDirect(),  tags[is]));
        sl.rCas = GetTH1(fResp, RespHist(StateCascade(), tags[is]));
        if (!sl.rDir || !sl.rCas) return;
        sl.rDir->Scale(1.0 / nDecDir);
        sl.rCas->Scale(1.0 / nDecCas);

        if (fBg) sl.bgShape = GetTH1(fBg, DataHist(tags[is]));
        if (!sl.bgShape) {
            // fixed double-expo shape evaluated at bin centres
            sl.bgShape = static_cast<TH1D*>(sl.rDir->Clone(Form("bgTruth_%s", tags[is].Data())));
            sl.bgShape->SetDirectory(nullptr);
            sl.bgShape->Reset();
            for (int b = 1; b <= sl.bgShape->GetNbinsX(); ++b)
                sl.bgShape->SetBinContent(b,
                    EvalOldFixedBackgroundShape(sl.bgShape->GetXaxis()->GetBinCenter(b)));
        }
        // placeholder toy-data hist (same binning)
        sl.data = static_cast<TH1D*>(sl.rDir->Clone(Form("toy_%s", tags[is].Data())));
        sl.data->SetDirectory(nullptr);
        S.push_back(sl);
    }

    gBinLo = S[0].rDir->GetXaxis()->FindBin(xminFit + 1e-6);
    gBinHi = S[0].rDir->GetXaxis()->FindBin(xmaxFit - 1e-6);

    for (auto& sl : S) {
        const double I = sl.bgShape->Integral(gBinLo, gBinHi);
        if (I > 0.) sl.bgShape->Scale(1.0 / I);
    }

    // truth expectations per slice
    std::vector<std::vector<double>> mu(S.size(),
        std::vector<double>(S[0].rDir->GetNbinsX() + 2, 0.));
    for (size_t is = 0; is < S.size(); ++is)
        for (int b = gBinLo; b <= gBinHi; ++b)
            mu[is][b] = NdirTrue * S[is].rDir->GetBinContent(b)
                      + NcasTrue * S[is].rCas->GetBinContent(b)
                      + bgCnts[is] * S[is].bgShape->GetBinContent(b);

    TRandom3 rng(seed);
    int nFail = 0;
    double sPd = 0., sPd2 = 0., sPc = 0., sPc2 = 0.;
    double sNd = 0., sNc = 0., sEd = 0., sEc = 0.;
    int nOK = 0;

    std::cout << "\n=== CLOSURE TEST: " << nToys << " toys,  truth N_dir="
              << NdirTrue << "  N_cas=" << NcasTrue << " ===\n";

    for (int t = 0; t < nToys; ++t) {
        for (size_t is = 0; is < S.size(); ++is) {
            S[is].data->Reset();
            for (int b = gBinLo; b <= gBinHi; ++b)
                S[is].data->SetBinContent(b, rng.Poisson(std::max(mu[is][b], 0.)));
        }
        FitOutcome r = RunCombinedFit(S, /*expoMode=*/false, false, 0);
        if (!r.ok || r.eNdir <= 0. || r.eNcas <= 0.) { ++nFail; continue; }
        const double pd = (r.Ndir - NdirTrue) / r.eNdir;
        const double pc = (r.Ncas - NcasTrue) / r.eNcas;
        sPd += pd;  sPd2 += pd * pd;
        sPc += pc;  sPc2 += pc * pc;
        sNd += r.Ndir; sNc += r.Ncas; sEd += r.eNdir; sEc += r.eNcas;
        ++nOK;
    }

    if (nOK == 0) { std::cerr << "[ERROR] all toys failed.\n"; return; }
    const double mPd = sPd / nOK, rPd = std::sqrt(std::max(sPd2 / nOK - mPd * mPd, 0.));
    const double mPc = sPc / nOK, rPc = std::sqrt(std::max(sPc2 / nOK - mPc * mPc, 0.));

    std::cout << "  successful toys : " << nOK << "  (failed " << nFail << ")\n"
              << std::fixed << std::setprecision(3)
              << "  pull N_direct   : mean " << mPd << "   RMS " << rPd << "\n"
              << "  pull N_cascade  : mean " << mPc << "   RMS " << rPc << "\n"
              << std::setprecision(1)
              << "  <N_dir fit> = " << sNd / nOK << "  <sigma> = " << sEd / nOK << "\n"
              << "  <N_cas fit> = " << sNc / nOK << "  <sigma> = " << sEc / nOK << "\n"
              << "  PASS criterion: |mean| <~ 0.1 and RMS ~ 1.\n"
              << std::setprecision(6);

    fResp->Close();
    if (fBg) fBg->Close();
}

// ═════════════════════════════════════════════════════════════════════════════
//  Backwards-compatible wrapper (old entry-point name / signature style).
// ═════════════════════════════════════════════════════════════════════════════
void fit_gamma_with_response(TString inputFile = "gamma_spectra_22O.root",
                             TString respFile  = "gamma_response_functions.root",
                             Double_t xminFit  = 0.5,
                             Double_t xmaxFit  = 10.0,
                             Int_t nSmooth     = 0)
{
    subtract_background(inputFile, respFile, "23O", "fit_populations_22O.root",
                        xminFit, xmaxFit, kTRUE, kTRUE, kFALSE, nSmooth, "");
}
