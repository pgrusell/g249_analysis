/** --------------------------------------------------------------------
 **
 **  fitExc.C
 **
 **  Fit of the 24O* excitation-energy spectrum from the 1n-breakup
 **  channel  25F(p,2p)24O* -> 23O + n  (G249).
 **
 **  Model:
 **     N(Eexc) = Sum_k A_k * [ BW_k(Erel) (x) R(Erel_tru -> Erel_rec) ]
 **               + A_bg * BG(Eexc)
 **
 **   * BW_k : Breit-Wigner with energy-dependent width
 **                Gamma(E) = Gamma0 * (E/E0)^(l+1/2)
 **            (low-energy neutron penetrability limit; l per resonance,
 **            set in gLorb[] below -- default l=2, d3/2 neutron decay).
 **            Set gEnergyDepWidth = kFALSE for a plain fixed-width BW.
 **
 **   * R    : response matrix built directly from the "ana" TTree of
 **            simFile (branches Erel, ErelTrue, GeV; -999 = no neutron).
 **            X = Erel_true (MeV, gNTrueBins bins over the true range,
 **            auto-detected from the tree and rounded outward to 0.1 MeV),
 **            Y = Eexc = Erel*1000+Sn with EXACTLY the data binning, so
 **            no rebinning is needed when folding. Each truth column is
 **            normalized by hTrueCount (generated ErelTrue histogram
 **            from the same tree), so the Erel-dependent (neutron)
 **            efficiency is folded in and the fitted A_k are
 **            efficiency-corrected yields. If useEfficiency is kFALSE
 **            each column is normalized to unity instead (pure line
 **            shape, flat efficiency assumed).
 **
 **   * BG   : non-resonant background, one of three modes (gBgMode):
 **              gBgMode == 0 : template built directly from the "ana"
 **                             TTree of bgFile as Erel*1000+Sn, with
 **                             the data binning; one free scale A_bg.
 **              gBgMode == 1 : parametric truth-level continuum
 **                             f_bg(Erel) = Erel^a * exp(-Erel/T),
 **                             folded through the same response matrix
 **                             R as the resonances (FoldContinuum);
 **                             free parameters A_bg, a, T.
 **              gBgMode == 2 : pure Maxwellian truth-level continuum
 **                             f_bg(Erel) = sqrt(Erel) * exp(-Erel/T),
 **                             folded through the same response matrix
 **                             R as the resonances (FoldMaxwell);
 **                             free parameters A_bg, T only (shape power
 **                             fixed at 1/2).
 **
 **  All Erel <-> Eexc conversions use Eexc = Erel + Sn (gSn, default
 **  4.2 MeV, i.e. exactly your KinTree->Draw expression). The response
 **  and background templates are filled directly onto the data/true
 **  binning from the tree, so no bin-overlap rebinning is needed
 **  anywhere.
 **
 **  Fit: binned Poisson likelihood (Minuit2/Migrad + Hesse).
 **
 **  NOTE ON THE FIT RANGE: the response matrix only covers the Erel_true
 **  range found in the simFile tree, i.e. Eexc < Sn + trueHi (trueHi is
 **  auto-detected and printed at run time). Do not extend fitHi beyond
 **  that; default range is [4.3, 13.8] MeV (also below S3n).
 **
 **  Usage:
 **    root -l 'fitExc.C()'                       // default paths below
 **    root -l 'fitExc.C("data.root","bg.root","sim.root")'
 **
 **  Tune number of resonances / initial values / l assignments in the
 **  "configuration" block right below.
 **
 **/

#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TAxis.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TLegend.h"
#include "TLine.h"
#include "TLatex.h"
#include "TStyle.h"
#include "TF1.h"
#include "TMath.h"
#include "TString.h"
#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "Math/Functor.h"
#include "TRandom.h"
#include <vector>
#include <functional>
#include <algorithm>
#include <iostream>
#include <iomanip>

// ----------------------------------------------------------------------------
// Configuration (edit here)
// ----------------------------------------------------------------------------
static const Int_t kMaxRes = 6;

static Int_t gNRes = 3; // number of BW resonances

static Int_t gLorb[kMaxRes] = {2, 1, 1, 0, 0, 0};

static Double_t gInitE[kMaxRes] = {0.60, 3.20, 6.20, 0., 0., 0.};
static Double_t gInitG[kMaxRes] = {0.05, 1.20, 1.40, 0., 0., 0.};

static Double_t gElo[kMaxRes] = {0.30, 2.50, 5.00, 0., 0., 0.};
static Double_t gEhi[kMaxRes] = {0.90, 3.90, 7.50, 0., 0., 0.};
static Double_t gGlim[2] = {0.01, 3.0};

// Per-resonance parameter fixing, to freeze a value (e.g. to a literature
// number) without touching the fit code below.
static Bool_t gFixGamma[kMaxRes] = {kTRUE, kFALSE, kFALSE, kFALSE, kFALSE, kFALSE};
static Bool_t gFixE[kMaxRes] = {kFALSE, kFALSE, kFALSE, kFALSE, kFALSE, kFALSE};

static Bool_t gEnergyDepWidth = kTRUE;

// Penetrability model for the energy-dependent width:
//   0 = low-energy approx  Gamma(E) = G0*(E/E0)^(l+1/2)  (old)
//   1 = exact neutron penetrability (closed-form Bessel):
//       Gamma(E) = G0 * P_l(E)/P_l(E0)
static Int_t gPenMode = 1;

// channel radius R = r0*(Af^(1/3) + 1) [fm] and masses [MeV], used when
// gPenMode == 1
static Double_t gR0Fm = 1.4;
static Double_t gAFrag = 23.;
static Double_t gMassNeutronMeV = 939.56542;
static Double_t gMassFragMeV = 23. * 931.494 + 14.6; // ~23O; only enters via
                                                     // the reduced mass, MeV precision is more than sufficient

// Sub-sampling of the BW inside each 0.1 MeV truth bin (narrow states!)
static Int_t gNSub = 10;

// Background mode: 0 = template simulada INCL (tree "ana" de bgFile)
//                  1 = continuo paramétrico E^a*exp(-E/T), doblado
//                  2 = maxwelliana pura sqrt(E)*exp(-E/T), doblada
//                  3 = event mixing (data-driven): template en Erel
//                      reconstruido, del tree "tMix" de gMixFile.
//                      NO se dobla por la matriz de respuesta
//                      (resolución/aceptancia ya incluidas por
//                      construcción).
static Int_t gBgMode = 3;

// Fichero de event mixing (gBgMode == 3), tree "tMix" con ramas
// Erel (MeV, ya convertido, sin centinela) y weight.
static const char *gMixFile =
    "/nucl_lustre/pablogrusell/g249/g249_analysis/final/macros/erel/23O_mixing_iterative.root";

// límites e inicialización del continuo paramétrico (gBgMode == 1)
static Double_t gBgAInit = 0.5, gBgALim[2] = {0.05, 3.0};
static Double_t gBgTInit = 3.0, gBgTLim[2] = {0.5, 12.0};

// límites e inicialización de la maxwelliana (gBgMode == 2)
static Double_t gBgMaxwellTInit = 3.0, gBgMaxwellTLim[2] = {0.5, 12.0};

// Toy MC (parametric bootstrap): refit gNToys pseudo-data realizations
// drawn from the nominal model to get MC-based uncertainties/bias/
// correlations. 0 = disabled. Can be overridden by the nToys argument
// of fitErel().
static Int_t gNToys = 500;
static Bool_t gToyFluctResp = kFALSE; // also fluctuate the response matrix
static UInt_t gToySeed = 12345;

// Threshold diagnostics (breakdown table + data/sim edge calibration),
// run right after the nominal fit. Set to kFALSE to skip them (e.g. during
// the ScanLPeak1() systematic scan, where only the summary table matters).
static Bool_t gThresholdDiag = kTRUE;

// l override for resonance 1 (gLorb[0]), used by ScanLPeak1() below;
// -1 = no override (use gLorb[0] as configured above).
static Int_t gLPeak1Override = -1;

// Optional s-wave threshold (virtual state) component:
//   f_v(E) = sqrt(E) / (E + Ev)^2   (scattering-length form)
// Added ON TOP of the gNRes BW resonances. Parameters: Av
// (amplitude, same meaning as the A_k) and Ev [MeV].
static Bool_t gVirtualState = kFALSE;
static Double_t gEvInit = 0.3;
static Double_t gEvLim[2] = {0.05, 2.0};
// Automatic significance test: if kTRUE, after the nominal fit
// re-run the fit with the virtual component disabled and print
// 2*DeltaNLL.
static Bool_t gVirtualSignif = kTRUE;

// ----------------------------------------------------------------------------
// Globals shared with the likelihood function
// ----------------------------------------------------------------------------
static Double_t gSn = 4.2;
static Double_t gFitLo = 4.2, gFitHi = 17.;

// Truth-axis binning for the response matrix (Erel_true, MeV); the range
// itself is auto-detected from the simFile tree at runtime (see fitErel()).
static Int_t gNTrueBins = 200;
static TH1D *gData = nullptr;  // data, Eexc binning
static TH1D *gBgExc = nullptr; // bg template, Eexc binning, unit integral in fit range
static TH2D *gRespN = nullptr; // normalized response (X=true Erel, Y=rec Erel)

// Raw (unnormalized) response counts and per-truth-column normalization
// denominator, saved when gRespN is built (section 2) so toy MC can
// re-fluctuate the response matrix (gToyFluctResp) without re-reading
// simFile. Index 1..nT, same as gRespN's X axis.
static TH2D *gRespRaw = nullptr;
static std::vector<Double_t> gRespDen;

// Summary results from the most recent fitErel() call, stashed here so
// ScanLPeak1() (a separate function, called after fitErel() returns) can
// read them off without fitErel() having to return a struct.
static Double_t gLastFval = 0.;
static Double_t gLastChi2 = 0.;
static Int_t gLastNdf = 0;
static Double_t gLastErel1 = 0., gLastA1 = 0.;
static Double_t gLastModelBin1 = 0., gLastDataBin1 = 0.;

// ----------------------------------------------------------------------------
// Dimensionless wavenumber rho = k*R, k = sqrt(2*mu*E)/hbar-c (E in MeV,
// R in fm), used by NeutronPenetrability below.
// ----------------------------------------------------------------------------
static Double_t NeutronRho(Double_t E)
{
    static const Double_t hbarc = 197.3269804; // MeV fm
    const Double_t mu = gMassNeutronMeV * gMassFragMeV /
                        (gMassNeutronMeV + gMassFragMeV);
    const Double_t R = gR0Fm * (TMath::Power(gAFrag, 1. / 3.) + 1.);
    return TMath::Sqrt(2. * mu * E) / hbarc * R;
}

// ----------------------------------------------------------------------------
// Exact single-neutron penetrability P_l(rho) (closed-form spherical Bessel
// functions), l = 0..3 -- covers all resonances used in this analysis.
// ----------------------------------------------------------------------------
static Double_t NeutronPenetrability(Double_t E, Int_t l)
{
    if (E <= 0.)
        return 0.;
    const Double_t r = NeutronRho(E);
    const Double_t r2 = r * r;
    switch (l)
    {
    case 0:
        return r;
    case 1:
        return r * r2 / (1. + r2);
    case 2:
        return r * r2 * r2 / (9. + 3. * r2 + r2 * r2);
    case 3:
        return r * r2 * r2 * r2 /
               (225. + 45. * r2 + 6. * r2 * r2 + r2 * r2 * r2);
    default:
        // fallback: low-energy power law (should not happen for l<=3,
        // which covers all our cases)
        return TMath::Power(r, 2 * l + 1);
    }
}

// ----------------------------------------------------------------------------
// Breit-Wigner line shape (not normalized; normalization done numerically)
// ----------------------------------------------------------------------------
static Double_t BWShape(Double_t E, Double_t E0, Double_t G0, Int_t l)
{
    if (E <= 0. || E0 <= 0. || G0 <= 0.)
        return 0.;
    Double_t G = G0;
    if (gEnergyDepWidth)
    {
        if (gPenMode == 1)
        {
            const Double_t p0 = NeutronPenetrability(E0, l);
            G = (p0 > 0.) ? G0 * NeutronPenetrability(E, l) / p0 : G0;
        }
        else
            G = G0 * TMath::Power(E / E0, l + 0.5);
    }
    const Double_t d = E - E0;
    return G / (d * d + 0.25 * G * G);
}

// ----------------------------------------------------------------------------
// Non-resonant continuum line shape (not normalized), used when
// gBgMode == 1: f(E) = E^a * exp(-E/T)
// ----------------------------------------------------------------------------
static Double_t ContinuumShape(Double_t E, Double_t a, Double_t T)
{
    if (E <= 0. || T <= 0.)
        return 0.;
    return TMath::Power(E, a) * TMath::Exp(-E / T);
}

// ----------------------------------------------------------------------------
// Pure Maxwellian continuum line shape (not normalized), used when
// gBgMode == 2: f(E) = sqrt(E) * exp(-E/T)  (shape power fixed at 1/2,
// i.e. ContinuumShape with a=0.5 -- only the amplitude and T are free)
// ----------------------------------------------------------------------------
static Double_t MaxwellShape(Double_t E, Double_t T)
{
    if (E <= 0. || T <= 0.)
        return 0.;
    return TMath::Sqrt(E) * TMath::Exp(-E / T);
}

// ----------------------------------------------------------------------------
// Fold a truth-level line shape (amplitude A = generated decays/events in
// the truth window covered by gRespN's X axis) through the response matrix
// and add onto 'out'. The Y axis of gRespN is already Eexc with the data
// binning, so the folded spectrum can be added bin-by-bin (index j)
// directly onto 'out'. Shared by FoldResonance and FoldContinuum below.
// ----------------------------------------------------------------------------
static void FoldShape(Double_t A, TH1D *out,
                      const std::function<Double_t(Double_t)> &shape)
{
    const TAxis *axT = gRespN->GetXaxis();
    const TAxis *axR = gRespN->GetYaxis();
    const Int_t nT = axT->GetNbins();
    const Int_t nR = axR->GetNbins();

    // shape weights on the truth axis, normalized to unity
    std::vector<Double_t> w(nT + 1, 0.);
    Double_t wsum = 0.;
    for (Int_t i = 1; i <= nT; ++i)
    {
        const Double_t lo = axT->GetBinLowEdge(i);
        const Double_t hi = axT->GetBinUpEdge(i);
        Double_t s = 0.;
        for (Int_t k = 0; k < gNSub; ++k)
            s += shape(lo + (k + 0.5) * (hi - lo) / gNSub);
        w[i] = s * (hi - lo) / gNSub;
        wsum += w[i];
    }
    if (wsum <= 0.)
        return;

    // Fold: reconstructed-Erel spectrum
    std::vector<Double_t> rec(nR + 1, 0.);
    for (Int_t i = 1; i <= nT; ++i)
    {
        const Double_t wi = w[i] / wsum;
        if (wi == 0.)
            continue;
        for (Int_t j = 1; j <= nR; ++j)
        {
            const Double_t r = gRespN->GetBinContent(i, j);
            if (r > 0.)
                rec[j] += wi * r;
        }
    }

    // gRespN's Y axis is already Eexc with the data binning: add directly
    for (Int_t j = 1; j <= nR; ++j)
        if (rec[j] != 0.)
            out->AddBinContent(j, A * rec[j]);
}

static void FoldResonance(Double_t A, Double_t E0, Double_t G0, Int_t l, TH1D *out)
{
    FoldShape(A, out, [E0, G0, l](Double_t E)
              { return BWShape(E, E0, G0, l); });
}

// ----------------------------------------------------------------------------
// Fold the parametric continuum E^a*exp(-E/T) (amplitude A = generated
// events in the truth window covered by gRespN's X axis) through the
// response matrix and add onto 'out'. Used when gBgMode == 1.
// ----------------------------------------------------------------------------
static void FoldContinuum(Double_t A, Double_t a, Double_t T, TH1D *out)
{
    FoldShape(A, out, [a, T](Double_t E)
              { return ContinuumShape(E, a, T); });
}

// ----------------------------------------------------------------------------
// Fold the pure Maxwellian sqrt(E)*exp(-E/T) (amplitude A = generated events
// in the truth window covered by gRespN's X axis) through the response
// matrix and add onto 'out'. Used when gBgMode == 2.
// ----------------------------------------------------------------------------
static void FoldMaxwell(Double_t A, Double_t T, TH1D *out)
{
    FoldShape(A, out, [T](Double_t E)
              { return MaxwellShape(E, T); });
}

// ----------------------------------------------------------------------------
// Optional s-wave threshold ("virtual state") line shape (not normalized):
// f_v(E) = sqrt(E) / (E + Ev)^2 -- the standard scattering-length form for
// an s-wave virtual state right at threshold. Confined to low Erel by
// construction (falls off as E^-3/2 at large E), so it does not compete
// with the non-resonant background away from threshold. Used when
// gVirtualState == kTRUE, on top of the gNRes BW resonances.
// ----------------------------------------------------------------------------
static Double_t VirtualShape(Double_t E, Double_t Ev)
{
    if (E <= 0. || Ev <= 0.)
        return 0.;
    const Double_t s = E + Ev;
    return TMath::Sqrt(E) / (s * s);
}

// ----------------------------------------------------------------------------
// Fold the virtual-state shape (amplitude A = generated events in the truth
// window covered by gRespN's X axis) through the response matrix and add
// onto 'out'. Used when gVirtualState == kTRUE.
// ----------------------------------------------------------------------------
static void FoldVirtual(Double_t A, Double_t Ev, TH1D *out)
{
    FoldShape(A, out, [Ev](Double_t E)
              { return VirtualShape(E, Ev); });
}

// ----------------------------------------------------------------------------
// Parameter-vector layout: p[3k]=A_k, p[3k+1]=E0_k, p[3k+2]=Gamma0_k for
// k=0..gNRes-1, then the background parameters (1-3 of them depending on
// gBgMode, starting at index 3*gNRes), then -- if gVirtualState -- the
// virtual-state amplitude/energy at Av=p[NBasePar()], Ev=p[NBasePar()+1].
// All inline nPar computations elsewhere use these two helpers instead of
// recomputing the ternary chain, so the virtual-state parameters are never
// missed.
// ----------------------------------------------------------------------------
static Int_t NBasePar() // parameter count without the virtual-state component
{
    return (gBgMode == 1)   ? 3 * gNRes + 3
           : (gBgMode == 2) ? 3 * gNRes + 2
                            : 3 * gNRes + 1;
}
static Int_t NPar() { return NBasePar() + (gVirtualState ? 2 : 0); }

// ----------------------------------------------------------------------------
// Build the full model on the data binning.
// Parameters: p[3k]=A_k, p[3k+1]=E0_k (Erel), p[3k+2]=Gamma0_k, p[3*nRes]=A_bg
// gBgMode == 1 adds two more: p[3*nRes+1]=a, p[3*nRes+2]=T (continuum shape)
// gBgMode == 2 adds one more:  p[3*nRes+1]=T (Maxwellian shape)
// gVirtualState adds two more, at the end: p[NBasePar()]=Av, p[NBasePar()+1]=Ev
// ----------------------------------------------------------------------------
static void BuildModel(const Double_t *p, TH1D *hM)
{
    hM->Reset();
    for (Int_t k = 0; k < gNRes; ++k)
        FoldResonance(p[3 * k], p[3 * k + 1], p[3 * k + 2], gLorb[k], hM);
    if (gBgMode == 1)
        FoldContinuum(p[3 * gNRes], p[3 * gNRes + 1], p[3 * gNRes + 2], hM);
    else if (gBgMode == 2)
        FoldMaxwell(p[3 * gNRes], p[3 * gNRes + 1], hM);
    else
        hM->Add(gBgExc, p[3 * gNRes]);
    if (gVirtualState)
        FoldVirtual(p[NBasePar()], p[NBasePar() + 1], hM);
}

// ----------------------------------------------------------------------------
// Binned Poisson negative log-likelihood (constant terms dropped)
// ----------------------------------------------------------------------------
static Double_t NLL(const Double_t *p)
{
    static TH1D *hM = nullptr;
    if (!hM)
    {
        hM = static_cast<TH1D *>(gData->Clone("hModel_tmp"));
        hM->SetDirectory(nullptr);
    }
    BuildModel(p, hM);

    const Int_t b1 = gData->FindFixBin(gFitLo + 1.e-9);
    const Int_t b2 = gData->FindFixBin(gFitHi - 1.e-9);
    Double_t nll = 0.;
    for (Int_t b = b1; b <= b2; ++b)
    {
        Double_t m = hM->GetBinContent(b);
        if (m < 1.e-9)
            m = 1.e-9;
        const Double_t n = gData->GetBinContent(b);
        nll += m - n * TMath::Log(m);
    }
    return nll;
}

// ----------------------------------------------------------------------------
// Book the nPar parameters and run Migrad (+ optionally Hesse) once on
// whatever 'gData' currently points to. If 'start' is non-empty, its values
// are used as starting points instead of the gInitE/gInitG/nDataFit-based
// defaults -- used to restart each toy MC fit from the nominal best-fit
// point. Returns the minimizer (ownership passed to the caller, who must
// delete it) and the Minimize() convergence flag in 'ok'. Returns nullptr
// (with ok = kFALSE) if Minuit2 is not available.
// ----------------------------------------------------------------------------
static ROOT::Math::Minimizer *
RunFit(const std::vector<Double_t> &start, Double_t nDataFit,
       Int_t printLevel, Int_t strategy, Bool_t &ok, Bool_t doHesse = kTRUE)
{
    const Int_t nPar = NPar();

    ROOT::Math::Minimizer *min =
        ROOT::Math::Factory::CreateMinimizer("Minuit2", "Migrad");
    if (!min)
    {
        std::cerr << "ERROR: Minuit2 not available\n";
        ok = kFALSE;
        return nullptr;
    }
    min->SetMaxFunctionCalls(200000);
    min->SetTolerance(0.1);
    min->SetStrategy(strategy);
    min->SetPrintLevel(printLevel);
    min->SetErrorDef(0.5); // NLL

    ROOT::Math::Functor fcn(&NLL, nPar);
    min->SetFunction(fcn);

    for (Int_t k = 0; k < gNRes; ++k)
    {
        const Double_t aInit = start.empty() ? 0.25 * nDataFit : start[3 * k];
        const Double_t eInit = start.empty() ? gInitE[k] : start[3 * k + 1];
        const Double_t gInit = start.empty() ? gInitG[k] : start[3 * k + 2];

        min->SetLimitedVariable(3 * k, Form("A%d", k + 1),
                                aInit, 0.01 * nDataFit,
                                0., 20. * nDataFit);
        if (gFixE[k])
            min->SetFixedVariable(3 * k + 1, Form("Erel%d", k + 1), gInitE[k]);
        else
            min->SetLimitedVariable(3 * k + 1, Form("Erel%d", k + 1),
                                    eInit, 0.02, gElo[k], gEhi[k]);
        if (gFixGamma[k])
            min->SetFixedVariable(3 * k + 2, Form("Gamma%d", k + 1), gInitG[k]);
        else
            min->SetLimitedVariable(3 * k + 2, Form("Gamma%d", k + 1),
                                    gInit, 0.02, gGlim[0], gGlim[1]);
    }

    const Double_t bgInit = start.empty() ? 0.3 * nDataFit : start[3 * gNRes];
    min->SetLimitedVariable(3 * gNRes, "Abg",
                            bgInit, 0.01 * nDataFit,
                            0., 20. * nDataFit);
    if (gBgMode == 1)
    {
        const Double_t aInit2 = start.empty() ? gBgAInit : start[3 * gNRes + 1];
        const Double_t tInit2 = start.empty() ? gBgTInit : start[3 * gNRes + 2];
        min->SetLimitedVariable(3 * gNRes + 1, "bgA", aInit2, 0.05,
                                gBgALim[0], gBgALim[1]);
        min->SetLimitedVariable(3 * gNRes + 2, "bgT", tInit2, 0.1,
                                gBgTLim[0], gBgTLim[1]);
    }
    else if (gBgMode == 2)
    {
        const Double_t tInit2 = start.empty() ? gBgMaxwellTInit : start[3 * gNRes + 1];
        min->SetLimitedVariable(3 * gNRes + 1, "bgT", tInit2, 0.1,
                                gBgMaxwellTLim[0], gBgMaxwellTLim[1]);
    }

    if (gVirtualState)
    {
        const Double_t avInit = start.empty() ? 0.05 * nDataFit : start[NBasePar()];
        const Double_t evInit = start.empty() ? gEvInit : start[NBasePar() + 1];
        min->SetLimitedVariable(NBasePar(), "Av",
                                avInit, 0.01 * nDataFit,
                                0., 20. * nDataFit);
        min->SetLimitedVariable(NBasePar() + 1, "Ev",
                                evInit, 0.02, gEvLim[0], gEvLim[1]);
    }

    ok = min->Minimize();
    if (doHesse)
        min->Hesse();
    return min;
}

// ----------------------------------------------------------------------------
// Fit an error-function edge in [Sn-0.5, Sn+0.9]:
//     f(x) = p0 * 0.5 * (1 + Erf((x-p1)/(sqrt(2)*p2)))
// Used to locate the reaction threshold in a given spectrum (data or
// simulated/reconstructed), for a data-vs-simulation calibration check.
// Seeds: p0 = mean content in [Sn+0.6, Sn+0.9] (plateau), p1 = Sn, p2 = 0.15.
// Prints tag, position +- error, sigma +- error; fills pos/posErr/sigma/
// sigmaErr for the caller.
// ----------------------------------------------------------------------------
static void FitEdge(TH1D *h, const char *tag, Double_t &pos, Double_t &posErr,
                    Double_t &sigma, Double_t &sigmaErr)
{
    pos = posErr = sigma = sigmaErr = 0.;
    if (!h)
        return;

    const Double_t lo = gSn - 0.5, hi = gSn + 0.9;
    auto *fEdge = new TF1(Form("fEdge_%s", tag),
                          "[0]*0.5*(1+TMath::Erf((x-[1])/(sqrt(2)*[2])))", lo, hi);

    const Int_t bPlateauLo = h->FindFixBin(gSn + 0.6 + 1.e-9);
    const Int_t bPlateauHi = h->FindFixBin(gSn + 0.9 - 1.e-9);
    const Double_t plateau =
        (bPlateauHi >= bPlateauLo)
            ? h->Integral(bPlateauLo, bPlateauHi) / (bPlateauHi - bPlateauLo + 1)
            : h->GetMaximum();

    fEdge->SetParameters(plateau, gSn, 0.15);
    h->Fit(fEdge, "QNRS");

    pos = fEdge->GetParameter(1);
    posErr = fEdge->GetParError(1);
    sigma = fEdge->GetParameter(2);
    sigmaErr = fEdge->GetParError(2);

    std::cout << "  [FitEdge] " << std::left << std::setw(12) << tag
              << std::right << "  pos = " << std::fixed << std::setprecision(4)
              << pos << " +- " << posErr << " MeV   sigma = " << sigma
              << " +- " << sigmaErr << " MeV\n";
}

// ----------------------------------------------------------------------------
// Main entry point
// ----------------------------------------------------------------------------
void fitErel(const char *dataFile =
                 "/nucl_lustre/pablogrusell/g249/g249_analysis/results/dataFiles/23O_analyzed.root",
             const char *bgFile =
                 "/nucl_lustre/pablogrusell/g249/g249_analysis/final/macros/erel/ana_results/full2_analysis.root",
             const char *simFile =
                 "/nucl_lustre/pablogrusell/g249/g249_analysis/final/macros/erel/ana_results/full2_analysis.root",
             const char *cut = "califa_opa > 1.25 && califa_opa < 1.65",
             Double_t Sn = 4.2,
             Double_t fitLo = 4.2,
             Double_t fitHi = 17,
             Bool_t useEfficiency = kTRUE,
             Bool_t smoothBg = kFALSE,
             const char *outFile = "fitExc_results.root",
             Int_t lPeak2 = -1,
             Int_t nToys = 0,
             Int_t printLevel = 1)
{
    gSn = Sn;
    gFitLo = fitLo;
    gFitHi = fitHi;

    if (lPeak2 >= 0)
    {
        gLorb[1] = lPeak2;
        std::cout << "[fitExc] lPeak2 override: gLorb[1] = " << gLorb[1]
                  << " (systematic-l variant)\n";
    }
    if (gLPeak1Override >= 0)
    {
        gLorb[0] = gLPeak1Override;
        std::cout << "[fitExc] lPeak1 override: gLorb[0] = " << gLorb[0]
                  << " (ScanLPeak1 systematic scan)\n";
    }

    // ------------------------------------------------------------------------
    // 1) Data: Eexc = Erel*1000 + Sn, same binning as your Draw command
    // ------------------------------------------------------------------------
    TFile *fD = TFile::Open(dataFile, "READ");
    if (!fD || fD->IsZombie())
    {
        std::cerr << "ERROR: cannot open " << dataFile << std::endl;
        return;
    }
    auto *tree = dynamic_cast<TTree *>(fD->Get("KinTree"));
    if (!tree)
    {
        std::cerr << "ERROR: no KinTree in " << dataFile << std::endl;
        return;
    }

    gData = new TH1D("hDataExc", "", 100, 0., 17.);
    tree->Draw(Form("Erel*1000.+%g>>hDataExc", gSn), cut, "goff");
    if (gData->GetEntries() == 0)
        std::cout << "[fitExc] WARNING: 0 entries after cut '" << cut << "'\n";
    gData->SetDirectory(nullptr);
    gData->SetTitle(Form("^{24}O* excitation energy;E_{exc} [MeV];counts / %.0f keV",
                         gData->GetBinWidth(1) * 1000.));
    fD->Close();

    // ------------------------------------------------------------------------
    // 2) Response matrix, built directly from the "ana" TTree and normalized
    //    per generated event (or per column). The true-Erel axis range is
    //    auto-detected from the tree itself.
    // ------------------------------------------------------------------------
    TFile *fS = TFile::Open(simFile, "READ");
    if (!fS || fS->IsZombie())
    {
        std::cerr << "ERROR: cannot open " << simFile << std::endl;
        return;
    }
    auto *simTree = dynamic_cast<TTree *>(fS->Get("ana"));
    if (!simTree)
    {
        std::cerr << "ERROR: no ana TTree in " << simFile << std::endl;
        return;
    }

    simTree->SetEstimate(simTree->GetEntries() + 1);
    const Long64_t nSel =
        simTree->Draw("ErelTrue*1000.", "ErelTrue > -998", "goff");
    if (nSel <= 0)
    {
        std::cerr << "ERROR: no ErelTrue > -998 entries in " << simFile << std::endl;
        return;
    }
    Double_t *v = simTree->GetV1();
    Double_t trueLo = TMath::MinElement(nSel, v);
    Double_t trueHi = TMath::MaxElement(nSel, v);
    trueLo = TMath::Floor(trueLo * 10.) / 10.;
    trueHi = TMath::Ceil(trueHi * 10.) / 10.;
    if (trueLo < 0.)
        trueLo = 0.;
    std::cout << "[fitExc] response truth range from tree: [" << trueLo
              << ", " << trueHi << "] MeV\n";

    gRespN = new TH2D("hRespN", "", gNTrueBins, trueLo, trueHi,
                      gData->GetNbinsX(), gData->GetXaxis()->GetXmin(),
                      gData->GetXaxis()->GetXmax());
    simTree->Draw(Form("Erel*1000.+%g : ErelTrue*1000. >> hRespN", gSn),
                  "Erel > -998 && ErelTrue > -998", "goff");
    gRespN->SetDirectory(nullptr);

    auto *hTrueCount = new TH1D("hTrueCount", "", gNTrueBins, trueLo, trueHi);
    simTree->Draw("ErelTrue*1000. >> hTrueCount", "ErelTrue > -998", "goff");
    hTrueCount->SetDirectory(nullptr);

    const TAxis *axT = gRespN->GetXaxis();
    const Int_t nT = axT->GetNbins();
    const Int_t nR = gRespN->GetNbinsY();

    // keep a copy of the raw (unnormalized) matrix + per-column denominators
    // for the toy MC response-fluctuation option (gToyFluctResp)
    gRespRaw = static_cast<TH2D *>(gRespN->Clone("hRespRaw"));
    gRespRaw->SetDirectory(nullptr);
    gRespDen.assign(nT + 1, 0.);

    Int_t nEmptyCols = 0;
    for (Int_t i = 1; i <= nT; ++i)
    {
        Double_t den; // normalization denominator for this truth column
        if (useEfficiency)
        {
            den = hTrueCount->GetBinContent(i);
        }
        else
        {
            den = 0.;
            for (Int_t j = 1; j <= nR; ++j)
                den += gRespN->GetBinContent(i, j);
        }
        gRespDen[i] = den;

        if (den <= 0.)
        {
            ++nEmptyCols;
            for (Int_t j = 1; j <= nR; ++j)
                gRespN->SetBinContent(i, j, 0.);
            continue;
        }
        for (Int_t j = 1; j <= nR; ++j)
            gRespN->SetBinContent(i, j, gRespN->GetBinContent(i, j) / den);
    }
    std::cout << "[fitExc] response: " << nT << " truth bins, "
              << nEmptyCols << " empty; normalization = "
              << (useEfficiency ? "per generated event (efficiency folded in)"
                                : "per column (flat efficiency assumed)")
              << std::endl;
    fS->Close();

    if (gFitHi > gSn + trueHi)
    {
        std::cout << "[fitExc] WARNING: fitHi=" << gFitHi
                  << " MeV exceeds the response-matrix coverage (Eexc < Sn+trueHi = "
                  << gSn + trueHi << " MeV). Clamping.\n";
        gFitHi = gSn + trueHi - 0.2;
    }

    const Double_t nDataFit =
        gData->Integral(gData->FindFixBin(gFitLo + 1.e-9),
                        gData->FindFixBin(gFitHi - 1.e-9));
    std::cout << "[fitExc] data counts in fit range [" << gFitLo << ","
              << gFitHi << "] MeV : " << nDataFit << std::endl;

    // ------------------------------------------------------------------------
    // 3) Background template, built directly from a TTree onto the data
    //    binning (only needed for gBgMode == 0 or 3; the parametric/
    //    Maxwellian continua need no template)
    // ------------------------------------------------------------------------
    if (gBgMode == 0)
    {
        TFile *fB = TFile::Open(bgFile, "READ");
        if (!fB || fB->IsZombie())
        {
            std::cerr << "ERROR: cannot open " << bgFile << std::endl;
            return;
        }
        auto *bgTree = dynamic_cast<TTree *>(fB->Get("ana"));
        if (!bgTree)
        {
            std::cerr << "ERROR: no ana TTree in " << bgFile << std::endl;
            return;
        }

        gBgExc = new TH1D("hBgExc", "", gData->GetNbinsX(),
                          gData->GetXaxis()->GetXmin(), gData->GetXaxis()->GetXmax());
        bgTree->Draw(Form("Erel*1000.+%g >> hBgExc", gSn), "Erel > -998", "goff");
        gBgExc->SetDirectory(nullptr);
        fB->Close();
        if (smoothBg)
            gBgExc->Smooth(2);

        // unit integral in the fit range -> A_bg = background counts in range
        const Double_t bgInt =
            gBgExc->Integral(gBgExc->FindFixBin(gFitLo + 1.e-9),
                             gBgExc->FindFixBin(gFitHi - 1.e-9));
        if (bgInt <= 0.)
        {
            std::cerr << "ERROR: background template empty in the fit range\n";
            return;
        }
        gBgExc->Scale(1. / bgInt);
    }
    else if (gBgMode == 3)
    {
        TFile *fMix = TFile::Open(gMixFile, "READ");
        if (!fMix || fMix->IsZombie())
        {
            std::cerr << "ERROR: cannot open " << gMixFile << std::endl;
            return;
        }
        auto *mixTree = dynamic_cast<TTree *>(fMix->Get("tMix"));
        if (!mixTree)
        {
            std::cerr << "ERROR: no tMix TTree in " << gMixFile << std::endl;
            return;
        }

        gBgExc = new TH1D("hBgExc", "", gData->GetNbinsX(),
                          gData->GetXaxis()->GetXmin(), gData->GetXaxis()->GetXmax());
        gBgExc->Sumw2();
        // tMix: Erel already in MeV (no *1000), no -999 sentinel; each
        // entry carries a "weight" branch used as the fill weight.
        mixTree->Draw(Form("Erel+%g >> hBgExc", gSn), "weight", "goff");
        gBgExc->SetDirectory(nullptr);
        const Long64_t nMix = mixTree->GetEntries();
        fMix->Close();
        if (smoothBg)
            gBgExc->Smooth(2);

        // unit integral in the fit range -> A_bg = background counts in range
        const Double_t bgInt =
            gBgExc->Integral(gBgExc->FindFixBin(gFitLo + 1.e-9),
                             gBgExc->FindFixBin(gFitHi - 1.e-9));
        if (bgInt <= 0.)
        {
            std::cerr << "ERROR: background template empty in the fit range\n";
            return;
        }
        gBgExc->Scale(1. / bgInt);
        std::cout << "[fitExc] background: event-mixing template (data-driven, "
                     "not folded), "
                  << nMix << " tMix entries used\n";
    }
    else if (gBgMode == 2)
    {
        std::cout << "[fitExc] background: pure Maxwellian sqrt(E) exp(-E/T), "
                     "folded through response\n";
    }
    else
    {
        std::cout << "[fitExc] background: parametric continuum E^a exp(-E/T), "
                     "folded through response\n";
    }

    // ------------------------------------------------------------------------
    // 4-5) Minimizer setup + nominal fit (via RunFit, shared with the toys)
    // ------------------------------------------------------------------------
    const Int_t nPar = NPar();

    std::cout << "\n[fitExc] Fitting " << gNRes
              << " BW resonance(s) + background, range [" << gFitLo << ", "
              << gFitHi << "] MeV ...\n";
    Bool_t ok = kFALSE;
    ROOT::Math::Minimizer *min =
        RunFit({}, nDataFit, printLevel, /*strategy=*/2, ok, /*doHesse=*/kTRUE);
    if (!min)
        return;

    const Double_t *p = min->X();
    const Double_t *e = min->Errors();

    // Goodness of fit (Pearson chi2, for reporting only)
    auto *hModel = static_cast<TH1D *>(gData->Clone("hModel"));
    hModel->SetDirectory(nullptr);
    BuildModel(p, hModel);
    Double_t chi2 = 0.;
    Int_t nBinsFit = 0;
    const Int_t b1 = gData->FindFixBin(gFitLo + 1.e-9);
    const Int_t b2 = gData->FindFixBin(gFitHi - 1.e-9);
    for (Int_t b = b1; b <= b2; ++b)
    {
        const Double_t m = hModel->GetBinContent(b);
        if (m > 0.)
        {
            const Double_t d = gData->GetBinContent(b) - m;
            chi2 += d * d / m;
            ++nBinsFit;
        }
    }
    const Int_t ndf = nBinsFit - nPar;

    // stash summary results for ScanLPeak1() (see end of file)
    gLastFval = min->MinValue();
    gLastChi2 = chi2;
    gLastNdf = ndf;
    gLastErel1 = p[1];
    gLastA1 = p[0];
    gLastModelBin1 = hModel->GetBinContent(b1);
    gLastDataBin1 = gData->GetBinContent(b1);

    // ------------------------------------------------------------------------
    // 6) Report
    // ------------------------------------------------------------------------
    std::cout << "\n===================== fitExc results =====================\n";
    std::cout << (ok ? "  Minimization converged.\n"
                     : "  WARNING: minimization did NOT converge!\n");
    std::cout << "  chi2/ndf (Pearson) = " << chi2 << " / " << ndf << " = "
              << (ndf > 0 ? chi2 / ndf : 0.) << "\n";
    std::cout << "  NLL at minimum      = " << std::fixed
              << std::setprecision(6) << min->MinValue() << "\n";
    if (gNRes >= 2)
        std::cout << "  l assignment used for peak 2: l = " << gLorb[1]
                  << (lPeak2 >= 0 ? " (override via lPeak2 argument)\n"
                                  : " (nominal)\n");
    std::cout << "  Penetrability model: "
              << (!gEnergyDepWidth
                      ? "n/a (fixed-width BW, gEnergyDepWidth=kFALSE)\n"
                  : gPenMode == 1
                      ? "exact neutron penetrability (closed-form Bessel)\n"
                      : "low-energy approx Gamma(E) = G0*(E/E0)^(l+1/2)\n");
    if (gEnergyDepWidth && gPenMode == 1)
    {
        const Double_t Rfm = gR0Fm * (TMath::Power(gAFrag, 1. / 3.) + 1.);
        std::cout << "    channel radius R = " << Rfm << " fm\n";
        for (Int_t k = 0; k < gNRes; ++k)
            std::cout << "    rho(E0_" << k + 1 << "=" << p[3 * k + 1]
                      << " MeV, l=" << gLorb[k]
                      << ") = " << NeutronRho(p[3 * k + 1]) << "\n";
    }
    std::cout << "\n";
    std::cout << std::fixed << std::setprecision(3);
    for (Int_t k = 0; k < gNRes; ++k)
    {
        std::cout << "  Resonance " << k + 1 << " (l=" << gLorb[k] << "):\n"
                  << "     Erel   = " << p[3 * k + 1] << " +- ";
        if (gFixE[k])
            std::cout << "(fixed)";
        else
            std::cout << e[3 * k + 1];
        std::cout << " MeV\n"
                  << "     Eexc   = " << p[3 * k + 1] + gSn << " +- ";
        if (gFixE[k])
            std::cout << "(fixed)";
        else
            std::cout << e[3 * k + 1];
        std::cout << " MeV   (Sn = " << gSn << ")\n"
                  << "     Gamma0 = " << p[3 * k + 2] << " +- ";
        if (gFixGamma[k])
            std::cout << "(fixed)";
        else
            std::cout << e[3 * k + 2];
        std::cout << " MeV\n"
                  << "     Yield  = " << p[3 * k] << " +- " << e[3 * k]
                  << (useEfficiency ? "  (efficiency-corrected decays, Erel 0-10 MeV)\n"
                                    : "  (NOT efficiency corrected)\n");
    }
    if (gVirtualState)
    {
        std::cout << "  Virtual state (s-wave threshold):\n"
                  << "     Ev    = " << p[NBasePar() + 1] << " +- "
                  << e[NBasePar() + 1] << " MeV\n"
                  << "     Yield = " << p[NBasePar()] << " +- " << e[NBasePar()]
                  << (useEfficiency ? "  (efficiency-corrected decays, Erel 0-10 MeV)\n"
                                    : "  (NOT efficiency corrected)\n");
    }
    if (gBgMode == 1)
    {
        std::cout << "  Continuum: A = " << p[3 * gNRes] << " +- " << e[3 * gNRes]
                  << " (generated events, eff-corrected)\n"
                  << "             a = " << p[3 * gNRes + 1] << " +- "
                  << e[3 * gNRes + 1] << "\n"
                  << "             T = " << p[3 * gNRes + 2] << " +- "
                  << e[3 * gNRes + 2] << " MeV\n";
    }
    else if (gBgMode == 2)
    {
        std::cout << "  Maxwellian: A = " << p[3 * gNRes] << " +- " << e[3 * gNRes]
                  << " (generated events, eff-corrected)\n"
                  << "              T = " << p[3 * gNRes + 1] << " +- "
                  << e[3 * gNRes + 1] << " MeV\n";
    }
    else
    {
        std::cout << "  Background counts in fit range = " << p[3 * gNRes]
                  << " +- " << e[3 * gNRes] << "\n";
        if (gBgMode == 3)
            std::cout << "  (event-mixing template)\n";
    }
    std::cout << "===========================================================\n\n";

    // ------------------------------------------------------------------------
    // 6a) Virtual-state significance test: nominal (with virtual state) fit
    //     vs a from-scratch refit with the virtual component disabled;
    //     2*DeltaNLL for the 2 extra degrees of freedom (Av, Ev).
    // ------------------------------------------------------------------------
    if (gVirtualState && gVirtualSignif)
    {
        const Double_t nllWith = min->MinValue();

        gVirtualState = kFALSE;
        Bool_t okNoV = kFALSE;
        ROOT::Math::Minimizer *minNoV =
            RunFit({}, nDataFit, /*printLevel=*/-1, /*strategy=*/2, okNoV,
                   /*doHesse=*/kFALSE);
        const Double_t nllWithout = minNoV ? minNoV->MinValue() : 0.;
        delete minNoV;
        gVirtualState = kTRUE;

        std::cout << "[fitExc] virtual-state significance:\n"
                  << "  NLL with    = " << std::fixed << std::setprecision(6)
                  << nllWith << "\n"
                  << "  NLL without = " << nllWithout << "\n"
                  << "  2*DeltaNLL  = " << 2. * (nllWithout - nllWith)
                  << "  (2 extra dof)\n\n";
    }

    // ------------------------------------------------------------------------
    // 6c) Threshold diagnostics (breakdown table + data/sim edge calibration
    //     + peak-1 l scan lives in ScanLPeak1() below). Always run right
    //     after the nominal fit, independent of the toy MC block below.
    // ------------------------------------------------------------------------
    if (gThresholdDiag)
    {
        // ---- 1) threshold breakdown table ----
        std::vector<TH1D *> diagComp(gNRes, nullptr);
        for (Int_t k = 0; k < gNRes; ++k)
        {
            diagComp[k] = static_cast<TH1D *>(gData->Clone(Form("hDiagRes%d", k + 1)));
            diagComp[k]->SetDirectory(nullptr);
            diagComp[k]->Reset();
            FoldResonance(p[3 * k], p[3 * k + 1], p[3 * k + 2], gLorb[k], diagComp[k]);
        }
        TH1D *diagVirt = nullptr;
        if (gVirtualState)
        {
            diagVirt = static_cast<TH1D *>(gData->Clone("hDiagVirt"));
            diagVirt->SetDirectory(nullptr);
            diagVirt->Reset();
            FoldVirtual(p[NBasePar()], p[NBasePar() + 1], diagVirt);
        }
        TH1D *diagBg = nullptr;
        if (gBgMode == 1)
        {
            diagBg = static_cast<TH1D *>(gData->Clone("hDiagBg"));
            diagBg->SetDirectory(nullptr);
            diagBg->Reset();
            FoldContinuum(p[3 * gNRes], p[3 * gNRes + 1], p[3 * gNRes + 2], diagBg);
        }
        else if (gBgMode == 2)
        {
            diagBg = static_cast<TH1D *>(gData->Clone("hDiagBg"));
            diagBg->SetDirectory(nullptr);
            diagBg->Reset();
            FoldMaxwell(p[3 * gNRes], p[3 * gNRes + 1], diagBg);
        }
        else
        {
            diagBg = static_cast<TH1D *>(gBgExc->Clone("hDiagBg"));
            diagBg->SetDirectory(nullptr);
            diagBg->Scale(p[3 * gNRes]);
        }

        std::cout << "---- threshold breakdown ----\n";
        std::cout << std::right << std::setw(7) << "Eexc" << std::setw(8) << "data";
        for (Int_t k = 0; k < gNRes; ++k)
            std::cout << std::setw(8) << Form("BW%d", k + 1);
        if (gVirtualState)
            std::cout << std::setw(8) << "virt";
        std::cout << std::setw(8) << "bg" << std::setw(8) << "model"
                  << std::setw(8) << "pull" << "\n";

        const Int_t nDiag = TMath::Min(5, b2 - b1 + 1);
        for (Int_t idx = 0; idx < nDiag; ++idx)
        {
            const Int_t b = b1 + idx;
            const Double_t eexc = gData->GetBinCenter(b);
            const Double_t data = gData->GetBinContent(b);
            const Double_t model = hModel->GetBinContent(b);
            const Double_t pull = (model > 0.) ? (data - model) / TMath::Sqrt(model) : 0.;

            std::cout << std::setw(7) << Form("%.2f", eexc)
                      << std::setw(8) << Form("%.1f", data);
            for (Int_t k = 0; k < gNRes; ++k)
                std::cout << std::setw(8) << Form("%.1f", diagComp[k]->GetBinContent(b));
            if (gVirtualState)
                std::cout << std::setw(8) << Form("%.1f", diagVirt->GetBinContent(b));
            std::cout << std::setw(8) << Form("%.1f", diagBg->GetBinContent(b))
                      << std::setw(8) << Form("%.1f", model)
                      << std::setw(8) << Form("%+.1f", pull) << "\n";
        }
        std::cout << "------------------------------\n\n";

        for (auto *h : diagComp)
            delete h;
        delete diagVirt;
        delete diagBg;

        // ---- 2) edge fit: data vs simulation (reconstructed) calibration ----
        Double_t posData = 0., errData = 0., sigData = 0., sigErrData = 0.;
        FitEdge(gData, "data", posData, errData, sigData, sigErrData);

        Double_t posSim = 0., errSim = 0., sigSim = 0., sigErrSim = 0.;
        TFile *fSedge = TFile::Open(simFile, "READ");
        if (!fSedge || fSedge->IsZombie())
        {
            std::cerr << "  WARNING: cannot reopen " << simFile
                      << " for the edge-fit diagnostic\n";
        }
        else
        {
            auto *simTreeEdge = dynamic_cast<TTree *>(fSedge->Get("ana"));
            if (!simTreeEdge)
            {
                std::cerr << "  WARNING: no ana TTree in " << simFile
                          << " for the edge-fit diagnostic\n";
            }
            else
            {
                auto *hSimRec = static_cast<TH1D *>(gData->Clone("hSimRec"));
                hSimRec->Reset();
                simTreeEdge->Draw(Form("Erel*1000.+%g >> hSimRec", gSn),
                                  "Erel > -998", "goff");
                hSimRec->SetDirectory(nullptr);
                FitEdge(hSimRec, "sim (reco)", posSim, errSim, sigSim, sigErrSim);
                delete hSimRec;
            }
            fSedge->Close();
        }

        const Double_t diffKeV = (posData - posSim) * 1000.;
        const Double_t diffErrKeV =
            TMath::Sqrt(errData * errData + errSim * errSim) * 1000.;
        std::cout << "  data-sim edge offset = " << std::fixed
                  << std::setprecision(1) << diffKeV << " +- " << diffErrKeV
                  << " keV (calibration offset)\n";
        if (diffErrKeV > 0. && TMath::Abs(diffKeV) > 2. * diffErrKeV)
            std::cout << "  WARNING: |data-sim edge offset| > 2 sigma -- check "
                         "the energy calibration between data and simulation!\n";
        std::cout << "\n";
    }

    // ------------------------------------------------------------------------
    // 6b) Toy MC (parametric bootstrap): refit nToysRun pseudo-data
    //     realizations drawn from the nominal model hModel, restarting each
    //     toy fit from the nominal best-fit point, to get MC-based
    //     uncertainties/bias/correlations. gToyFluctResp additionally
    //     re-fluctuates the response matrix (Poisson on the raw migration
    //     counts, fixed per-column denominator) for each toy.
    // ------------------------------------------------------------------------
    const Int_t nToysRun = (nToys >= 0) ? nToys : gNToys;

    std::vector<TString> parNames(nPar);
    for (Int_t i = 0; i < nPar; ++i)
        parNames[i] = min->VariableName(i);
    std::vector<Int_t> freeIdx;
    for (Int_t i = 0; i < nPar; ++i)
        if (!min->IsFixedVariable(i))
            freeIdx.push_back(i);
    const Int_t nFree = (Int_t)freeIdx.size();

    TTree *toyTree = nullptr;
    std::vector<TH1D *> toyMarginals;
    TH2D *toyCorrHist = nullptr;

    if (nToysRun > 0)
    {
        TH1D *dataReal = gData;
        const std::vector<Double_t> startNom(p, p + nPar);
        std::vector<std::vector<Double_t>> toyVals(nFree);

        toyTree = new TTree("toys", "Toy MC fit results");
        std::vector<Double_t> branchVal(nPar, 0.);
        Int_t toyStatus = 0;
        for (Int_t i = 0; i < nPar; ++i)
            toyTree->Branch(parNames[i], &branchVal[i]);
        toyTree->Branch("status", &toyStatus);

        // backup of the (already normalized) response matrix, restored
        // after each toy if gToyFluctResp perturbs it in place
        auto *respBackup = static_cast<TH2D *>(gRespN->Clone("hRespBackup"));
        respBackup->SetDirectory(nullptr);

        gRandom->SetSeed(gToySeed);

        Int_t nConverged = 0;
        for (Int_t t = 1; t <= nToysRun; ++t)
        {
            auto *toyData = static_cast<TH1D *>(dataReal->Clone("hToyData"));
            toyData->SetDirectory(nullptr);
            toyData->Reset();
            for (Int_t b = b1; b <= b2; ++b)
                toyData->SetBinContent(b, gRandom->Poisson(hModel->GetBinContent(b)));

            if (gToyFluctResp)
            {
                const Int_t nTr = gRespN->GetNbinsX();
                const Int_t nRr = gRespN->GetNbinsY();
                for (Int_t i = 1; i <= nTr; ++i)
                {
                    const Double_t den = gRespDen[i];
                    if (den <= 0.)
                        continue;
                    for (Int_t j = 1; j <= nRr; ++j)
                    {
                        const Double_t raw = gRespRaw->GetBinContent(i, j);
                        const Double_t fluct = (raw > 0.) ? gRandom->Poisson(raw) : 0.;
                        gRespN->SetBinContent(i, j, fluct / den);
                    }
                }
            }

            gData = toyData;
            Bool_t okToy = kFALSE;
            ROOT::Math::Minimizer *minToy =
                RunFit(startNom, nDataFit, /*printLevel=*/-1, /*strategy=*/1,
                       okToy, /*doHesse=*/kFALSE);

            if (gToyFluctResp)
                for (Int_t i = 1; i <= respBackup->GetNbinsX(); ++i)
                    for (Int_t j = 1; j <= respBackup->GetNbinsY(); ++j)
                        gRespN->SetBinContent(i, j, respBackup->GetBinContent(i, j));

            if (minToy)
            {
                const Int_t statusToy = minToy->Status();
                if (okToy && statusToy == 0)
                {
                    const Double_t *pToy = minToy->X();
                    for (Int_t i = 0; i < nPar; ++i)
                        branchVal[i] = pToy[i];
                    toyStatus = statusToy;
                    toyTree->Fill();
                    for (Int_t f = 0; f < nFree; ++f)
                        toyVals[f].push_back(pToy[freeIdx[f]]);
                    ++nConverged;
                }
                delete minToy;
            }
            delete toyData;

            if (t % 50 == 0)
                std::cout << "[fitExc] toys: " << t << " / " << nToysRun << "\n";
        }

        gData = dataReal;
        delete respBackup;

        std::cout << "[fitExc] toys converged: " << nConverged << " / "
                  << nToysRun << "\n";
        if (nConverged < 0.95 * nToysRun)
            std::cout << "[fitExc] WARNING: less than 95% of toys converged!\n";

        if (nConverged > 0)
        {
            std::cout << "\n----------------- toy MC summary -----------------\n";
            std::cout << std::left << std::setw(10) << "param"
                      << std::right << std::setw(12) << "nominal"
                      << std::setw(12) << "median"
                      << std::setw(10) << "-sigma"
                      << std::setw(10) << "+sigma"
                      << std::setw(10) << "bias" << "\n";
            std::cout << std::fixed << std::setprecision(4);

            for (Int_t f = 0; f < nFree; ++f)
            {
                const Int_t idx = freeIdx[f];
                std::vector<Double_t> vals = toyVals[f];
                const Int_t n = (Int_t)vals.size();

                Double_t vlo = *std::min_element(vals.begin(), vals.end());
                Double_t vhi = *std::max_element(vals.begin(), vals.end());
                const Double_t span = (vhi > vlo) ? (vhi - vlo) : 1.;
                vlo -= 0.1 * span;
                vhi += 0.1 * span;

                auto *hMarg = new TH1D(Form("toy_%s", parNames[idx].Data()),
                                       Form(";%s;toys", parNames[idx].Data()),
                                       60, vlo, vhi);
                hMarg->SetDirectory(nullptr);
                for (Double_t val : vals)
                    hMarg->Fill(val);
                toyMarginals.push_back(hMarg);

                Double_t prob[3] = {0.16, 0.50, 0.84};
                Double_t quant[3] = {0., 0., 0.};
                TMath::Quantiles(n, 3, vals.data(), quant, prob, kFALSE);

                std::cout << std::left << std::setw(10) << parNames[idx]
                          << std::right
                          << std::setw(12) << p[idx]
                          << std::setw(12) << quant[1]
                          << std::setw(10) << (quant[1] - quant[0])
                          << std::setw(10) << (quant[2] - quant[1])
                          << std::setw(10) << (quant[1] - p[idx]) << "\n";
            }

            // Pearson correlation matrix of the toy parameter values
            std::vector<Double_t> mean(nFree, 0.), sigma(nFree, 0.);
            for (Int_t f = 0; f < nFree; ++f)
            {
                Double_t s = 0.;
                for (Double_t x : toyVals[f])
                    s += x;
                mean[f] = s / toyVals[f].size();
                Double_t s2 = 0.;
                for (Double_t x : toyVals[f])
                    s2 += (x - mean[f]) * (x - mean[f]);
                sigma[f] = TMath::Sqrt(s2 / toyVals[f].size());
            }

            std::vector<std::vector<Double_t>> corr(nFree, std::vector<Double_t>(nFree, 0.));
            Double_t maxAbsRho = -1.;
            Int_t bi = -1, bj = -1;
            for (Int_t fi = 0; fi < nFree; ++fi)
            {
                for (Int_t fj = 0; fj < nFree; ++fj)
                {
                    if (sigma[fi] <= 0. || sigma[fj] <= 0.)
                    {
                        corr[fi][fj] = (fi == fj) ? 1. : 0.;
                        continue;
                    }
                    Double_t cov = 0.;
                    for (size_t kk = 0; kk < toyVals[fi].size(); ++kk)
                        cov += (toyVals[fi][kk] - mean[fi]) * (toyVals[fj][kk] - mean[fj]);
                    cov /= toyVals[fi].size();
                    corr[fi][fj] = cov / (sigma[fi] * sigma[fj]);
                    if (fi != fj && TMath::Abs(corr[fi][fj]) > maxAbsRho)
                    {
                        maxAbsRho = TMath::Abs(corr[fi][fj]);
                        bi = fi;
                        bj = fj;
                    }
                }
            }

            std::cout << "\n  toy correlation matrix:\n       ";
            for (Int_t f = 0; f < nFree; ++f)
                std::cout << std::setw(8) << parNames[freeIdx[f]];
            std::cout << "\n";
            std::cout << std::setprecision(2);
            for (Int_t fi = 0; fi < nFree; ++fi)
            {
                std::cout << std::setw(7) << parNames[freeIdx[fi]];
                for (Int_t fj = 0; fj < nFree; ++fj)
                    std::cout << std::setw(8) << corr[fi][fj];
                std::cout << "\n";
            }
            std::cout << "---------------------------------------------------\n\n";

            if (bi >= 0 && bj >= 0)
            {
                toyCorrHist = new TH2D(
                    Form("toy_corr_%s_%s", parNames[freeIdx[bi]].Data(),
                         parNames[freeIdx[bj]].Data()),
                    Form(";%s;%s", parNames[freeIdx[bi]].Data(),
                         parNames[freeIdx[bj]].Data()),
                    toyMarginals[bi]->GetNbinsX(),
                    toyMarginals[bi]->GetXaxis()->GetXmin(),
                    toyMarginals[bi]->GetXaxis()->GetXmax(),
                    toyMarginals[bj]->GetNbinsX(),
                    toyMarginals[bj]->GetXaxis()->GetXmin(),
                    toyMarginals[bj]->GetXaxis()->GetXmax());
                toyCorrHist->SetDirectory(nullptr);
                for (size_t kk = 0; kk < toyVals[bi].size(); ++kk)
                    toyCorrHist->Fill(toyVals[bi][kk], toyVals[bj][kk]);
            }
        }
    }

    // ------------------------------------------------------------------------
    // 7) Draw -- two pads: main fit (top) + pull panel (bottom)
    // ------------------------------------------------------------------------
    gStyle->SetOptStat(0);
    gStyle->SetCanvasPreferGL();
    auto *c = new TCanvas("cFitExc", "24O* excitation energy fit", 950, 800);

    auto *padTop = new TPad("padTop", "padTop", 0., 0.30, 1., 1.00);
    padTop->SetLeftMargin(0.12);
    padTop->SetBottomMargin(0.02);
    padTop->SetTopMargin(0.08);
    padTop->Draw();

    auto *padBot = new TPad("padBot", "padBot", 0., 0.00, 1., 0.30);
    padBot->SetLeftMargin(0.12);
    padBot->SetTopMargin(0.02);
    padBot->SetBottomMargin(0.35);
    padBot->Draw();

    padTop->cd();

    gData->SetMarkerStyle(20);
    gData->SetMarkerSize(0.8);
    gData->SetLineColor(kBlack);
    gData->GetXaxis()->SetLabelSize(0.);
    gData->GetXaxis()->SetTitleSize(0.);
    gData->Draw("E1");

    // individual resonances (each on top of the background, for the eye)
    const Int_t resCol[kMaxRes] = {kOrange + 1, kMagenta + 1, kCyan + 2,
                                   kGreen + 2, kRed + 1, kBlue + 1};

    std::vector<TH1D *> comps;
    for (Int_t k = 0; k < gNRes; ++k)
    {
        auto *hk = static_cast<TH1D *>(gData->Clone(Form("hRes%d", k + 1)));
        hk->SetDirectory(nullptr);
        hk->Reset();
        FoldResonance(p[3 * k], p[3 * k + 1], p[3 * k + 2], gLorb[k], hk);
        // hk->Add(hBgComp);
        hk->SetLineColor(resCol[k]);
        hk->SetFillColorAlpha(resCol[k], 0.3);
        hk->SetLineWidth(2);
        hk->SetLineStyle(1);
        hk->Draw("HIST SAME");
        comps.push_back(hk);
    }

    // virtual-state component (drawn on its own, no background stacked)
    TH1D *hVirt = nullptr;
    if (gVirtualState)
    {
        hVirt = static_cast<TH1D *>(gData->Clone("hVirtualComponent"));
        hVirt->SetDirectory(nullptr);
        hVirt->Reset();
        FoldVirtual(p[NBasePar()], p[NBasePar() + 1], hVirt);
        hVirt->SetLineColor(kRed + 2);
        hVirt->SetFillColorAlpha(kRed + 2, 0.3);
        hVirt->SetLineWidth(2);
        hVirt->SetLineStyle(1);
        hVirt->Draw("HIST SAME");
    }

    hModel->SetLineColor(kRed);
    hModel->SetLineWidth(3);
    hModel->Draw("HIST SAME");
    gData->Draw("E1 SAME");

    // background component
    TH1D *hBgComp = nullptr;
    if (gBgMode == 1)
    {
        hBgComp = static_cast<TH1D *>(gData->Clone("hBgComp"));
        hBgComp->SetDirectory(nullptr);
        hBgComp->Reset();
        FoldContinuum(p[3 * gNRes], p[3 * gNRes + 1], p[3 * gNRes + 2], hBgComp);
    }
    else if (gBgMode == 2)
    {
        hBgComp = static_cast<TH1D *>(gData->Clone("hBgComp"));
        hBgComp->SetDirectory(nullptr);
        hBgComp->Reset();
        FoldMaxwell(p[3 * gNRes], p[3 * gNRes + 1], hBgComp);
    }
    else
    {
        hBgComp = static_cast<TH1D *>(gBgExc->Clone("hBgComp"));
        hBgComp->SetDirectory(nullptr);
        hBgComp->Scale(p[3 * gNRes]);
    }
    hBgComp->SetLineColor(kGray + 2);
    hBgComp->SetLineStyle(2);
    hBgComp->SetLineWidth(2);
    hBgComp->SetFillColorAlpha(kGray, 0.35);
    hBgComp->Draw("HIST SAME");

    auto *leg = new TLegend(0.58, 0.60, 0.93, 0.92);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->AddEntry(gData, "data", "pe");
    leg->AddEntry(hModel, "total fit", "l");
    for (Int_t k = 0; k < gNRes; ++k)
        leg->AddEntry(comps[k],
                      Form("E_{exc}=%.2f, #Gamma_{0}=%.2f MeV",
                           p[3 * k + 1] + gSn, p[3 * k + 2]),
                      "l");
    if (hVirt)
        leg->AddEntry(hVirt,
                      Form("virtual s-wave, E_{v}=%.2f MeV", p[NBasePar() + 1]),
                      "l");
    leg->AddEntry(hBgComp,
                  gBgMode == 1   ? "continuum E^{a}e^{-E/T}"
                  : gBgMode == 2 ? "Maxwellian #sqrt{E}e^{-E/T}"
                  : gBgMode == 3 ? "mixed-event bg"
                                 : "non-resonant bg",
                  "lf");
    leg->Draw();

    for (Double_t x : {gFitLo, gFitHi})
    {
        auto *l = new TLine(x, 0., x, gData->GetMaximum() * 0.35);
        l->SetLineStyle(3);
        l->SetLineColor(kBlack);
        l->Draw();
    }

    auto *txtChi2 = new TLatex();
    txtChi2->SetNDC();
    txtChi2->SetTextSize(0.035);
    txtChi2->DrawLatex(0.15, 0.93,
                       Form("#chi^{2}/ndf = %.1f / %d", chi2, ndf));

    // ------------------------------------------------------------------------
    // 7b) Pull panel: (data - model)/sqrt(model), only within [gFitLo,gFitHi]
    // ------------------------------------------------------------------------
    auto *hPulls = static_cast<TH1D *>(gData->Clone("hPulls"));
    hPulls->SetDirectory(nullptr);
    hPulls->Reset();
    hPulls->SetTitle(";E_{exc} [MeV];pull");

    Double_t maxAbsPull = 0.;
    for (Int_t b = b1; b <= b2; ++b)
    {
        const Double_t m = hModel->GetBinContent(b);
        if (m > 0.)
        {
            const Double_t pull = (gData->GetBinContent(b) - m) / TMath::Sqrt(m);
            hPulls->SetBinContent(b, pull);
            hPulls->SetBinError(b, 0.);
            if (TMath::Abs(pull) > maxAbsPull)
                maxAbsPull = TMath::Abs(pull);
        }
    }
    Double_t pullRange = maxAbsPull * 1.2;
    if (pullRange < 3.)
        pullRange = 3.;

    padBot->cd();
    hPulls->SetMarkerStyle(20);
    hPulls->SetMarkerSize(0.7);
    hPulls->SetMarkerColor(kBlack);
    hPulls->SetLineColor(kBlack);
    hPulls->GetYaxis()->SetRangeUser(-pullRange, pullRange);
    hPulls->GetXaxis()->SetTitle("E_{exc} [MeV]");
    hPulls->GetYaxis()->SetTitle("pull");
    hPulls->GetXaxis()->SetLabelSize(0.10);
    hPulls->GetXaxis()->SetTitleSize(0.10);
    hPulls->GetXaxis()->SetTitleOffset(1.3);
    hPulls->GetYaxis()->SetLabelSize(0.10);
    hPulls->GetYaxis()->SetTitleSize(0.10);
    hPulls->GetYaxis()->SetTitleOffset(0.5);
    hPulls->GetYaxis()->SetNdivisions(505);
    hPulls->Draw("P");

    auto *lZero = new TLine(0., 0., 17., 0.);
    lZero->SetLineColor(kGray + 1);
    lZero->SetLineStyle(2);
    lZero->Draw();

    auto *lPlus2 = new TLine(0., 1., 17., 1.);
    lPlus2->SetLineColor(kGray);
    lPlus2->SetLineStyle(3);
    lPlus2->Draw();

    auto *lMinus2 = new TLine(0., -1., 17., -1.);
    lMinus2->SetLineColor(kGray);
    lMinus2->SetLineStyle(3);
    lMinus2->Draw();

    hPulls->Draw("P SAME");

    c->cd();
    c->Update();

    // ------------------------------------------------------------------------
    // 8) Save
    // ------------------------------------------------------------------------
    TFile *fout = TFile::Open(outFile, "RECREATE");
    gData->Write("hData");
    hModel->Write("hModelTotal");
    hBgComp->Write("hBgComponent");
    for (size_t k = 0; k < comps.size(); ++k)
        comps[k]->Write(Form("hResonance%zu_plusBg", k + 1));
    if (hVirt)
        hVirt->Write("hVirtualComponent");
    gRespN->Write("hRespNorm");
    hPulls->Write("hPulls");
    c->Write("cFitExc");
    if (toyTree)
        toyTree->Write();
    for (auto *h : toyMarginals)
        h->Write();
    if (toyCorrHist)
        toyCorrHist->Write();
    fout->Close();
    std::cout << "[fitExc] output written to " << outFile << std::endl;
}

// ----------------------------------------------------------------------------
// ScanLPeak1: systematic scan of the l value assigned to resonance 1
// (gLorb[0], default l=2 -- see the header comment above BWShape). Runs the
// full nominal fitErel() fit three times, with gLorb[0] = 0, 1, 2, toys
// disabled and Migrad printLevel = -1 (only the final summary table below is
// printed). Same default arguments as fitErel() where applicable.
//
// Implementation note: inputs (data/response/background) are NOT loaded
// once and shared -- that would require splitting fitErel() into a
// LoadInputs() + Fit() pair, which is invasive given how much of fitErel()
// (background mode branching, response-matrix truth-range detection, toy MC,
// drawing) is entangled with the loading step. Instead this simply calls
// fitErel() three times, varying gLPeak1Override (checked near the top of
// fitErel(), same pattern as the existing lPeak2 argument) -- each call
// reloads its own inputs, which is slower but simple and safe.
// ----------------------------------------------------------------------------
void ScanLPeak1(const char *dataFile =
                    "/nucl_lustre/pablogrusell/g249/g249_analysis/results/dataFiles/23O_analyzed.root",
                const char *bgFile =
                    "/nucl_lustre/pablogrusell/g249/g249_analysis/final/macros/erel/ana_results/full2_analysis.root",
                const char *simFile =
                    "/nucl_lustre/pablogrusell/g249/g249_analysis/final/macros/erel/ana_results/full2_analysis.root",
                const char *cut = "califa_opa > 1.25 && califa_opa < 1.65",
                Double_t Sn = 4.2,
                Double_t fitLo = 4.2,
                Double_t fitHi = 17,
                Bool_t useEfficiency = kTRUE,
                Bool_t smoothBg = kFALSE)
{
    const Int_t lVals[3] = {0, 1, 2};
    Double_t fval[3], chi2v[3], erel1[3], a1[3], modelBin1[3], dataBin1[3];

    const Bool_t diagBackup = gThresholdDiag;
    gThresholdDiag = kFALSE; // keep the scan output to just the summary table

    for (Int_t i = 0; i < 3; ++i)
    {
        gLPeak1Override = lVals[i];
        fitErel(dataFile, bgFile, simFile, cut, Sn, fitLo, fitHi, useEfficiency,
                smoothBg, Form("fitExc_results_lpeak1_%d.root", lVals[i]),
                /*lPeak2=*/-1, /*nToys=*/0, /*printLevel=*/-1);
        fval[i] = gLastFval;
        chi2v[i] = gLastChi2;
        erel1[i] = gLastErel1;
        a1[i] = gLastA1;
        modelBin1[i] = gLastModelBin1;
        dataBin1[i] = gLastDataBin1;
    }

    gLPeak1Override = -1; // restore nominal (unoverridden) behavior
    gThresholdDiag = diagBackup;

    const Double_t fvalBest = *std::min_element(fval, fval + 3);

    std::cout << "\n============= ScanLPeak1: l(peak 1) systematic scan =============\n";
    std::cout << std::left << std::setw(4) << "l"
              << std::right << std::setw(12) << "FVAL"
              << std::setw(14) << "2*DeltaNLL"
              << std::setw(14) << "chi2(Pearson)"
              << std::setw(10) << "Erel1"
              << std::setw(12) << "A1"
              << std::setw(14) << "model(bin1)"
              << std::setw(12) << "data(bin1)" << "\n";
    std::cout << std::fixed << std::setprecision(4);
    for (Int_t i = 0; i < 3; ++i)
    {
        std::cout << std::left << std::setw(4) << lVals[i]
                  << std::right << std::setw(12) << fval[i]
                  << std::setw(14) << 2. * (fval[i] - fvalBest)
                  << std::setw(14) << chi2v[i]
                  << std::setw(10) << erel1[i]
                  << std::setw(12) << a1[i]
                  << std::setw(14) << modelBin1[i]
                  << std::setw(12) << dataBin1[i] << "\n";
    }
    std::cout << "===================================================================\n\n";
}