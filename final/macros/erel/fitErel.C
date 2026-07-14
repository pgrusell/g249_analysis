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
 **   * R    : response matrix hErelRecVsTrue from anaDelta.C
 **            (X = Erel_true 0..10 MeV, Y = Erel_rec -5..10 MeV).
 **            Each truth column is normalized per GENERATED event using
 **            hErelTrue from the same file, so the Erel-dependent
 **            (neutron) efficiency is folded in and the fitted A_k are
 **            efficiency-corrected yields. If hErelTrue is absent the
 **            columns are normalized to unity (pure line shape, flat
 **            efficiency assumed).
 **
 **   * BG   : non-resonant background template = hErel from the
 **            bg simulation (reconstructed Erel), shifted by Sn and
 **            re-binned onto the data binning; one free scale A_bg.
 **
 **  All Erel <-> Eexc conversions use Eexc = Erel + Sn (gSn, default
 **  4.2 MeV, i.e. exactly your KinTree->Draw expression). Bin-overlap
 **  weighting is used everywhere the binnings differ (0.15 MeV rec bins
 **  vs 0.2 MeV data bins), so no interpolation artefacts.
 **
 **  Fit: binned Poisson likelihood (Minuit2/Migrad + Hesse).
 **
 **  NOTE ON THE FIT RANGE: the response matrix only covers
 **  Erel_true < 10 MeV, i.e. Eexc < Sn + 10. Do not extend fitHi
 **  beyond that; default range is [4.3, 13.8] MeV (also below S3n).
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
#include "TMath.h"
#include "TString.h"
#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "Math/Functor.h"
#include <vector>
#include <iostream>
#include <iomanip>

// ----------------------------------------------------------------------------
// Configuration (edit here)
// ----------------------------------------------------------------------------
static const Int_t kMaxRes = 6;

static Int_t gNRes = 3; // number of BW resonances

// Orbital angular momentum of the emitted neutron, per resonance
// (enters the energy-dependent width). 24O* -> 23O(1/2+) + n:
//
// MODELO DE 3 RESONANCIAS. Asignacion nominal: pico 1 (4.8 MeV) l=2,
// pico 2 (7.4 MeV) l=2, pico 3 (10.4 MeV) l=1. El pico 2 contiene un
// doblete no resuelto (l=1 y l=2): el fit se repite con gLorb[1]=1 y la
// diferencia se reporta como incertidumbre sistematica (ver argumento
// lPeak2 de fitErel()).

static Int_t gLorb[kMaxRes] = {2, 1, 1, 0, 0, 0};

static Double_t gInitE[kMaxRes] = {0.60, 3.20, 6.20, 0., 0., 0.};
static Double_t gInitG[kMaxRes] = {0.05, 1.20, 1.40, 0., 0., 0.};

static Double_t gElo[kMaxRes] = {0.30, 2.50, 5.00, 0., 0., 0.};
static Double_t gEhi[kMaxRes] = {0.90, 3.90, 7.50, 0., 0., 0.};
static Double_t gGlim[2] = {0.01, 3.0};

// Per-resonance parameter fixing, to freeze a value (e.g. to a literature
// number) without touching the fit code below.
static Bool_t gFixGamma[kMaxRes] = {kFALSE, kFALSE, kFALSE, kFALSE, kFALSE, kFALSE};
static Bool_t gFixE[kMaxRes] = {kFALSE, kFALSE, kFALSE, kFALSE, kFALSE, kFALSE};

static Bool_t gEnergyDepWidth = kTRUE;

// Sub-sampling of the BW inside each 0.1 MeV truth bin (narrow states!)
static Int_t gNSub = 10;

// ----------------------------------------------------------------------------
// Globals shared with the likelihood function
// ----------------------------------------------------------------------------
static Double_t gSn = 4.2;
static Double_t gFitLo = 4.3, gFitHi = 20.;
static TH1D *gData = nullptr;  // data, Eexc binning
static TH1D *gBgExc = nullptr; // bg template, Eexc binning, unit integral in fit range
static TH2D *gRespN = nullptr; // normalized response (X=true Erel, Y=rec Erel)

// ----------------------------------------------------------------------------
// Breit-Wigner line shape (not normalized; normalization done numerically)
// ----------------------------------------------------------------------------
static Double_t BWShape(Double_t E, Double_t E0, Double_t G0, Int_t l)
{
    if (E <= 0. || E0 <= 0. || G0 <= 0.)
        return 0.;
    const Double_t G =
        gEnergyDepWidth ? G0 * TMath::Power(E / E0, l + 0.5) : G0;
    const Double_t d = E - E0;
    return G / (d * d + 0.25 * G * G);
}

// ----------------------------------------------------------------------------
// Add the interval [lo,hi] with total content c into 'out', distributing
// the content over the overlapping bins proportionally to the overlap.
// ----------------------------------------------------------------------------
static void AddOverlap(TH1D *out, Double_t lo, Double_t hi, Double_t c)
{
    if (c == 0. || hi <= lo)
        return;
    const TAxis *ax = out->GetXaxis();
    Int_t b1 = ax->FindFixBin(lo + 1.e-9);
    Int_t b2 = ax->FindFixBin(hi - 1.e-9);
    if (b1 < 1)
        b1 = 1;
    if (b2 > out->GetNbinsX())
        b2 = out->GetNbinsX();
    for (Int_t b = b1; b <= b2; ++b)
    {
        const Double_t ov =
            TMath::Min(hi, ax->GetBinUpEdge(b)) - TMath::Max(lo, ax->GetBinLowEdge(b));
        if (ov > 0.)
            out->AddBinContent(b, c * ov / (hi - lo));
    }
}

// ----------------------------------------------------------------------------
// Fold one BW resonance (amplitude A = generated decays in the truth window
// 0..10 MeV) through the response matrix, shift Erel_rec -> Eexc, and add
// onto 'out' (which has the data binning).
// ----------------------------------------------------------------------------
static void FoldResonance(Double_t A, Double_t E0, Double_t G0, Int_t l, TH1D *out)
{
    const TAxis *axT = gRespN->GetXaxis();
    const TAxis *axR = gRespN->GetYaxis();
    const Int_t nT = axT->GetNbins();
    const Int_t nR = axR->GetNbins();

    // BW weights on the truth axis, normalized to unity
    std::vector<Double_t> w(nT + 1, 0.);
    Double_t wsum = 0.;
    for (Int_t i = 1; i <= nT; ++i)
    {
        const Double_t lo = axT->GetBinLowEdge(i);
        const Double_t hi = axT->GetBinUpEdge(i);
        Double_t s = 0.;
        for (Int_t k = 0; k < gNSub; ++k)
            s += BWShape(lo + (k + 0.5) * (hi - lo) / gNSub, E0, G0, l);
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

    // Shift to Eexc and rebin onto the data axis
    for (Int_t j = 1; j <= nR; ++j)
        if (rec[j] != 0.)
            AddOverlap(out,
                       axR->GetBinLowEdge(j) + gSn,
                       axR->GetBinUpEdge(j) + gSn,
                       A * rec[j]);
}

// ----------------------------------------------------------------------------
// Build the full model on the data binning.
// Parameters: p[3k]=A_k, p[3k+1]=E0_k (Erel), p[3k+2]=Gamma0_k, p[3*nRes]=A_bg
// ----------------------------------------------------------------------------
static void BuildModel(const Double_t *p, TH1D *hM)
{
    hM->Reset();
    for (Int_t k = 0; k < gNRes; ++k)
        FoldResonance(p[3 * k], p[3 * k + 1], p[3 * k + 2], gLorb[k], hM);
    hM->Add(gBgExc, p[3 * gNRes]);
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
// Main entry point
// ----------------------------------------------------------------------------
void fitErel(const char *dataFile =
                 "/nucl_lustre/pablogrusell/g249/g249_analysis/results/dataFiles/23O_analyzed.root",
             const char *bgFile =
                 "/nucl_lustre/pablogrusell/g249/g249_analysis/final/macros/erel/ana_results/bg_analysis.root",
             const char *simFile =
                 "/nucl_lustre/pablogrusell/g249/g249_analysis/final/macros/erel/ana_results/full_analysis.root",
             const char *cut = "califa_opa > 1.25 && califa_opa < 1.65",
             Double_t Sn = 4.2,
             Double_t fitLo = 4.3,
             Double_t fitHi = 20,
             Bool_t useEfficiency = kTRUE,
             Bool_t smoothBg = kTRUE,
             const char *outFile = "fitExc_results.root",
             Int_t lPeak2 = -1)
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

    if (gFitHi > gSn + 20.)
    {
        std::cout << "[fitExc] WARNING: fitHi=" << gFitHi
                  << " MeV exceeds the response-matrix coverage (Eexc < Sn+10 = "
                  << gSn + 10. << " MeV). Clamping.\n";
        gFitHi = gSn + 10. - 0.2;
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

    gData = new TH1D("hDataExc", "", 100, 0., 20.);
    tree->Draw(Form("Erel*1000.+%g>>hDataExc", gSn), cut, "goff");
    if (gData->GetEntries() == 0)
        std::cout << "[fitExc] WARNING: 0 entries after cut '" << cut << "'\n";
    gData->SetDirectory(nullptr);
    gData->SetTitle(Form("^{24}O* excitation energy;E_{exc} [MeV];counts / %.0f keV",
                         gData->GetBinWidth(1) * 1000.));
    fD->Close();

    const Double_t nDataFit =
        gData->Integral(gData->FindFixBin(gFitLo + 1.e-9),
                        gData->FindFixBin(gFitHi - 1.e-9));
    std::cout << "[fitExc] data counts in fit range [" << gFitLo << ","
              << gFitHi << "] MeV : " << nDataFit << std::endl;

    // ------------------------------------------------------------------------
    // 2) Response matrix, normalized per generated event (or per column)
    // ------------------------------------------------------------------------
    TFile *fS = TFile::Open(simFile, "READ");
    if (!fS || fS->IsZombie())
    {
        std::cerr << "ERROR: cannot open " << simFile << std::endl;
        return;
    }
    auto *hResp = dynamic_cast<TH2D *>(fS->Get("hErelRecVsTrue"));
    if (!hResp)
    {
        std::cerr << "ERROR: no hErelRecVsTrue in " << simFile << std::endl;
        return;
    }
    auto *hTrue = dynamic_cast<TH1D *>(fS->Get("hErelTrue")); // may be null

    gRespN = static_cast<TH2D *>(hResp->Clone("hRespNorm"));
    gRespN->SetDirectory(nullptr);

    const TAxis *axT = gRespN->GetXaxis();
    const Int_t nT = axT->GetNbins();
    const Int_t nR = gRespN->GetNbinsY();
    const Double_t dEtrue = axT->GetBinWidth(1);

    Int_t nEmptyCols = 0;
    for (Int_t i = 1; i <= nT; ++i)
    {
        Double_t colsum = 0.;
        for (Int_t j = 1; j <= nR; ++j)
            colsum += gRespN->GetBinContent(i, j);

        Double_t den = colsum; // default: unit-normalized columns
        if (useEfficiency && hTrue)
        {
            // expected generated events in this truth bin, from hErelTrue
            // (different binning -> use the interpolated density)
            const Double_t Ec = axT->GetBinCenter(i);
            const Double_t dens = hTrue->Interpolate(Ec) / hTrue->GetBinWidth(1);
            den = dens * dEtrue;
        }

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
              << ((useEfficiency && hTrue) ? "per generated event (efficiency folded in)"
                                           : "per column (flat efficiency assumed)")
              << std::endl;
    if (useEfficiency && !hTrue)
        std::cout << "[fitExc] WARNING: hErelTrue not found -- fitted yields "
                     "are NOT efficiency corrected.\n";
    fS->Close();

    // ------------------------------------------------------------------------
    // 3) Background template: hErel (rec Erel) -> Eexc, data binning
    // ------------------------------------------------------------------------
    TFile *fB = TFile::Open(bgFile, "READ");
    if (!fB || fB->IsZombie())
    {
        std::cerr << "ERROR: cannot open " << bgFile << std::endl;
        return;
    }
    auto *hBgRel = dynamic_cast<TH1D *>(fB->Get("hErel"));
    if (!hBgRel)
    {
        std::cerr << "ERROR: no hErel in " << bgFile << std::endl;
        return;
    }
    auto *hBg = static_cast<TH1D *>(hBgRel->Clone("hBgRelC"));
    hBg->SetDirectory(nullptr);
    fB->Close();
    if (smoothBg)
        hBg->Smooth(2);

    gBgExc = static_cast<TH1D *>(gData->Clone("hBgExc"));
    gBgExc->SetDirectory(nullptr);
    gBgExc->Reset();
    for (Int_t j = 1; j <= hBg->GetNbinsX(); ++j)
        AddOverlap(gBgExc,
                   hBg->GetXaxis()->GetBinLowEdge(j) + gSn,
                   hBg->GetXaxis()->GetBinUpEdge(j) + gSn,
                   hBg->GetBinContent(j));
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

    // ------------------------------------------------------------------------
    // 4) Minimizer setup
    // ------------------------------------------------------------------------
    const Int_t nPar = 3 * gNRes + 1;
    ROOT::Math::Minimizer *min =
        ROOT::Math::Factory::CreateMinimizer("Minuit2", "Migrad");
    if (!min)
    {
        std::cerr << "ERROR: Minuit2 not available\n";
        return;
    }
    min->SetMaxFunctionCalls(200000);
    min->SetTolerance(0.1);
    min->SetStrategy(2);
    min->SetPrintLevel(1);
    min->SetErrorDef(0.5); // NLL

    ROOT::Math::Functor fcn(&NLL, nPar);
    min->SetFunction(fcn);

    for (Int_t k = 0; k < gNRes; ++k)
    {
        min->SetLimitedVariable(3 * k, Form("A%d", k + 1),
                                0.25 * nDataFit, 0.01 * nDataFit,
                                0., 20. * nDataFit);
        if (gFixE[k])
            min->SetFixedVariable(3 * k + 1, Form("Erel%d", k + 1), gInitE[k]);
        else
            min->SetLimitedVariable(3 * k + 1, Form("Erel%d", k + 1),
                                    gInitE[k], 0.02, gElo[k], gEhi[k]);
        if (gFixGamma[k])
            min->SetFixedVariable(3 * k + 2, Form("Gamma%d", k + 1), gInitG[k]);
        else
            min->SetLimitedVariable(3 * k + 2, Form("Gamma%d", k + 1),
                                    gInitG[k], 0.02, gGlim[0], gGlim[1]);
    }
    min->SetLimitedVariable(3 * gNRes, "Abg",
                            0.3 * nDataFit, 0.01 * nDataFit,
                            0., 5. * nDataFit);

    // ------------------------------------------------------------------------
    // 5) Fit
    // ------------------------------------------------------------------------
    std::cout << "\n[fitExc] Fitting " << gNRes
              << " BW resonance(s) + background, range [" << gFitLo << ", "
              << gFitHi << "] MeV ...\n";
    const Bool_t ok = min->Minimize();
    min->Hesse();

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
    std::cout << "  Background counts in fit range = " << p[3 * gNRes]
              << " +- " << e[3 * gNRes] << "\n";
    std::cout << "===========================================================\n\n";

    // ------------------------------------------------------------------------
    // 7) Draw -- two pads: main fit (top) + pull panel (bottom)
    // ------------------------------------------------------------------------
    gStyle->SetOptStat(0);
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

    // background component
    auto *hBgComp = static_cast<TH1D *>(gBgExc->Clone("hBgComp"));
    hBgComp->SetDirectory(nullptr);
    hBgComp->Scale(p[3 * gNRes]);
    hBgComp->SetLineColor(kGray + 2);
    hBgComp->SetLineStyle(2);
    hBgComp->SetLineWidth(2);
    hBgComp->SetFillColorAlpha(kGray, 0.35);
    hBgComp->Draw("HIST SAME");

    // individual resonances (each on top of the background, for the eye)
    const Int_t resCol[kMaxRes] = {kOrange + 1, kRed + 1, kBlue + 1,
                                   kGreen + 2, kMagenta + 1, kCyan + 2};
    std::vector<TH1D *> comps;
    for (Int_t k = 0; k < gNRes; ++k)
    {
        auto *hk = static_cast<TH1D *>(gData->Clone(Form("hRes%d", k + 1)));
        hk->SetDirectory(nullptr);
        hk->Reset();
        FoldResonance(p[3 * k], p[3 * k + 1], p[3 * k + 2], gLorb[k], hk);
        hk->Add(hBgComp);
        hk->SetLineColor(resCol[k]);
        hk->SetLineWidth(2);
        hk->SetLineStyle(7);
        hk->Draw("HIST SAME");
        comps.push_back(hk);
    }

    hModel->SetLineColor(kRed);
    hModel->SetLineWidth(3);
    hModel->Draw("HIST SAME");
    gData->Draw("E1 SAME");

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
    leg->AddEntry(hBgComp, "non-resonant bg", "lf");
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
    hPulls->Draw("E1");

    auto *lZero = new TLine(0., 0., 20., 0.);
    lZero->SetLineColor(kGray + 1);
    lZero->SetLineStyle(2);
    lZero->Draw();

    auto *lPlus2 = new TLine(0., 2., 20., 2.);
    lPlus2->SetLineColor(kGray);
    lPlus2->SetLineStyle(3);
    lPlus2->Draw();

    auto *lMinus2 = new TLine(0., -2., 20., -2.);
    lMinus2->SetLineColor(kGray);
    lMinus2->SetLineStyle(3);
    lMinus2->Draw();

    hPulls->Draw("E1 SAME");

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
    gRespN->Write("hRespNorm");
    hPulls->Write("hPulls");
    c->Write("cFitExc");
    fout->Close();
    std::cout << "[fitExc] output written to " << outFile << std::endl;
}