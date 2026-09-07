// ============================================================================
//  CONFIGURATION
// ============================================================================

static const Int_t kMaxRes = 6; // storage size of the per-resonance arrays

// ---- Resonances ------------------------------------------------------------
static Int_t gNRes = 3; // number of Breit-Wigner resonances

// Angular momentum of every resonance
static Int_t gLorb[kMaxRes] = {2, 1, 1, 0, 0, 0};

// Starting energies and gammas for the resonances
static Double_t gInitE[kMaxRes] = {0.60, 3.20, 6.20, 0., 0., 0.};
static Double_t gInitG[kMaxRes] = {0.05, 1.20, 1.40, 0., 0., 0.};

// Limits in the values of the resonances energies
static Double_t gElo[kMaxRes] = {0.30, 2.50, 5.00, 0., 0., 0.};
static Double_t gEhi[kMaxRes] = {0.90, 3.90, 7.50, 0., 0., 0.};

// Limits in the gammas
static Double_t gGlim[2] = {0.01, 5.0};

// Freeze different parameters
static Bool_t gFixE[kMaxRes] = {kFALSE, kFALSE, kFALSE, kFALSE, kFALSE, kFALSE};
static Bool_t gFixGamma[kMaxRes] = {true, kFALSE, kFALSE, kFALSE, kFALSE, kFALSE};

// ---- Penetrability ---------------------------------------------------------
// Channel radius
static Double_t gR0Fm = 1.3;

static Bool_t gUseShift = kTRUE;

static Double_t gAFrag = 23.;                        // fragment mass number (23O)
static Double_t gMassNeutronMeV = 939.56542;         // MeV
static Double_t gMassFragMeV = 23. * 931.494 + 14.6; // MeV; only the reduced
                                                     // mass matters, so MeV
                                                     // precision is plenty

// ---- Numerics --------------------------------------------------------------
// Sub samples of the covariance matrix (not used now)
static Int_t gNSub = 1;
static Int_t gNTrueBins = 200;

// ---- Toy MC ----------------------------------------------------------------
// just for uncertainties
static Int_t gNToys = 500;
static UInt_t gToySeed = 12345;

// ---- Input files -----------------------------------------------------------
static const char *gDataFileDef =
    "/nucl_lustre/pablogrusell/g249/g249_analysis/results/dataFiles/23O_analyzed.root";
static const char *gSimFileDef =
    "/nucl_lustre/pablogrusell/g249/g249_analysis/final/macros/erel/ana_results/full2_analysis.root";
static const char *gMixFileDef =
    "/nucl_lustre/pablogrusell/g249/g249_analysis/final/macros/erel/23O_mixing_iterative.root";

// ---- Data histogram binning (Eexc [MeV]) -----------------------------------
static Int_t gNDataBins = 100;
static Double_t gDataLo = 0.;
static Double_t gDataHi = 20.;

// ============================================================================
//  GLOBAL STATE shared with the likelihood function
//  (Minuit calls NLL() through a plain function pointer, so the ingredients
//   have to live somewhere it can reach.)
// ============================================================================

static Double_t gSn = 4.2;     // neutron separation energy [MeV]
static Double_t gFitLo = 4.2;  // fit range [MeV], in Eexc
static Double_t gFitHi = 14.0; //
static TH1D *gData = nullptr;  // data, Eexc binning
static TH1D *gBgExc = nullptr; // mixed-event template, unit integral in range
static TH2D *gRespN = nullptr; // normalised response (X = Erel_true, Y = Eexc)

// ============================================================================
//  LINE SHAPES
// ============================================================================

/// Dimensionless wave number rho = k*R for the n + fragment channel.
/// E is the relative energy in MeV; the result is dimensionless.
static Double_t NeutronRho(Double_t E)
{
    static const Double_t hbarc = 197.3269804; // MeV fm
    const Double_t mu = gMassNeutronMeV * gMassFragMeV /
                        (gMassNeutronMeV + gMassFragMeV);
    const Double_t R = gR0Fm * (TMath::Power(gAFrag, 1. / 3.) + 1.);
    return TMath::Sqrt(2. * mu * E) / hbarc * R;
}

/// Exact single-neutron penetrability P_l(E) (closed forms of the spherical
/// Bessel functions; no Coulomb term because the emitted particle is neutral).
/// For rho << 1 these reduce to rho^(2l+1), i.e. the familiar E^(l+1/2) law,
/// but they saturate for rho >~ 1 -- precisely the regime of the states
/// fitted here, so the exact form must be used.
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
        // Not reached for the l values used here; fall back on the
        // low-energy power law rather than returning nonsense.
        return TMath::Power(r, 2 * l + 1);
    }
}

/// Exact single-neutron shift function S_l(E) (closed forms of the spherical
/// Bessel functions), the real-part counterpart of NeutronPenetrability()'s
/// imaginary part in the R-matrix channel quantity (Lane & Thomas, Rev. Mod.
/// Phys. 30, 257 (1958)). The denominators are EXACTLY those of
/// NeutronPenetrability() for the same l -- keep the two in sync if either
/// is ever changed.
static Double_t NeutronShift(Double_t E, Int_t l)
{
    if (E <= 0.)
        return 0.;
    const Double_t r = NeutronRho(E);
    const Double_t r2 = r * r;
    switch (l)
    {
    case 0:
        return 0.;
    case 1:
        return -1. / (1. + r2);
    case 2:
        return -(18. + 3. * r2) / (9. + 3. * r2 + r2 * r2);
    case 3:
        return -(675. + 90. * r2 + 6. * r2 * r2) /
               (225. + 45. * r2 + 6. * r2 * r2 + r2 * r2 * r2);
    default:
        // Not reached for the l values used here; S_l -> 0 as rho -> 0 for
        // any l, so this is the sane low-energy fallback.
        return 0.;
    }
}

/// Breit-Wigner with a penetrability-driven width and its R-matrix level
/// shift, normalised so that the width equals Gamma0 at the resonance energy
/// and Delta(E0) = 0 by construction. The overall scale is irrelevant:
/// FoldShape() normalises the discretised shape to unit area and the
/// amplitude is a free parameter. gUseShift = kFALSE drops Delta (Delta = 0),
/// i.e. the line shape used before the shift was added.
static Double_t BWShape(Double_t E, Double_t E0, Double_t G0, Int_t l)
{
    if (E <= 0. || E0 <= 0. || G0 <= 0.)
        return 0.;
    const Double_t p0 = NeutronPenetrability(E0, l);
    const Double_t G = (p0 > 0.) ? G0 * NeutronPenetrability(E, l) / p0 : G0;
    const Double_t Delta = (gUseShift && p0 > 0.)
                               ? G0 * (NeutronShift(E0, l) - NeutronShift(E, l)) /
                                     (2. * p0)
                               : 0.;
    const Double_t d = E - E0 - Delta;
    return G / (d * d + 0.25 * G * G);
}

// ============================================================================
//  FOLDING
// ============================================================================

/// Fold a truth-level line shape through the response matrix and ADD the
/// result onto 'out' (which must have the data binning).
///
/// Three steps:
///   1. discretise the shape on the truth axis, integrating each bin with
///      gNSub sub-samples, and normalise the weights to sum to one;
///   2. multiply by the response matrix, rec_j = SUM_i w_i * R_ij -- each
///      truth bin has its OWN reconstructed distribution (resolution, bias,
///      tails and efficiency all measured, nothing assumed);
///   3. add A * rec_j onto bin j of 'out'. No rebinning: the Y axis of
///      gRespN was built with the data binning on purpose.
///
/// 'A' is the number of generated decays inside the truth window; what
/// survives into 'out' is A times the efficiency, because the columns of
/// gRespN are normalised per generated event.
static void FoldShape(Double_t A, TH1D *out,
                      const std::function<Double_t(Double_t)> &shape)
{
    const TAxis *axT = gRespN->GetXaxis();
    const Int_t nT = axT->GetNbins();
    const Int_t nR = gRespN->GetNbinsY();

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

    for (Int_t j = 1; j <= nR; ++j)
        if (rec[j] != 0.)
            out->AddBinContent(j, A * rec[j]);
}

/// Convenience wrapper: fold one Breit-Wigner resonance onto 'out'.
static void FoldResonance(Double_t A, Double_t E0, Double_t G0, Int_t l,
                          TH1D *out)
{
    FoldShape(A, out, [E0, G0, l](Double_t E)
              { return BWShape(E, E0, G0, l); });
}

// ============================================================================
//  MODEL AND LIKELIHOOD
// ============================================================================

/// Number of fit parameters. Layout:
///     p[3k+0] = A_k       amplitude of resonance k (efficiency-corrected)
///     p[3k+1] = Erel_k    resonance energy [MeV]
///     p[3k+2] = Gamma0_k  width at resonance [MeV]
///     p[3*gNRes] = A_bg   background counts inside the fit range
static Int_t NPar() { return 3 * gNRes + 1; }

/// Build the complete model on the data binning.
static void BuildModel(const Double_t *p, TH1D *hM)
{
    hM->Reset();
    for (Int_t k = 0; k < gNRes; ++k)
        FoldResonance(p[3 * k], p[3 * k + 1], p[3 * k + 2], gLorb[k], hM);
    hM->Add(gBgExc, p[3 * gNRes]); // mixed-event template: already
                                   // reconstructed, so never folded
}

/// Binned Poisson negative log-likelihood over the fit range, with the
/// parameter-independent ln(n!) term dropped.
static Double_t NLL(const Double_t *p)
{
    static TH1D *hM = nullptr;
    if (!hM)
    {
        hM = static_cast<TH1D *>(gData->Clone("hModel_scratch"));
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
            m = 1.e-9; // guard against log(0)
        nll += m - gData->GetBinContent(b) * TMath::Log(m);
    }
    return nll;
}

// ============================================================================
//  FITTING
// ============================================================================

/// Book the parameters and run Migrad (and optionally Hesse) on whatever
/// gData currently points to.
///
/// @param start      if non-empty, used as starting values instead of the
///                   configured defaults -- the toys restart from the nominal
///                   best fit, which is both faster and more stable.
/// @param nDataFit   data counts in range; sets the amplitude scale.
/// @param ok         set to the Migrad convergence flag.
/// @return           the minimizer; the CALLER owns it and must delete it.
static ROOT::Math::Minimizer *RunFit(const std::vector<Double_t> &start,
                                     Double_t nDataFit, Int_t printLevel,
                                     Int_t strategy, Bool_t &ok,
                                     Bool_t doHesse = kTRUE)
{
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
    min->SetErrorDef(0.5); // 0.5 = NLL (1.0 would be chi2)

    ROOT::Math::Functor fcn(&NLL, NPar());
    min->SetFunction(fcn);

    for (Int_t k = 0; k < gNRes; ++k)
    {
        const Double_t aInit = start.empty() ? 0.25 * nDataFit : start[3 * k];
        const Double_t eInit = start.empty() ? gInitE[k] : start[3 * k + 1];
        const Double_t gInit = start.empty() ? gInitG[k] : start[3 * k + 2];

        min->SetLimitedVariable(3 * k, Form("A%d", k + 1), aInit,
                                0.01 * nDataFit, 0., 20. * nDataFit);

        if (gFixE[k])
            min->SetFixedVariable(3 * k + 1, Form("Erel%d", k + 1), gInitE[k]);
        else
            min->SetLimitedVariable(3 * k + 1, Form("Erel%d", k + 1), eInit,
                                    0.02, gElo[k], gEhi[k]);

        if (gFixGamma[k])
            min->SetFixedVariable(3 * k + 2, Form("Gamma%d", k + 1), gInitG[k]);
        else
            min->SetLimitedVariable(3 * k + 2, Form("Gamma%d", k + 1), gInit,
                                    0.02, gGlim[0], gGlim[1]);
    }

    const Double_t bgInit = start.empty() ? 0.3 * nDataFit : start[3 * gNRes];
    min->SetLimitedVariable(3 * gNRes, "Abg", bgInit, 0.01 * nDataFit, 0.,
                            5. * nDataFit);

    ok = min->Minimize();
    if (doHesse)
        min->Hesse();
    return min;
}

// ============================================================================
//  INPUT LOADING
// ============================================================================

/// Read the data spectrum into gData (Eexc = Erel*1000 + Sn, in MeV).
/// The histogram is booked explicitly and filled by name, because
/// TTree::Draw's ">>h(n,lo,hi)" syntax would create a TH1F instead.
/// SetDirectory(nullptr) detaches it so it survives the file being closed.
static Bool_t LoadData(const char *dataFile, const char *cut)
{
    TFile *f = TFile::Open(dataFile, "READ");
    if (!f || f->IsZombie())
    {
        std::cerr << "ERROR: cannot open " << dataFile << "\n";
        return kFALSE;
    }
    auto *tree = dynamic_cast<TTree *>(f->Get("KinTree"));
    if (!tree)
    {
        std::cerr << "ERROR: no KinTree in " << dataFile << "\n";
        f->Close();
        return kFALSE;
    }

    gData = new TH1D("hDataExc", "", gNDataBins, gDataLo, gDataHi);
    tree->Draw(Form("Erel*1000.+%g>>hDataExc", gSn), cut, "goff");
    gData->SetDirectory(nullptr);
    gData->SetTitle(Form("^{24}O* excitation energy;E_{exc} [MeV];"
                         "counts / %.0f keV",
                         gData->GetBinWidth(1) * 1000.));
    f->Close();

    if (gData->GetEntries() == 0)
    {
        std::cerr << "ERROR: 0 entries after cut '" << cut << "'\n";
        return kFALSE;
    }
    std::cout << "[fitExc] data: " << gData->GetEntries()
              << " entries after cut\n";
    return kTRUE;
}

/// Build gRespN from the "ana" TTree of the signal simulation.
///
/// The truth range is taken from the tree itself (excluding the -999
/// sentinel) so that changing the simulation does not require touching the
/// code. 'trueHi' is returned to the caller, which uses it to cap the fit
/// range: the model is undefined -- not zero -- beyond Sn + trueHi.
///
/// Each column i is divided by the number of generated events with that
/// ErelTrue, hence SUM_j R_ij = efficiency(E_i) < 1 and the fitted
/// amplitudes are efficiency-corrected.
static Bool_t LoadResponse(const char *simFile, Double_t &trueHi)
{
    TFile *f = TFile::Open(simFile, "READ");
    if (!f || f->IsZombie())
    {
        std::cerr << "ERROR: cannot open " << simFile << "\n";
        return kFALSE;
    }
    auto *tree = dynamic_cast<TTree *>(f->Get("ana"));
    if (!tree)
    {
        std::cerr << "ERROR: no 'ana' TTree in " << simFile << "\n";
        f->Close();
        return kFALSE;
    }

    // --- truth range, auto-detected -----------------------------------------
    tree->SetEstimate(tree->GetEntries() + 1);
    const Long64_t nSel = tree->Draw("ErelTrue*1000.", "ErelTrue > -998", "goff");
    if (nSel <= 0)
    {
        std::cerr << "ERROR: no ErelTrue > -998 entries in " << simFile << "\n";
        f->Close();
        return kFALSE;
    }
    Double_t *v = tree->GetV1();
    Double_t trueLo = TMath::Floor(TMath::MinElement(nSel, v) * 10.) / 10.;
    trueHi = TMath::Ceil(TMath::MaxElement(nSel, v) * 10.) / 10.;
    if (trueLo < 0.)
        trueLo = 0.;
    std::cout << "[fitExc] response truth range from tree: [" << trueLo << ", "
              << trueHi << "] MeV\n";

    // --- migration counts and generated counts ------------------------------
    // Note TTree::Draw's "y:x" convention: Eexc goes on Y, ErelTrue on X.
    gRespN = new TH2D("hRespN", "", gNTrueBins, trueLo, trueHi, gNDataBins,
                      gDataLo, gDataHi);
    tree->Draw(Form("Erel*1000.+%g : ErelTrue*1000. >> hRespN", gSn),
               "Erel > -998 && ErelTrue > -998", "goff");
    gRespN->SetDirectory(nullptr);

    auto *hGen = new TH1D("hTrueCount", "", gNTrueBins, trueLo, trueHi);
    tree->Draw("ErelTrue*1000. >> hTrueCount", "ErelTrue > -998", "goff");
    hGen->SetDirectory(nullptr);
    f->Close();

    // --- normalise column by column -----------------------------------------
    const Int_t nT = gRespN->GetNbinsX();
    const Int_t nR = gRespN->GetNbinsY();
    Int_t nEmpty = 0;
    Double_t minGen = 1.e30;
    for (Int_t i = 1; i <= nT; ++i)
    {
        const Double_t den = hGen->GetBinContent(i);
        if (den <= 0.)
        {
            ++nEmpty;
            for (Int_t j = 1; j <= nR; ++j)
                gRespN->SetBinContent(i, j, 0.);
            continue;
        }
        minGen = TMath::Min(minGen, den);
        for (Int_t j = 1; j <= nR; ++j)
            gRespN->SetBinContent(i, j, gRespN->GetBinContent(i, j) / den);
    }
    delete hGen;

    std::cout << "[fitExc] response: " << nT << " truth bins, " << nEmpty
              << " empty, >= " << (Int_t)minGen
              << " generated events per filled column\n";
    if (minGen < 50.)
        std::cout << "[fitExc] WARNING: some truth columns are poorly "
                     "populated; consider fewer gNTrueBins or more statistics\n";
    return kTRUE;
}

/// Build gBgExc from the "tMix" TTree of the event-mixing file.
///
/// Beware of the unit conventions, which differ from the "ana" tree:
/// tMix::Erel is ALREADY in MeV and has no -999 sentinel, and every entry
/// carries a "weight" (the iterative-mixing weight times the normalisation)
/// which MUST be used -- the numeric selection argument of TTree::Draw acts
/// as a per-entry weight.
///
/// The template is normalised to unit integral inside the fit range, so A_bg
/// reads directly as "background counts in range". No smoothing is applied:
/// the ~10^7 weighted pairs make the statistical fluctuations negligible, and
/// TH1::Smooth would wash out the hard kinematic edge at Erel = 0.
static Bool_t LoadMixedBackground(const char *mixFile)
{
    TFile *f = TFile::Open(mixFile, "READ");
    if (!f || f->IsZombie())
    {
        std::cerr << "ERROR: cannot open " << mixFile << "\n";
        return kFALSE;
    }
    auto *tree = dynamic_cast<TTree *>(f->Get("tMix"));
    if (!tree)
    {
        std::cerr << "ERROR: no 'tMix' TTree in " << mixFile << "\n";
        f->Close();
        return kFALSE;
    }

    gBgExc = new TH1D("hBgExc", "", gNDataBins, gDataLo, gDataHi);
    gBgExc->Sumw2();
    tree->Draw(Form("Erel+%g >> hBgExc", gSn), "weight", "goff");
    gBgExc->SetDirectory(nullptr);
    const Long64_t nMix = tree->GetEntries();
    f->Close();

    const Double_t integ = gBgExc->Integral(gBgExc->FindFixBin(gFitLo + 1.e-9),
                                            gBgExc->FindFixBin(gFitHi - 1.e-9));
    if (integ <= 0.)
    {
        std::cerr << "ERROR: mixed-event template is empty in the fit range\n";
        return kFALSE;
    }
    gBgExc->Scale(1. / integ);

    std::cout << "[fitExc] background: event-mixing template (data-driven, "
                 "not folded), "
              << nMix << " weighted pairs\n";
    return kTRUE;
}

// ============================================================================
//  MAIN ENTRY POINT
// ============================================================================

/// @param dataFile  ROOT file with the "KinTree" of the analysed data
/// @param simFile   ROOT file with the "ana" tree of the signal simulation
/// @param mixFile   ROOT file with the "tMix" tree of the event mixing
/// @param cut       selection applied to the data tree
/// @param Sn        neutron separation energy [MeV]; Eexc = Erel + Sn
/// @param fitLo     lower edge of the fit range [MeV] (Eexc)
/// @param fitHi     upper edge; capped at Sn + response coverage
/// @param nToys     number of toy experiments; < 0 uses gNToys, 0 disables
/// @param outFile   output ROOT file
void fitExc(const char *dataFile = gDataFileDef,
            const char *simFile = gSimFileDef,
            const char *mixFile = gMixFileDef,
            const char *cut = "califa_opa > 1.25 && califa_opa < 1.65",
            Double_t Sn = 4.2,
            Double_t fitLo = 4.2,
            Double_t fitHi = 20,
            Int_t nToys = 0,
            const char *outFile = "fitExc_results.root")
{
    gSn = Sn;
    gFitLo = fitLo;
    gFitHi = fitHi;

    // ------------------------------------------------------------------
    // 1) Inputs
    // ------------------------------------------------------------------
    if (!LoadData(dataFile, cut))
        return;

    Double_t trueHi = 0.;
    if (!LoadResponse(simFile, trueHi))
        return;

    // The response matrix defines how far the model exists at all. Cap the
    // fit range BEFORE the background template is normalised, so that A_bg
    // refers to the same range the likelihood uses.
    if (gFitHi > gSn + trueHi)
    {
        std::cout << "[fitExc] WARNING: fitHi = " << gFitHi
                  << " MeV exceeds the response coverage (Sn + trueHi = "
                  << gSn + trueHi << " MeV). Clamping.\n";
        gFitHi = gSn + trueHi - gData->GetBinWidth(1);
    }

    if (!LoadMixedBackground(mixFile))
        return;

    const Int_t b1 = gData->FindFixBin(gFitLo + 1.e-9);
    const Int_t b2 = gData->FindFixBin(gFitHi - 1.e-9);
    const Double_t nDataFit = gData->Integral(b1, b2);
    std::cout << "[fitExc] fit range [" << gFitLo << ", " << gFitHi
              << "] MeV -> " << nDataFit << " counts in " << (b2 - b1 + 1)
              << " bins\n";

    // ------------------------------------------------------------------
    // 2) Nominal fit
    // ------------------------------------------------------------------
    std::cout << "\n[fitExc] Fitting " << gNRes
              << " BW resonance(s) + mixed-event background ...\n";
    Bool_t ok = kFALSE;
    ROOT::Math::Minimizer *min =
        RunFit({}, nDataFit, /*printLevel=*/1, /*strategy=*/2, ok);
    if (!min)
        return;

    const Int_t nPar = NPar();
    const Double_t *p = min->X();
    const Double_t *e = min->Errors();

    // Best-fit model, and the Pearson chi2 as a goodness-of-fit indicator
    // (NOT used in the minimisation).
    auto *hModel = static_cast<TH1D *>(gData->Clone("hModel"));
    hModel->SetDirectory(nullptr);
    BuildModel(p, hModel);

    Double_t chi2 = 0.;
    Int_t nBinsFit = 0;
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

    // ------------------------------------------------------------------
    // 3) Report
    // ------------------------------------------------------------------
    const Double_t Rfm = gR0Fm * (TMath::Power(gAFrag, 1. / 3.) + 1.);

    std::cout << "\n===================== fitExc results =====================\n";
    std::cout << (ok ? "  Minimization converged.\n"
                     : "  WARNING: minimization did NOT converge!\n");
    // CovMatrixStatus: 3 = accurate, 2 = forced positive definite,
    // 1 = approximate, 0 = not available. Below 3, prefer the toy errors.
    std::cout << "  CovMatrixStatus     = " << min->CovMatrixStatus()
              << (min->CovMatrixStatus() < 3
                      ? "  (< 3: quote the toy MC errors, not Hesse)\n"
                      : "  (accurate)\n");
    std::cout << "  chi2/ndf (Pearson)  = " << std::fixed << std::setprecision(2)
              << chi2 << " / " << ndf << " = " << (ndf > 0 ? chi2 / ndf : 0.)
              << "\n";
    std::cout << "  NLL at minimum      = " << std::setprecision(6)
              << min->MinValue() << "\n";
    std::cout << "  channel radius R    = " << std::setprecision(3) << Rfm
              << " fm (r0 = " << gR0Fm << " fm)\n\n";

    // Truth-axis range of the response matrix, used below for the numerical
    // peak scan (same range FoldShape() integrates the line shape over).
    const Double_t truthLo = gRespN->GetXaxis()->GetXmin();
    const Double_t truthHi = gRespN->GetXaxis()->GetXmax();

    for (Int_t k = 0; k < gNRes; ++k)
    {
        std::cout << "  Resonance " << k + 1 << "  (l = " << gLorb[k]
                  << ",  rho(E0) = " << NeutronRho(p[3 * k + 1]) << ")\n";
        std::cout << "     Erel   = " << p[3 * k + 1] << " +- ";
        if (gFixE[k])
            std::cout << "(fixed)";
        else
            std::cout << e[3 * k + 1];
        std::cout << " MeV\n";
        std::cout << "     Eexc   = " << p[3 * k + 1] + gSn << " +- ";
        if (gFixE[k])
            std::cout << "(fixed)";
        else
            std::cout << e[3 * k + 1];
        std::cout << " MeV   (Sn = " << gSn << ")\n";
        std::cout << "     Gamma0 = " << p[3 * k + 2] << " +- ";
        if (gFixGamma[k])
            std::cout << "(fixed)";
        else
            std::cout << e[3 * k + 2];
        std::cout << " MeV\n";
        std::cout << "     Yield  = " << p[3 * k] << " +- " << e[3 * k]
                  << "   (efficiency-corrected decays)\n";

        // ---- level-shift diagnostics ----------------------------------
        // Erel_k is the R-matrix POLE (Delta(E0) = 0 by construction), not
        // the observed peak -- see the gInitE note and the MODEL header.
        // Report S_l(E0) and Delta(E) at the fit-range edges to gauge the
        // size of the effect, and the numerically-found peak for comparison.
        const Double_t E0 = p[3 * k + 1];
        const Double_t G0 = p[3 * k + 2];
        const Int_t l = gLorb[k];
        const Double_t p0 = NeutronPenetrability(E0, l);
        const Double_t sl0 = NeutronShift(E0, l);
        std::cout << "     S_l(E0)= " << sl0 << "\n";

        auto deltaAt = [&](Double_t E) -> Double_t
        {
            if (p0 <= 0. || E <= 0.)
                return 0.;
            return G0 * (sl0 - NeutronShift(E, l)) / (2. * p0);
        };
        const Double_t ErelLo = gFitLo - gSn;
        const Double_t ErelHi = gFitHi - gSn;
        std::cout << "     Delta(E_fitLo=" << ErelLo << ")=" << deltaAt(ErelLo)
                  << " MeV,  Delta(E_fitHi=" << ErelHi << ")=" << deltaAt(ErelHi)
                  << " MeV   (gUseShift = " << (gUseShift ? "kTRUE" : "kFALSE")
                  << ")\n";

        Double_t erelPeak = E0;
        Double_t bwMax = -1.;
        const Int_t nScan = 2000;
        for (Int_t s = 0; s <= nScan; ++s)
        {
            const Double_t E = truthLo + (truthHi - truthLo) * s / nScan;
            const Double_t bw = BWShape(E, E0, G0, l);
            if (bw > bwMax)
            {
                bwMax = bw;
                erelPeak = E;
            }
        }
        std::cout << "     pole Erel = " << E0 << " MeV,  peak Erel = "
                  << erelPeak << " MeV   (numerical scan over truth range ["
                  << truthLo << ", " << truthHi << "] MeV)\n";
    }
    std::cout << "  Background (mixed events) = " << p[3 * gNRes] << " +- "
              << e[3 * gNRes] << " counts in range\n";
    std::cout << "===========================================================\n\n";

    // ------------------------------------------------------------------
    // 4) Toy MC (parametric bootstrap)
    //
    //    Each toy: draw Poisson pseudo-data from the best-fit model, point
    //    gData at it, refit from the nominal minimum, record the parameters.
    //    Migrad strategy 1 and no Hesse -- the toys measure the spread
    //    themselves, so the per-toy covariance is not needed.
    // ------------------------------------------------------------------
    const Int_t nToysRun = (nToys >= 0) ? nToys : gNToys;

    // Names and indices of the FREE parameters (fixed ones have no spread).
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

    if (nToysRun > 0)
    {
        TH1D *dataReal = gData; // keep the real data safe
        const std::vector<Double_t> startNom(p, p + nPar);
        std::vector<std::vector<Double_t>> toyVals(nFree);

        toyTree = new TTree("toys", "Toy MC fit results");
        std::vector<Double_t> branchVal(nPar, 0.);
        Int_t toyStatus = 0;
        for (Int_t i = 0; i < nPar; ++i)
            toyTree->Branch(parNames[i], &branchVal[i]);
        toyTree->Branch("status", &toyStatus);

        gRandom->SetSeed(gToySeed);
        Int_t nConverged = 0;

        for (Int_t t = 1; t <= nToysRun; ++t)
        {
            auto *toyData = static_cast<TH1D *>(dataReal->Clone("hToyData"));
            toyData->SetDirectory(nullptr);
            toyData->Reset();
            for (Int_t b = b1; b <= b2; ++b)
                toyData->SetBinContent(
                    b, gRandom->Poisson(hModel->GetBinContent(b)));

            gData = toyData;
            Bool_t okToy = kFALSE;
            ROOT::Math::Minimizer *minToy =
                RunFit(startNom, nDataFit, /*printLevel=*/-1, /*strategy=*/1,
                       okToy, /*doHesse=*/kFALSE);

            if (minToy)
            {
                if (okToy && minToy->Status() == 0)
                {
                    const Double_t *pToy = minToy->X();
                    for (Int_t i = 0; i < nPar; ++i)
                        branchVal[i] = pToy[i];
                    toyStatus = minToy->Status();
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

        gData = dataReal; // ALWAYS restore before anything else uses gData

        std::cout << "[fitExc] toys converged: " << nConverged << " / "
                  << nToysRun << "\n";
        if (nConverged < 0.95 * nToysRun)
            std::cout << "[fitExc] WARNING: fewer than 95% of the toys "
                         "converged -- suspect a degenerate parameter\n";

        if (nConverged > 0)
        {
            // ---- marginals, quantiles, bias --------------------------------
            // The central 68% interval runs from the 16% to the 84% quantile
            // (NOT 32-68%); quoting the median -/+ those distances captures
            // any asymmetry, which Hesse cannot.
            std::cout << "\n------------------- toy MC summary "
                         "-------------------\n";
            std::cout << std::left << std::setw(10) << "param" << std::right
                      << std::setw(12) << "nominal" << std::setw(12) << "median"
                      << std::setw(10) << "-sigma" << std::setw(10) << "+sigma"
                      << std::setw(10) << "bias" << std::setw(10) << "Hesse"
                      << "\n";
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
                          << std::right << std::setw(12) << p[idx]
                          << std::setw(12) << quant[1] << std::setw(10)
                          << (quant[1] - quant[0]) << std::setw(10)
                          << (quant[2] - quant[1]) << std::setw(10)
                          << (quant[1] - p[idx]) << std::setw(10) << e[idx]
                          << "\n";
            }

            // ---- correlation matrix ----------------------------------------
            // Pearson correlations of the toy parameter values. Strong
            // correlations are not a bug in themselves (an amplitude and its
            // width are naturally correlated), but they inflate the effective
            // uncertainty of the pair and should be quoted alongside it.
            std::vector<Double_t> mean(nFree, 0.), sig(nFree, 0.);
            for (Int_t f = 0; f < nFree; ++f)
            {
                const Int_t n = (Int_t)toyVals[f].size();
                Double_t s = 0.;
                for (Double_t x : toyVals[f])
                    s += x;
                mean[f] = s / n;
                Double_t s2 = 0.;
                for (Double_t x : toyVals[f])
                    s2 += (x - mean[f]) * (x - mean[f]);
                sig[f] = TMath::Sqrt(s2 / n);
            }

            std::cout << "\n  toy correlation matrix:\n         ";
            for (Int_t f = 0; f < nFree; ++f)
                std::cout << std::setw(8) << parNames[freeIdx[f]];
            std::cout << "\n"
                      << std::setprecision(2);
            for (Int_t fi = 0; fi < nFree; ++fi)
            {
                std::cout << std::setw(9) << parNames[freeIdx[fi]];
                for (Int_t fj = 0; fj < nFree; ++fj)
                {
                    Double_t rho = (fi == fj) ? 1. : 0.;
                    if (sig[fi] > 0. && sig[fj] > 0.)
                    {
                        Double_t cov = 0.;
                        const Int_t n = (Int_t)toyVals[fi].size();
                        for (Int_t kk = 0; kk < n; ++kk)
                            cov += (toyVals[fi][kk] - mean[fi]) *
                                   (toyVals[fj][kk] - mean[fj]);
                        rho = cov / n / (sig[fi] * sig[fj]);
                    }
                    std::cout << std::setw(8) << rho;
                }
                std::cout << "\n";
            }
            std::cout << "--------------------------------------------------"
                         "------\n\n";
        }
    }

    // ------------------------------------------------------------------
    // 5) Plot: fit on top, pulls below
    // ------------------------------------------------------------------
    gStyle->SetCanvasPreferGL();
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

    // ---- top pad -----------------------------------------------------------
    padTop->cd();
    gData->SetMarkerStyle(20);
    gData->SetMarkerSize(0.8);
    gData->SetLineColor(kBlack);
    gData->GetXaxis()->SetLabelSize(0.); // the X axis belongs to the pull pad
    gData->GetXaxis()->SetTitleSize(0.);
    gData->Draw("E1");

    // Individual components are drawn BARE (from zero), not stacked on the
    // background: the total is the red curve, so the coloured areas are not
    // additive by eye where they overlap.
    const Int_t resCol[kMaxRes] = {kOrange + 1, kMagenta + 1, kCyan + 2,
                                   kGreen + 2, kRed + 1, kBlue + 1};
    std::vector<TH1D *> comps;
    for (Int_t k = 0; k < gNRes; ++k)
    {
        auto *hk = static_cast<TH1D *>(gData->Clone(Form("hRes%d", k + 1)));
        hk->SetDirectory(nullptr);
        hk->Reset();
        FoldResonance(p[3 * k], p[3 * k + 1], p[3 * k + 2], gLorb[k], hk);
        hk->SetLineColor(resCol[k]);
        hk->SetFillColorAlpha(resCol[k], 0.3);
        hk->SetLineWidth(2);
        hk->Draw("HIST SAME");
        comps.push_back(hk);
    }

    auto *hBgComp = static_cast<TH1D *>(gBgExc->Clone("hBgComp"));
    hBgComp->SetDirectory(nullptr);
    hBgComp->Scale(p[3 * gNRes]);
    hBgComp->SetLineColor(kGray + 2);
    hBgComp->SetLineStyle(2);
    hBgComp->SetLineWidth(2);
    hBgComp->SetFillColorAlpha(kGray, 0.35);
    hBgComp->Draw("HIST SAME");

    hModel->SetLineColor(kRed);
    hModel->SetLineWidth(3);
    hModel->Draw("HIST SAME");
    gData->Draw("E1 SAME");

    auto *leg = new TLegend(0.58, 0.58, 0.93, 0.92);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->AddEntry(gData, "data", "pe");
    leg->AddEntry(hModel, "total fit", "l");
    for (Int_t k = 0; k < gNRes; ++k)
        leg->AddEntry(comps[k],
                      Form("E_{exc}=%.2f, #Gamma_{0}=%.2f MeV (l=%d)",
                           p[3 * k + 1] + gSn, p[3 * k + 2], gLorb[k]),
                      "lf");
    leg->AddEntry(hBgComp, "mixed-event bg", "lf");
    leg->Draw();

    // Fit-range markers
    for (Double_t x : {gFitLo, gFitHi})
    {
        auto *l = new TLine(x, 0., x, gData->GetMaximum() * 0.35);
        l->SetLineStyle(3);
        l->SetLineColor(kBlack);
        l->Draw();
    }

    auto *txt = new TLatex();
    txt->SetNDC();
    txt->SetTextSize(0.035);
    txt->DrawLatex(0.15, 0.93, Form("#chi^{2}/ndf = %.1f / %d", chi2, ndf));

    // ---- bottom pad: pulls -------------------------------------------------
    // (data - model)/sqrt(model), ONLY for the bins that entered the
    // likelihood. A TGraph is used rather than a histogram so that bins
    // outside the fit range produce no marker at all (a TH1 drawn with "P"
    // would show a point at zero for every empty bin). If the model is
    // correct these should scatter like a standard normal: no structure,
    // ~68% within +-1. Coherent runs of same-sign pulls indicate a model
    // deficiency; an isolated large one smells like a data artefact.
    auto *gPulls = new TGraph();
    gPulls->SetName("gPulls");
    Double_t maxAbsPull = 0.;
    for (Int_t b = b1; b <= b2; ++b)
    {
        const Double_t m = hModel->GetBinContent(b);
        if (m <= 0.)
            continue;
        const Double_t pull = (gData->GetBinContent(b) - m) / TMath::Sqrt(m);
        gPulls->SetPoint(gPulls->GetN(), gData->GetBinCenter(b), pull);
        maxAbsPull = TMath::Max(maxAbsPull, TMath::Abs(pull));
    }
    const Double_t pullRange = TMath::Max(3., maxAbsPull * 1.2);

    padBot->cd();
    // An empty frame fixes the axes, so the X range still matches the top pad
    // exactly while the graph supplies only the points inside the fit range.
    auto *hFrame = new TH1D("hPullFrame", ";E_{exc} [MeV];pull", gNDataBins,
                            gDataLo, gDataHi);
    hFrame->SetDirectory(nullptr);
    hFrame->GetYaxis()->SetRangeUser(-pullRange, pullRange);
    hFrame->GetXaxis()->SetLabelSize(0.10);
    hFrame->GetXaxis()->SetTitleSize(0.10);
    hFrame->GetXaxis()->SetTitleOffset(1.3);
    hFrame->GetYaxis()->SetLabelSize(0.10);
    hFrame->GetYaxis()->SetTitleSize(0.10);
    hFrame->GetYaxis()->SetTitleOffset(0.5);
    hFrame->GetYaxis()->SetNdivisions(505);
    hFrame->Draw("AXIS");

    for (Double_t y : {-1., 0., 1.})
    {
        auto *l = new TLine(gDataLo, y, gDataHi, y);
        l->SetLineColor(y == 0. ? kGray + 1 : kGray);
        l->SetLineStyle(y == 0. ? 2 : 3);
        l->Draw();
    }

    gPulls->SetMarkerStyle(20);
    gPulls->SetMarkerSize(0.7);
    gPulls->SetMarkerColor(kBlack);
    gPulls->Draw("P SAME");

    c->cd();
    c->Update();

    // ------------------------------------------------------------------
    // 6) Save
    // ------------------------------------------------------------------
    TFile *fout = TFile::Open(outFile, "RECREATE");
    gData->Write("hData");
    hModel->Write("hModelTotal");
    hBgComp->Write("hBgComponent");
    for (size_t k = 0; k < comps.size(); ++k)
        comps[k]->Write(Form("hResonance%zu", k + 1));
    gRespN->Write("hRespNorm");
    gBgExc->Write("hBgTemplate");
    gPulls->Write("gPulls");
    c->Write("cFitExc");
    if (toyTree)
        toyTree->Write();
    for (auto *h : toyMarginals)
        h->Write();
    fout->Close();

    delete min;
    std::cout << "[fitExc] output written to " << outFile << "\n";
}