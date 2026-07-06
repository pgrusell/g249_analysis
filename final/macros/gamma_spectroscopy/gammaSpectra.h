#pragma once

#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TF2.h>
#include <TVector3.h>
#include <TMath.h>
#include <TString.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TParameter.h>
#include <iostream>
#include <vector>
#include <cmath>
#include <string>
#include <fstream>
#include <sstream>

// Offsets loaded from a settings .txt file (same format as 23O1n.txt / 24O.txt).
// Line order:  fx_frag  fy_frag  dx_neu  dy_neu  beta_match
// For Doppler correction we only use fragOffsetX/Y and betaMatch.
struct OffsetParams
{
    double fragOffsetX = 0.0; // line 0: fx offset (added to -px/pz direction)
    double fragOffsetY = 0.0; // line 1: fy offset (added to -py/pz direction)
    double betaMatch   = 0.0; // line 4: additive correction to beta_frag
};

// Read offsets from $repopath/final/settings/<fileName>.
// Returns a zeroed struct (no correction) if the file cannot be opened.
inline OffsetParams readOffsets(const std::string& fileName)
{
    OffsetParams p;
    const char* repopath = getenv("repopath");
    if (!repopath) {
        std::cerr << "[readOffsets] $repopath not set — applying no corrections.\n";
        return p;
    }
    std::string path = std::string(repopath) + "/final/settings/" + fileName;
    std::ifstream in(path);
    if (!in.is_open()) {
        std::cerr << "[readOffsets] Cannot open " << path << " — applying no corrections.\n";
        return p;
    }
    std::vector<double> vals;
    std::string line;
    while (std::getline(in, line)) {
        if (line.empty() || line[0] == '#') continue;
        vals.push_back(std::atof(line.c_str()));
    }
    if (vals.size() >= 1) p.fragOffsetX = vals[0];
    if (vals.size() >= 2) p.fragOffsetY = vals[1];
    if (vals.size() >= 5) p.betaMatch   = vals[4];
    std::cout << "[readOffsets] Loaded from " << path << "\n"
              << "  fragOffsetX = " << p.fragOffsetX << "\n"
              << "  fragOffsetY = " << p.fragOffsetY << "\n"
              << "  betaMatch   = " << p.betaMatch   << "\n";
    return p;
}

struct Fit2DParams
{
    double amplitude = 0.0;
    double muAoQ     = 0.0;
    double sigmaAoQ  = 0.0;
    double muZ       = 0.0;
    double sigmaZ    = 0.0;
    bool   valid     = false;
};

inline bool insideEllipse(double aoq, double z,
                           double muAoQ, double sigmaAoQ,
                           double muZ,   double sigmaZ,
                           double k)
{
    if (sigmaAoQ <= 0.0 || sigmaZ <= 0.0) return false;
    double da = (aoq - muAoQ) / sigmaAoQ;
    double dz = (z   - muZ)   / sigmaZ;
    return (da * da + dz * dz) <= k * k;
}

// Doppler-correct a gamma energy (in any units) given the fragment direction
// (px,py,pz need not be normalised) and beta.
inline double DopplerCorrect(double E,
                              double theta_g, double phi_g,
                              double px_frag, double py_frag, double pz_frag,
                              double beta)
{
    if (beta <= 0.0 || beta >= 1.0) return E;
    const double gamma_fac = 1.0 / std::sqrt(1.0 - beta * beta);
    TVector3 gdir;
    gdir.SetMagThetaPhi(1.0, theta_g, phi_g);
    TVector3 fdir(px_frag, py_frag, pz_frag);
    if (fdir.Mag() <= 0.0) fdir.SetXYZ(0.0, 0.0, 1.0);
    fdir = fdir.Unit();
    const double cos_theta = gdir.Dot(fdir);
    return E * gamma_fac * (1.0 - beta * cos_theta);
}

// ─────────────────────────────────────────────────────────────────────────────
//  EXACT crystal-position Doppler correction using R3BCalifaGeometry.
//
//  R3BCalifaGeometry::GetAngles(crystalID) returns a TVector3 that IS the
//  Cartesian (x,y,z) of the crystal CENTRE in the global lab frame (cm) — NOT
//  angles, despite the name (it transforms the crystal local origin to master
//  coordinates via the loaded TGeo geometry). See R3BCalifaGeometry.cxx.
//
//  Therefore, with the FOOT-reconstructed vertex as the emission point, the
//  gamma direction is EXACTLY (crystal_centre - vertex) — no CALIFA-radius
//  assumption, and CALIFA's non-cylindrical forward geometry is handled
//  correctly because each crystal's true position is used.
//
//  Call InitCalifaGeometry(2025) ONCE before using DopplerCorrectCrystal.
// ─────────────────────────────────────────────────────────────────────────────
#include "R3BCalifaGeometry.h"

inline R3BCalifaGeometry* InitCalifaGeometry(int version = 2025)
{
    auto* geo = R3BCalifaGeometry::Instance();
    if (geo) geo->Init(version);   // 2025 = g249
    return geo;
}

// Exact Doppler correction: gamma direction from vertex to the true crystal
// centre (looked up by crystal ID). Falls back to the supplied (theta,phi)
// origin-direction if the geometry/ID is unavailable (geo==nullptr or id<0).
inline double DopplerCorrectCrystal(double E,
                                    int    crystalID,
                                    double theta_g, double phi_g, // fallback only
                                    double px_frag, double py_frag, double pz_frag,
                                    double beta,
                                    double vx, double vy, double vz,
                                    R3BCalifaGeometry* geo)
{
    if (beta <= 0.0 || beta >= 1.0) return E;
    const double gamma_fac = 1.0 / std::sqrt(1.0 - beta * beta);

    TVector3 crystal;
    bool haveCrys = false;
    if (geo && crystalID >= 0) {
        crystal = geo->GetAngles(crystalID);   // <-- Cartesian crystal centre [cm]
        if (crystal.Mag() > 0.0) haveCrys = true;
    }

    TVector3 gdir;
    if (haveCrys) {
        gdir = crystal - TVector3(vx, vy, vz); // exact, no radius assumption
    } else {
        // fallback: origin-based direction from (theta,phi) minus vertex offset
        TVector3 c; c.SetMagThetaPhi(30.0, theta_g, phi_g);
        gdir = c - TVector3(vx, vy, vz);
    }
    if (gdir.Mag() <= 0.0) return E;
    gdir = gdir.Unit();

    TVector3 fdir(px_frag, py_frag, pz_frag);
    if (fdir.Mag() <= 0.0) fdir.SetXYZ(0.0, 0.0, 1.0);
    fdir = fdir.Unit();

    const double cos_theta = gdir.Dot(fdir);
    return E * gamma_fac * (1.0 - beta * cos_theta);
}

// Fit a 2D Gaussian (no correlation term) to a TH2F within [aoqMin,aoqMax] x [zMin,zMax].
// Returns Fit2DParams with valid=true on success.
inline Fit2DParams Fit2DGaussian(TH2F* h,
                                  double aoqMin, double aoqMax,
                                  double zMin,   double zMax,
                                  int uid = 0)
{
    Fit2DParams p;
    if (!h || h->GetEntries() < 10) return p;

    TString fname = Form("fGammaPID_%d", uid);
    auto* f2 = new TF2(fname,
                        "[0]*exp(-0.5*((x-[1])/[2])^2 - 0.5*((y-[3])/[4])^2)",
                        aoqMin, aoqMax, zMin, zMax);

    f2->SetParameters(h->GetMaximum(),
                      0.5 * (aoqMin + aoqMax), 0.01,
                      0.5 * (zMin + zMax),     0.15);
    f2->SetParLimits(0, 0.0, 1.0e9);
    f2->SetParLimits(1, aoqMin, aoqMax);
    f2->SetParLimits(2, 0.001, 0.05);
    f2->SetParLimits(3, zMin, zMax);
    f2->SetParLimits(4, 0.02, 0.5);

    // Q=quiet, I=use integral of bin, R=respect range, 0=don't draw
    int status = h->Fit(f2, "QIRO");

    p.amplitude = f2->GetParameter(0);
    p.muAoQ     = f2->GetParameter(1);
    p.sigmaAoQ  = f2->GetParameter(2);
    p.muZ       = f2->GetParameter(3);
    p.sigmaZ    = f2->GetParameter(4);

    // Mark as valid only when fit converged and parameters are physically sensible
    p.valid = (status == 0 || status == 4000) // Minuit convergence codes
           && p.sigmaAoQ > 0.001
           && p.sigmaZ   > 0.01
           && p.muAoQ > aoqMin && p.muAoQ < aoqMax
           && p.muZ   > zMin   && p.muZ   < zMax;

    delete f2;
    return p;
}

// ═════════════════════════════════════════════════════════════════════════════
//  GammaCfg — SINGLE SOURCE OF TRUTH shared by
//     gammaSpectra.C           (data spectra)
//     califaGammaResponseSim.C (simulated response templates)
//     subtract_background.C    (simultaneous population fit)
//
//  The multiplicity-resolved population fit is a bin-by-bin Poisson likelihood,
//  so DATA and TEMPLATE histograms must have IDENTICAL binning, IDENTICAL
//  cluster acceptance and IDENTICAL multiplicity definition. Everything that
//  defines those lives here and ONLY here.
// ═════════════════════════════════════════════════════════════════════════════
namespace GammaCfg
{
    // ---- Doppler-corrected energy histogram binning (50 keV / bin) ----------
    constexpr int    NBINS_E   = 200;
    constexpr double E_MIN_MEV = 0.0;
    constexpr double E_MAX_MEV = 10.0;

    // ---- Gamma theta histogram range (degrees), 1 deg/bin --------------------
    constexpr int    NBINS_TH   = 100;
    constexpr double TH_MIN_DEG = 0.0;
    constexpr double TH_MAX_DEG = 100.0;

    // ---- Analysis-level cluster acceptance on the LAB energy -----------------
    // MUST mirror the eventFilter save window (CALIFA_GAMMA_EMIN / EMAX):
    // multiplicity M is DEFINED as the number of clusters passing this cut.
    constexpr double ECL_LAB_MIN_MEV = 0.10;   // 100 keV
    constexpr double ECL_LAB_MAX_MEV = 20.0;   // 20 MeV

    inline bool PassCluster(double eLabMeV)
    {
        return (eLabMeV > ECL_LAB_MIN_MEV) && (eLabMeV < ECL_LAB_MAX_MEV);
    }

    // ---- Multiplicity slice bookkeeping ---------------------------------------
    //  exact slices : M == m,  m = 1 .. MULT_MAX
    //  mleq  slices : M <= m,  m = 1 .. MULT_MAX
    //  mgeq  slices : M >= m,  m = 2 .. MULT_MAX   (no upper cap on M)
    // The simultaneous fit uses the DISJOINT partition { mult1, mult2, mgeq3 }.
    constexpr int MULT_MAX = 5;

    // ---- Gamma-gamma coincidence gates on the DOPPLER-CORRECTED energy -------
    // +-2 sigma around the measured 22O line positions
    // (sigma(1383) ~ 0.196 MeV, sigma(3199) ~ 0.452 MeV from the response fits)
    constexpr double GATE_1383_LO = 1.00, GATE_1383_HI = 1.77;
    constexpr double GATE_3199_LO = 2.30, GATE_3199_HI = 4.10;

    inline bool In1383Gate(double eCorrMeV)
    {
        return (eCorrMeV > GATE_1383_LO) && (eCorrMeV < GATE_1383_HI);
    }
    inline bool In3199Gate(double eCorrMeV)
    {
        return (eCorrMeV > GATE_3199_LO) && (eCorrMeV < GATE_3199_HI);
    }

    // ---- Canonical slice tags --------------------------------------------------
    inline TString SliceAll()        { return TString("all"); }
    inline TString SliceExact(int m) { return TString(Form("mult%d", m)); }
    inline TString SliceLeq(int m)   { return TString(Form("mleq%d", m)); }
    inline TString SliceGeq(int m)   { return TString(Form("mgeq%d", m)); }

    // ---- Canonical histogram / parameter names ---------------------------------
    // DATA (produced by gammaSpectra.C, consumed by subtract_background.C):
    //   Doppler-corrected energy spectrum of slice <s>:  h1_gamma_E_corr_<s>
    inline TString DataHist(const TString& slice) { return "h1_gamma_E_corr_" + slice; }

    // RESPONSE (produced by califaGammaResponseSim.C):
    //   hResp_<state>_<slice>  +  TParameter<double>  nDecays_<state>
    inline TString RespHist(const TString& state, const TString& slice)
    {
        return "hResp_" + state + "_" + slice;
    }
    inline TString NDecaysName(const TString& state) { return "nDecays_" + state; }
    inline TString RespMultHist(const TString& state) { return "hRespMult_" + state; }

    // State tags for the two fitted 22O components (higher-lying states are
    // deliberately NOT fitted — their feeding thresholds sit at/above 11.8 MeV
    // in 24O*, right below S3n, so they are expected to be negligible).
    inline TString StateDirect()  { return TString("direct_3199");    } // 2+ -> 0+
    inline TString StateCascade() { return TString("casc_1383_3199"); } // 3+ -> 2+ -> 0+

    // ─────────────────────────────────────────────────────────────────────────
    //  ARMEL-METHOD CONFIGURATION (single-line response functions)
    //
    //  Following Kamenyero's thesis Annexe 1, the population extraction uses:
    //    * SINGLE-gamma response templates (one per transition, NOT cascades)
    //    * photopeak areas fitted at each multiplicity (M==1 primary, M==2 check)
    //    * the per-line photopeak EFFICIENCY eps(E) from the simulation
    //  and then his direct/cascade equations (implemented in subtract_background.C).
    //
    //  Each transition is identified by an integer keV label; the response
    //  template and its efficiency are stored under names derived from it.
    // ─────────────────────────────────────────────────────────────────────────

    // All single transitions in 22O to build response functions for.
    // (energy in keV).  1383 and 3199 are the two used by the direct/cascade
    // decomposition; the rest are available for fitting higher-lying feeders
    // if/when they are seen.
    inline std::vector<int> AllTransitionsKeV()
    {
        return { 3199,   // 2+(3199)    -> 0+     (direct ground-state transition)
                 1383,   // 3+(4582)    -> 2+     (cascade upper rung; unique tag)
                 1710,   // (4909)      -> 2+
                 2600,   // (5800)      -> 2+
                 3310,   // (6509)      -> 2+
                 3710,   // (6936)      -> 2+
                 1218,   // (5800)      -> 3+
                 2354 }; // (6936)      -> 3+
    }

    // The two transitions used for the direct-vs-cascade decomposition.
    inline int TransDirectKeV()   { return 3199; } // shared: direct 2+ AND cascade lower rung
    inline int TransCascadeKeV()  { return 1383; } // unique cascade tag (3+ -> 2+)

    // Response-template name for a single transition, e.g. hResp_line_3199
    inline TString RespLineHist(int keV) { return TString(Form("hResp_line_%d", keV)); }

    // 2D Doppler-corrected-energy vs CALIFA theta for a transition (angular QA:
    // a correct correction gives a FLAT band at the line energy vs theta).
    inline TString RespEvsThetaHist(int keV) { return TString(Form("hResp_EvsTheta_%d", keV)); }

    // Photopeak efficiency of a transition, stored as a TParameter<double>,
    // eps = (photopeak counts in +-2sigma) / (thrown decays).  Name: eps_3199
    inline TString EffName(int keV) { return TString(Form("eps_%d", keV)); }

    // Efficiency-vs-energy curve (percent) written by the response builder.
    inline TString EffCurveGraph() { return TString("gEfficiencyPercent"); }
    inline TString EffCurveHist()  { return TString("hEfficiencyPercent"); }

    // Fit window (+-, in MeV) around a line's nominal energy when extracting its
    // photopeak area from data and when integrating the template efficiency.
    // Scales with energy to follow the ~const fractional resolution.
    inline double LineHalfWindowMeV(int keV)
    {
        const double E = keV / 1000.0;
        const double frac = 0.12;          // ~+-2.5 sigma at 9.4% resolution
        const double floorMeV = 0.20;
        return std::max(frac * E, floorMeV);
    }

    // Book-keeping parameters written by gammaSpectra.C:
    //   nFragEvents      : events after PID+OPA(+NeuLAND) cuts, NO gamma condition
    //                      (the denominator for populations)
    //   nCoinc_1383_3199 : photopeak-photopeak coincidence pairs (cascade tag)
    inline TString NFragName()  { return TString("nFragEvents"); }
    inline TString NCoincName() { return TString("nCoinc_1383_3199"); }
}
