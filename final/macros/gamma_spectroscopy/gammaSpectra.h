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
