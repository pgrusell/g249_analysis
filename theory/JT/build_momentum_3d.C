// build_momentum_3d.C
//
// Reads two text files containing 1D momentum distributions:
//   - sigt file: f_T(p_x) (transverse Cartesian component; p_y follows the same)
//   - sigl file: f_L(p_z) (longitudinal component)
// and builds a 3D histogram h3(p_x, p_y, p_z) under the assumption that the
// three Cartesian components are independent:
//
//      f(p_x, p_y, p_z) = f_T(p_x) * f_T(p_y) * f_L(p_z)
//
// This preserves the expected symmetries (azimuthal symmetry around z, and
// the symmetry p_x -> -p_x, etc.) but introduces no extra correlations.
//
// Usage in ROOT:
//   root -l -b -q 'build_momentum_3d.C("sigt.txt","sigl.txt","momentum3d.root")'
//
// Each input file is two columns: momentum value, density value.

#include <TFile.h>
#include <TH1D.h>
#include <TH3D.h>
#include <TString.h>
#include <fstream>
#include <iostream>
#include <vector>

// --- Tunables -------------------------------------------------------------
static const int NBINS = 200;      // bins per axis in the 3D histogram
static const double PMIN = -500.0; // axis lower edge (MeV/c, presumably)
static const double PMAX = 500.0;  // axis upper edge
// --------------------------------------------------------------------------

// Load a two-column text file (p, density) into parallel vectors.
bool LoadProfile(const TString &path,
                 std::vector<double> &p,
                 std::vector<double> &f)
{
    std::ifstream in(path.Data());
    if (!in)
    {
        std::cerr << "ERROR: cannot open " << path << std::endl;
        return false;
    }
    double a, b;
    while (in >> a >> b)
    {
        p.push_back(a);
        f.push_back(b);
    }
    std::cout << "Loaded " << p.size() << " points from " << path << std::endl;
    return !p.empty();
}

// Linear interpolation of the profile (p_i, f_i) at an arbitrary x.
// Assumes p_i is sorted ascending and uniformly spaced.
double Interp(const std::vector<double> &p,
              const std::vector<double> &f,
              double x)
{
    const int N = (int)p.size();
    if (N < 2)
        return 0.0;
    if (x <= p.front() || x >= p.back())
        return 0.0;

    const double dp = p[1] - p[0];
    const double t = (x - p[0]) / dp;
    int i = (int)t;
    if (i < 0)
        i = 0;
    if (i > N - 2)
        i = N - 2;
    const double frac = t - i;
    return f[i] * (1.0 - frac) + f[i + 1] * frac;
}

// Build a TH1D with the same binning as the 3D axes from a profile vector.
// Each bin gets the interpolated density at its center.
TH1D *MakeProfileHist(const char *name, const char *title,
                      const std::vector<double> &p,
                      const std::vector<double> &f)
{
    TH1D *h = new TH1D(name, title, NBINS, PMIN, PMAX);
    for (int i = 1; i <= NBINS; ++i)
    {
        double xc = h->GetBinCenter(i);
        h->SetBinContent(i, Interp(p, f, xc));
    }
    return h;
}

void build_momentum_3d(const char *sigtFile = "/nucl_lustre/pablogrusell/g249/g249_analysis/theory/JT/25F/sigt_d52-gs.txt",
                       const char *siglFile = "/nucl_lustre/pablogrusell/g249/g249_analysis/theory/JT/25F/sigl_d52-gs.txt",
                       const char *outFile = "momentum3d_d52.root")
{
    // 1. Read the two input files.
    std::vector<double> pT, fT, pL, fL;
    if (!LoadProfile(sigtFile, pT, fT))
        return;
    if (!LoadProfile(siglFile, pL, fL))
        return;

    // 2. Build 1D histograms on the target binning (for storage / cross-check).
    TH1D *hT = MakeProfileHist("hT", "Transverse component;p_{x} (or p_{y});f_{T}",
                               pT, fT);
    TH1D *hL = MakeProfileHist("hL", "Longitudinal component;p_{z};f_{L}",
                               pL, fL);

    // 3. Build the 3D histogram as the outer product f_T(px) * f_T(py) * f_L(pz).
    TH3D *h3 = new TH3D("h3_pxpypz",
                        "3D momentum;p_{x};p_{y};p_{z}",
                        NBINS, PMIN, PMAX,
                        NBINS, PMIN, PMAX,
                        NBINS, PMIN, PMAX);

    // Cache 1D bin contents to avoid repeated lookups inside the triple loop.
    std::vector<double> fx(NBINS + 2, 0.0), fy(NBINS + 2, 0.0), fz(NBINS + 2, 0.0);
    for (int i = 1; i <= NBINS; ++i)
    {
        fx[i] = hT->GetBinContent(i);
        fy[i] = hT->GetBinContent(i); // p_y uses the same transverse profile
        fz[i] = hL->GetBinContent(i);
    }

    double sum = 0.0;
    for (int iz = 1; iz <= NBINS; ++iz)
    {
        const double vz = fz[iz];
        if (vz == 0.0)
            continue;
        for (int iy = 1; iy <= NBINS; ++iy)
        {
            const double vyz = fy[iy] * vz;
            if (vyz == 0.0)
                continue;
            for (int ix = 1; ix <= NBINS; ++ix)
            {
                const double v = fx[ix] * vyz;
                if (v == 0.0)
                    continue;
                h3->SetBinContent(ix, iy, iz, v);
                sum += v;
            }
        }
    }
    std::cout << "Filled 3D histogram; raw integral (sum of bin contents) = "
              << sum << std::endl;

    // 4. Optional: normalise the 3D histogram so it integrates to 1
    //    (Integral() multiplied by the bin volume).
    const double binVol = h3->GetXaxis()->GetBinWidth(1) * h3->GetYaxis()->GetBinWidth(1) * h3->GetZaxis()->GetBinWidth(1);
    const double integral = h3->Integral() * binVol;
    std::cout << "Integral (with bin volume) before normalisation = "
              << integral << std::endl;
    if (integral > 0.0)
    {
        h3->Scale(1.0 / integral);
        std::cout << "Scaled to unit integral." << std::endl;
    }

    // 5. Write everything to a ROOT file.
    TFile fout(outFile, "RECREATE");
    hT->Write();
    hL->Write();
    h3->Write();
    fout.Close();
    std::cout << "Wrote " << outFile << std::endl;
}
