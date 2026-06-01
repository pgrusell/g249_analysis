// ===========================================================================
// plotTheoryOverlay.C
//
// Standalone diagnostic macro: read BOTH theoretical momentum-distribution
// file formats, normalize each to unit AREA, and overlay them on one canvas.
//
// Purpose: visually verify that the multi-column reader produces a sensible
// template shape (and how it compares to the legacy 2-column distribution),
// independently of the RooFit machinery.
//
//   format 0 : legacy two-column      ->  Qi  value
//   format 1 : JT multi-column (5 col) ->  Qi  dSdQz  dSdQt  dSdQy  dSdQ
//              (the fit uses the dS/dQ_y column = column index 3)
//
// Both files use a possibly NON-uniform Q grid (the multi-col file has a
// double-width gap across Q=0), so the histogram is built with VARIABLE bin
// edges derived from the actual Q_i values. Each point sits at its bin
// center; edges are the midpoints between neighbours.
//
// Usage:
//   root -l 'plotTheoryOverlay.C("legacy.txt", 0, "multicol.txt", 1)'
// or just edit the defaults below and run:
//   root -l plotTheoryOverlay.C
// ===========================================================================

#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <iostream>

#include "TH1F.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TStyle.h"

// ---------------------------------------------------------------------------
// Read one column ('valCol', 0-based) of momentum value vs Q_i (column 0)
// from a whitespace-separated text file. Lines that do not parse as
// 'nColsExpected' doubles are skipped (handles headers / labels / units /
// blank lines). Returns a histogram on VARIABLE bin edges built from the
// actual Q_i grid (robust to non-uniform spacing and gaps).
// ---------------------------------------------------------------------------
TH1F *readDistToHist(const std::string &txtFile,
                     int nColsExpected, // 2 for legacy, 5 for multi-col
                     int valCol,        // which column holds the value
                     const std::string &histName)
{
    std::vector<double> Qi, val;
    std::ifstream f(txtFile.c_str());
    if (!f.is_open())
    {
        std::cerr << "[readDistToHist] Cannot open " << txtFile << "\n";
        return nullptr;
    }

    std::string line;
    while (std::getline(f, line))
    {
        std::istringstream ss(line);
        std::vector<double> cols;
        double v;
        while (ss >> v)
            cols.push_back(v);
        // Keep only lines that have exactly the expected number of numeric
        // columns -> skips header / label / units / blank lines.
        if ((int)cols.size() != nColsExpected)
            continue;
        Qi.push_back(cols[0]);
        val.push_back(cols[valCol]);
    }

    const int nPts = (int)Qi.size();
    if (nPts < 2)
    {
        std::cerr << "[readDistToHist] Not enough points in " << txtFile << "\n";
        return nullptr;
    }

    // Variable bin edges from the actual grid (handles the gap at Q=0).
    std::vector<double> edges(nPts + 1);
    edges[0] = Qi[0] - 0.5 * (Qi[1] - Qi[0]);
    for (int i = 1; i < nPts; i++)
        edges[i] = 0.5 * (Qi[i - 1] + Qi[i]);
    edges[nPts] = Qi[nPts - 1] + 0.5 * (Qi[nPts - 1] - Qi[nPts - 2]);

    auto *h = new TH1F(histName.c_str(), histName.c_str(), nPts, edges.data());
    h->SetDirectory(nullptr);
    for (int i = 0; i < nPts; i++)
        h->SetBinContent(i + 1, val[i]);

    std::cout << "[readDistToHist] " << txtFile << "\n"
              << "    points=" << nPts
              << "  Qi=[" << Qi.front() << ", " << Qi.back() << "]"
              << "  rawIntegral(width)=" << h->Integral("width") << "\n";
    return h;
}

// ---------------------------------------------------------------------------
// Normalize a histogram to unit AREA (integral including bin widths).
// ---------------------------------------------------------------------------
void normalizeArea(TH1F *h)
{
    if (!h)
        return;
    const double area = h->Integral("width");
    if (area > 0)
        h->Scale(1.0 / area);
}

// ---------------------------------------------------------------------------
// Main.
// ---------------------------------------------------------------------------
void plotTheoryOverlay(
    std::string fileLegacy = "/nucl_lustre/pablogrusell/g249/g249_analysis/theory/JT/25F/sigt_p12-3.txt",
    int fmtLegacy = 0,
    std::string fileMultiCol = "/nucl_lustre/pablogrusell/g249/g249_analysis/theory/CB/momdis_1p12_19_785.txt",
    int fmtMultiCol = 1,
    double xView = 600.0) // x-axis half-range for the plot [MeV/c]
{
    gStyle->SetOptStat(0);

    // --- read both, choosing columns by format ---------------------------
    //   format 0 (legacy)    : 2 columns, value in column 1
    //   format 1 (multi-col) : 5 columns, dS/dQ_y in column 3
    TH1F *hA = (fmtLegacy == 1)
                   ? readDistToHist(fileLegacy, 5, 3, "hLegacy")
                   : readDistToHist(fileLegacy, 2, 1, "hLegacy");

    TH1F *hB = (fmtMultiCol == 1)
                   ? readDistToHist(fileMultiCol, 5, 3, "hMultiCol")
                   : readDistToHist(fileMultiCol, 2, 1, "hMultiCol");

    if (!hA || !hB)
    {
        std::cerr << "[plotTheoryOverlay] One of the files failed to load. Abort.\n";
        return;
    }

    // --- normalize each to unit area -------------------------------------
    normalizeArea(hA);
    normalizeArea(hB);

    // --- style -----------------------------------------------------------
    hA->SetLineColor(kBlue + 1);
    hA->SetLineWidth(2);
    hB->SetLineColor(kRed + 1);
    hB->SetLineWidth(2);

    hA->SetTitle("Normalized theory momentum distributions;p [MeV/c];normalized dS/dQ (unit area)");

    const double ymax = 1.25 * std::max(hA->GetMaximum(), hB->GetMaximum());
    hA->GetYaxis()->SetRangeUser(0.0, ymax);
    hA->GetXaxis()->SetRangeUser(-xView, xView);

    // --- draw ------------------------------------------------------------
    auto *c = new TCanvas("c_theory_overlay",
                          "Normalized theory distributions", 900, 650);
    hA->Draw("HIST");
    hB->Draw("HIST SAME");

    auto *leg = new TLegend(0.62, 0.74, 0.89, 0.88);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->AddEntry(hA, Form("Tostevin", fmtLegacy), "l");
    leg->AddEntry(hB, Form("Bertulani", fmtMultiCol), "l");
    leg->Draw();

    c->Update();

    // --- quick numeric comparison ---------------------------------------
    std::cout << "\n=== normalized peak heights (unit area) ===\n"
              << "  legacy     peak = " << hA->GetMaximum()
              << "  at Q = " << hA->GetBinCenter(hA->GetMaximumBin()) << " MeV/c\n"
              << "  multi-col  peak = " << hB->GetMaximum()
              << "  at Q = " << hB->GetBinCenter(hB->GetMaximumBin()) << " MeV/c\n";

    // FWHM-ish via half-max crossing on each side
    auto halfWidth = [](TH1F *h) -> double
    {
        const double half = h->GetMaximum() / 2.0;
        int lo = -1, hi = -1;
        for (int i = 1; i <= h->GetNbinsX(); i++)
            if (h->GetBinContent(i) >= half)
            {
                lo = i;
                break;
            }
        for (int i = h->GetNbinsX(); i >= 1; i--)
            if (h->GetBinContent(i) >= half)
            {
                hi = i;
                break;
            }
        if (lo < 0 || hi < 0)
            return -1;
        return h->GetBinCenter(hi) - h->GetBinCenter(lo);
    };
    std::cout << "  legacy     FWHM ~ " << halfWidth(hA) << " MeV/c\n"
              << "  multi-col  FWHM ~ " << halfWidth(hB) << " MeV/c\n";
}