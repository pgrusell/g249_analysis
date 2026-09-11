// ═══════════════════════════════════════════════════════════════════════════
//  neutronEventDisplay.C
//
//  Interactive event display for the NeuLAND hits stored by
//  final/macros/data_reduction/eventFilter.C  (tree "FilterDataTree").
//
//  Six projections of the neutron hit coordinates are shown per event:
//      y vs x   (front view)      t vs x
//      x vs z   (top view)        t vs y
//      y vs z   (side view)       t vs z
//  plus a 3D view and an information panel. The panel gives the relative
//  energy this neutron contributes, computed exactly as in
//  final/macros/momentum_distributions/dataAnalysis.cpp with the offsets of
//  final/settings/23O1n.txt.
//
//  A second canvas holds the Erel spectrum of the whole selection, built three
//  times: taking the first hit in z (what the tree and the display use),
//  taking the first hit in time, and keeping only the events where both
//  criteria pick the same hit.
//
//  Only quasi-free scattering events are shown by default, that is
//  1.25 < califa_opa < 1.65 rad. Use setQFS(false) to see all of them.
//
//  Usage:
//      root -l 'neutronEventDisplay.C'                  // default 23O file
//      root -l 'neutronEventDisplay.C("myfile.root")'
//
//  Then, at the ROOT prompt (or with the buttons of the control panel):
//      nextEvent()          prevEvent()        randomEvent()
//      showEvent(entry)     printEvent()       saveEvent("name.png")
//      setMinHits(2)        applyCut("beta_neu>0.6")
//      setQFS(false)        setRange("t", 55, 80)
//      loadOffsets("23O1n.txt")                 playEvents(20, 0.5)
//      nedHelp()
//
//  With the mouse over the canvas: n = next, p = previous, r = random,
//  d = dump the hit table, s = save the canvas, h = help.
// ═══════════════════════════════════════════════════════════════════════════

#include <iostream>
#include <iomanip>
#include <fstream>
#include <cstdlib>
#include <vector>
#include <string>
#include <cmath>
#include <algorithm>

#include "TROOT.h"
#include "TSystem.h"
#include "TFile.h"
#include "TTree.h"
#include "TTreeFormula.h"
#include "TStyle.h"
#include "TColor.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TH1.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TMath.h"
#include "TMarker.h"
#include "TLine.h"
#include "TBox.h"
#include "TLatex.h"
#include "TPaveText.h"
#include "TLegend.h"
#include "TPolyMarker3D.h"
#include "TRandom3.h"
#include "TString.h"

#include "TGClient.h"
#include "TGFrame.h"
#include "TGLayout.h"
#include "TGButton.h"
#include "TGLabel.h"
#include "TGNumberEntry.h"

// ─── Configuration ──────────────────────────────────────────────────────────

static const char *kNedDefaultFile =
    "/nucl_lustre/pablogrusell/g249/g249_analysis/results/dataFiles/"
    "data_22O_with_neutrons.root";
static const char *kNedTreeName = "FilterDataTree";

// quasi-free scattering selection on the CALIFA opening angle [rad]
static constexpr double kNedOpaMin = 1.25;
static constexpr double kNedOpaMax = 1.65;

static constexpr double kNedMn = 0.939565;  // neutron mass [GeV]
static constexpr double kNedC = 29.9792458; // speed of light [cm/ns]

static const char *kNedDefaultOffsets = "23O1n.txt";
static const char *kNedSettingsDir =
    "/nucl_lustre/pablogrusell/g249/g249_analysis/final/settings/";

// Variable indices used by the panel definitions: 0=x, 1=y, 2=z, 3=t
static const int kNedU[6] = {0, 2, 2, 0, 1, 2}; // horizontal axis of each panel
static const int kNedV[6] = {1, 0, 1, 3, 3, 3}; // vertical axis of each panel

static const char *kNedVarName[4] = {"x [cm]", "y [cm]", "z [cm]", "t [ns]"};

static const char *kNedPanelTitle[6] = {
    "Front view (beam into page): y vs x",
    "Top view: x vs z",
    "Side view: y vs z",
    "Time vs x",
    "Time vs y",
    "Time vs z"};

// ─── State ──────────────────────────────────────────────────────────────────

static TFile *gNedFile = nullptr;
static TTree *gNedTree = nullptr;
static TCanvas *gNedCanvas = nullptr;

// branch buffers (all NeuLAND hits of the event)
static std::vector<double> *gNedX = nullptr;
static std::vector<double> *gNedY = nullptr;
static std::vector<double> *gNedZ = nullptr;
static std::vector<double> *gNedT = nullptr;
static std::vector<int> *gNedPaddle = nullptr;

// branch buffers (first hit = the one used as "the neutron", plus kinematics)
static double gNedXf = 0, gNedYf = 0, gNedZf = 0, gNedTf = 0;
static double gNedBetaNeu = 0, gNedPNeu = 0;
static double gNedPxNeu = 0, gNedPyNeu = 0, gNedPzNeu = 0;
static double gNedAoQ = 0, gNedZfrag = 0, gNedBetaFrag = 0, gNedOpa = 0;
static double gNedMFrag = 0;
static double gNedPxFrag = 0, gNedPyFrag = 0, gNedPzFrag = 0;

// Erel calibration: transverse direction offsets and the beta matching value,
// read from final/settings/<file> (fragment x, fragment y, neutron x,
// neutron y, beta matching)
static double gNedFragOff[2] = {0., 0.};
static double gNedNeuOff[2] = {0., 0.};
static double gNedBetaMatch = 0.;
static TString gNedOffsetsFile = "";

// selection
static std::vector<Long64_t> gNedList; // tree entries passing the selection
static int gNedPos = 0;                // current position inside gNedList
static Long64_t gNedEntry = -1;        // current tree entry
static bool gNedInList = false;        // is that entry part of the selection?
static int gNedMinHits = 1;
static TString gNedCut = "";
static bool gNedQFS = true; // keep only quasi-free scattering events

// axis limits for x, y, z, t
static double gNedLo[4] = {0, 0, 0, 0};
static double gNedHi[4] = {0, 0, 0, 0};

// drawing objects that survive between events
static TH2F *gNedFrame[6] = {nullptr, nullptr, nullptr, nullptr, nullptr, nullptr};
static TH3F *gNedFrame3D = nullptr;

// Erel spectrum canvas
static TCanvas *gNedErelCanvas = nullptr;
static TH1F *gNedErelZ = nullptr;    // first hit in z (the one used by the display)
static TH1F *gNedErelT = nullptr;    // first hit in time
static TH1F *gNedErelSame = nullptr; // events where both are the same hit
static int gNedErelBins = 100;
static double gNedErelLo = -5.;
static double gNedErelHi = 10.;

// GUI
static TGMainFrame *gNedGui = nullptr;
static TGNumberEntry *gNedGotoBox = nullptr;
static TGLabel *gNedStatus = nullptr;

// forward declarations (needed by the interpreted GUI slots)
void nextEvent();
void prevEvent();
void randomEvent();
void showEvent(Long64_t entry);
void printEvent();
void saveEvent(const char *fileName = nullptr);
static void nedMakeFrames();
static void nedDrawEvent();
static void nedFillErelSpectrum();
void erelSpectrum(int nbins = 100, double lo = -5., double hi = 10.);
void nedHelp();

// ─── Small helpers ──────────────────────────────────────────────────────────

/// Drop the current file and leave the display in a well defined idle state.
static void nedUnload()
{
    if (gNedFile)
        gNedFile->Close();
    gNedFile = nullptr;
    gNedTree = nullptr;
    gNedList.clear();
    gNedEntry = -1;
    gNedInList = false;
}

static bool nedReady()
{
    if (!gNedTree)
    {
        std::cerr << "[ERROR] No tree loaded. Run neutronEventDisplay() first.\n";
        return false;
    }
    return true;
}

/// Colour of a hit according to its time (blue = early, red = late).
static Int_t nedTimeColor(double t)
{
    const int n = TColor::GetNumberOfColors();
    if (n <= 0)
        return kAzure + 1;

    double span = gNedHi[3] - gNedLo[3];
    double f = (span > 0) ? (t - gNedLo[3]) / span : 0.5;
    f = std::max(0.0, std::min(1.0, f));
    return TColor::GetColorPalette((int)(f * (n - 1)));
}

/// Value of variable `v` (0=x,1=y,2=z,3=t) for one hit.
static double nedVar(int v, double x, double y, double z, double t)
{
    switch (v)
    {
    case 0:
        return x;
    case 1:
        return y;
    case 2:
        return z;
    default:
        return t;
    }
}

static void nedPadRange(int p, double &uLo, double &uHi, double &vLo, double &vHi)
{
    uLo = gNedLo[kNedU[p]];
    uHi = gNedHi[kNedU[p]];
    vLo = gNedLo[kNedV[p]];
    vHi = gNedHi[kNedV[p]];
}

static void nedSetStatus(const TString &text)
{
    if (gNedStatus)
    {
        gNedStatus->SetText(text.Data());
        if (gNedGui)
            gNedGui->Layout();
    }
}

// ─── Tree access ────────────────────────────────────────────────────────────

static void nedSetBranches()
{
    gNedTree->SetBranchAddress("n_x", &gNedX);
    gNedTree->SetBranchAddress("n_y", &gNedY);
    gNedTree->SetBranchAddress("n_z", &gNedZ);
    gNedTree->SetBranchAddress("n_t", &gNedT);
    gNedTree->SetBranchAddress("n_paddle", &gNedPaddle);

    gNedTree->SetBranchAddress("x_neu_hit", &gNedXf);
    gNedTree->SetBranchAddress("y_neu_hit", &gNedYf);
    gNedTree->SetBranchAddress("z_neu_hit", &gNedZf);
    gNedTree->SetBranchAddress("tof_neuland", &gNedTf);

    gNedTree->SetBranchAddress("beta_neu", &gNedBetaNeu);
    gNedTree->SetBranchAddress("p_neu", &gNedPNeu);
    gNedTree->SetBranchAddress("px_neu", &gNedPxNeu);
    gNedTree->SetBranchAddress("py_neu", &gNedPyNeu);
    gNedTree->SetBranchAddress("pz_neu", &gNedPzNeu);

    gNedTree->SetBranchAddress("M_frag", &gNedMFrag);
    gNedTree->SetBranchAddress("px_frag", &gNedPxFrag);
    gNedTree->SetBranchAddress("py_frag", &gNedPyFrag);
    gNedTree->SetBranchAddress("pz_frag", &gNedPzFrag);

    gNedTree->SetBranchAddress("AoQ_frag", &gNedAoQ);
    gNedTree->SetBranchAddress("Z_frag_est", &gNedZfrag);
    gNedTree->SetBranchAddress("beta_frag", &gNedBetaFrag);
    gNedTree->SetBranchAddress("califa_opa", &gNedOpa);
}

/// Robust global x/y/z/t ranges: 0.5 % - 99.5 % quantiles of all hits plus a
/// margin. Quantiles are used because a few hits carry meaningless times
/// (tens of thousands of ns) that would otherwise squeeze every axis.
static void nedComputeRanges()
{
    std::vector<double> val[4];
    const Long64_t n = gNedTree->GetEntries();
    for (Long64_t i = 0; i < n; ++i)
    {
        gNedTree->GetEntry(i);
        for (size_t h = 0; h < gNedX->size(); ++h)
        {
            val[0].push_back((*gNedX)[h]);
            val[1].push_back((*gNedY)[h]);
            val[2].push_back((*gNedZ)[h]);
            val[3].push_back((*gNedT)[h]);
        }
    }

    for (int k = 0; k < 4; ++k)
    {
        if (val[k].empty())
        {
            gNedLo[k] = 0;
            gNedHi[k] = 1;
            continue;
        }

        std::sort(val[k].begin(), val[k].end());
        const size_t last = val[k].size() - 1;
        const double lo = val[k][(size_t)(0.005 * last)];
        const double hi = val[k][(size_t)(0.995 * last)];
        const double margin = 0.08 * (hi - lo) + 1e-6;
        gNedLo[k] = lo - margin;
        gNedHi[k] = hi + margin;
    }
}

// ─── Relative energy ────────────────────────────────────────────────────────

/// Read the offsets used by dataAnalysis.cpp: fragment x, fragment y,
/// neutron x, neutron y, beta matching (one value per line).
void loadOffsets(const char *offFile = kNedDefaultOffsets)
{
    const char *repo = getenv("repopath");
    const TString path = repo ? TString::Format("%s/final/settings/%s", repo, offFile)
                              : TString::Format("%s%s", kNedSettingsDir, offFile);

    std::ifstream in(path.Data());
    if (!in.is_open())
    {
        std::cerr << "[WARN] Cannot open " << path
                  << ", Erel will use zero offsets\n";
        return;
    }

    std::vector<double> v;
    std::string line;
    while (std::getline(in, line))
        if (!line.empty())
            v.push_back(std::atof(line.c_str()));

    if (v.size() < 5)
    {
        std::cerr << "[WARN] " << path << " has " << v.size()
                  << " values instead of 5, Erel will use zero offsets\n";
        return;
    }

    gNedFragOff[0] = v[0];
    gNedFragOff[1] = v[1];
    gNedNeuOff[0] = v[2];
    gNedNeuOff[1] = v[3];
    gNedBetaMatch = v[4];
    gNedOffsetsFile = offFile;

    std::cout << "[OK] Erel offsets from " << path << " : frag ("
              << gNedFragOff[0] << ", " << gNedFragOff[1] << ")  neu ("
              << gNedNeuOff[0] << ", " << gNedNeuOff[1] << ")  betaMatch "
              << gNedBetaMatch << "\n";

    if (gNedEntry >= 0)
        nedDrawEvent();
}

/// Relative energy [GeV] of the fragment plus a neutron of the given velocity
/// and momentum, same expression as DataAnalysis::getData().
static double nedErelFrom(double betaN, double pxNeu, double pyNeu, double pzNeu)
{
    const double betaF = gNedBetaFrag + gNedBetaMatch;

    if (betaF <= 0. || betaF >= 1. || betaN <= 0. || betaN >= 1. ||
        pzNeu == 0. || gNedPzFrag == 0. || gNedMFrag <= 0.)
        return -999.;

    const double dxNeu = pxNeu / pzNeu - gNedNeuOff[0];
    const double dyNeu = pyNeu / pzNeu - gNedNeuOff[1];
    const double dxFrag = gNedPxFrag / gNedPzFrag - gNedFragOff[0];
    const double dyFrag = gNedPyFrag / gNedPzFrag - gNedFragOff[1];

    const double cosAng =
        (dxNeu * dxFrag + dyNeu * dyFrag + 1.0) /
        (std::sqrt(dxNeu * dxNeu + dyNeu * dyNeu + 1.0) *
         std::sqrt(dxFrag * dxFrag + dyFrag * dyFrag + 1.0));

    const double gammaN = 1.0 / std::sqrt(1.0 - betaN * betaN);
    const double gammaF = 1.0 / std::sqrt(1.0 - betaF * betaF);
    const double mF = gNedMFrag;

    return std::sqrt(mF * mF + kNedMn * kNedMn +
                     2.0 * gammaN * gammaF * mF * kNedMn *
                         (1.0 - betaN * betaF * cosAng)) -
           mF - kNedMn;
}

/// Relative energy [GeV] of the event as displayed: the neutron branches of
/// the tree, that is the first hit in z.
static double nedErel()
{
    return nedErelFrom(gNedBetaNeu, gNedPxNeu, gNedPyNeu, gNedPzNeu);
}

/// Relative energy [GeV] using one NeuLAND hit, with the neutron velocity and
/// momentum rebuilt from its position and time of flight as in eventFilter.C.
static double nedErelFromHit(double x, double y, double z, double t)
{
    const double L = std::sqrt(x * x + y * y + z * z);
    if (L <= 0. || t <= 0.)
        return -999.;

    const double beta = L / (kNedC * t);
    if (beta <= 0. || beta >= 1.)
        return -999.;

    const double p = kNedMn * beta / std::sqrt(1.0 - beta * beta);
    return nedErelFrom(beta, p * x / L, p * y / L, p * z / L);
}

/// Index in the hit vectors of the hit the tree kept as the neutron (the first
/// one in z). Times carry no smearing, so it is matched by time of flight.
static int nedFirstHitIndex()
{
    int best = -1;
    double bestDz = 1e30;

    for (size_t h = 0; h < gNedT->size(); ++h)
    {
        if (std::abs((*gNedT)[h] - gNedTf) > 1e-9)
            continue;
        const double dz = std::abs((*gNedZ)[h] - gNedZf);
        if (dz < bestDz)
        {
            bestDz = dz;
            best = (int)h;
        }
    }
    return best;
}

/// Index of the earliest hit in time.
static int nedFirstTimeIndex()
{
    int best = -1;
    double tMin = 1e30;

    for (size_t h = 0; h < gNedT->size(); ++h)
        if ((*gNedT)[h] < tMin)
        {
            tMin = (*gNedT)[h];
            best = (int)h;
        }
    return best;
}

// ─── Selection ──────────────────────────────────────────────────────────────

/// Rebuild the list of entries passing the QFS window, the minimum number of
/// hits and the user cut.
static void nedRescan()
{
    if (!nedReady())
        return;

    gNedList.clear();

    TTreeFormula *form = nullptr;
    if (!gNedCut.IsNull())
    {
        form = new TTreeFormula("nedCut", gNedCut.Data(), gNedTree);
        if (form->GetNdim() == 0)
        {
            std::cerr << "[ERROR] Invalid cut: " << gNedCut << " (ignored)\n";
            delete form;
            form = nullptr;
            gNedCut = "";
        }
    }

    const Long64_t n = gNedTree->GetEntries();
    for (Long64_t i = 0; i < n; ++i)
    {
        gNedTree->GetEntry(i);
        if ((int)gNedX->size() < gNedMinHits)
            continue;
        if (gNedQFS && (gNedOpa <= kNedOpaMin || gNedOpa >= kNedOpaMax))
            continue;
        if (form && form->EvalInstance(0) == 0)
            continue;

        gNedList.push_back(i);
    }

    delete form;

    if (gNedErelCanvas)
        erelSpectrum(gNedErelBins, gNedErelLo, gNedErelHi);

    std::cout << "[OK] " << gNedList.size() << " / " << n
              << " events pass the selection (minHits = " << gNedMinHits;
    if (gNedQFS)
        std::cout << ", QFS " << kNedOpaMin << " < califa_opa < " << kNedOpaMax;
    if (!gNedCut.IsNull())
        std::cout << ", cut = \"" << gNedCut << "\"";
    std::cout << ")\n";

    gNedPos = 0;
}

// ─── Canvas construction ────────────────────────────────────────────────────

/// (Re)create the axis frames from the current gNedLo / gNedHi ranges.
static void nedMakeFrames()
{
    for (int p = 0; p < 6; ++p)
    {
        double uLo, uHi, vLo, vHi;
        nedPadRange(p, uLo, uHi, vLo, vHi);

        const TString title = TString::Format("%s;%s;%s",
                                              kNedPanelTitle[p],
                                              kNedVarName[kNedU[p]],
                                              kNedVarName[kNedV[p]]);

        delete gNedFrame[p];
        gNedFrame[p] = new TH2F(Form("nedFrame%d", p), title,
                                1, uLo, uHi, 1, vLo, vHi);

        for (TAxis *a : {gNedFrame[p]->GetXaxis(), gNedFrame[p]->GetYaxis()})
        {
            a->SetTitleSize(0.05);
            a->SetLabelSize(0.04);
        }
        gNedFrame[p]->GetXaxis()->SetTitleOffset(1.05);
        gNedFrame[p]->GetYaxis()->SetTitleOffset(1.35);
    }

    delete gNedFrame3D;
    gNedFrame3D = new TH3F("nedFrame3D", "3D view;z [cm];x [cm];y [cm]",
                           1, gNedLo[2], gNedHi[2],
                           1, gNedLo[0], gNedHi[0],
                           1, gNedLo[1], gNedHi[1]);
    for (TAxis *a : {gNedFrame3D->GetXaxis(), gNedFrame3D->GetYaxis(),
                     gNedFrame3D->GetZaxis()})
    {
        a->SetTitleSize(0.05);
        a->SetLabelSize(0.035);
        a->SetTitleOffset(1.8);
    }
}

static void nedBuildCanvas()
{
    gStyle->SetOptStat(0);
    gStyle->SetPalette(kBird);
    gStyle->SetPadTickX(1);
    gStyle->SetPadTickY(1);
    gStyle->SetTitleFontSize(0.06);

    nedMakeFrames();

    if (gNedCanvas)
        return;

    gNedCanvas = new TCanvas("cNeutronEventDisplay",
                             "NeuLAND neutron event display", 1500, 800);
    gNedCanvas->Divide(4, 2, 0.001, 0.001);

    for (int p = 1; p <= 8; ++p)
    {
        TVirtualPad *pad = gNedCanvas->cd(p);
        pad->SetLeftMargin(0.15);
        pad->SetRightMargin(0.04);
        pad->SetBottomMargin(0.13);
        pad->SetTopMargin(0.10);
    }

    if (TPad *pad3d = dynamic_cast<TPad *>(gNedCanvas->cd(8)))
    {
        pad3d->SetTheta(25.);
        pad3d->SetPhi(-50.);
    }
}

// ─── Drawing ────────────────────────────────────────────────────────────────

/// Colour scale (time) drawn at the bottom of the information pad.
static void nedDrawTimeScale()
{
    const int nBox = 40;
    const double x0 = 0.10, x1 = 0.92, y0 = 0.06, y1 = 0.12;

    for (int i = 0; i < nBox; ++i)
    {
        const double f = (i + 0.5) / nBox;
        const double t = gNedLo[3] + f * (gNedHi[3] - gNedLo[3]);
        TBox *b = new TBox(x0 + (x1 - x0) * i / nBox, y0,
                           x0 + (x1 - x0) * (i + 1) / nBox, y1);
        b->SetFillColor(nedTimeColor(t));
        b->SetBit(kCanDelete);
        b->Draw();
    }

    TLatex *lo = new TLatex(x0, y1 + 0.02, Form("%.0f ns", gNedLo[3]));
    TLatex *hi = new TLatex(x1, y1 + 0.02, Form("%.0f ns", gNedHi[3]));
    TLatex *mid = new TLatex(0.5 * (x0 + x1), y0 - 0.05, "hit colour = time");
    for (TLatex *l : {lo, hi, mid})
    {
        l->SetTextSize(0.05);
        l->SetBit(kCanDelete);
    }
    hi->SetTextAlign(31);
    mid->SetTextAlign(21);
    lo->Draw();
    hi->Draw();
    mid->Draw();
}

static void nedDrawInfo()
{
    TPaveText *pt = new TPaveText(0.02, 0.20, 0.98, 0.97, "NDC");
    pt->SetBit(kCanDelete);
    pt->SetFillColor(0);
    pt->SetBorderSize(0);
    pt->SetTextAlign(12);
    pt->SetTextFont(42);
    pt->SetTextSize(0.055);

    if (gNedInList)
        pt->AddText(Form("#bf{entry %lld}   (%d / %d selected)",
                         gNedEntry, gNedPos + 1, (int)gNedList.size()));
    else
        pt->AddText(Form("#bf{entry %lld}   (not in the current selection)",
                         gNedEntry));
    pt->AddText(Form("NeuLAND hits: %d", (int)gNedX->size()));
    pt->AddText("");
    pt->AddText("#bf{First hit (used as neutron)}");
    pt->AddText(Form("x = %.1f   y = %.1f   z = %.1f cm", gNedXf, gNedYf, gNedZf));
    pt->AddText(Form("t = %.2f ns   #beta = %.4f", gNedTf, gNedBetaNeu));
    pt->AddText(Form("p = %.3f GeV/c   (p_{z} = %.3f)", gNedPNeu, gNedPzNeu));
    pt->AddText("");

    const double erel = nedErel();
    if (erel > -900.)
        pt->AddText(Form("#bf{E_{rel} = %.3f MeV}", erel * 1000.));
    else
        pt->AddText("#bf{E_{rel} = not available}");

    pt->AddText("");
    pt->AddText("#bf{Fragment}");
    pt->AddText(Form("A/Q = %.3f   Z = %.2f   #beta = %.4f",
                     gNedAoQ, gNedZfrag, gNedBetaFrag));
    pt->AddText(Form("CALIFA opa = %.3f rad (%.1f deg)",
                     gNedOpa, gNedOpa * TMath::RadToDeg()));
    if (gNedQFS)
        pt->AddText(Form("#it{QFS window %.2f - %.2f rad}", kNedOpaMin, kNedOpaMax));
    pt->Draw();

    nedDrawTimeScale();
}

static void nedDraw3D()
{
    gNedFrame3D->Draw();

    const int nHits = (int)gNedX->size();
    if (nHits > 0)
    {
        TPolyMarker3D *pm = new TPolyMarker3D(nHits);
        pm->SetBit(kCanDelete);
        for (int h = 0; h < nHits; ++h)
            pm->SetPoint(h, (*gNedZ)[h], (*gNedX)[h], (*gNedY)[h]);
        pm->SetMarkerStyle(20);
        pm->SetMarkerSize(1.1);
        pm->SetMarkerColor(kAzure + 2);
        pm->Draw();
    }

    TPolyMarker3D *first = new TPolyMarker3D(1);
    first->SetBit(kCanDelete);
    first->SetPoint(0, gNedZf, gNedXf, gNedYf);
    first->SetMarkerStyle(29);
    first->SetMarkerSize(2.2);
    first->SetMarkerColor(kRed + 1);
    first->Draw();
}

/// Draw the currently loaded entry.
static void nedDrawEvent()
{
    if (!gNedCanvas)
        return;

    for (int p = 0; p < 6; ++p)
    {
        TVirtualPad *pad = gNedCanvas->cd(p < 3 ? p + 1 : p + 2);
        pad->Clear();
        pad->SetGrid();

        gNedFrame[p]->Draw("AXIS");

        // projected neutron trajectory from the target (0,0,0)
        if (p < 3)
        {
            TLine *tr = new TLine(nedVar(kNedU[p], 0, 0, 0, 0),
                                  nedVar(kNedV[p], 0, 0, 0, 0),
                                  nedVar(kNedU[p], gNedXf, gNedYf, gNedZf, gNedTf),
                                  nedVar(kNedV[p], gNedXf, gNedYf, gNedZf, gNedTf));
            tr->SetLineColor(kRed + 1);
            tr->SetLineStyle(2);
            tr->SetBit(kCanDelete);
            tr->Draw();
        }

        // all hits of the event, coloured by time
        for (size_t h = 0; h < gNedX->size(); ++h)
        {
            const double x = (*gNedX)[h], y = (*gNedY)[h];
            const double z = (*gNedZ)[h], t = (*gNedT)[h];
            TMarker *m = new TMarker(nedVar(kNedU[p], x, y, z, t),
                                     nedVar(kNedV[p], x, y, z, t), 20);
            m->SetMarkerColor(nedTimeColor(t));
            m->SetMarkerSize(1.6);
            m->SetBit(kCanDelete);
            m->Draw();
        }

        // first hit highlighted
        TMarker *f = new TMarker(nedVar(kNedU[p], gNedXf, gNedYf, gNedZf, gNedTf),
                                 nedVar(kNedV[p], gNedXf, gNedYf, gNedZf, gNedTf), 29);
        f->SetMarkerColor(kRed + 1);
        f->SetMarkerSize(2.6);
        f->SetBit(kCanDelete);
        f->Draw();

        pad->Modified();
    }

    TVirtualPad *info = gNedCanvas->cd(4);
    info->Clear();
    nedDrawInfo();
    info->Modified();

    TVirtualPad *pad3d = gNedCanvas->cd(8);
    pad3d->Clear();
    nedDraw3D();
    pad3d->Modified();

    gNedCanvas->Update();

    nedSetStatus(gNedInList
                     ? TString::Format("entry %lld   %d / %d selected   hits: %d",
                                       gNedEntry, gNedPos + 1,
                                       (int)gNedList.size(), (int)gNedX->size())
                     : TString::Format("entry %lld   not in the selection   hits: %d",
                                       gNedEntry, (int)gNedX->size()));
}

// ─── Erel spectrum canvas ───────────────────────────────────────────────────

/// Refill the three Erel histograms from the events of the current selection.
static void nedFillErelSpectrum()
{
    if (!nedReady() || !gNedErelZ)
        return;

    gNedErelZ->Reset();
    gNedErelT->Reset();
    gNedErelSame->Reset();

    const Long64_t keep = gNedEntry;
    int nSame = 0;

    for (size_t i = 0; i < gNedList.size(); ++i)
    {
        gNedTree->GetEntry(gNedList[i]);

        // (1) neutron of the tree: first hit in z
        const double erelZ = nedErel();
        if (erelZ > -900.)
            gNedErelZ->Fill(erelZ * 1000.);

        // (2) earliest hit in time
        const int idxT = nedFirstTimeIndex();
        if (idxT >= 0)
        {
            const double erelT = nedErelFromHit((*gNedX)[idxT], (*gNedY)[idxT],
                                                (*gNedZ)[idxT], (*gNedT)[idxT]);
            if (erelT > -900.)
                gNedErelT->Fill(erelT * 1000.);
        }

        // (3) events where both criteria pick the same hit
        if (idxT >= 0 && idxT == nedFirstHitIndex())
        {
            ++nSame;
            if (erelZ > -900.)
                gNedErelSame->Fill(erelZ * 1000.);
        }
    }

    std::cout << "[OK] Erel spectrum over " << gNedList.size() << " events: "
              << "first in z <Erel> = " << gNedErelZ->GetMean() << " MeV, "
              << "first in time <Erel> = " << gNedErelT->GetMean() << " MeV, "
              << "same hit in " << nSame << " events ("
              << (gNedList.empty() ? 0. : 100. * nSame / gNedList.size())
              << " %)\n";

    if (keep >= 0)
        gNedTree->GetEntry(keep);
}

/// Second canvas: Erel of the current selection for the three ways of picking
/// the neutron hit.
void erelSpectrum(int nbins, double lo, double hi)
{
    if (!nedReady())
        return;

    gNedErelBins = nbins;
    gNedErelLo = lo;
    gNedErelHi = hi;

    delete gNedErelZ;
    delete gNedErelT;
    delete gNedErelSame;

    gNedErelZ = new TH1F("nedErelZ", "E_{rel};E_{rel} [MeV];counts",
                         nbins, lo, hi);
    gNedErelT = new TH1F("nedErelT", "E_{rel};E_{rel} [MeV];counts",
                         nbins, lo, hi);
    gNedErelSame = new TH1F("nedErelSame", "E_{rel};E_{rel} [MeV];counts",
                            nbins, lo, hi);

    gNedErelZ->SetLineColor(kBlack);
    gNedErelT->SetLineColor(kRed + 1);
    gNedErelSame->SetLineColor(kAzure + 2);
    for (TH1F *h : {gNedErelZ, gNedErelT, gNedErelSame})
    {
        h->SetLineWidth(2);
        h->SetStats(0);
    }

    nedFillErelSpectrum();

    if (!gNedErelCanvas)
    {
        gNedErelCanvas = new TCanvas("cNeutronErel",
                                     "Erel of the selected events", 800, 600);
        gNedErelCanvas->SetLeftMargin(0.13);
        gNedErelCanvas->SetTicks(1, 1);
    }

    gNedErelCanvas->cd();
    gNedErelCanvas->Clear();

    gNedErelZ->SetMaximum(1.25 * std::max(gNedErelZ->GetMaximum(),
                                          gNedErelT->GetMaximum()));
    gNedErelZ->Draw("HIST");
    gNedErelT->Draw("HIST SAME");
    gNedErelSame->Draw("HIST SAME");

    TLegend *leg = new TLegend(0.46, 0.66, 0.96, 0.88);
    leg->SetBit(kCanDelete);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetTextSize(0.030);
    leg->SetHeader("hit used as the neutron");
    leg->AddEntry(gNedErelZ,
                  Form("first in z: %d ev, %.2f MeV",
                       (int)gNedErelZ->GetEntries(), gNedErelZ->GetMean()),
                  "l");
    leg->AddEntry(gNedErelT,
                  Form("first in time: %d ev, %.2f MeV",
                       (int)gNedErelT->GetEntries(), gNedErelT->GetMean()),
                  "l");
    leg->AddEntry(gNedErelSame,
                  Form("same hit: %d ev, %.2f MeV",
                       (int)gNedErelSame->GetEntries(), gNedErelSame->GetMean()),
                  "l");
    leg->Draw();

    TLatex *sel = new TLatex(0.16, 0.86,
                             Form("#it{%s}#it{, minHits = %d%s}",
                                  gNedQFS ? Form("QFS %.2f - %.2f rad", kNedOpaMin, kNedOpaMax)
                                          : "no QFS window",
                                  gNedMinHits,
                                  gNedCut.IsNull() ? "" : Form(", %s", gNedCut.Data())));
    sel->SetBit(kCanDelete);
    sel->SetNDC();
    sel->SetTextSize(0.028);
    sel->Draw();

    gNedErelCanvas->Update();
}

// ─── Navigation (also used as GUI slots) ────────────────────────────────────

/// Show a given tree entry, whether or not it passes the current selection.
void showEvent(Long64_t entry)
{
    if (!nedReady())
        return;
    if (entry < 0 || entry >= gNedTree->GetEntries())
    {
        std::cerr << "[ERROR] Entry " << entry << " out of range [0, "
                  << gNedTree->GetEntries() - 1 << "]\n";
        return;
    }

    gNedTree->GetEntry(entry);
    gNedEntry = entry;

    auto it = std::find(gNedList.begin(), gNedList.end(), entry);
    gNedInList = (it != gNedList.end());
    if (gNedInList)
        gNedPos = (int)(it - gNedList.begin());

    nedDrawEvent();
}

/// Show the i-th event of the current selection.
static void nedShowIndex(int i)
{
    if (!nedReady())
        return;
    if (gNedList.empty())
    {
        std::cerr << "[ERROR] No event passes the current selection.\n";
        return;
    }
    gNedPos = ((i % (int)gNedList.size()) + (int)gNedList.size()) % (int)gNedList.size();
    showEvent(gNedList[gNedPos]);
}

void nextEvent() { nedShowIndex(gNedPos + 1); }
void prevEvent() { nedShowIndex(gNedPos - 1); }

void randomEvent()
{
    if (gNedList.empty())
        return;
    static TRandom3 rng(0);
    nedShowIndex((int)(rng.Uniform(0, gNedList.size())));
}

void printEvent()
{
    if (!nedReady() || gNedEntry < 0)
        return;

    std::cout << "\n─── entry " << gNedEntry << " ─ " << gNedX->size()
              << " NeuLAND hits ───\n"
              << std::setw(6) << "hit" << std::setw(9) << "paddle"
              << std::setw(11) << "x [cm]" << std::setw(11) << "y [cm]"
              << std::setw(12) << "z [cm]" << std::setw(11) << "t [ns]" << "\n";

    for (size_t h = 0; h < gNedX->size(); ++h)
        std::cout << std::fixed << std::setprecision(2)
                  << std::setw(6) << h
                  << std::setw(9) << (*gNedPaddle)[h]
                  << std::setw(11) << (*gNedX)[h]
                  << std::setw(11) << (*gNedY)[h]
                  << std::setw(12) << (*gNedZ)[h]
                  << std::setw(11) << (*gNedT)[h] << "\n";

    std::cout << "first hit: x = " << gNedXf << "  y = " << gNedYf
              << "  z = " << gNedZf << "  t = " << gNedTf
              << "  beta = " << gNedBetaNeu << "  p = " << gNedPNeu << " GeV/c\n\n";
}

void saveEvent(const char *fileName)
{
    if (!gNedCanvas)
        return;
    const TString name = fileName ? TString(fileName)
                                  : TString::Format("neutron_event_%lld.png", gNedEntry);
    gNedCanvas->SaveAs(name);
}

// ─── Selection control ──────────────────────────────────────────────────────

void setMinHits(int n)
{
    gNedMinHits = std::max(1, n);
    nedRescan();
    nedShowIndex(0);
}

void applyCut(const char *cut)
{
    gNedCut = (cut ? cut : "");
    nedRescan();
    nedShowIndex(0);
}

/// Manual axis range for one coordinate: var is "x", "y", "z" or "t".
void setRange(const char *var, double lo, double hi)
{
    const TString v(var);
    int k = -1;
    if (v == "x")
        k = 0;
    else if (v == "y")
        k = 1;
    else if (v == "z")
        k = 2;
    else if (v == "t")
        k = 3;

    if (k < 0 || lo >= hi)
    {
        std::cerr << "[ERROR] setRange(\"x\"|\"y\"|\"z\"|\"t\", lo, hi) with lo < hi\n";
        return;
    }

    gNedLo[k] = lo;
    gNedHi[k] = hi;
    nedMakeFrames();
    if (gNedEntry >= 0)
        nedDrawEvent();
}

/// Keep only quasi-free scattering events (1.25 < califa_opa < 1.65 rad).
void setQFS(bool on)
{
    gNedQFS = on;
    nedRescan();
    nedShowIndex(0);
}

void playEvents(int nEvents = 20, double seconds = 1.0)
{
    for (int i = 0; i < nEvents; ++i)
    {
        nextEvent();
        gSystem->ProcessEvents();
        gSystem->Sleep((ULong_t)(seconds * 1000));
    }
}

// ─── GUI ────────────────────────────────────────────────────────────────────

void nedGuiGoto()
{
    if (gNedGotoBox)
        showEvent((Long64_t)gNedGotoBox->GetNumberEntry()->GetIntNumber());
}

void nedGuiSave() { saveEvent(nullptr); }

void nedCloseGui()
{
    if (gNedGui)
    {
        gNedGui->CloseWindow();
        gNedGui = nullptr;
        gNedGotoBox = nullptr;
        gNedStatus = nullptr;
    }
}

/// Keyboard shortcuts while the mouse is over the canvas.
void nedCanvasEvent(Int_t event, Int_t px, Int_t /*py*/, TObject * /*sel*/)
{
    if (event != kKeyPress)
        return;

    switch (px)
    {
    case 'n':
        nextEvent();
        break;
    case 'p':
        prevEvent();
        break;
    case 'r':
        randomEvent();
        break;
    case 'd':
        printEvent();
        break;
    case 's':
        saveEvent(nullptr);
        break;
    case 'h':
        nedHelp();
        break;
    default:
        break;
    }
}

static void nedBuildGui()
{
    if (gROOT->IsBatch() || gNedGui)
        return;

    gNedGui = new TGMainFrame(gClient->GetRoot(), 520, 90);
    gNedGui->SetWindowName("Neutron event display");

    TGHorizontalFrame *row = new TGHorizontalFrame(gNedGui, 520, 40);

    TGTextButton *bPrev = new TGTextButton(row, "  << Prev  ");
    TGTextButton *bNext = new TGTextButton(row, "  Next >>  ");
    TGTextButton *bRand = new TGTextButton(row, " Random ");
    TGTextButton *bPrint = new TGTextButton(row, " Print ");
    TGTextButton *bSave = new TGTextButton(row, " Save ");

    bPrev->Connect("Clicked()", 0, 0, "prevEvent()");
    bNext->Connect("Clicked()", 0, 0, "nextEvent()");
    bRand->Connect("Clicked()", 0, 0, "randomEvent()");
    bPrint->Connect("Clicked()", 0, 0, "printEvent()");
    bSave->Connect("Clicked()", 0, 0, "nedGuiSave()");

    for (TGTextButton *b : {bPrev, bNext, bRand, bPrint, bSave})
        row->AddFrame(b, new TGLayoutHints(kLHintsLeft | kLHintsCenterY, 3, 3, 3, 3));

    row->AddFrame(new TGLabel(row, "entry:"),
                  new TGLayoutHints(kLHintsLeft | kLHintsCenterY, 8, 3, 3, 3));
    gNedGotoBox = new TGNumberEntry(row, 0, 8, -1,
                                    TGNumberFormat::kNESInteger,
                                    TGNumberFormat::kNEANonNegative);
    row->AddFrame(gNedGotoBox,
                  new TGLayoutHints(kLHintsLeft | kLHintsCenterY, 3, 3, 3, 3));

    TGTextButton *bGo = new TGTextButton(row, " Go ");
    bGo->Connect("Clicked()", 0, 0, "nedGuiGoto()");
    row->AddFrame(bGo, new TGLayoutHints(kLHintsLeft | kLHintsCenterY, 3, 3, 3, 3));

    gNedGui->AddFrame(row, new TGLayoutHints(kLHintsExpandX));

    gNedStatus = new TGLabel(gNedGui, "                                        ");
    gNedGui->AddFrame(gNedStatus,
                      new TGLayoutHints(kLHintsLeft | kLHintsExpandX, 6, 6, 2, 4));

    gNedGui->Connect("CloseWindow()", 0, 0, "nedCloseGui()");
    gNedGui->MapSubwindows();
    gNedGui->Resize(gNedGui->GetDefaultSize());
    gNedGui->MapWindow();
}

// ─── Help ───────────────────────────────────────────────────────────────────

void nedHelp()
{
    std::cout << "\n─── neutron event display ───────────────────────────────\n"
              << "  nextEvent()              next selected event\n"
              << "  prevEvent()              previous selected event\n"
              << "  randomEvent()            random selected event\n"
              << "  showEvent(entry)         jump to a given tree entry\n"
              << "  printEvent()             hit table of the current event\n"
              << "  saveEvent(\"name.png\")    save the canvas\n"
              << "  setMinHits(n)            keep events with >= n hits\n"
              << "  applyCut(\"beta_neu>0.6\") TTree cut on scalar branches\n"
              << "  setQFS(false)            drop the 1.25 < opa < 1.65 window\n"
              << "  setRange(\"t\", 50, 80)     axis range of x, y, z or t\n"
              << "  loadOffsets(\"23O1n.txt\")  offsets used for E_rel\n"
              << "  erelSpectrum(100,-5,10)  redraw the E_rel canvas (bins, range)\n"
              << "  playEvents(20, 0.5)      slideshow: n events, seconds each\n"
              << "  canvas keys: n/p/r/d/s/h\n"
              << "─────────────────────────────────────────────────────────\n\n";
}

// ═══════════════════════════════════════════════════════════════════════════
//  MAIN
// ═══════════════════════════════════════════════════════════════════════════

void neutronEventDisplay(const char *fileName = kNedDefaultFile,
                         int minHits = 1,
                         const char *cut = "",
                         bool qfsOnly = true)
{
    TH1::AddDirectory(kFALSE);
    nedUnload();

    if (gNedOffsetsFile.IsNull())
        loadOffsets();

    gNedFile = TFile::Open(fileName, "READ");
    if (!gNedFile || gNedFile->IsZombie())
    {
        std::cerr << "[ERROR] Cannot open " << fileName << "\n";
        return;
    }

    gNedTree = dynamic_cast<TTree *>(gNedFile->Get(kNedTreeName));
    if (!gNedTree)
    {
        std::cerr << "[ERROR] Tree " << kNedTreeName << " not found in "
                  << fileName << "\n";
        nedUnload();
        return;
    }
    if (!gNedTree->GetBranch("n_x"))
    {
        std::cerr << "[ERROR] No NeuLAND hit branches (n_x, n_y, n_z, n_t) in "
                  << fileName << ". Was it produced with hasNeutrons = true?\n";
        nedUnload();
        return;
    }

    std::cout << "[OK] " << fileName << " : " << gNedTree->GetEntries()
              << " entries\n";

    nedSetBranches();
    nedComputeRanges();
    nedBuildCanvas();

    gNedMinHits = std::max(1, minHits);
    gNedCut = cut ? cut : "";
    gNedQFS = qfsOnly;
    nedRescan();

    nedBuildGui();
    if (!gROOT->IsBatch())
        gNedCanvas->Connect("ProcessedEvent(Int_t,Int_t,Int_t,TObject*)", 0, 0,
                            "nedCanvasEvent(Int_t,Int_t,Int_t,TObject*)");

    nedShowIndex(0);
    erelSpectrum();
    nedHelp();
}
