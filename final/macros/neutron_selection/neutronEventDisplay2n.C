// ═══════════════════════════════════════════════════════════════════════════
//  neutronEventDisplay2n.C
//
//  Event display and Erel spectra for the 2n breakup channel, the counterpart
//  of neutronEventDisplay.C. It reads the tree written by
//  final/macros/data_reduction/eventFilter.C (tree "FilterDataTree") for the
//  25F -> 22O + 2n selection.
//
//  The NeuLAND hits of each event are grouped into clusters by causality, with
//  causalClustering.h (see the long comment there): two hits are joined when,
//  within the detector resolution, one particle emitted at the earlier hit
//  could have produced the later one, either a scattered neutron / recoil no
//  faster than the incident neutron, or a prompt gamma travelling at c. The
//  clusters are grown by single linkage under that relation. The resolutions
//  (setResolutions), the width of the window (setCausalSigma), the gamma
//  branch (setGammaBranch, setBeamBetaWindow) and the strict-to-first-hit
//  systematic (setStrictToFirst) are all tunable; causalAcceptance() prints
//  the n-n relative energy below which a genuine pair is merged.
//
//  The old clustering, a fixed sphere (two hits closer than dR in space and dT
//  in time, 50 cm and 5 ns by default, setClusterCuts), is kept for comparison
//  and selected with setUseCausal(false). It fails in three ways at once: it
//  merges genuine pairs closer than dR, it splits real scattering chains
//  slower than dR / dT, and it lets a gamma emitted at the first interaction
//  through as a second neutron.
//
//  Two clusters alone do not mean two neutrons: one neutron that scatters and
//  interacts again further on also leaves two groups of hits, and then the
//  second neutron was simply missed. With the sphere this was handled by a
//  causality test between the clusters (ned2TwoNeutrons): the apparent
//  velocity between their earliest hits, dr / (c dt), had to exceed the
//  neutron velocity. With causal clustering that test is subsumed and it is
//  OFF by default (g2NBetaRelCut = 0): two clusters that survive the causal
//  clustering are, by construction, not connectable by any particle through
//  ANY pair of their hits, which is strictly stronger than comparing only
//  their two earliest hits. setCausalityCut(-1) switches the old test back on
//  (on top of either clustering), which is what the sphere selection used.
//
//  Every cluster must also look like a neutron (causal::neutronLike): the
//  velocity implied by its first hit in time, beta0 = |r| / (c t), has to lie
//  inside a window, 0.55 - 0.90 by default (setNeutronBetaWindow, lo >= hi
//  switches it off). Without it a cluster of gammas (beta0 ~ 1) enters as a
//  neutron and pushes Erel up.
//
//  The exact selection used before the causal clustering is recovered with
//      setUseCausal(false); setCausalityCut(-1); setNeutronBetaWindow(0, 0);
//
//  The surviving events with exactly two clusters are the 2n candidates. Each
//  cluster gives one neutron, and Erel is built from the invariant mass of
//  fragment + n1 + n2 in two ways: taking the first hit in z of each cluster,
//  and taking the first hit in time of each cluster. Both numbers appear in
//  the panel of the display and as spectra of the second canvas, which shows a
//  pad per definition and holds a
//  third one with the events where the two criteria pick the same hit in both
//  clusters, where the two definitions give the same energy, and a fourth one
//  that keeps, event by event, whichever of the two choices gives the smaller
//  relative energy.
//
//  Events with a single cluster, that is one neutron seen and one missed, are
//  outside the selection from the start, as are events with more clusters.
//
//  Six projections of the hit coordinates are shown per event:
//      y vs x   (front view)      t vs x
//      x vs z   (top view)        t vs y
//      y vs z   (side view)       t vs z
//  plus a 3D view and an information panel. Hits are coloured by cluster, a
//  filled star marks the first hit in z of a cluster and an open square the
//  first hit in time.
//
//  Only quasi-free scattering events are shown by default, that is
//  1.25 < califa_opa < 1.65 rad. Use setQFS2n(false) to see all of them.
//
//  Usage:
//      root -l 'neutronEventDisplay2n.C'
//      root -l 'neutronEventDisplay2n.C("myfile.root")'
//
//  Then, at the ROOT prompt (or with the buttons of the control panel):
//      nextEvent2n()        prevEvent2n()      randomEvent2n()
//      showEvent2n(entry)   printEvent2n()     saveEvent2n("name.png")
//      setUseCausal(true)   setCausalSigma(2)  setResolutions(4, 0.15)
//      setGammaBranch(true) setBeamBetaWindow(0.70, 0.80)
//      setStrictToFirst(false)                 causalAcceptance(1300, 0.75)
//      setNeutronBetaWindow(0.55, 0.90)        (lo >= hi switches it off)
//      setClusterCuts(50, 5)  (sphere only)    setRequiredClusters(2)
//      setCausalityCut(0)   (<0 event neutron velocity, 0 off, >0 fixed)
//      setMinHits2n(4)      applyCut2n("beta_neu>0.6")
//          (setMinHits2n asks for that many hits in every cluster)
//      setQFS2n(false)      setRange2n("t", 55, 80)
//      loadOffsets2n("23O1n.txt")              erelSpectrum2n(100, -5, 20)
//      ned2nHelp()
//
//  With the mouse over the canvas: n = next, p = previous, r = random,
//  d = dump the hit table, s = save the canvas, h = help.
//
//  The names all carry the 2n suffix so that this macro and the 1n one can be
//  loaded in the same ROOT session.
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
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TMath.h"
#include "TMarker.h"
#include "TLine.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TPaveText.h"
#include "TPolyMarker3D.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include "TRandom3.h"
#include "TString.h"

#include "TGClient.h"
#include "TGFrame.h"
#include "TGLayout.h"
#include "TGButton.h"
#include "TGLabel.h"
#include "TGNumberEntry.h"

#include "causalClustering.h"

// ─── Configuration ──────────────────────────────────────────────────────────

static const char *k2NDefaultFile =
    "/nucl_lustre/pablogrusell/g249/g249_analysis/results/dataFiles/"
    "data_22O_with_neutrons.root";
static const char *k2NTreeName = "FilterDataTree";

// quasi-free scattering selection on the CALIFA opening angle [rad]
static constexpr double k2NOpaMin = 1.25;
static constexpr double k2NOpaMax = 1.65;

static constexpr double k2NMn = 0.939565;  // neutron mass [GeV]
static constexpr double k2NC = 29.9792458; // speed of light [cm/ns]

static const char *k2NDefaultOffsets = "22O2n.txt";
static const char *k2NSettingsDir =
    "/nucl_lustre/pablogrusell/g249/g249_analysis/final/settings/";

// Variable indices used by the panel definitions: 0=x, 1=y, 2=z, 3=t
static const int k2NU[6] = {0, 2, 2, 0, 1, 2}; // horizontal axis of each panel
static const int k2NV[6] = {1, 0, 1, 3, 3, 3}; // vertical axis of each panel

static const char *k2NVarName[4] = {"x [cm]", "y [cm]", "z [cm]", "t [ns]"};

static const char *k2NPanelTitle[6] = {
    "Front view (beam into page): y vs x",
    "Top view: x vs z",
    "Side view: y vs z",
    "Time vs x",
    "Time vs y",
    "Time vs z"};

static const int k2NClusterColor[8] = {kAzure + 2, kRed + 1, kGreen + 2,
                                       kMagenta + 1, kOrange + 7, kCyan + 2,
                                       kYellow + 2, kGray + 2};

// ─── State ──────────────────────────────────────────────────────────────────

static TFile *g2NFile = nullptr;
static TTree *g2NTree = nullptr;
static TCanvas *g2NCanvas = nullptr;

// branch buffers (all NeuLAND hits of the event)
static std::vector<double> *g2NX = nullptr;
static std::vector<double> *g2NY = nullptr;
static std::vector<double> *g2NZ = nullptr;
static std::vector<double> *g2NT = nullptr;
static std::vector<int> *g2NPaddle = nullptr;

// branch buffers (fragment, and the single neutron kept by the filter)
static double g2NXf = 0, g2NYf = 0, g2NZf = 0, g2NTf = 0;
static double g2NBetaNeu = 0, g2NPNeu = 0;
static double g2NAoQ = 0, g2NZfrag = 0, g2NBetaFrag = 0, g2NOpa = 0;
static double g2NMFrag = 0;
static double g2NPxFrag = 0, g2NPyFrag = 0, g2NPzFrag = 0;

// Erel calibration, read from final/settings/<file> (fragment x, fragment y,
// neutron x, neutron y, beta matching)
static double g2NFragOff[2] = {0., 0.};
static double g2NNeuOff[2] = {0., 0.};
static double g2NBetaMatch = 0.;
static TString g2NOffsetsFile = "";

// clustering
static causal::Config g2NCausal; // resoluciones y ventana de causalidad
static bool g2NUseCausal = true; // false vuelve a la esfera

// sphere clustering, only used with g2NUseCausal = false
static double g2NdR = 50.; // maximum distance between two hits of a cluster [cm]
static double g2NdT = 5.;  // maximum time difference between them [ns]

// causality test between the two clusters: a negative value compares the
// apparent velocity with the neutron velocity of the event itself, a positive
// one with that fixed value, and zero switches the test off. Off by default:
// the causal clustering already guarantees that no particle connects any pair
// of hits of two different clusters (see the header of the macro)
static double g2NBetaRelCut = 0.;

// cluster quality: window on beta0, the velocity implied by the first hit in
// time of every cluster (causal::neutronLike). lo >= hi switches it off
static double g2NNeuBetaLo = 0.55;
static double g2NNeuBetaHi = 0.90;

// selection
static std::vector<Long64_t> g2NList; // tree entries passing the selection
static int g2NPos = 0;                // current position inside g2NList
static Long64_t g2NEntry = -1;        // current tree entry
static bool g2NInList = false;        // is that entry part of the selection?
static int g2NMinHits = 1;            // minimum number of hits in every cluster
static int g2NReqClusters = 2;        // 0 keeps every cluster multiplicity
static TString g2NCut = "";
static bool g2NQFS = true;

// axis limits for x, y, z, t
static double g2NLo[4] = {0, 0, 0, 0};
static double g2NHi[4] = {0, 0, 0, 0};

// drawing objects that survive between events
static TH2F *g2NFrame[6] = {nullptr, nullptr, nullptr, nullptr, nullptr, nullptr};
static TH3F *g2NFrame3D = nullptr;

// Erel spectrum canvas
static TCanvas *g2NErelCanvas = nullptr;
static TH1F *g2NErelPos = nullptr;  // first hit in z of each cluster
static TH1F *g2NErelTime = nullptr; // first hit in time of each cluster
static TH1F *g2NErelSame = nullptr; // events where both criteria agree
static TH1F *g2NErelMin = nullptr;  // the smaller of the two energies
static int g2NErelBins = 100;
static double g2NErelLo = -5.;
static double g2NErelHi = 20.;

// GUI
static TGMainFrame *g2NGui = nullptr;
static TGNumberEntry *g2NGotoBox = nullptr;
static TGLabel *g2NStatus = nullptr;

// forward declarations (needed by the interpreted GUI slots)
void nextEvent2n();
void prevEvent2n();
void randomEvent2n();
void showEvent2n(Long64_t entry);
void printEvent2n();
void saveEvent2n(const char *fileName = nullptr);
static void ned2MakeFrames();
static void ned2DrawEvent();
void erelSpectrum2n(int nbins = 100, double lo = -5., double hi = 20.);
void ned2nHelp();

/// One group of hits, with the two candidate "first" hits of the group.
struct N2Cluster
{
    std::vector<int> hits;
    int firstZ = -1; // hit with the smallest z
    int firstT = -1; // earliest hit
    double zMin = 0.;
    double tMin = 0.;
};

static double ned2HitBeta(int hit);

// ─── Small helpers ──────────────────────────────────────────────────────────

/// Drop the current file and leave the display in a well defined idle state.
static void ned2Unload()
{
    if (g2NFile)
        g2NFile->Close();
    g2NFile = nullptr;
    g2NTree = nullptr;
    g2NList.clear();
    g2NEntry = -1;
    g2NInList = false;
}

static bool ned2Ready()
{
    if (!g2NTree)
    {
        std::cerr << "[ERROR] No tree loaded. Run neutronEventDisplay2n() first.\n";
        return false;
    }
    return true;
}

/// Value of variable `v` (0=x, 1=y, 2=z, 3=t) for one hit.
static double ned2Var(int v, double x, double y, double z, double t)
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

static void ned2PadRange(int p, double &uLo, double &uHi, double &vLo, double &vHi)
{
    uLo = g2NLo[k2NU[p]];
    uHi = g2NHi[k2NU[p]];
    vLo = g2NLo[k2NV[p]];
    vHi = g2NHi[k2NV[p]];
}

static Int_t ned2ClusterColor(int i)
{
    return k2NClusterColor[i % 8];
}

static void ned2SetStatus(const TString &text)
{
    if (g2NStatus)
    {
        g2NStatus->SetText(text.Data());
        if (g2NGui)
            g2NGui->Layout();
    }
}

// ─── Tree access ────────────────────────────────────────────────────────────

static void ned2SetBranches()
{
    g2NTree->SetBranchAddress("n_x", &g2NX);
    g2NTree->SetBranchAddress("n_y", &g2NY);
    g2NTree->SetBranchAddress("n_z", &g2NZ);
    g2NTree->SetBranchAddress("n_t", &g2NT);
    g2NTree->SetBranchAddress("n_paddle", &g2NPaddle);

    g2NTree->SetBranchAddress("x_neu_hit", &g2NXf);
    g2NTree->SetBranchAddress("y_neu_hit", &g2NYf);
    g2NTree->SetBranchAddress("z_neu_hit", &g2NZf);
    g2NTree->SetBranchAddress("tof_neuland", &g2NTf);

    g2NTree->SetBranchAddress("beta_neu", &g2NBetaNeu);
    g2NTree->SetBranchAddress("p_neu", &g2NPNeu);

    g2NTree->SetBranchAddress("M_frag", &g2NMFrag);
    g2NTree->SetBranchAddress("px_frag", &g2NPxFrag);
    g2NTree->SetBranchAddress("py_frag", &g2NPyFrag);
    g2NTree->SetBranchAddress("pz_frag", &g2NPzFrag);

    g2NTree->SetBranchAddress("AoQ_frag", &g2NAoQ);
    g2NTree->SetBranchAddress("Z_frag_est", &g2NZfrag);
    g2NTree->SetBranchAddress("beta_frag", &g2NBetaFrag);
    g2NTree->SetBranchAddress("califa_opa", &g2NOpa);
}

/// Robust global x/y/z/t ranges: 0.5 % - 99.5 % quantiles of all hits plus a
/// margin, because a few hits carry meaningless times.
static void ned2ComputeRanges()
{
    std::vector<double> val[4];
    const Long64_t n = g2NTree->GetEntries();
    for (Long64_t i = 0; i < n; ++i)
    {
        g2NTree->GetEntry(i);
        for (size_t h = 0; h < g2NX->size(); ++h)
        {
            val[0].push_back((*g2NX)[h]);
            val[1].push_back((*g2NY)[h]);
            val[2].push_back((*g2NZ)[h]);
            val[3].push_back((*g2NT)[h]);
        }
    }

    for (int k = 0; k < 4; ++k)
    {
        if (val[k].empty())
        {
            g2NLo[k] = 0;
            g2NHi[k] = 1;
            continue;
        }

        std::sort(val[k].begin(), val[k].end());
        const size_t last = val[k].size() - 1;
        const double lo = val[k][(size_t)(0.005 * last)];
        const double hi = val[k][(size_t)(0.995 * last)];
        const double margin = 0.08 * (hi - lo) + 1e-6;
        g2NLo[k] = lo - margin;
        g2NHi[k] = hi + margin;
    }
}

// ─── Clustering ─────────────────────────────────────────────────────────────

/// Single linkage clustering of the hits of the current event: two hits are
/// joined when they are within g2NdR in space and g2NdT in time. Clusters come
/// back ordered by the z of their first hit.
static std::vector<N2Cluster> ned2FindClustersSphere()
{
    std::vector<N2Cluster> clusters;
    const int n = (int)g2NX->size();
    if (n == 0)
        return clusters;

    std::vector<bool> used(n, false);

    for (int seed = 0; seed < n; ++seed)
    {
        if (used[seed])
            continue;

        N2Cluster c;
        std::vector<int> stack{seed};
        used[seed] = true;

        while (!stack.empty())
        {
            const int a = stack.back();
            stack.pop_back();
            c.hits.push_back(a);

            for (int b = 0; b < n; ++b)
            {
                if (used[b])
                    continue;

                const double dx = (*g2NX)[a] - (*g2NX)[b];
                const double dy = (*g2NY)[a] - (*g2NY)[b];
                const double dz = (*g2NZ)[a] - (*g2NZ)[b];

                if (std::sqrt(dx * dx + dy * dy + dz * dz) < g2NdR &&
                    std::abs((*g2NT)[a] - (*g2NT)[b]) < g2NdT)
                {
                    used[b] = true;
                    stack.push_back(b);
                }
            }
        }

        for (int h : c.hits)
        {
            if (c.firstZ < 0 || (*g2NZ)[h] < c.zMin)
            {
                c.firstZ = h;
                c.zMin = (*g2NZ)[h];
            }
            if (c.firstT < 0 || (*g2NT)[h] < c.tMin)
            {
                c.firstT = h;
                c.tMin = (*g2NT)[h];
            }
        }

        clusters.push_back(c);
    }

    std::sort(clusters.begin(), clusters.end(),
              [](const N2Cluster &a, const N2Cluster &b)
              { return a.zMin < b.zMin; });

    return clusters;
}

/// Clusters of the current event, causal (causalClustering.h) or sphere
/// depending on g2NUseCausal. Causal clusters come back ordered by the time of
/// their first interaction, sphere clusters by z: Erel does not care, but the
/// cluster labels of the display do.
static std::vector<N2Cluster> ned2FindClusters()
{
    if (!g2NUseCausal)
        return ned2FindClustersSphere();

    const std::vector<causal::Cluster> cc =
        causal::findClusters(*g2NX, *g2NY, *g2NZ, *g2NT, g2NCausal);

    std::vector<N2Cluster> out;
    out.reserve(cc.size());
    for (const causal::Cluster &c : cc)
    {
        N2Cluster n;
        n.hits = c.hits;
        n.firstZ = c.firstZ;
        n.firstT = c.firstT;
        n.zMin = c.zMin;
        n.tMin = c.tMin;
        out.push_back(n);
    }
    return out;
}

/// Every cluster of the event has to reach the required number of hits, so
/// that no neutron is built on a group too small to place an interaction
/// point. An event without clusters never qualifies.
static bool ned2ClustersBigEnough(const std::vector<N2Cluster> &cl)
{
    if (cl.empty())
        return false;

    for (const N2Cluster &c : cl)
        if ((int)c.hits.size() < g2NMinHits)
            return false;

    return true;
}

static bool ned2NeutronBetaCutOn() { return g2NNeuBetaLo < g2NNeuBetaHi; }

/// Does one cluster look like a neutron? The velocity implied by its first hit
/// in time has to fall inside [g2NNeuBetaLo, g2NNeuBetaHi]: a gamma cluster
/// sits at beta0 ~ 1, noise anywhere. Minimum hits are handled apart.
static bool ned2NeutronLike(const N2Cluster &c)
{
    if (!ned2NeutronBetaCutOn())
        return true;

    causal::Cluster cc;
    cc.hits = c.hits;
    cc.firstT = c.firstT;
    cc.beta0 = ned2HitBeta(c.firstT);
    return causal::neutronLike(cc, g2NCausal, 1, g2NNeuBetaLo, g2NNeuBetaHi);
}

/// Every cluster of the event has to look like a neutron.
static bool ned2ClustersNeutronLike(const std::vector<N2Cluster> &cl)
{
    for (const N2Cluster &c : cl)
        if (!ned2NeutronLike(c))
            return false;
    return true;
}

// ─── Are the two clusters two neutrons? ─────────────────────────────────────

/// Velocity of the neutron that would have produced a hit, from the flight
/// path and the time of flight.
static double ned2HitBeta(int hit)
{
    if (hit < 0 || hit >= (int)g2NX->size())
        return -1.;

    const double x = (*g2NX)[hit], y = (*g2NY)[hit], z = (*g2NZ)[hit];
    const double t = (*g2NT)[hit];
    if (t <= 0.)
        return -1.;

    return std::sqrt(x * x + y * y + z * z) / (k2NC * t);
}

/// Apparent velocity between the earliest hits of the two clusters. A single
/// particle can only connect them when this stays below its own velocity.
static double ned2BetaRel(const std::vector<N2Cluster> &cl)
{
    if (cl.size() != 2)
        return -1.;

    const int a = cl[0].firstT, b = cl[1].firstT;
    const double dx = (*g2NX)[a] - (*g2NX)[b];
    const double dy = (*g2NY)[a] - (*g2NY)[b];
    const double dz = (*g2NZ)[a] - (*g2NZ)[b];
    const double dt = std::abs((*g2NT)[a] - (*g2NT)[b]);

    if (dt <= 0.) // same time, far apart: two particles for sure
        return 1e9;

    return std::sqrt(dx * dx + dy * dy + dz * dz) / (k2NC * dt);
}

/// Velocity the two clusters have to beat to count as two neutrons.
static double ned2BetaRelThreshold(const std::vector<N2Cluster> &cl)
{
    if (g2NBetaRelCut > 0.)
        return g2NBetaRelCut;

    // the neutron cannot scatter and then move faster than it arrived
    return std::max(ned2HitBeta(cl[0].firstT), ned2HitBeta(cl[1].firstT));
}

/// True when no single neutron could have made both clusters.
static bool ned2TwoNeutrons(const std::vector<N2Cluster> &cl)
{
    if (cl.size() != 2)
        return false;
    if (g2NBetaRelCut == 0.)
        return true;

    const double thr = ned2BetaRelThreshold(cl);
    if (thr <= 0.) // a hit with no usable time, the test cannot be made
        return false;

    return ned2BetaRel(cl) > thr;
}

// ─── Relative energy ────────────────────────────────────────────────────────

/// Read the offsets used by dataAnalysis.cpp: fragment x, fragment y,
/// neutron x, neutron y, beta matching (one value per line).
void loadOffsets2n(const char *offFile = k2NDefaultOffsets)
{
    const char *repo = getenv("repopath");
    const TString path = repo ? TString::Format("%s/final/settings/%s", repo, offFile)
                              : TString::Format("%s%s", k2NSettingsDir, offFile);

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

    g2NFragOff[0] = v[0];
    g2NFragOff[1] = v[1];
    g2NNeuOff[0] = v[2];
    g2NNeuOff[1] = v[3];
    g2NBetaMatch = v[4];
    g2NOffsetsFile = offFile;

    std::cout << "[OK] Erel offsets from " << path << " : frag ("
              << g2NFragOff[0] << ", " << g2NFragOff[1] << ")  neu ("
              << g2NNeuOff[0] << ", " << g2NNeuOff[1] << ")  betaMatch "
              << g2NBetaMatch << "\n";

    if (g2NEntry >= 0)
        ned2DrawEvent();
}

/// Four momentum of the fragment, with the transverse direction offsets and
/// the beta matching of the settings file.
static bool ned2FragmentP4(TLorentzVector &p4)
{
    const double betaF = g2NBetaFrag + g2NBetaMatch;
    if (betaF <= 0. || betaF >= 1. || g2NPzFrag == 0. || g2NMFrag <= 0.)
        return false;

    TVector3 dir(g2NPxFrag / g2NPzFrag - g2NFragOff[0],
                 g2NPyFrag / g2NPzFrag - g2NFragOff[1], 1.0);
    dir = dir.Unit();

    const double gammaF = 1.0 / std::sqrt(1.0 - betaF * betaF);
    p4.SetVectM(dir * (gammaF * g2NMFrag * betaF), g2NMFrag);
    return true;
}

/// Four momentum of a neutron from the position and the time of flight of one
/// hit, as the velocity is built in eventFilter.C.
static bool ned2NeutronP4(int hit, TLorentzVector &p4)
{
    if (hit < 0 || hit >= (int)g2NX->size())
        return false;

    const double x = (*g2NX)[hit], y = (*g2NY)[hit];
    const double z = (*g2NZ)[hit], t = (*g2NT)[hit];

    const double L = std::sqrt(x * x + y * y + z * z);
    if (L <= 0. || t <= 0. || z == 0.)
        return false;

    const double beta = L / (k2NC * t);
    if (beta <= 0. || beta >= 1.)
        return false;

    TVector3 dir(x / z - g2NNeuOff[0], y / z - g2NNeuOff[1], 1.0);
    dir = dir.Unit();

    const double gamma = 1.0 / std::sqrt(1.0 - beta * beta);
    p4.SetVectM(dir * (gamma * k2NMn * beta), k2NMn);
    return true;
}

/// Relative energy [GeV] of fragment + the neutrons given by `hits`, from the
/// invariant mass of the system. For one neutron this is the same expression
/// as DataAnalysis::getData().
static double ned2Erel(const std::vector<int> &hits)
{
    TLorentzVector sys;
    if (!ned2FragmentP4(sys))
        return -999.;

    for (int h : hits)
    {
        TLorentzVector n;
        if (!ned2NeutronP4(h, n))
            return -999.;
        sys += n;
    }

    return sys.M() - g2NMFrag - hits.size() * k2NMn;
}

/// Erel of a two cluster event taking the first hit in z of each cluster.
static double ned2ErelByPosition(const std::vector<N2Cluster> &cl)
{
    if (cl.size() != 2)
        return -999.;
    return ned2Erel({cl[0].firstZ, cl[1].firstZ});
}

/// Erel of a two cluster event taking the first hit in time of each cluster.
static double ned2ErelByTime(const std::vector<N2Cluster> &cl)
{
    if (cl.size() != 2)
        return -999.;
    return ned2Erel({cl[0].firstT, cl[1].firstT});
}

// ─── Selection ──────────────────────────────────────────────────────────────

/// Rebuild the list of entries passing the QFS window, the minimum number of
/// hits, the required number of clusters and the user cut.
static void ned2Rescan()
{
    if (!ned2Ready())
        return;

    g2NList.clear();

    TTreeFormula *form = nullptr;
    if (!g2NCut.IsNull())
    {
        form = new TTreeFormula("ned2Cut", g2NCut.Data(), g2NTree);
        if (form->GetNdim() == 0)
        {
            std::cerr << "[ERROR] Invalid cut: " << g2NCut << " (ignored)\n";
            delete form;
            form = nullptr;
            g2NCut = "";
        }
    }

    // denominators: events passing QFS and the user cut, with at least one hit,
    // split by cluster multiplicity 0, 1, 2, 3, > 3. The second row counts only
    // events whose clusters all pass the quality cuts (hits and beta0)
    Long64_t nPre = 0, nHit = 0;
    Long64_t nMult[5] = {0, 0, 0, 0, 0};
    Long64_t nMultGood[5] = {0, 0, 0, 0, 0};

    const Long64_t n = g2NTree->GetEntries();
    for (Long64_t i = 0; i < n; ++i)
    {
        g2NTree->GetEntry(i);
        if (g2NQFS && (g2NOpa <= k2NOpaMin || g2NOpa >= k2NOpaMax))
            continue;
        if (form && form->EvalInstance(0) == 0)
            continue;

        const std::vector<N2Cluster> cl = ned2FindClusters();

        ++nPre;
        if (!g2NX->empty())
            ++nHit;
        const int m = std::min((int)cl.size(), 4);
        ++nMult[m];
        if (!cl.empty() && ned2ClustersBigEnough(cl) && ned2ClustersNeutronLike(cl))
            ++nMultGood[m];

        if (g2NReqClusters > 0 && (int)cl.size() != g2NReqClusters)
            continue;
        if (!ned2ClustersBigEnough(cl))
            continue;
        if (!ned2ClustersNeutronLike(cl))
            continue;
        // a pair a single neutron could have made means the other one was
        // missed, so the event is not a 2n candidate
        if (g2NReqClusters == 2 && !ned2TwoNeutrons(cl))
            continue;

        g2NList.push_back(i);
    }

    delete form;

    if (g2NErelCanvas)
        erelSpectrum2n(g2NErelBins, g2NErelLo, g2NErelHi);

    std::cout << "[OK] " << g2NList.size() << " / " << n
              << " events pass the selection (>= " << g2NMinHits
              << " hits per cluster";
    if (g2NReqClusters > 0)
        std::cout << ", " << g2NReqClusters << " clusters";
    if (ned2NeutronBetaCutOn())
        std::cout << ", " << g2NNeuBetaLo << " < beta0 < " << g2NNeuBetaHi;
    if (g2NReqClusters == 2 && g2NBetaRelCut != 0.)
        std::cout << ", causality "
                  << (g2NBetaRelCut > 0.
                          ? Form("beta_rel > %.2f", g2NBetaRelCut)
                          : "beta_rel > neutron velocity");
    if (g2NUseCausal)
        std::cout << Form(", causal nSigma = %.2f, sigmaR = %.1f cm, sigmaT = %.3f ns,"
                          " gamma branch %s",
                          g2NCausal.nSigma, g2NCausal.sigmaR, g2NCausal.sigmaT,
                          g2NCausal.gammaBranch ? "on" : "off");
    else
        std::cout << ", sphere dR = " << g2NdR << " cm, dT = " << g2NdT << " ns";
    if (g2NQFS)
        std::cout << ", QFS " << k2NOpaMin << " < califa_opa < " << k2NOpaMax;
    if (!g2NCut.IsNull())
        std::cout << ", cut = \"" << g2NCut << "\"";
    std::cout << ")\n";

    std::cout << "     " << nPre << " events after QFS/cut, " << nHit
              << " with >= 1 hit. Cluster multiplicity:\n"
              << Form("     %-22s %8s %8s %8s %8s %8s\n",
                      "", "0", "1", "2", "3", ">3")
              << Form("     %-22s %8lld %8lld %8lld %8lld %8lld\n", "all",
                      nMult[0], nMult[1], nMult[2], nMult[3], nMult[4])
              << Form("     %-22s %8lld %8lld %8lld %8lld %8lld\n",
                      "all clusters pass qual.",
                      nMultGood[0], nMultGood[1], nMultGood[2], nMultGood[3],
                      nMultGood[4]);

    g2NPos = 0;
}

// ─── Canvas construction ────────────────────────────────────────────────────

/// (Re)create the axis frames from the current g2NLo / g2NHi ranges.
static void ned2MakeFrames()
{
    for (int p = 0; p < 6; ++p)
    {
        double uLo, uHi, vLo, vHi;
        ned2PadRange(p, uLo, uHi, vLo, vHi);

        const TString title = TString::Format("%s;%s;%s",
                                              k2NPanelTitle[p],
                                              k2NVarName[k2NU[p]],
                                              k2NVarName[k2NV[p]]);

        delete g2NFrame[p];
        g2NFrame[p] = new TH2F(Form("ned2Frame%d", p), title,
                               1, uLo, uHi, 1, vLo, vHi);

        for (TAxis *a : {g2NFrame[p]->GetXaxis(), g2NFrame[p]->GetYaxis()})
        {
            a->SetTitleSize(0.05);
            a->SetLabelSize(0.04);
        }
        g2NFrame[p]->GetXaxis()->SetTitleOffset(1.05);
        g2NFrame[p]->GetYaxis()->SetTitleOffset(1.35);
    }

    delete g2NFrame3D;
    g2NFrame3D = new TH3F("ned2Frame3D", "3D view;z [cm];x [cm];y [cm]",
                          1, g2NLo[2], g2NHi[2],
                          1, g2NLo[0], g2NHi[0],
                          1, g2NLo[1], g2NHi[1]);
    for (TAxis *a : {g2NFrame3D->GetXaxis(), g2NFrame3D->GetYaxis(),
                     g2NFrame3D->GetZaxis()})
    {
        a->SetTitleSize(0.05);
        a->SetLabelSize(0.035);
        a->SetTitleOffset(1.8);
    }
}

static void ned2BuildCanvas()
{
    gStyle->SetOptStat(0);
    gStyle->SetPadTickX(1);
    gStyle->SetPadTickY(1);
    gStyle->SetTitleFontSize(0.06);

    ned2MakeFrames();

    if (g2NCanvas)
        return;

    g2NCanvas = new TCanvas("cNeutron2nEventDisplay",
                            "NeuLAND 2n event display", 1500, 800);
    g2NCanvas->Divide(4, 2, 0.001, 0.001);

    for (int p = 1; p <= 8; ++p)
    {
        TVirtualPad *pad = g2NCanvas->cd(p);
        pad->SetLeftMargin(0.15);
        pad->SetRightMargin(0.04);
        pad->SetBottomMargin(0.13);
        pad->SetTopMargin(0.10);
    }

    if (TPad *pad3d = dynamic_cast<TPad *>(g2NCanvas->cd(8)))
    {
        pad3d->SetTheta(25.);
        pad3d->SetPhi(-50.);
    }
}

// ─── Drawing ────────────────────────────────────────────────────────────────

static void ned2DrawInfo(const std::vector<N2Cluster> &cl)
{
    TPaveText *pt = new TPaveText(0.02, 0.02, 0.98, 0.97, "NDC");
    pt->SetBit(kCanDelete);
    pt->SetFillColor(0);
    pt->SetBorderSize(0);
    pt->SetTextAlign(12);
    pt->SetTextFont(42);
    pt->SetTextSize(0.048);

    if (g2NInList)
        pt->AddText(Form("#bf{entry %lld}   (%d / %d selected)",
                         g2NEntry, g2NPos + 1, (int)g2NList.size()));
    else
        pt->AddText(Form("#bf{entry %lld}   (not in the current selection)",
                         g2NEntry));
    pt->AddText(Form("hits: %d      clusters: %d",
                     (int)g2NX->size(), (int)cl.size()));
    pt->AddText("");

    const int nShow = std::min((int)cl.size(), 4);
    for (int i = 0; i < nShow; ++i)
        pt->AddText(Form("#color[%d]{#bullet} cluster %d: %d hits, "
                         "z_{1} = %.0f, t_{1} = %.1f ns",
                         ned2ClusterColor(i), i, (int)cl[i].hits.size(),
                         cl[i].zMin, cl[i].tMin));
    if ((int)cl.size() > nShow)
        pt->AddText(Form("... and %d more clusters", (int)cl.size() - nShow));

    pt->AddText("");

    if (cl.size() == 2)
    {
        const double bRel = ned2BetaRel(cl);
        const double bThr = ned2BetaRelThreshold(cl);
        const bool two = ned2TwoNeutrons(cl);

        pt->AddText(Form("#beta_{rel} = %s   (%s)",
                         bRel > 1e8 ? "#infty" : Form("%.2f", bRel),
                         g2NBetaRelCut == 0.
                             ? "test off"
                             : TString::Format("threshold %.2f", bThr).Data()));
        pt->AddText(two ? "#color[418]{two neutrons}"
                        : "#color[633]{one neutron scattered, partner missed}");

        const double ePos = ned2ErelByPosition(cl);
        const double eTime = ned2ErelByTime(cl);
        pt->AddText(ePos > -900. ? Form("#bf{E_{rel} (first in z) = %.3f MeV}",
                                        ePos * 1000.)
                                 : "#bf{E_{rel} (first in z) = n/a}");
        pt->AddText(eTime > -900. ? Form("#bf{E_{rel} (first in t) = %.3f MeV}",
                                         eTime * 1000.)
                                  : "#bf{E_{rel} (first in t) = n/a}");
        if (cl[0].firstZ == cl[0].firstT && cl[1].firstZ == cl[1].firstT)
            pt->AddText("#it{both criteria pick the same hits}");
    }
    else
    {
        pt->AddText("#bf{E_{rel}: needs exactly two clusters}");
    }

    pt->AddText("");
    pt->AddText("#bf{Fragment}");
    pt->AddText(Form("A/Q = %.3f   Z = %.2f   #beta = %.4f",
                     g2NAoQ, g2NZfrag, g2NBetaFrag));
    pt->AddText(Form("CALIFA opa = %.3f rad (%.1f deg)",
                     g2NOpa, g2NOpa * TMath::RadToDeg()));
    if (g2NQFS)
        pt->AddText(Form("#it{QFS window %.2f - %.2f rad}", k2NOpaMin, k2NOpaMax));
    if (g2NUseCausal)
        pt->AddText(Form("#it{causal clusters: %.1f#sigma, #sigma_{R} = %.1f cm, "
                         "#sigma_{t} = %.2f ns, #geq %d hits}",
                         g2NCausal.nSigma, g2NCausal.sigmaR, g2NCausal.sigmaT,
                         g2NMinHits));
    else
        pt->AddText(Form("#it{clusters: dR < %.0f cm, dt < %.0f ns, #geq %d hits}",
                         g2NdR, g2NdT, g2NMinHits));
    pt->AddText("#it{star: first in z, square: first in time}");
    pt->Draw();
}

static void ned2Draw3D(const std::vector<N2Cluster> &cl)
{
    g2NFrame3D->Draw();

    for (size_t i = 0; i < cl.size(); ++i)
    {
        TPolyMarker3D *pm = new TPolyMarker3D((int)cl[i].hits.size());
        pm->SetBit(kCanDelete);
        for (size_t h = 0; h < cl[i].hits.size(); ++h)
        {
            const int j = cl[i].hits[h];
            pm->SetPoint((int)h, (*g2NZ)[j], (*g2NX)[j], (*g2NY)[j]);
        }
        pm->SetMarkerStyle(20);
        pm->SetMarkerSize(1.1);
        pm->SetMarkerColor(ned2ClusterColor((int)i));
        pm->Draw();

        TPolyMarker3D *first = new TPolyMarker3D(1);
        first->SetBit(kCanDelete);
        const int j = cl[i].firstZ;
        first->SetPoint(0, (*g2NZ)[j], (*g2NX)[j], (*g2NY)[j]);
        first->SetMarkerStyle(29);
        first->SetMarkerSize(2.2);
        first->SetMarkerColor(ned2ClusterColor((int)i));
        first->Draw();
    }
}

/// Draw the currently loaded entry.
static void ned2DrawEvent()
{
    if (!g2NCanvas)
        return;

    const std::vector<N2Cluster> cl = ned2FindClusters();

    for (int p = 0; p < 6; ++p)
    {
        TVirtualPad *pad = g2NCanvas->cd(p < 3 ? p + 1 : p + 2);
        pad->Clear();
        pad->SetGrid();

        g2NFrame[p]->Draw("AXIS");

        for (size_t i = 0; i < cl.size(); ++i)
        {
            const Int_t col = ned2ClusterColor((int)i);

            // every hit of the cluster
            for (int j : cl[i].hits)
            {
                TMarker *m = new TMarker(
                    ned2Var(k2NU[p], (*g2NX)[j], (*g2NY)[j], (*g2NZ)[j], (*g2NT)[j]),
                    ned2Var(k2NV[p], (*g2NX)[j], (*g2NY)[j], (*g2NZ)[j], (*g2NT)[j]), 20);
                m->SetMarkerColor(col);
                m->SetMarkerSize(1.4);
                m->SetBit(kCanDelete);
                m->Draw();
            }

            // first hit in z: filled star
            const int jz = cl[i].firstZ;
            TMarker *mz = new TMarker(
                ned2Var(k2NU[p], (*g2NX)[jz], (*g2NY)[jz], (*g2NZ)[jz], (*g2NT)[jz]),
                ned2Var(k2NV[p], (*g2NX)[jz], (*g2NY)[jz], (*g2NZ)[jz], (*g2NT)[jz]), 29);
            mz->SetMarkerColor(col);
            mz->SetMarkerSize(2.6);
            mz->SetBit(kCanDelete);
            mz->Draw();

            // first hit in time: open square on top
            const int jt = cl[i].firstT;
            TMarker *mt = new TMarker(
                ned2Var(k2NU[p], (*g2NX)[jt], (*g2NY)[jt], (*g2NZ)[jt], (*g2NT)[jt]),
                ned2Var(k2NV[p], (*g2NX)[jt], (*g2NY)[jt], (*g2NZ)[jt], (*g2NT)[jt]), 25);
            mt->SetMarkerColor(kBlack);
            mt->SetMarkerSize(2.2);
            mt->SetBit(kCanDelete);
            mt->Draw();

            // line from the target to the first hit in z
            if (p < 3)
            {
                TLine *tr = new TLine(
                    ned2Var(k2NU[p], 0, 0, 0, 0),
                    ned2Var(k2NV[p], 0, 0, 0, 0),
                    ned2Var(k2NU[p], (*g2NX)[jz], (*g2NY)[jz], (*g2NZ)[jz], (*g2NT)[jz]),
                    ned2Var(k2NV[p], (*g2NX)[jz], (*g2NY)[jz], (*g2NZ)[jz], (*g2NT)[jz]));
                tr->SetLineColor(col);
                tr->SetLineStyle(2);
                tr->SetBit(kCanDelete);
                tr->Draw();
            }
        }

        pad->Modified();
    }

    TVirtualPad *info = g2NCanvas->cd(4);
    info->Clear();
    ned2DrawInfo(cl);
    info->Modified();

    TVirtualPad *pad3d = g2NCanvas->cd(8);
    pad3d->Clear();
    ned2Draw3D(cl);
    pad3d->Modified();

    g2NCanvas->Update();

    ned2SetStatus(g2NInList
                      ? TString::Format("entry %lld   %d / %d selected   "
                                        "hits: %d   clusters: %d",
                                        g2NEntry, g2NPos + 1, (int)g2NList.size(),
                                        (int)g2NX->size(), (int)cl.size())
                      : TString::Format("entry %lld   not in the selection   "
                                        "hits: %d   clusters: %d",
                                        g2NEntry, (int)g2NX->size(), (int)cl.size()));
}

// ─── Erel spectrum canvas ───────────────────────────────────────────────────

/// Refill the two Erel histograms from the two cluster events of the selection.
static void ned2FillErelSpectrum()
{
    if (!ned2Ready() || !g2NErelPos)
        return;

    g2NErelPos->Reset();
    g2NErelTime->Reset();
    g2NErelSame->Reset();
    g2NErelMin->Reset();

    const Long64_t keep = g2NEntry;
    int nTwo = 0, nSameHits = 0;

    for (size_t i = 0; i < g2NList.size(); ++i)
    {
        g2NTree->GetEntry(g2NList[i]);

        const std::vector<N2Cluster> cl = ned2FindClusters();
        if (cl.size() != 2 || !ned2TwoNeutrons(cl))
            continue;
        ++nTwo;

        const bool sameHits = (cl[0].firstZ == cl[0].firstT &&
                               cl[1].firstZ == cl[1].firstT);
        if (sameHits)
            ++nSameHits;

        const double ePos = ned2ErelByPosition(cl);
        if (ePos > -900.)
        {
            g2NErelPos->Fill(ePos * 1000.);
            // both criteria pick the same hits, so the two values coincide
            if (sameHits)
                g2NErelSame->Fill(ePos * 1000.);
        }

        const double eTime = ned2ErelByTime(cl);
        if (eTime > -900.)
            g2NErelTime->Fill(eTime * 1000.);

        // take, event by event, whichever hit gives the smaller energy
        if (ePos > -900. && eTime > -900.)
            g2NErelMin->Fill(std::min(ePos, eTime) * 1000.);
    }

    std::cout << "[OK] Erel over " << nTwo << " 2n candidates: "
              << "first in z <Erel> = " << g2NErelPos->GetMean() << " MeV, "
              << "first in time <Erel> = " << g2NErelTime->GetMean() << " MeV, "
              << "both criteria agree in " << nSameHits << " events ("
              << (nTwo ? 100. * nSameHits / nTwo : 0.) << " %)\n";

    if (keep >= 0)
        g2NTree->GetEntry(keep);
}

/// Second canvas: Erel of the two cluster events, taking the first hit in z
/// and the first hit in time of each cluster.
void erelSpectrum2n(int nbins, double lo, double hi)
{
    if (!ned2Ready())
        return;

    g2NErelBins = nbins;
    g2NErelLo = lo;
    g2NErelHi = hi;

    delete g2NErelPos;
    delete g2NErelTime;
    delete g2NErelSame;
    delete g2NErelMin;

    g2NErelPos = new TH1F("ned2ErelPos",
                          "E_{rel} (2n), first hit in z of each cluster"
                          ";E_{rel} [MeV];counts",
                          nbins, lo, hi);
    g2NErelTime = new TH1F("ned2ErelTime",
                           "E_{rel} (2n), first hit in time of each cluster"
                           ";E_{rel} [MeV];counts",
                           nbins, lo, hi);
    g2NErelSame = new TH1F("ned2ErelSame",
                           "E_{rel} (2n), same hit in z and time in both clusters"
                           ";E_{rel} [MeV];counts",
                           nbins, lo, hi);
    g2NErelMin = new TH1F("ned2ErelMin",
                          "E_{rel} (2n), whichever hit gives less energy"
                          ";E_{rel} [MeV];counts",
                          nbins, lo, hi);

    g2NErelPos->SetLineColor(kBlack);
    g2NErelTime->SetLineColor(kRed + 1);
    g2NErelSame->SetLineColor(kAzure + 2);
    g2NErelMin->SetLineColor(kGreen + 2);
    for (TH1F *h : {g2NErelPos, g2NErelTime, g2NErelSame, g2NErelMin})
    {
        h->SetLineWidth(2);
        h->SetStats(0);
    }

    ned2FillErelSpectrum();

    if (!g2NErelCanvas)
        g2NErelCanvas = new TCanvas("cNeutron2nErel",
                                    "Erel of the 2n candidates", 1100, 800);

    g2NErelCanvas->cd();
    g2NErelCanvas->Clear();
    g2NErelCanvas->Divide(2, 2, 0.001, 0.001);

    TH1F *hist[4] = {g2NErelPos, g2NErelTime, g2NErelSame, g2NErelMin};

    const TString clustering =
        g2NUseCausal
            ? TString::Format("causal %.1f#sigma (#sigma_{R} %.1f cm, #sigma_{t} %.2f ns, "
                              "#gamma %s%s)",
                              g2NCausal.nSigma, g2NCausal.sigmaR, g2NCausal.sigmaT,
                              g2NCausal.gammaBranch ? "on" : "off",
                              g2NCausal.strictToFirst ? ", strict" : "")
            : TString::Format("dR < %.0f cm, dt < %.0f ns", g2NdR, g2NdT);
    const TString betaWin =
        ned2NeutronBetaCutOn()
            ? TString::Format("%.2f < #beta_{0} < %.2f, ", g2NNeuBetaLo, g2NNeuBetaHi)
            : TString("");

    const TString selection =
        Form("#it{%s, %s, #geq %d hits/cluster, %s%s%s}",
             g2NQFS ? Form("QFS %.2f - %.2f rad", k2NOpaMin, k2NOpaMax)
                    : "no QFS window",
             clustering.Data(), g2NMinHits, betaWin.Data(),
             g2NBetaRelCut == 0.
                 ? "no causality cut"
                 : (g2NBetaRelCut > 0.
                        ? Form("#beta_{rel} > %.2f", g2NBetaRelCut)
                        : "#beta_{rel} > #beta_{n}"),
             g2NCut.IsNull() ? "" : Form(", %s", g2NCut.Data()));

    for (int i = 0; i < 4; ++i)
    {
        TVirtualPad *pad = g2NErelCanvas->cd(i + 1);
        pad->SetLeftMargin(0.13);
        pad->SetRightMargin(0.04);
        pad->SetBottomMargin(0.12);
        pad->SetTopMargin(0.10);
        pad->SetTicks(1, 1);

        hist[i]->SetMaximum(1.30 * hist[i]->GetMaximum());
        hist[i]->Draw("HIST");

        TLatex *num = new TLatex(0.17, 0.84,
                                 Form("%d events, mean %.2f MeV",
                                      (int)hist[i]->GetEntries(),
                                      hist[i]->GetMean()));
        num->SetBit(kCanDelete);
        num->SetNDC();
        num->SetTextSize(0.042);
        num->Draw();

        if (i == 0)
        {
            TLatex *sel = new TLatex(0.17, 0.78, selection);
            sel->SetBit(kCanDelete);
            sel->SetNDC();
            sel->SetTextSize(0.034);
            sel->Draw();
        }

        pad->Modified();
    }

    g2NErelCanvas->Update();
}

// ─── Navigation (also used as GUI slots) ────────────────────────────────────

/// Show a given tree entry, whether or not it passes the current selection.
void showEvent2n(Long64_t entry)
{
    if (!ned2Ready())
        return;
    if (entry < 0 || entry >= g2NTree->GetEntries())
    {
        std::cerr << "[ERROR] Entry " << entry << " out of range [0, "
                  << g2NTree->GetEntries() - 1 << "]\n";
        return;
    }

    g2NTree->GetEntry(entry);
    g2NEntry = entry;

    auto it = std::find(g2NList.begin(), g2NList.end(), entry);
    g2NInList = (it != g2NList.end());
    if (g2NInList)
        g2NPos = (int)(it - g2NList.begin());

    ned2DrawEvent();
}

/// Show the i-th event of the current selection.
static void ned2ShowIndex(int i)
{
    if (!ned2Ready())
        return;
    if (g2NList.empty())
    {
        std::cerr << "[ERROR] No event passes the current selection.\n";
        return;
    }
    g2NPos = ((i % (int)g2NList.size()) + (int)g2NList.size()) % (int)g2NList.size();
    showEvent2n(g2NList[g2NPos]);
}

void nextEvent2n() { ned2ShowIndex(g2NPos + 1); }
void prevEvent2n() { ned2ShowIndex(g2NPos - 1); }

void randomEvent2n()
{
    if (g2NList.empty())
        return;
    static TRandom3 rng(0);
    ned2ShowIndex((int)(rng.Uniform(0, g2NList.size())));
}

void printEvent2n()
{
    if (!ned2Ready() || g2NEntry < 0)
        return;

    const std::vector<N2Cluster> cl = ned2FindClusters();

    std::cout << "\n─── entry " << g2NEntry << " ─ " << g2NX->size()
              << " NeuLAND hits ─ " << cl.size() << " clusters ───\n"
              << std::setw(6) << "hit" << std::setw(9) << "cluster"
              << std::setw(9) << "paddle"
              << std::setw(11) << "x [cm]" << std::setw(11) << "y [cm]"
              << std::setw(12) << "z [cm]" << std::setw(11) << "t [ns]"
              << "   flag\n";

    std::vector<int> label(g2NX->size(), -1);
    for (size_t i = 0; i < cl.size(); ++i)
        for (int j : cl[i].hits)
            label[j] = (int)i;

    for (size_t h = 0; h < g2NX->size(); ++h)
    {
        TString flag;
        for (size_t i = 0; i < cl.size(); ++i)
        {
            if (cl[i].firstZ == (int)h)
                flag += " firstZ";
            if (cl[i].firstT == (int)h)
                flag += " firstT";
        }

        std::cout << std::fixed << std::setprecision(2)
                  << std::setw(6) << h
                  << std::setw(9) << label[h]
                  << std::setw(9) << (*g2NPaddle)[h]
                  << std::setw(11) << (*g2NX)[h]
                  << std::setw(11) << (*g2NY)[h]
                  << std::setw(12) << (*g2NZ)[h]
                  << std::setw(11) << (*g2NT)[h]
                  << "  " << flag << "\n";
    }

    if (cl.size() == 2)
        std::cout << "beta_rel = " << ned2BetaRel(cl) << " (threshold "
                  << ned2BetaRelThreshold(cl) << "), "
                  << (ned2TwoNeutrons(cl)
                          ? "two neutrons"
                          : "one neutron scattered, partner missed")
                  << "\n";

    if (cl.size() == 2)
        std::cout << "Erel (first in z) = " << ned2ErelByPosition(cl) * 1000.
                  << " MeV, Erel (first in time) = "
                  << ned2ErelByTime(cl) * 1000. << " MeV\n\n";
    else
        std::cout << "Erel needs exactly two clusters, this event has "
                  << cl.size() << "\n\n";
}

void saveEvent2n(const char *fileName)
{
    if (!g2NCanvas)
        return;
    const TString name = fileName ? TString(fileName)
                                  : TString::Format("neutron2n_event_%lld.png", g2NEntry);
    g2NCanvas->SaveAs(name);
}

// ─── Selection control ──────────────────────────────────────────────────────

/// Distance [cm] and time [ns] that join two hits into the same cluster.
void setClusterCuts(double dR, double dT)
{
    if (dR <= 0. || dT <= 0.)
    {
        std::cerr << "[ERROR] setClusterCuts(dR, dT) needs positive values\n";
        return;
    }
    g2NdR = dR;
    g2NdT = dT;
    ned2Rescan();
    ned2ShowIndex(0);
}

/// Threshold on the apparent velocity between the two clusters. A negative
/// value compares it with the neutron velocity of the event, so that pairs a
/// single neutron could have made are rejected. Zero, the default, switches
/// the test off (the causal clustering subsumes it), a positive value is used
/// as a fixed threshold.
void setCausalityCut(double betaMin = -1.)
{
    g2NBetaRelCut = betaMin;
    ned2Rescan();
    ned2ShowIndex(0);
}

/// Causal clustering (true, default) or the old fixed sphere (false).
void setUseCausal(bool on)
{
    g2NUseCausal = on;
    ned2Rescan();
    ned2ShowIndex(0);
}

/// Width of the causality window in standard deviations (g2NCausal.nSigma).
void setCausalSigma(double n)
{
    if (n <= 0.)
    {
        std::cerr << "[ERROR] setCausalSigma(n) needs n > 0\n";
        return;
    }
    g2NCausal.nSigma = n;
    ned2Rescan();
    ned2ShowIndex(0);
}

/// One-hit position [cm] and time [ns] resolutions of the causal clustering.
void setResolutions(double sigmaR, double sigmaT)
{
    if (sigmaR <= 0. || sigmaT <= 0.)
    {
        std::cerr << "[ERROR] setResolutions(sigmaR, sigmaT) needs positive values\n";
        return;
    }
    g2NCausal.sigmaR = sigmaR;
    g2NCausal.sigmaT = sigmaT;
    ned2Rescan();
    ned2ShowIndex(0);
}

/// Follow prompt gammas emitted at an interaction point into the same cluster.
void setGammaBranch(bool on)
{
    g2NCausal.gammaBranch = on;
    ned2Rescan();
    ned2ShowIndex(0);
}

/// Beam velocity window of the gamma branch guard: a later hit whose apparent
/// velocity from the target lies inside it is never taken as a gamma.
/// lo >= hi switches the guard off.
void setBeamBetaWindow(double lo, double hi)
{
    g2NCausal.beamBetaLo = lo;
    g2NCausal.beamBetaHi = hi;
    ned2Rescan();
    ned2ShowIndex(0);
}

/// Link only against the earliest hit of a cluster (merging systematic).
void setStrictToFirst(bool on)
{
    g2NCausal.strictToFirst = on;
    ned2Rescan();
    ned2ShowIndex(0);
}

/// Window on beta0, the velocity implied by the first hit in time of every
/// cluster. lo >= hi switches the cut off.
void setNeutronBetaWindow(double lo, double hi)
{
    g2NNeuBetaLo = lo;
    g2NNeuBetaHi = hi;
    ned2Rescan();
    ned2ShowIndex(0);
}

/// n-n relative energy below which the causal clustering merges a genuine
/// pair, for a flight path L [cm] and a neutron velocity beta, compared with
/// the current sphere.
void causalAcceptance(double L = 1300., double beta = 0.75)
{
    causal::printAcceptance(g2NCausal, L, beta, g2NdR);
}

/// Number of clusters an event must have to enter the selection, 0 for any.
void setRequiredClusters(int n)
{
    g2NReqClusters = std::max(0, n);
    ned2Rescan();
    ned2ShowIndex(0);
}

/// Minimum number of hits every cluster of the event must have.
void setMinHits2n(int n)
{
    g2NMinHits = std::max(1, n);
    ned2Rescan();
    ned2ShowIndex(0);
}

void applyCut2n(const char *cut)
{
    g2NCut = (cut ? cut : "");
    ned2Rescan();
    ned2ShowIndex(0);
}

/// Keep only quasi-free scattering events (1.25 < califa_opa < 1.65 rad).
void setQFS2n(bool on)
{
    g2NQFS = on;
    ned2Rescan();
    ned2ShowIndex(0);
}

/// Manual axis range for one coordinate: var is "x", "y", "z" or "t".
void setRange2n(const char *var, double lo, double hi)
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
        std::cerr << "[ERROR] setRange2n(\"x\"|\"y\"|\"z\"|\"t\", lo, hi) with lo < hi\n";
        return;
    }

    g2NLo[k] = lo;
    g2NHi[k] = hi;
    ned2MakeFrames();
    if (g2NEntry >= 0)
        ned2DrawEvent();
}

void playEvents2n(int nEvents = 20, double seconds = 1.0)
{
    for (int i = 0; i < nEvents; ++i)
    {
        nextEvent2n();
        gSystem->ProcessEvents();
        gSystem->Sleep((ULong_t)(seconds * 1000));
    }
}

// ─── GUI ────────────────────────────────────────────────────────────────────

void ned2GuiGoto()
{
    if (g2NGotoBox)
        showEvent2n((Long64_t)g2NGotoBox->GetNumberEntry()->GetIntNumber());
}

void ned2GuiSave() { saveEvent2n(nullptr); }

void ned2CloseGui()
{
    if (g2NGui)
    {
        g2NGui->CloseWindow();
        g2NGui = nullptr;
        g2NGotoBox = nullptr;
        g2NStatus = nullptr;
    }
}

/// Keyboard shortcuts while the mouse is over the canvas.
void ned2CanvasEvent(Int_t event, Int_t px, Int_t /*py*/, TObject * /*sel*/)
{
    if (event != kKeyPress)
        return;

    switch (px)
    {
    case 'n':
        nextEvent2n();
        break;
    case 'p':
        prevEvent2n();
        break;
    case 'r':
        randomEvent2n();
        break;
    case 'd':
        printEvent2n();
        break;
    case 's':
        saveEvent2n(nullptr);
        break;
    case 'h':
        ned2nHelp();
        break;
    default:
        break;
    }
}

static void ned2BuildGui()
{
    if (gROOT->IsBatch() || g2NGui)
        return;

    g2NGui = new TGMainFrame(gClient->GetRoot(), 520, 90);
    g2NGui->SetWindowName("Neutron 2n event display");

    TGHorizontalFrame *row = new TGHorizontalFrame(g2NGui, 520, 40);

    TGTextButton *bPrev = new TGTextButton(row, "  << Prev  ");
    TGTextButton *bNext = new TGTextButton(row, "  Next >>  ");
    TGTextButton *bRand = new TGTextButton(row, " Random ");
    TGTextButton *bPrint = new TGTextButton(row, " Print ");
    TGTextButton *bSave = new TGTextButton(row, " Save ");

    bPrev->Connect("Clicked()", 0, 0, "prevEvent2n()");
    bNext->Connect("Clicked()", 0, 0, "nextEvent2n()");
    bRand->Connect("Clicked()", 0, 0, "randomEvent2n()");
    bPrint->Connect("Clicked()", 0, 0, "printEvent2n()");
    bSave->Connect("Clicked()", 0, 0, "ned2GuiSave()");

    for (TGTextButton *b : {bPrev, bNext, bRand, bPrint, bSave})
        row->AddFrame(b, new TGLayoutHints(kLHintsLeft | kLHintsCenterY, 3, 3, 3, 3));

    row->AddFrame(new TGLabel(row, "entry:"),
                  new TGLayoutHints(kLHintsLeft | kLHintsCenterY, 8, 3, 3, 3));
    g2NGotoBox = new TGNumberEntry(row, 0, 8, -1,
                                   TGNumberFormat::kNESInteger,
                                   TGNumberFormat::kNEANonNegative);
    row->AddFrame(g2NGotoBox,
                  new TGLayoutHints(kLHintsLeft | kLHintsCenterY, 3, 3, 3, 3));

    TGTextButton *bGo = new TGTextButton(row, " Go ");
    bGo->Connect("Clicked()", 0, 0, "ned2GuiGoto()");
    row->AddFrame(bGo, new TGLayoutHints(kLHintsLeft | kLHintsCenterY, 3, 3, 3, 3));

    g2NGui->AddFrame(row, new TGLayoutHints(kLHintsExpandX));

    g2NStatus = new TGLabel(g2NGui, "                                        ");
    g2NGui->AddFrame(g2NStatus,
                     new TGLayoutHints(kLHintsLeft | kLHintsExpandX, 6, 6, 2, 4));

    g2NGui->Connect("CloseWindow()", 0, 0, "ned2CloseGui()");
    g2NGui->MapSubwindows();
    g2NGui->Resize(g2NGui->GetDefaultSize());
    g2NGui->MapWindow();
}

// ─── Help ───────────────────────────────────────────────────────────────────

void ned2nHelp()
{
    std::cout << "\n─── 2n neutron event display ────────────────────────────\n"
              << "  nextEvent2n()            next selected event\n"
              << "  prevEvent2n()            previous selected event\n"
              << "  randomEvent2n()          random selected event\n"
              << "  showEvent2n(entry)       jump to a given tree entry\n"
              << "  printEvent2n()           hit table with cluster labels\n"
              << "  saveEvent2n(\"name.png\")  save the display canvas\n"
              << "  setUseCausal(true)       causal clustering, false = sphere\n"
              << "  setCausalSigma(2)        causality window [sigma]\n"
              << "  setResolutions(4, 0.15)  hit sigma_R [cm], sigma_t [ns]\n"
              << "  setGammaBranch(true)     join prompt gammas to the cluster\n"
              << "  setBeamBetaWindow(.7,.8) gamma branch guard, lo>=hi off\n"
              << "  setStrictToFirst(false)  link only to the first hit\n"
              << "  causalAcceptance(1300, 0.75)  merged Erel(nn) vs sphere\n"
              << "  setNeutronBetaWindow(.55,.9) beta0 of every cluster,\n"
              << "                           lo >= hi switches it off\n"
              << "  setClusterCuts(50, 5)    sphere dR [cm] and dt [ns]\n"
              << "  setRequiredClusters(2)   cluster multiplicity, 0 for any\n"
              << "  setCausalityCut(0)       beta_rel threshold: <0 neutron\n"
              << "                           velocity, 0 off, >0 fixed value\n"
              << "  old sphere selection:    setUseCausal(false);\n"
              << "                           setCausalityCut(-1);\n"
              << "                           setNeutronBetaWindow(0, 0)\n"
              << "  setMinHits2n(n)          every cluster needs >= n hits\n"
              << "  applyCut2n(\"beta_neu>0.6\") TTree cut on scalar branches\n"
              << "  setQFS2n(false)          drop the 1.25 < opa < 1.65 window\n"
              << "  setRange2n(\"t\", 50, 80)  axis range of x, y, z or t\n"
              << "  loadOffsets2n(\"23O1n.txt\") offsets used for E_rel\n"
              << "  erelSpectrum2n(100,-5,20)  redraw the E_rel canvas\n"
              << "  playEvents2n(20, 0.5)    slideshow: n events, seconds each\n"
              << "  canvas keys: n/p/r/d/s/h\n"
              << "─────────────────────────────────────────────────────────\n\n";
}

// ═══════════════════════════════════════════════════════════════════════════
//  MAIN
// ═══════════════════════════════════════════════════════════════════════════

void neutronEventDisplay2n(const char *fileName = k2NDefaultFile,
                           int minHits = 1,
                           const char *cut = "",
                           bool qfsOnly = true)
{
    TH1::AddDirectory(kFALSE);
    ned2Unload();

    if (g2NOffsetsFile.IsNull())
        loadOffsets2n();

    g2NFile = TFile::Open(fileName, "READ");
    if (!g2NFile || g2NFile->IsZombie())
    {
        std::cerr << "[ERROR] Cannot open " << fileName << "\n";
        return;
    }

    g2NTree = dynamic_cast<TTree *>(g2NFile->Get(k2NTreeName));
    if (!g2NTree)
    {
        std::cerr << "[ERROR] Tree " << k2NTreeName << " not found in "
                  << fileName << "\n";
        ned2Unload();
        return;
    }
    if (!g2NTree->GetBranch("n_x"))
    {
        std::cerr << "[ERROR] No NeuLAND hit branches (n_x, n_y, n_z, n_t) in "
                  << fileName << ". Was it produced with hasNeutrons = true?\n";
        ned2Unload();
        return;
    }

    std::cout << "[OK] " << fileName << " : " << g2NTree->GetEntries()
              << " entries\n";

    ned2SetBranches();
    ned2ComputeRanges();
    ned2BuildCanvas();

    g2NMinHits = std::max(1, minHits);
    g2NCut = cut ? cut : "";
    g2NQFS = qfsOnly;
    ned2Rescan();

    ned2BuildGui();
    if (!gROOT->IsBatch())
        g2NCanvas->Connect("ProcessedEvent(Int_t,Int_t,Int_t,TObject*)", 0, 0,
                           "ned2CanvasEvent(Int_t,Int_t,Int_t,TObject*)");

    ned2ShowIndex(0);
    erelSpectrum2n();
    ned2nHelp();
}

void erel1nCheck(double sn = 4.2, int nbins = 100, double lo = 0., double hi = 20.)
{
    if (!ned2Ready())
        return;
    TH1F *h = new TH1F("ned2Check1n", "E_{exc} (1n);E_{exc} [MeV];counts", nbins, lo, hi);
    const Long64_t keep = g2NEntry;
    for (Long64_t i = 0; i < g2NTree->GetEntries(); ++i)
    {
        g2NTree->GetEntry(i);
        if (g2NQFS && (g2NOpa <= k2NOpaMin || g2NOpa >= k2NOpaMax))
            continue;
        TLorentzVector sys, n;
        if (!ned2FragmentP4(sys) || g2NTf <= 0. || g2NZf == 0.)
            continue;
        const double L = std::sqrt(g2NXf * g2NXf + g2NYf * g2NYf + g2NZf * g2NZf);
        const double beta = L / (k2NC * g2NTf);
        if (beta <= 0. || beta >= 1.)
            continue;
        TVector3 dir(g2NXf / g2NZf - g2NNeuOff[0], g2NYf / g2NZf - g2NNeuOff[1], 1.);
        n.SetVectM(dir.Unit() * (k2NMn * beta / std::sqrt(1. - beta * beta)), k2NMn);
        sys += n;
        h->Fill((sys.M() - g2NMFrag - k2NMn) * 1000. + sn);
    }
    if (keep >= 0)
        g2NTree->GetEntry(keep);
    new TCanvas("cCheck1n", "check 1n", 800, 600);
    h->Draw("HIST");
}
