#include <array>
#include <string>
#include <vector>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <algorithm>
#include <functional>
#include <cmath>

#include "TFile.h"
#include "TSystem.h"
#include "TString.h"
#include "TVector3.h"
#include "TVector2.h"
#include "TLorentzVector.h"
#include "TRandom3.h"
#include "TMath.h"
#include "TClonesArray.h"
#include "ROOT/RDataFrame.hxx"

// ─── Physical constants ─────────────────────────────────────────────────────
static constexpr double AMU_GeV = 0.93149410242;
static constexpr double Z_FRAG = 8.;
static constexpr double SEL_Z_MIN = 7.6;
static constexpr double SEL_Z_MAX = 8.5;
static constexpr double MN_GeV = 0.939565420;
static constexpr double C_CM_PER_NS = 29.9792458;
static constexpr double FIB33_OFF = 25.;
static constexpr double FIB31_OFF = -25.;

// ─── 25F incoming graphical cut (polygon from TCutG) ────────────────────────
static const std::vector<std::pair<double, double>> INCOMING_25F_POLYGON = {
    {2.77227, 8.55258},
    {2.77403, 8.36228},
    {2.77587, 8.24598},
    {2.77849, 8.21426},
    {2.78194, 8.31470},
    {2.78276, 8.38871},
    {2.78512, 8.69531},
    {2.78587, 9.09178},
    {2.78490, 9.36138},
    {2.78295, 9.53582},
    {2.78055, 9.77370},
    {2.77834, 9.95344},
    {2.77647, 10.10670},
    {2.77534, 10.14900},
    {2.77272, 10.11730},
    {2.77059, 9.95344},
    {2.76924, 9.71027},
    {2.76905, 9.46182},
    {2.76939, 9.19750},
    {2.77014, 8.92791},
    {2.77100, 8.74289},
    {2.77182, 8.58959},
    {2.77227, 8.55258}};

// ─── Structs ────────────────────────────────────────────────────────────────

struct NFirst
{
        int idx;
        double tns;
        TVector3 pos;
};

/// All PID / kinematic cuts for a single reaction channel.
struct ReactionConfig
{
        double aoqOutMin = 0;
        double aoqOutMax = 0;
        double massAMU = 0;
        std::string outName;
        bool hasNeutrons = false;
        bool isUnreacted = false;

        // Optional elliptical outgoing-PID cut
        bool useEllipse = false;
        double ell_mu_AoQ = 0;
        double ell_sig_AoQ = 0;
        double ell_mu_Z = 0;
        double ell_sig_Z = 0;
        double ell_k = 0;
};

/// CALIFA top-2 cluster result
struct Top2Clusters
{
        double e1, th1, ph1, e2, th2, ph2;
        bool good;
};

// ─── Free helpers ───────────────────────────────────────────────────────────

/// Thread-safe point-in-polygon test (ray-casting algorithm).
static bool isInsidePolygon(double x, double y,
                            const std::vector<std::pair<double, double>> &poly)
{
        bool inside = false;
        const int n = (int)poly.size();
        for (int i = 0, j = n - 1; i < n; j = i++)
        {
                double xi = poly[i].first, yi = poly[i].second;
                double xj = poly[j].first, yj = poly[j].second;
                if (((yi > y) != (yj > y)) &&
                    (x < (xj - xi) * (y - yi) / (yj - yi) + xi))
                        inside = !inside;
        }
        return inside;
}

static bool insideEllipse(double aoq, double z,
                          double mu_AoQ, double sig_AoQ,
                          double mu_Z, double sig_Z, double k)
{
        double d2 = std::pow((aoq - mu_AoQ) / sig_AoQ, 2) +
                    std::pow((z - mu_Z) / sig_Z, 2);
        return d2 < k * k;
}

static inline void trim_inplace(std::string &s)
{
        auto notSpace = [](unsigned char c)
        { return !std::isspace(c); };
        s.erase(s.begin(), std::find_if(s.begin(), s.end(), notSpace));
        s.erase(std::find_if(s.rbegin(), s.rend(), notSpace).base(), s.end());
}

static inline void strip_quotes_inplace(std::string &s)
{
        if (s.size() >= 2 &&
            ((s.front() == '"' && s.back() == '"') ||
             (s.front() == '\'' && s.back() == '\'')))
                s = s.substr(1, s.size() - 2);
}

static std::vector<std::string> read_file_list(const std::string &txtPath)
{
        std::ifstream in(txtPath);
        if (!in.is_open())
                throw std::runtime_error("Cannot open file list: " + txtPath);

        std::vector<std::string> files;
        files.reserve(256);
        std::string line;

        while (std::getline(in, line))
        {
                auto hashPos = line.find('#');
                if (hashPos != std::string::npos)
                        line = line.substr(0, hashPos);
                trim_inplace(line);
                if (line.empty())
                        continue;
                strip_quotes_inplace(line);
                trim_inplace(line);
                if (!line.empty())
                        files.push_back(line);
        }

        if (files.empty())
                throw std::runtime_error("File list is empty: " + txtPath);
        return files;
}

/// Incoming 25F selection using the graphical cut polygon (thread-safe).
static bool isGoodIncoming(double AoQ, double Z)
{
        return isInsidePolygon(AoQ, Z, INCOMING_25F_POLYGON);
}

// ─── Reaction configuration factory ────────────────────────────────────────

static ReactionConfig makeReactionConfig(const TString &reaction)
{
        ReactionConfig cfg;

        if (reaction == "25F23O")
        {
                cfg = {2.785, 2.88,
                       23.015696686 - 8 * 0.00511,
                       "data_23O", true, false};
        }
        else if (reaction == "25F22O")
        {
                cfg = {2.67, 2.755,
                       22.009965744 - 8.0 * 0.00511,
                       "data_22O", true, false};
        }
        else if (reaction == "25F24O")
        {
                cfg = {2.89, 3.02,
                       24.019861000 - 8 * 0.00511,
                       "data_24O", false, false,
                       true, 2.95766, 0.0221273, 8.06447, 0.201233, 3.};
        }
        else if (reaction == "25Fp2p")
        {
                std::cout << "[WARN] Testing only - not for kinematic extraction.\n";
                cfg = {2.3, 3.1,
                       24.019861000 - 8 * 0.00511,
                       "data_25Fp2p", false, false};
        }
        else if (reaction == "25F25F")
        {
                cfg = {2.71, 2.77,
                       1.0 /* dummy */,
                       "data_25F", false, true,
                       true, 2.740, 0.00988, 8.976, 0.1973, 3.};
        }

        return cfg;
}

// ─── Reusable RDataFrame column builders ────────────────────────────────────

/// FRS incoming columns + filter (now uses graphical cut polygon)
static ROOT::RDF::RNode defineFrsIncoming(ROOT::RDF::RNode node)
{
        return node
            .Define("in_AoQ", [](TClonesArray &f)
                    { return ((R3BFrsData *)f.UncheckedAt(0))->GetAq(); }, {"FrsData"})
            .Define("in_Z", [](TClonesArray &f)
                    { return ((R3BFrsData *)f.UncheckedAt(0))->GetZ(); }, {"FrsData"})
            .Define("beta_proj", [](TClonesArray &f)
                    { return ((R3BFrsData *)f.UncheckedAt(0))->GetBeta(); }, {"FrsData"})
            .Filter([](double aoq, double z)
                    { return isGoodIncoming(aoq, z); },
                    {"in_AoQ", "in_Z"}, "incoming 25F graphical cut");
}

/// Outgoing fragment PID from MDF
static ROOT::RDF::RNode defineFragmentPID(ROOT::RDF::RNode node)
{
        return node
            .Define("frag", [](TClonesArray &m)
                    { return (R3BTrackingParticle *)m.UncheckedAt(0); }, {"FragmentMDFTrack"})
            .Define("AoQ_frag", [](R3BTrackingParticle *t)
                    { return (double)t->GetMass(); }, {"frag"})
            .Define("Z_frag_est", [](R3BTrackingParticle *t)
                    { return (double)t->GetCharge(); }, {"frag"});
}

/// FOOT charge (foot_z_1..8) and position (PosFoot5..8) using GetHitIndexByName
static ROOT::RDF::RNode defineFootColumns(ROOT::RDF::RNode node)
{
        // Charge array: feet 1-4 from IncomingTrackFoot, feet 5-8 from OutgoingTrackFoot
        node = node.Define("foot_z_arr",
                           [](TClonesArray &infoot, TClonesArray &outfoot, TClonesArray &hitFoot)
                           {
                std::array<double, 8> z;
                z.fill(-1.0);

                // Feet 1-4 from IncomingTrackFoot
                if (infoot.GetEntriesFast() > 0)
                {
                    auto *trk = static_cast<R3BTrackingParticle *>(infoot.UncheckedAt(0));
                    for (int i = 0; i < 4; ++i)
                    {
                        int idx = trk->GetHitIndexByName(Form("foot%d", i + 1));
                        if (idx >= 0 && idx < hitFoot.GetEntriesFast())
                        {
                            auto *h = static_cast<R3BFootHitData *>(hitFoot.UncheckedAt(idx));
                            z[i] = h->GetZCharge();
                        }
                    }
                }

                // Feet 5-8 from OutgoingTrackFoot
                if (outfoot.GetEntriesFast() > 0)
                {
                    auto *trk = static_cast<R3BTrackingParticle *>(outfoot.UncheckedAt(0));
                    for (int i = 4; i < 8; ++i)
                    {
                        int idx = trk->GetHitIndexByName(Form("foot%d", i + 1));
                        if (idx >= 0 && idx < hitFoot.GetEntriesFast())
                        {
                            auto *h = static_cast<R3BFootHitData *>(hitFoot.UncheckedAt(idx));
                            z[i] = h->GetZCharge();
                        }
                    }
                }

                return z; }, {"IncomingTrackFoot", "OutgoingTrackFoot", "FootHitData"});

        for (int i = 0; i < 8; ++i)
        {
                const int idx = i;
                node = node.Define(Form("foot_z_%d", i + 1),
                                   [idx](const std::array<double, 8> &z)
                                   { return z[idx]; },
                                   {"foot_z_arr"});
        }

        // Positions for feet 5-8 from OutgoingTrackFoot
        node = node.Define("foot_pos_arr",
                           [](TClonesArray &outfoot, TClonesArray &hitFoot)
                           {
                std::array<double, 8> pos;
                pos.fill(-999.0);

                if (outfoot.GetEntriesFast() > 0)
                {
                    auto *trk = static_cast<R3BTrackingParticle *>(outfoot.UncheckedAt(0));
                    for (int i = 4; i < 8; ++i)
                    {
                        int idx = trk->GetHitIndexByName(Form("foot%d", i + 1));
                        if (idx >= 0 && idx < hitFoot.GetEntriesFast())
                        {
                            auto *h = static_cast<R3BFootHitData *>(hitFoot.UncheckedAt(idx));
                            pos[i] = h->GetPos();
                        }
                    }
                }

                return pos; }, {"OutgoingTrackFoot", "FootHitData"});

        for (int i = 4; i < 8; ++i)
        {
                const int idx = i;
                node = node.Define(Form("PosFoot%d", i + 1),
                                   [idx](const std::array<double, 8> &p)
                                   { return p[idx]; },
                                   {"foot_pos_arr"});
        }

        return node;
}

/// Generic best-hit extractor for a single fiber detector using GetHitIndexByName.
static ROOT::RDF::RNode defineSingleFiber(ROOT::RDF::RNode node,
                                          const std::string &hitBranch,
                                          const std::string &prefix,
                                          const std::string &elossName,
                                          double xOffset,
                                          const std::string &trackBranch,
                                          const std::string &hitName)
{
        const std::string arrCol = prefix + "_best";

        node = node.Define(arrCol,
                           [xOffset, hitName](TClonesArray &trks, TClonesArray &hits)
                           {
                std::array<double, 3> res = {-999.0, -999.0, -999.0};

                if (trks.GetEntriesFast() > 0)
                {
                    auto *trk = static_cast<R3BTrackingParticle *>(trks.UncheckedAt(0));
                    int idx = trk->GetHitIndexByName(hitName.c_str());
                    if (idx >= 0 && idx < hits.GetEntriesFast())
                    {
                        auto *hit = static_cast<R3BFiberMAPMTHitData *>(hits.UncheckedAt(idx));
                        res = {hit->GetX() + xOffset, hit->GetY(), hit->GetEloss()};
                    }
                }

                return res; }, {trackBranch, hitBranch});

        node = node.Define(prefix + "X",
                           [](const std::array<double, 3> &a)
                           { return a[0]; }, {arrCol})
                   .Define(prefix + "Y",
                           [](const std::array<double, 3> &a)
                           { return a[1]; }, {arrCol})
                   .Define(elossName,
                           [](const std::array<double, 3> &a)
                           { return a[2]; }, {arrCol});
        return node;
}

/// All four fiber detectors + ToFD
static ROOT::RDF::RNode defineFibersAndTofd(ROOT::RDF::RNode node)
{
        // Fi30 and Fi31: incoming fibers -> IncomingTrackFoot
        node = defineSingleFiber(node, "Fi30Hit", "fib30", "ElossFib30", 0.0,
                                 "IncomingTrackFoot", "Fi30Hit");
        node = defineSingleFiber(node, "Fi31Hit", "fib31", "ElossFib31", FIB31_OFF,
                                 "IncomingTrackFoot", "Fi31Hit");

        // Fi32 and Fi33: outgoing fibers -> OutgoingTrackFoot
        node = defineSingleFiber(node, "Fi32Hit", "fib32", "ElossFib32", 0.0,
                                 "OutgoingTrackFoot", "Fi32Hit");
        node = defineSingleFiber(node, "Fi33Hit", "fib33", "ElossFib33", FIB33_OFF,
                                 "OutgoingTrackFoot", "Fi33Hit");

        node = node.Define("tofdX", [](TClonesArray &t)
                           {
        for (int i = 0; i < t.GetEntriesFast(); ++i) {
            auto *h = (R3BTofdHitData *)t.UncheckedAt(i);
            if (h->GetDetId() == 1) return (double)h->GetX();
        }
        return -999.0; }, {"TofdHit"});

        return node;
}

/// Outgoing start position from OutgoingTrackFoot
static ROOT::RDF::RNode defineOutgoingStartPos(ROOT::RDF::RNode node)
{
        return node
            .Define("out_startpos", [](TClonesArray &ot)
                    {
            std::array<double, 3> res = {-999.0, -999.0, -999.0};
            if (ot.GetEntriesFast() > 0) {
                auto *trk = (R3BTrackingParticle *)ot.UncheckedAt(0);
                TVector3 pos = trk->GetStartPosition();
                res = {pos.X(), pos.Y(), pos.Z()};
            }
            return res; }, {"OutgoingTrackFoot"})
            .Define("out_startX", [](const std::array<double, 3> &p)
                    { return p[0]; }, {"out_startpos"})
            .Define("out_startY", [](const std::array<double, 3> &p)
                    { return p[1]; }, {"out_startpos"})
            .Define("out_startZ", [](const std::array<double, 3> &p)
                    { return p[2]; }, {"out_startpos"});
}

// ─── CALIFA ─────────────────────────────────────────────────────────────────

static Top2Clusters findTop2Clusters(TClonesArray &clu,
                                     double EminClu = 20e3,
                                     double EminP = 20e3)
{
        Top2Clusters c{-1, -999, -999, -1, -999, -999, false};

        for (int i = 0; i < clu.GetEntriesFast(); ++i)
        {
                auto *hit = (R3BCalifaClusterData *)clu.UncheckedAt(i);
                double E = hit->GetEnergy();
                if (E < EminClu)
                        continue;
                double th = hit->GetTheta(), ph = hit->GetPhi();

                if (E > c.e1)
                {
                        c.e2 = c.e1;
                        c.th2 = c.th1;
                        c.ph2 = c.ph1;
                        c.e1 = E;
                        c.th1 = th;
                        c.ph1 = ph;
                }
                else if (E > c.e2)
                {
                        c.e2 = E;
                        c.th2 = th;
                        c.ph2 = ph;
                }
        }

        if (c.e1 < EminP || c.e2 < EminP)
                return c;

        double dphiDeg = std::abs(TVector2::Phi_mpi_pi(c.ph2 - c.ph1)) * TMath::RadToDeg();
        // if (std::abs(dphiDeg - 180.0) > 30.0)
        //         return c;

        if (std::abs(dphiDeg - 180.0) > 60.0)
                return c;

        c.good = true;
        return c;
}

static ROOT::RDF::RNode defineCalifaColumns(ROOT::RDF::RNode node)
{
        node = node
                   .Define("califa_top2", [](TClonesArray &clu)
                           { return findTop2Clusters(clu); },
                           {"CalifaClusterData"})
                   .Define("califa_opa", [](const Top2Clusters &c)
                           {
            if (!c.good) return -999.0;
            TVector3 v1; v1.SetMagThetaPhi(1.0, c.th1, c.ph1);
            TVector3 v2; v2.SetMagThetaPhi(1.0, c.th2, c.ph2);
            return v1.Angle(v2); }, {"califa_top2"})
                   .Define("califa_swap", []()
                           {
            static thread_local TRandom3 rng(0);
            return rng.Rndm() < 0.5; })
                   .Define("califa_theta_L", [](bool sw, const Top2Clusters &c)
                           { return (!c.good) ? -999.0 : (sw ? c.th2 : c.th1); }, {"califa_swap", "califa_top2"})
                   .Define("califa_phi_L", [](bool sw, const Top2Clusters &c)
                           { return (!c.good) ? -999.0 : (sw ? c.ph2 : c.ph1); }, {"califa_swap", "califa_top2"})
                   .Define("califa_theta_R", [](bool sw, const Top2Clusters &c)
                           { return (!c.good) ? -999.0 : (sw ? c.th1 : c.th2); }, {"califa_swap", "califa_top2"})
                   .Define("califa_phi_R", [](bool sw, const Top2Clusters &c)
                           { return (!c.good) ? -999.0 : (sw ? c.ph1 : c.ph2); }, {"califa_swap", "califa_top2"});

        return node.Filter([](double opa)
                           { return opa > -990.0; },
                           {"califa_opa"}, "good CALIFA event");
}

// ─── Neutrons ───────────────────────────────────────────────────────────────

static ROOT::RDF::RNode defineNeutronColumns(ROOT::RDF::RNode node)
{
        return node
            .Define("n_first", [](TClonesArray &nl)
                    {
            int idx = -1;  double bestZ = 1e99, bestT = -1;
            TVector3 bestPos(0, 0, 0);
            for (int i = 0; i < nl.GetEntriesFast(); ++i)
            {
                auto *h = static_cast<R3BNeulandHit *>(nl.UncheckedAt(i));
                double z = h->GetPosition().Z(), t = h->GetT();
                if (t <= 76.0 && z < bestZ)
                {
                    bestZ = z;  bestT = t;  idx = i;
                    TVector3 pos = h->GetPosition();
                    int paddle = h->GetPaddle();
                    static thread_local TRandom3 rng(0);
                    if (((paddle / 50) % 2) == 0)
                        bestPos.SetXYZ(pos.X(),
                                       pos.Y() + rng.Uniform(-2.5, 2.5),
                                       pos.Z() + rng.Uniform(-2.5, 2.5));
                    else
                        bestPos.SetXYZ(pos.X() + rng.Uniform(-2.5, 2.5),
                                       pos.Y(),
                                       pos.Z() + rng.Uniform(-2.5, 2.5));
                }
            }
            return NFirst{idx, bestT, bestPos}; }, {"NeulandHits"})
            .Filter([](const NFirst &nf)
                    { return nf.idx >= 0; }, {"n_first"})
            .Define("beta_neu", [](const NFirst &nf)
                    {
            double L = nf.pos.Mag(), t = nf.tns;
            return (L > 0 && t > 0) ? L / (C_CM_PER_NS * t) : -1.0; }, {"n_first"})
            .Filter([](double b)
                    { return b > 0 && b < 1; }, {"beta_neu"})
            .Define("n_dir", [](const NFirst &nf)
                    { return nf.pos.Mag() > 0 ? nf.pos.Unit() : TVector3(0, 0, 1); }, {"n_first"})
            .Define("p_neu", [](double bn)
                    {
            double g = 1.0 / std::sqrt(1.0 - bn * bn);
            return g * MN_GeV * bn; }, {"beta_neu"})
            .Define("P4_neu", [](double p, const TVector3 &dir, double bn)
                    {
            double g = 1.0 / std::sqrt(1.0 - bn * bn);
            return TLorentzVector(dir * p, g * MN_GeV); }, {"p_neu", "n_dir", "beta_neu"})
            .Define("px_neu", [](const TLorentzVector &N)
                    { return N.Px(); }, {"P4_neu"})
            .Define("py_neu", [](const TLorentzVector &N)
                    { return N.Py(); }, {"P4_neu"})
            .Define("pz_neu", [](const TLorentzVector &N)
                    { return N.Pz(); }, {"P4_neu"})
            .Define("x_neu_hit", [](const NFirst &nf)
                    { return nf.pos.X(); }, {"n_first"})
            .Define("y_neu_hit", [](const NFirst &nf)
                    { return nf.pos.Y(); }, {"n_first"})
            .Define("z_neu_hit", [](const NFirst &nf)
                    { return nf.pos.Z(); }, {"n_first"})
            .Define("tof_neuland", [](const NFirst &nf)
                    { return nf.tns; }, {"n_first"});
}

// ─── Output column lists ────────────────────────────────────────────────────

static std::vector<std::string> detectorColumns()
{
        return {
            "foot_z_1", "foot_z_2", "foot_z_3", "foot_z_4",
            "foot_z_5", "foot_z_6", "foot_z_7", "foot_z_8",
            "PosFoot5", "PosFoot6", "PosFoot7", "PosFoot8",
            "fib30X", "fib30Y", "ElossFib30",
            "fib31X", "fib31Y", "ElossFib31",
            "fib32X", "fib32Y", "ElossFib32",
            "fib33X", "fib33Y", "ElossFib33",
            "tofdX",
            "out_startX", "out_startY", "out_startZ"};
}

static std::vector<std::string> buildOutputColumns(bool hasNeutrons)
{
        std::vector<std::string> cols = {
            "AoQ_frag", "Z_frag_est", "A_frag", "M_frag",
            "beta_frag", "p_frag"};

        if (hasNeutrons)
                cols.insert(cols.end(), {"beta_neu", "p_neu"});

        cols.insert(cols.end(), {"califa_opa",
                                 "califa_theta_L", "califa_phi_L",
                                 "califa_theta_R", "califa_phi_R",
                                 "px_frag", "py_frag", "pz_frag",
                                 "beta_proj"});

        if (hasNeutrons)
                cols.insert(cols.end(), {"px_neu", "py_neu", "pz_neu",
                                         "x_neu_hit", "y_neu_hit", "z_neu_hit", "tof_neuland"});

        cols.insert(cols.end(), {"px_in", "py_in", "pz_in"});

        auto det = detectorColumns();
        cols.insert(cols.end(), det.begin(), det.end());
        return cols;
}

static std::vector<std::string> buildUnreactedColumns()
{
        std::vector<std::string> cols = {
            "in_AoQ", "in_Z", "AoQ_frag", "Z_frag_est", "beta_proj"};
        auto det = detectorColumns();
        cols.insert(cols.end(), det.begin(), det.end());
        return cols;
}

// ═══════════════════════════════════════════════════════════════════════════
//  MAIN
// ═══════════════════════════════════════════════════════════════════════════

void eventFilter(std::string setting = "",
                 TString reaction = "",
                 bool test = false,
                 bool append = false)
{
        ROOT::EnableImplicitMT(16);

        // ── Load file list ──────────────────────────────────────────────────
        const std::string listTxt =
            std::string(getenv("repopath")) + "/final/settings/" + setting;

        std::vector<std::string> files;
        try
        {
                files = read_file_list(listTxt);
        }
        catch (const std::exception &e)
        {
                std::cerr << "[ERROR] " << e.what() << "\n";
                return;
        }
        std::cout << "[OK] Loaded " << files.size()
                  << " ROOT files from " << listTxt << "\n";

        // ── Reaction configuration ──────────────────────────────────────────
        ReactionConfig cfg = makeReactionConfig(reaction);

        if (cfg.outName.empty())
        {
                std::cerr << "[ERROR] Unknown/invalid reaction: " << reaction << "\n";
                return;
        }

        std::string outFileName = cfg.outName + (test ? "_test_60opa" : "") + ".root";
        double M_FRAG_GeV = cfg.massAMU * AMU_GeV;

        // ── Open output file ────────────────────────────────────────────────
        std::string outFile =
            std::string(getenv("repopath")) + "/results/final/" + outFileName;
        bool exists = (gSystem->AccessPathName(outFile.c_str()) == kFALSE);
        TFile *fout = new TFile(outFile.c_str(),
                                (append && exists) ? "UPDATE" : "RECREATE");

        // ── Base filters ────────────────────────────────────────────────────
        ROOT::RDataFrame df("evt", files);

        ROOT::RDF::RNode base = df
                                    .Filter([](TClonesArray &f)
                                            { return f.GetEntriesFast() > 0; }, {"FrsData"})
                                    .Filter([](TClonesArray &m)
                                            { return m.GetEntriesFast() > 0; }, {"FragmentMDFTrack"})
                                    .Filter([](TClonesArray &i)
                                            { return i.GetEntriesFast() > 0; }, {"IncomingTrackFoot"});

        if (!cfg.isUnreacted)
                base = base.Filter([](TClonesArray &c)
                                   { return c.GetEntriesFast() >= 2; },
                                   {"CalifaClusterData"});

        ROOT::RDF::RNode df_cut = base;
        if (cfg.hasNeutrons)
                df_cut = df_cut.Filter([](TClonesArray &n)
                                       { return n.GetEntriesFast() > 0; },
                                       {"NeulandHits"});

        // ── FRS incoming (graphical cut) + fragment PID ─────────────────────
        auto df_frs = defineFrsIncoming(df_cut);
        auto df_frag = defineFragmentPID(df_frs);

        // ═════════════════  UNREACTED 25F PATH  ══════════════════════════════
        if (cfg.isUnreacted)
        {
                double unr_mu_AoQ = cfg.ell_mu_AoQ;
                double unr_sig_AoQ = cfg.ell_sig_AoQ;
                double unr_mu_Z = cfg.ell_mu_Z;
                double unr_sig_Z = cfg.ell_sig_Z;
                double unr_k = cfg.ell_k;

                auto df_sel = df_frag
                                  .Filter([=](double aoq, double z)
                                          { return insideEllipse(aoq, z,
                                                                 unr_mu_AoQ, unr_sig_AoQ, unr_mu_Z, unr_sig_Z, unr_k); }, {"AoQ_frag", "Z_frag_est"}, "outgoing 25F (unreacted)");

                auto df_det = defineFibersAndTofd(defineFootColumns(df_sel));
                df_det = defineOutgoingStartPos(df_det);

                df_det.Snapshot("FilterDataTree", outFile, buildUnreactedColumns());
                std::cout << "\n[OK] TTree (unreacted 25F) saved in: " << outFile << "\n";
                if (fout)
                        fout->Close();
                return;
        }

        // ═════════════════  STANDARD REACTION PATH  ══════════════════════════

        // Full fragment kinematics
        auto df_kin = df_frag
                          .Define("A_frag", [](double aoq, double z)
                                  { return aoq * z; },
                                  {"AoQ_frag", "Z_frag_est"})
                          .Define("P_frag_PoQ_vec", [](R3BTrackingParticle *t)
                                  { return t->GetStartMomentum(); }, {"frag"})
                          .Define("P_over_Q_frag", [](const TVector3 &p)
                                  { return p.Mag(); },
                                  {"P_frag_PoQ_vec"})
                          .Define("P_frag_dir", [](const TVector3 &p)
                                  { return p.Mag() > 0 ? p.Unit() : TVector3(0, 0, 1); }, {"P_frag_PoQ_vec"})
                          .Define("P_in_vec", [](TClonesArray &it)
                                  { return ((R3BTrackingParticle *)it.UncheckedAt(0))->GetStartMomentum(); }, {"IncomingTrackFoot"})
                          .Define("px_in", [](const TVector3 &p)
                                  { return p.X(); }, {"P_in_vec"})
                          .Define("py_in", [](const TVector3 &p)
                                  { return p.Y(); }, {"P_in_vec"})
                          .Define("pz_in", [](const TVector3 &p)
                                  { return p.Z(); }, {"P_in_vec"})
                          .Define("Z_frag", []()
                                  { return Z_FRAG; }, {})
                          .Define("M_frag", [M_FRAG_GeV]()
                                  { return M_FRAG_GeV; }, {})
                          .Define("p_frag", [](double pq, double Z)
                                  { return (pq > 0 && Z > 0) ? pq * Z : 0.0; }, {"P_over_Q_frag", "Z_frag"})
                          .Define("beta_frag", [](double p, double M)
                                  {
            if (p <= 0 || M <= 0) return 0.0;
            return p / std::sqrt(p * p + M * M); }, {"p_frag", "M_frag"})
                          .Define("P4_frag", [](double p, double Mf, const TVector3 &dir)
                                  {
            if (p <= 0) return TLorentzVector(0, 0, 0, 0);
            return TLorentzVector(dir * p, std::sqrt(p * p + Mf * Mf)); }, {"p_frag", "M_frag", "P_frag_dir"});

        // ── Outgoing PID cut ────────────────────────────────────────────────
        auto df_filtered = df_kin.Filter(
            [&cfg](double aoq, double z)
            {
                    if (cfg.useEllipse)
                            return insideEllipse(aoq, z,
                                                 cfg.ell_mu_AoQ, cfg.ell_sig_AoQ,
                                                 cfg.ell_mu_Z, cfg.ell_sig_Z, cfg.ell_k);
                    return (aoq >= cfg.aoqOutMin && aoq <= cfg.aoqOutMax &&
                            z >= SEL_Z_MIN && z <= SEL_Z_MAX);
            },
            {"AoQ_frag", "Z_frag_est"});

        // ── Neutrons (conditional) ──────────────────────────────────────────
        ROOT::RDF::RNode df_post_n = df_filtered;
        if (cfg.hasNeutrons)
                df_post_n = defineNeutronColumns(df_post_n);

        // ── CALIFA + detectors + fragment momentum ──────────────────────────
        auto df_cal = defineCalifaColumns(df_post_n);
        auto df_det = defineFibersAndTofd(defineFootColumns(df_cal));

        auto df_out = df_det
                          .Define("px_frag", [](const TLorentzVector &F)
                                  { return F.Px(); }, {"P4_frag"})
                          .Define("py_frag", [](const TLorentzVector &F)
                                  { return F.Py(); }, {"P4_frag"})
                          .Define("pz_frag", [](const TLorentzVector &F)
                                  { return F.Pz(); }, {"P4_frag"});

        df_out = defineOutgoingStartPos(df_out);

        // ── Snapshot ────────────────────────────────────────────────────────
        df_out.Snapshot("FilterDataTree", outFile,
                        buildOutputColumns(cfg.hasNeutrons));

        std::cout << "\n[OK] TTree saved in: " << outFile << "\n";
        if (fout)
                fout->Close();
}