// ===========================================================================
//  studyRandomizeNeuland.C
// ---------------------------------------------------------------------------
//  Study of the effect of the NeuLAND hit-position randomization on Erel.
//
//  Input : ROOT file produced by the *modified* eventFilter, i.e. the one that
//          stores the CENTRAL bar position (x_neu_hit, y_neu_hit, z_neu_hit),
//          the raw ToF (tof_neuland) and the paddle id (paddle_neu).
//
//  What it does
//  ------------
//   1) Recomputes Erel  (a) with NO randomization  and
//                       (b) with 4 independent randomizations,
//      and overlays the 5 spectra in the same canvas (+ ratio pad).
//
//   2) Runs <nReplicas> independent randomizations, and builds the Erel
//      spectrum with the randomization uncertainty attached:
//         content(bin) = mean over replicas
//         error(bin)   = RMS over replicas  (systematic band)
//      It also produces the per-event <Erel> spectrum and the distribution of
//      the per-event sigma(Erel) induced by the randomization.
//
//  Usage
//  -----
//    root -l 'studyRandomizeNeuland.C+("/nucl_lustre/pablogrusell/g249/g249_analysis/results/dataFiles/data_23O_no_rand.root")'
//    root -l 'studyRandomizeNeuland.C+("...root","23O1n.txt",200)'   // more replicas
//    root -l 'studyRandomizeNeuland.C+("...root","23O1n.txt",100,100,0,10,200000)' // quick test
//
//  The Erel formula, the beta-matching shift and the transverse-momentum
//  offsets are exactly the ones used in DataAnalysis::getData().
// ===========================================================================

#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <limits>
#include <string>
#include <vector>

#include "TCanvas.h"
#include "TFile.h"
#include "TH1D.h"
#include "TLatex.h"
#include "TLeaf.h"
#include "TLegend.h"
#include "TMath.h"
#include "TPad.h"
#include "TRandom3.h"
#include "TStyle.h"
#include "TTree.h"
#include "TVector3.h"

// ---------------------------------------------------------------------------
//  Constants / knobs
// ---------------------------------------------------------------------------
namespace nrand
{
    static constexpr double kC_cm_ns = 29.9792458; // speed of light [cm/ns]
    static constexpr double kMn_GeV = 0.939565420; // neutron mass   [GeV]
    static constexpr double kErelShiftMeV = 4.2;

    // Half width of the uniform smearing, in cm.
    //  - kHalfBar   : along the bar (the coordinate that is NOT measured by the
    //                 time difference of the two PMTs) -> half of the bar width
    //  - kHalfDepth : along z (bar thickness)
    static constexpr double kHalfBar = 2.5;
    static constexpr double kHalfDepth = 2.5;

    // The original eventFilter used  ((paddle / 50) % 2) == 0  with 1-based
    // paddle ids, which puts paddle 50 in the "wrong" plane group.
    // Keep it false to reproduce EXACTLY the old (randomized) production;
    // set it to true to use the (more correct) ((paddle-1)/50) % 2.
    static bool kFixPaddleGrouping = false;

    struct Offsets
    {
        double fragX = 0., fragY = 0.; // fFragmentPOffsets[0], [1]
        double neuX = 0., neuY = 0.;   // fNeutronPOffsets[0], [1]
        double dBeta = 0.;             // fBetaMatchValue
    };

    // Everything needed, per event, to recompute Erel for any smearing.
    struct Evt
    {
        double x, y, z;        // central NeuLAND hit position [cm]
        double tof;            // NeuLAND ToF [ns] (already offset-corrected)
        int paddle;            // NeuLAND paddle id
        double betaFrag;       // beta_frag + beta-match shift
        double gammaFrag;      // 1/sqrt(1-beta^2)
        double mFrag;          // M_frag [GeV]
        double fxFrag, fyFrag; // px/pz - offset , py/pz - offset
    };

    // -----------------------------------------------------------------------
    inline TVector3 smearPos(const Evt &e, TRandom3 &rng)
    {
        const int plane = kFixPaddleGrouping ? (e.paddle - 1) / 50 : e.paddle / 50;
        if ((plane % 2) == 0) // horizontal bar: x is measured, smear y and z
            return TVector3(e.x,
                            e.y + rng.Uniform(-kHalfBar, kHalfBar),
                            e.z + rng.Uniform(-kHalfDepth, kHalfDepth));
        // vertical bar: y is measured, smear x and z
        return TVector3(e.x + rng.Uniform(-kHalfBar, kHalfBar),
                        e.y,
                        e.z + rng.Uniform(-kHalfDepth, kHalfDepth));
    }

    // Erel [MeV] for a given neutron position; returns NaN if unphysical.
    inline double erelMeV(const Evt &e, const TVector3 &pos, const Offsets &off)
    {
        const double L = pos.Mag();
        if (L <= 0. || e.tof <= 0.)
            return std::numeric_limits<double>::quiet_NaN();

        const double betaN = L / (kC_cm_ns * e.tof);
        if (betaN <= 0. || betaN >= 1.)
            return std::numeric_limits<double>::quiet_NaN();

        if (pos.Z() == 0.)
            return std::numeric_limits<double>::quiet_NaN();

        // Same construction as DataAnalysis::getData(): direction cosines built
        // from the px/pz , py/pz ratios (identical to x/z , y/z) minus offsets.
        const double dx = pos.X() / pos.Z() - off.neuX;
        const double dy = pos.Y() / pos.Z() - off.neuY;

        const double cosAng =
            (dx * e.fxFrag + dy * e.fyFrag + 1.0) /
            (std::sqrt(dx * dx + dy * dy + 1.0) *
             std::sqrt(e.fxFrag * e.fxFrag + e.fyFrag * e.fyFrag + 1.0));

        const double gammaN = 1.0 / std::sqrt(1.0 - betaN * betaN);
        const double mf = e.mFrag;

        const double erel =
            std::sqrt(mf * mf + kMn_GeV * kMn_GeV +
                      2.0 * gammaN * e.gammaFrag * mf * kMn_GeV *
                          (1.0 - betaN * e.betaFrag * cosAng)) -
            mf - kMn_GeV;

        return erel * 1000.; // GeV -> MeV
    }

    // -----------------------------------------------------------------------
    inline Offsets readOffsets(const char *offFile)
    {
        Offsets o;
        if (!offFile || std::string(offFile).empty())
        {
            std::cout << "[offsets] none requested -> all offsets = 0\n";
            return o;
        }

        const char *repo = getenv("repopath");
        std::string path = offFile;
        if (repo && path.find('/') == std::string::npos)
            path = std::string(repo) + "/final/settings/" + offFile;

        std::ifstream in(path.c_str());
        if (!in.is_open())
        {
            std::cerr << "[offsets] WARNING: cannot open " << path
                      << " -> all offsets = 0\n";
            return o;
        }

        std::vector<double> v;
        std::string line;
        while (std::getline(in, line))
        {
            if (line.empty())
                continue;
            v.push_back(std::atof(line.c_str()));
        }
        if (v.size() < 5)
        {
            std::cerr << "[offsets] WARNING: " << path << " has only " << v.size()
                      << " values (5 expected) -> all offsets = 0\n";
            return o;
        }

        o.fragX = v[0];
        o.fragY = v[1];
        o.neuX = v[2];
        o.neuY = v[3];
        o.dBeta = v[4];

        std::cout << "[offsets] " << path << "\n"
                  << "          frag  (" << o.fragX << ", " << o.fragY << ")\n"
                  << "          neu   (" << o.neuX << ", " << o.neuY << ")\n"
                  << "          dBeta  " << o.dBeta << "\n";
        return o;
    }

    // Style helper for the overlaid spectra.
    inline void styleHist(TH1D *h, int color, int style, int width = 2)
    {
        h->SetLineColor(color);
        h->SetLineWidth(width);
        h->SetLineStyle(style);
        h->SetMarkerColor(color);
        h->SetStats(kFALSE);
        h->SetTitle("");
    }
} // namespace nrand

// ===========================================================================
//  MAIN
// ===========================================================================
void studyRandomizeNeuland(
    const char *inFile = "/nucl_lustre/pablogrusell/g249/g249_analysis/results/dataFiles/data_23O_no_rand.root",
    const char *offFile = "23O1n.txt",
    int nReplicas = 100,
    int nBins = 70,
    double eMin = 0.,
    double eMax = 17.,
    Long64_t maxEntries = -1,
    const char *outBase = "studyRandomizeNeuland")
{
    using namespace nrand;

    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);
    gStyle->SetPadTickX(1);
    gStyle->SetPadTickY(1);

    // -----------------------------------------------------------------------
    // 0) Offsets + input tree
    // -----------------------------------------------------------------------
    const Offsets off = readOffsets(offFile);

    TFile *f = TFile::Open(inFile, "READ");
    if (!f || f->IsZombie())
    {
        std::cerr << "[studyRandomizeNeuland] ERROR: cannot open " << inFile << "\n";
        return;
    }

    TTree *t = dynamic_cast<TTree *>(f->Get("FilterDataTree"));
    if (!t)
    {
        std::cerr << "[studyRandomizeNeuland] ERROR: no TTree \"FilterDataTree\" in "
                  << inFile << "\n";
        f->Close();
        return;
    }

    if (!t->GetBranch("paddle_neu"))
    {
        std::cerr << "[studyRandomizeNeuland] ERROR: branch \"paddle_neu\" not found.\n"
                  << "  This macro needs the output of the *modified* eventFilter\n"
                  << "  (central positions + paddle id, no randomization applied).\n";
        f->Close();
        return;
    }

    Double_t x_neu, y_neu, z_neu, tof_neu;
    Double_t beta_frag, M_frag, px_frag, py_frag, pz_frag, califa_opa;

    t->SetBranchAddress("x_neu_hit", &x_neu);
    t->SetBranchAddress("y_neu_hit", &y_neu);
    t->SetBranchAddress("z_neu_hit", &z_neu);
    t->SetBranchAddress("tof_neuland", &tof_neu);
    t->SetBranchAddress("beta_frag", &beta_frag);
    t->SetBranchAddress("M_frag", &M_frag);
    t->SetBranchAddress("px_frag", &px_frag);
    t->SetBranchAddress("py_frag", &py_frag);
    t->SetBranchAddress("pz_frag", &pz_frag);
    t->SetBranchAddress("califa_opa", &califa_opa);

    // paddle_neu is written as Int_t by RDataFrame, but stay tolerant.
    Int_t paddle_i = -1;
    Double_t paddle_d = -1.;
    bool paddleIsInt = true;
    {
        TLeaf *lf = t->GetLeaf("paddle_neu");
        const TString tn = lf ? lf->GetTypeName() : "Int_t";
        paddleIsInt = !(tn == "Double_t" || tn == "Float_t" || tn == "Double32_t");
        if (paddleIsInt)
            t->SetBranchAddress("paddle_neu", &paddle_i);
        else
            t->SetBranchAddress("paddle_neu", &paddle_d);
    }

    const Long64_t nAll = t->GetEntries();
    const Long64_t nUse = (maxEntries > 0 && maxEntries < nAll) ? maxEntries : nAll;
    std::cout << "[input] " << inFile << "\n"
              << "        entries in tree : " << nAll << "\n"
              << "        entries used    : " << nUse << "\n";

    // -----------------------------------------------------------------------
    // 1) One single pass over the tree -> cache everything needed
    // -----------------------------------------------------------------------
    std::vector<Evt> evts;
    evts.reserve(nUse);

    Long64_t nBadFrag = 0;
    for (Long64_t i = 0; i < nUse; ++i)
    {
        t->GetEntry(i);

        if (califa_opa <= 1.25 || califa_opa >= 1.65)
            continue;

        if (pz_frag == 0. || M_frag <= 0.)
        {
            ++nBadFrag;
            continue;
        }

        Evt e;
        e.x = x_neu;
        e.y = y_neu;
        e.z = z_neu;
        e.tof = tof_neu;
        e.paddle = paddleIsInt ? paddle_i : (int)paddle_d;

        e.betaFrag = beta_frag + off.dBeta;
        if (e.betaFrag <= 0. || e.betaFrag >= 1.)
        {
            ++nBadFrag;
            continue;
        }
        e.gammaFrag = 1.0 / std::sqrt(1.0 - e.betaFrag * e.betaFrag);
        e.mFrag = M_frag;
        e.fxFrag = px_frag / pz_frag - off.fragX;
        e.fyFrag = py_frag / pz_frag - off.fragY;

        evts.push_back(e);
    }
    f->Close();

    const size_t nEvt = evts.size();
    std::cout << "[input] usable events    : " << nEvt
              << "   (rejected: " << nBadFrag << ")\n\n";
    if (nEvt == 0)
        return;

    // -----------------------------------------------------------------------
    // 2) PART 1 : no randomization + 4 randomizations, overlaid
    // -----------------------------------------------------------------------
    const int kShow = 4;

    TH1D *hNoRand = new TH1D("hErel_norand", ";E_{rel} [MeV];counts / bin",
                             nBins, eMin, eMax);
    hNoRand->Sumw2();

    std::vector<TH1D *> hRand(kShow, nullptr);
    for (int k = 0; k < kShow; ++k)
    {
        hRand[k] = new TH1D(Form("hErel_rand%d", k + 1),
                            ";E_{rel} [MeV];counts / bin", nBins, eMin, eMax);
        hRand[k]->Sumw2();
    }

    // no randomization
    for (size_t i = 0; i < nEvt; ++i)
    {
        const Evt &e = evts[i];
        const double er = erelMeV(e, TVector3(e.x, e.y, e.z), off);
        if (std::isfinite(er))
            hNoRand->Fill(er + kErelShiftMeV);
    }

    // 4 independent randomizations
    for (int k = 0; k < kShow; ++k)
    {
        TRandom3 rng(20250901 + 1000 * (k + 1)); // fixed, reproducible seeds
        for (size_t i = 0; i < nEvt; ++i)
        {
            const Evt &e = evts[i];
            const double er = erelMeV(e, smearPos(e, rng), off);
            if (std::isfinite(er))
                hRand[k]->Fill(er + kErelShiftMeV);
        }
    }

    // ---- report ----
    auto peakOf = [](TH1D *h)
    { return h->GetXaxis()->GetBinCenter(h->GetMaximumBin()); };

    std::cout << "----------------------------------------------------------\n"
              << " spectrum        entries      mean[MeV]   rms[MeV]  peak[MeV]\n"
              << "----------------------------------------------------------\n";
    std::cout << Form(" no randomiz. %10.0f %11.4f %10.4f %10.4f\n",
                      hNoRand->GetEntries(), hNoRand->GetMean(),
                      hNoRand->GetRMS(), peakOf(hNoRand));
    for (int k = 0; k < kShow; ++k)
        std::cout << Form(" random #%d    %10.0f %11.4f %10.4f %10.4f\n", k + 1,
                          hRand[k]->GetEntries(), hRand[k]->GetMean(),
                          hRand[k]->GetRMS(), peakOf(hRand[k]));
    std::cout << "----------------------------------------------------------\n\n";

    // ---- draw ----
    const int col[kShow] = {kRed + 1, kBlue + 1, kGreen + 2, kMagenta + 1};

    TCanvas *c1 = new TCanvas("c1_randomizations",
                              "Erel: effect of the NeuLAND randomization", 950, 800);
    TPad *p1 = new TPad("p1", "", 0., 0.32, 1., 1.);
    TPad *p2 = new TPad("p2", "", 0., 0.00, 1., 0.32);
    p1->SetBottomMargin(0.02);
    p1->SetLeftMargin(0.12);
    p2->SetTopMargin(0.03);
    p2->SetBottomMargin(0.32);
    p2->SetLeftMargin(0.12);
    p2->SetGridy();
    p1->Draw();
    p2->Draw();

    p1->cd();
    styleHist(hNoRand, kBlack, 1, 3);
    hNoRand->GetYaxis()->SetTitleOffset(1.25);
    hNoRand->GetXaxis()->SetLabelSize(0.);
    hNoRand->SetMaximum(1.25 * hNoRand->GetMaximum());
    hNoRand->Draw("HIST");
    for (int k = 0; k < kShow; ++k)
    {
        styleHist(hRand[k], col[k], 1, 2);
        hRand[k]->Draw("HIST SAME");
    }
    hNoRand->Draw("HIST SAME");

    TLegend *leg = new TLegend(0.58, 0.60, 0.90, 0.89);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->AddEntry(hNoRand, Form("no randomization (#mu=%.3f)", hNoRand->GetMean()), "l");
    for (int k = 0; k < kShow; ++k)
        leg->AddEntry(hRand[k], Form("randomization #%d (#mu=%.3f)", k + 1, hRand[k]->GetMean()), "l");
    leg->Draw();

    TLatex tl;
    tl.SetNDC();
    tl.SetTextSize(0.040);
    tl.DrawLatex(0.15, 0.85, "^{23}O + n   E_{rel}");

    // ratio pad: randomized / non randomized
    p2->cd();
    TH1D *hRef = (TH1D *)hNoRand->Clone("hRef_forRatio");
    for (int k = 0; k < kShow; ++k)
    {
        TH1D *r = (TH1D *)hRand[k]->Clone(Form("hRatio%d", k + 1));
        r->Divide(hRef);
        styleHist(r, col[k], 1, 2);
        r->SetTitle("");
        r->GetYaxis()->SetTitle("rand / no rand");
        r->GetYaxis()->SetNdivisions(505);
        r->GetYaxis()->SetTitleSize(0.11);
        r->GetYaxis()->SetTitleOffset(0.45);
        r->GetYaxis()->SetLabelSize(0.09);
        r->GetXaxis()->SetTitle("E_{rel} [MeV]");
        r->GetXaxis()->SetTitleSize(0.13);
        r->GetXaxis()->SetTitleOffset(1.0);
        r->GetXaxis()->SetLabelSize(0.10);
        r->SetMinimum(0.5);
        r->SetMaximum(1.5);
        r->Draw(k == 0 ? "HIST" : "HIST SAME");
    }
    c1->cd();
    c1->Update();

    // -----------------------------------------------------------------------
    // 3) PART 2 : Erel with the randomization uncertainty
    // -----------------------------------------------------------------------
    if (nReplicas < 2)
        nReplicas = 2;
    std::cout << "[part 2] running " << nReplicas << " randomization replicas ...\n";

    // bin-by-bin accumulators over the replicas
    std::vector<double> sum(nBins + 2, 0.), sum2(nBins + 2, 0.);

    // per-event accumulators (mean and sigma of Erel for each event)
    std::vector<double> evSum(nEvt, 0.), evSum2(nEvt, 0.);
    std::vector<int> evN(nEvt, 0);

    TH1D *hTmp = new TH1D("hTmp_replica", "", nBins, eMin, eMax);
    hTmp->SetDirectory(nullptr);

    for (int r = 0; r < nReplicas; ++r)
    {
        TRandom3 rng(987654321u + 7919u * (unsigned)r);
        hTmp->Reset();

        for (size_t i = 0; i < nEvt; ++i)
        {
            const Evt &e = evts[i];
            const double er = erelMeV(e, smearPos(e, rng), off);
            if (!std::isfinite(er))
                continue;
            const double shiftedEr = er + kErelShiftMeV;
            hTmp->Fill(shiftedEr);
            evSum[i] += shiftedEr;
            evSum2[i] += shiftedEr * shiftedEr;
            evN[i] += 1;
        }

        for (int b = 0; b <= nBins + 1; ++b)
        {
            const double c = hTmp->GetBinContent(b);
            sum[b] += c;
            sum2[b] += c * c;
        }

        if ((r + 1) % 10 == 0 || r + 1 == nReplicas)
            std::cout << "         replica " << r + 1 << " / " << nReplicas << "\r"
                      << std::flush;
    }
    std::cout << "\n";
    delete hTmp;

    // mean spectrum with the randomization (systematic) uncertainty
    TH1D *hSys = new TH1D("hErel_mean_sys",
                          ";E_{rel} [MeV];counts / bin", nBins, eMin, eMax);
    TH1D *hStat = new TH1D("hErel_mean_stat",
                           ";E_{rel} [MeV];counts / bin", nBins, eMin, eMax);
    TH1D *hRel = new TH1D("hErel_relunc",
                          ";E_{rel} [MeV];#sigma_{rand} / counts", nBins, eMin, eMax);

    for (int b = 0; b <= nBins + 1; ++b)
    {
        const double mean = sum[b] / nReplicas;
        double var = sum2[b] / nReplicas - mean * mean;
        if (var < 0.)
            var = 0.;
        const double sig = std::sqrt(var); // spread induced by the randomization

        hSys->SetBinContent(b, mean);
        hSys->SetBinError(b, sig);

        hStat->SetBinContent(b, mean);
        hStat->SetBinError(b, std::sqrt(std::max(mean, 0.)));

        if (b >= 1 && b <= nBins)
            hRel->SetBinContent(b, mean > 0. ? sig / mean : 0.);
    }

    TH1D *hTotal = (TH1D *)hSys->Clone("hErel_total_uncertainty");
    hTotal->SetTitle(";E_{x} [MeV];counts / bin");
    TH1D *hSysUnc = (TH1D *)hSys->Clone("hErel_systematic_uncertainty");
    TH1D *hStatUnc = (TH1D *)hSys->Clone("hErel_statistical_uncertainty");
    TH1D *hTotalUnc = (TH1D *)hSys->Clone("hErel_total_uncertainty_per_bin");
    for (int b = 0; b <= nBins + 1; ++b)
    {
        const double sys = hSys->GetBinError(b);
        const double stat = hStat->GetBinError(b);
        const double total = std::sqrt(sys * sys + stat * stat);
        hTotal->SetBinError(b, total);
        hSysUnc->SetBinContent(b, sys);
        hStatUnc->SetBinContent(b, stat);
        hTotalUnc->SetBinContent(b, total);
        hSysUnc->SetBinError(b, 0.);
        hStatUnc->SetBinError(b, 0.);
        hTotalUnc->SetBinError(b, 0.);
    }

    // per-event mean Erel and per-event sigma
    TH1D *hEvtMean = new TH1D("hErel_evtmean",
                              ";#LTE_{rel}#GT_{rand} [MeV];counts / bin",
                              nBins, eMin, eMax);
    TH1D *hEvtSig = new TH1D("hErel_evtsigma",
                             ";#sigma_{rand}(E_{rel}) per event [MeV];events",
                             120, 0., 3.);
    double meanSigma = 0.;
    Long64_t nSigma = 0;
    for (size_t i = 0; i < nEvt; ++i)
    {
        if (evN[i] < 2)
            continue;
        const double m = evSum[i] / evN[i];
        double v = evSum2[i] / evN[i] - m * m;
        if (v < 0.)
            v = 0.;
        const double s = std::sqrt(v);
        hEvtMean->Fill(m);
        hEvtSig->Fill(s);
        meanSigma += s;
        ++nSigma;
    }
    if (nSigma)
        meanSigma /= nSigma;

    std::cout << "[part 2] mean per-event sigma(Erel) from randomization : "
              << Form("%.4f MeV", meanSigma) << "\n"
              << "[part 2] most probable per-event sigma                 : "
              << Form("%.4f MeV\n\n",
                      hEvtSig->GetEntries() > 0
                          ? hEvtSig->GetXaxis()->GetBinCenter(
                                hEvtSig->GetMaximumBin())
                          : 0.);

    // ---- draw ----
    TCanvas *c2 = new TCanvas("c2_erel_uncertainty",
                              "Erel with randomization uncertainty", 950, 800);
    TPad *q1 = new TPad("q1", "", 0., 0.32, 1., 1.);
    TPad *q2 = new TPad("q2", "", 0., 0.00, 1., 0.32);
    q1->SetBottomMargin(0.02);
    q1->SetLeftMargin(0.12);
    q2->SetTopMargin(0.03);
    q2->SetBottomMargin(0.32);
    q2->SetLeftMargin(0.12);
    q2->SetGridy();
    q1->Draw();
    q2->Draw();

    q1->cd();
    TH1D *hBand = (TH1D *)hSys->Clone("hErel_band");
    hBand->SetFillColorAlpha(kAzure + 1, 0.45);
    hBand->SetFillStyle(1001);
    hBand->SetLineColor(kAzure + 2);
    hBand->SetMarkerStyle(0);
    hBand->SetStats(kFALSE);
    hBand->SetTitle("");
    hBand->GetYaxis()->SetTitleOffset(1.25);
    hBand->GetXaxis()->SetLabelSize(0.);
    hBand->SetMaximum(1.25 * hBand->GetMaximum());
    hBand->Draw("E2");

    TH1D *hPts = (TH1D *)hStat->Clone("hErel_points");
    hPts->SetLineColor(kAzure + 2);
    hPts->SetMarkerColor(kAzure + 2);
    hPts->SetMarkerStyle(20);
    hPts->SetMarkerSize(0.7);
    hPts->SetStats(kFALSE);
    hPts->Draw("E1 SAME");

    TH1D *hNoRand2 = (TH1D *)hNoRand->Clone("hErel_norand_c2");
    styleHist(hNoRand2, kBlack, 2, 2);
    hNoRand2->Draw("HIST SAME");

    TLegend *leg2 = new TLegend(0.50, 0.66, 0.90, 0.89);
    leg2->SetBorderSize(0);
    leg2->SetFillStyle(0);
    leg2->AddEntry(hPts, Form("#LTE_{rel}#GT over %d randomizations", nReplicas), "lep");
    leg2->AddEntry(hBand, "randomization uncertainty (RMS)", "f");
    leg2->AddEntry(hNoRand2, "no randomization", "l");
    leg2->Draw();

    TLatex tl2;
    tl2.SetNDC();
    tl2.SetTextSize(0.038);
    tl2.DrawLatex(0.15, 0.85, Form("#LT#sigma_{rand}(E_{rel})#GT_{event} = %.3f MeV", meanSigma));

    q2->cd();
    hRel->SetLineColor(kAzure + 2);
    hRel->SetLineWidth(2);
    hRel->SetStats(kFALSE);
    hRel->SetTitle("");
    hRel->GetYaxis()->SetTitle("#sigma_{rand}/counts");
    hRel->GetYaxis()->SetNdivisions(505);
    hRel->GetYaxis()->SetTitleSize(0.11);
    hRel->GetYaxis()->SetTitleOffset(0.45);
    hRel->GetYaxis()->SetLabelSize(0.09);
    hRel->GetXaxis()->SetTitle("E_{rel} [MeV]");
    hRel->GetXaxis()->SetTitleSize(0.13);
    hRel->GetXaxis()->SetTitleOffset(1.0);
    hRel->GetXaxis()->SetLabelSize(0.10);
    hRel->SetMinimum(0.);
    hRel->Draw("HIST");
    c2->cd();
    c2->Update();

    TCanvas *c4 = new TCanvas("c4_total_uncertainty",
                              "Excitation energy with total uncertainty", 950, 700);
    c4->Divide(1, 2);

    c4->cd(1);
    hTotal->SetLineColor(kAzure + 2);
    hTotal->SetLineWidth(2);
    hTotal->SetMarkerColor(kAzure + 2);
    hTotal->SetMarkerStyle(20);
    hTotal->SetMarkerSize(0.7);
    hTotal->SetStats(kFALSE);
    hTotal->SetTitle("");
    hTotal->GetYaxis()->SetTitleOffset(1.25);
    hTotal->Draw("E1");

    c4->cd(2);
    hStat->SetLineColor(kAzure + 2);
    hStat->SetLineWidth(2);
    hStat->SetMarkerColor(kAzure + 2);
    hStat->SetMarkerStyle(20);
    hStat->SetMarkerSize(0.7);
    hStat->SetStats(kFALSE);
    hStat->SetTitle("Poisson statistical uncertainty only;E_{x} [MeV];counts / bin");
    hStat->Draw("E1");
    c4->Update();

    // per-event quantities
    TCanvas *c3 = new TCanvas("c3_perevent", "Per-event randomization spread",
                              1100, 480);
    c3->Divide(2, 1);

    c3->cd(1);
    styleHist(hEvtMean, kOrange + 7, 1, 2);
    hEvtMean->GetYaxis()->SetTitleOffset(1.35);
    hEvtMean->Draw("HIST");
    TH1D *hNoRand3 = (TH1D *)hNoRand->Clone("hErel_norand_c3");
    styleHist(hNoRand3, kBlack, 2, 2);
    hNoRand3->Draw("HIST SAME");
    TLegend *leg3 = new TLegend(0.45, 0.72, 0.89, 0.88);
    leg3->SetBorderSize(0);
    leg3->SetFillStyle(0);
    leg3->AddEntry(hEvtMean, "event-by-event #LTE_{rel}#GT", "l");
    leg3->AddEntry(hNoRand3, "no randomization", "l");
    leg3->Draw();

    c3->cd(2);
    styleHist(hEvtSig, kViolet + 1, 1, 2);
    hEvtSig->GetYaxis()->SetTitleOffset(1.35);
    hEvtSig->Draw("HIST");
    c3->Update();

    // -----------------------------------------------------------------------
    // 4) Save
    // -----------------------------------------------------------------------
    const TString rootOut = TString(outBase) + ".root";
    TFile *fout = new TFile(rootOut, "RECREATE");
    hNoRand->Write();
    for (int k = 0; k < kShow; ++k)
        hRand[k]->Write();
    hSys->Write();
    hStat->Write();
    hRel->Write();
    hTotal->Write();
    hSysUnc->Write();
    hStatUnc->Write();
    hTotalUnc->Write();
    hEvtMean->Write();
    hEvtSig->Write();
    c1->Write();
    c2->Write();
    c3->Write();
    c4->Write();
    fout->Close();

    c1->SaveAs(TString(outBase) + "_randomizations.pdf");
    c1->SaveAs(TString(outBase) + "_randomizations.png");
    c2->SaveAs(TString(outBase) + "_uncertainty.pdf");
    c2->SaveAs(TString(outBase) + "_uncertainty.png");
    c3->SaveAs(TString(outBase) + "_perevent.pdf");
    c3->SaveAs(TString(outBase) + "_perevent.png");
    c4->SaveAs(TString(outBase) + "_total_uncertainty.pdf");
    c4->SaveAs(TString(outBase) + "_total_uncertainty.png");

    std::cout << "[output] histograms and canvases saved in " << rootOut << "\n";
}