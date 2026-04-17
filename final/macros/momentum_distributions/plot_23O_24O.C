#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TKey.h>
#include <TStyle.h>
#include <iostream>
#include <string>

// ---------------------------------------------------------------
TTree *GetTreeFromFile(TFile *f)
{
    if (!f || f->IsZombie())
        return nullptr;

    TTree *t = nullptr;

    f->GetObject("KinTree", t);
    if (t)
        return t;

    TIter next(f->GetListOfKeys());
    TKey *k;
    while ((k = (TKey *)next()))
    {
        if (std::string(k->GetClassName()) == "TTree")
            return (TTree *)k->ReadObj();
    }
    return nullptr;
}

// ---------------------------------------------------------------
void Normalize(TH1D *h)
{
    if (!h)
        return;
    double I = h->Integral();
    if (I > 0)
        h->Scale(1.0 / I);
}

// ---------------------------------------------------------------
void plot_23O_24O()
{
    gStyle->SetOptStat(0);
    gStyle->SetCanvasPreferGL();

    const int colFrag = kOrange + 1; // 24O
    const int col1n = kCyan + 1;     // 23O+n

    std::string repopath = getenv("repopath") ? getenv("repopath") : ".";

    std::string file23 = repopath + "/results/final/23O_analyzed_test.root";
    std::string file24 = repopath + "/results/final/24O_analyzed_test.root";

    TFile *f23 = TFile::Open(file23.c_str(), "READ");
    TFile *f24 = TFile::Open(file24.c_str(), "READ");
    if (!f23 || f23->IsZombie() || !f24 || f24->IsZombie())
    {
        std::cerr << "[plot_23O_24O] Error loading files:\n"
                  << "  " << file23 << "\n"
                  << "  " << file24 << "\n";
        return;
    }

    TTree *t23 = GetTreeFromFile(f23);
    TTree *t24 = GetTreeFromFile(f24);
    if (!t23 || !t24)
    {
        std::cerr << "[plot_23O_24O] Error: missing TTrees (expected KinTree)\n";
        return;
    }

    // (p,2p) in CALIFA
    TString cond = "califa_opa > 1.2 && califa_opa < 1.8";

    TH1D *hErel23 = new TH1D("hErel23",
                             "E_{rel} (^{23}O+n);E_{rel} [MeV];A.U.",
                             80, 0, 10);

    t23->Draw("Erel*1000>>hErel23", cond, "goff");

    hErel23->SetLineColor(col1n);
    hErel23->SetFillColorAlpha(col1n, 0.3);
    hErel23->SetLineWidth(2);

    TCanvas *c1 = new TCanvas("c1", "Erel 23O+n", 800, 600);
    hErel23->Draw("HIST");

    // ===========================================================
    // 1) px, py in rotated
    // ===========================================================
    TH1D *hPx23 = new TH1D("hPx23", "^{23}O+n (CM rot);p_{x} [MeV/c];A.U.",
                           100, -0.5, 0.500);
    TH1D *hPy23 = new TH1D("hPy23", "^{23}O+n (CM rot);p_{y} [MeV/c];A.U.",
                           100, -0.500, 0.500);

    t23->Draw("px_rf_rot>>hPx23", cond, "goff");
    t23->Draw("py_rf_rot>>hPy23", cond, "goff");

    Normalize(hPx23);
    Normalize(hPy23);

    hPx23->SetLineColor(col1n);
    hPx23->SetFillColorAlpha(col1n, 0.3);
    hPx23->SetLineWidth(2);

    hPy23->SetLineColor(col1n);
    hPy23->SetFillColorAlpha(col1n, 0.3);
    hPy23->SetLineWidth(2);

    TCanvas *c2 = new TCanvas("c2", "px, py 23O+n (CM rot)", 800, 600);
    hPx23->Draw("HIST");
    hPy23->Draw("HIST SAME");

    TLegend *leg2 = new TLegend(0.6, 0.75, 0.88, 0.88);
    leg2->AddEntry(hPx23, "p_{x}  (^{23}O+n)", "lf");
    leg2->AddEntry(hPy23, "p_{y}  (^{23}O+n)", "lf");
    leg2->Draw();

    // ===========================================================
    // 3) px, py en CM rotado: 24O (fragment-only)
    // ===========================================================
    TH1D *hPx24 = new TH1D("hPx24", "^{24}O (CM rot);p_{x} [MeV/c];A.U.",
                           100, -0.500, 0.500);
    TH1D *hPy24 = new TH1D("hPy24", "^{24}O (CM rot);p_{y} [MeV/c];A.U.",
                           100, -0.500, 0.500);

    t24->Draw("px_rf_rot>>hPx24", cond, "goff");
    t24->Draw("py_rf_rot>>hPy24", cond, "goff");

    Normalize(hPx24);
    Normalize(hPy24);

    hPx24->SetLineColor(colFrag);
    hPx24->SetFillColorAlpha(colFrag, 0.3);
    hPx24->SetLineWidth(2);

    hPy24->SetLineColor(colFrag);
    hPy24->SetFillColorAlpha(colFrag, 0.3);
    hPy24->SetLineWidth(2);

    TCanvas *c3 = new TCanvas("c3", "px, py 24O (CM rot)", 800, 600);
    hPx24->Draw("HIST");
    hPy24->Draw("HIST SAME");

    TLegend *leg3 = new TLegend(0.6, 0.75, 0.88, 0.88);
    leg3->AddEntry(hPx24, "p_{x} (^{24}O)", "lf");
    leg3->AddEntry(hPy24, "p_{y} (^{24}O)", "lf");
    leg3->Draw();

    // ===========================================================
    // 4) Comparación px
    // ===========================================================
    TH1D *hPx24_v = new TH1D("hPx24_v", ";p_{x} [MeV/c];A.U.",
                             35, -0.500 * 1000, 0.500 * 1000);
    TH1D *hPx23_v = new TH1D("hPx23_v", ";p_{x} [MeV/c];A.U.",
                             35, -0.500 * 1000, 0.500 * 1000);

    t24->Draw("px_rf_rot*1000>>hPx24_v", cond, "goff");
    t23->Draw("px_rf_rot*1000>>hPx23_v", cond, "goff");

    Normalize(hPx24_v);
    Normalize(hPx23_v);

    hPx24_v->SetLineColor(colFrag);
    hPx24_v->SetFillColorAlpha(colFrag, 0.3);
    hPx24_v->SetLineWidth(2);

    hPx23_v->SetLineColor(col1n);
    hPx23_v->SetFillColorAlpha(col1n, 0.3);
    hPx23_v->SetLineWidth(2);

    TCanvas *c4 = new TCanvas("c4", "px comparison (CM rot)", 800, 600);
    hPx24_v->Draw("HIST");
    hPx23_v->Draw("HIST SAME");

    TLegend *leg4 = new TLegend(0.55, 0.75, 0.88, 0.88);
    leg4->AddEntry(hPx24_v, "^{24}O", "lf");
    leg4->AddEntry(hPx23_v, "^{23}O + n", "lf");
    leg4->Draw();

    // ===========================================================
    // 5) Comparación py
    // ===========================================================
    TH1D *hPy24_v = new TH1D("hPy24_v", ";p_{y} [MeV/c];A.U.",
                             35, -0.500 * 1000, 0.500 * 1000);
    TH1D *hPy23_v = new TH1D("hPy23_v", "p_{y} (CM rot);p_{y} [MeV/c];A.U.",
                             35, -0.500 * 1000, 0.500 * 1000);

    t24->Draw("py_rf_rot*1000>>hPy24_v", cond, "goff");
    t23->Draw("py_rf_rot*1000>>hPy23_v", cond, "goff");

    Normalize(hPy24_v);
    Normalize(hPy23_v);

    hPy24_v->SetLineColor(colFrag);
    hPy24_v->SetFillColorAlpha(colFrag, 0.3);
    hPy24_v->SetLineWidth(2);

    hPy23_v->SetLineColor(col1n);
    hPy23_v->SetFillColorAlpha(col1n, 0.3);
    hPy23_v->SetLineWidth(2);

    TCanvas *c5 = new TCanvas("c5", "py comparison (CM rot)", 800, 600);
    hPy24_v->Draw("HIST");
    hPy23_v->Draw("HIST SAME");

    TLegend *leg5 = new TLegend(0.55, 0.75, 0.88, 0.88);
    leg5->AddEntry(hPy24_v, "^{24}O", "lf");
    leg5->AddEntry(hPy23_v, "^{23}O + n", "lf");
    leg5->Draw();

    // ===========================================================
    // 6) Comparación pT (YA guardado como pT)
    // ===========================================================
    TH1D *hPt24_v = new TH1D("hPt24_v", "p_{T} = #sqrt{p_{x}^{2}+p_{y}^{2}};p_{T} [MeV/c];A.U.",
                             40, 0, 0.400 * 1000);
    TH1D *hPt23_v = new TH1D("hPt23_v", "sqrt{p_x^2+p_y^2};p_{T} [MeV/c];A.U.",
                             40, 0, 0.400 * 1000);

    t24->Draw("pT * 1000>>hPt24_v", cond, "goff");
    t23->Draw("pT * 1000>>hPt23_v", cond, "goff");

    Normalize(hPt24_v);
    Normalize(hPt23_v);

    hPt24_v->SetLineColor(colFrag);
    hPt24_v->SetFillColorAlpha(colFrag, 0.3);
    hPt24_v->SetLineWidth(2);

    hPt23_v->SetLineColor(col1n);
    hPt23_v->SetFillColorAlpha(col1n, 0.3);
    hPt23_v->SetLineWidth(2);

    TCanvas *c6 = new TCanvas("c6", "pT comparison (CM rot)", 800, 600);
    hPt24_v->Draw("HIST");
    hPt23_v->Draw("HIST SAME");

    TLegend *leg6 = new TLegend(0.55, 0.75, 0.88, 0.88);
    leg6->AddEntry(hPt24_v, "^{24}O", "lf");
    leg6->AddEntry(hPt23_v, "^{23}O + n", "lf");
    leg6->Draw();

    // ===========================================================
    // 7) Comparación pL (YA guardado como pL)
    // ===========================================================
    TH1D *hPl24_v = new TH1D("hPl24_v", "p_{L} (CM );p_{L} [MeV/c];A.U.",
                             50, -1.000, 1.000);
    TH1D *hPl23_v = new TH1D("hPl23_v", "p_{L} (CM );p_{L} [MeV/c];A.U.",
                             50, -1.000, 1.000);

    t24->Draw("pz_sys_CM>>hPl24_v", cond, "goff");
    t23->Draw("pz_sys_CM>>hPl23_v", cond, "goff");

    Normalize(hPl24_v);
    Normalize(hPl23_v);

    hPl24_v->SetLineColor(colFrag);
    hPl24_v->SetFillColorAlpha(colFrag, 0.3);
    hPl24_v->SetLineWidth(2);

    hPl23_v->SetLineColor(col1n);
    hPl23_v->SetFillColorAlpha(col1n, 0.3);
    hPl23_v->SetLineWidth(2);

    TCanvas *c7 = new TCanvas("c7", "pL comparison (CM rot)", 800, 600);
    hPl24_v->Draw("HIST");
    hPl23_v->Draw("HIST SAME");

    TLegend *leg7 = new TLegend(0.55, 0.75, 0.88, 0.88);
    leg7->AddEntry(hPl24_v, "^{24}O", "lf");
    leg7->AddEntry(hPl23_v, "^{23}O + n", "lf");
    leg7->Draw();

    // ===========================================================
    // 8) Correlaciones 2D (px vs py) con corte en |pL| (antes era pz_cm_*)
    // ===========================================================
    TH2D *hCorr24 = new TH2D("hCorr24",
                             "^{24}O (CM rot);p_{x} [MeV/c];p_{y} [MeV/c]",
                             50, -0.500, 0.500, 50, -0.500, 0.500);

    t24->Draw("py_rf_rot:px_rf_rot>>hCorr24",
              cond, "goff");

    TCanvas *c8 = new TCanvas("c8", "Corr 24O (CM rot)", 800, 600);
    hCorr24->Draw("COLZ");

    TH2D *hCorr23 = new TH2D("hCorr23",
                             "^{23}O+n (CM rot);p_{x} [MeV/c];p_{y} [MeV/c]",
                             50, -0.500, 0.500, 50, -0.500, 0.500);

    t23->Draw("py_rf_rot:px_rf_rot>>hCorr23",
              cond, "goff");

    TCanvas *c9 = new TCanvas("c9", "Corr 23O+n (CM rot)", 800, 600);
    hCorr23->Draw("COLZ");

    // ===========================================================
    // 9) pT por regiones de Erel (solo 23O+n)
    // ===========================================================
    TH1D *hPt23_p1 = new TH1D("hPt23_p1", "E_{rel}<1.8;p_{T} [MeV/c];A.U.",
                              50, 0, 0.300);
    TH1D *hPt23_p2 = new TH1D("hPt23_p2", "2.1<E_{rel}<4;p_{T} [MeV/c];A.U.",
                              50, 0, 0.300);
    TH1D *hPt23_p3 = new TH1D("hPt23_p3", "5<E_{rel}<7;p_{T} [MeV/c];A.U.",
                              50, 0, 0.300);

    t23->Draw("pT>>hPt23_p1", cond + " && Erel*1000 < 1.8", "goff");
    t23->Draw("pT>>hPt23_p2", cond + " && Erel*1000 > 2.1 && Erel*1000 < 4.0", "goff");
    t23->Draw("pT>>hPt23_p3", cond + " && Erel*1000 > 5.0 && Erel*1000 < 7.0", "goff");

    Normalize(hPt23_p1);
    Normalize(hPt23_p2);
    Normalize(hPt23_p3);

    hPt23_p1->SetLineColor(kCyan + 1);
    hPt23_p1->SetFillColorAlpha(kCyan + 1, 0.3);
    hPt23_p1->SetLineWidth(2);

    hPt23_p2->SetLineColor(kOrange + 1);
    hPt23_p2->SetFillColorAlpha(kOrange + 1, 0.3);
    hPt23_p2->SetLineWidth(2);

    hPt23_p3->SetLineColor(kGreen + 1);
    hPt23_p3->SetFillColorAlpha(kGreen + 1, 0.3);
    hPt23_p3->SetLineWidth(2);

    TCanvas *c10 = new TCanvas("c10", "pT peaks (23O+n)", 800, 600);
    hPt23_p1->Draw("HIST");
    hPt23_p2->Draw("HIST SAME");
    hPt23_p3->Draw("HIST SAME");

    TLegend *leg10 = new TLegend(0.55, 0.72, 0.88, 0.88);
    leg10->AddEntry(hPt23_p1, "E_{rel}<1.8", "lf");
    leg10->AddEntry(hPt23_p2, "2.1<E_{rel}<4", "lf");
    leg10->AddEntry(hPt23_p3, "5<E_{rel}<7", "lf");
    leg10->Draw();

    // ===========================================================
    // 10) Comparaciones por picos de Erel:
    //     - px_rf_rot: 24O (sin corte) vs 23O+n (corte en el pico)
    //     - pT:        24O (sin corte) vs 23O+n (corte en el pico)
    // ===========================================================

    TString cut_p1 = cond + " && Erel*1000 < 1.8";
    TString cut_p2 = cond + " && Erel*1000 > 2.1 && Erel*1000 < 4.0";
    TString cut_p3 = cond + " && Erel*1000 > 5.0 && Erel*1000 < 7.0";

    // ---------------------------
    // A) px_rf_rot (3 canvases)
    // ---------------------------
    const int nBinsPx = 40;
    const double pxMin = -0.5;
    const double pxMax = 0.5;

    // Pico 1
    TH1D *hPx24_p1 = new TH1D("hPx24_p1", "p_{x}; p_{x} [GeV/c]; A.U.", nBinsPx, pxMin, pxMax);
    TH1D *hPx23_p1 = new TH1D("hPx23_p1", "p_{x}; p_{x} [GeV/c]; A.U.", nBinsPx, pxMin, pxMax);

    t24->Draw("px_rf_rot>>hPx24_p1", cond, "goff");   // 24O: sin corte Erel
    t23->Draw("px_rf_rot>>hPx23_p1", cut_p1, "goff"); // 23O+n: corte pico 1

    Normalize(hPx24_p1);
    Normalize(hPx23_p1);

    hPx24_p1->SetLineColor(colFrag);
    hPx24_p1->SetFillColorAlpha(colFrag, 0.3);
    hPx24_p1->SetLineWidth(2);

    hPx23_p1->SetLineColor(col1n);
    hPx23_p1->SetFillColorAlpha(col1n, 0.3);
    hPx23_p1->SetLineWidth(2);

    TCanvas *c11 = new TCanvas("c11", "px_rf_rot peak 1", 800, 600);
    hPx24_p1->Draw("HIST");
    hPx23_p1->Draw("HIST SAME");

    TLegend *leg11 = new TLegend(0.55, 0.75, 0.88, 0.88);
    leg11->AddEntry(hPx24_p1, "^{24}O (all)", "lf");
    leg11->AddEntry(hPx23_p1, "^{23}O+n (E_{rel} peak 1)", "lf");
    leg11->Draw();

    // Pico 2
    TH1D *hPx24_p2 = new TH1D("hPx24_p2", "; p_{x} [GeV/c]; A.U.", nBinsPx, pxMin, pxMax);
    TH1D *hPx23_p2 = new TH1D("hPx23_p2", "; p_{x} [GeV/c]; A.U.", nBinsPx, pxMin, pxMax);

    t24->Draw("px_rf_rot>>hPx24_p2", cond, "goff");
    t23->Draw("px_rf_rot>>hPx23_p2", cut_p2, "goff");

    Normalize(hPx24_p2);
    Normalize(hPx23_p2);

    hPx24_p2->SetLineColor(colFrag);
    hPx24_p2->SetFillColorAlpha(colFrag, 0.3);
    hPx24_p2->SetLineWidth(2);

    hPx23_p2->SetLineColor(col1n);
    hPx23_p2->SetFillColorAlpha(col1n, 0.3);
    hPx23_p2->SetLineWidth(2);

    TCanvas *c12 = new TCanvas("c12", "px_rf_rot peak 2", 800, 600);
    hPx24_p2->Draw("HIST");
    hPx23_p2->Draw("HIST SAME");

    TLegend *leg12 = new TLegend(0.55, 0.75, 0.88, 0.88);
    leg12->AddEntry(hPx24_p2, "^{24}O (all)", "lf");
    leg12->AddEntry(hPx23_p2, "^{23}O+n (E_{rel} peak 2)", "lf");
    leg12->Draw();

    // Pico 3
    TH1D *hPx24_p3 = new TH1D("hPx24_p3", "; p_{x} [GeV/c]; A.U.", nBinsPx, pxMin, pxMax);
    TH1D *hPx23_p3 = new TH1D("hPx23_p3", "; p_{x} [GeV/c]; A.U.", nBinsPx, pxMin, pxMax);

    t24->Draw("px_rf_rot>>hPx24_p3", cond, "goff");
    t23->Draw("px_rf_rot>>hPx23_p3", cut_p3, "goff");

    Normalize(hPx24_p3);
    Normalize(hPx23_p3);

    hPx24_p3->SetLineColor(colFrag);
    hPx24_p3->SetFillColorAlpha(colFrag, 0.3);
    hPx24_p3->SetLineWidth(2);

    hPx23_p3->SetLineColor(col1n);
    hPx23_p3->SetFillColorAlpha(col1n, 0.3);
    hPx23_p3->SetLineWidth(2);

    TCanvas *c13 = new TCanvas("c13", "px_rf_rot peak 3", 800, 600);
    hPx24_p3->Draw("HIST");
    hPx23_p3->Draw("HIST SAME");

    TLegend *leg13 = new TLegend(0.55, 0.75, 0.88, 0.88);
    leg13->AddEntry(hPx24_p3, "^{24}O (all)", "lf");
    leg13->AddEntry(hPx23_p3, "^{23}O+n (E_{rel} peak 3)", "lf");
    leg13->Draw();

    // ---------------------------
    // B) pT (3 canvases)
    // ---------------------------
    const int nBinsPt = 40;
    const double ptMin = 0.0;
    const double ptMax = 400.0;

    // Pico 1
    TH1D *hPt24_comp_p1 = new TH1D("hPt24_comp_p1", "p_{T} = #sqrt{p_{x}^{2}+p_{y}^{2}}; p_{T} [MeV/c]; A.U.", 20, ptMin, ptMax);
    TH1D *hPt23_comp_p1 = new TH1D("hPt23_comp_p1", "p_{T} = #sqrt{p_{x}^{2}+p_{y}^{2}}; p_{T} [MeV/c]; A.U.", 20, ptMin, ptMax);

    t24->Draw("pT*1000>>hPt24_comp_p1", cond, "goff");
    t23->Draw("pT*1000>>hPt23_comp_p1", cut_p1, "goff");

    Normalize(hPt24_comp_p1);
    Normalize(hPt23_comp_p1);

    hPt24_comp_p1->SetLineColor(colFrag);
    hPt24_comp_p1->SetFillColorAlpha(colFrag, 0.3);
    hPt24_comp_p1->SetLineWidth(2);

    hPt23_comp_p1->SetLineColor(col1n);
    hPt23_comp_p1->SetFillColorAlpha(col1n, 0.3);
    hPt23_comp_p1->SetLineWidth(2);

    TCanvas *c14 = new TCanvas("c14", "pT peak 1", 800, 600);
    c14->SetFillColor(TColor::GetColor(249, 249, 249));
    hPt24_comp_p1->Draw("HIST");
    hPt23_comp_p1->Draw("HIST SAME");
    hPt24_comp_p1->Draw("E1 SAME");
    hPt23_comp_p1->Draw("E1 SAME");

    TLegend *leg14 = new TLegend(0.55, 0.75, 0.88, 0.88);
    leg14->AddEntry(hPt24_comp_p1, "^{24}O (all)", "lf");
    leg14->AddEntry(hPt23_comp_p1, "^{23}O+n (E_{rel} peak 1)", "lf");
    leg14->Draw();

    // Pico 2
    TH1D *hPt24_comp_p2 = new TH1D("hPt24_comp_p2", "p_{T} = #sqrt{p_{x}^{2}+p_{y}^{2}}; p_{T} [MeV/c]; A.U.", nBinsPt, ptMin, ptMax);
    TH1D *hPt23_comp_p2 = new TH1D("hPt23_comp_p2", "p_{T} = #sqrt{p_{x}^{2}+p_{y}^{2}}; p_{T} [MeV/c]; A.U.", nBinsPt, ptMin, ptMax);

    t24->Draw("pT*1000>>hPt24_comp_p2", cond, "goff");
    t23->Draw("pT*1000>>hPt23_comp_p2", cut_p2, "goff");

    Normalize(hPt24_comp_p2);
    Normalize(hPt23_comp_p2);

    hPt24_comp_p2->SetLineColor(colFrag);
    hPt24_comp_p2->SetFillColorAlpha(colFrag, 0.3);
    hPt24_comp_p2->SetLineWidth(2);

    hPt23_comp_p2->SetLineColor(col1n);
    hPt23_comp_p2->SetFillColorAlpha(col1n, 0.3);
    hPt23_comp_p2->SetLineWidth(2);

    TCanvas *c15 = new TCanvas("c15", "pT peak 2", 800, 600);
    c15->SetFillColor(TColor::GetColor(249, 249, 249));
    hPt24_comp_p2->Draw("HIST");
    hPt23_comp_p2->Draw("HIST SAME");
    hPt24_comp_p2->Draw("E1 SAME");
    hPt23_comp_p2->Draw("E1 SAME");

    TLegend *leg15 = new TLegend(0.55, 0.75, 0.88, 0.88);
    leg15->AddEntry(hPt24_comp_p2, "^{24}O (all)", "lf");
    leg15->AddEntry(hPt23_comp_p2, "^{23}O+n (E_{rel} peak 2)", "lf");
    leg15->Draw();

    // Pico 3
    TH1D *hPt24_comp_p3 = new TH1D("hPt24_comp_p3", "p_{T} = #sqrt{p_{x}^{2}+p_{y}^{2}}; p_{T} [MeV/c]; A.U.", nBinsPt, ptMin, ptMax);
    TH1D *hPt23_comp_p3 = new TH1D("hPt23_comp_p3", "p_{T} = #sqrt{p_{x}^{2}+p_{y}^{2}}; p_{T} [MeV/c]; A.U.", nBinsPt, ptMin, ptMax);

    t24->Draw("pT*1000>>hPt24_comp_p3", cond, "goff");
    t23->Draw("pT*1000>>hPt23_comp_p3", cut_p3, "goff");

    Normalize(hPt24_comp_p3);
    Normalize(hPt23_comp_p3);

    hPt24_comp_p3->SetLineColor(colFrag);
    hPt24_comp_p3->SetFillColorAlpha(colFrag, 0.3);
    hPt24_comp_p3->SetLineWidth(2);

    hPt23_comp_p3->SetLineColor(col1n);
    hPt23_comp_p3->SetFillColorAlpha(col1n, 0.3);
    hPt23_comp_p3->SetLineWidth(2);

    TCanvas *c16 = new TCanvas("c16", "pT peak 3", 800, 600);
    c16->SetFillColor(TColor::GetColor(249, 249, 249));
    hPt24_comp_p3->Draw("HIST");
    hPt23_comp_p3->Draw("HIST SAME");
    hPt24_comp_p3->Draw("E1 SAME");
    hPt23_comp_p3->Draw("E1 SAME");

    TLegend *leg16 = new TLegend(0.55, 0.75, 0.88, 0.88);
    leg16->AddEntry(hPt24_comp_p3, "^{24}O (all)", "lf");
    leg16->AddEntry(hPt23_comp_p3, "^{23}O+n (E_{rel} peak 3)", "lf");
    leg16->Draw();

    std::cout << "[plot_23O_24O] DONE (KinTree)" << std::endl;
}