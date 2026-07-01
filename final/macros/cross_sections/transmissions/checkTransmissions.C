void checkTransmissions()
{
    gStyle->SetOptStat(0);
    gStyle->SetCanvasPreferGL();
    Int_t ciBg = TColor::GetColor("#f9f9f9");

    const char *file24 = "/nucl_lustre/pablogrusell/g249/g249_analysis/results/final/data_24O.root";
    const char *file25 = "/nucl_lustre/pablogrusell/g249/g249_analysis/results/final/data_25F.root";

    // Unreacted cuts
    auto *unreacted_25f_cut = new TCutG("unreacted_25f_cut", 23);
    unreacted_25f_cut->SetPoint(0, 2.74329, 8.48631);
    unreacted_25f_cut->SetPoint(1, 2.75374, 8.55453);
    unreacted_25f_cut->SetPoint(2, 2.7594, 8.66688);
    unreacted_25f_cut->SetPoint(3, 2.76376, 8.80331);
    unreacted_25f_cut->SetPoint(4, 2.76637, 8.92368);
    unreacted_25f_cut->SetPoint(5, 2.7655, 9.04005);
    unreacted_25f_cut->SetPoint(6, 2.76289, 9.15642);
    unreacted_25f_cut->SetPoint(7, 2.75853, 9.26877);
    unreacted_25f_cut->SetPoint(8, 2.75418, 9.3731);
    unreacted_25f_cut->SetPoint(9, 2.74939, 9.4373);
    unreacted_25f_cut->SetPoint(10, 2.7446, 9.48545);
    unreacted_25f_cut->SetPoint(11, 2.73937, 9.48545);
    unreacted_25f_cut->SetPoint(12, 2.73153, 9.45736);
    unreacted_25f_cut->SetPoint(13, 2.72413, 9.36507);
    unreacted_25f_cut->SetPoint(14, 2.71934, 9.23667);
    unreacted_25f_cut->SetPoint(15, 2.71542, 9.08018);
    unreacted_25f_cut->SetPoint(16, 2.71454, 8.9598);
    unreacted_25f_cut->SetPoint(17, 2.71454, 8.83139);
    unreacted_25f_cut->SetPoint(18, 2.72064, 8.65885);
    unreacted_25f_cut->SetPoint(19, 2.72979, 8.53045);
    unreacted_25f_cut->SetPoint(20, 2.7385, 8.47427);
    unreacted_25f_cut->SetPoint(21, 2.74285, 8.47026);
    unreacted_25f_cut->SetPoint(22, 2.74329, 8.48631);
    unreacted_25f_cut->SetVarX("AoQ_frag");
    unreacted_25f_cut->SetVarY("Z_frag_est");

    // 24O cuts
    auto *reacted_23o_cutg = new TCutG("reacted_23o_cutg", 13);
    reacted_23o_cutg->SetVarX("AoQ_frag");
    reacted_23o_cutg->SetVarY("Z_frag_est");
    reacted_23o_cutg->SetTitle("Graph");
    reacted_23o_cutg->SetFillStyle(1000);
    reacted_23o_cutg->SetPoint(0, 2.78967, 8.38836);
    reacted_23o_cutg->SetPoint(1, 2.78589, 8.14281);
    reacted_23o_cutg->SetPoint(2, 2.794, 7.85634);
    reacted_23o_cutg->SetPoint(3, 2.81074, 7.68753);
    reacted_23o_cutg->SetPoint(4, 2.85233, 7.67218);
    reacted_23o_cutg->SetPoint(5, 2.87663, 7.83588);
    reacted_23o_cutg->SetPoint(6, 2.87933, 8.11212);
    reacted_23o_cutg->SetPoint(7, 2.87231, 8.30139);
    reacted_23o_cutg->SetPoint(8, 2.85503, 8.5009);
    reacted_23o_cutg->SetPoint(9, 2.83612, 8.56229);
    reacted_23o_cutg->SetPoint(10, 2.80966, 8.52648);
    reacted_23o_cutg->SetPoint(11, 2.79616, 8.44463);
    reacted_23o_cutg->SetPoint(12, 2.78967, 8.38836);
    reacted_23o_cutg->SetVarX("AoQ_frag");
    reacted_23o_cutg->SetVarY("Z_frag_est");

    TFile *f24 = TFile::Open(file24);
    TFile *f25 = TFile::Open(file25);
    if (!f24 || f24->IsZombie())
    {
        printf("Failed to open %s\n", file24);
        return;
    }
    if (!f25 || f25->IsZombie())
    {
        printf("Failed to open %s\n", file25);
        return;
    }

    TTree *t24 = (TTree *)f24->Get("FilterDataTree");
    TTree *t25 = (TTree *)f25->Get("FilterDataTree");
    if (!t24 || !t25)
    {
        printf("Could not find FilterDataTree in one of the files\n");
        return;
    }

    // 3sigma elliptical identification cuts (same parameters as calculateCrossSections.C)
    // 24O data: Z = 8.06, AoQ = 2.958, + opening-angle cut on CALIFA
    const TString cut24 = "califa_opa>1.25&&califa_opa<1.65 && reacted_23o_cutg";
    // 25F data: Z = 8.97, AoQ = 2.740
    const TString cut25 = "unreacted_25f_cut";

    // --- Canvas 1a: Foots charge ---
    TCanvas *cFootsCharge = new TCanvas("cFootsCharge", "Foots Charge: 24O vs 25F", 800, 700);
    cFootsCharge->SetFillColor(ciBg);
    cFootsCharge->Divide(2, 2);

    // --- Canvas 1b: Foots position ---
    TCanvas *cFootsPos = new TCanvas("cFootsPos", "Foots Position: 24O vs 25F", 800, 700);
    cFootsPos->SetFillColor(ciBg);
    cFootsPos->Divide(2, 2);

    // --- Canvas 2: Fibers 32 & 30 ---
    TCanvas *cFib3230 = new TCanvas("cFib3230", "Fibers 32 & 30: 24O vs 25F", 1000, 500);
    cFib3230->SetFillColor(ciBg);
    cFib3230->Divide(2, 1);

    // --- Canvas 3: Fibers 33 & 31 ---
    TCanvas *cFib3331 = new TCanvas("cFib3331", "Fibers 33 & 31: 24O vs 25F", 700, 500);
    cFib3331->SetFillColor(ciBg);

    // --- Canvas 4: ToFD ---
    TCanvas *cToFD = new TCanvas("cToFD", "ToFD: 24O vs 25F", 700, 500);
    cToFD->SetFillColor(ciBg);

    // Helper: apply background color to all pads of a canvas
    auto setBgColor = [&](TCanvas *cvs)
    {
        for (int i = 0; i <= 20; ++i)
        {
            TPad *p = (TPad *)cvs->GetPad(i);
            if (p)
                p->SetFillColor(ciBg);
        }
    };
    setBgColor(cFootsCharge);
    setBgColor(cFootsPos);
    setBgColor(cFib3230);
    setBgColor(cFib3331);
    setBgColor(cToFD);

    // helper lambda to draw two histograms overlayed from two trees
    auto drawOverlay = [&](TCanvas *cvs, int padN, TTree *ta, TTree *tb, const char *varA, const char *nameA, const char *title, int colorA, int colorB, int nbins = 100, double xlo = 0, double xhi = -1, const char *cutA = "", const char *cutB = "")
    {
        cvs->cd(padN);
        TString hA = TString::Format("%s_%s", nameA, "24");
        TString hB = TString::Format("%s_%s", nameA, "25");
        if (xhi > xlo)
        {
            ta->Draw(Form("%s>>%s(%d,%g,%g)", varA, hA.Data(), nbins, xlo, xhi), cutA, "goff");
            tb->Draw(Form("%s>>%s(%d,%g,%g)", varA, hB.Data(), nbins, xlo, xhi), cutB, "goff");
        }
        else
        {
            ta->Draw(Form("%s>>%s", varA, hA.Data()), cutA, "goff");
            tb->Draw(Form("%s>>%s", varA, hB.Data()), cutB, "goff");
        }
        TH1 *ha = (TH1 *)gDirectory->Get(hA.Data());
        TH1 *hb = (TH1 *)gDirectory->Get(hB.Data());
        if (!ha || !hb)
        {
            return;
        }
        // normalize to unit area for shape comparison
        double ia = ha->Integral();
        if (ia > 0)
            ha->Scale(1.0 / ia);
        double ib = hb->Integral();
        if (ib > 0)
            hb->Scale(1.0 / ib);

        // fill with alpha for visual overlap and keep line for border
        ha->SetLineColor(colorA);
        ha->SetLineWidth(2);
        ha->SetFillColorAlpha(colorA, 0.3);
        ha->SetFillStyle(1001);

        hb->SetLineColor(colorB);
        hb->SetLineWidth(2);
        hb->SetFillColorAlpha(colorB, 0.3);
        hb->SetFillStyle(1001);

        ha->SetTitle(title);
        ha->Draw("HIST");
        hb->Draw("HIST SAME");
        TLegend *leg = new TLegend(0.65, 0.75, 0.88, 0.88);
        leg->SetBorderSize(0);
        leg->AddEntry(ha, "24O", "f");
        leg->AddEntry(hb, "25F", "f");
        leg->Draw();
    };

    // Colors: 24O orange, 25F blue
    const int col24 = kOrange + 7;
    const int col25 = kCyan - 3;

    // 1D charge histograms for foots 5..8 (variables foot_z_5..foot_z_8)
    for (int f = 5; f <= 8; ++f)
    {
        TString var = TString::Format("foot_z_%d", f);
        TString title = TString::Format("Charge foot %d;%s;Counts", f, var.Data());
        drawOverlay(cFootsCharge, f - 4, t24, t25, var.Data(), Form("charge_f%d", f), title.Data(), col24, col25, 200, 0, 15, cut24.Data(), cut25.Data());
    }

    // 1D position histograms for foots 5..8 (PosFoot5..PosFoot8)
    for (int f = 5; f <= 8; ++f)
    {
        TString var = TString::Format("PosFoot%d", f);
        TString title = TString::Format("Position foot %d;%s [mm];Counts", f, var.Data());
        drawOverlay(cFootsPos, f - 4, t24, t25, var.Data(), Form("pos_f%d", f), title.Data(), col24, col25, 200, -50, 50, cut24.Data(), cut25.Data());
    }

    // fiber 32 (fib32X)
    drawOverlay(cFib3230, 1, t24, t25, "fib32X", "fib32X", "Fiber 32 X;fib32X [cm];Counts", col24, col25, 160, -80, 80, cut24.Data(), cut25.Data());

    // fiber 30 (fib30Y)
    drawOverlay(cFib3230, 2, t24, t25, "fib30Y", "fib30Y", "Fiber 30 (fib30Y);fib30Y [cm];Counts", col24, col25, 160, -80, 80, cut24.Data(), cut25.Data());

    // Overlap positions for fibers 33 and 31 using fib31X and fib33X on same pad
    cFib3331->cd();
    TString h313 = "fib31X_24";
    TString h333 = "fib33X_24";
    // create from 24O — explicit range [-80, 80]
    t24->Draw("fib31X>>fib31X_24(160,-80,80)", cut24.Data(), "goff");
    t24->Draw("fib33X>>fib33X_24(160,-80,80)", cut24.Data(), "goff");
    TH1 *h31_24 = (TH1 *)gDirectory->Get("fib31X_24");
    TH1 *h33_24 = (TH1 *)gDirectory->Get("fib33X_24");
    // create from 25F — explicit range [-80, 80]
    t25->Draw("fib31X>>fib31X_25(160,-80,80)", cut25.Data(), "goff");
    t25->Draw("fib33X>>fib33X_25(160,-80,80)", cut25.Data(), "goff");
    TH1 *h31_25 = (TH1 *)gDirectory->Get("fib31X_25");
    TH1 *h33_25 = (TH1 *)gDirectory->Get("fib33X_25");
    // Normalize fib31/fib33: each isotope independently, but both sides share the same scale
    {
        // 24O: use max integral of fib31 vs fib33 so relative heights are preserved
        double i31_24 = h31_24 ? h31_24->Integral() : 0;
        double i33_24 = h33_24 ? h33_24->Integral() : 0;
        double imax24 = std::max(i31_24, i33_24);
        if (imax24 > 0)
        {
            if (h31_24)
                h31_24->Scale(1.0 / imax24);
            if (h33_24)
                h33_24->Scale(1.0 / imax24);
        }
        // 25F: same approach independently
        double i31_25 = h31_25 ? h31_25->Integral() : 0;
        double i33_25 = h33_25 ? h33_25->Integral() : 0;
        double imax25 = std::max(i31_25, i33_25);
        if (imax25 > 0)
        {
            if (h31_25)
                h31_25->Scale(1.0 / imax25);
            if (h33_25)
                h33_25->Scale(1.0 / imax25);
        }
    }
    // Compute global Y-axis max so no histogram overflows the frame
    {
        double ymax = 0;
        if (h31_24)
            ymax = std::max(ymax, h31_24->GetMaximum());
        if (h33_24)
            ymax = std::max(ymax, h33_24->GetMaximum());
        if (h31_25)
            ymax = std::max(ymax, h31_25->GetMaximum());
        if (h33_25)
            ymax = std::max(ymax, h33_25->GetMaximum());
        ymax *= 1.1; // 10 % headroom
        if (h31_24)
            h31_24->SetMaximum(ymax);
    }
    // Style and draw — make sure we are on the right canvas
    cFib3331->cd();
    if (h31_24)
    {
        h31_24->SetLineColor(col24);
        h31_24->SetLineWidth(2);
        h31_24->SetFillColorAlpha(col24, 0.3);
        h31_24->SetFillStyle(1001);
        h31_24->SetTitle("Fiber positions: fib31X & fib33X;X position [cm];Normalized Counts");
        h31_24->Draw("HIST");
    }
    if (h33_24)
    {
        h33_24->SetLineColor(col24);
        h33_24->SetLineStyle(2);
        h33_24->SetLineWidth(2);
        h33_24->SetFillColorAlpha(col24, 0.3);
        h33_24->SetFillStyle(1001);
        h33_24->Draw("HIST SAME");
    }
    if (h31_25)
    {
        h31_25->SetLineColor(col25);
        h31_25->SetLineWidth(2);
        h31_25->SetLineStyle(3);
        h31_25->SetFillColorAlpha(col25, 0.3);
        h31_25->SetFillStyle(1001);
        h31_25->Draw("HIST SAME");
    }
    if (h33_25)
    {
        h33_25->SetLineColor(col25);
        h33_25->SetLineWidth(2);
        h33_25->SetLineStyle(4);
        h33_25->SetFillColorAlpha(col25, 0.3);
        h33_25->SetFillStyle(1001);
        h33_25->Draw("HIST SAME");
    }
    {
        TLegend *leg = new TLegend(0.55, 0.65, 0.88, 0.88);
        leg->SetBorderSize(0);
        if (h31_24)
            leg->AddEntry(h31_24, "24O fib31X", "f");
        if (h33_24)
            leg->AddEntry(h33_24, "24O fib33X", "f");
        if (h31_25)
            leg->AddEntry(h31_25, "25F fib31X", "f");
        if (h33_25)
            leg->AddEntry(h33_25, "25F fib33X", "f");
        leg->Draw();
    }

    // X position in ToFD (tofdX)
    drawOverlay(cToFD, 0, t24, t25, "tofdX", "tofdX", "ToFD X Position;tofdX [cm];Counts", col24, col25, 200, -150, 150, cut24.Data(), cut25.Data());

    // ==================== Eloss vs Position for each fiber ====================
    // 4 fibers × 2 isotopes → 4 columns × 2 rows
    TCanvas *cFibEloss = new TCanvas("cFibEloss", "Fiber Eloss vs Position: 24O (top) vs 25F (bottom)", 1600, 700);
    cFibEloss->SetFillColor(ciBg);
    cFibEloss->Divide(4, 2);
    setBgColor(cFibEloss);

    // {fiber label, position variable, eloss variable, x-axis label, xlo, xhi}
    struct FibInfo
    {
        const char *label;
        const char *posVar;
        const char *elossVar;
        const char *xLabel;
        double xlo, xhi;
    };
    FibInfo fibs[4] = {
        {"Fib30", "fib30Y", "ElossFib30", "fib30Y [cm]", -80, 80},
        {"Fib31", "fib31X", "ElossFib31", "fib31X [cm]", -80, 80},
        {"Fib33", "fib33X", "ElossFib33", "fib33X [cm]", -80, 80},
        {"Fib32", "fib32X", "ElossFib32", "fib32X [cm]", -80, 80},
    };

    for (int i = 0; i < 4; ++i)
    {
        // 24O → top row (pads 1..4)
        cFibEloss->cd(i + 1);
        TString drawExpr24 = TString::Format("%s:%s>>hFibEl24_%d(100,%g,%g,100,0,30)",
                                             fibs[i].elossVar, fibs[i].posVar,
                                             i, fibs[i].xlo, fibs[i].xhi);
        t24->Draw(drawExpr24.Data(), cut24.Data(), "COLZ");
        TH2 *hfe24 = (TH2 *)gDirectory->Get(Form("hFibEl24_%d", i));
        if (hfe24)
            hfe24->SetTitle(Form("24O %s;%s;Eloss", fibs[i].label, fibs[i].xLabel));

        // 25F → bottom row (pads 5..8)
        cFibEloss->cd(i + 5);
        TString drawExpr25 = TString::Format("%s:%s>>hFibEl25_%d(100,%g,%g,100,0,30)",
                                             fibs[i].elossVar, fibs[i].posVar,
                                             i, fibs[i].xlo, fibs[i].xhi);
        t25->Draw(drawExpr25.Data(), cut25.Data(), "COLZ");
        TH2 *hfe25 = (TH2 *)gDirectory->Get(Form("hFibEl25_%d", i));
        if (hfe25)
            hfe25->SetTitle(Form("25F %s;%s;Eloss", fibs[i].label, fibs[i].xLabel));
    }

    // ==================== Foot charge vs position ====================
    // Foots 5..8, 2 isotopes → 4 columns × 2 rows
    TCanvas *cFootZvsPos = new TCanvas("cFootZvsPos", "Foot Charge vs Position: 24O (top) vs 25F (bottom)", 1600, 700);
    cFootZvsPos->SetFillColor(ciBg);
    cFootZvsPos->Divide(4, 2);
    setBgColor(cFootZvsPos);

    for (int f = 5; f <= 8; ++f)
    {
        int col = f - 4; // pad column 1..4

        // 24O → top row
        cFootZvsPos->cd(col);
        TString expr24 = TString::Format("foot_z_%d:PosFoot%d>>hFootZP24_%d(641,-48,48,300,0,15)", f, f, f);
        t24->Draw(expr24.Data(), cut24.Data(), "COLZ");
        TH2 *hfp24 = (TH2 *)gDirectory->Get(Form("hFootZP24_%d", f));
        if (hfp24)
            hfp24->SetTitle(Form("24O Foot %d;PosFoot%d [mm];Charge", f, f));

        // 25F → bottom row
        cFootZvsPos->cd(col + 4);
        TString expr25 = TString::Format("foot_z_%d:PosFoot%d>>hFootZP25_%d(641,-48,48,300,0,15)", f, f, f);
        t25->Draw(expr25.Data(), cut25.Data(), "COLZ");
        TH2 *hfp25 = (TH2 *)gDirectory->Get(Form("hFootZP25_%d", f));
        if (hfp25)
            hfp25->SetTitle(Form("25F Foot %d;PosFoot%d [mm];Charge", f, f));
    }

    // ==================== tofdX vs Z_frag_est ====================
    TCanvas *cToFDvsZ = new TCanvas("cToFDvsZ", "tofdX vs Z_frag_est: 24O (left) vs 25F (right)", 1200, 600);
    cToFDvsZ->SetFillColor(ciBg);
    cToFDvsZ->Divide(2, 1);
    setBgColor(cToFDvsZ);

    // 24O
    cToFDvsZ->cd(1);
    t24->Draw("tofdX:Z_frag_est>>hToFDZ24(100,0,15,200,-150,150)", cut24.Data(), "COLZ");
    TH2 *hToFDZ24 = (TH2 *)gDirectory->Get("hToFDZ24");
    if (hToFDZ24)
        hToFDZ24->SetTitle("24O;Z_{Plane 1};X [cm]");

    // 25F
    cToFDvsZ->cd(2);
    t25->Draw("tofdX:Z_frag_est>>hToFDZ25(100,0,15,200,-150,150)", cut25.Data(), "COLZ");
    TH2 *hToFDZ25 = (TH2 *)gDirectory->Get("hToFDZ25");
    if (hToFDZ25)
        hToFDZ25->SetTitle("25F;Z_{Plane 1};X [cm]");

    cFootsCharge->Update();
    cFootsPos->Update();
    cFib3230->Update();
    cFib3331->Modified();
    cFib3331->Update();
    cToFD->Update();
    cFibEloss->Update();
    cFootZvsPos->Update();
    cToFDvsZ->Update();

    // leave files open for interactive inspection
}
