// erel_corrected.C
// Corrects the Erel spectrum bin-by-bin using the measured efficiency TGraph.

void erel_corrected()
{
    // =====================================================================
    //  1. Build the efficiency TGraph from measured data (values in %)
    // =====================================================================

    Double_t eff_x[] = {
        0.025, 0.075, 0.125, 0.175, 0.225, 0.275, 0.325, 0.375, 0.425, 0.475,
        0.525, 0.575, 0.625, 0.675, 0.725, 0.775, 0.825, 0.875, 0.925, 0.975,
        1.025, 1.075, 1.125, 1.175, 1.225, 1.275, 1.325, 1.375, 1.425, 1.475,
        1.525, 1.575, 1.625, 1.675, 1.725, 1.775, 1.825, 1.875, 1.925, 1.975,
        2.025, 2.075, 2.125, 2.175, 2.225, 2.275, 2.325, 2.375, 2.425, 2.475,
        2.525, 2.575, 2.625, 2.675, 2.725, 2.775, 2.825, 2.875, 2.925, 2.975,
        3.025, 3.075, 3.125, 3.175, 3.225, 3.275, 3.325, 3.375, 3.425, 3.475,
        3.525, 3.575, 3.625, 3.675, 3.725, 3.775, 3.825, 3.875, 3.925, 3.975,
        4.025, 4.075, 4.125, 4.175, 4.225, 4.275, 4.325, 4.375, 4.425, 4.475,
        4.525, 4.575, 4.625, 4.675, 4.725, 4.775, 4.825, 4.875, 4.925, 4.975,
        5.025, 5.075, 5.125, 5.175, 5.225, 5.275, 5.325, 5.375, 5.425, 5.475,
        5.525, 5.575, 5.625, 5.675, 5.725, 5.775, 5.825, 5.875, 5.925, 5.975,
        6.025, 6.075, 6.125, 6.175, 6.225, 6.275, 6.325, 6.375, 6.425, 6.475,
        6.525, 6.575, 6.625, 6.675, 6.725, 6.775, 6.825, 6.875, 6.925, 6.975,
        7.025, 7.075, 7.125, 7.175, 7.225, 7.275, 7.325, 7.375, 7.425, 7.475,
        7.525, 7.575, 7.625, 7.675, 7.725, 7.775, 7.825, 7.875, 7.925, 7.975,
        8.025, 8.075, 8.125, 8.175, 8.225, 8.275, 8.325, 8.375, 8.425, 8.475,
        8.525, 8.575, 8.625, 8.675, 8.725, 8.775, 8.825, 8.875, 8.925, 8.975,
        9.025, 9.075, 9.125, 9.175, 9.225, 9.275, 9.325, 9.375, 9.425, 9.475,
        9.525, 9.575, 9.625, 9.675, 9.725, 9.775, 9.825, 9.875, 9.925, 9.975,
        10.025, 10.075, 10.125, 10.175, 10.225, 10.275, 10.325, 10.375, 10.425, 10.475,
        10.525, 10.575, 10.625, 10.675, 10.725, 10.775, 10.825, 10.875, 10.925, 10.975,
        11.025, 11.075, 11.125, 11.175, 11.225, 11.275, 11.325, 11.375, 11.425, 11.475,
        11.525, 11.575, 11.625, 11.675, 11.725, 11.775, 11.825, 11.875, 11.925, 11.975};

    Double_t eff_y[] = {
        79.8624, 77.3101, 78.4156, 77.0169, 78.2946, 78.4749, 76.8986, 77.6381, 78.0087, 79.2948,
        79.8064, 76.8142, 78.6929, 78.3796, 77.7376, 79.2093, 77.1689, 78.3365, 79.9085, 79.1383,
        76.7857, 77.155, 78.0368, 77.9184, 79.1841, 78.9003, 77.8128, 77.1454, 78.2729, 77.2467,
        79.5308, 79.8041, 77.9653, 78.6969, 79.6656, 79.3945, 79.6013, 76.9658, 77.5932, 77.8737,
        78.5182, 78.2298, 75.7215, 78.3213, 78.4886, 78.7562, 77.6419, 78.6172, 76.9302, 77.3876,
        79.486, 78.9788, 78.4387, 79.03, 77.9215, 78.4076, 77.5641, 76.9476, 79.0905, 79.5181,
        78.3747, 77.4301, 77.623, 77.4179, 79.1741, 77.4059, 77.3043, 78.0657, 76.8411, 77.5894,
        78.3161, 78.1852, 79.8518, 78.2254, 79.1058, 77.8548, 79.4602, 77.0698, 78.2648, 79.3658,
        79.1827, 78.6848, 77.3376, 77.6631, 77.1242, 77.0805, 76.1568, 76.5566, 75.9284, 79.2275,
        77.9825, 78.3448, 76.8484, 77.1824, 77.2263, 76.1571, 75.0, 75.1295, 74.6777, 75.2549,
        75.0228, 75.3352, 73.0699, 75.6614, 73.7651, 72.6862, 71.4684, 74.4259, 72.0127, 73.0946,
        72.719, 71.9266, 69.4737, 72.3521, 68.9334, 70.1892, 69.4628, 69.9599, 70.0047, 67.2048,
        66.5461, 66.9867, 65.4529, 69.0884, 65.4702, 65.9266, 65.5337, 65.1579, 62.7248, 63.8543,
        59.3126, 62.2827, 60.9309, 61.7174, 59.9825, 59.7816, 57.9824, 59.7915, 57.8852, 56.5137,
        58.041, 57.6659, 56.7262, 58.3371, 55.8879, 55.3122, 55.8355, 52.5352, 54.6089, 53.3147,
        55.5455, 51.4259, 53.8143, 54.2582, 53.067, 52.8851, 50.947, 53.7431, 49.2049, 49.6324,
        50.0909, 49.4831, 49.0187, 49.6938, 47.869, 49.7963, 49.0824, 46.701, 47.4884, 46.8313,
        45.3729, 47.5091, 48.892, 47.4532, 44.6404, 44.8856, 46.0706, 44.0525, 41.6784, 42.0657,
        41.3617, 44.1441, 43.247, 41.5753, 39.3925, 42.0276, 40.0879, 40.4408, 37.4008, 36.3776,
        35.8338, 36.654, 36.2894, 36.7823, 35.3055, 35.5403, 37.5899, 46.4151, 49.6088, 48.8506,
        58.2418, 57.6023, 56.8421, 58.6667, 60.2094, 68.3871, 71.5328, 73.0769, 77.7778, 79.1045,
        67.2414, 56.25, 70.9677, 95.0, 50.0, 76.9231, 92.3077, 87.5, 100.0, 100.0,
        85.7143, 60.0, 0.0, 100.0, 50.0, 100.0, 0.0, 100.0, 100.0, 100.0,
        0.0};

    // Filter out NaN / zero-efficiency points and keep only the reliable region
    // (the data becomes very noisy with poor statistics above ~12 MeV)
    const int nRaw = 231; // entries listed above (up to 11.975 MeV + a few beyond)
    std::vector<double> vx, vy;

    for (int i = 0; i < nRaw; i++)
    {
        // Skip points with 0% efficiency or clearly unreliable (stats too low)
        if (eff_y[i] > 0. && eff_y[i] <= 100. && eff_x[i] <= 12.0)
        {
            vx.push_back(eff_x[i]);
            vy.push_back(eff_y[i] / 100.); // convert % -> fraction
        }
    }

    auto *gEff = new TGraph(vx.size(), vx.data(), vy.data());
    gEff->SetName("gEff");
    gEff->SetTitle("Neutron detection efficiency;E_{rel} [MeV];#varepsilon");

    // =====================================================================
    //  2. Data files & graphical cuts  (same as inclusive_xs.C)
    // =====================================================================

    std::string ox23WithNeulandFilePath =
        "/nucl_lustre/pablogrusell/g249/g249_analysis/results/dataFiles/23O_analyzed.root";

    auto *reacted_23o_cutg = new TCutG("reacted_23o_cutg", 13);
    reacted_23o_cutg->SetVarX("AoQ_frag");
    reacted_23o_cutg->SetVarY("Z_frag_est");
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

    double opaMin = 1.25;
    double opaMax = 1.65;
    TString condOpa = Form("califa_opa > %f && califa_opa < %f", opaMin, opaMax);
    TString fullCut = TString("reacted_23o_cutg && Erel*1000 < 10. && ") + condOpa;

    // =====================================================================
    //  3. Fill the raw Erel histogram
    // =====================================================================

    auto *fIn = new TFile(ox23WithNeulandFilePath.c_str(), "READ");
    if (!fIn || fIn->IsZombie())
    {
        std::cerr << "ERROR: cannot open " << ox23WithNeulandFilePath << std::endl;
        return;
    }
    auto *tree = static_cast<TTree *>(fIn->Get("KinTree"));

    const int nBins = 70;
    const double erelMin = 0.;
    const double erelMax = 10.1;

    auto *hRaw = new TH1F("hRaw", "E_{rel} sin corregir;E_{rel} [MeV];Counts / bin",
                          nBins, erelMin, erelMax);
    hRaw->Sumw2();

    tree->Draw("Erel*1000>>hRaw", fullCut, "goff");

    // =====================================================================
    //  4. Build the efficiency-corrected histogram  (bin-by-bin)
    // =====================================================================

    auto *hCorr = static_cast<TH1F *>(hRaw->Clone("hCorr"));
    hCorr->SetTitle("E_{rel} corregido por eficiencia;E_{rel} [MeV];Counts / #varepsilon");
    hCorr->Reset();

    for (int ib = 1; ib <= hRaw->GetNbinsX(); ib++)
    {
        double energy = hRaw->GetBinCenter(ib);
        double counts = hRaw->GetBinContent(ib);
        double err = hRaw->GetBinError(ib);

        // Interpolate efficiency from the TGraph
        double eff = gEff->Eval(energy);

        // Safety: skip bins outside the reliable efficiency range or with eff <= 0
        if (eff <= 0.01 || energy > 12.0)
        {
            hCorr->SetBinContent(ib, 0);
            hCorr->SetBinError(ib, 0);
            continue;
        }

        hCorr->SetBinContent(ib, counts / eff);
        hCorr->SetBinError(ib, err / eff);
    }

    // =====================================================================
    //  5. Draw
    // =====================================================================

    gStyle->SetOptStat(0);

    auto *c = new TCanvas("cErel", "Erel corrected", 1200, 800);
    c->Divide(2, 2);

    // --- Pad 1: Efficiency graph ---
    c->cd(1);
    gEff->SetLineColor(kTeal + 7);
    gEff->SetLineWidth(2);
    gEff->SetMarkerStyle(20);
    gEff->SetMarkerSize(0.4);
    gEff->SetMarkerColor(kTeal + 7);
    gEff->Draw("APL");

    // --- Pad 2: Raw Erel ---
    c->cd(2);
    hRaw->SetLineColor(kAzure + 1);
    hRaw->SetLineWidth(2);
    hRaw->GetXaxis()->SetRangeUser(0, 12);
    hRaw->Draw("HIST");

    // --- Pad 3: Corrected Erel ---
    c->cd(3);
    hCorr->SetLineColor(kRed + 1);
    hCorr->SetLineWidth(2);
    hCorr->GetXaxis()->SetRangeUser(0, 12);
    hCorr->Draw("HIST");

    // --- Pad 4: Overlay (normalised shapes) ---
    c->cd(4);
    auto *hRawNorm = static_cast<TH1F *>(hRaw->Clone("hRawNorm"));
    auto *hCorrNorm = static_cast<TH1F *>(hCorr->Clone("hCorrNorm"));

    if (hRawNorm->Integral() > 0)
        hRawNorm->Scale(1. / hRawNorm->Integral());
    if (hCorrNorm->Integral() > 0)
        hCorrNorm->Scale(1. / hCorrNorm->Integral());

    hRawNorm->SetTitle("Comparacion normalizada;E_{rel} [MeV];u.a.");
    hRawNorm->GetXaxis()->SetRangeUser(0, 12);
    hRawNorm->SetLineColor(kAzure + 1);
    hCorrNorm->SetLineColor(kRed + 1);

    hRawNorm->Draw("HIST E");
    hCorrNorm->Draw("HIST E SAME");

    auto *leg = new TLegend(0.55, 0.75, 0.88, 0.88);
    leg->AddEntry(hRawNorm, "Sin corregir", "l");
    leg->AddEntry(hCorrNorm, "Corregido", "l");
    leg->Draw();

    c->Update();
    c->SaveAs("erel_corrected.pdf");

    std::cout << "\n=== Done. Output saved to erel_corrected.pdf ===\n"
              << std::endl;
}