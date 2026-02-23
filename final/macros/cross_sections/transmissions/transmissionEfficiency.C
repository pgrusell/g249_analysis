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

// 3-sigma elliptical cut for outgoing 24O identification (from eventFilter)
static bool insideEllipse24O(double aoq, double z)
{
    constexpr double mu_AoQ = 2.95766;
    constexpr double sig_AoQ = 0.0221273;
    constexpr double mu_Z = 8.06447;
    constexpr double sig_Z = 0.201233;
    constexpr double k = 3.0;
    double d2 = std::pow((aoq - mu_AoQ) / sig_AoQ, 2) +
                std::pow((z - mu_Z) / sig_Z, 2);
    return d2 < k * k;
}

// 3-sigma elliptical cut for unreacted 25F identification (from eventFilter)
static bool insideEllipse25F(double aoq, double z)
{
    constexpr double mu_AoQ = 2.740;
    constexpr double sig_AoQ = 0.00988;
    constexpr double mu_Z = 8.976;
    constexpr double sig_Z = 0.1973;
    constexpr double k = 3.0;
    double d2 = std::pow((aoq - mu_AoQ) / sig_AoQ, 2) +
                std::pow((z - mu_Z) / sig_Z, 2);
    return d2 < k * k;
}

void transmissionEfficiency()
{
    auto *ch = new TChain("evt");
    TString path = "/nucl_lustre/pablogrusell/g249/root_files/unpackedData/";
    ch->AddFile(path + "g249_all_det_offline_0001_20260224_103705_160_8.root");

    // Basic data
    auto *frsDataBranch = new TClonesArray("R3BFrsData");
    auto *footHitDataBranch = new TClonesArray("R3BFootHitData");
    auto *tofdHitDataBranch = new TClonesArray("R3BTofdHitData");
    auto *fib30DataBranch = new TClonesArray("R3BBunchedFiberHitData");
    auto *fib31DataBranch = new TClonesArray("R3BBunchedFiberHitData");
    auto *fib33DataBranch = new TClonesArray("R3BBunchedFiberHitData");
    auto *fib32DataBranch = new TClonesArray("R3BBunchedFiberHitData");

    // Higher levels
    auto *incomingFootDataBranch = new TClonesArray("R3BTrackingParticle");
    auto *outgoingFootDataBranch = new TClonesArray("R3BTrackingParticle");
    auto *mdfDataBranch = new TClonesArray("R3BTrackingParticle");

    ch->SetBranchAddress("FrsData", &frsDataBranch);
    ch->SetBranchAddress("FootHitData", &footHitDataBranch);
    ch->SetBranchAddress("TofdHit", &tofdHitDataBranch);
    ch->SetBranchAddress("Fi30Hit", &fib30DataBranch);
    ch->SetBranchAddress("Fi31Hit", &fib31DataBranch);
    ch->SetBranchAddress("Fi32Hit", &fib32DataBranch);
    ch->SetBranchAddress("Fi33Hit", &fib33DataBranch);

    ch->SetBranchAddress("IncomingTrackFoot", &incomingFootDataBranch);
    ch->SetBranchAddress("OutgoingTrackFoot", &outgoingFootDataBranch);
    ch->SetBranchAddress("FragmentMDFTrack", &mdfDataBranch);

    int nEntries = ch->GetEntries();

    // ---- Cut counters ----
    int nTotal = nEntries;
    int nFrsOk = 0;          // 1) Exactly 1 FRS entry
    int nIncoming25F = 0;    // 2) Inside 25F polygon
    int nUpstreamFoot = 0;   // 3) FOOT 1-4 valid
    int nChargeFoot5_8 = 0;  // 4) Charge cut on FOOT 5-8
    int nDownstreamFoot = 0; // 5) FOOT 5-8 valid
    int nFibers = 0;         // 6) Fiber detectors
    int nTofdPlane1 = 0;     // 7) ToFD plane 1 present
    int nTofdCharge = 0;     // 8) ToFD charge cut

    // ---- 2D histograms: position vs Eloss (filled after charge cut on FOOT 5-8) ----
    // All events passing charge cut
    auto *hFib30 = new TH2D("hFib30", "Fib30: Y vs Eloss;fib30Y [cm];Eloss", 160, -80, 80, 100, 0, 30);
    auto *hFib31 = new TH2D("hFib31", "Fib31: X vs Eloss;fib31X [cm];Eloss", 160, -80, 80, 100, 0, 30);
    auto *hFib32 = new TH2D("hFib32", "Fib32: X vs Eloss;fib32X [cm];Eloss", 160, -80, 80, 100, 0, 30);
    auto *hFib33 = new TH2D("hFib33", "Fib33: X vs Eloss;fib33X [cm];Eloss", 160, -80, 80, 100, 0, 30);
    auto *hTofdXvsZ = new TH2D("hTofdXvsZ", "ToFD: X vs Z;Z_{Plane 1};X [cm]", 300, 0, 15, 300, -60, 60);

    // 24O MDF-selected events (3#sigma ellipse on AoQ_frag vs Z_frag)
    auto *hFib30_24O = new TH2D("hFib30_24O", "^{24}O Fib30: Y vs Eloss;fib30Y [cm];Eloss", 160, -80, 80, 100, 0, 30);
    auto *hFib31_24O = new TH2D("hFib31_24O", "^{24}O Fib31: X vs Eloss;fib31X [cm];Eloss", 160, -80, 80, 100, 0, 30);
    auto *hFib32_24O = new TH2D("hFib32_24O", "^{24}O Fib32: X vs Eloss;fib32X [cm];Eloss", 160, -80, 80, 100, 0, 30);
    auto *hFib33_24O = new TH2D("hFib33_24O", "^{24}O Fib33: X vs Eloss;fib33X [cm];Eloss", 160, -80, 80, 100, 0, 30);
    auto *hTofdXvsZ_24O = new TH2D("hTofdXvsZ_24O", "^{24}O ToFD: X vs Z;Z_{Plane 1};X [cm]", 300, 0, 15, 300, -60, 60);
    auto *hAoQvsZ_MDF = new TH2D("hAoQvsZ_MDF", "MDF: AoQ vs Z (after charge cut);AoQ;Z", 200, 2.5, 3.2, 200, 5, 12);

    // ---- ToFD charge vs X-position approach (from slides) ----
    // Z=8 in FOOTs: all events, then Z=8 + MDF exists
    auto *hTofdZvsX_Z8 = new TH2D("hTofdZvsX_Z8",
                                  "ToFD: Z vs X (FOOT Z=8);X [cm];Z_{ToFD}", 300, -60, 60, 300, 0, 15);
    auto *hTofdZvsX_Z8_MDF = new TH2D("hTofdZvsX_Z8_MDF",
                                      "ToFD: Z vs X (FOOT Z=8 + MDF exists);X [cm];Z_{ToFD}", 300, -60, 60, 300, 0, 15);

    // Z=9 in FOOTs: all events, then Z=9 + MDF exists
    auto *hTofdZvsX_Z9 = new TH2D("hTofdZvsX_Z9",
                                  "ToFD: Z vs X (FOOT Z=9);X [cm];Z_{ToFD}", 300, -60, 60, 300, 0, 15);
    auto *hTofdZvsX_Z9_MDF = new TH2D("hTofdZvsX_Z9_MDF",
                                      "ToFD: Z vs X (FOOT Z=9 + MDF exists);X [cm];Z_{ToFD}", 300, -60, 60, 300, 0, 15);

    int n24O = 0;
    int n25F_mdf = 0;

    for (auto i = 0; i < nEntries; i++)
    {
        ch->GetEntry(i);

        ///////// Incoming part //////////////

        // 1) FRS Data
        if (frsDataBranch->GetEntriesFast() != 1)
            continue;
        nFrsOk++;

        auto *frsData = static_cast<R3BFrsData *>(frsDataBranch->At(0));
        double aoq_In = frsData->GetAq();
        double z_In = frsData->GetZ();

        // 2) 25F as incoming
        if (!isInsidePolygon(aoq_In, z_In, INCOMING_25F_POLYGON))
            continue;
        nIncoming25F++;

        // FOOT Data
        std::vector<bool> isValidFoot(8);
        std::vector<double> charge(8);
        for (int j = 0; j < footHitDataBranch->GetEntriesFast(); j++)
        {
            auto *hit = static_cast<R3BFootHitData *>(footHitDataBranch->At(j));
            isValidFoot[hit->GetDetId() - 1] = true;
            charge[hit->GetDetId() - 1] = hit->GetZCharge();
        }

        // 3) Upstream FOOT (1-4)
        if (!isValidFoot[0] || !isValidFoot[1] || !isValidFoot[2] || !isValidFoot[3])
            continue;
        nUpstreamFoot++;

        ///////// MDF requirement /////////////

        double chMinO = 7.5; // Z=8 (Oxygen)
        double chMaxO = 8.5;
        double chMinF = 8.5; // Z=9 (Fluorine)
        double chMaxF = 9.5;

        // ---- Check if all downstream FOOTs are consistent with Z=8 or Z=9 ----
        bool isZ8 = (charge[4] >= chMinO && charge[4] <= chMaxO &&
                     charge[5] >= chMinO && charge[5] <= chMaxO &&
                     charge[6] >= chMinO && charge[6] <= chMaxO &&
                     charge[7] >= chMinO && charge[7] <= chMaxO);

        bool isZ9 = (charge[4] >= chMinF && charge[4] <= chMaxF &&
                     charge[5] >= chMinF && charge[5] <= chMaxF &&
                     charge[6] >= chMinF && charge[6] <= chMaxF &&
                     charge[7] >= chMinF && charge[7] <= chMaxF);

        // ---- Z=8 path: fill ToFD charge vs X, then check MDF for 24O ----
        if (isZ8)
        {
            nChargeFoot5_8++;

            // Fill generic 2D histograms (fibers + ToFD)
            for (int j = 0; j < fib30DataBranch->GetEntriesFast(); j++)
            {
                auto *hit = static_cast<R3BBunchedFiberHitData *>(fib30DataBranch->At(j));
                hFib30->Fill(hit->GetY(), hit->GetEloss());
            }
            for (int j = 0; j < fib31DataBranch->GetEntriesFast(); j++)
            {
                auto *hit = static_cast<R3BBunchedFiberHitData *>(fib31DataBranch->At(j));
                hFib31->Fill(hit->GetX(), hit->GetEloss());
            }
            for (int j = 0; j < fib32DataBranch->GetEntriesFast(); j++)
            {
                auto *hit = static_cast<R3BBunchedFiberHitData *>(fib32DataBranch->At(j));
                hFib32->Fill(hit->GetX(), hit->GetEloss());
            }
            for (int j = 0; j < fib33DataBranch->GetEntriesFast(); j++)
            {
                auto *hit = static_cast<R3BBunchedFiberHitData *>(fib33DataBranch->At(j));
                hFib33->Fill(hit->GetX(), hit->GetEloss());
            }
            for (int j = 0; j < tofdHitDataBranch->GetEntriesFast(); j++)
            {
                auto *hit = static_cast<R3BTofdHitData *>(tofdHitDataBranch->At(j));
                if (hit->GetDetId() == 1)
                    hTofdXvsZ->Fill(hit->GetEloss(), hit->GetX());
            }

            // ToFD charge vs X for all Z=8 FOOT events
            for (int j = 0; j < tofdHitDataBranch->GetEntriesFast(); j++)
            {
                auto *hit = static_cast<R3BTofdHitData *>(tofdHitDataBranch->At(j));
                if (hit->GetDetId() == 1)
                    hTofdZvsX_Z8->Fill(hit->GetX(), hit->GetEloss());
            }

            // MDF identification
            if (mdfDataBranch->GetEntriesFast() > 0)
            {
                auto *mdf = static_cast<R3BTrackingParticle *>(mdfDataBranch->At(0));
                double aoqFrag = mdf->GetMass();
                double zFrag = mdf->GetCharge();

                hAoQvsZ_MDF->Fill(aoqFrag, zFrag);

                // Z=8 + MDF exists (no PID selection)
                for (int j = 0; j < tofdHitDataBranch->GetEntriesFast(); j++)
                {
                    auto *hit = static_cast<R3BTofdHitData *>(tofdHitDataBranch->At(j));
                    if (hit->GetDetId() == 1)
                        hTofdZvsX_Z8_MDF->Fill(hit->GetX(), hit->GetEloss());
                }

                // 24O ellipse selection (for fiber/ToFD 2D plots)
                if (insideEllipse24O(aoqFrag, zFrag))
                {
                    n24O++;

                    for (int j = 0; j < fib30DataBranch->GetEntriesFast(); j++)
                    {
                        auto *hit = static_cast<R3BBunchedFiberHitData *>(fib30DataBranch->At(j));
                        hFib30_24O->Fill(hit->GetY(), hit->GetEloss());
                    }
                    for (int j = 0; j < fib31DataBranch->GetEntriesFast(); j++)
                    {
                        auto *hit = static_cast<R3BBunchedFiberHitData *>(fib31DataBranch->At(j));
                        hFib31_24O->Fill(hit->GetX(), hit->GetEloss());
                    }
                    for (int j = 0; j < fib32DataBranch->GetEntriesFast(); j++)
                    {
                        auto *hit = static_cast<R3BBunchedFiberHitData *>(fib32DataBranch->At(j));
                        hFib32_24O->Fill(hit->GetX(), hit->GetEloss());
                    }
                    for (int j = 0; j < fib33DataBranch->GetEntriesFast(); j++)
                    {
                        auto *hit = static_cast<R3BBunchedFiberHitData *>(fib33DataBranch->At(j));
                        hFib33_24O->Fill(hit->GetX(), hit->GetEloss());
                    }
                    for (int j = 0; j < tofdHitDataBranch->GetEntriesFast(); j++)
                    {
                        auto *hit = static_cast<R3BTofdHitData *>(tofdHitDataBranch->At(j));
                        if (hit->GetDetId() == 1)
                            hTofdXvsZ_24O->Fill(hit->GetEloss(), hit->GetX());
                    }
                }
            }
        }

        // ---- Z=9 path: fill ToFD charge vs X, then check if MDF exists ----
        if (isZ9)
        {
            // ToFD charge vs X for all Z=9 FOOT events
            for (int j = 0; j < tofdHitDataBranch->GetEntriesFast(); j++)
            {
                auto *hit = static_cast<R3BTofdHitData *>(tofdHitDataBranch->At(j));
                if (hit->GetDetId() == 1)
                    hTofdZvsX_Z9->Fill(hit->GetX(), hit->GetEloss());
            }

            // Z=9 + MDF exists (no PID selection)
            if (mdfDataBranch->GetEntriesFast() > 0)
            {
                n25F_mdf++;

                for (int j = 0; j < tofdHitDataBranch->GetEntriesFast(); j++)
                {
                    auto *hit = static_cast<R3BTofdHitData *>(tofdHitDataBranch->At(j));
                    if (hit->GetDetId() == 1)
                        hTofdZvsX_Z9_MDF->Fill(hit->GetX(), hit->GetEloss());
                }
            }
        }

        // Continue cut flow only for Z=8
        if (!isZ8)
            continue;

        // 5) All downstream FOOTs valid
        if (!isValidFoot[4] || !isValidFoot[5] || !isValidFoot[6] || !isValidFoot[7])
            continue;
        nDownstreamFoot++;

        // 6) All fibers
        if (fib30DataBranch->GetEntriesFast() == 0)
            continue;
        if (fib32DataBranch->GetEntriesFast() == 0)
            continue;
        if (fib33DataBranch->GetEntriesFast() == 0 && fib31DataBranch->GetEntriesFast() == 0)
            continue;
        nFibers++;

        // 7) ToFD first plane present
        bool isGoodTofd = false;
        double chargeTofD = 0;
        for (int j = 0; j < tofdHitDataBranch->GetEntriesFast(); j++)
        {
            auto *hit = static_cast<R3BTofdHitData *>(tofdHitDataBranch->At(j));
            if (hit->GetDetId() == 1)
            {
                isGoodTofd = true;
                chargeTofD = hit->GetEloss();
            }
        }

        if (!isGoodTofd)
            continue;
        nTofdPlane1++;

        // 8) ToFD charge cut
        if (chargeTofD < 7.5 || chargeTofD > 8.5)
            continue;
        nTofdCharge++;
    }

    // ---- Print summary ----
    auto pct = [&](int n)
    { return (nTotal > 0) ? 100.0 * n / nTotal : 0.0; };
    std::cout << "=== Cut flow ===" << std::endl;
    std::cout << Form("Total entries:          %d (%.2f%%)", nTotal, pct(nTotal)) << std::endl;
    std::cout << Form("1) FRS == 1:            %d (%.2f%%)", nFrsOk, pct(nFrsOk)) << std::endl;
    std::cout << Form("2) 25F polygon:         %d (%.2f%%)", nIncoming25F, pct(nIncoming25F)) << std::endl;
    std::cout << Form("3) Upstream FOOT 1-4:   %d (%.2f%%)", nUpstreamFoot, pct(nUpstreamFoot)) << std::endl;
    std::cout << Form("4) Charge FOOT 5-8:     %d (%.2f%%)", nChargeFoot5_8, pct(nChargeFoot5_8)) << std::endl;
    std::cout << Form("5) Downstream FOOT 5-8: %d (%.2f%%)", nDownstreamFoot, pct(nDownstreamFoot)) << std::endl;
    std::cout << Form("6) Fibers:              %d (%.2f%%)", nFibers, pct(nFibers)) << std::endl;
    std::cout << Form("7) ToFD plane 1:        %d (%.2f%%)", nTofdPlane1, pct(nTofdPlane1)) << std::endl;
    std::cout << Form("8) ToFD charge:         %d (%.2f%%)", nTofdCharge, pct(nTofdCharge)) << std::endl;

    // ---- Bar chart histogram ----
    const int nCuts = 9;
    const char *cutLabels[nCuts] = {
        "Total",
        "FRS==1",
        "^{25}F polygon",
        "FOOT 1-4",
        "Z(FOOT 5-8)",
        "FOOT 5-8",
        "Fibers",
        "ToFD plane 1",
        "ToFD Z cut"};
    int cutValues[nCuts] = {
        nTotal, nFrsOk, nIncoming25F, nUpstreamFoot,
        nChargeFoot5_8, nDownstreamFoot, nFibers, nTofdPlane1, nTofdCharge};

    auto *c1 = new TCanvas("c1", "Cut Flow - Transmission Efficiency", 1200, 600);
    c1->SetLogy();
    c1->SetGridy();
    c1->SetBottomMargin(0.18);

    auto *hCutFlow = new TH1D("hCutFlow", "Cut flow: ^{25}F transmission efficiency;Selection cut;Entries", nCuts, 0, nCuts);
    hCutFlow->SetStats(0);

    for (int k = 0; k < nCuts; k++)
    {
        hCutFlow->SetBinContent(k + 1, cutValues[k]);
        hCutFlow->GetXaxis()->SetBinLabel(k + 1, cutLabels[k]);
    }

    hCutFlow->SetFillColor(kAzure + 1);
    hCutFlow->SetFillStyle(1001);
    hCutFlow->SetLineColor(kAzure + 3);
    hCutFlow->SetLineWidth(2);
    hCutFlow->GetXaxis()->SetLabelSize(0.045);
    hCutFlow->GetXaxis()->SetLabelOffset(0.01);
    hCutFlow->GetYaxis()->SetLabelSize(0.04);
    hCutFlow->GetYaxis()->SetTitleSize(0.045);
    hCutFlow->SetMinimum(0.5);

    hCutFlow->Draw("BAR");

    // Add text labels on top of each bar (percentage w.r.t. total)
    auto *tex = new TLatex();
    tex->SetTextSize(0.03);
    tex->SetTextAlign(21); // center-bottom
    tex->SetTextFont(42);
    for (int k = 0; k < nCuts; k++)
    {
        double xpos = hCutFlow->GetBinCenter(k + 1);
        double ypos = cutValues[k] * 1.3; // slightly above bar in log scale
        double pct = (nTotal > 0) ? 100.0 * cutValues[k] / nTotal : 0.0;
        tex->DrawLatex(xpos, ypos, Form("%.2f%%", pct));
    }

    c1->Update();
    c1->SaveAs("cutflow_transmissionEfficiency.png");
    c1->SaveAs("cutflow_transmissionEfficiency.pdf");

    std::cout << "\nPlot saved as cutflow_transmissionEfficiency.png/.pdf" << std::endl;

    // ==================== 2D plots: Fiber Eloss vs Position (after charge cut) ====================
    gStyle->SetOptStat(0);

    auto *cFibEloss = new TCanvas("cFibEloss", "Fiber Eloss vs Position (after Z cut on FOOT 5-8)", 1600, 400);
    cFibEloss->Divide(4, 1);

    cFibEloss->cd(1);
    hFib30->Draw("COLZ");

    cFibEloss->cd(2);
    hFib31->Draw("COLZ");

    cFibEloss->cd(3);
    hFib33->Draw("COLZ");

    cFibEloss->cd(4);
    hFib32->Draw("COLZ");

    cFibEloss->Update();
    cFibEloss->SaveAs("fibEloss_afterChargeCut.png");
    cFibEloss->SaveAs("fibEloss_afterChargeCut.pdf");

    // ==================== 2D plot: ToFD X vs Z (after charge cut) ====================
    auto *cToFD = new TCanvas("cToFD", "ToFD X vs Z (after Z cut on FOOT 5-8)", 700, 500);
    hTofdXvsZ->Draw("COLZ");
    cToFD->Update();
    cToFD->SaveAs("tofdXvsZ_afterChargeCut.png");
    cToFD->SaveAs("tofdXvsZ_afterChargeCut.pdf");

    std::cout << "2D plots saved as fibEloss_afterChargeCut and tofdXvsZ_afterChargeCut (.png/.pdf)" << std::endl;
    std::cout << Form("\n24O events (3sigma ellipse on MDF): %d (%.2f%% of total)", n24O,
                      (nTotal > 0) ? 100.0 * n24O / nTotal : 0.0)
              << std::endl;

    // ==================== MDF AoQ vs Z (after charge cut) ====================
    auto *cMDF = new TCanvas("cMDF", "MDF: AoQ vs Z (after Z cut on FOOT 5-8)", 800, 600);
    hAoQvsZ_MDF->Draw("COLZ");
    // Draw 24O ellipse
    auto *ell24O = new TEllipse(2.95766, 8.06447, 3.0 * 0.0221273, 3.0 * 0.201233);
    ell24O->SetFillStyle(0);
    ell24O->SetLineColor(kRed);
    ell24O->SetLineWidth(2);
    ell24O->Draw("SAME");
    // Draw 25F ellipse
    auto *ell25F = new TEllipse(2.740, 8.976, 3.0 * 0.00988, 3.0 * 0.1973);
    ell25F->SetFillStyle(0);
    ell25F->SetLineColor(kBlue);
    ell25F->SetLineWidth(2);
    ell25F->Draw("SAME");
    auto *legMDF = new TLegend(0.65, 0.75, 0.88, 0.88);
    legMDF->SetBorderSize(0);
    legMDF->AddEntry(ell24O, "^{24}O (3#sigma)", "l");
    legMDF->AddEntry(ell25F, "^{25}F (3#sigma)", "l");
    legMDF->Draw();
    cMDF->Update();
    cMDF->SaveAs("mdf_AoQvsZ_afterChargeCut.png");
    cMDF->SaveAs("mdf_AoQvsZ_afterChargeCut.pdf");

    // ==================== 24O: Fiber Eloss vs Position ====================
    auto *cFibEloss24O = new TCanvas("cFibEloss24O", "24O Fiber Eloss vs Position (3#sigma MDF cut)", 1600, 400);
    cFibEloss24O->Divide(4, 1);

    cFibEloss24O->cd(1);
    hFib30_24O->Draw("COLZ");

    cFibEloss24O->cd(2);
    hFib31_24O->Draw("COLZ");

    cFibEloss24O->cd(3);
    hFib33_24O->Draw("COLZ");

    cFibEloss24O->cd(4);
    hFib32_24O->Draw("COLZ");

    cFibEloss24O->Update();
    cFibEloss24O->SaveAs("fibEloss_24O_mdfCut.png");
    cFibEloss24O->SaveAs("fibEloss_24O_mdfCut.pdf");

    // ==================== 24O: ToFD X vs Z ====================
    auto *cToFD24O = new TCanvas("cToFD24O", "24O ToFD X vs Z (3#sigma MDF cut)", 700, 500);
    hTofdXvsZ_24O->Draw("COLZ");
    cToFD24O->Update();
    cToFD24O->SaveAs("tofdXvsZ_24O_mdfCut.png");
    cToFD24O->SaveAs("tofdXvsZ_24O_mdfCut.pdf");

    std::cout << "24O plots saved as fibEloss_24O_mdfCut and tofdXvsZ_24O_mdfCut (.png/.pdf)" << std::endl;

    // ==================== Comparison: ToFD charge vs X-position ====================
    // Z=8: FOOT selection vs FOOT + MDF exists
    auto *cCompZ8 = new TCanvas("cCompZ8", "Z=8: ToFD charge vs X — FOOT vs FOOT+MDF", 1400, 600);
    cCompZ8->Divide(2, 1);

    cCompZ8->cd(1);
    hTofdZvsX_Z8->Draw("COLZ");
    auto *texZ8a = new TLatex();
    texZ8a->SetNDC();
    texZ8a->SetTextSize(0.04);
    texZ8a->DrawLatex(0.15, 0.92, Form("FOOT Z=8: %d entries", (int)hTofdZvsX_Z8->GetEntries()));

    cCompZ8->cd(2);
    hTofdZvsX_Z8_MDF->Draw("COLZ");
    auto *texZ8b = new TLatex();
    texZ8b->SetNDC();
    texZ8b->SetTextSize(0.04);
    texZ8b->DrawLatex(0.15, 0.92, Form("FOOT Z=8 + MDF: %d entries", (int)hTofdZvsX_Z8_MDF->GetEntries()));

    cCompZ8->Update();
    cCompZ8->SaveAs("tofd_chargeVsX_Z8_comparison.png");
    cCompZ8->SaveAs("tofd_chargeVsX_Z8_comparison.pdf");

    // Z=9: FOOT selection vs FOOT + MDF exists
    auto *cCompZ9 = new TCanvas("cCompZ9", "Z=9: ToFD charge vs X — FOOT vs FOOT+MDF", 1400, 600);
    cCompZ9->Divide(2, 1);

    cCompZ9->cd(1);
    hTofdZvsX_Z9->Draw("COLZ");
    auto *texZ9a = new TLatex();
    texZ9a->SetNDC();
    texZ9a->SetTextSize(0.04);
    texZ9a->DrawLatex(0.15, 0.92, Form("FOOT Z=9: %d entries", (int)hTofdZvsX_Z9->GetEntries()));

    cCompZ9->cd(2);
    hTofdZvsX_Z9_MDF->Draw("COLZ");
    auto *texZ9b = new TLatex();
    texZ9b->SetNDC();
    texZ9b->SetTextSize(0.04);
    texZ9b->DrawLatex(0.15, 0.92, Form("FOOT Z=9 + MDF: %d entries", (int)hTofdZvsX_Z9_MDF->GetEntries()));

    cCompZ9->Update();
    cCompZ9->SaveAs("tofd_chargeVsX_Z9_comparison.png");
    cCompZ9->SaveAs("tofd_chargeVsX_Z9_comparison.pdf");

    // ---- Entry comparison summary ----
    double entriesZ8 = hTofdZvsX_Z8->GetEntries();
    double entriesZ8_MDF = hTofdZvsX_Z8_MDF->GetEntries();
    double entriesZ9 = hTofdZvsX_Z9->GetEntries();
    double entriesZ9_MDF = hTofdZvsX_Z9_MDF->GetEntries();

    std::cout << "\n==================== Entry Comparison ====================" << std::endl;
    std::cout << Form("Z=8 FOOT only:       %.0f entries", entriesZ8) << std::endl;
    std::cout << Form("Z=8 FOOT + MDF:      %.0f entries (%.2f%% of FOOT Z=8)",
                      entriesZ8_MDF, (entriesZ8 > 0) ? 100.0 * entriesZ8_MDF / entriesZ8 : 0.0)
              << std::endl;
    std::cout << Form("   -> Lost w/o MDF:  %.0f entries (%.2f%%)",
                      entriesZ8 - entriesZ8_MDF,
                      (entriesZ8 > 0) ? 100.0 * (entriesZ8 - entriesZ8_MDF) / entriesZ8 : 0.0)
              << std::endl;
    std::cout << Form("Z=9 FOOT only:       %.0f entries", entriesZ9) << std::endl;
    std::cout << Form("Z=9 FOOT + MDF:      %.0f entries (%.2f%% of FOOT Z=9)",
                      entriesZ9_MDF, (entriesZ9 > 0) ? 100.0 * entriesZ9_MDF / entriesZ9 : 0.0)
              << std::endl;
    std::cout << Form("   -> Lost w/o MDF:  %.0f entries (%.2f%%)",
                      entriesZ9 - entriesZ9_MDF,
                      (entriesZ9 > 0) ? 100.0 * (entriesZ9 - entriesZ9_MDF) / entriesZ9 : 0.0)
              << std::endl;
    std::cout << "==========================================================" << std::endl;
}