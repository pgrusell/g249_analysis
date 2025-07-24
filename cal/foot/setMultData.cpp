// This macro analyzes the data from a losHit, footHit and footCal data file, to fill the following histograms:
// :etaFootAsicChar: eta vs energy. One histogram per Foot, asic and charge (8 and 9)
// :etaFootAsic: eta vs energy. One histogram per Foot, asic
// :etaFootAsicMult: eta vs energy. One histogram per Foot, asic and multiplicity
// :posFootChar: energy vs position. One histogram per FOOT and charge in LOS
// :posFootMult: energy vs position. One histogram per FOOT and multiplicity
// :chargeEnergyFootAsic: charge correlation FOOT vs LOS. One histogram per FOOT and ASIC
// :chargeEnergyFootAsicMult: charge correlation FOOT vs LOS. One histogram per FOOT, ASIC and multiplicity
// :hChargeFull: Energy spectra in LOS

// Includes
#include "R3BFootCalData.h"
#include "R3BFootHitData.h"
#include "R3BLosHitData.h"
#include "TClonesArray.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TString.h"
#include "TTree.h"

#include <cmath>
#include <cstdio>
#include <iostream>
#include <vector>

int main(int argc, char **argv)
{
    TString inFilePath =
        "/nucl_lustre/pablogrusell/g249/root_files/etaCorr/g249_data_incoming_online_20250724_095255.root";
    TString outFilePath = "/nucl_lustre/pablogrusell/g249/root_files/etaCorr/foot_etaplot_mult_v5.root";

    if (argc > 1)
        inFilePath = argv[1];
    if (argc > 2)
        outFilePath = argv[2];

    const int nFoots = 8;
    const int nAsics = 10;
    const int nCharges = 2;
    const int nStrips = 640;
    const int maxMult = 4;

    // Open input ROOT file
    auto *inFile = new TFile(inFilePath, "READ");
    if (!inFile || inFile->IsZombie())
    {
        std::cerr << "Error: Cannot open input file " << inFilePath << std::endl;
        return 1;
    }

    auto *evt = (TTree *)inFile->Get("evt");
    if (!evt)
    {
        std::cerr << "Error: Tree 'evt' not found in file." << std::endl;
        return 1;
    }

    // Branches
    evt->SetBranchStatus("*", 0);
    evt->SetBranchStatus("FootHitData*", 1);
    evt->SetBranchStatus("FootCalData*", 1);
    evt->SetBranchStatus("LosHit*", 1);

    auto *footCalDataArray = new TClonesArray("R3BFootCalData");
    evt->SetBranchAddress("FootCalData", &footCalDataArray);

    auto *footDataArray = new TClonesArray("R3BFootHitData");
    evt->SetBranchAddress("FootHitData", &footDataArray);

    auto *losDataArray = new TClonesArray("R3BLosHitData");
    evt->SetBranchAddress("LosHit", &losDataArray);

    // Histograms
    std::vector<std::vector<std::vector<TH2F *>>> etaFootAsicChar(
        nFoots, std::vector<std::vector<TH2F *>>(nAsics, std::vector<TH2F *>(nCharges, nullptr)));
    std::vector<std::vector<TH2F *>> posFootChar(nFoots, std::vector<TH2F *>(nCharges, nullptr));
    std::vector<std::vector<TH2F *>> etaFootAsic(nFoots, std::vector<TH2F *>(nAsics, nullptr));
    std::vector<std::vector<std::vector<TH2F *>>> etaFootAsicMult(
        nFoots, std::vector<std::vector<TH2F *>>(nAsics, std::vector<TH2F *>(maxMult, nullptr)));
    std::vector<TH2F *> posFoot(nFoots, nullptr);
    std::vector<std::vector<TH2F *>> chargeEnergyFootAsic(nFoots, std::vector<TH2F *>(nAsics, nullptr));
    std::vector<TH2F *> stripEnergyFoot(nFoots, nullptr);
    std::vector<std::vector<std::vector<TH2F *>>> posFootMult(
        nFoots, std::vector<std::vector<TH2F *>>(maxMult, std::vector<TH2F *>(1, nullptr)));
    std::vector<std::vector<std::vector<TH2F *>>> chargeEnergyFootAsicMult(
        nFoots, std::vector<std::vector<TH2F *>>(nAsics, std::vector<TH2F *>(maxMult, nullptr)));
    TH1F *hChargeFull = new TH1F("charge_full", "Charge full spectrum", 150, 0, 15);

    for (int iFoot = 0; iFoot < nFoots; ++iFoot)
    {
        auto namePos = Form("pos_foot_%i", iFoot);
        auto titlePos = Form("Energy vs Position for Foot %i", iFoot);
        posFoot[iFoot] = new TH2F(namePos, titlePos, 640, 0, 640, 1000, 0, 5000);

        auto nameStrip = Form("strip_vs_energy_foot_%i", iFoot);
        auto titleStrip = Form("Strip vs Energy for Foot %i", iFoot);
        stripEnergyFoot[iFoot] = new TH2F(nameStrip, titleStrip, nStrips, 0, nStrips, 1000, 0, 5000);

        for (int m = 0; m < maxMult; ++m)
        {
            auto namePosM = Form("pos_foot_%i_mult_%i", iFoot, m + 1);
            auto titlePosM = Form("Energy vs Position for Foot %i, mult %i", iFoot, m + 1);
            posFootMult[iFoot][m][0] = new TH2F(namePosM, titlePosM, 640, 0, 640, 1000, 0, 5000);
        }

        for (int iAsic = 0; iAsic < nAsics; ++iAsic)
        {
            auto nameEta = Form("eta_foot_%i_asic_%i", iFoot, iAsic);
            auto titleEta = Form("Eta vs Energy for Foot %i, ASIC %i", iFoot, iAsic);
            etaFootAsic[iFoot][iAsic] = new TH2F(nameEta, titleEta, 100, 0, 1, 1000, 0, 5000);

            for (int m = 0; m < maxMult; ++m)
            {
                auto nameEtaM = Form("eta_foot_%i_asic_%i_mult_%i", iFoot, iAsic, m + 1);
                auto titleEtaM = Form("Eta vs Energy for Foot %i, ASIC %i, mult %i", iFoot, iAsic, m + 1);
                etaFootAsicMult[iFoot][iAsic][m] = new TH2F(nameEtaM, titleEtaM, 100, 0, 1, 1000, 0, 5000);
            }

            auto nameCharge = Form("charge_full_vs_energy_foot_%i_asic_%i", iFoot, iAsic);
            auto titleCharge = Form("Charge full vs Energy for Foot %i, ASIC %i", iFoot, iAsic);
            chargeEnergyFootAsic[iFoot][iAsic] = new TH2F(nameCharge, titleCharge, 200, 0, 15, 1000, 0, 5000);

            for (int m = 0; m < maxMult; ++m)
            {
                auto nameChargeM = Form("charge_full_vs_energy_foot_%i_asic_%i_mult_%i", iFoot, iAsic, m + 1);
                auto titleChargeM = Form("Charge full vs Energy for Foot %i, ASIC %i, mult %i", iFoot, iAsic, m + 1);
                chargeEnergyFootAsicMult[iFoot][iAsic][m] =
                    new TH2F(nameChargeM, titleChargeM, 200, 0, 15, 1000, 0, 5000);
            }
        }

        for (int iCharge = 0; iCharge < nCharges; ++iCharge)
        {
            auto namePosC = Form("pos_foot_%i_char_%i", iFoot, iCharge);
            auto titlePosC = Form("Energy vs Position for Foot %i, Charge %i", iFoot, iCharge);
            posFootChar[iFoot][iCharge] = new TH2F(namePosC, titlePosC, 640, 0, 640, 1000, 0, 5000);

            for (int iAsic = 0; iAsic < nAsics; ++iAsic)
            {
                auto nameEtaC = Form("eta_foot_%i_asic_%i_char_%i", iFoot, iAsic, iCharge);
                auto titleEtaC = Form("Eta vs Energy for Foot %i, ASIC %i, Charge %i", iFoot, iAsic, iCharge);
                etaFootAsicChar[iFoot][iAsic][iCharge] = new TH2F(nameEtaC, titleEtaC, 100, 0, 1, 1000, 0, 5000);
            }
        }
    }

    // Event loop
    int nEvents = evt->GetEntries();
    for (int iEvent = 0; iEvent < nEvents; ++iEvent)
    {
        if (iEvent % 10000 == 0 || iEvent == nEvents - 1)
        {
            float percent = 100.0f * iEvent / nEvents;
            printf("\rProcessing events... %.1f%%", percent);
            fflush(stdout);
        }

        evt->GetEntry(iEvent);
        int nHitFoot = footDataArray->GetEntries();
        int nHitLos = losDataArray->GetEntries();
        int nHitCal = footCalDataArray->GetEntries();

        if (nHitLos != 1)
            continue;

        auto *hitLos = (R3BLosHitData *)losDataArray->At(0);
        double charge_full = hitLos->GetZ();
        int charge = -1;

        if ((charge_full < 7.89 + 0.26) && (charge_full > 7.89 - 0.26))
            charge = 0;
        else if ((charge_full < 8.68 + 0.31) && (charge_full > 8.68 - 0.31))
            charge = 1;

        hChargeFull->Fill(charge_full);

        for (int iHit = 0; iHit < nHitFoot; ++iHit)
        {
            auto *hitFoot = (R3BFootHitData *)footDataArray->At(iHit);
            int footId = hitFoot->GetDetId() - 1;
            double eta = hitFoot->GetEta();
            double energy = hitFoot->GetEnergy();
            double pos = hitFoot->GetPos();
            int mult = hitFoot->GetMulStrip();

            pos += 50;
            pos *= 6.4;
            int asic = pos / 64;

            if (footId < 0 || footId >= nFoots || asic < 0 || asic >= nAsics)
                continue;

            etaFootAsic[footId][asic]->Fill(eta, energy);

            if (mult >= 1 && mult <= maxMult)
            {
                etaFootAsicMult[footId][asic][mult - 1]->Fill(eta, energy);
                posFootMult[footId][mult - 1][0]->Fill(pos, energy);
                chargeEnergyFootAsicMult[footId][asic][mult - 1]->Fill(charge_full, energy);
            }

            if (mult <= 1)
            {
                posFoot[footId]->Fill(pos, energy);
                chargeEnergyFootAsic[footId][asic]->Fill(charge_full, energy);
            }

            if (charge >= 0 && charge < nCharges)
            {
                etaFootAsicChar[footId][asic][charge]->Fill(eta, energy);

                if (mult <= 1)
                    posFootChar[footId][charge]->Fill(pos, energy);
            }
        }

        for (int i = 0; i < nHitCal; ++i)
        {
            auto *cal = (R3BFootCalData *)footCalDataArray->At(i);
            int footId = cal->GetDetId() - 1;
            int strip = cal->GetStripId();
            double energy = cal->GetEnergy();

            if (footId >= 0 && footId < nFoots && strip >= 0 && strip < nStrips)
                stripEnergyFoot[footId]->Fill(strip, energy);
        }
    }

    // Write output
    auto *fileOut = new TFile(outFilePath, "RECREATE");

    for (auto &footVec : etaFootAsicChar)
        for (auto &asicVec : footVec)
            for (auto &hist : asicVec)
                if (hist)
                    hist->Write();

    for (auto &footVec : posFootChar)
        for (auto &hist : footVec)
            if (hist)
                hist->Write();

    for (auto &footVec : etaFootAsic)
        for (auto &hist : footVec)
            if (hist)
                hist->Write();

    for (auto &footVec : etaFootAsicMult)
        for (auto &asicVec : footVec)
            for (auto *hist : asicVec)
                if (hist)
                    hist->Write();

    for (auto *hist : posFoot)
        if (hist)
            hist->Write();

    for (auto &footVec : chargeEnergyFootAsic)
        for (auto *hist : footVec)
            if (hist)
                hist->Write();

    for (auto *hist : stripEnergyFoot)
        if (hist)
            hist->Write();

    for (auto &footVec : posFootMult)
        for (auto &multVec : footVec)
            for (auto *hist : multVec)
                if (hist)
                    hist->Write();

    for (auto &footVec : chargeEnergyFootAsicMult)
        for (auto &asicVec : footVec)
            for (auto *hist : asicVec)
                if (hist)
                    hist->Write();

    if (hChargeFull)
        hChargeFull->Write();

    fileOut->Close();
    std::cout << "\nHistograms saved to " << outFilePath << std::endl;

    return 0;
}