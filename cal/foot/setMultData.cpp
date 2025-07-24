// This macro analyzes the data from a losHit, footHit and footCal data file, to fill the following histograms:
// :etaFootAsicChar: eta vs energy. One histogram per Foot, asic and charge (8 and 9)
// :etaFootAsic: eta vs energy. One histogram per Foot, asic
// :etaFootAsicMult: eta vs energy. One histogram per Foot, asic and multiplicity
// :posFootChar: energy vs position. One histogram per FOOT and charge in LOS
// :posFootMult: energy vs position. One histogram per FOOT and multiplicity
// :chargeEnergyFootAsic: charge correlation FOOT vs LOS. One histogram per FOOT and ASIC
// :chargeEnergyFootAsicMult: charge correlation FOOT vs LOS. One histogram per FOOT, ASIC and multiplicity
// :hChargeFull: Energy spectra in LOS
// :etaFootAsicMultChar: eta vs energy. One histogram per Foot, Asic, Mult, and Charge
// :energyFootAsicMultChar: energy in FOOT. One 1D histogram per Foot, Asic, Mult, and Charge
// :energyVsEventFootAsicMultChar: energy vs event number. One histogram 2D per Foot, Asic, Mult, and Charge

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
#include <map>

int main(int argc, char **argv)
{
    TString inFilePath = "/nucl_lustre/pablogrusell/g249/g249_analysis/root_files/etaCorrection/g249_data_incoming_online_20250728_183008.root";
    TString outFilePath = "/nucl_lustre/pablogrusell/g249/g249_analysis/root_files/etaCorrection/foot_etaplot_mult2.root";

    if (argc > 1)
        inFilePath = argv[1];
    if (argc > 2)
        outFilePath = argv[2];

    const int nFoots = 8;
    const int nAsics = 5;
    const int maxMult = 4;
    const int nStrips = 640;
    const std::vector<int> chargesZ = {8, 9};
    const int nCharges = chargesZ.size();

    std::map<int, int> zToIndex = {{8, 0}, {9, 1}};
    std::map<int, std::pair<float, float>> zRanges = {
        {8, {7.4, 8.3}},
        {9, {8.32, 9.2}}};

    int nProcessedEvents = 0;

    auto *inFile = new TFile(inFilePath, "READ");
    if (!inFile || inFile->IsZombie())
        return 1;

    auto *evt = (TTree *)inFile->Get("evt");
    if (!evt)
        return 1;

    int nEvents = evt->GetEntries() / 10;

    evt->SetBranchStatus("*", 0);
    evt->SetBranchStatus("FootHitData*", 1);
    evt->SetBranchStatus("FootCalData*", 1);
    evt->SetBranchStatus("LosHit*", 1);

    evt->SetCacheSize(100 * 1024 * 1024);
    evt->AddBranchToCache("FootHitData", true);
    evt->AddBranchToCache("FootCalData", true);
    evt->AddBranchToCache("LosHit", true);

    auto *footCalDataArray = new TClonesArray("R3BFootCalData");
    auto *footDataArray = new TClonesArray("R3BFootHitData");
    auto *losDataArray = new TClonesArray("R3BLosHitData");

    evt->SetBranchAddress("FootCalData", &footCalDataArray);
    evt->SetBranchAddress("FootHitData", &footDataArray);
    evt->SetBranchAddress("LosHit", &losDataArray);

    std::vector<std::vector<std::vector<TH2F *>>> etaFootAsicChar(nFoots, std::vector<std::vector<TH2F *>>(nAsics, std::vector<TH2F *>(nCharges)));
    std::vector<std::vector<TH2F *>> etaFootAsic(nFoots, std::vector<TH2F *>(nAsics));
    std::vector<std::vector<std::vector<TH2F *>>> etaFootAsicMult(nFoots, std::vector<std::vector<TH2F *>>(nAsics, std::vector<TH2F *>(maxMult)));
    std::vector<std::vector<std::vector<std::vector<TH2F *>>>> etaFootAsicMultChar(nFoots, std::vector<std::vector<std::vector<TH2F *>>>(nAsics, std::vector<std::vector<TH2F *>>(maxMult, std::vector<TH2F *>(nCharges))));
    std::vector<std::vector<std::vector<std::vector<TH1F *>>>> energyFootAsicMultChar(nFoots, std::vector<std::vector<std::vector<TH1F *>>>(nAsics, std::vector<std::vector<TH1F *>>(maxMult, std::vector<TH1F *>(nCharges))));
    std::vector<std::vector<std::vector<std::vector<TH2F *>>>> energyVsEventFootAsicMultChar(
        nFoots, std::vector<std::vector<std::vector<TH2F *>>>(
                    nAsics, std::vector<std::vector<TH2F *>>(
                                maxMult, std::vector<TH2F *>(nCharges))));
    std::vector<std::vector<TH2F *>> posFootChar(nFoots, std::vector<TH2F *>(nCharges));
    std::vector<TH2F *> posFoot(nFoots);
    std::vector<std::vector<std::vector<TH2F *>>> posFootMult(nFoots, std::vector<std::vector<TH2F *>>(maxMult, std::vector<TH2F *>(1)));
    std::vector<std::vector<TH2F *>> chargeEnergyFootAsic(nFoots, std::vector<TH2F *>(nAsics));
    std::vector<std::vector<std::vector<TH2F *>>> chargeEnergyFootAsicMult(nFoots, std::vector<std::vector<TH2F *>>(nAsics, std::vector<TH2F *>(maxMult)));
    std::vector<TH2F *> stripEnergyFoot(nFoots);
    TH1F *hChargeFull = new TH1F("charge_full", "Charge full spectrum", 150, 0, 15);

    for (int f = 0; f < nFoots; ++f)
    {
        posFoot[f] = new TH2F(Form("pos_foot_%d", f), "", 640, 0, 640, 200, 0, 15);
        stripEnergyFoot[f] = new TH2F(Form("strip_vs_energy_foot_%d", f), "", nStrips, 0, nStrips, 200, 0, 5000);
        for (int m = 0; m < maxMult; ++m)
            posFootMult[f][m][0] = new TH2F(Form("pos_foot_%d_mult_%d", f, m + 1), "", 640, 0, 640, 200, 0, 15);
        for (int a = 0; a < nAsics; ++a)
        {
            etaFootAsic[f][a] = new TH2F(Form("eta_foot_%d_asic_%d", f, a), "", 100, 0, 1, 200, 0, 5000);
            chargeEnergyFootAsic[f][a] = new TH2F(Form("charge_full_vs_energy_foot_%d_asic_%d", f, a), "", 200, 0, 15, 200, 0, 15);
            for (int m = 0; m < maxMult; ++m)
            {
                etaFootAsicMult[f][a][m] = new TH2F(Form("eta_foot_%d_asic_%d_mult_%d", f, a, m + 1), "", 100, 0, 1, 200, 0, 5000);
                chargeEnergyFootAsicMult[f][a][m] = new TH2F(Form("charge_full_vs_energy_foot_%d_asic_%d_mult_%d", f, a, m + 1), "", 200, 0, 15, 200, 0, 15);
                for (int c = 0; c < nCharges; ++c)
                {
                    etaFootAsicMultChar[f][a][m][c] = new TH2F(Form("eta_foot_%d_asic_%d_mult_%d_charge_%d", f, a, m + 1, chargesZ[c]), "", 100, 0, 1, 200, 0, 15);
                    energyFootAsicMultChar[f][a][m][c] = new TH1F(Form("energy_foot_%d_asic_%d_mult_%d_charge_%d", f, a, m + 1, chargesZ[c]), "", 200, 0, 15);
                    energyVsEventFootAsicMultChar[f][a][m][c] = new TH2F(
                        Form("energy_vs_event_foot_%d_asic_%d_mult_%d_charge_%d", f, a, m + 1, chargesZ[c]),
                        "Energy vs Event;Event number;Energy [a.u.]", nEvents, 0, nEvents, 200, 0, 15);
                }
            }
            for (int c = 0; c < nCharges; ++c)
                etaFootAsicChar[f][a][c] = new TH2F(Form("eta_foot_%d_asic_%d_char_%d", f, a, chargesZ[c]), "", 100, 0, 1, 200, 0, 5000);
        }
        for (int c = 0; c < nCharges; ++c)
            posFootChar[f][c] = new TH2F(Form("pos_foot_%d_char_%d", f, chargesZ[c]), "", 640, 0, 640, 200, 0, 15);
    }

    for (int iEvent = 0; iEvent < nEvents; ++iEvent)
    {
        if (iEvent % 10000 == 0 || iEvent == nEvents - 1)
        {
            float percent = 100.0f * iEvent / nEvents;
            printf("\rProcessing events... %.1f%%", percent);
            fflush(stdout);
        }

        evt->GetEntry(iEvent);
        if (losDataArray->GetEntries() != 1)
            continue;

        auto *hitLos = (R3BLosHitData *)losDataArray->At(0);
        double charge_full = hitLos->GetZ();
        int chargeIdx = -1;
        for (auto &[z, range] : zRanges)
        {
            if (charge_full >= range.first && charge_full <= range.second)
            {
                chargeIdx = zToIndex[z];
                break;
            }
        }

        hChargeFull->Fill(charge_full);
        nProcessedEvents++;

        for (int iHit = 0; iHit < footDataArray->GetEntries(); ++iHit)
        {
            auto *hitFoot = (R3BFootHitData *)footDataArray->At(iHit);
            int foot = hitFoot->GetDetId() - 1;
            double eta = hitFoot->GetEta();
            double energy = hitFoot->GetEnergy();
            double footCharge = hitFoot->GetZCharge();
            double pos = (hitFoot->GetPos() + 50) * 6.4;
            int asic = pos / 64 - 3;

            int mult = hitFoot->GetMulStrip();

            if (foot < 0 || foot >= nFoots || asic < 0 || asic >= nAsics || mult < 1 || mult > maxMult)
                continue;

            etaFootAsic[foot][asic]->Fill(eta, energy);
            etaFootAsicMult[foot][asic][mult - 1]->Fill(eta, energy);
            posFootMult[foot][mult - 1][0]->Fill(pos, footCharge);
            chargeEnergyFootAsicMult[foot][asic][mult - 1]->Fill(charge_full, footCharge);

            if (chargeIdx >= 0)
            {
                etaFootAsicChar[foot][asic][chargeIdx]->Fill(eta, energy);
                etaFootAsicMultChar[foot][asic][mult - 1][chargeIdx]->Fill(eta, energy);
                energyFootAsicMultChar[foot][asic][mult - 1][chargeIdx]->Fill(footCharge);
                energyVsEventFootAsicMultChar[foot][asic][mult - 1][chargeIdx]->Fill(iEvent, footCharge);
                if (mult == 1)
                    posFootChar[foot][chargeIdx]->Fill(pos, footCharge);
            }

            posFoot[foot]->Fill(pos, footCharge);
            chargeEnergyFootAsic[foot][asic]->Fill(charge_full, footCharge);
        }

        for (int i = 0; i < footCalDataArray->GetEntries(); ++i)
        {
            auto *cal = (R3BFootCalData *)footCalDataArray->At(i);
            int foot = cal->GetDetId() - 1;
            int strip = cal->GetStripId();
            if (foot >= 0 && foot < nFoots && strip >= 0 && strip < nStrips)
                stripEnergyFoot[foot]->Fill(strip, cal->GetEnergy());
        }
    }

    auto *outFile = new TFile(outFilePath, "RECREATE");

    auto write = [](auto &container)
    {
        for (auto &v1 : container)
            for (auto &v2 : v1)
                for (auto &v3 : v2)
                    for (auto &obj : v3)
                        if (obj)
                            obj->Write();
    };

    write(etaFootAsicMultChar);
    write(energyFootAsicMultChar);
    write(energyVsEventFootAsicMultChar);

    for (auto &v : etaFootAsic)
        for (auto *h : v)
            if (h)
                h->Write();

    for (auto &v : etaFootAsicMult)
        for (auto &m : v)
            for (auto *h : m)
                if (h)
                    h->Write();

    for (auto &v : etaFootAsicChar)
        for (auto &m : v)
            for (auto *h : m)
                if (h)
                    h->Write();

    for (auto &v : posFootChar)
        for (auto *h : v)
            if (h)
                h->Write();

    for (auto &v : posFootMult)
        for (auto &m : v)
            for (auto *h : m)
                if (h)
                    h->Write();

    for (auto *h : posFoot)
        if (h)
            h->Write();

    for (auto &v : chargeEnergyFootAsic)
        for (auto *h : v)
            if (h)
                h->Write();

    for (auto &v : chargeEnergyFootAsicMult)
        for (auto &m : v)
            for (auto *h : m)
                if (h)
                    h->Write();

    for (auto *h : stripEnergyFoot)
        if (h)
            h->Write();

    if (hChargeFull)
        hChargeFull->Write();

    outFile->Close();

    std::cout << "\nProcessed events: " << nProcessedEvents << std::endl;
    std::cout << "Histograms saved to " << outFilePath << std::endl;

    return 0;
}