#pragma once

#include <vector>
#include <array>
#include <utility>
#include <iostream>
#include <iomanip>
#include <cmath>

#include "TString.h"
#include "TFile.h"
#include "TChain.h"
#include "TROOT.h"
#include "TTreeReader.h"
#include "TTreeReaderArray.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TCanvas.h"

#include "../../utils/plotStyles.h"

class footAna
{
public:
    // === constructor para un solo fichero ===
    explicit footAna(const std::string &file) : fDataPaths({TString(file.c_str())})
    {
        initIOAndHists();
        buildChain();
    }

    // === constructor para varios ficheros ===
    explicit footAna(const std::vector<TString> &files) : fDataPaths(files)
    {
        initIOAndHists();
        buildChain();
    }

    ~footAna()
    {
        if (fResults)
        {
            fResults->Write();
            fResults->Close();
            delete fResults;
            fResults = nullptr;
        }

        if (fChain)
        {
            delete fChain;
            fChain = nullptr;
        }
    }

    void fillHistograms()
    {
        if (!fChain || fChain->GetEntries() == 0)
            return;

        TTreeReader reader(fChain);
        TTreeReaderArray<UChar_t> fDetId(reader, "FootHitData.fDetId");
        TTreeReaderArray<UShort_t> fMulStrip(reader, "FootHitData.fMulStrip");
        TTreeReaderArray<double> fZCharge(reader, "FootHitData.fZCharge");
        TTreeReaderArray<double> fEta(reader, "FootHitData.fEta");
        TTreeReaderArray<double> losZ(reader, "LosHit.fZ");
        TTreeReaderArray<double> tofdZ(reader, "TofdHit.fZCharge");
        TTreeReaderArray<UShort_t> footNbHit(reader, "FootHitData.fNbHit");
        TTreeReaderArray<double> footPos(reader, "FootHitData.fPos"); // cm
        TTreeReaderArray<UInt_t> planeTofd(reader, "TofdHit.fPlaneId");
        TTreeReaderArray<UInt_t> barTofd(reader, "TofdHit.fBarId");
        TTreeReaderArray<double> frsAq(reader, "FrsData.fAq");
        TTreeReaderArray<double> frsZ(reader, "FrsData.fZ");

        const Long64_t nTot = fChain->GetEntries();
        Long64_t iEv = 0;
        const Long64_t step = std::max<Long64_t>(1, nTot / 100);

        while (reader.Next())
        {
            ++iEv;
            if (iEv % step == 0 || iEv == nTot)
            {
                const int pct = static_cast<int>((100.0 * iEv) / std::max<Long64_t>(1, nTot) + 0.5);
                std::cout << "\rProcesing " << std::setw(3) << pct
                          << "% (" << iEv << "/" << nTot << ")" << std::flush;
            }

            if (frsAq.GetSize() != 1)
                continue;

            if (losZ.GetSize() != 1)
                continue;

            /*
            if (planeTofd.GetSize() != 4)
                continue;

            bool isGoodTofd = true;

            for (int i = 0; i < 4; i++)
                if (planeTofd[i] != i + 1)
                    isGoodTofd = false;

            if (!isGoodTofd)
                continue;
            */

            const double zLos = losZ[0];
            const int los_bin = losChargeBin(zLos);

            if (frsAq.GetSize() != 1)
                continue;
            if (frsZ.GetSize() != 1)
                continue;

            IncomingId->Fill(frsAq[0], frsZ[0]);

            /*

            if (fDetId.GetSize() != 8)
                continue;

            bool isGoodFoot = true;
            for (int i = 0; i < 8; i++)
                if (fDetId[i] != i + 1)
                {
                    isGoodFoot = false;
                    std::cout << iEv << std::endl;
                }


            if (!isGoodFoot)
                continue;

            */

            // --- loop FOOT hits ---

            if ((fZCharge[4] > 8.4) || (fZCharge[4] < 8))
                continue;

            for (int i = 0; i < fDetId.GetSize(); i++)
            {
                if (fMulStrip[i] <= 1)
                    continue;

                const int det = static_cast<int>(fDetId[i]);
                if (det < 1 || det > fNbDets)
                    continue;

                const double x_mm = footPos[i] * 10.0; // cm -> mm
                const double strips = (x_mm / 10. + 48.) / 96. * 640.;
                const int stripVal = static_cast<int>(strips);

                // ASIC value
                int asic = static_cast<int>(std::floor(strips / 64.0));
                if (asic < 0)
                    asic = 0;
                if (asic >= kNbAsics)
                    asic = kNbAsics - 1;

                chargeFootVsLos[det - 1]->Fill(fZCharge[i], zLos);
                energyVsEta[det - 1]->Fill(fEta[i], fZCharge[i]);
                FootpositionVsCharge[det - 1]->Fill(x_mm, fZCharge[i]);
                energyVsEtaAsic[det - 1][asic]->Fill(fEta[i], fZCharge[i]);
                if (los_bin >= 0)
                    energyVsEtaAsicLos[det - 1][asic][los_bin]->Fill(fEta[i], fZCharge[i]);

                for (int j = 0; j < planeTofd.GetSize(); j++)
                {
                    const int plane = static_cast<int>(planeTofd[j]); // 1..4
                    if (plane >= 1 && plane <= 4)
                    {
                        FootVsTofdCharge[det - 1][plane - 1]->Fill(tofdZ[j], fZCharge[i]);
                        if (tofdChargePlane[plane - 1])
                            tofdChargePlane[plane - 1]->Fill(tofdZ[j]);

                        if (tofdZ[j] > 8.33)
                            std::cout
                                << "[WARN] Wrong value? " << tofdZ[j] << " at evt = " << iEv << std::endl;

                        const int bar = static_cast<int>(barTofd[j]);
                        if (FootChargeVsLosBar[det - 1][plane - 1])
                        {
                            FootChargeVsLosBar[det - 1][plane - 1]->Fill(bar, fZCharge[i]);
                            // FootChargeVsLosBar[det - 1][plane - 1]->Fill(bar, tofdZ[j]);
                        }
                    }
                }
            }

            // correlaciones 1-2, 3-4, 5-6, 7-8
            for (int i = 0; i < fDetId.GetSize() - 1; i++)
            {
                const int deti = static_cast<int>(fDetId[i]);

                for (int j = i + 1; j < fDetId.GetSize(); j++)
                {
                    const int detj = static_cast<int>(fDetId[j]);

                    if (detj - deti != 1)
                        continue;
                    if ((detj - 1) / 2 != (deti - 1) / 2)
                        continue;

                    if (fMulStrip[i] > 0 && fMulStrip[j] > 0)
                        chargeCorrelationFoots[(deti - 1) / 2]->Fill(fZCharge[i], fZCharge[j]);
                }
            }
        }

        // --- escribir al TFile ---
        fResults->cd();

        for (auto &h : chargeFootVsLos)
            h->Write();
        for (auto &h : energyVsEta)
            h->Write();
        for (auto &h : chargeCorrelationFoots)
            h->Write();
        for (auto &h : FootpositionVsCharge)
            h->Write();

        for (auto &vec : FootVsTofdCharge)
            for (auto &h : vec)
                h->Write();

        for (auto &detvec : FootChargeVsLosBar)
            for (auto *h : detvec)
                if (h)
                    h->Write();

        for (auto &vec : energyVsEtaAsic)
            for (auto *h : vec)
                h->Write();

        for (auto &detvec : energyVsEtaAsicLos)
            for (auto &asicvec : detvec)
                for (auto *h : asicvec)
                    h->Write();

        for (auto *h : tofdChargePlane)
            if (h)
                h->Write();
    }

    void saveFigures(TString name = "results_")
    {
        double rightMarg = 0.15;

        auto draw_and_save = [&](TH2D *h)
        {
            setOpenGL();
            auto *c = new TCanvas("c", "", 800, 600);
            setCanvasStyle(c, rightMarg);
            setHistogramStyle(h, h->GetXaxis()->GetTitle(), h->GetYaxis()->GetTitle());
            h->Draw("zcol");
            gPad->SetLogz();
            c->SaveAs(TString(resultsPath + name + h->GetName()) + ".pdf");
            delete c;
        };

        auto draw_and_save_1d = [&](TH1D *h)
        {
            setOpenGL();
            auto *c = new TCanvas("c1d", "", 800, 600);
            setCanvasStyle(c, 0.05);
            setHistogramStyle(h, "Z_{ToFD}", "Counts #", kCyan + 3);
            h->Draw();
            c->SaveAs(TString(resultsPath + name + h->GetName()) + ".pdf");
            delete c;
        };

        for (auto &h : chargeFootVsLos)
            draw_and_save(h);
        draw_and_save(IncomingId);
        for (auto &h : chargeCorrelationFoots)
            draw_and_save(h);
        for (auto &h : FootpositionVsCharge)
            draw_and_save(h);
        for (auto &vec : FootVsTofdCharge)
            for (auto &h : vec)
                draw_and_save(h);

        for (auto &detvec : FootChargeVsLosBar)
            for (auto *h : detvec)
                if (h)
                    draw_and_save(h);

        for (auto *h : tofdChargePlane)
            if (h)
                draw_and_save_1d(h);
    }

private:
    void initIOAndHists()
    {
        resultsPath = static_cast<TString>(getenv("repopath")) + "/results/cal/";
        fResults = new TFile(resultsPath + "resultsFCuts.root", "RECREATE");

        // Parámetros de binning
        const int nbx = 300; // bins eje X (Z)
        const double xmin = 0., xmax = 15.;
        const int nby = 600; // bins eje Y (Z)
        const double ymin = 0., ymax = 15.;

        // Número de paddles por plano (ajusta si conoces el real)
        const int nBarsGuess = 100;
        const double barMin = 0.5, barMax = nBarsGuess + 0.5; // centrado en enteros

        chargeFootVsLos.reserve(fNbDets);
        energyVsEta.reserve(fNbDets);
        FootpositionVsCharge.reserve(fNbDets);
        FootVsTofdCharge.resize(fNbDets);
        FootChargeVsLosBar.resize(fNbDets); // NUEVO
        energyVsEtaAsic.resize(fNbDets);
        energyVsEtaAsicLos.resize(fNbDets);

        IncomingId = new TH2D("h2_IncomingId",
                              "Incoming ID (FRS);A/Q_{FRS};Z_{LOS}",
                              1000, 2.5, 3.0, 1000, 4.0, 15.0);

        for (int p = 0; p < 4; ++p)
        {
            tofdChargePlane[p] = new TH1D(
                Form("h1_TofdCharge_plane%d", p + 1),
                Form("ToFD Plane %d;Z_{ToFD};Counts", p + 1),
                nbx, xmin, xmax);
        }

        // --- Por detector ---
        for (int det = 1; det <= fNbDets; det++)
        {
            chargeFootVsLos.push_back(new TH2D(
                Form("h2_FootVsLos_z_det%d", det),
                Form("Det %d;Z_{FOOT};Z_{LOS}", det),
                nbx, xmin, xmax, nby, ymin, ymax));

            energyVsEta.push_back(new TH2D(
                Form("h2_EnergyVsEta_det%d", det),
                Form("Det %d;#eta;Z_{FOOT}", det),
                nbx, 0, 1, nbx, xmin, xmax));

            FootpositionVsCharge.push_back(new TH2D(
                Form("h2_ChargeVsPos_det%d", det),
                Form("Det %d;x [mm];Z_{FOOT}", det),
                nbx, -320, 320, nbx, xmin, xmax));

            // 4 planos ToFD: Z_ToFD vs Z_FOOT
            FootVsTofdCharge[det - 1].resize(4);
            for (int j = 0; j < 4; j++)
            {
                FootVsTofdCharge[det - 1][j] = new TH2D(
                    Form("h2_FootVsTofd_det%d_plane%d", det, j + 1),
                    Form("Det %d - Plane %d;Z_{ToFD %d};Z_{FOOT %d}",
                         det, j + 1, j + 1, det),
                    nbx, xmin, xmax, nbx, xmin, xmax);
            }

            // NUEVO: barId (paddle) vs Z_FOOT por plano
            FootChargeVsLosBar[det - 1].resize(4);
            for (int j = 0; j < 4; ++j)
            {
                FootChargeVsLosBar[det - 1][j] = new TH2D(
                    Form("h2_FootChargeVsLosBar_det%d_plane%d", det, j + 1),
                    Form("Det %d vs LOS plane %d;Paddle (barId);Z_{FOOT %d}",
                         det, j + 1, det),
                    nBarsGuess, barMin, barMax, nbx, xmin, xmax);
            }

            // Por ASIC
            energyVsEtaAsic[det - 1].resize(kNbAsics, nullptr);
            for (int a = 0; a < kNbAsics; ++a)
            {
                energyVsEtaAsic[det - 1][a] = new TH2D(
                    Form("h2_EnergyVsEta_det%d_asic%d", det, a),
                    Form("Det %d, ASIC %d;#eta;Z_{FOOT}", det, a),
                    nbx, 0, 1, nbx, xmin, xmax);
            }

            // Por ASIC y bin de LOS
            energyVsEtaAsicLos[det - 1].resize(kNbAsics);
            for (int a = 0; a < kNbAsics; ++a)
            {
                energyVsEtaAsicLos[det - 1][a].resize(kNbLosBins, nullptr);
                for (int b = 0; b < kNbLosBins; ++b)
                {
                    const int losLabel = fLosBinLabels[b];
                    energyVsEtaAsicLos[det - 1][a][b] = new TH2D(
                        Form("h2_EnergyVsEta_det%d_asic%d_LOS%d", det, a, losLabel),
                        Form("Det %d, ASIC %d, LOS %d;#eta;Z_{FOOT}", det, a, losLabel),
                        nbx, 0, 1, nbx, xmin, xmax);
                }
            }
        }

        // Correlaciones (1-2, 3-4, 5-6, 7-8)
        for (int det = 1; det < fNbDets; det += 2)
        {
            chargeCorrelationFoots.push_back(new TH2D(
                Form("h2_ChargeCorrelation_det%d_det%d", det, det + 1),
                Form("Dets %d and %d;Z_{FOOT %d};Z_{FOOT %d}",
                     det, det + 1, det, det + 1),
                nbx, xmin, xmax, nbx, xmin, xmax));
        }
    }

    void buildChain()
    {
        fChain = new TChain("evt");
        for (const auto &p : fDataPaths)
            fChain->Add(p);
    }

    int losChargeBin(double z) const
    {
        for (int i = 0; i < kNbLosBins; ++i)
            if (z >= fLosBins[i].first && z <= fLosBins[i].second)
                return i;
        return -1;
    }

private:
    // --- configuración ---
    const int fNbDets = 8;
    static constexpr int kNbAsics = 11;
    static constexpr int kNbLosBins = 6;

    const std::array<std::pair<double, double>, kNbLosBins> fLosBins{
        std::make_pair(5.4, 6.3),
        std::make_pair(6.6, 7.4),
        std::make_pair(7.5, 8.5),
        std::make_pair(8.5, 9.5),
        std::make_pair(9.6, 10.2),
        std::make_pair(10.6, 11.4)};
    const std::array<int, kNbLosBins> fLosBinLabels{6, 7, 8, 9, 10, 11};

    // --- IO ---
    std::vector<TString> fDataPaths;
    TChain *fChain = nullptr;
    TFile *fResults = nullptr;
    TString resultsPath;

    // --- histos ---
    TH2D *IncomingId = nullptr;

    std::vector<TH2D *> chargeFootVsLos;                              // [det]
    std::vector<TH2D *> energyVsEta;                                  // [det]
    std::vector<TH2D *> chargeCorrelationFoots;                       // [(1-2),(3-4),(5-6),(7-8)]
    std::vector<std::vector<TH2D *>> FootVsTofdCharge;                // [det][plane]
    std::vector<TH2D *> FootpositionVsCharge;                         // [det]
    std::vector<std::vector<TH2D *>> energyVsEtaAsic;                 // [det][asic]
    std::vector<std::vector<std::vector<TH2D *>>> energyVsEtaAsicLos; // [det][asic][losBin]

    // NUEVO: barId (paddle) vs Z_FOOT por plano
    std::vector<std::vector<TH2D *>> FootChargeVsLosBar; // [det][plane]

    // ToFD charge por plano
    std::array<TH1D *, 4> tofdChargePlane{}; // [plane 0..3]
};