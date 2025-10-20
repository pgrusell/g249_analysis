#include <vector>
#include <array>   // NEW
#include <utility> // NEW
#include "../../utils/plotStyles.h"

class footAna
{
public:
    // === constructor for one file ===
    footAna(const std::string &file) : fDataPaths({TString(file.c_str())})
    {
        initIOAndHists();
        buildChain();
    }

    // === constructor for many files ===
    footAna(const std::vector<TString> &files) : fDataPaths(files)
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
        TTreeReaderArray<double> footPos(reader, "FootHitData.fPos"); // cm
        TTreeReaderArray<UInt_t> planeTofd(reader, "TofdHit.fPlaneId");

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

            if (losZ.GetSize() != 1)
                continue;

            const double zLos = losZ[0];
            const int los_bin = losChargeBin(zLos);

            if (fDetId.GetSize() == 0)
                continue;

            for (int i = 0; i < fDetId.GetSize(); i++)
            {
                if (fMulStrip[i] <= 0)
                    continue;

                const int det = static_cast<int>(fDetId[i]);
                if (det < 1 || det > fNbDets)
                    continue;

                const double x_mm = footPos[i] * 10.0; // cm -> mm
                const double strips = (x_mm / 10. + 48.) / 96. * 640;

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

                // ToFD
                for (int j = 0; j < planeTofd.GetSize(); j++)
                {
                    const int plane = static_cast<int>(planeTofd[j]);
                    if (plane >= 1 && plane <= 4)
                        FootVsTofdCharge[det - 1][plane - 1]->Fill(tofdZ[j], fZCharge[i]);
                }
            }

            // correlations 1-2, 3-4, 5-6, 7-8
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

        // write
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
        for (auto &vec : energyVsEtaAsic)
            for (auto *h : vec)
                h->Write();
        for (auto &detvec : energyVsEtaAsicLos)
            for (auto &asicvec : detvec)
                for (auto *h : asicvec)
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
            c->SaveAs(TString(resultsPath + name + h->GetName()) + ".pdf");
            delete c;
        };

        for (auto &h : chargeFootVsLos)
            draw_and_save(h);
        for (auto &h : energyVsEta)
        {
            setOpenGL();
            auto *c = new TCanvas("c", "", 800, 600);
            setCanvasStyle(c, rightMarg);
            setHistogramStyle(h, "#eta", "Z_{FOOT}");
            h->Draw("zcol");
            gPad->SetLogz();
            c->SaveAs(TString(resultsPath + name + h->GetName()) + ".pdf");
            delete c;
        }
        for (auto &h : chargeCorrelationFoots)
            draw_and_save(h);
        for (auto &h : FootpositionVsCharge)
            draw_and_save(h);
        for (auto &vec : FootVsTofdCharge)
            for (auto &h : vec)
                draw_and_save(h);
        for (auto &vec : energyVsEtaAsic)
            for (auto *h : vec)
            {
                setOpenGL();
                auto *c = new TCanvas("c", "", 800, 600);
                setCanvasStyle(c, rightMarg);
                setHistogramStyle(h, "#eta", "Z_{FOOT}");
                h->Draw("zcol");
                gPad->SetLogz();
                c->SaveAs(TString(resultsPath + name + h->GetName()) + ".pdf");
                delete c;
            }

        for (auto &detvec : energyVsEtaAsicLos)
            for (auto &asicvec : detvec)
                for (auto *h : asicvec)
                {
                    setOpenGL();
                    auto *c = new TCanvas("c", "", 800, 600);
                    setCanvasStyle(c, rightMarg);
                    setHistogramStyle(h, "#eta", "Z_{FOOT}");
                    h->Draw("zcol");
                    gPad->SetLogz();
                    c->SaveAs(TString(resultsPath + name + h->GetName()) + ".pdf");
                    delete c;
                }
    }

private:
    void initIOAndHists()
    {
        resultsPath = static_cast<TString>(getenv("repopath")) + "/results/cal/";
        fResults = new TFile(resultsPath + "resultsWCorrv2.root", "RECREATE");

        // Histograms
        const int nbx = 300;
        // const double xmin = 0., xmax = 5000.;
        const double xmin = 0., xmax = 15.;
        const int nby = 600;
        // const double ymin = 0., ymax = 5000.;
        const double ymin = 0., ymax = 15.;

        chargeFootVsLos.reserve(fNbDets);
        energyVsEta.reserve(fNbDets);
        FootpositionVsCharge.reserve(fNbDets);
        FootVsTofdCharge.resize(fNbDets);
        energyVsEtaAsic.resize(fNbDets);
        energyVsEtaAsicLos.resize(fNbDets);

        for (int det = 1; det <= fNbDets; det++)
        {
            chargeFootVsLos.push_back(new TH2D(Form("h2_FootVsLos_z_det%d", det),
                                               Form("Det %d;Z_{FOOT};Z_{LOS}", det),
                                               nbx, xmin, xmax, nby, ymin, ymax));

            energyVsEta.push_back(new TH2D(Form("h2_EnergyVsEta_det%d", det),
                                           Form("Det %d;#eta;Z_{FOOT}", det),
                                           nbx, 0, 1, nbx, xmin, xmax));

            FootpositionVsCharge.push_back(new TH2D(Form("h2_ChargeVsPos_det%d", det),
                                                    Form("Det %d;x [mm];Z_{FOOT}", det),
                                                    nbx, -320, 320, nbx, xmin, xmax));

            // 4 tofd planes
            for (int j = 0; j < 4; j++)
            {
                FootVsTofdCharge[det - 1].push_back(new TH2D(
                    Form("h2_FootVsTofd_det%d_plane%d", det, j + 1),
                    Form("Det %d - Plane %d;Z_{ToFD %d};Z_{FOOT %d}", det, j + 1, j + 1, det),
                    nbx, xmin, xmax, nbx, xmin, xmax));
            }

            energyVsEtaAsic[det - 1].resize(kNbAsics, nullptr);
            for (int a = 0; a < kNbAsics; ++a)
            {
                energyVsEtaAsic[det - 1][a] = new TH2D(
                    Form("h2_EnergyVsEta_det%d_asic%d", det, a),
                    Form("Det %d, ASIC %d;#eta;Z_{FOOT}", det, a),
                    nbx, 0, 1, nbx, xmin, xmax);
            }

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

        // correlations (1-2, 3-4, 5-6, 7-8)
        for (int det = 1; det < fNbDets; det += 2)
        {
            chargeCorrelationFoots.push_back(new TH2D(
                Form("h2_ChargeCorrelation_det%d_det%d", det, det + 1),
                Form("Dets %d and %d;Z_{FOOT %d};Z_{FOOT %d}", det, det + 1, det, det + 1),
                nbx, xmin, xmax, nbx, xmin, xmax));
        }
    }

    void buildChain()
    {
        fChain = new TChain("evt");
        for (const auto &p : fDataPaths)
        {
            fChain->Add(p);
        }
    }

    // ---------- NEW: helpers ----------

    int losChargeBin(double z) const
    {
        for (int i = 0; i < kNbLosBins; ++i)
        {
            if (z >= fLosBins[i].first && z <= fLosBins[i].second)
                return i;
        }
        return -1;
    }

private:
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

    std::vector<TString> fDataPaths;
    TChain *fChain = nullptr;

    TFile *fResults = nullptr;
    TString resultsPath;

    std::vector<TH2D *> chargeFootVsLos;
    std::vector<TH2D *> energyVsEta;
    std::vector<TH2D *> chargeCorrelationFoots;
    std::vector<std::vector<TH2D *>> FootVsTofdCharge;
    std::vector<TH2D *> FootpositionVsCharge;
    std::vector<std::vector<TH2D *>> energyVsEtaAsic;                 // [det][asic]
    std::vector<std::vector<std::vector<TH2D *>>> energyVsEtaAsicLos; // [det][asic][losBin]
};