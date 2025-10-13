#include <vector>
#include "../../utils/plotStyles.h"

class footAna
{
public:
    footAna(const std::string &file) : fDataPath(file)
    {

        resultsPath = static_cast<TString>(getenv("repopath")) + "/results/cal/";

        fResults = new TFile(resultsPath + "results.root", "RECREATE");
        fDataFile = TFile::Open(fDataPath, "READ");

        // Histograms
        const int nbx = 120;
        const double xmin = 5., xmax = 14.;
        const int nby = 600;
        const double ymin = 5., ymax = 14.;

        for (int det = 1; det <= fNbDets; det++)
        {
            chargeFootVsLos.push_back(new TH2D(Form("h2_FootVsLos_z_det%d", det),
                                               Form("Det %d;Z_{FOOT};Z_{LOS}", det),
                                               nbx, xmin, xmax, nby, ymin, ymax));

            energyVsEta.push_back(new TH2D(Form("h2_EnergyVsEta_det%d", det),
                                           Form("Det %d;#eta;Z_{FOOT}", det),
                                           nbx, 0, 1, nbx, xmin, xmax));
        }

        for (int det = 1; det < fNbDets; det += 2)
        {
            chargeCorrelationFoots.push_back(new TH2D(Form("h2_ChargeCorrelation_det%d_det%d", det, det + 1),
                                                      Form("Dets %d and %d;Z_{FOOT %d};Z_{FOOT %d}", det, det + 1, det, det + 1),
                                                      nbx, xmin, xmax, nbx, xmin, xmax));
        }
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
    }

    void fillHistograms()
    {
        TTree *tree = (TTree *)fDataFile->Get("evt");

        TTreeReader reader(tree);
        TTreeReaderArray<UChar_t> fDetId(reader, "FootHitData.fDetId");
        TTreeReaderArray<UShort_t> fMulStrip(reader, "FootHitData.fMulStrip");
        TTreeReaderArray<double> fZCharge(reader, "FootHitData.fZCharge");
        TTreeReaderArray<double> fEta(reader, "FootHitData.fEta");
        TTreeReaderArray<double> losZ(reader, "LosHit.fZ");

        // Event loop
        while (reader.Next())
        {
            if (losZ.GetSize() != 1)
                continue;

            double zLos = losZ[0];

            if (fDetId.GetSize() == 0)
                continue;

            for (int i = 0; i < fDetId.GetSize(); i++)
            {
                if (fMulStrip[i] > 0)
                {
                    int det = fDetId[i];
                    if (det >= 1 && det <= 8)
                    {
                        chargeFootVsLos[det - 1]->Fill(fZCharge[i], zLos);
                        energyVsEta[det - 1]->Fill(fEta[i], fZCharge[i]);
                    }
                }
            }

            // Fill correlations

            for (int i = 0; i < fDetId.GetSize() - 1; i++)
            {
                int deti = fDetId[i];

                for (int j = i + 1; j < fDetId.GetSize(); j++)
                {

                    int detj = fDetId[j];

                    if (detj - deti != 1)
                        continue;

                    if ((detj - 1) / 2 != (deti - 1) / 2)
                        continue;

                    chargeCorrelationFoots[(deti - 1) / 2]->Fill(fZCharge[i], fZCharge[j]);
                }
            }
        }

        fResults->cd();

        for (auto &h : chargeFootVsLos)
            h->Write();
        for (auto &h : energyVsEta)
            h->Write();
        for (auto &h : chargeCorrelationFoots)
            h->Write();
    }

    void saveFigures(TString name = "results_")
    {

        double rightMarg = 0.15;

        for (auto &h : chargeFootVsLos)
        {
            setOpenGL();
            auto *c = new TCanvas("c", "", 800, 600);
            setCanvasStyle(c, rightMarg);
            setHistogramStyle(h, "Z_{FOOT}", "Z_{LOS}");
            h->Draw("zcol");
            c->SaveAs(TString(resultsPath + name + h->GetName()) + ".pdf");
            delete c;
        }

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
        {
            setOpenGL();
            auto *c = new TCanvas("c", "", 800, 600);
            setCanvasStyle(c, rightMarg);
            setHistogramStyle(h, h->GetXaxis()->GetTitle(), h->GetYaxis()->GetTitle());
            h->Draw("zcol");
            c->SaveAs(TString(resultsPath + name + h->GetName()) + ".pdf");
            delete c;
        }
    }

private:
    const int fNbDets = 8;

    TString fDataPath;
    TFile *fResults = nullptr;
    TFile *fDataFile = nullptr;

    TString resultsPath;

    std::vector<TH2D *> chargeFootVsLos;
    std::vector<TH2D *> energyVsEta;
    std::vector<TH2D *> chargeCorrelationFoots;
};