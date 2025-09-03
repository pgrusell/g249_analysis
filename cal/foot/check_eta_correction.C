#include <TFile.h>
#include <TTree.h>
#include <TTreeReader.h>
#include <TTreeReaderArray.h>
#include <TH2D.h>
#include <vector>
#include <string>
#include <limits>

void check_eta_correction(TString inFile = "g249_data_incoming_online_20250728_183008.root")
{
    TString repopath = gSystem->Getenv("repopath");

    // Open file
    TFile *fin = TFile::Open(repopath + "/data/" + inFile);
    TTree *tree = (TTree *)fin->Get("evt");

    // Read branches
    TTreeReader reader(tree);
    TTreeReaderArray<UChar_t> fDetId(reader, "FootHitData.fDetId");
    TTreeReaderArray<UShort_t> fMulStrip(reader, "FootHitData.fMulStrip");
    TTreeReaderArray<double> fZCharge(reader, "FootHitData.fZCharge");
    TTreeReaderArray<double> fEta(reader, "FootHitData.fEta");
    TTreeReaderArray<double> losZ(reader, "LosHit.fZ");

    // Histograms
    const int nbx = 120;
    const double xmin = 5., xmax = 14.;
    const int nby = 600;
    const double ymin = 5., ymax = 14.;

    std::vector<TH2D *> chargeFootVsLos;
    std::vector<TH2D *> energyVsEta;
    for (int det = 1; det <= 8; det++)
    {
        chargeFootVsLos.push_back(new TH2D(Form("h2_FootVsLos_z_det%d", det),
                                           Form("Det %d;Z_{FOOT};Z_{LOS}", det),
                                           nbx, xmin, xmax, nby, ymin, ymax));

        energyVsEta.push_back(new TH2D(Form("h2_EnergyVsEta_det%d", det),
                                       Form("Det %d;#eta;Z_{FOOT}", det),
                                       nbx, 0, 1, nbx, xmin, xmax));
    }

    // event loop
    while (reader.Next())
    {
        if (losZ.GetSize() != 1)
            continue;
        double zLos = losZ[0];

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
    }

    // save
    TFile out(repopath + "/results/cal/chargeId_perDet.root", "RECREATE");
    for (auto h : chargeFootVsLos)
        h->Write();
    for (auto h : energyVsEta)
        h->Write();
    out.Close();

    fin->Close();
}