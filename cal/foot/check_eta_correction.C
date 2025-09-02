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
    if (!fin || fin->IsZombie())
    {
        Error("check_eta_correction", "No pude abrir el fichero %s", inFile.Data());
        return;
    }

    TTree *tree = (TTree *)fin->Get("evt");
    if (!tree)
    {
        Error("check_eta_correction", "No existe el TTree 'evt'.");
        return;
    }

    // Read branches
    TTreeReader reader(tree);
    TTreeReaderArray<UChar_t> fDetId(reader, "FootHitData.fDetId");
    TTreeReaderArray<UShort_t> fMulStrip(reader, "FootHitData.fMulStrip");
    TTreeReaderArray<double> fZCharge(reader, "FootHitData.fZCharge");
    TTreeReaderArray<double> losZ(reader, "LosHit.fZ");

    // Histograms
    const int nbx = 120;
    const double xmin = 5., xmax = 14.;
    const int nby = 600;
    const double ymin = 5., ymax = 14.;

    std::vector<TH2D *> hists;
    for (int det = 1; det <= 8; det++)
    {
        hists.push_back(new TH2D(Form("h2_FootVsLos_z_det%d", det),
                                 Form("Det %d;Z_{FOOT};Z_{LOS}", det),
                                 nbx, xmin, xmax, nby, ymin, ymax));
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
                    hists[det - 1]->Fill(fZCharge[i], zLos);
            }
        }
    }

    // save
    TFile out(repopath + "/results/cal/chargeId_perDet.root", "RECREATE");
    for (auto h : hists)
        h->Write();
    out.Close();

    fin->Close();
}