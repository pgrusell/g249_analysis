// This class will:

// 1) Calculate the offsets from son TH2D to align the tof
// 2) Apply the offsets to show that it has worked

// * If the TH2D with ToF vs Paddle does not exist it can be retrived
//   from the corresponding rootfile with neuland hit data.

#pragma once

#include <fstream>

#include "../../utils/plotStyles.h"

class neulandAlignment
{
public:
    neulandAlignment(TString hitsFile = "", std::vector<TString> dataFile = {})
    {

        // Absolute path of the file
        TString histFilePath = static_cast<TString>(getenv("repopath")) + "/results/cal/" + hitsFile;

        // If the dataFile contains the histograms get it
        if (hitsFile != "")
        {
            auto *histDataFile = new TFile(histFilePath, "READ");
            auto *h = static_cast<TH2D *>(histDataFile->Get("hTofVsPad"));

            // Clean the histogram
            int threshold = 40;

            TH2D *hTofVsPaddlesNotCleaned = (TH2D *)h->RebinY(1, "hTofVsPaddlesRebbined");

            delete h;

            hTofVsPaddles = (TH2D *)hTofVsPaddlesNotCleaned->Clone("hTofVsPaddles");

            for (int ix = 1; ix <= hTofVsPaddles->GetNbinsX(); ix++)
            {
                for (int iy = 1; iy <= hTofVsPaddles->GetNbinsY(); iy++)
                {
                    double content = hTofVsPaddles->GetBinContent(ix, iy);
                    if (content < threshold)
                        hTofVsPaddles->SetBinContent(ix, iy, 0);
                }
            }
        }

        // Call the method to retrieve it from the data
        else
        {
            buildHistogram(dataFile, false);
        }
    }

    // This method will build the ToF vs Paddle histogram from the dataFile
    // If the offsets have already been calculated we can apply them
    void buildHistogram(std::vector<TString> dataFileAbs, bool withOffsets = false, bool both = false)
    {

        TString absPath = static_cast<TString>(getenv("repopath")) + "/data/";

        // Chain all files
        TChain chain("evt");
        for (const auto &f : dataFileAbs)
        {
            chain.Add(f);
        }

        // Reader
        TTreeReader reader(&chain);
        TTreeReaderArray<R3BNeulandHit> neuland(reader, "NeulandHits");

        int ievt = 0;
        const double Lref = 1557.0; // cm
        const double c = 29.979;    // cm/ns
        auto *hTofVsPad = new TH2D("hTofVsPad", "ToF vs Paddle;Paddle;ToF (ns)",
                                   1300, 0.0, 1300.0,
                                   200, 46.0, 55.0);

        while (reader.Next())
        {
            if ((++ievt % 100000) == 0)
            {
                std::cout << ievt << '\n';
            }

            for (auto const &neu : neuland)
            {
                // const double tof_corr = neu.GetT() - (neu.GetPosition().Mag() - Lref) / c - withOffsets * fTofOffsets[neu.GetPaddle() - 1];
                const double v_light = neu.GetPosition().Mag() / neu.GetT();
                const double tof_corr = Lref / v_light - withOffsets * fTofOffsets[neu.GetPaddle() - 1];
                hTofVsPad->Fill(neu.GetPaddle() - 1, tof_corr);
            }
        }

        withOffsets ? hTofVsPaddlesCorr = hTofVsPad : hTofVsPaddles = hTofVsPad;
        TString titHist = withOffsets ?

                                      auto *f = new TFile("histogram_built.root", "RECREATE");
        hTofVsPaddlesCorr->Write();
        delete f;

        if (both)
        {
        }
    }

    // This method takes the corrected histogram from rootFile
    void setCorrectedFromRoot(TString histoPath)
    {

        auto *histDataFile = new TFile(histoPath, "READ");
        auto *h = static_cast<TH2D *>(histDataFile->Get("hTofVsPad"));
        hTofVsPaddlesCorr = h;
    }

    // This method can be used to initialize the offsets
    void calculateOffsets()
    {

        TH1D *profile;

        // Iterate over all the bins and perform a Y projection
        for (int i = 0; i < fNPaddles; i++)
        {
            profile = hTofVsPaddles->ProjectionY(Form("h_proj_%i", i + 1), i + 1, i + 1);
            calculateOffsetFromProjection(profile);
            delete profile;
        }

        std::ofstream out("offsets_v1.txt");

        for (int i = 0; i < fNPaddles; i++)
        {
            out << i + 1 << " " << fTofOffsets[i] << " " << fFitPars[i][0] << " " << fFitPars[i][1] << std::endl;
        }

        out.close();
    }

    void setOffsetsFromTxt(std::string inputFile = "offsets_v1.txt")
    {
        std::ifstream input(inputFile);
        std::string line;
        fTofOffsets.resize(fNPaddles);

        while (std::getline(input, line))
        {
            std::istringstream ss(line);
            std::vector<double> fileVals;
            double val;

            while (ss >> val)
            {
                fileVals.push_back(val);
            }

            fTofOffsets[fileVals[0] - 1] = fileVals[1];
        }
    }

    // Method to fit the particular profile
    void calculateOffsetFromProjection(TH1D *profile)
    {
        if (profile->Integral() < 1)
        {
            fTofOffsets.push_back(0.);
            fFitPars.push_back(std::vector<double>{0., 0.});
            return;
        }

        double center = profile->GetBinCenter(profile->GetMaximumBin());
        double min = center - 0.5;
        double max = center + 0.5;

        auto f1 = new TF1("f_gauss", "gaus", 0, 0);
        f1->SetRange(min, max);

        profile->Fit(f1, "QR");
        double offset = f1->GetParameter(1) - 51.95;

        // profile->Draw();
        fTofOffsets.push_back(offset);
        fFitPars.push_back(std::vector<double>{f1->GetParameter(1), f1->GetParameter(2)});
    }

    // Method to check the alignment
    void checkAlignment(int nProfiles = 5, bool corrected = false)
    {

        // Plot different profiles to check the alignment
        TH2 *hToCheck = corrected ? (TH2 *)hTofVsPaddlesCorr : (TH2 *)hTofVsPaddles;
        TString figTit = corrected ? "figCorr.png" : "figUncorr.png";

        setOpenGL();
        auto c = new TCanvas("c");
        setCanvasStyle(c);

        std::vector<std::pair<TH1D *, double>> profiles;
        profiles.reserve(nProfiles);
        const int nY = fNPaddles;

        for (int i = 0; i < nProfiles; ++i)
        {
            const int num = 100 * i;
            TH1D *prof = hToCheck->ProjectionY(Form("h_prof_%d", i), num, num);
            profiles.emplace_back(prof, prof->GetMaximum());
        }

        std::sort(profiles.begin(), profiles.end(),
                  [](const auto &a, const auto &b)
                  { return a.second > b.second; });

        auto colors = makeViridisColors(profiles.size());

        for (auto i = 0; i < colors.size(); i++)
        {
            setHistogramStyle(profiles[i].first, "ToF [ns]", "# Counts", colors[i]);

            if (i == 0)
                profiles[i].first->Draw();
            else
                profiles[i].first->Draw("SAME");
        }

        c->SaveAs(figTit);
        // delete c;
    }

private:
    const int fNPaddles = 1300;
    TH2D *hTofVsPaddles = nullptr;
    TH2D *hTofVsPaddlesCorr = nullptr;
    std::vector<double> fTofOffsets;
    std::vector<std::vector<double>> fFitPars;
};
