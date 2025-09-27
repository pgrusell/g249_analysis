// This class will:

// 1) Calculate the offsets from son TH2D to align the tof
// 2) Apply the offsets to show that it has worked

// * If the TH2D with ToF vs Paddle does not exist it can be retrived
//   from the corresponding rootfile with neuland hit data.

#pragma once

#include <fstream>

#include "../../utils/plotStyles.h"
#include "neuland_alignment_analysis.hh"

class neulandAlignment
{
public:
    neulandAlignment(TString hitsFile = "", std::vector<TString> dataFile = {})
    {

        // Absolute path of the file
        TString histFilePath = static_cast<TString>(getenv("repopath")) + "/results/cal/" + hitsFile;

        fResults = new TFile("resultsTest.root", "RECREATE");

        // If the dataFile contains the histograms get it
        if (hitsFile != "")
        {

            auto *histDataFile = new TFile(histFilePath, "READ");
            auto *h = static_cast<TH2D *>(histDataFile->Get("hTofVsPad"));

            // Clean the histogram
            int threshold = 0;
            TH2D *hTofVsPaddlesNotCleaned = (TH2D *)h->RebinY(1, "hTofVsPaddlesRebbined");
            neulandAlignmentAnalysis::cleanHistogram(hTofVsPaddlesNotCleaned, hTofVsPaddles, "hTofVsPaddles", threshold);

            hTofVsPaddles->SetDirectory(nullptr);
            delete histDataFile;
        }
        // Call the method to retrieve it from the data
        else
        {
            buildHistogram(dataFile, false);
        }
    }

    ~neulandAlignment()
    {
        if (fResults)
        {
            fResults->Write();
            fResults->Close();
            delete fResults;
            fResults = nullptr;
        }
    }

    void buildHistogram(std::vector<TString> dataFileAbs, bool withOffsets = false, int TofOrSpeed = 0)
    {

        if (withOffsets && fTofOffsets.size() == 0)
        {
            std::cerr << "[FATAL] Offsets have not been initialized\n";
            return;
        }

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
        const double timeTarget = 7.423887;

        double binTofMin = -50.;
        double binTofMax = 50.;
        double nBinsTof = 2200;

        if (withOffsets)
        {
            binTofMin = 46.0;
            binTofMax = 55.0;
            nBinsTof = 200;
        }

        /*
        auto *hTofVsPad = new TH2D("hTofVsPad", "ToF vs Paddle;Paddle;ToF (ns)",
                                   1300, -0.5, 1299.5,
                                   200, 46.0, 55.0);
        */

        auto *hTofVsPad = new TH2D("hTofVsPad", "ToF vs Paddle;Paddle;ToF (ns)",
                                   1300, -0.5, 1299.5,
                                   nBinsTof, binTofMin, binTofMax);

        while (reader.Next())
        {
            if ((++ievt % 100000) == 0)
            {
                std::cout << 100.0 * ievt / chain.GetEntries() << " %\n";
            }

            for (auto const &neu : neuland)
            {

                double dToFOffset = 0;

                if (withOffsets)
                {
                    dToFOffset = (std::abs(fTofOffsets[neu.GetPaddle() - 1]) < 2.5) ? fTofOffsets[neu.GetPaddle() - 1] : 0;
                }

                const double timeDiff = neu.GetT() - timeTarget - neu.GetPosition().Mag() / c;

                if (!std::isnan(timeDiff) && neu.GetT() > 0)
                    hTofVsPad->Fill(neu.GetPaddle() - 1, timeDiff - dToFOffset);
            }
        }

        withOffsets ? hTofVsPaddlesCorr = hTofVsPad : hTofVsPaddles = hTofVsPad;
        TString titHist = withOffsets ? "hCorrected" : "hUncorrected";

        auto *hToSave = hTofVsPad->Clone(titHist);
        fResults->cd();
        hToSave->Write();
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

        std::ofstream out("offsets_v2.txt");

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

        int binMin = profile->FindBin(-9);
        int binMax = profile->FindBin(-7.6);

        int maxBin = binMin;
        double maxContent = profile->GetBinContent(binMin);

        for (int b = binMin + 1; b <= binMax; ++b)
        {
            double content = profile->GetBinContent(b);
            if (content > maxContent)
            {
                maxContent = content;
                maxBin = b;
            }
        }

        double center = profile->GetBinCenter(maxBin);

        // double center = profile->GetBinCenter(profile->GetMaximumBin());
        double min = center - 0.5;
        double max = center + 0.5;

        auto f1 = new TF1("f_gauss", "gaus", 0, 0);
        f1->SetRange(min, max);

        profile->Fit(f1, "QR");
        double offset = f1->GetParameter(1);

        // profile->Draw();
        fTofOffsets.push_back(offset);
        fFitPars.push_back(std::vector<double>{f1->GetParameter(1), f1->GetParameter(2)});
    }

private:
    const int fNPaddles = 1300;
    TH2D *hTofVsPaddles = nullptr;
    TH2D *hTofVsPaddlesCorr = nullptr;
    TH2D *hOffsets = nullptr;
    std::vector<double> fTofOffsets;
    std::vector<std::vector<double>> fFitPars;
    TFile *fResults = nullptr;
};
