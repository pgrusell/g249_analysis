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
    neulandAlignment(TString hitsFile = "", std::vector<TString> dataFile = {}, int TofOrDt = 1)
    {

        // Absolute path of the file
        TString histFilePath = static_cast<TString>(getenv("repopath")) + "/results/cal/" + hitsFile;

        fResults = new TFile(static_cast<TString>(getenv("repopath")) + "/results/cal/resultsTest.root", "UPDATE");

        // If the dataFile contains the histograms get it
        if (hitsFile != "")
        {

            // auto *histDataFile = new TFile(histFilePath, "READ");

            TString histName = TofOrDt ? "hUncorrectedToF" : "hUncorrecteddT";
            auto *h = static_cast<TH2D *>(fResults->Get(histName));

            // Clean the histogram
            int threshold = 0;
            TH2D *hTofVsPaddlesNotCleaned = (TH2D *)h->RebinY(1, "hTofVsPaddlesRebbined");
            neulandAlignmentAnalysis::cleanHistogram(hTofVsPaddlesNotCleaned, hTofVsPaddles, "hTofVsPaddles", threshold);

            hTofVsPaddles->SetDirectory(nullptr);
            // delete histDataFile;
        }
        // Call the method to retrieve it from the data
        else
        {
            buildHistogram(dataFile, false, TofOrDt);
        }

        hToFVsEnergy = new TH2D("hToFVsEnergy", "hToFVsEnergy:E[chn]:ToF[ns]", 1000, 0, 200, 1000, 40, 100);
        hXYPeak1 = new TH2D("hXYPeak1", "hXYPeak1;X [cm];Y [cm]", 1000, -250, 100, 1000, -250, 100);
        hXYPeak2 = new TH2D("hXYPeak2", "hXYPeak2;X [cm];Y [cm]", 1000, -250, 100, 1000, -250, 100);
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

    void buildHistogram(std::vector<TString> dataFileAbs, bool withOffsets = false, int TofOrDt = 1, int full = 0)
    {

        // Four options:
        // - dT = time_hit - time_target - L_hit / c ( + time_offset)
        // - ToF = [time_hit - time_target (+ time_offset)] / L_hit * L_0

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
        TTreeReaderArray<R3BLosHitData> los(reader, "LosHit");
        int ievt = 0;

        // If not offsets, set them to zero
        if (!withOffsets)
        {
            fTofOffsets.clear();
            fTofOffsets.resize(fNPaddles);
        }

        double binTofMin, binTofMax;
        int nBinsTof;
        TString histogramTitle;
        TString units;

        if (TofOrDt)
        {
            if (!full)
            {
                binTofMin = 0.;
                binTofMax = 500;
                nBinsTof = 11000;
                histogramTitle = "hToF";
                units = "ToF vs Paddle;Paddle;ToF (ns)";
            }

            else
            {
                binTofMin = 46.0;
                binTofMax = 55.0;
                nBinsTof = 200;
                histogramTitle = "hToF_full";
                units = "ToF vs Paddle;Paddle;ToF (ns)";
            }
        }
        else
        {
            binTofMin = -50.;
            binTofMax = 50.;
            nBinsTof = 2200;
            histogramTitle = "hdeltaTime";
            units = "#Deltat vs Paddle;Paddle;#Deltat (ns)";
        }

        if (withOffsets)
            histogramTitle += "Corrected";

        histogramTitle += "VsPaddle";

        auto *hTofVsPad = new TH2D(histogramTitle, units,
                                   1300, -0.5, 1299.5,
                                   nBinsTof, binTofMin, binTofMax);

        // Calculate dT or ToF
        const double Lref = 1557.0; // cm
        const double c = 29.979;    // cm/ns
        double timeTarget = 7.423887;

        if (!withOffsets)
            timeTarget = 0.;

        auto computeTime = [=](double timeHit, double distHit, double timeOffset)
        {
            const double timeCorr = timeHit - timeTarget - timeOffset;

            if (TofOrDt)
                return timeCorr / distHit * Lref;
            else
                return timeCorr - distHit / c;
        };

        // Iterate over events
        while (reader.Next())
        {
            if ((++ievt % 100000) == 0)
            {
                std::cout << 100.0 * ievt / chain.GetEntries() << " %\n";
            }

            // Only 1 hit in LOS
            if (los.GetSize() != 1)
                continue;

            /*
            for (auto const &neu : neuland)
            {

                const double time = computeTime(neu.GetT(), neu.GetPosition().Mag(), fTofOffsets[neu.GetPaddle() - 1]);

                if (!std::isnan(time) && neu.GetT() > 0)
                {

                    const double energy = neu.GetE();
                    double xPos = neu.GetPosition().x();
                    double yPos = neu.GetPosition().y();

                    // std::cout << energy << " " << xPos << " " << yPos << std::endl;

                    // Identify the coordinate measured by the bar
                    if (2 * xPos == static_cast<int>(2 * xPos))
                        xPos = gRandom->Uniform(xPos - 2.5, xPos + 2.5);

                    if (2 * yPos == static_cast<int>(2 * yPos))
                        yPos = gRandom->Uniform(yPos - 2.5, yPos + 2.5);

                    hTofVsPad->Fill(neu.GetPaddle() - 1, time);
                    hToFVsEnergy->Fill(energy, time);

                    if ((time > 51.5) && (time < 52.4))
                        hXYPeak1->Fill(xPos, yPos);

                    if ((time > 52.4) && (time < 52.9))
                        hXYPeak2->Fill(xPos, yPos);
                }
            }
            */

            if (neuland.GetSize() == 0)
                continue;

            auto neu = neuland[0];
            const double time = computeTime(neu.GetT(), neu.GetPosition().Mag(), fTofOffsets[neu.GetPaddle() - 1]);

            if (!std::isnan(time) && neu.GetT() > 0)
            {

                const double energy = neu.GetE();
                double xPos = neu.GetPosition().x();
                double yPos = neu.GetPosition().y();

                // std::cout << energy << " " << xPos << " " << yPos << std::endl;

                // Identify the coordinate measured by the bar
                if (2 * xPos == static_cast<int>(2 * xPos))
                    xPos = gRandom->Uniform(xPos - 2.5, xPos + 2.5);

                if (2 * yPos == static_cast<int>(2 * yPos))
                    yPos = gRandom->Uniform(yPos - 2.5, yPos + 2.5);

                if (energy < 6.)
                    continue;

                hTofVsPad->Fill(neu.GetPaddle() - 1, time);
                hToFVsEnergy->Fill(energy, time);

                if ((time > 51.5) && (time < 52.4))
                    hXYPeak1->Fill(xPos, yPos);

                if ((time > 52.4) && (time < 52.9))
                    hXYPeak2->Fill(xPos, yPos);
            }
        }

        withOffsets ? hTofVsPaddlesCorr = hTofVsPad : hTofVsPaddles = hTofVsPad;
        TString titHist = withOffsets ? "hCorrected" : "hUncorrected";

        if (TofOrDt)
        {
            if (full)
                titHist += "_full";

            titHist += "ToF";
        }
        else
            titHist += "dT";

        auto *hToSave = hTofVsPad->Clone(titHist);
        fResults->cd();
        hToSave->Write("", TObject::kOverwrite);
        hTofVsPad->Write();
        hToFVsEnergy->Write();
        hXYPeak1->Write();
        hXYPeak2->Write();
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

        fTofOffsets.clear();

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
    TH2D *hToFVsEnergy = nullptr;
    TH2D *hXYPeak1 = nullptr;
    TH2D *hXYPeak2 = nullptr;
    std::vector<double> fTofOffsets;
    std::vector<std::vector<double>> fFitPars;
    TFile *fResults = nullptr;
};
