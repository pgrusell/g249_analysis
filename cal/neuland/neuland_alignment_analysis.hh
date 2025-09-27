#pragma once

#include "../../utils/plotStyles.h"

class neulandAlignmentAnalysis
{
public:
    neulandAlignmentAnalysis(TString resultsFile = "")
    {

        TString resultsPath = static_cast<TString>(getenv("repopath")) + "/results/cal/";

        TString histFilePath = resultsPath + "/results/cal/" + resultsFile;
        auto *histDataFile = new TFile(resultsFile, "READ");

        hUncorrected = static_cast<TH2D *>(histDataFile->Get("hUncorrected"));
        hCorrected = static_cast<TH2D *>(histDataFile->Get("hCorrected"));

        cleanHistogram(hUncorrected, hUncorrectedCleaned, "hUncorrectedCleaned");
        cleanHistogram(hCorrected, hCorrectedCleaned, "hCorrectedCleaned");

        // save ToF vs paddle histograms in pdf files
        printCanvases();

        // save the ToF projection over Y axis
        totalProfile();

        // plot profiles for uncorrected
        checkAlignment(hUncorrected, resultsPath + "single_prof_uncorr.pdf", 10);

        // plot profiles for corrected
        checkAlignment(hCorrected, resultsPath + "single_prof_corr.pdf", 10);
    }

    static void cleanHistogram(TH2D *&h, TH2D *&hClean, TString name, int threshold = 30)
    {
        hClean = (TH2D *)h->Clone(name);

        for (int ix = 1; ix <= hClean->GetNbinsX(); ix++)
        {
            for (int iy = 1; iy <= hClean->GetNbinsY(); iy++)
            {
                double content = hClean->GetBinContent(ix, iy);
                if (content < threshold)
                    hClean->SetBinContent(ix, iy, 0);
            }
        }
    }

    void printCanvases()
    {
        setOpenGL();

        // hUncorrected
        auto *c1 = new TCanvas("c1", "hUncorrected", 800, 600);
        setCanvasStyle(c1);
        setHistogramStyle(hUncorrected, "Paddle #", "ToF [ns]");
        hUncorrected->Draw("COLZ");
        c1->SaveAs("hUncorrected.pdf");
        delete c1;

        // hCorrected
        auto *c2 = new TCanvas("c2", "hCorrected", 800, 600);
        setCanvasStyle(c2);
        setHistogramStyle(hCorrected, "Paddle #", "ToF [ns]");
        hCorrected->Draw("COLZ");
        c2->SaveAs("hCorrected.pdf");
        delete c2;

        // hUncorrectedCleaned
        auto *c3 = new TCanvas("c3", "hUncorrectedCleaned", 800, 600);
        setCanvasStyle(c3);
        setHistogramStyle(hUncorrectedCleaned, "Paddle #", "ToF [ns]");
        hUncorrectedCleaned->Draw("COLZ");
        c3->SaveAs("hUncorrectedCleaned.pdf");
        delete c3;

        // hCorrectedCleaned
        auto *c4 = new TCanvas("c4", "hCorrectedCleaned", 800, 600);
        setCanvasStyle(c4);
        setHistogramStyle(hCorrectedCleaned, "Paddle #", "ToF [ns]");
        hCorrectedCleaned->Draw("COLZ");
        c4->SaveAs("hCorrectedCleaned.pdf");
        delete c4;
    }

    void totalProfile()
    {

        auto *prCorrected = hCorrected->ProjectionY("prCorr");
        auto *prUncorrected = hUncorrected->ProjectionY("prUncorr");

        setHistogramStyle(prCorrected, "ToF [ns]", "# counts", kCyan - 6);
        setHistogramStyle(prUncorrected, "ToF [ns]", "# counts", kOrange - 3);

        auto *c5 = new TCanvas("c5", "hProfTotBoth", 800, 600);
        setCanvasStyle(c5);
        prCorrected->Draw();
        prUncorrected->Draw("same");
        c5->SaveAs("hProfTot.pdf");
        delete c5;
    }

    void checkAlignment(TH2D *&hInit, TString title, int nProfiles = 5)
    {

        TH2D *h = (TH2D *)hInit->RebinY(5, "h");

        setOpenGL();
        auto c6 = new TCanvas("c6");
        setCanvasStyle(c6);

        std::vector<std::pair<TH1D *, double>> profiles;
        profiles.reserve(nProfiles);
        const int nY = fNPaddles;

        for (int i = 0; i < nProfiles; ++i)
        {
            const int num = 100 * i + 10;
            TH1D *prof = h->ProjectionY(Form("h_prof_%d", i), num, num);
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

        c6->SaveAs(title);
        delete c6;
        delete h;
    }

private:
    const int fNPaddles = 1300;

    TH2D *hUncorrected = nullptr;
    TH2D *hCorrected = nullptr;

    TH2D *hUncorrectedCleaned = nullptr;
    TH2D *hCorrectedCleaned = nullptr;
};