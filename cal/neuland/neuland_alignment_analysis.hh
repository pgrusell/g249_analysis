#pragma once

#include "../../utils/plotStyles.h"

class neulandAlignmentAnalysis
{
public:
    neulandAlignmentAnalysis(TString resultsFile = "")
    {

        TString resultsPath = static_cast<TString>(getenv("repopath")) + "/results/cal/";

        TString histFilePath = resultsPath + resultsFile;
        std::cout << histFilePath << std::endl;
        auto *histDataFile = new TFile(histFilePath, "READ");

        hUncorrecteddT = static_cast<TH2D *>(histDataFile->Get("hUncorrecteddT"));
        hCorrecteddT = static_cast<TH2D *>(histDataFile->Get("hCorrecteddT"));

        hUncorrectedToF = static_cast<TH2D *>(histDataFile->Get("hUncorrectedToF"));
        hCorrectedToF = static_cast<TH2D *>(histDataFile->Get("hCorrectedToF"));

        cleanHistogram(hUncorrecteddT, hUncorrectedCleaned, "hUncorrectedCleaned");
        cleanHistogram(hCorrecteddT, hCorrectedCleaned, "hCorrectedCleaned");

        // save ToF vs paddle histograms in pdf files
        printCanvases(resultsPath);

        // save the ToF projection over Y axis
        totalProfile(resultsPath);

        // plot profiles for uncorrected
        checkAlignment(hUncorrecteddT, resultsPath + "single_prof_uncorr_dt", 8, -10, -6, "#Deltat [ns]", true);

        // plot profiles for corrected
        checkAlignment(hCorrecteddT, resultsPath + "single_prof_corr_dt.pdf", 8, -2, -12, "#Deltat [ns]", false);

        // plot profiles for uncorrected
        checkAlignment(hUncorrectedToF, resultsPath + "single_prof_uncorr_tof", 10, 46, 54, "ToF [ns]");

        // plot profiles for corrected
        checkAlignment(hCorrectedToF, resultsPath + "single_prof_corr_tof", 10, 46, 54, "ToF [ns]");
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

    static void printCanvas(TH2D *&h, TString title, TString xTit, TString yTit, TString path, double ymin = 0, double ymax = 0)
    {

        auto *c = new TCanvas("c", title, 800, 600);
        setCanvasStyle(c, 0.15);
        setHistogramStyle(h, xTit, yTit, 0, 0, ymin, ymax);
        h->Draw("COLZ");
        c->SaveAs(path + title + ".pdf");
        delete c;
    }

    void printCanvases(TString path)
    {
        setOpenGL();

        printCanvas(hUncorrecteddT, "hUncorrecteddT", "Paddle #", "#Deltat [ns]", path, -10, -6);
        printCanvas(hCorrecteddT, "hCorrecteddT", "Paddle #", "#Deltat [ns]", path, -2, 2);
        printCanvas(hUncorrectedCleaned, "hUncorrectedCleaned", "Paddle #", "#Deltat [ns]", path, -10, -6);
        printCanvas(hCorrectedCleaned, "hCorrectedCleaned", "Paddle #", "#Deltat [ns]", path, -2, 2);
        printCanvas(hUncorrectedToF, "hUncorrectedToF", "Paddle #", "ToF [ns]", path);
        printCanvas(hCorrectedToF, "hCorrectedToF", "Paddle #", "ToF [ns]", path);
    }

    void totalProfile(TString path)
    {
        auto *prCorrected = hCorrectedToF->ProjectionY("prCorr");
        auto *prUncorrected = hUncorrectedToF->ProjectionY("prUncorr");

        setHistogramStyle(prCorrected, "ToF [ns]", "# counts", kCyan - 6);
        setHistogramStyle(prUncorrected, "ToF [ns]", "# counts", kOrange - 3);

        auto *c5 = new TCanvas("c5", "hProfTotBoth", 800, 600);
        setCanvasStyle(c5);

        prCorrected->Draw();
        prUncorrected->Draw("same");

        auto *leg = new TLegend(0.15, 0.75, 0.25, 0.9);
        leg->SetTextSize(0.04);
        leg->SetBorderSize(0);
        leg->SetFillStyle(0);
        leg->AddEntry(prCorrected, "new params (FWHM = 0.6 ns)", "l");
        leg->AddEntry(prUncorrected, "old params (FWHM = 0.8 ns)", "l");
        leg->Draw();

        c5->SaveAs(path + "hProfTot.pdf");

        delete c5;
    }

    void checkAlignment(TH2D *&hInit, TString title, int nProfiles = 5, double xmin = -10, double xmax = -6, TString lab = "#Deltat [ns]", bool single = false)
    {

        TString titleOne = title + "One.pdf";
        title += ".pdf";

        TH2D *h = (TH2D *)hInit->RebinY(2, "h");

        setOpenGL();
        auto c6 = new TCanvas("c6", "", 800, 600);
        setCanvasStyle(c6);

        std::vector<std::pair<TH1D *, double>> profiles;
        profiles.reserve(nProfiles);
        const int nY = fNPaddles;

        for (int i = 0; i < nProfiles; ++i)
        {
            const int num = 100 * i + 10;
            TH1D *prof = h->ProjectionY(Form("h_prof_%d", i), num, num);
            int binLow = prof->FindBin(xmin);
            int binHigh = prof->FindBin(xmax);

            double maxInRange = 0.0;
            for (int b = binLow; b <= binHigh; ++b)
            {
                double content = prof->GetBinContent(b);
                if (content > maxInRange)
                    maxInRange = content;
            }

            profiles.emplace_back(prof, maxInRange);
        }

        std::sort(profiles.begin(), profiles.end(),
                  [](const auto &a, const auto &b)
                  { return a.second > b.second; });

        auto colors = makeViridisColors(profiles.size());
        auto legend = new TLegend(0.8 - 0.03, 0.7 - 0.03, 1. - 0.03, 1. - 0.03);

        for (auto i = 0; i < colors.size(); i++)
        {
            int paddleNum = 100 * i + 10;
            TString entry = Form("Paddle %d", paddleNum);
            setHistogramStyle(profiles[i].first, lab, "# Counts", colors[i], xmin, xmax);

            if (i == 0)
                profiles[i].first->Draw();
            else
                profiles[i].first->Draw("SAME");

            legend->AddEntry(profiles[i].first, entry, "l");
        }

        legend->SetLineWidth(2);

        legend->Draw();
        c6->SaveAs(title);

        delete c6;
        delete h;
        delete legend;

        if (single)
        {
            // Plot a fit over one of the profiles
            auto c7 = new TCanvas("c7", "", 800, 600);
            setCanvasStyle(c7);
            int k = 2;
            profiles[k].first->Draw();

            int binMin = profiles[k].first->FindBin(xmin);
            int binMax = profiles[k].first->FindBin(xmax);

            int maxBin = binMin;
            double maxContent = profiles[k].first->GetBinContent(binMin);

            for (int b = binMin + 1; b <= binMax; ++b)
            {
                double content = profiles[k].first->GetBinContent(b);
                if (content > maxContent)
                {
                    maxContent = content;
                    maxBin = b;
                }
            }

            double center = profiles[k].first->GetBinCenter(maxBin);

            // double center = profile->GetBinCenter(profile->GetMaximumBin());
            double min = center - 0.3;
            double max = center + 0.3;

            auto f1 = new TF1("f_gauss", "gaus", 0, 0);
            f1->SetRange(min, max);

            profiles[k].first->Fit(f1, "R");

            auto leg = new TLegend(0.8, 0.85, 0.97, 0.95);

            TString muStr = Form("#mu = %.2f ns", f1->GetParameter(1));
            TString sigmaStr = Form("#sigma = %.2f ns", f1->GetParameter(2));

            leg->AddEntry((TObject *)0, muStr, "");    // (TObject*)0 = texto sin símbolo
            leg->AddEntry((TObject *)0, sigmaStr, ""); // segunda línea
            leg->Draw();

            c7->SaveAs(titleOne);

            delete c7;
        }
    }

private:
    const int fNPaddles = 1300;

    TH2D *hUncorrecteddT = nullptr;
    TH2D *hCorrecteddT = nullptr;

    TH2D *hUncorrectedToF = nullptr;
    TH2D *hCorrectedToF = nullptr;

    TH2D *hUncorrectedCleaned = nullptr;
    TH2D *hCorrectedCleaned = nullptr;
};