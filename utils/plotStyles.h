#pragma once

#include <TStyle.h>

inline void setOpenGL()
{
    gStyle->SetCanvasPreferGL(kTRUE);
}

inline void setHistogramStyle(TH1 *h, TString xTit, TString yTit, int color = kCyan, double xmin = 0, double xmax = 0)
{

    h->SetTitle("");
    h->SetFillColor(color - 3);
    h->SetFillStyle(1001);
    h->SetFillColorAlpha(color - 3, 0.3);
    h->SetLineColor(color);
    h->SetLineWidth(2);

    h->GetXaxis()
        ->SetTitleOffset(1.1);
    h->GetYaxis()->SetTitleOffset(1.8);

    h->GetXaxis()->SetTitle(xTit);
    h->GetYaxis()->SetTitle(yTit);

    h->GetXaxis()->SetAxisColor(kBlack);
    h->GetYaxis()->SetAxisColor(kBlack);

    if (xmin != 0 && xmax != 0)
        h->GetXaxis()->SetRangeUser(xmin, xmax);

    h->SetStats(0);
}

inline void setHistogramStyle(TH2 *h, TString xTit, TString yTit, double xmin = 0, double xmax = 0, double ymin = 0, double ymax = 0)
{
    gStyle->SetPalette(kViridis);
    h->SetTitle("");
    h->GetXaxis()
        ->SetTitleOffset(1.1);
    h->GetYaxis()->SetTitleOffset(1.2);

    h->GetXaxis()->SetTitle(xTit);
    h->GetYaxis()->SetTitle(yTit);

    h->GetXaxis()->SetAxisColor(kBlack);
    h->GetYaxis()->SetAxisColor(kBlack);

    if (xmin != 0 && xmax != 0)
        h->GetXaxis()->SetRangeUser(xmin, xmax);

    if (ymin != 0 && ymax != 0)
        h->GetYaxis()->SetRangeUser(ymin, ymax);

    h->SetStats(0);
}

inline void setCanvasStyle(TCanvas *c, double rightMargin = 0.05)
{
    gStyle->SetPadRightMargin(rightMargin);
    gStyle->SetPadLeftMargin(0.15);
    gStyle->SetCanvasPreferGL(kTRUE);
    gStyle->SetOptStat(0);
    c->SetLeftMargin(0.12);
    c->SetBottomMargin(0.12);
    c->SetTopMargin(0.08);
    c->SetFrameLineColor(kBlack);
    c->SetTicks(1, 1);

    c->SetFillStyle(4000);
    c->SetFillColor(0);
    c->SetFrameFillStyle(4000);
    c->SetFrameFillColor(0);
}

inline std::vector<Color_t> makeViridisColors(int n)
{
    std::vector<Color_t> colors;
    colors.reserve(n);

    gStyle->SetPalette(kViridis);

    for (int i = 0; i < n; ++i)
    {
        int idx = static_cast<int>(255.0 * i / (n - 1));
        colors.push_back(TColor::GetColorPalette(idx));
    }
    return colors;
}