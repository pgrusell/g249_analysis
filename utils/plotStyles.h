#pragma once

#include <TStyle.h>

inline void setOpenGL()
{
    gStyle->SetCanvasPreferGL(kTRUE);
}

inline void setHistogramStyle(TH1 *h, TString xTit, TString yTit, int color = kCyan)
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
}

inline void setHistogramStyle(TH2 *h, TString xTit, TString yTit)
{
    h->SetTitle("");
    h->GetXaxis()
        ->SetTitleOffset(1.1);
    h->GetYaxis()->SetTitleOffset(1.2);

    h->GetXaxis()->SetTitle(xTit);
    h->GetYaxis()->SetTitle(yTit);

    h->GetXaxis()->SetAxisColor(kBlack);
    h->GetYaxis()->SetAxisColor(kBlack);
}

inline void setCanvasStyle(TCanvas *c)
{
    gStyle->SetPadRightMargin(0.05);
    gStyle->SetPadLeftMargin(0.15);
    // gStyle->SetPadRightMargin(0.15);
    gStyle->SetCanvasPreferGL(kTRUE);
    gStyle->SetOptStat(0);
    c->SetLeftMargin(0.12);
    c->SetBottomMargin(0.12);
    c->SetTopMargin(0.08);
    c->SetFrameLineColor(kBlack);
    c->SetTicks(1, 1);
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