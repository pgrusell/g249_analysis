#pragma once

#include <TStyle.h>

inline void setHistogramStyle(TH1 *h, TString xTit, TString yTit, int color = kCyan)
{

    h->SetFillColor(color - 3);
    h->SetFillColorAlpha(color - 3, 0.3);
    h->SetLineColor(kBlack);
    h->SetLineWidth(2);

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
    gStyle->SetCanvasPreferGL(kTRUE);
    c->SetLeftMargin(0.12);
    c->SetBottomMargin(0.12);
    c->SetTopMargin(0.08);
    c->SetRightMargin(0.05);
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