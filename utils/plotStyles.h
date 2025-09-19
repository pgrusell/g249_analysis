#pragma once

inline void setHistogramStyle(TH1F *h, TString xTit, TString yTit)
{
    h->SetFillColor(kCyan - 3);
    h->SetFillColorAlpha(kCyan - 3, 0.3);
    h->SetLineColor(kCyan + 2);
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
    c->SetLeftMargin(0.12);
    c->SetBottomMargin(0.12);
    c->SetTopMargin(0.08);
    c->SetRightMargin(0.05);
    c->SetFrameLineColor(kBlack);
    c->SetTicks(1, 1);
}