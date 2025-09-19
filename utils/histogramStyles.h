#pragma once

inline void setHistogramStyle(TH1F *h, TString xTit, TString yTit)
{
    h->SetFillColor(kCyan - 3);
    h->SetFillColorAlpha(kCyan - 3, 0.3);
    h->SetLineColor(kCyan + 2);
    h->SetLineWidth(2);

    h->GetXaxis()
        ->SetTitleOffset(1.1);
    h->GetXaxis()->SetTitleOffset(1.2);

    h->GetXaxis()->SetTitle(xTit);
    h->GetYaxis()->SetTitle(yTit);

    h->GetXaxis()->SetAxisColor(kBlack);
    h->GetYaxis()->SetAxisColor(kBlack);
}