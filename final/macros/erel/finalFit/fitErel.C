#include "Math/Factory.h"
#include "Math/Functor.h"
#include "Math/Minimizer.h"

double Penetrability(double E, int l)
{
    const double R = 5.4;    // 1.4*(23^(1/3)+1) fm
    const double mu = 900.4; // 23/24 * m_n, MeV

    if (E <= 0.)
        return 0.;

    double rho = std::sqrt(2. * mu * E) / 197.327 * R;
    double r2 = rho * rho;

    switch (l)
    {
    case 0:
        return rho;
    case 1:
        return rho * r2 / (1. + r2);
    case 2:
        return rho * r2 * r2 / (9. + 3. * r2 + r2 * r2);
    default:
        return rho;
    }
}

double BW(double E, double Er, double Gr, int l)
{
    if (E <= 0.)
        return 0.;

    double G = Gr * Penetrability(E, l) / Penetrability(Er, l);
    return (G / (2. * TMath::Pi())) / ((E - Er) * (E - Er) + 0.25 * G * G);
}

struct FoldedNLLs
{
    std::vector<std::vector<double>> R; // R[j][i] = P(rec j | gen i)
    std::vector<double> Egen, dEgen;
    std::vector<double> data, bg;
    std::vector<int> lvals = {2, 0, 2, 1}; // 2+, 1+, res(3), res(6)
    int nGen = 0, nRec = 0;
    int nRes = 4;
    int jMax = -1; // ultimo bin de rec usado en el fit (-1 = todos)

    double Model(double E, const double *p) const
    {
        double v = 0.;
        for (int k = 0; k < nRes; k++)
            v += p[3 * k] * BW(E, p[3 * k + 1], p[3 * k + 2], lvals[k]);
        return v;
    }

    double Predict(int j, const double *p) const
    {
        double nu = p[3 * nRes] * bg[j];
        for (int i = 0; i < nGen; i++)
            nu += R[j][i] * Model(Egen[i], p) * dEgen[i];
        return nu;
    }

    double operator()(const double *p) const
    {
        int jEnd = (jMax > 0) ? jMax : nRec;
        double nll = 0.;
        for (int j = 0; j < jEnd; j++)
        {
            double nu = Predict(j, p);
            if (nu <= 0.)
                nu = 1e-12;
            nll += nu - data[j] * std::log(nu);
        }
        return 2. * nll;
    }
};

void fitErel()
{
    int nBins = 80;
    double minRec = 0, maxRec = 14;
    int nGenBins = 500;
    double minGen = 0, maxGen = 25;

    TString cut = "califa_opa > 1.25 && califa_opa < 1.65";

    //////// Datos
    auto *f = new TFile("/nucl_lustre/pablogrusell/g249/g249_analysis/results/final/23O_analyzed_final.root");
    auto *tr = static_cast<TTree *>(f->Get("KinTree"));

    auto *hErelRec = new TH1F("hErelRec", "", nBins, minRec, maxRec);
    tr->Draw("Erel*1000>>hErelRec", cut, "goff");

    //////// Mirar la region baja con binado fino ANTES de ajustar nada
    auto *hLow = new TH1F("hLow", ";E_{rel} [MeV];counts / 50 keV", 40, 0, 2);
    tr->Draw("Erel*1000>>hLow", cut, "goff");

    auto *cLow = new TCanvas("cLow", "region baja", 700, 500);
    hLow->SetMarkerStyle(20);
    hLow->Draw("E");

    //////// Matriz de respuesta
    auto *fSim = new TFile("/nucl_lustre/pablogrusell/g249/g249_analysis/final/macros/erel/ana_results/full2_analysis.root");
    auto *trSim = static_cast<TTree *>(fSim->Get("ana"));

    auto *hResp = new TH2F("hResp", "", nGenBins, minGen, maxGen, nBins, minRec, maxRec);
    hResp->Sumw2();
    trSim->Draw("Erel*1000:ErelTrue*1000>>hResp", "", "goff");

    auto *hErelGen = new TH1F("hErelGen", "", nGenBins, minGen, maxGen);
    trSim->Draw("ErelTrue*1000>>hErelGen", "", "goff");

    for (int i = 1; i <= hResp->GetNbinsX(); i++)
    {
        double norm = hErelGen->GetBinContent(i);
        if (norm <= 0)
            continue;

        for (int j = 1; j <= hResp->GetNbinsY(); j++)
        {
            hResp->SetBinContent(i, j, hResp->GetBinContent(i, j) / norm);
            hResp->SetBinError(i, j, hResp->GetBinError(i, j) / norm);
        }
    }

    //////// Fondo no resonante
    auto *fBg = new TFile("/nucl_lustre/pablogrusell/g249/g249_analysis/final/macros/erel/23O_mixing_iterative.root");
    auto *trBg = static_cast<TTree *>(fBg->Get("tMix"));

    auto *hBg = new TH1F("hBg", "", nBins, minRec, maxRec);
    trBg->Draw("Erel>>hBg", "weight", "goff");
    hBg->Scale(1. / hBg->Integral());

    //////// Cachear
    FoldedNLL fcn;
    fcn.nGen = hResp->GetNbinsX();
    fcn.nRec = hResp->GetNbinsY();

    fcn.R.assign(fcn.nRec, std::vector<double>(fcn.nGen, 0.));
    for (int j = 1; j <= fcn.nRec; j++)
        for (int i = 1; i <= fcn.nGen; i++)
            fcn.R[j - 1][i - 1] = hResp->GetBinContent(i, j);

    for (int i = 1; i <= fcn.nGen; i++)
    {
        fcn.Egen.push_back(hResp->GetXaxis()->GetBinCenter(i));
        fcn.dEgen.push_back(hResp->GetXaxis()->GetBinWidth(i));
    }

    for (int j = 1; j <= fcn.nRec; j++)
    {
        fcn.data.push_back(hErelRec->GetBinContent(j));
        fcn.bg.push_back(hBg->GetBinContent(j));
    }

    // fcn.jMax = hErelRec->FindBin(8.0);   // descomentar para restringir a 0-8

    //////// Comprobacion: que anchura permite la penetrabilidad
    printf("\n--- penetrabilidad ---\n");
    for (double E : {0.3, 0.6, 1.0, 1.3, 3.1, 6.2})
        printf("E = %.2f : P(l=0) = %.4f   P(l=1) = %.4f   P(l=2) = %.5f\n",
               E, Penetrability(E, 0), Penetrability(E, 1), Penetrability(E, 2));
    printf("\n");

    //////// Fit
    int nPar = 3 * fcn.nRes + 1;

    auto *min = ROOT::Math::Factory::CreateMinimizer("Minuit2", "Migrad");
    min->SetErrorDef(1.);
    min->SetPrintLevel(1);
    min->SetStrategy(2);

    ROOT::Math::Functor fun(fcn, nPar);
    min->SetFunction(fun);

    // 2+ : l=2, obligado por reglas de seleccion. Estrecho.
    min->SetLimitedVariable(0, "A_2p", 300., 10., 0., 1e6);
    min->SetLimitedVariable(1, "E_2p", 0.65, 0.005, 0.25, 1.05);
    min->SetLimitedVariable(2, "G_2p", 0.05, 0.002, 0.001, 0.40);

    // 1+ : l=0 domina (sin barrera)
    min->SetLimitedVariable(3, "A_1p", 300., 10., 0., 1e6);
    min->SetLimitedVariable(4, "E_1p", 1.25, 0.005, 0.85, 1.90);
    min->SetLimitedVariable(5, "G_1p", 0.10, 0.005, 0.001, 0.70);

    // resonancia ~3 MeV
    min->SetLimitedVariable(6, "A_r3", 1300., 10., 0., 1e6);
    min->SetLimitedVariable(7, "E_r3", 3.10, 0.01, 2.2, 4.2);
    min->SetLimitedVariable(8, "G_r3", 0.80, 0.01, 0.02, 3.0);

    // resonancia ~6 MeV
    min->SetLimitedVariable(9, "A_r6", 2500., 10., 0., 1e6);
    min->SetLimitedVariable(10, "E_r6", 6.20, 0.01, 5.0, 7.5);
    min->SetLimitedVariable(11, "G_r6", 1.20, 0.01, 0.02, 5.0);

    min->SetLimitedVariable(12, "Abg", 800., 10., 0., 1e6);

    min->Minimize();
    min->Hesse();

    const double *par = min->X();
    const double *err = min->Errors();

    printf("\n--- resultados ---\n");
    const char *names[4] = {"2+ ", "1+ ", "r3 ", "r6 "};
    for (int k = 0; k < fcn.nRes; k++)
        printf("%s (l=%d): Er = %.3f +- %.3f   G = %.4f +- %.4f   A = %.0f +- %.0f\n",
               names[k], fcn.lvals[k],
               par[3 * k + 1], err[3 * k + 1],
               par[3 * k + 2], err[3 * k + 2],
               par[3 * k], err[3 * k]);
    printf("Abg = %.1f +- %.1f  (%.1f%% de %.0f eventos)\n",
           par[12], err[12], 100. * par[12] / hErelRec->Integral(), hErelRec->Integral());

    printf("\n--- correlaciones con Abg ---\n");
    for (int k = 0; k < fcn.nRes; k++)
        printf("corr(A_%s, Abg) = %+.3f\n", names[k], min->Correlation(3 * k, 12));

    //////// Dibujar
    auto *c1 = new TCanvas("c1", "fit", 800, 600);

    auto *hFit = static_cast<TH1F *>(hErelRec->Clone("hFit"));
    hFit->Reset();
    for (int j = 0; j < fcn.nRec; j++)
        hFit->SetBinContent(j + 1, fcn.Predict(j, par));

    auto *hBgFit = static_cast<TH1F *>(hBg->Clone("hBgFit"));
    hBgFit->Scale(par[12]);

    hErelRec->SetMarkerStyle(20);
    hErelRec->Draw("E");
    hFit->SetLineColor(kRed);
    hFit->SetLineWidth(2);
    hFit->Draw("hist same");
    hBgFit->SetLineColor(kGreen + 2);
    hBgFit->SetLineWidth(2);
    hBgFit->Draw("hist same");

    int cols[4] = {kMagenta, kOrange + 7, kBlue, kCyan + 2};
    for (int k = 0; k < fcn.nRes; k++)
    {
        std::vector<double> pk(nPar, 0.);
        pk[3 * k] = par[3 * k];
        pk[3 * k + 1] = par[3 * k + 1];
        pk[3 * k + 2] = par[3 * k + 2];

        auto *hk = static_cast<TH1F *>(hErelRec->Clone(Form("hRes%d", k)));
        hk->Reset();
        for (int j = 0; j < fcn.nRec; j++)
            hk->SetBinContent(j + 1, fcn.Predict(j, pk.data()));

        hk->SetLineColor(cols[k]);
        hk->SetLineStyle(2);
        hk->Draw("hist same");
    }
}