// eventMixingIterative.C
//   Fondo no resonante por EVENT-MIXING ITERATIVO CON PESOS (metodo de A. Revel, tesis).
//   Corrige la "correlacion residual" <C>(p) que sobrevive al mixing estandar cuando
//   la correlacion es fuerte (ec. 2.10-2.13 de la tesis).
//
//   Idea:
//     - Se construye el fondo mezclado sobre TODOS los pares i!=j (O(N^2)).
//     - C(x) = datos / fondo  ->  peso por particula w_i = 1/<C>(p_i)
//     - Se rehace el fondo pesando cada par virtual por w_i*w_j, y se itera.
//
//   Uso:  root -l 'eventMixingIterative.C+("/ruta/data_23O.root")'
#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <vector>
#include <cmath>
#include <iostream>

void eventMixing(const char *inName =
                     "/nucl_lustre/pablogrusell/g249/g249_analysis/results/final/data_23O.root",
                 const char *outName = "23O_mixing_iterative.root",
                 int nIter = 25,                                    // iteraciones del algoritmo (Fig 2.3 usa ~10)
                 int nInnerIter = 3,                                // sub-iteraciones para <C>(p_i) (ec. 2.13)
                 int nbins = 150, double xlo = -1, double xhi = 20) // MeV, como fErel
{
    // --- constantes y offsets (23O1n.txt) ------------------------------------
    const double m_neut = 0.939565; // GeV

    const double fFragOffX = -0.0026286137;
    const double fFragOffY = -0.016158624;
    const double fNeuOffX = -0.0012187827;
    const double fNeuOffY = -0.018529642;
    const double fBetaMatch = -0.00140415;

    const double fxOff = -fFragOffX, fyOff = -fFragOffY;
    const double dxOff = -fNeuOffX, dyOff = -fNeuOffY;

    // --- leer arbol ----------------------------------------------------------
    TFile *f = TFile::Open(inName, "READ");
    if (!f || f->IsZombie())
    {
        std::cerr << "ERROR: no se pudo abrir " << inName << std::endl;
        return;
    }
    TTree *t = dynamic_cast<TTree *>(f->Get("FilterDataTree"));
    if (!t)
    {
        std::cerr << "ERROR: no se encontro FilterDataTree.\n";
        f->Close();
        return;
    }

    Double_t M_frag, beta_frag, beta_neu;
    Double_t px_frag, py_frag, pz_frag;
    Double_t px_neu, py_neu, pz_neu;

    t->SetBranchAddress("M_frag", &M_frag);
    t->SetBranchAddress("beta_frag", &beta_frag);
    t->SetBranchAddress("px_frag", &px_frag);
    t->SetBranchAddress("py_frag", &py_frag);
    t->SetBranchAddress("pz_frag", &pz_frag);
    t->SetBranchAddress("beta_neu", &beta_neu);
    t->SetBranchAddress("px_neu", &px_neu);
    t->SetBranchAddress("py_neu", &py_neu);
    t->SetBranchAddress("pz_neu", &pz_neu);

    struct FragKin
    {
        double beta, M, fx, fy;
    };
    struct NeuKin
    {
        double beta, dx, dy;
    };
    std::vector<FragKin> frags;
    std::vector<NeuKin> neus;

    const Long64_t Nall = t->GetEntries();
    frags.reserve(Nall);
    neus.reserve(Nall);

    for (Long64_t i = 0; i < Nall; ++i)
    {
        t->GetEntry(i);
        FragKin fk;
        fk.beta = beta_frag + fBetaMatch;
        fk.M = M_frag;
        fk.fx = px_frag / pz_frag + fxOff;
        fk.fy = py_frag / pz_frag + fyOff;
        NeuKin nk;
        nk.beta = beta_neu;
        nk.dx = px_neu / pz_neu + dxOff;
        nk.dy = py_neu / pz_neu + dyOff;
        frags.push_back(fk);
        neus.push_back(nk);
    }
    f->Close();

    const int N = (int)frags.size();
    if (N < 2)
    {
        std::cerr << "ERROR: pocos eventos (" << N << ").\n";
        return;
    }
    std::cout << "Eventos: " << N << "   (pares virtuales ~ " << (double)N * (N - 1) << ")\n";

    // --- Erel(fragmento i, neutron j): misma formula que getData() -----------
    auto erelOf = [&](int i, int j) -> double
    {
        const FragKin &fk = frags[i];
        const NeuKin &nk = neus[j];
        double cos_ang =
            (nk.dx * fk.fx + nk.dy * fk.fy + 1.0) /
            (std::sqrt(nk.dx * nk.dx + nk.dy * nk.dy + 1.0) *
             std::sqrt(fk.fx * fk.fx + fk.fy * fk.fy + 1.0));
        double g_neu = 1.0 / std::sqrt(1.0 - nk.beta * nk.beta);
        double g_frag = 1.0 / std::sqrt(1.0 - fk.beta * fk.beta);
        double m_f = fk.M;
        double Erel = std::sqrt(m_f * m_f + m_neut * m_neut +
                                2.0 * g_neu * g_frag * m_f * m_neut *
                                    (1.0 - nk.beta * fk.beta * cos_ang)) -
                      m_f - m_neut;
        return Erel * 1000.0; // MeV
    };

    // --- histograma de datos (mismo evento) ----------------------------------
    TH1F *hReal = new TH1F("hReal", "Erel data;E_{rel} [MeV];counts", nbins, xlo, xhi);
    hReal->Sumw2();
    for (int i = 0; i < N; ++i)
        hReal->Fill(erelOf(i, i));
    const double sReal = hReal->Integral();

    // --- utilidades de binning para C(x) -------------------------------------
    const double bw = (xhi - xlo) / nbins;
    auto binOf = [&](double x) -> int
    {
        int b = (int)((x - xlo) / bw);
        if (b < 0 || b >= nbins)
            return -1;
        return b;
    };

    // Precalcular Erel_ij? -> N^2 doubles puede ser mucho (3090^2*8B ~ 76 MB, OK).
    // Para N mas grandes, quitar el cache y recalcular en el bucle.
    std::vector<float> erelCache((size_t)N * N);
    for (int i = 0; i < N; ++i)
        for (int j = 0; j < N; ++j)
            erelCache[(size_t)i * N + j] = (float)erelOf(i, j);
    auto ERij = [&](int i, int j) -> double
    { return erelCache[(size_t)i * N + j]; };

    // C(x): funcion de correlacion por bin. Inicializada a 1.
    std::vector<double> C(nbins, 1.0);

    // pesos por particula: wFrag[i] = 1/<C>(frag_i), wNeu[j] = 1/<C>(neu_j)
    std::vector<double> wFrag(N, 1.0), wNeu(N, 1.0);

    // helper: valor de C en la Erel del par (i,j)
    auto Cij = [&](int i, int j) -> double
    {
        int b = binOf(ERij(i, j));
        return (b < 0) ? 1.0 : C[b];
    };

    TH1F *hMix = new TH1F("hMix", "Erel mixed (weighted);E_{rel} [MeV];counts", nbins, xlo, xhi);
    hMix->Sumw2();

    // guardaremos la evolucion para inspeccion
    std::vector<TH1F *> snaps;

    for (int it = 0; it < nIter; ++it)
    {
        // (1) Actualizar pesos por particula resolviendo el sistema acoplado (ec 2.13)
        //     <C>(p_i) = 1/(N-1) * sum_{j!=i} C(x_ij) * w_partner_j
        //     Para el fragmento i, el partner es el neutron j (peso wNeu[j]).
        //     Para el neutron j, el partner es el fragmento i (peso wFrag[i]).
        //     Iteramos nInnerIter veces el sistema acoplado.
        if (it == 0)
        {
            std::fill(wFrag.begin(), wFrag.end(), 1.0);
            std::fill(wNeu.begin(), wNeu.end(), 1.0);
        }
        else
        {
            for (int inner = 0; inner < nInnerIter; ++inner)
            {
                // <C> para cada fragmento i (partner neutron j, peso wNeu[j])
                std::vector<double> avgFrag(N, 0.0), avgNeu(N, 0.0);
                for (int i = 0; i < N; ++i)
                {
                    double acc = 0.0, wsum = 0.0;
                    for (int j = 0; j < N; ++j)
                    {
                        if (j == i)
                            continue;
                        acc += Cij(i, j) * wNeu[j];
                        wsum += wNeu[j];
                    }
                    avgFrag[i] = (wsum > 0) ? acc / wsum : 1.0;
                }
                // <C> para cada neutron j (partner fragmento i, peso wFrag[i])
                for (int j = 0; j < N; ++j)
                {
                    double acc = 0.0, wsum = 0.0;
                    for (int i = 0; i < N; ++i)
                    {
                        if (i == j)
                            continue;
                        acc += Cij(i, j) * wFrag[i];
                        wsum += wFrag[i];
                    }
                    avgNeu[j] = (wsum > 0) ? acc / wsum : 1.0;
                }
                for (int i = 0; i < N; ++i)
                    wFrag[i] = (avgFrag[i] > 1e-9) ? 1.0 / avgFrag[i] : 1.0;
                for (int j = 0; j < N; ++j)
                    wNeu[j] = (avgNeu[j] > 1e-9) ? 1.0 / avgNeu[j] : 1.0;
            }
        }

        // (2) Construir el fondo mezclado pesado: todos los pares i!=j, peso wFrag[i]*wNeu[j]
        hMix->Reset();
        for (int i = 0; i < N; ++i)
            for (int j = 0; j < N; ++j)
            {
                if (j == i)
                    continue;
                hMix->Fill(ERij(i, j), wFrag[i] * wNeu[j]);
            }

        // (3) Normalizar el fondo al area de los datos y actualizar C(x) = datos/fondo
        double sMix = hMix->Integral();
        if (sMix <= 0)
        {
            std::cerr << "AVISO: fondo vacio en iter " << it << "\n";
            break;
        }
        double norm = sReal / sMix;

        for (int b = 1; b <= nbins; ++b)
        {
            double d = hReal->GetBinContent(b);
            double m = hMix->GetBinContent(b) * norm;
            C[b - 1] = (m > 0) ? d / m : 1.0;
        }

        // snapshot del fondo (normalizado) para ver convergencia
        TH1F *snap = (TH1F *)hMix->Clone(Form("hMix_it%d", it));
        snap->Scale(norm);
        snap->SetDirectory(nullptr);
        snaps.push_back(snap);

        std::cout << "  iter " << it
                  << "  norm=" << norm
                  << "  C medio=" << [&]
        { double s=0; for(double c:C) s+=c; return s/nbins; }()
                  << std::endl;
    }

    // fondo final normalizado
    double sMix = hMix->Integral();
    TH1F *hMixFinal = (TH1F *)hMix->Clone("hMixNonResonant");
    hMixFinal->SetTitle("Non-resonant background (iterative mixing);E_{rel} [MeV];counts");
    if (sMix > 0)
        hMixFinal->Scale(sReal / sMix);

    // template area=1 (pdf) para el fit
    TH1F *hTemplate = (TH1F *)hMixFinal->Clone("hMixTemplate");
    hTemplate->SetTitle("Non-resonant template (area=1);E_{rel} [MeV];pdf");
    double area = hTemplate->Integral("width");
    if (area > 0)
        hTemplate->Scale(1.0 / area);

    // funcion de correlacion como histograma
    TH1F *hC = new TH1F("hCorr", "Correlation C(x)=data/bkg;E_{rel} [MeV];C", nbins, xlo, xhi);
    for (int b = 1; b <= nbins; ++b)
        hC->SetBinContent(b, C[b - 1]);

    // --- dibujar -------------------------------------------------------------
    TCanvas *c1 = new TCanvas("c1", "data vs non-resonant bkg", 900, 600);
    hReal->SetLineColor(kBlack);
    hReal->SetLineWidth(2);
    hMixFinal->SetLineColor(kRed);
    hMixFinal->SetLineWidth(2);
    hReal->Draw("hist");
    hMixFinal->Draw("hist same");
    auto *leg = new TLegend(0.6, 0.7, 0.88, 0.88);
    leg->AddEntry(hReal, "data", "l");
    leg->AddEntry(hMixFinal, "non-resonant (iter.)", "l");
    leg->Draw();

    TCanvas *c2 = new TCanvas("c2", "iteration convergence", 900, 600);
    for (size_t k = 0; k < snaps.size(); ++k)
    {
        snaps[k]->SetLineColor(k == 0 ? kBlack : (k + 1));
        snaps[k]->SetLineWidth(2);
        snaps[k]->Draw(k == 0 ? "hist" : "hist same");
    }

    TCanvas *c3 = new TCanvas("c3", "correlation function", 900, 600);
    hC->SetLineColor(kBlue);
    hC->SetLineWidth(2);
    hC->Draw("hist");

    // --- guardar -------------------------------------------------------------
    TFile *fout = new TFile(outName, "RECREATE");

    // arbol con los pares reales (mismo evento) que llenan hReal (peso=1)
    Float_t brErelReal;
    Double_t brWeightReal;
    TTree *tReal = new TTree("tReal", "Pares reales (mismo evento) que llenan hReal");
    tReal->Branch("Erel", &brErelReal, "Erel/F");
    tReal->Branch("weight", &brWeightReal, "weight/D");
    for (int i = 0; i < N; ++i)
    {
        brErelReal = (Float_t)erelOf(i, i);
        brWeightReal = 1.0;
        tReal->Fill();
    }

    // arbol con los pares mezclados pesados (i!=j) que llenan hMixFinal
    // (usa los pesos wFrag/wNeu de la ultima iteracion y el mismo factor
    //  de normalizacion sReal/sMix aplicado a hMixFinal)
    const double normFinal = (sMix > 0) ? sReal / sMix : 1.0;

    Float_t brErelMix;
    Double_t brWeightMix;
    TTree *tMix = new TTree("tMix", "Pares mezclados (i!=j) que llenan hMixNonResonant");
    tMix->Branch("Erel", &brErelMix, "Erel/F");
    tMix->Branch("weight", &brWeightMix, "weight/D");
    for (int i = 0; i < N; ++i)
        for (int j = 0; j < N; ++j)
        {
            if (j == i)
                continue;
            brErelMix = (Float_t)ERij(i, j);
            brWeightMix = wFrag[i] * wNeu[j] * normFinal;
            tMix->Fill();
        }

    hReal->Write();
    hMixFinal->Write(); // fondo no resonante (cuentas, normalizado a datos)
    hTemplate->Write(); // template area=1 para el fit
    hC->Write();        // funcion de correlacion final
    tReal->Write();     // eventos (Erel, weight) que reconstruyen hReal
    tMix->Write();      // eventos (Erel, weight) que reconstruyen hMixFinal
    for (auto *s : snaps)
        s->Write();
    fout->Close();

    std::cout << "Guardado en: " << outName
              << "  (hReal, hMixNonResonant, hMixTemplate, hCorr, hMix_itN, tReal, tMix)\n";
}