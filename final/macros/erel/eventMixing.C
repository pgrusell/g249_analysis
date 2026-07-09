#include <TH1.h>
#include <TRotation.h>
#include <TVector3.h>

namespace
{
    // Masa del neutron en GeV (mismo valor que DataAnalysis::m_neut)
    constexpr double m_neut = 0.939565;

    struct Offsets
    {
        double fragP[2] = {0., 0.}; // offsets angulares px/pz, py/pz del fragmento
        double neuP[2] = {0., 0.};  // offsets angulares px/pz, py/pz del neutron
        double betaMatch = 0.;      // correccion aditiva de beta del fragmento
        bool ok = false;
    };

    // Cinematica minima que necesitamos para el Erel, cacheada por evento.
    // Se guardan tanto las componentes de laboratorio (px,py,pz) como su version
    // rotada evento-a-evento al marco del proyectil de ESE evento (pxR,pyR,pzR),
    // ver rotateToBeamFrame() y la nota fisica junto a computeCosAngle().
    struct Frag
    {
        double M, beta;
        double px, py, pz;
        double pxR, pyR, pzR;
    };
    struct Neu
    {
        double beta;
        double px, py, pz;
        double pxR, pyR, pzR;
    };

    // Misma convencion de lectura que DataAnalysis::setOffsetsFromTxt
    Offsets loadOffsets(const std::string &offFile)
    {
        Offsets o;

        const char *repo = getenv("repopath");
        if (!repo)
        {
            std::cerr << "[eventMixing] ERROR: la variable de entorno 'repopath' no esta definida\n";
            return o;
        }

        const std::string txtPath = std::string(repo) + "/final/settings/" + offFile;

        std::ifstream in(txtPath);
        if (!in)
        {
            std::cerr << "[eventMixing] ERROR: no se puede abrir el fichero de offsets " << txtPath << "\n";
            return o;
        }

        std::vector<double> v;
        std::string line;
        while (std::getline(in, line))
            v.push_back(std::atof(line.c_str()));

        if (v.size() < 5)
        {
            std::cerr << "[eventMixing] ERROR: el fichero de offsets tiene < 5 lineas\n";
            return o;
        }

        o.fragP[0] = v[0];
        o.fragP[1] = v[1];
        o.neuP[0] = v[2];
        o.neuP[1] = v[3];
        o.betaMatch = v[4];
        o.ok = true;
        return o;
    }

    // Rota un vector de momento al marco en el que la direccion del haz (por evento)
    // define el eje Z, con la misma convencion que DataAnalysis::getData(): phi/theta
    // de la direccion del haz, R.RotateZ(-phi) seguido de R.RotateY(-theta).
    TVector3 rotateToBeamFrame(const TVector3 &v, const TVector3 &beamDir)
    {
        const double phi = beamDir.Phi();
        const double theta = beamDir.Theta();

        TRotation R;
        R.RotateZ(-phi);
        R.RotateY(-theta);

        TVector3 out(v);
        out.Transform(R);
        return out;
    }

    // Coseno del angulo de apertura neutron-fragmento (misma formula que
    // DataAnalysis::getData(), proyectando las pendientes transversales px/pz, py/pz).
    //
    // Nota fisica sobre 'rotateToProjFrame':
    // En un evento real, el neutron y el fragmento comparten la direccion del haz
    // incidente DE ESE EVENTO (rama px_in/py_in/pz_in), que fluctua evento a evento
    // por la optica/reconstruccion del haz. Esa inclinacion comun se cancela en el
    // angulo relativo neutron-fragmento cuando ambos vienen del MISMO evento. Si en
    // vez de eso usamos las pendientes de LABORATORIO (con offsets constantes,
    // promediados sobre todo el run) y mezclamos eventos (event mixing), el neutron
    // trae la inclinacion de haz del evento A y el fragmento la del evento B: esa
    // diferencia YA NO se cancela y se suma como un angulo de apertura espurio,
    // inflando el Erel del fondo mezclado (maximo desplazado a 5-6 MeV con cola
    // hasta 20 MeV en vez de picar cerca del umbral).
    //
    // Rotando cada evento a su propio marco del proyectil ANTES de cachear/mezclar
    // (misma maquinaria que la rotacion a projectile-RF de DataAnalysis::getData():
    // TRotation con RotateZ(-phi)/RotateY(-theta) segun la direccion del haz de ESE
    // evento) se elimina esa componente espuria: tanto el neutron como el fragmento
    // quedan expresados en un marco donde el haz de SU evento ya es el eje Z, y el
    // angulo relativo que sobrevive al mezclar es el fisicamente comparable.
    //
    // 'applyLabOffsets' controla si ademas se siguen restando los offsets angulares
    // constantes (fNeutronPOffsets/fFragmentPOffsets). Esos offsets se calibraron
    // originalmente para corregir un desalineamiento medio en el marco de LABORATORIO;
    // en el marco rotado podrian estar corrigiendo (parcial o totalmente) lo mismo que
    // ya corrige la rotacion evento-a-evento, y aplicarlos ahi podria sobre-corregir.
    // Se deja como opcion para comparar ambos casos en vez de asumir cual es correcto.
    double computeCosAngle(const Frag &fr, const Neu &nu, const Offsets &off,
                            bool rotateToProjFrame, bool applyLabOffsets)
    {
        double npx, npy, npz;
        double fpx, fpy, fpz;

        if (rotateToProjFrame)
        {
            npx = nu.pxR;
            npy = nu.pyR;
            npz = nu.pzR;
            fpx = fr.pxR;
            fpy = fr.pyR;
            fpz = fr.pzR;
        }
        else
        {
            npx = nu.px;
            npy = nu.py;
            npz = nu.pz;
            fpx = fr.px;
            fpy = fr.py;
            fpz = fr.pz;
        }

        const double offNx = applyLabOffsets ? off.neuP[0] : 0.0;
        const double offNy = applyLabOffsets ? off.neuP[1] : 0.0;
        const double offFx = applyLabOffsets ? off.fragP[0] : 0.0;
        const double offFy = applyLabOffsets ? off.fragP[1] : 0.0;

        const double dx_neu = npx / npz - offNx;
        const double dy_neu = npy / npz - offNy;

        const double fx_frag = fpx / fpz - offFx;
        const double fy_frag = fpy / fpz - offFy;

        return (dx_neu * fx_frag + dy_neu * fy_frag + 1.0) /
               (std::sqrt(dx_neu * dx_neu + dy_neu * dy_neu + 1.0) *
                std::sqrt(fx_frag * fx_frag + fy_frag * fy_frag + 1.0));
    }

    // Misma formula de masa invariante que en DataAnalysis::getData().
    // El fragmento (fr) y el neutron (nu) pueden provenir de eventos DISTINTOS
    // (event mixing) o del MISMO evento (Erel correlacionado).
    double computeErel(const Frag &fr, const Neu &nu, const Offsets &off,
                        bool rotateToProjFrame, bool applyLabOffsets)
    {
        // Offset de beta aplicado al fragmento (igual que beta_frag += fBetaMatchValue)
        const double beta_frag = fr.beta + off.betaMatch;

        const double cos_ang = computeCosAngle(fr, nu, off, rotateToProjFrame, applyLabOffsets);

        const double gamma_neu = 1.0 / std::sqrt(1.0 - nu.beta * nu.beta);
        const double gamma_frag = 1.0 / std::sqrt(1.0 - beta_frag * beta_frag);

        const double m_f = fr.M;

        const double Erel =
            std::sqrt(m_f * m_f + m_neut * m_neut +
                      2.0 * gamma_neu * gamma_frag * m_f * m_neut *
                          (1.0 - nu.beta * beta_frag * cos_ang)) -
            m_f - m_neut;

        return Erel;
    }

    // Angulo de apertura en grados, a partir del cos_ang usado en el Erel (con clamp
    // por seguridad numerica antes de acos).
    double openingAngleDeg(const Frag &fr, const Neu &nu, const Offsets &off,
                           bool rotateToProjFrame, bool applyLabOffsets)
    {
        double c = computeCosAngle(fr, nu, off, rotateToProjFrame, applyLabOffsets);
        if (c > 1.0)
            c = 1.0;
        if (c < -1.0)
            c = -1.0;
        return TMath::ACos(c) * TMath::RadToDeg();
    }
}

void eventMixing(
    TString dataFilePath = "/nucl_lustre/pablogrusell/g249/g249_analysis/results/dataFiles/data_23O.root",
    TString offFile = "23O1n.txt",
    TString outFileName = "23O_background_mixed.root",
    int nCounts = 6000,
    unsigned int seed = -1,
    int nbins = 100, double elo = 0, double ehi = 20.,
    bool rotateToProjFrame = true, // rotar cada evento a su propio marco del haz antes de cachear/mezclar
    bool applyLabOffsets = true)   // aplicar ademas los offsets angulares constantes (ver nota en computeCosAngle)
{
    // ---------------------------------------------------------------------
    // 1) Offsets (momento + beta), misma convencion que la clase
    // ---------------------------------------------------------------------
    Offsets off = loadOffsets(offFile.Data());
    if (!off.ok)
        return;

    // ---------------------------------------------------------------------
    // 2) Abrir fichero y TTree
    // ---------------------------------------------------------------------
    TFile *f = TFile::Open(dataFilePath, "READ");
    if (!f || f->IsZombie())
    {
        std::cerr << "[eventMixing] ERROR: no se pudo abrir " << dataFilePath << "\n";
        return;
    }

    TTree *t = dynamic_cast<TTree *>(f->Get("FilterDataTree"));
    if (!t)
    {
        std::cerr << "[eventMixing] ERROR: no se encontro el TTree \"FilterDataTree\"\n";
        f->Close();
        return;
    }

    Double_t M_frag, beta_frag, px_frag, py_frag, pz_frag;
    Double_t beta_neu, px_neu, py_neu, pz_neu;
    Double_t califa_opa;
    Double_t px_in, py_in, pz_in; // direccion del haz incidente, por evento

    // Solo activamos las ramas que usamos (mas rapido)
    t->SetBranchStatus("*", 0);
    const char *used[] = {"M_frag", "beta_frag", "px_frag", "py_frag", "pz_frag",
                          "beta_neu", "px_neu", "py_neu", "pz_neu", "califa_opa",
                          "px_in", "py_in", "pz_in"};
    for (auto br : used)
        t->SetBranchStatus(br, 1);

    t->SetBranchAddress("M_frag", &M_frag);
    t->SetBranchAddress("beta_frag", &beta_frag);
    t->SetBranchAddress("px_frag", &px_frag);
    t->SetBranchAddress("py_frag", &py_frag);
    t->SetBranchAddress("pz_frag", &pz_frag);
    t->SetBranchAddress("beta_neu", &beta_neu);
    t->SetBranchAddress("px_neu", &px_neu);
    t->SetBranchAddress("py_neu", &py_neu);
    t->SetBranchAddress("pz_neu", &pz_neu);
    t->SetBranchAddress("califa_opa", &califa_opa);
    t->SetBranchAddress("px_in", &px_in);
    t->SetBranchAddress("py_in", &py_in);
    t->SetBranchAddress("pz_in", &pz_in);

    // ---------------------------------------------------------------------
    // 3) Cachear la cinematica evento a evento (frag + neu del MISMO evento).
    //    Guardamos ambos alineados por indice para poder exigir luego iF != iN.
    //    Si rotateToProjFrame, tambien cacheamos las componentes rotadas al
    //    marco del haz de ESE evento (ver rotateToBeamFrame()).
    // ---------------------------------------------------------------------
    std::vector<Frag> frags;
    std::vector<Neu> neus;

    // Distribuciones de control de las variables de entrada
    TH1F *hBetaFrag = new TH1F("hBetaFrag", "#beta del fragmento (QFS, cacheados);#beta_{frag};Cuentas", 100, 0., 1.);
    TH1F *hBetaNeu = new TH1F("hBetaNeu", "#beta del neutron (QFS, cacheados);#beta_{neu};Cuentas", 100, 0., 1.);
    TH1F *hMFrag = new TH1F("hMFrag", "Masa del fragmento (QFS, cacheados);M_{frag} [GeV];Cuentas", 100, 0., 30.);
    hBetaFrag->SetCanExtend(TH1::kXaxis);
    hBetaNeu->SetCanExtend(TH1::kXaxis);
    hMFrag->SetCanExtend(TH1::kXaxis);

    const Long64_t nentries = t->GetEntries();
    frags.reserve(nentries);
    neus.reserve(nentries);

    for (Long64_t i = 0; i < nentries; ++i)
    {
        t->GetEntry(i);

        // Filtro minimo de validez para evitar NaN/inf
        const bool okFrag = (pz_frag != 0.) && (M_frag > 0.) &&
                            (beta_frag > 0.) && (beta_frag < 1.);
        const bool okNeu = (pz_neu != 0.) &&
                           (beta_neu > 0.) && (beta_neu < 1.);
        // Solo eventos QFS (quasi-free scattering), definidos por el angulo de apertura de CALIFA
        const bool okOpa = (califa_opa > 1.25) && (califa_opa < 1.65);

        if (!okFrag || !okNeu || !okOpa)
            continue;

        Frag fr;
        fr.M = M_frag;
        fr.beta = beta_frag;
        fr.px = px_frag;
        fr.py = py_frag;
        fr.pz = pz_frag;
        fr.pxR = fr.px;
        fr.pyR = fr.py;
        fr.pzR = fr.pz;

        Neu nu;
        nu.beta = beta_neu;
        nu.px = px_neu;
        nu.py = py_neu;
        nu.pz = pz_neu;
        nu.pxR = nu.px;
        nu.pyR = nu.py;
        nu.pzR = nu.pz;

        if (rotateToProjFrame)
        {
            const TVector3 beamDir(px_in, py_in, pz_in);
            if (beamDir.Mag2() <= 0.)
                continue; // sin direccion de haz valida para este evento, no se puede rotar

            const TVector3 vFrag = rotateToBeamFrame(TVector3(px_frag, py_frag, pz_frag), beamDir);
            const TVector3 vNeu = rotateToBeamFrame(TVector3(px_neu, py_neu, pz_neu), beamDir);

            fr.pxR = vFrag.X();
            fr.pyR = vFrag.Y();
            fr.pzR = vFrag.Z();
            nu.pxR = vNeu.X();
            nu.pyR = vNeu.Y();
            nu.pzR = vNeu.Z();
        }

        frags.push_back(fr);
        neus.push_back(nu);

        hBetaFrag->Fill(beta_frag);
        hBetaNeu->Fill(beta_neu);
        hMFrag->Fill(M_frag);
    }

    const size_t N = frags.size();
    std::cout << "[eventMixing] Eventos validos cacheados: " << N << "\n";
    if (N < 2)
    {
        std::cerr << "[eventMixing] ERROR: no hay suficientes eventos para mezclar\n";
        f->Close();
        return;
    }

    // ---------------------------------------------------------------------
    // 3b) Erel correcto: fragmento y neutron del MISMO evento (sin mixing).
    //     Tambien llenamos el angulo de apertura para pares del MISMO evento,
    //     como referencia para comparar con el de los pares MEZCLADOS.
    // ---------------------------------------------------------------------
    TH1F *hErel = new TH1F(
        "hErel",
        "E_{rel} correlacionado (mismo evento);E_{rel} [MeV];Cuentas",
        nbins, elo, ehi);
    hErel->SetLineColor(kBlue + 1);
    hErel->SetLineWidth(2);

    TH1F *hOpaSame = new TH1F(
        "hOpaSame",
        "Angulo de apertura neutron-fragmento;#theta_{n-frag} [deg];Cuentas",
        100, 0., 5.);
    hOpaSame->SetCanExtend(TH1::kXaxis);
    hOpaSame->SetLineColor(kBlue + 1);
    hOpaSame->SetLineWidth(2);

    for (size_t i = 0; i < N; ++i)
    {
        const double Erel = computeErel(frags[i], neus[i], off, rotateToProjFrame, applyLabOffsets);
        if (!std::isfinite(Erel))
            continue;

        hErel->Fill(Erel * 1000.0); // GeV -> MeV
        hOpaSame->Fill(openingAngleDeg(frags[i], neus[i], off, rotateToProjFrame, applyLabOffsets));
    }

    std::cout << "[eventMixing] Erel correcto (mismo evento) llenado con "
              << hErel->GetEntries() << " cuentas\n";

    // ---------------------------------------------------------------------
    // 4) Event mixing: neutron del evento iN + fragmento del evento iF (iF != iN),
    //    hasta acumular nCounts entradas en el histograma.
    // ---------------------------------------------------------------------
    TH1F *hMix = new TH1F(
        "hErelMixed",
        "Fondo no resonante (event mixing);E_{rel} [MeV];Cuentas",
        nbins, elo, ehi);
    hMix->SetLineColor(kRed + 1);
    hMix->SetLineWidth(2);

    TH1F *hOpaMixed = new TH1F(
        "hOpaMixed",
        "Angulo de apertura neutron-fragmento;#theta_{n-frag} [deg];Cuentas",
        100, 0., 5.);
    hOpaMixed->SetCanExtend(TH1::kXaxis);
    hOpaMixed->SetLineColor(kRed + 1);
    hOpaMixed->SetLineWidth(2);

    TRandom3 rng(seed);

    const Long64_t maxAttempts = 500LL * nCounts + 1000000LL; // salvaguarda anti-bucle infinito
    Long64_t attempts = 0;
    Long64_t filled = 0;

    while (hMix->GetEntries() < nCounts && attempts < maxAttempts)
    {
        ++attempts;

        const size_t iN = rng.Integer(N); // evento del que tomamos el NEUTRON
        size_t iF = rng.Integer(N);       // evento del que tomamos el FRAGMENTO
        if (iF == iN)
            continue; // deben ser eventos distintos -> descorrelacion

        const double Erel = computeErel(frags[iF], neus[iN], off, rotateToProjFrame, applyLabOffsets);
        if (!std::isfinite(Erel))
            continue;

        hMix->Fill(Erel * 1000.0); // GeV -> MeV (igual que getData)
        hOpaMixed->Fill(openingAngleDeg(frags[iF], neus[iN], off, rotateToProjFrame, applyLabOffsets));
        ++filled;
    }

    std::cout << "[eventMixing] Pares mezclados aceptados: " << filled
              << "  (intentos: " << attempts << ")\n";
    if (hMix->GetEntries() < nCounts)
        std::cerr << "[eventMixing] AVISO: solo se alcanzaron " << hMix->GetEntries()
                  << " cuentas antes del limite de intentos\n";

    // ---------------------------------------------------------------------
    // 5) Guardar y dibujar
    // ---------------------------------------------------------------------
    const std::string outPath =
        std::string(getenv("repopath")) + "/results/final/" + outFileName.Data();

    TFile *fout = new TFile(outPath.c_str(), "RECREATE");
    hMix->Write();
    hErel->Write();
    hOpaSame->Write();
    hOpaMixed->Write();
    hBetaFrag->Write();
    hBetaNeu->Write();
    hMFrag->Write();
    fout->Close();
    std::cout << "[eventMixing] Histogramas guardados en: " << outPath << "\n";

    TCanvas *c = new TCanvas("cMix", "Event mixing background", 800, 600);
    hErel->Draw("HIST");
    hMix->Draw("HIST SAME");
    c->Update();

    TCanvas *cOpa = new TCanvas("cOpa", "Angulo de apertura: mismo evento vs mezclado", 800, 600);
    hOpaSame->Draw("HIST");
    hOpaMixed->Draw("HIST SAME");
    cOpa->Update();

    TCanvas *cInputs = new TCanvas("cInputs", "Distribuciones de entrada (QFS)", 1200, 400);
    cInputs->Divide(3, 1);
    cInputs->cd(1);
    hBetaFrag->Draw("HIST");
    cInputs->cd(2);
    hBetaNeu->Draw("HIST");
    cInputs->cd(3);
    hMFrag->Draw("HIST");
    cInputs->Update();

    f->Close();
}
