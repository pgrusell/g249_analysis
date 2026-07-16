#include <vector>
#include <cmath>
#include <iostream>

void getNR(const char *inName = "/nucl_lustre/g249/sim/p_F25_630_v1.root",
           const char *outName = "p_F25_630_v1_p2p_24O_to_23On.root")
{
    // --- Canal: (p,2p) 25F -> 24O* -> 23O + n ---------------------------------
    const Short_t Z_frag = 8, A_frag = 23, S_frag = 0; // fragmento frio final: 23O
    const Short_t Z_rem = 8, A_rem = 24, S_rem = 0;    // remanente caliente: 24O
    const Int_t PDG_NEUTRON = 2112;
    const Int_t PDG_PROTON = 2212;

    // --- Fichero de entrada ---------------------------------------------------
    TFile *fin = TFile::Open(inName, "READ");
    if (!fin || fin->IsZombie())
    {
        std::cerr << "ERROR: no se pudo abrir " << inName << std::endl;
        return;
    }
    TTree *tin = nullptr;
    fin->GetObject("et", tin);
    if (!tin)
    {
        std::cerr << "ERROR: no se encontro el arbol 'et'" << std::endl;
        fin->Close();
        return;
    }

    // --- Ramas de PARTICULAS emitidas (estado final, indexadas por nParticles) -
    const Int_t MAXP = 5000;
    Short_t nParticles = 0;
    Short_t A[MAXP], Z[MAXP], S[MAXP], origin[MAXP];
    Int_t PDGCode[MAXP];
    Float_t px[MAXP], py[MAXP], pz[MAXP], EKin[MAXP];

    tin->SetBranchAddress("nParticles", &nParticles);
    tin->SetBranchAddress("A", A);
    tin->SetBranchAddress("Z", Z);
    tin->SetBranchAddress("S", S);
    tin->SetBranchAddress("PDGCode", PDGCode);
    tin->SetBranchAddress("px", px);
    tin->SetBranchAddress("py", py);
    tin->SetBranchAddress("pz", pz);
    tin->SetBranchAddress("EKin", EKin);
    tin->SetBranchAddress("origin", origin); // origen de cada particula (cascada/evapor/...)

    // --- Ramas de REMANENTES (prefragmento caliente, antes de la desexcitacion) -
    const Int_t MAXREM = 256;
    Short_t nRemnants = 0;
    Short_t ARem[MAXREM], ZRem[MAXREM], SRem[MAXREM];

    tin->SetBranchAddress("nRemnants", &nRemnants);
    tin->SetBranchAddress("ARem", ARem);
    tin->SetBranchAddress("ZRem", ZRem);
    tin->SetBranchAddress("SRem", SRem);

    // --- Fichero de salida ----------------------------------------------------
    TFile *fout = new TFile(outName, "RECREATE");
    if (!fout || fout->IsZombie())
    {
        std::cerr << "ERROR: no se pudo crear " << outName << std::endl;
        fin->Close();
        return;
    }

    TTree *tout = tin->CloneTree(0);

    std::vector<float> erelFN; // Erel(23O + n) por cada neutron
    tout->Branch("erelFN", &erelFN);
    std::vector<short> nOrigin; // origin del neutron, ALINEADO con erelFN
    tout->Branch("nOrigin", &nOrigin);
    float erelFNall = -1.f;
    tout->Branch("erelFNall", &erelFNall);
    Int_t nNeutrons = 0, nProtons = 0;
    tout->Branch("nNeutrons", &nNeutrons, "nNeutrons/I");
    tout->Branch("nProtons", &nProtons, "nProtons/I");

    // --- Helper: cuadrivector a partir de EKin y (px,py,pz) -------------------
    auto fourMom = [](float ekin, float ppx, float ppy, float ppz, double &E, double &m)
    {
        double p2 = (double)ppx * ppx + (double)ppy * ppy + (double)ppz * ppz;
        m = (p2 - (double)ekin * ekin) / (2.0 * (double)ekin);
        E = (double)ekin + m;
    };

    const Long64_t nEntries = tin->GetEntries();
    Long64_t nKept = 0;
    std::cout << "Procesando " << nEntries << " eventos..." << std::endl;

    for (Long64_t i = 0; i < nEntries; ++i)
    {
        tin->GetEntry(i);

        // (A) SELECCION (p,2p)24O: el remanente de la cascada debe ser 24O.
        bool has24Orem = false;
        for (Short_t r = 0; r < nRemnants; ++r)
            if (ZRem[r] == Z_rem && ARem[r] == A_rem && SRem[r] == S_rem)
            {
                has24Orem = true;
                break;
            }
        if (!has24Orem)
            continue;

        // (B) el fragmento frio final debe ser 23O.
        Int_t idxFrag = -1;
        for (Short_t k = 0; k < nParticles; ++k)
            if (Z[k] == Z_frag && A[k] == A_frag && S[k] == S_frag)
            {
                idxFrag = k;
                break;
            }
        if (idxFrag < 0)
            continue;
        if (EKin[idxFrag] <= 0.f)
            continue;

        // --- cuadrivector del 23O ---
        double Ef, mf;
        fourMom(EKin[idxFrag], px[idxFrag], py[idxFrag], pz[idxFrag], Ef, mf);
        const double Pfx = px[idxFrag], Pfy = py[idxFrag], Pfz = pz[idxFrag];

        // --- recorrer neutrones y protones ---
        erelFN.clear();
        nOrigin.clear();
        double Esum = Ef, Pxs = Pfx, Pys = Pfy, Pzs = Pfz, mNsum = 0.0;
        nNeutrons = 0;
        nProtons = 0;

        for (Short_t k = 0; k < nParticles; ++k)
        {
            if (PDGCode[k] == PDG_PROTON)
            {
                ++nProtons;
                continue;
            }
            if (k == idxFrag)
                continue;
            if (PDGCode[k] != PDG_NEUTRON)
                continue;
            if (EKin[k] <= 0.f)
                continue;

            double En, mn;
            fourMom(EKin[k], px[k], py[k], pz[k], En, mn);

            double Et = Ef + En;
            double Pxt = Pfx + px[k], Pyt = Pfy + py[k], Pzt = Pfz + pz[k];
            double Minv = std::sqrt(Et * Et - (Pxt * Pxt + Pyt * Pyt + Pzt * Pzt));
            erelFN.push_back((float)(Minv - mf - mn));
            nOrigin.push_back(origin[k]);

            Esum += En;
            Pxs += px[k];
            Pys += py[k];
            Pzs += pz[k];
            mNsum += mn;
            ++nNeutrons;
        }

        if (nNeutrons > 0)
        {
            double Minv = std::sqrt(Esum * Esum - (Pxs * Pxs + Pys * Pys + Pzs * Pzs));
            erelFNall = (float)(Minv - mf - mNsum);
        }
        else
            erelFNall = -1.f;

        tout->Fill();
        ++nKept;

        if (i > 0 && (i % 500000) == 0)
            std::cout << "  " << i << " / " << nEntries << "   (seleccionados: " << nKept << ")" << std::endl;
    }

    fout->cd();
    tout->Write();

    std::cout << "-----------------------------------------------------------\n";
    std::cout << "Eventos (p,2p) 24O -> 23O + n: " << nKept << " de " << nEntries << std::endl;
    std::cout << "Guardado en: " << outName << "  (ramas nuevas: erelFN[], nOrigin[], erelFNall, nNeutrons, nProtons)" << std::endl;

    fout->Close();
    fin->Close();
}