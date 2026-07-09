void getNR(const char *inName = "/nucl_lustre/g249/sim/p_F25_630_v1.root",
           const char *outName = "p_F25_630_v1_23O.root")
{
    // --- Nucleo objetivo: 23O -------------------------------------------------
    const Short_t Z_target = 8;  // oxigeno
    const Short_t A_target = 23; // numero masico
    const Short_t S_target = 0;  // extraneza nula (23O es un nucleo normal)

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
        std::cerr << "ERROR: no se encontro el arbol 'et' en el fichero" << std::endl;
        fin->Close();
        return;
    }

    // --- Variables de seleccion (se enlazan ANTES de clonar) ------------------
    // nRemnants es muy pequeno; 256 es de sobra.
    const Int_t MAXREM = 256;
    Short_t nRemnants = 0;
    Short_t ARem[MAXREM];
    Short_t ZRem[MAXREM];
    Short_t SRem[MAXREM];

    tin->SetBranchAddress("nParticles", &nRemnants);
    tin->SetBranchAddress("A", ARem);
    tin->SetBranchAddress("Z", ZRem);
    tin->SetBranchAddress("S", SRem);

    // --- Fichero de salida (crearlo antes de CloneTree) -----------------------
    TFile *fout = new TFile(outName, "RECREATE");
    if (!fout || fout->IsZombie())
    {
        std::cerr << "ERROR: no se pudo crear " << outName << std::endl;
        fin->Close();
        return;
    }

    // Clona la estructura del arbol con 0 entradas. Las ramas de seleccion
    // comparten los buffers definidos arriba; el resto usa buffers propios.
    TTree *tout = tin->CloneTree(0);

    // --- Bucle sobre los eventos ----------------------------------------------
    const Long64_t nEntries = tin->GetEntries();
    Long64_t nKept = 0;

    std::cout << "Procesando " << nEntries << " eventos..." << std::endl;

    for (Long64_t i = 0; i < nEntries; ++i)
    {
        tin->GetEntry(i);

        bool keep = false;
        Short_t nr = nRemnants;
        for (Short_t r = 0; r < nr; ++r)
        {
            if (ZRem[r] == Z_target && ARem[r] == A_target && SRem[r] == S_target)
            {
                keep = true;
                break;
            }
        }

        if (keep)
        {
            tout->Fill();
            ++nKept;
        }

        if (i > 0 && (i % 500000) == 0)
            std::cout << "  " << i << " / " << nEntries
                      << "   (seleccionados hasta ahora: " << nKept << ")" << std::endl;
    }

    // --- Guardar y cerrar -----------------------------------------------------
    fout->cd();
    tout->Write();

    std::cout << "-----------------------------------------------------------\n";
    std::cout << "Eventos con remanente 23O (Z=8, A=23, S=0): "
              << nKept << " de " << nEntries << std::endl;
    std::cout << "Arbol filtrado guardado en: " << outName << std::endl;

    fout->Close();
    fin->Close();
}