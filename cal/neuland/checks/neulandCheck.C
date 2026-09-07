void neulandCheck(std::string inFilePath = "/nucl_lustre/pablogrusell/g249/root_files/unpackedData/g249_all_det_offline_0001_20260223_132131_160_0.root",
                  std::string outFilePath = "out_reduced.root")
{
    auto *f = TFile::Open(inFilePath.c_str(), "READ");
    if (f == nullptr || f->IsZombie())
    {
        std::cerr << "ERROR: no se puede abrir el fichero de entrada: " << inFilePath << std::endl;
        gSystem->Exit(1);
    }

    auto *evt = dynamic_cast<TTree *>(f->Get("evt"));
    if (evt == nullptr)
    {
        std::cerr << "ERROR: no hay TTree 'evt' en " << inFilePath << std::endl;
        f->Close();
        gSystem->Exit(2);
    }

    auto *neulandData = new TClonesArray("R3BNeulandHit");
    evt->SetBranchAddress("NeulandHits", &neulandData);

    auto *frsData = new TClonesArray("R3BFrsData");
    evt->SetBranchAddress("FrsData", &frsData);

    auto *outFile = new TFile(outFilePath.c_str(), "RECREATE");
    if (outFile == nullptr || outFile->IsZombie())
    {
        std::cerr << "ERROR: no se puede crear el fichero de salida: " << outFilePath << std::endl;
        gSystem->Exit(3);
    }
    outFile->cd();

    auto *outTree = new TTree("outtree", "");

    double xneu, yneu, zneu, tneu;
    double tof;
    int paddle;

    outTree->Branch("xneu", &xneu);
    outTree->Branch("yneu", &yneu);
    outTree->Branch("zneu", &zneu);
    outTree->Branch("tneu", &tneu);
    outTree->Branch("paddle", &paddle);
    outTree->Branch("tof", &tof);

    const Long64_t nEntries = evt->GetEntries();
    std::cout << "Entradas: " << nEntries << std::endl;

    // 25F cut
    TCutG *cutg_incoming = new TCutG("cut_25f_incoming", 14);
    cutg_incoming->SetPoint(0, 8.23352, 2.7854);
    cutg_incoming->SetPoint(1, 8.12607, 2.77784);
    cutg_incoming->SetPoint(2, 8.32665, 2.77145);
    cutg_incoming->SetPoint(3, 8.72779, 2.76862);
    cutg_incoming->SetPoint(4, 9.24355, 2.76626);
    cutg_incoming->SetPoint(5, 9.80229, 2.76744);
    cutg_incoming->SetPoint(6, 10.2321, 2.77004);
    cutg_incoming->SetPoint(7, 10.3324, 2.77547);
    cutg_incoming->SetPoint(8, 9.9384, 2.7828);
    cutg_incoming->SetPoint(9, 9.51576, 2.78753);
    cutg_incoming->SetPoint(10, 9.04298, 2.79084);
    cutg_incoming->SetPoint(11, 8.67049, 2.78942);
    cutg_incoming->SetPoint(12, 8.31948, 2.78611);
    cutg_incoming->SetPoint(13, 8.23352, 2.7854);

    for (Long64_t i = 0; i < nEntries; i++)
    {
        evt->GetEntry(i);

        if (i % 100000 == 0)
            std::cout << i << "/" << nEntries << std::endl;

        // Solo 25F
        if (frsData->GetEntries() < 1)
            continue;

        auto *frshit = static_cast<R3BFrsData *>(frsData->At(0));

        const double Z = frshit->GetZ();
        const double Aq = frshit->GetAq();

        if (!cutg_incoming->IsInside(Z, Aq))
            continue;

        // Primer hit en NeuLAND (en posicion)
        const int nNeulandHits = neulandData->GetEntries();
        if (nNeulandHits == 0)
            continue;

        zneu = 99999;
        xneu = 99999;
        yneu = 99999;
        tneu = 99999;
        paddle = 9999;

        for (int j = 0; j < nNeulandHits; j++)
        {
            auto *neuhit = static_cast<R3BNeulandHit *>(neulandData->At(j));

            const double zNew = neuhit->GetPosition().Z();

            if (zNew < zneu)
            {
                zneu = zNew;
                xneu = neuhit->GetPosition().X();
                yneu = neuhit->GetPosition().Y();
                tneu = neuhit->GetT();
                paddle = neuhit->GetPaddle();
            }
        }

        tof = tneu / (TMath::Sqrt(zneu * zneu + xneu * xneu + yneu * yneu)) * 1557.0;

        outTree->Fill();
    }

    outFile->cd();
    outTree->Write();
    outFile->Close();

    f->Close();

    std::cout << "Escrito: " << outFilePath << std::endl;
}