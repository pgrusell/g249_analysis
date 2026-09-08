void neulandCheck(std::string inFilePath = "/nucl_lustre/pablogrusell/g249/root_files/unpackedData/g249_all_det_offline_0001_20260223_132131_160_0.root",
                  std::string outFilePath = "out_reduced.root", bool offsetsVal = true)
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

    // tof offsets

    std::vector<double> offsets(1300);

    if (offsetsVal)
    {
        offsets[19 - 1] = 2.1280;
        offsets[30 - 1] = 1.3686;
        offsets[62 - 1] = -3.3266;
        offsets[72 - 1] = -3.3507;
        offsets[79 - 1] = -0.6977;
        offsets[88 - 1] = -3.1225;
        offsets[115 - 1] = 4.6680;
        offsets[166 - 1] = 2.5123;
        offsets[167 - 1] = -2.5644;
        offsets[168 - 1] = 4.5560;
        offsets[201 - 1] = 5.1455;
        offsets[233 - 1] = 2.9769;
        offsets[280 - 1] = 3.8662;
        offsets[295 - 1] = 2.7434;
        offsets[339 - 1] = 3.2547;
        offsets[371 - 1] = -3.4727;
        offsets[393 - 1] = 6.0279;
        offsets[456 - 1] = -2.7744;
        offsets[599 - 1] = 3.3694;
        offsets[622 - 1] = 5.0871;
        offsets[663 - 1] = 3.1136;
        offsets[765 - 1] = 7.2081;
        offsets[767 - 1] = 2.5812;
        offsets[875 - 1] = 7.2277;
        offsets[930 - 1] = 2.7249;
        offsets[934 - 1] = 1.8959;
        offsets[981 - 1] = 0.7595;
        offsets[1028 - 1] = 2.0740;
        offsets[1162 - 1] = 4.9011;
        offsets[1165 - 1] = 2.4801;
        offsets[1191 - 1] = 1.9692;
        offsets[1227 - 1] = -1.6468;
        offsets[1273 - 1] = 2.5166;
    }

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
                paddle = neuhit->GetPaddle();
                tneu = neuhit->GetT(); //+ offsets[paddle - 1];
            }
        }

        // tof = tneu / (TMath::Sqrt(zneu * zneu + xneu * xneu + yneu * yneu)) * 1557.0;
        double d = TMath::Sqrt(xneu * xneu + yneu * yneu + zneu * zneu);
        tof = tneu / d * 1557.0 + offsets[paddle];
        tneu += 1557. / d * offsets[paddle];

        outTree->Fill();
    }

    outFile->cd();
    outTree->Write();
    outFile->Close();

    f->Close();

    std::cout << "Escrito: " << outFilePath << std::endl;
}