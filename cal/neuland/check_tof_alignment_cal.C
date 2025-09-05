void check_tof_alignment_cal(
    TString absPath = "",
    std::vector<TString> inFiles = {
        "g249_all_det_offline_20250828_210751.root",
        "g249_all_det_offline_20250828_210751_1.root",
        "g249_all_det_offline_20250828_210751_2.root"},
    TString fileOutName = "tof_alignment_losNeuland.root", TString calFile = "tfine_run112to117_2.txt")
{
    // Paths
    TString repopath = gSystem->Getenv("repopath");
    if (absPath.IsNull() || absPath == "")
        absPath = repopath + "/data/";
    if (absPath.IsNull() || absPath == "")
    {
        std::cerr << "[WARN] $repopath not defined, using '.'\n";
        absPath = ".";
    }

    // Open cal file and store the pedestals
    int iPaddle;
    double iPedestal;
    std::vector<double> pedestals(1300);

    std::ifstream fin(repopath + "/data/" + calFile);
    while (fin >> iPaddle >> iPedestal)
    {
        std::cout << iPaddle << " " << iPedestal << "\n";
        pedestals[iPaddle - 1] = iPedestal;
    }

    // Chain all files
    TChain chain("evt");
    for (const auto &f : inFiles)
    {
        TString full = absPath + "/" + f;
        if (gSystem->AccessPathName(full))
        {
            std::cerr << "[WARN] not found: " << full << "\n";
            continue;
        }
        chain.Add(full);
    }
    if (chain.GetEntries() == 0)
    {
        std::cerr << "[ERROR] evt is empty.\n";
        return;
    }

    // Reader
    TTreeReader reader(&chain);

    // NEULAND
    TTreeReaderArray<int> neuPaddle(reader, "NeulandHits.fPaddle");
    TTreeReaderArray<double> neuT(reader, "NeulandHits.fT");            // ns
    TTreeReaderArray<TVector3> neuPos(reader, "NeulandHits.fPosition"); // cm

    // LOS
    TTreeReaderArray<double> losZ(reader, "LosHit.fZ"); // double

    const double Lref = 1557.0; // cm
    const double c = 29.979;    // cm/ns

    TH2D hTofVsPaddles("hTofVsPaddles", "ToF vs Paddle;Paddle;ToF (ns)",
                       1300, 0.0, 1300.0,
                       1000, 46.0, 55.0);

    Long64_t iev = 0;
    while (reader.Next())
    {
        if ((iev++ % 100000) == 0)
            std::cout << iev << '\n';

        if (losZ.GetSize() != 1)
            continue;
        int nNeu = std::min({neuPaddle.GetSize(), neuT.GetSize(), neuPos.GetSize()});
        if (nNeu <= 0)
            continue;

        for (int i = 0; i < nNeu; ++i)
        {
            const double tof_corr = neuT[i] - (neuPos[i].Mag() - Lref) / c;
            const double tof_diff = 51.98 - pedestals[neuPaddle[i] - 1];
            hTofVsPaddles.Fill(neuPaddle[i], tof_corr + tof_diff);
        }
    }

    // Save
    TString outDir = repopath + "/results/cal";
    gSystem->mkdir(outDir, kTRUE);

    TFile fout(outDir + "/" + fileOutName, "RECREATE");
    hTofVsPaddles.Write();
    fout.Close();
}