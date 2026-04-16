void plot(TString fileName = "25F/sigt_d52-gs.txt")
{
    std::vector<double> binCenters;
    std::vector<double> binIntegral;

    std::ifstream infile(fileName.Data());
    if (!infile.is_open())
    {
        std::cerr << "Error opening file: " << fileName << std::endl;
        return;
    }

    std::string line;
    while (std::getline(infile, line))
    {
        std::istringstream iss(line);
        double binCent, binVal;

        if (!(iss >> binCent >> binVal))
            continue;

        binCenters.push_back(binCent);
        binIntegral.push_back(binVal);
    }

    const int nBins = binCenters.size();
    if (nBins < 2)
    {
        std::cerr << "Not enough bins read from file." << std::endl;
        return;
    }

    std::cout << "nBins = " << nBins << std::endl;

    const double binWidth = binCenters[1] - binCenters[0];
    const double maxBin = binCenters[nBins - 1] + binWidth / 2.;

    std::cout << "maxBin = " << maxBin << std::endl;

    auto *mom = new TH1F("mom", ";p [MeV/c];probability", nBins, -maxBin, maxBin);

    for (int i = 0; i < nBins; i++)
    {
        mom->SetBinContent(i + 1, binIntegral[i]);
    }

    mom->Draw();
}