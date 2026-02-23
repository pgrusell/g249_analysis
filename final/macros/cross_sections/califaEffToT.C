constexpr int kNCrystals = 2544;
constexpr int kPetalsBarrel = 1953;
constexpr int kPetalsEndcap = 2433;
constexpr int kOutputCrystals = 2455;
constexpr int kProbeOffset = kNCrystals; // PR partner = crId + 2544
constexpr int kEntriesDivisor = 100;

// Helpers

std::string repoPath()
{
    const char *env = getenv("repopath");
    if (!env)
        throw std::runtime_error("Environment variable 'repopath' is not set");
    return std::string(env);
}

std::vector<std::string> readFileList(const std::string &fileName = "25Ftest.txt")
{
    const std::string fullPath = repoPath() + "/final/settings/" + fileName;
    std::cout << "Reading file list: " << fullPath << std::endl;

    std::ifstream ifs(fullPath);
    if (!ifs.is_open())
        throw std::runtime_error("Cannot open file list: " + fullPath);

    std::vector<std::string> paths;
    std::string line;
    while (std::getline(ifs, line))
    {
        line.erase(std::remove(line.begin(), line.end(), '"'), line.end());
        if (!line.empty())
            paths.push_back(line);
    }
    return paths;
}

bool isValidAngle(const TVector3 &v)
{
    if (v.Mag() < 1e-6)
        return false;
    double theta = v.Theta() * TMath::RadToDeg();
    double phi = v.Phi();
    return std::isfinite(theta) && std::isfinite(phi);
}

double getThetaDeg(int crystalId)
{
    TVector3 angles = R3BCalifaGeometry::Instance()->GetAngles(crystalId);
    return angles.Theta() * TMath::RadToDeg();
}

void drawBoundaryLine(double x, double yMin, double yMax, Color_t color)
{
    auto *line = new TLine(x, yMin, x, yMax);
    line->SetLineColor(color);
    line->SetLineStyle(2);
    line->Draw();
}

// Data processing

struct EfficiencyData
{
    int nTotal = 0;
    int nRec = 0;
    std::vector<int> nTotalPerCrystal;
    std::vector<int> nRecPerCrystal;
    std::vector<double> thetaDeg;

    EfficiencyData() : nTotalPerCrystal(kNCrystals, 0),
                       nRecPerCrystal(kNCrystals, 0),
                       thetaDeg(kNCrystals, 0.0) {}

    double efficiency() const
    {
        return (nTotal > 0) ? 1.0 - static_cast<double>(nRec) / nTotal : 0.0;
    }
};

EfficiencyData computeEfficiency(TChain *evt)
{
    EfficiencyData data;

    // Precompute theta per crystal
    for (int i = 0; i < kNCrystals; i++)
        data.thetaDeg[i] = getThetaDeg(i + 1);

    auto *califaArray = new TClonesArray("R3BCalifaCrystalCalData");
    evt->SetBranchAddress("CalifaCrystalCalData", &califaArray);

    const Long64_t nEntries = evt->GetEntries() / kEntriesDivisor;

    for (Long64_t i = 0; i < nEntries; i++)
    {
        evt->GetEntry(i);
        const int nHits = califaArray->GetEntries();

        // Collect hit data for this event
        std::vector<int> crId(nHits);
        std::vector<double> calEnergy(nHits);
        std::vector<double> totEnergy(nHits);

        for (int j = 0; j < nHits; j++)
        {
            auto *hit = static_cast<R3BCalifaCrystalCalData *>(califaArray->At(j));
            crId[j] = hit->GetCrystalId();
            calEnergy[j] = hit->GetCalEnergy();
            totEnergy[j] = hit->GetToTEnergy();
        }

        // Build a set of crystal IDs for fast PR-partner lookup
        std::set<int> idSet(crId.begin(), crId.end());

        for (int j = 0; j < nHits; j++)
        {
            if (crId[j] > kNCrystals)
                continue; // only first half
            if (idSet.count(crId[j] + kProbeOffset) == 0)
                continue; // no PR partner
            if (calEnergy[j] > 0)
                continue; // calEnergy must be 0

            const int idx = crId[j] - 1;
            data.nTotal++;
            data.nTotalPerCrystal[idx]++;

            if (totEnergy[j] > 0)
                continue; // totEnergy must be 0
            data.nRec++;
            data.nRecPerCrystal[idx]++;
        }
    }

    std::cout << "Efficiency = " << data.efficiency() << std::endl;
    std::cout << "Entries processed: " << nEntries << std::endl;

    return data;
}

TEfficiency *makeEfficiency(const EfficiencyData &data)
{
    auto *hTotal = new TH1D("hTotal", "Total per crystal", kNCrystals, 0.5, kNCrystals + 0.5);
    auto *hPassed = new TH1D("hPassed", "Passed per crystal", kNCrystals, 0.5, kNCrystals + 0.5);

    for (int i = 0; i < kNCrystals; i++)
    {
        hTotal->SetBinContent(i + 1, data.nTotalPerCrystal[i]);
        hPassed->SetBinContent(i + 1, data.nTotalPerCrystal[i] - data.nRecPerCrystal[i]);
    }

    auto *eff = new TEfficiency(*hPassed, *hTotal);
    eff->SetTitle("ToT Efficiency per Crystal;Crystal ID;Efficiency");
    eff->SetMarkerStyle(20);
    eff->SetMarkerSize(0.5);
    return eff;
}

TGraph *makeThetaGraph()
{
    std::vector<double> gX, gY;
    for (int i = 0; i < kNCrystals; i++)
    {
        TVector3 angles = R3BCalifaGeometry::Instance()->GetAngles(i + 1);
        if (!isValidAngle(angles))
            continue;

        gX.push_back(i + 1);
        gY.push_back(angles.Theta() * TMath::RadToDeg());
    }

    auto *gr = new TGraph(gX.size(), gX.data(), gY.data());
    gr->SetTitle("#theta per crystal (geometry);Crystal ID;#theta (deg)");
    gr->SetMarkerStyle(20);
    gr->SetMarkerSize(0.3);
    gr->SetMarkerColor(kAzure + 2);
    return gr;
}

TH2D *makeCorrelationHist(const EfficiencyData &data, TEfficiency *eff)
{
    auto *h = new TH2D("hCorr", "ToT Efficiency vs #theta;#theta (deg);Efficiency",
                       180, 0., 180., 100, 0., 1.);
    h->SetStats(0);

    for (int i = 0; i < kNCrystals; i++)
    {
        if (data.nTotalPerCrystal[i] == 0)
            continue;

        TVector3 angles = R3BCalifaGeometry::Instance()->GetAngles(i + 1);
        if (!isValidAngle(angles))
            continue;

        double theta = angles.Theta() * TMath::RadToDeg();
        double effVal = eff->GetEfficiency(i + 1);
        if (!std::isfinite(effVal))
            continue;

        h->Fill(theta, effVal);
    }
    return h;
}

void drawCanvasEfficiency(TEfficiency *eff)
{
    auto *c = new TCanvas("c1", "ToT Efficiency per Crystal", 1200, 600);
    eff->Draw("AP");
    gPad->Update();

    if (auto *gr = eff->GetPaintedGraph())
    {
        gr->GetXaxis()->SetRangeUser(0.5, kNCrystals + 0.5);
        gr->GetYaxis()->SetRangeUser(0.7, 1.03);
        gPad->Update();
    }

    drawBoundaryLine(kPetalsBarrel, 0.7, 1.03, kRed);
    drawBoundaryLine(kPetalsEndcap, 0.7, 1.03, kBlue);
    c->Update();
}

void drawCanvasTheta(TGraph *gr)
{
    auto *c = new TCanvas("c2", "#theta per Crystal", 1200, 600);
    gr->Draw("AP");
    gPad->Update();

    gr->GetXaxis()->SetRangeUser(0.5, kNCrystals + 0.5);
    gr->GetYaxis()->SetRangeUser(0., 180.);
    gPad->Update();

    drawBoundaryLine(kPetalsBarrel, 0., 180., kRed);
    drawBoundaryLine(kPetalsEndcap, 0., 180., kBlue);
    c->Update();
}

TF1 *drawCanvasCorrelation(TH2D *hCorr)
{
    auto *c = new TCanvas("c3", "ToT Efficiency vs #theta", 800, 600);
    gStyle->SetPalette(kBird);
    hCorr->Draw("COLZ");
    gPad->Update();

    auto *prof = hCorr->ProfileX("profTotEff");
    auto *func = new TF1("funcEtotEff", "[0]*(1-[1]*exp(-[2]*x))");
    prof->Fit(func);
    prof->Draw();
    c->Update();

    return func;
}

void writeParameterFile(const EfficiencyData &data, TF1 *func)
{
    const int nEntries = kOutputCrystals * 2;
    const std::string outPath = repoPath() + "/final/macros/cross_sections/califaToTEffPar.txt";

    std::ofstream ofs(outPath);
    if (!ofs.is_open())
        throw std::runtime_error("Cannot open output file: " + outPath);

    ofs << "califaToTEffPar:  Float_t \\" << std::endl;

    for (int i = 0; i < nEntries; i++)
    {
        const int crId = i + 1;
        double val = 1.0;

        if (crId <= kPetalsBarrel)
        {
            double theta = data.thetaDeg[crId - 1];
            if (std::isfinite(theta) && theta > 1e-6)
                val = func->Eval(theta);
        }

        if (i % 10 == 0)
            ofs << "  ";

        ofs << Form("%.6f", val);

        bool isLineEnd = ((i + 1) % 10 == 0);
        bool isLast = (i == nEntries - 1);

        if (isLineEnd)
            ofs << (isLast ? "\n" : "  \\\n");
        else
            ofs << " ";
    }

    ofs.close();
    std::cout << "Written califaToTEffPar to: " << outPath << std::endl;
}

void califaEffToT()
{
    // Build TChain from file list
    std::vector<std::string> files = readFileList();
    for (const auto &f : files)
        std::cout << f << std::endl;

    auto *evt = new TChain("evt");
    for (const auto &f : files)
        evt->AddFile(f.c_str());

    // Init geometry
    R3BCalifaGeometry::Instance()->Init(2025);

    // Compute efficiency
    EfficiencyData data = computeEfficiency(evt);

    // Build plot objects
    TEfficiency *eff = makeEfficiency(data);
    TGraph *grTheta = makeThetaGraph();
    TH2D *hCorr = makeCorrelationHist(data, eff);

    // Draw canvases
    drawCanvasEfficiency(eff);
    drawCanvasTheta(grTheta);
    TF1 *fitFunc = drawCanvasCorrelation(hCorr);

    // Write output parameter file
    writeParameterFile(data, fitFunc);
}