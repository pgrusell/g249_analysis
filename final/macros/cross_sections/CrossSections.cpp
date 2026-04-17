#include "CrossSections.hpp"

CrossSections::CrossSections() {}
CrossSections::~CrossSections() {}

// ------
// Tree loading
// ------

TreeData CrossSections::LoadTreeData(TTree *tree, bool loadOpa)
{
    TreeData d;
    const Long64_t n = tree->GetEntries();
    d.Z_frag_est.resize(n);
    d.AoQ_frag.resize(n);
    if (loadOpa)
        d.califa_opa.resize(n);

    tree->SetBranchStatus("*", 0);
    tree->SetBranchStatus("Z_frag_est", 1);
    tree->SetBranchStatus("AoQ_frag", 1);

    double z = 0, aoq = 0, opa = 0;
    tree->SetBranchAddress("Z_frag_est", &z);
    tree->SetBranchAddress("AoQ_frag", &aoq);
    if (loadOpa)
    {
        tree->SetBranchStatus("califa_opa", 1);
        tree->SetBranchAddress("califa_opa", &opa);
    }

    for (Long64_t i = 0; i < n; ++i)
    {
        tree->GetEntry(i);
        d.Z_frag_est[i] = z;
        d.AoQ_frag[i] = aoq;
        if (loadOpa)
            d.califa_opa[i] = opa;
    }

    tree->SetBranchStatus("*", 1);
    tree->ResetBranchAddresses();
    return d;
}

void CrossSections::SetTrees(TString fileNameFragment, TString fileNameUnreacted,
                             TString treeName)
{
    TFile *fFrag = TFile::Open(fileNameFragment, "READ");
    if (!fFrag || fFrag->IsZombie())
    {
        std::cerr << "[ERROR] Cannot open fragment file: " << fileNameFragment << "\n";
        return;
    }
    TTree *tFrag = static_cast<TTree *>(fFrag->Get(treeName));
    if (!tFrag)
    {
        std::cerr << "[ERROR] Tree '" << treeName << "' not found in " << fileNameFragment << "\n";
        fFrag->Close();
        delete fFrag;
        return;
    }

    TFile *fUnr = TFile::Open(fileNameUnreacted, "READ");
    if (!fUnr || fUnr->IsZombie())
    {
        std::cerr << "[ERROR] Cannot open unreacted file: " << fileNameUnreacted << "\n";
        fFrag->Close();
        delete fFrag;
        return;
    }
    TTree *tUnr = static_cast<TTree *>(fUnr->Get(treeName));
    if (!tUnr)
    {
        std::cerr << "[ERROR] Tree '" << treeName << "' not found in " << fileNameUnreacted << "\n";
        fFrag->Close();
        delete fFrag;
        fUnr->Close();
        delete fUnr;
        return;
    }

    std::cout << "[INFO] Loading fragment (outgoing) tree...\n";
    fDataFragment = LoadTreeData(tFrag, /*loadOpa=*/true);
    std::cout << "[INFO] Loading unreacted (incoming) tree...\n";
    fDataUnreacted = LoadTreeData(tUnr, /*loadOpa=*/false);
    fTreesLoaded = true;

    std::cout << "[INFO] Trees loaded. fragment=" << fDataFragment.Z_frag_est.size()
              << " events, unreacted=" << fDataUnreacted.Z_frag_est.size() << " events.\n";
}

// -------
// 2D Gaussian fit (from TChain)
// -------

Fit2DParams CrossSections::Fit2D(TString label,
                                 double AoQ_min, double AoQ_max,
                                 double Z_min, double Z_max,
                                 std::vector<TString> files,
                                 TString treeName, int nBins)
{
    if (fFitCache.count(label))
        return fFitCache[label];

    static int callCount = 0;
    ++callCount;

    auto *chain = new TChain(treeName);
    for (auto &f : files)
        chain->Add(f);

    TString hname = Form("hFit2D_%d", callCount);
    auto *h = new TH2F(hname, "PID;AoQ;Z", nBins, 2, 3.3, nBins, 7.5, 10.5);

    auto *cTmp = new TCanvas(Form("cFit2D_%d", callCount));
    chain->Draw(Form("Z_frag_est:AoQ_frag>>%s", hname.Data()), "", "colz");

    TString fname = Form("fFit2D_%d", callCount);
    auto *f2 = new TF2(fname,
                       "[0]*exp(-0.5*((x-[1])/[2])^2 -0.5*((y-[3])/[4])^2)",
                       AoQ_min, AoQ_max, Z_min, Z_max);

    double mu_AoQ = 0.5 * (AoQ_min + AoQ_max);
    double mu_Z = 0.5 * (Z_min + Z_max);
    f2->SetParameters(h->GetMaximum(), mu_AoQ, 0.010, mu_Z, 0.15);
    f2->SetParLimits(0, 0, 1e9);
    f2->SetParLimits(1, AoQ_min, AoQ_max);
    f2->SetParLimits(2, 0.001, 0.05);
    f2->SetParLimits(3, Z_min, Z_max);
    f2->SetParLimits(4, 0.02, 0.5);

    h->Fit(f2, "IR");

    Fit2DParams p;
    p.amplitude = f2->GetParameter(0);
    p.muAoQ = f2->GetParameter(1);
    p.sigmaAoQ = f2->GetParameter(2);
    p.muZ = f2->GetParameter(3);
    p.sigmaZ = f2->GetParameter(4);
    p.errAmplitude = f2->GetParError(0);
    p.errMuAoQ = f2->GetParError(1);
    p.errSigmaAoQ = f2->GetParError(2);
    p.errMuZ = f2->GetParError(3);
    p.errSigmaZ = f2->GetParError(4);

    fFitCache[label] = p;

    if (fVerbose)
    {
        std::cout << "[Fit2D] label = " << label << "\n";
        std::cout << "  mu_AoQ   = " << p.muAoQ << " +/- " << p.errMuAoQ << "\n";
        std::cout << "  sigma_AoQ= " << p.sigmaAoQ << " +/- " << p.errSigmaAoQ << "\n";
        std::cout << "  mu_Z     = " << p.muZ << " +/- " << p.errMuZ << "\n";
        std::cout << "  sigma_Z  = " << p.sigmaZ << " +/- " << p.errSigmaZ << "\n";
    }

    delete cTmp;
    delete f2;
    delete chain;

    return p;
}

// -------
// 2D Gaussian fit (from existing TH2F)
// -------

Fit2DParams CrossSections::Fit2DFromHisto(TH2F *h,
                                          double AoQ_min, double AoQ_max,
                                          double Z_min, double Z_max,
                                          int uniqueId)
{
    Fit2DParams p;

    if (!h || h->GetEntries() == 0)
        return p;

    TString fname = Form("fFit2DH_%d", uniqueId);
    auto *f2 = new TF2(fname,
                       "[0]*exp(-0.5*((x-[1])/[2])^2 -0.5*((y-[3])/[4])^2)",
                       AoQ_min, AoQ_max, Z_min, Z_max);

    double mu_AoQ = 0.5 * (AoQ_min + AoQ_max);
    double mu_Z = 0.5 * (Z_min + Z_max);
    f2->SetParameters(h->GetMaximum(), mu_AoQ, 0.010, mu_Z, 0.15);
    f2->SetParLimits(0, 0, 1e9);
    f2->SetParLimits(1, AoQ_min, AoQ_max);
    f2->SetParLimits(2, 0.001, 0.05);
    f2->SetParLimits(3, Z_min, Z_max);
    f2->SetParLimits(4, 0.02, 0.5);

    h->Fit(f2, "QIR0"); // Q=quiet, I=use integral, R=range, 0=do not draw

    p.amplitude = f2->GetParameter(0);
    p.muAoQ = f2->GetParameter(1);
    p.sigmaAoQ = f2->GetParameter(2);
    p.muZ = f2->GetParameter(3);
    p.sigmaZ = f2->GetParameter(4);
    p.errAmplitude = f2->GetParError(0);
    p.errMuAoQ = f2->GetParError(1);
    p.errSigmaAoQ = f2->GetParError(2);
    p.errMuZ = f2->GetParError(3);
    p.errSigmaZ = f2->GetParError(4);

    delete f2;
    return p;
}

double CrossSections::RunFit2DForNucleus(TString nuc, double k)
{
    std::vector<TString> files;
    double aoqmin, aoqmax, zmin, zmax;
    int nBins;

    if (nuc == "25F")
    {
        files = {"/nucl_lustre/pablogrusell/g249/g249_analysis/results/final/data_25F.root"};
        aoqmin = 2.71;
        aoqmax = 2.77;
        zmin = 8.3;
        zmax = 9.5;
        nBins = 2500;
    }
    else if (nuc == "24O")
    {
        files = {"/nucl_lustre/pablogrusell/g249/g249_analysis/results/final/data_24O_test.root"};
        aoqmin = 2.91;
        aoqmax = 3.0;
        zmin = 7.6;
        zmax = 8.5;
        nBins = 250;
    }
    else if (nuc == "23O")
    {
        files = {"/nucl_lustre/pablogrusell/g249/g249_analysis/results/final/data_25Fp2p_test.root"};
        aoqmin = 2.78;
        aoqmax = 2.88;
        zmin = 7.6;
        zmax = 8.5;
        nBins = 250;
    }
    else
    {
        std::cerr << "[ERROR] Unknown nucleus: " << nuc << "\n";
        return 0.0;
    }

    Fit2DParams pars = Fit2D(nuc, aoqmin, aoqmax, zmin, zmax, files, "FilterDataTree", nBins);

    if (nuc == "24O")
    {
        fMuZ_Fragment = pars.muZ;
        fMuAoQ_Fragment = pars.muAoQ;
        fSigmaZ_Fragment = pars.sigmaZ;
        fSigmaAoQ_Fragment = pars.sigmaAoQ;
    }
    else if (nuc == "25F")
    {
        fMuZ_Unreacted = pars.muZ;
        fMuAoQ_Unreacted = pars.muAoQ;
        fSigmaZ_Unreacted = pars.sigmaZ;
        fSigmaAoQ_Unreacted = pars.sigmaAoQ;
    }

    double eff = 1.0 - std::exp(-0.5 * k * k);

    if (fVerbose)
    {
        std::cout << "[RunFit2DForNucleus] " << nuc << " k=" << k << "\n";
        std::cout << "  Efficiency = " << eff << "\n";
        std::cout << "  mu_AoQ     = " << pars.muAoQ << "\n";
        std::cout << "  sigma_AoQ  = " << pars.sigmaAoQ << "\n";
        std::cout << "  mu_Z       = " << pars.muZ << "\n";
        std::cout << "  sigma_Z    = " << pars.sigmaZ << "\n";
    }

    return eff;
}

// ----
// Fit cache I/O
// ----

void CrossSections::SaveFitCache(TString fileName) const
{
    std::ofstream out(fileName.Data());
    if (!out.is_open())
    {
        std::cerr << "[ERROR] Cannot open file for writing: " << fileName << "\n";
        return;
    }

    out << "# CrossSections fit cache\n";
    out << "# label  amplitude  muAoQ  sigmaAoQ  muZ  sigmaZ"
        << "  errAmplitude  errMuAoQ  errSigmaAoQ  errMuZ  errSigmaZ\n";

    for (const auto &entry : fFitCache)
    {
        const auto &p = entry.second;
        out << entry.first << "  "
            << p.amplitude << "  " << p.muAoQ << "  " << p.sigmaAoQ << "  "
            << p.muZ << "  " << p.sigmaZ << "  "
            << p.errAmplitude << "  " << p.errMuAoQ << "  " << p.errSigmaAoQ << "  "
            << p.errMuZ << "  " << p.errSigmaZ << "\n";
    }

    out.close();
    std::cout << "[INFO] Fit cache saved to " << fileName << " (" << fFitCache.size() << " entries)\n";
}

void CrossSections::LoadFitCache(TString fileName)
{
    std::ifstream in(fileName.Data());
    if (!in.is_open())
    {
        std::cerr << "[ERROR] Cannot open file for reading: " << fileName << "\n";
        return;
    }

    fFitCache.clear();
    std::string line;
    int count = 0;

    while (std::getline(in, line))
    {
        if (line.empty() || line[0] == '#')
            continue;

        std::istringstream iss(line);
        std::string label;
        Fit2DParams p;

        iss >> label >> p.amplitude >> p.muAoQ >> p.sigmaAoQ >> p.muZ >> p.sigmaZ >> p.errAmplitude >> p.errMuAoQ >> p.errSigmaAoQ >> p.errMuZ >> p.errSigmaZ;

        if (iss.fail())
        {
            std::cerr << "[WARN] Skipping malformed line: " << line << "\n";
            continue;
        }

        fFitCache[label.c_str()] = p;
        ++count;
    }

    in.close();
    std::cout << "[INFO] Fit cache loaded from " << fileName << " (" << count << " entries)\n";

    if (fFitCache.count("24O"))
    {
        const auto &p = fFitCache["24O"];
        fMuZ_Fragment = p.muZ;
        fMuAoQ_Fragment = p.muAoQ;
        fSigmaZ_Fragment = p.sigmaZ;
        fSigmaAoQ_Fragment = p.sigmaAoQ;
    }
    if (fFitCache.count("25F"))
    {
        const auto &p = fFitCache["25F"];
        fMuZ_Unreacted = p.muZ;
        fMuAoQ_Unreacted = p.muAoQ;
        fSigmaZ_Unreacted = p.sigmaZ;
        fSigmaAoQ_Unreacted = p.sigmaAoQ;
    }
}

const Fit2DParams &CrossSections::GetCachedFit(TString label) const
{
    auto it = fFitCache.find(label);
    if (it == fFitCache.end())
        throw std::runtime_error(Form("Fit2DParams not found for label '%s'", label.Data()));
    return it->second;
}

// ------
// OPA cut + signal extraction
// ------

OpaCutResults CrossSections::CalculateOpaLimits(const TreeData &data, double eff,
                                                double muZ, double muAoQ,
                                                double sigmaZ, double sigmaAoQ,
                                                double kEllipse, double kRoi)
{
    TH1F *h = BuildOpaHistogram(data, muZ, muAoQ, sigmaZ, sigmaAoQ, kEllipse, fOpaCallCount++);

    if (h->GetEntries() == 0)
    {
        if (fVerbose)
            std::cout << "[OPA] No entries after ellipse cut.\n";
        delete h;
        OpaCutResults res;
        res.kRoi = kRoi;
        return res;
    }

    OpaCutResults res = FitOpaHistogram(h, kRoi, false);
    return res;
}

double CrossSections::CalculateOutgoing(double eff, const OpaCutResults &opa)
{
    return opa.totalSignalFromRoi / eff;
}

IncomingFitResults CrossSections::CalculateIncoming(const TreeData &data,
                                                    double muZ, double muAoQ,
                                                    double sigmaZ, double sigmaAoQ,
                                                    double k, double purity)
{
    IncomingFitResults res;

    const double eff = 1.0 - std::exp(-0.5 * k * k);
    res.efficiency = eff;
    res.purity = purity;

    const double dZ2 = (k * sigmaZ) * (k * sigmaZ);
    const double dAoQ2 = (k * sigmaAoQ) * (k * sigmaAoQ);

    double counts = 0.0;
    const std::size_t n = data.Z_frag_est.size();
    for (std::size_t i = 0; i < n; ++i)
    {
        const double dz = data.Z_frag_est[i] - muZ;
        const double daoq = data.AoQ_frag[i] - muAoQ;
        if (dz * dz / dZ2 + daoq * daoq / dAoQ2 <= 1.0)
            counts += 1.0;
    }
    res.countsInRoi = counts;

    const double signalInRoi = counts * purity;
    res.yield = (eff > 0.0) ? signalInRoi / eff : 0.0;

    if (fVerbose)
    {
        std::cout << "\n=== INCOMING RESULTS (k=" << k << ") ===\n";
        std::cout << "Counts in ROI : " << counts << "\n";
        std::cout << "Purity        : " << purity << "\n";
        std::cout << "Signal in ROI : " << signalInRoi << "\n";
        std::cout << "Efficiency    : " << eff << "\n";
        std::cout << "Yield (total) : " << res.yield << "\n";
        std::cout << "==========================================\n\n";
    }

    return res;
}

// -----
// Cross-section computation
// -----

double CrossSections::ComputeCrossSectionWithParams(double eff25F, double eff24O,
                                                    double k1, double k2, double kOpa)
{
    auto opa = CalculateOpaLimits(fDataFragment, eff25F,
                                  fMuZ_Fragment, fMuAoQ_Fragment,
                                  fSigmaZ_Fragment, fSigmaAoQ_Fragment,
                                  k1, kOpa);
    double outgoing = CalculateOutgoing(eff24O, opa);

    double incoming = CalculateIncoming(fDataUnreacted,
                                        fMuZ_Unreacted, fMuAoQ_Unreacted,
                                        fSigmaZ_Unreacted, fSigmaAoQ_Unreacted,
                                        k2)
                          .yield;

    double nproj = incoming / (1. - fReactionRate);
    return -TMath::Log(1. - outgoing / nproj / fEffFrag / fEffToTeffp2p) * kMt / kRho / kZTarget / kNa * 1e24 * 1e3;
}

double CrossSections::ComputeCrossSection()
{
    if (!fTreesLoaded)
    {
        std::cerr << "[ERROR] Trees not loaded. Call SetTrees() first.\n";
        return 0.0;
    }

    double effFragment = 1.0 - std::exp(-0.5 * fKEllipseFragment * fKEllipseFragment);
    double effUnreacted = 1.0 - std::exp(-0.5 * fKEllipseUnreacted * fKEllipseUnreacted);

    double xs = ComputeCrossSectionWithParams(effFragment, effUnreacted,
                                              fKEllipseFragment, fKEllipseUnreacted, fKOpa);

    if (fVerbose)
        std::cout << "[RESULT] Cross section = " << xs << " mb\n";

    return xs;
}

void CrossSections::PlotNbOfYieldsVsK()
{
    if (!fTreesLoaded)
    {
        std::cerr << "[ERROR] Trees not loaded.\n";
        return;
    }

    const int nPoints = 50;
    const double kMin = 0.5, kMax = 4.0;
    const double kStep = (kMax - kMin) / (nPoints - 1);

    const std::vector<double> kOpaValues = {2.0, 2.5, 3.0, 3.5, 4.0};
    const std::vector<Color_t> colors = {kBlue, kRed, kGreen + 2, kMagenta, kOrange + 1};

    std::vector<std::vector<double>> kVals(kOpaValues.size(), std::vector<double>(nPoints));
    std::vector<std::vector<double>> yieldOut(kOpaValues.size(), std::vector<double>(nPoints));
    std::vector<std::vector<double>> yieldIn(kOpaValues.size(), std::vector<double>(nPoints));

    for (int iOpa = 0; iOpa < (int)kOpaValues.size(); ++iOpa)
    {
        const double kOpa = kOpaValues[iOpa];
        for (int i = 0; i < nPoints; ++i)
        {
            const double k = kMin + i * kStep;
            const double eff = 1.0 - std::exp(-0.5 * k * k);
            kVals[iOpa][i] = k;

            auto opa = CalculateOpaLimits(fDataFragment, eff,
                                          fMuZ_Fragment, fMuAoQ_Fragment,
                                          fSigmaZ_Fragment, fSigmaAoQ_Fragment,
                                          k, kOpa);
            yieldOut[iOpa][i] = CalculateOutgoing(eff, opa);

            if (iOpa == 0)
            {
                auto resIn = CalculateIncoming(fDataUnreacted,
                                               fMuZ_Unreacted, fMuAoQ_Unreacted,
                                               fSigmaZ_Unreacted, fSigmaAoQ_Unreacted,
                                               k);
                yieldIn[0][i] = resIn.yield;
            }
        }
    }

    auto *cIn = new TCanvas("cIncoming", "Incoming vs k", 900, 650);
    cIn->SetGrid();
    auto *grIn = new TGraph(nPoints, kVals[0].data(), yieldIn[0].data());
    grIn->SetTitle("Incoming (^{25}F) vs k;k (ellipse cut);Number of yields");
    grIn->SetMarkerStyle(20);
    grIn->SetMarkerSize(0.8);
    grIn->SetMarkerColor(kBlue);
    grIn->SetLineColor(kBlue);
    grIn->SetLineWidth(2);
    grIn->Draw("APL");
    auto *legIn = new TLegend(0.60, 0.75, 0.88, 0.88);
    legIn->SetBorderSize(1);
    legIn->AddEntry(grIn, "Incoming (^{25}F)", "lp");
    legIn->Draw();
    cIn->Update();

    auto *cOut = new TCanvas("cOutgoing", "Outgoing vs k", 900, 650);
    cOut->SetGrid();
    auto *mg = new TMultiGraph();
    auto *legOut = new TLegend(0.15, 0.60, 0.45, 0.88);
    legOut->SetBorderSize(1);

    for (int iOpa = 0; iOpa < (int)kOpaValues.size(); ++iOpa)
    {
        auto *gr = new TGraph(nPoints, kVals[iOpa].data(), yieldOut[iOpa].data());
        gr->SetMarkerStyle(20 + iOpa);
        gr->SetMarkerSize(0.8);
        gr->SetMarkerColor(colors[iOpa]);
        gr->SetLineColor(colors[iOpa]);
        gr->SetLineWidth(2);
        mg->Add(gr, "PL");
        legOut->AddEntry(gr, Form("k_{opa} = %.1f", kOpaValues[iOpa]), "lp");
    }

    mg->SetTitle("Outgoing (^{24}O) vs k;k (ellipse cut);Number of yields");
    mg->Draw("A");
    legOut->Draw();
    cOut->Update();
}

void CrossSections::PlotOutgoingVsKopa()
{
    if (!fTreesLoaded)
    {
        std::cerr << "[ERROR] Trees not loaded.\n";
        return;
    }

    const int nPoints = 30;
    const double kopaMin = 0.5, kopaMax = 4.5;
    const double kopaStep = (kopaMax - kopaMin) / (nPoints - 1);

    const double kFixed = 2.5;
    const double eff = 1.0 - std::exp(-0.5 * kFixed * kFixed);

    std::vector<double> kopaVals(nPoints), yOut(nPoints);

    for (int i = 0; i < nPoints; ++i)
    {
        const double kOpa = kopaMin + i * kopaStep;
        kopaVals[i] = kOpa;
        auto opa = CalculateOpaLimits(fDataFragment, eff,
                                      fMuZ_Fragment, fMuAoQ_Fragment,
                                      fSigmaZ_Fragment, fSigmaAoQ_Fragment,
                                      kFixed, kOpa);
        yOut[i] = CalculateOutgoing(eff, opa);
    }

    auto *c = new TCanvas("cOutVsKopa", "Outgoing vs k_{opa}", 900, 650);
    c->SetGrid();
    auto *gr = new TGraph(nPoints, kopaVals.data(), yOut.data());
    gr->SetTitle("Outgoing (^{24}O) vs k_{opa} (k_{PID}=2.5);k_{opa};Number of yields");
    gr->SetMarkerStyle(20);
    gr->SetMarkerSize(0.8);
    gr->SetMarkerColor(kRed);
    gr->SetLineColor(kRed);
    gr->SetLineWidth(2);
    gr->Draw("APL");
    auto *leg = new TLegend(0.60, 0.75, 0.88, 0.88);
    leg->SetBorderSize(1);
    leg->AddEntry(gr, "Outgoing, k_{PID}=2.5", "lp");
    leg->Draw();
    c->Update();
}

TH1F *CrossSections::BuildOpaHistogram(const TreeData &data,
                                       double muZ, double muAoQ,
                                       double sigmaZ, double sigmaAoQ,
                                       double kEllipse, int uniqueId) const
{
    const double dZ2 = (kEllipse * sigmaZ) * (kEllipse * sigmaZ);
    const double dAoQ2 = (kEllipse * sigmaAoQ) * (kEllipse * sigmaAoQ);

    TString hname = Form("hOpaBuilt_%d", uniqueId);
    auto *h = new TH1F(hname, "", 125, 0, 3);
    h->SetDirectory(nullptr);

    const std::size_t n = data.Z_frag_est.size();
    for (std::size_t i = 0; i < n; ++i)
    {
        const double dz = data.Z_frag_est[i] - muZ;
        const double daoq = data.AoQ_frag[i] - muAoQ;
        if (dz * dz / dZ2 + daoq * daoq / dAoQ2 <= 1.0)
            h->Fill(data.califa_opa[i]);
    }
    return h;
}

OpaCutResults CrossSections::FitOpaHistogram(TH1F *h, double kRoi, bool storeFit)
{
    ++fOpaCallCount;
    OpaCutResults res;
    res.kRoi = kRoi;

    if (h->GetEntries() == 0)
        return res;

    const double center = h->GetBinCenter(h->GetMaximumBin());

    TString fname = Form("funcOpaFit_%d", fOpaCallCount);
    auto *func = new TF1(fname,
                         "[0] + [1]*x + [2]*exp(-0.5*((x-[3])/[4])^2) + [5]*x*x",
                         0.4, 2.2);
    func->SetParameters(10.0, 0.0, 1000.0, center, 0.08, 0.0);
    func->SetParLimits(2, 0.0, 1e9);
    func->SetParLimits(3, 0.4, 2.5);
    func->SetParLimits(4, 0.01, 0.5);

    int fitStatus = h->Fit(func, storeFit ? "RQI" : "RQNI");

    res.chi2u = (double)func->GetChisquare() / (double)func->GetNDF();

    const double p0 = func->GetParameter(0);
    const double p1 = func->GetParameter(1);
    const double A = func->GetParameter(2);
    const double mu = func->GetParameter(3);
    const double sg = std::abs(func->GetParameter(4));
    const double p2 = func->GetParameter(5);

    res.mu = mu;
    res.sigma = sg;
    res.roiMin = mu - kRoi * sg;
    res.roiMax = mu + kRoi * sg;

    if (res.roiMin < h->GetXaxis()->GetXmin())
        res.roiMin = h->GetXaxis()->GetXmin();
    if (res.roiMax > h->GetXaxis()->GetXmax())
        res.roiMax = h->GetXaxis()->GetXmax();

    TString signame = Form("sigOpaFit_%d", fOpaCallCount);
    TString bkgname = Form("bkgOpaFit_%d", fOpaCallCount);

    auto *sig = new TF1(signame, "[0]*exp(-0.5*((x-[1])/[2])^2)", 0.0, 3.0);
    sig->SetParameters(A, mu, sg);

    auto *bkg = new TF1(bkgname, "[0] + [1]*x + [2]*x*x", 0.0, 3.0);
    bkg->SetParameters(p0, p1, p2);

    const double binWidth = h->GetXaxis()->GetBinWidth(1);

    res.signalInRoi = sig->Integral(res.roiMin, res.roiMax) / binWidth;
    res.backgroundInRoi = bkg->Integral(res.roiMin, res.roiMax) / binWidth;
    if (res.backgroundInRoi < 0.0)
        res.backgroundInRoi = 0.0;

    const int bin1 = h->FindBin(res.roiMin);
    const int bin2 = h->FindBin(res.roiMax);
    res.eventsRoi = h->Integral(bin1, bin2);
    res.signalTotal = A * sg * std::sqrt(2.0 * TMath::Pi()) / binWidth;

    res.efficiency = (res.signalTotal > 0.0) ? res.signalInRoi / res.signalTotal : 0.0;

    const double denom = res.signalInRoi + res.backgroundInRoi;
    res.purity = (denom > 0.0) ? res.signalInRoi / denom : 0.0;

    res.selectedSignal = res.eventsRoi * res.purity;
    res.totalSignalFromRoi = (res.efficiency > 0.0) ? res.selectedSignal / res.efficiency : 0.0;

    if (fVerbose)
    {
        std::cout << "\n=== OPA CUT RESULTS ===\n";
        std::cout << "Fit status: " << fitStatus << "\n";
        std::cout << "ROI: [" << res.roiMin << ", " << res.roiMax << "]\n";
        std::cout << "mu = " << res.mu << ", sigma = " << res.sigma << "\n";
        std::cout << "Events in ROI (data): " << res.eventsRoi << "\n";
        std::cout << "Signal in ROI (fit): " << res.signalInRoi << "\n";
        std::cout << "Background in ROI (fit): " << res.backgroundInRoi << "\n";
        std::cout << "Signal total from Gaussian fit: " << res.signalTotal << "\n";
        std::cout << "Efficiency (ROI / total Gaussian): " << res.efficiency << "\n";
        std::cout << "Purity (ROI): " << res.purity << "\n";
        std::cout << "Selected signal = " << res.selectedSignal << "\n";
        std::cout << "Total signal extrapolated = " << res.totalSignalFromRoi << "\n";
        std::cout << "=======================\n\n";
    }

    if (storeFit)
    {
        sig->SetLineColor(kGreen + 2);
        sig->SetLineStyle(2);
        bkg->SetLineColor(kMagenta);
        bkg->SetLineStyle(2);
        h->GetListOfFunctions()->Add(sig);
        h->GetListOfFunctions()->Add(bkg);
    }
    else
    {
        delete func;
        delete sig;
        delete bkg;
    }

    return res;
}

// ----
// Systematic uncertainty on the cross section
// ----

TH1F *CrossSections::ComputeCrossSectionSystematics(int N,
                                                    double kEllipseMin, double kEllipseMax,
                                                    double kOpaMin, double kOpaMax,
                                                    double effP2PMin, double effP2PMax)
{
    if (!fTreesLoaded)
    {
        std::cerr << "[ERROR] Trees not loaded. Call SetTrees() first.\n";
        return nullptr;
    }

    // Save the original value so we can restore it at the end
    const double origEffP2P = fEffToTeffp2p;
    const bool savedVerbose = fVerbose;
    fVerbose = false;

    auto *hXs = new TH1F("hXsSyst", ";#sigma [mb];counts", 500, 2.7, 3.3);
    hXs->SetDirectory(nullptr);

    // Also keep track of the sampled parameters vs xs for correlation plots
    auto *hKellVsXs = new TH2F("hKellVsXs", ";k_{ellipse};#sigma [mb]",
                               100, kEllipseMin, kEllipseMax, 200, 0, 10);
    hKellVsXs->SetDirectory(nullptr);

    auto *hKopaVsXs = new TH2F("hKopaVsXs", ";k_{opa};#sigma [mb]",
                               100, kOpaMin, kOpaMax, 200, 0, 10);
    hKopaVsXs->SetDirectory(nullptr);

    auto *hEffVsXs = new TH2F("hEffVsXs", ";#varepsilon_{ToT #times p2p};#sigma [mb]",
                              100, effP2PMin, effP2PMax, 200, 0, 10);
    hEffVsXs->SetDirectory(nullptr);

    std::cout << "[ComputeCrossSectionSystematics] Running " << N << " iterations...\n";

    for (int i = 0; i < N; ++i)
    {
        if (i % 1000 == 0)
            std::cout << "  iteration " << i << "/" << N << "\n";

        // 1) Sample parameters uniformly
        const double kEll1 = gRandom->Uniform(kEllipseMin, kEllipseMax);
        const double kEll2 = gRandom->Uniform(kEllipseMin, kEllipseMax);
        const double kOpa = gRandom->Uniform(kOpaMin, kOpaMax);
        const double effP2P = gRandom->Uniform(effP2PMin, effP2PMax);

        // 2) Derived efficiencies
        const double effEll1 = 1.0 - std::exp(-0.5 * kEll1 * kEll1);
        const double effEll2 = 1.0 - std::exp(-0.5 * kEll2 * kEll2);

        // 3) Temporarily set fEffToTeffp2p to the sampled value
        fEffToTeffp2p = effP2P;

        // 4) Outgoing: OPA cut with the sampled kEll (PID ellipse) and kOpa
        auto opa = CalculateOpaLimits(fDataFragment, effEll1,
                                      fMuZ_Fragment, fMuAoQ_Fragment,
                                      fSigmaZ_Fragment, fSigmaAoQ_Fragment,
                                      kEll1, kOpa);
        double outgoing = CalculateOutgoing(effEll1, opa);

        // 5) Incoming: ellipse count with the sampled kEll
        double incoming = CalculateIncoming(fDataUnreacted,
                                            fMuZ_Unreacted, fMuAoQ_Unreacted,
                                            fSigmaZ_Unreacted, fSigmaAoQ_Unreacted,
                                            kEll2)
                              .yield;

        // 6) Cross-section formula (same as ComputeCrossSectionWithParams)
        double nproj = incoming / (1. - fReactionRate);
        double xs = -TMath::Log(1. - outgoing / nproj / fEffFrag / fEffToTeffp2p) * kMt / kRho / kZTarget / kNa * 1e24 * 1e3;

        hXs->Fill(xs);
        hKellVsXs->Fill(kEll1, xs);
        hKopaVsXs->Fill(kOpa, xs);
        hEffVsXs->Fill(effP2P, xs);
    }

    // Restore original value
    fEffToTeffp2p = origEffP2P;
    fVerbose = savedVerbose;

    // 1D distribution
    auto *c1 = new TCanvas("cXsSyst", "Cross-section systematic distribution", 900, 650);
    hXs->SetLineColor(kBlue);
    hXs->SetLineWidth(2);
    hXs->Draw();

    double mean = hXs->GetMean();
    double rms = hXs->GetRMS();

    auto *leg = new TLegend(0.55, 0.72, 0.88, 0.88);
    leg->AddEntry(hXs, Form("N = %d", N), "l");
    leg->AddEntry((TObject *)nullptr, Form("Mean = %.3f mb", mean), "");
    leg->AddEntry((TObject *)nullptr, Form("RMS  = %.3f mb (%.1f%%)", rms, rms / mean * 100), "");
    leg->Draw();
    c1->Update();

    // 2D correlation plots
    auto *c2 = new TCanvas("cXsSystCorr", "XS vs sampled parameters", 1200, 400);
    c2->Divide(3, 1);

    c2->cd(1);
    hKellVsXs->Draw("colz");

    c2->cd(2);
    hKopaVsXs->Draw("colz");

    c2->cd(3);
    hEffVsXs->Draw("colz");

    c2->Update();

    std::cout << "\n[ComputeCrossSectionSystematics] Done.\n";
    std::cout << "  Mean  = " << mean << " mb\n";
    std::cout << "  RMS   = " << rms << " mb\n";
    std::cout << "  Rel.  = " << rms / mean * 100 << " %\n\n";

    return hXs;
}

// ---
// Analytical cross-section with error propagation
// ---

CrossSectionResult CrossSections::ComputeCrossSectionAnalytical(double dNp2p,
                                                                double dNproj,
                                                                double dEffP2P,
                                                                double dEffFrag)
{
    if (!fTreesLoaded)
    {
        std::cerr << "[ERROR] Trees not loaded. Call SetTrees() first.\n";
        return {};
    }

    CrossSectionResult res;

    const double effFragment = 1.0 - std::exp(-0.5 * fKEllipseFragment * fKEllipseFragment);

    auto opa = CalculateOpaLimits(fDataFragment, effFragment,
                                  fMuZ_Fragment, fMuAoQ_Fragment,
                                  fSigmaZ_Fragment, fSigmaAoQ_Fragment,
                                  fKEllipseFragment, fKOpa);
    res.Np2p = CalculateOutgoing(effFragment, opa);

    const double effUnreacted = 1.0 - std::exp(-0.5 * fKEllipseUnreacted * fKEllipseUnreacted);

    auto incoming = CalculateIncoming(fDataUnreacted,
                                      fMuZ_Unreacted, fMuAoQ_Unreacted,
                                      fSigmaZ_Unreacted, fSigmaAoQ_Unreacted,
                                      fKEllipseUnreacted);
    double yieldIncoming = incoming.yield;
    res.Nproj = yieldIncoming / (1.0 - fReactionRate);

    const double C = kMt / kRho / kZTarget / kNa * 1e24 * 1e3; // prefactor [mb]
    const double denom = res.Nproj * fEffFrag * fEffToTeffp2p;
    const double f = res.Np2p / denom; // ratio
    res.xs = -TMath::Log(1.0 - f) * C;

    //   Let g = Np2p / (Nproj * effFrag * effP2P)
    //   xs = -C * ln(1 - g)
    //   d(xs)/d(x_i) = C / (1 - g) * d(g)/d(x_i)
    //
    //   d(g)/d(Np2p)    =  1 / (Nproj * effFrag * effP2P)
    //   d(g)/d(Nproj)   = -Np2p / (Nproj^2 * effFrag * effP2P)
    //   d(g)/d(effP2P)  = -Np2p / (Nproj * effFrag * effP2P^2)
    //   d(g)/d(effFrag) = -Np2p / (Nproj * effFrag^2 * effP2P)

    const double omf = 1.0 - f;

    const double dxs_dNp2p = C / (omf * res.Nproj * fEffFrag * fEffToTeffp2p);
    const double dxs_dNproj = C * res.Np2p / (omf * res.Nproj * res.Nproj * fEffFrag * fEffToTeffp2p);
    const double dxs_dEff = C * res.Np2p / (omf * res.Nproj * fEffFrag * fEffToTeffp2p * fEffToTeffp2p);
    const double dxs_dEffFrag = C * res.Np2p / (omf * res.Nproj * fEffFrag * fEffFrag * fEffToTeffp2p);

    res.dNp2p = dNp2p;
    res.dNproj = dNproj;
    res.dEffP2P = dEffP2P;
    res.dEffFrag = dEffFrag;

    res.contNp2p = std::abs(dxs_dNp2p * dNp2p);
    res.contNproj = std::abs(dxs_dNproj * dNproj);
    res.contEff = std::abs(dxs_dEff * dEffP2P);
    res.contEffFrag = std::abs(dxs_dEffFrag * dEffFrag);

    res.dxs = std::sqrt(res.contNp2p * res.contNp2p +
                        res.contNproj * res.contNproj +
                        res.contEff * res.contEff +
                        res.contEffFrag * res.contEffFrag);

    std::cout << "  ANALYTICAL CROSS SECTION + ERROR PROPAGATION\n";

    std::cout << "  Nominal parameters:\n";
    std::cout << "    kEllipse (fragment)  = " << fKEllipseFragment << "\n";
    std::cout << "    kEllipse (unreacted) = " << fKEllipseUnreacted << "\n";
    std::cout << "    kOpa                 = " << fKOpa << "\n";
    std::cout << "    effToTeffp2p         = " << fEffToTeffp2p << "\n";
    std::cout << "    effFrag              = " << fEffFrag << "\n";
    std::cout << "    reactionRate         = " << fReactionRate << "\n\n";

    std::cout << "  Computed quantities:\n";
    std::cout << "    Np2p  (outgoing)     = " << res.Np2p << "\n";
    std::cout << "    Nproj (projectiles)  = " << res.Nproj << "\n";
    std::cout << "    ratio f              = " << f << "\n\n";

    std::cout << "  Cross section:\n";
    std::cout << "    sigma = " << res.xs << " +/- " << res.dxs << " mb\n";
    std::cout << "    rel.  = " << res.dxs / res.xs * 100 << " %\n\n";

    std::cout << "  Error budget:\n";
    std::cout << "    from Np2p     : " << res.contNp2p << " mb  ("
              << res.contNp2p / res.xs * 100 << " %)\n";
    std::cout << "    from Nproj    : " << res.contNproj << " mb  ("
              << res.contNproj / res.xs * 100 << " %)\n";
    std::cout << "    from effP2P   : " << res.contEff << " mb  ("
              << res.contEff / res.xs * 100 << " %)\n";
    std::cout << "    from effFrag  : " << res.contEffFrag << " mb  ("
              << res.contEffFrag / res.xs * 100 << " %)\n";

    return res;
}