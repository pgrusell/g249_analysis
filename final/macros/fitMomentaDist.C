// fitMomentaDist.C
// ROOT macro: load many theory momentum distributions (.txt), resample them to a common uniform binning,
// and fit an experimental momentum histogram as a linear combination of N chosen theory components.
//
// Usage (in ROOT):
//   .L fitMomentaDist.C+
//   fitMomentaDist(50, -400, 400, "Qy", {"1d52_Sp25p509"}); // 1 component
//   fitMomentaDist(40, -400, 400, "Qy", {"1d52_Sp30p930","2s12_Sp30p930","1p12_Sp30p930"}); // 3 comp
//
// Notes:
// - Units: if exp variable is in GeV/c, multiply by 1000 in the TTree::Draw expression.
// - Model: counts(x) = sum_i a_i * shape_i(x), with a_i >= 0, where shape_i are area-normalized.
// - You can choose any subset of components to fit via the vector<string> comps.
//

// ------------------------------------------------------------
// Your helpers (kept, slightly hardened)
// ------------------------------------------------------------
TTree *GetTreeFromFile(TFile *f)
{
    if (!f || f->IsZombie())
        return nullptr;

    TTree *t = nullptr;

    f->GetObject("KinTree", t);
    if (t)
        return t;

    TIter next(f->GetListOfKeys());
    TKey *k;
    while ((k = (TKey *)next()))
    {
        if (std::string(k->GetClassName()) == "TTree")
            return (TTree *)k->ReadObj();
    }
    return nullptr;
}

void Normalize(TH1D *h)
{
    if (!h)
        return;
    double I = h->Integral("width"); // more correct when comparing shapes of differential distributions
    if (I > 0)
        h->Scale(1.0 / I);
}

struct MomDistsHists
{
    TH1D *hQz = nullptr;
    TH1D *hQt = nullptr;
    TH1D *hQy = nullptr;
    TH1D *hQ = nullptr;
};

static std::vector<double> BuildEdgesFromCenters(const std::vector<double> &centers)
{
    if (centers.size() < 2)
        throw std::runtime_error("Need at least 2 Qi points to build bin edges.");

    std::vector<double> c = centers;
    std::sort(c.begin(), c.end());

    const size_t n = c.size();
    std::vector<double> edges(n + 1, 0.0);

    for (size_t i = 1; i < n; ++i)
        edges[i] = 0.5 * (c[i - 1] + c[i]);

    const double firstStep = c[1] - c[0];
    const double lastStep = c[n - 1] - c[n - 2];

    edges[0] = c[0] - 0.5 * firstStep;
    edges[n] = c[n - 1] + 0.5 * lastStep;

    return edges;
}

// Reads txt (Qi, dS/dQz, dS/dQt, dS/dQy, dS/dQ)
MomDistsHists ReadMomentumDistributionsTxt(const std::string &filepath,
                                           const std::string &nameTag = "theory",
                                           const std::string &titleTag = "theory")
{
    std::ifstream in(filepath);
    if (!in.is_open())
        throw std::runtime_error("Cannot open file: " + filepath);

    std::vector<double> Qi, vQz, vQt, vQy, vQ;

    std::string line;
    while (std::getline(in, line))
    {
        if (line.empty())
            continue;

        bool hasAlpha = false;
        for (char ch : line)
        {
            if ((ch >= 'a' && ch <= 'z') || (ch >= 'A' && ch <= 'Z') || ch == '(' || ch == ')')
            {
                if (ch != 'E' && ch != 'e')
                {
                    hasAlpha = true;
                    break;
                }
            }
        }
        if (hasAlpha)
            continue;

        std::istringstream ss(line);

        double x, a, b, c, d;
        if (!(ss >> x >> a >> b >> c >> d))
            continue;

        Qi.push_back(x);
        vQz.push_back(a);
        vQt.push_back(b);
        vQy.push_back(c);
        vQ.push_back(d);
    }

    if (Qi.size() < 2)
        throw std::runtime_error("Parsed <2 data rows from: " + filepath);

    std::vector<size_t> idx(Qi.size());
    for (size_t i = 0; i < idx.size(); ++i)
        idx[i] = i;
    std::sort(idx.begin(), idx.end(), [&](size_t i, size_t j)
              { return Qi[i] < Qi[j]; });

    std::vector<double> QiS, QzS, QtS, QyS, QS;
    QiS.reserve(Qi.size());
    QzS.reserve(Qi.size());
    QtS.reserve(Qi.size());
    QyS.reserve(Qi.size());
    QS.reserve(Qi.size());

    for (auto i : idx)
    {
        QiS.push_back(Qi[i]);
        QzS.push_back(vQz[i]);
        QtS.push_back(vQt[i]);
        QyS.push_back(vQy[i]);
        QS.push_back(vQ[i]);
    }

    std::vector<double> edges = BuildEdgesFromCenters(QiS);
    const int nbins = static_cast<int>(QiS.size());

    MomDistsHists out;

    out.hQz = new TH1D(Form("hQz_%s", nameTag.c_str()),
                       Form("dS/dQz (%s);Q_{z} (MeV/c);dS/dQ_{z} (mb/(MeV/c))", titleTag.c_str()),
                       nbins, edges.data());

    out.hQt = new TH1D(Form("hQt_%s", nameTag.c_str()),
                       Form("dS/dQt (%s);Q_{t} (MeV/c);dS/dQ_{t} (mb/(MeV/c))", titleTag.c_str()),
                       nbins, edges.data());

    out.hQy = new TH1D(Form("hQy_%s", nameTag.c_str()),
                       Form("dS/dQy (%s);Q_{y} (MeV/c);dS/dQ_{y} (mb/(MeV/c))", titleTag.c_str()),
                       nbins, edges.data());

    out.hQ = new TH1D(Form("hQ_%s", nameTag.c_str()),
                      Form("dS/dQ (%s);Q (MeV/c);dS/dQ (mb/(MeV/c))", titleTag.c_str()),
                      nbins, edges.data());

    for (int i = 0; i < nbins; ++i)
    {
        const int bin = i + 1;
        out.hQz->SetBinContent(bin, QzS[i]);
        out.hQt->SetBinContent(bin, QtS[i]);
        out.hQy->SetBinContent(bin, QyS[i]);
        out.hQ->SetBinContent(bin, QS[i]);
    }

    out.hQz->SetDirectory(nullptr);
    out.hQt->SetDirectory(nullptr);
    out.hQy->SetDirectory(nullptr);
    out.hQ->SetDirectory(nullptr);

    return out;
}

// --- Resample theory to uniform binning by averaging within each bin ---
static TH1D *ResampleTheoryToUniformBinning(const TH1D *hThIn,
                                            const char *name,
                                            int nbins, double xmin, double xmax,
                                            int nSamplePerBin = 30)
{
    if (!hThIn)
        return nullptr;

    TH1D *h = new TH1D(name, hThIn->GetTitle(), nbins, xmin, xmax);
    h->Sumw2(false);

    for (int ib = 1; ib <= nbins; ++ib)
    {
        const double lo = h->GetXaxis()->GetBinLowEdge(ib);
        const double hi = h->GetXaxis()->GetBinUpEdge(ib);

        double sum = 0.0;
        for (int k = 0; k < nSamplePerBin; ++k)
        {
            const double x = lo + (k + 0.5) * (hi - lo) / nSamplePerBin;
            sum += hThIn->Interpolate(x);
        }
        const double avg = sum / nSamplePerBin;

        h->SetBinContent(ib, avg);
        h->SetBinError(ib, 0.0);
    }

    h->SetDirectory(nullptr);
    return h;
}

// ------------------------------------------------------------
// Multi-component machinery
// ------------------------------------------------------------
static TH1D *GetComponentHisto(const MomDistsHists &m, const std::string &which)
{
    if (which == "Qz")
        return m.hQz;
    if (which == "Qt")
        return m.hQt;
    if (which == "Qy")
        return m.hQy;
    if (which == "Q")
        return m.hQ;
    return nullptr;
}

// key: short label, value: filepath
static std::vector<std::pair<std::string, std::string>> GetTheoryFileList(const std::string &theoryDir)
{
    return {
        {"1d52_Sp25p509", theoryDir + "/25F(p,2p)-1d52-Sp-25.509.txt"},

        {"1d52_Sp29p850", theoryDir + "/25F(p,2p)-1d52-Sp-29.850.txt"},
        {"2s12_Sp29p850", theoryDir + "/25F(p,2p)-2s12-Sp-29.850.txt"},

        {"1d52_Sp30p930", theoryDir + "/25F(p,2p)-1d52-Sp-30.930.txt"},
        {"2s12_Sp30p930", theoryDir + "/25F(p,2p)-2s12-Sp-30.930.txt"},
        {"1p12_Sp30p930", theoryDir + "/25F(p,2p)-1p12-Sp-30.930.txt"},

        {"1d52_Sp33p000", theoryDir + "/25F(p,2p)-1d52-Sp-33.000.txt"},
        {"2s12_Sp33p000", theoryDir + "/25F(p,2p)-2s12-Sp-33.000.txt"},
        {"1p12_Sp33p000", theoryDir + "/25F(p,2p)-1p12-Sp-33.000.txt"},
        {"1p32_Sp33p000", theoryDir + "/25F(p,2p)-1p32-Sp-33.000.txt"},

        {"1d52_Sp35p392", theoryDir + "/25F(p,2p)-1d52-Sp-35.392.txt"},
        {"2s12_Sp35p392", theoryDir + "/25F(p,2p)-2s12-Sp-35.392.txt"},
        {"1p12_Sp35p392", theoryDir + "/25F(p,2p)-1p12-Sp-35.392.txt"},
        {"1p32_Sp35p392", theoryDir + "/25F(p,2p)-1p32-Sp-35.392.txt"},
    };
}

// Loads all txts, keeps original binning
static std::map<std::string, MomDistsHists> LoadAllTheoryRaw(const std::string &theoryDir)
{
    std::map<std::string, MomDistsHists> out;
    for (const auto &kv : GetTheoryFileList(theoryDir))
    {
        const std::string &key = kv.first;
        const std::string &path = kv.second;
        out[key] = ReadMomentumDistributionsTxt(path, key, key);
    }
    return out;
}

// Resample all theory components to common uniform binning and normalize to shape=1
static std::map<std::string, TH1D *> LoadAllTheoryResampledShapes(const std::string &theoryDir,
                                                                  const std::string &which,
                                                                  int nbins, double xmin, double xmax,
                                                                  int nSamplePerBin = 30)
{
    auto raw = LoadAllTheoryRaw(theoryDir);

    std::map<std::string, TH1D *> out;
    for (auto &kv : raw)
    {
        const std::string &key = kv.first;
        MomDistsHists &m = kv.second;

        TH1D *hComp = GetComponentHisto(m, which);
        if (!hComp)
            throw std::runtime_error("Invalid 'which' in LoadAllTheoryResampledShapes: " + which);

        TH1D *h = ResampleTheoryToUniformBinning(hComp, Form("h_%s_%s", which.c_str(), key.c_str()),
                                                 nbins, xmin, xmax, nSamplePerBin);
        Normalize(h);
        out[key] = h;
    }
    return out;
}

// Fit exp with N components: model(x) = sum_i a_i * shape_i(x), a_i >= 0
static TF1 *FitMultiComponent(TH1D *hExp,
                              const std::vector<std::string> &comps,
                              const std::map<std::string, TH1D *> &shapes,
                              double xmin, double xmax,
                              const char *fname = "fMulti")
{
    if (!hExp)
        return nullptr;
    if (comps.empty())
        throw std::runtime_error("FitMultiComponent: empty comps");

    std::vector<const TH1D *> templ;
    templ.reserve(comps.size());
    for (const auto &k : comps)
    {
        auto it = shapes.find(k);
        if (it == shapes.end() || !it->second)
            throw std::runtime_error("FitMultiComponent: missing shape for component: " + k);
        templ.push_back(it->second);
    }

    const int npar = (int)templ.size();

    TF1 *f = new TF1(fname, [templ](double *x, double *p)
                     {
                         double y = 0.0;
                         for (int i = 0; i < (int)templ.size(); ++i)
                             y += p[i] * templ[i]->Interpolate(x[0]);
                         return y; }, xmin, xmax, npar);

    for (int i = 0; i < npar; ++i)
    {
        f->SetParName(i, Form("a%d", i));
        f->SetParameter(i, hExp->Integral() / npar);
        f->SetParLimits(i, 0.0, 1e18);
    }

    hExp->Fit(f, "R0"); // silent + range only
    return f;
}

// Build model histogram from fit params (in same binning as hExp)
static TH1D *BuildModelFromFit(const char *name,
                               const std::vector<std::string> &comps,
                               const std::map<std::string, TH1D *> &shapes,
                               const TF1 *f)
{
    if (comps.empty())
        return nullptr;
    auto it0 = shapes.find(comps[0]);
    if (it0 == shapes.end())
        return nullptr;

    TH1D *hModel = (TH1D *)it0->second->Clone(name);
    hModel->SetDirectory(nullptr);
    hModel->Reset("ICES");

    for (size_t i = 0; i < comps.size(); ++i)
    {
        auto it = shapes.find(comps[i]);
        if (it == shapes.end() || !it->second)
            continue;
        const double a = f ? f->GetParameter((int)i) : 0.0;
        hModel->Add(it->second, a);
    }
    return hModel;
}

// ------------------------------------------------------------
// Main macro
// ------------------------------------------------------------
// You choose:
//  - nbinsResample: binning for BOTH exp & theory in [-400,400] (or user range)
//  - whichTheory: "Qy" / "Qz" / "Qt" / "Q"
//  - compsToFit: subset of keys (see list below) with any size N
//
// Available keys:
//  1d52_Sp25p509
//  1d52_Sp29p850, 2s12_Sp29p850
//  1d52_Sp30p930, 2s12_Sp30p930, 1p12_Sp30p930
//  1d52_Sp33p000, 2s12_Sp33p000, 1p12_Sp33p000, 1p32_Sp33p000
//  1d52_Sp35p392, 2s12_Sp35p392, 1p12_Sp35p392, 1p32_Sp35p392
//
void fitMomentaDist(int nbinsResample = 50,
                    double xmin = -400.0,
                    double xmax = +400.0,
                    std::string whichTheory = "Qy",
                    std::vector<std::string> compsToFit = {"1d52_Sp25p509"},
                    std::string expVarExpr = "(1000.0*px_rf_rot -2.7798553)", // change to your variable; keep *1000 if GeV/c
                    TString cond = "califa_opa > 1.2 && califa_opa < 1.4")
{
    gStyle->SetOptStat(0);

    std::string repopath = getenv("repopath") ? getenv("repopath") : ".";

    std::string file23 = repopath + "/results/final/23O_analyzed_test.root";
    std::string file24 = repopath + "/results/final/24O_analyzed_test.root";

    TFile *f23 = TFile::Open(file23.c_str(), "READ");
    TFile *f24 = TFile::Open(file24.c_str(), "READ");
    if (!f23 || f23->IsZombie() || !f24 || f24->IsZombie())
    {
        std::cerr << "[fitMomentaDist] Error loading files:\n"
                  << "  " << file23 << "\n"
                  << "  " << file24 << "\n";
        return;
    }

    TTree *t23 = GetTreeFromFile(f23);
    TTree *t24 = GetTreeFromFile(f24);
    if (!t23 || !t24)
    {
        std::cerr << "[fitMomentaDist] Error: missing TTrees (expected KinTree)\n";
        return;
    }

    // ---- Load & resample ALL theory shapes to the binning YOU chose ----
    std::string theoryDir = "../../theory";
    std::map<std::string, TH1D *> theoryShapes;
    try
    {
        theoryShapes = LoadAllTheoryResampledShapes(theoryDir, whichTheory, nbinsResample, xmin, xmax, 40);
    }
    catch (const std::exception &e)
    {
        std::cerr << "[fitMomentaDist] Theory loading failed: " << e.what() << "\n";
        return;
    }

    // ---- Build EXP histogram with EXACT SAME binning ----
    TH1D *hExp = new TH1D("hExp", Form("^{24}O: %s;momentum (MeV/c);Counts", expVarExpr.c_str()),
                          nbinsResample, xmin, xmax);
    hExp->Sumw2();

    // Fill
    t24->Draw(Form("%s>>hExp", expVarExpr.c_str()), cond, "goff");

    // quick sanity print
    std::cout << "[fitMomentaDist] Exp entries=" << hExp->GetEntries()
              << "  integral=" << hExp->Integral() << "\n";

    // ---- Fit N components ----
    TF1 *f = nullptr;
    try
    {
        f = FitMultiComponent(hExp, compsToFit, theoryShapes, xmin, xmax, "fMulti");
    }
    catch (const std::exception &e)
    {
        std::cerr << "[fitMomentaDist] Fit failed: " << e.what() << "\n";
        return;
    }

    // ---- Build model histogram from fitted coefficients ----
    TH1D *hModel = BuildModelFromFit("hModel", compsToFit, theoryShapes, f);
    if (hModel)
    {
        hModel->SetLineColor(kRed);
        hModel->SetLineWidth(2);
    }

    // ---- Draw ----
    TCanvas *c = new TCanvas("c_fit", "Multi-component fit", 950, 720);
    hExp->SetMarkerStyle(20);
    hExp->Draw("E1");
    if (hModel)
        hModel->Draw("HIST SAME");

    // legend + print coefficients
    TLegend *leg = new TLegend(0.60, 0.65, 0.88, 0.88);
    leg->SetBorderSize(0);
    leg->AddEntry(hExp, "Exp", "pe");
    if (hModel)
        leg->AddEntry(hModel, "Fit model", "l");
    leg->Draw();

    if (f)
    {
        std::cout << "\n[fitMomentaDist] chi2/ndf = " << f->GetChisquare() << "/" << f->GetNDF() << "\n";
        for (size_t i = 0; i < compsToFit.size(); ++i)
        {
            std::cout << "  a[" << i << "]  " << compsToFit[i] << " = "
                      << f->GetParameter((int)i) << " +- " << f->GetParError((int)i) << "\n";
        }
        std::cout << std::endl;
    }

    // Note: theoryShapes histograms are heap-allocated; ROOT will keep them alive for the session.
    // If you want to clean explicitly, delete them here (after you're done using them).
}