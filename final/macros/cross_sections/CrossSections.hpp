#pragma once

// -----
// Data structures
// -----

struct TreeData
{
    std::vector<double> Z_frag_est;
    std::vector<double> AoQ_frag;
    std::vector<double> califa_opa; // empty if not loaded (incoming tree)
};

struct OpaCutResults
{
    double roiMin = 0.0;
    double roiMax = 0.0;

    double mu = 0.0;
    double sigma = 0.0;
    double kRoi = 0.0;

    double eventsRoi = 0.0;
    double signalInRoi = 0.0;
    double backgroundInRoi = 0.0;
    double signalTotal = 0.0;

    double efficiency = 0.0;
    double purity = 0.0;

    double selectedSignal = 0.0;
    double totalSignalFromRoi = 0.0;

    double chi2u = 0.;
};

struct IncomingFitResults
{
    double yield = 0.0;
    double countsInRoi = 0.0;
    double efficiency = 0.0;
    double purity = 0.0;
};

/// Parameters from the 2D Gaussian fit (5 params + 5 errors).
struct Fit2DParams
{
    double amplitude = 0.0;
    double muAoQ = 0.0;
    double sigmaAoQ = 0.0;
    double muZ = 0.0;
    double sigmaZ = 0.0;

    double errAmplitude = 0.0;
    double errMuAoQ = 0.0;
    double errSigmaAoQ = 0.0;
    double errMuZ = 0.0;
    double errSigmaZ = 0.0;

    std::vector<double> toVector() const
    {
        return {amplitude, muAoQ, sigmaAoQ, muZ, sigmaZ,
                errAmplitude, errMuAoQ, errSigmaAoQ, errMuZ, errSigmaZ};
    }

    static Fit2DParams fromVector(const std::vector<double> &v)
    {
        Fit2DParams p;
        if (v.size() < 10)
            return p;
        p.amplitude = v[0];
        p.muAoQ = v[1];
        p.sigmaAoQ = v[2];
        p.muZ = v[3];
        p.sigmaZ = v[4];
        p.errAmplitude = v[5];
        p.errMuAoQ = v[6];
        p.errSigmaAoQ = v[7];
        p.errMuZ = v[8];
        p.errSigmaZ = v[9];
        return p;
    }
};

/// Result of the analytical cross-section computation with propagated errors.
struct CrossSectionResult
{
    double xs = 0.0;  // central value [mb]
    double dxs = 0.0; // total uncertainty [mb]

    double Np2p = 0.0;  // outgoing yield (computed)
    double Nproj = 0.0; // number of projectiles (computed)

    double dNp2p = 0.0;    // input uncertainty on Np2p
    double dNproj = 0.0;   // input uncertainty on Nproj
    double dEffP2P = 0.0;  // input uncertainty on effToTeffp2p
    double dEffFrag = 0.0; // input uncertainty on effFrag

    double contNp2p = 0.0;    // contribution to dxs from Np2p  [mb]
    double contNproj = 0.0;   // contribution to dxs from Nproj [mb]
    double contEff = 0.0;     // contribution to dxs from eff   [mb]
    double contEffFrag = 0.0; // contribution to dxs from effFrag [mb]
};

// ====
// CrossSections class
// ====

class CrossSections
{
public:
    CrossSections();
    ~CrossSections();

    // --- Configuration setters

    void SetEfficiencyp2pAndTot(double val) { fEffToTeffp2p = val; }
    void SetReactionRate(double val) { fReactionRate = val; }
    void SetEffFrag(double val) { fEffFrag = val; }

    void SetKEllipseFragment(double val) { fKEllipseFragment = val; }
    void SetKEllipseUnreacted(double val) { fKEllipseUnreacted = val; }
    void SetKOpa(double val) { fKOpa = val; }

    void SetVerbose(bool v) { fVerbose = v; }

    // --- Tree loading

    void SetTrees(TString fileNameFragment, TString fileNameUnreacted,
                  TString treeName = "FilterDataTree");

    // --- 2D Gaussian fit

    Fit2DParams Fit2D(TString label,
                      double AoQ_min, double AoQ_max,
                      double Z_min, double Z_max,
                      std::vector<TString> files,
                      TString treeName = "FilterDataTree", int nBins = 2500);

    /// Fit a 2D Gaussian directly on an existing TH2F histogram.
    Fit2DParams Fit2DFromHisto(TH2F *h,
                               double AoQ_min, double AoQ_max,
                               double Z_min, double Z_max,
                               int uniqueId = 0);

    double RunFit2DForNucleus(TString nuc, double k);

    // --- Fit cache I/O

    void SaveFitCache(TString fileName) const;
    void LoadFitCache(TString fileName);
    bool HasCachedFit(TString label) const { return fFitCache.count(label) > 0; }
    const Fit2DParams &GetCachedFit(TString label) const;

    // ---- Physics computation

    double ComputeCrossSection();

    /// Compute the cross section at nominal k values and propagate
    /// uncertainties analytically from dNp2p, dNproj, dEffP2P and dEffFrag.
    /// The central values of Np2p and Nproj are computed internally
    /// from the loaded data using the configured k parameters.
    CrossSectionResult ComputeCrossSectionAnalytical(double dNp2p = 83.0,
                                                     double dNproj = 2200.0,
                                                     double dEffP2P = 0.0003,
                                                     double dEffFrag = 0.0);

    OpaCutResults CalculateOpaLimits(const TreeData &data, double eff,
                                     double muZ, double muAoQ,
                                     double sigmaZ, double sigmaAoQ,
                                     double kEllipse, double kRoi);

    double CalculateOutgoing(double eff, const OpaCutResults &opa);

    IncomingFitResults CalculateIncoming(const TreeData &data,
                                         double muZ, double muAoQ,
                                         double sigmaZ, double sigmaAoQ,
                                         double k, double purity = 1.0);

    // ── Systematic plots ─────────────────────────────────────────────────────

    void PlotNbOfYieldsVsK();
    void PlotOutgoingVsKopa();

    /// Compute the cross-section distribution from systematic variations.
    /// For each of N iterations the following are sampled uniformly:
    ///   - k_ellipse  in [kEllipseMin, kEllipseMax]  (same k for fragment & unreacted)
    ///   - k_opa      in [kOpaMin, kOpaMax]
    ///   - effToTeffp2p in [effP2PMin, effP2PMax]
    /// Returns the histogram of cross-section values (caller owns it).
    TH1F *ComputeCrossSectionSystematics(int N = 10000,
                                         double kEllipseMin = 2.0, double kEllipseMax = 3.0,
                                         double kOpaMin = 2.5, double kOpaMax = 3.5,
                                         double effP2PMin = 0.566313, double effP2PMax = 0.568788);

private:
    // --- Tree data
    TreeData fDataFragment;
    TreeData fDataUnreacted;
    bool fTreesLoaded = false;

    // --- Physical constants
    static constexpr double kMt = 1.000784;
    static constexpr double kRho = 0.0715;
    static constexpr double kZTarget = 5;
    static constexpr double kNa = 6.022E23;

    // --- Configurable parameters
    double fKEllipseFragment = 3.0;
    double fKEllipseUnreacted = 3.0;
    double fKOpa = 3.0;

    double fEffToTeffp2p = 0.566897;
    double fEffFrag = 0.9853;
    double fReactionRate = 0.0024;

    bool fVerbose = false;

    // --- Fit parameters
    double fMuZ_Fragment = 8.06109;
    double fMuAoQ_Fragment = 2.95752;
    double fSigmaZ_Fragment = 0.1;
    double fSigmaAoQ_Fragment = 0.01;

    double fMuZ_Unreacted = 8.97888;
    double fMuAoQ_Unreacted = 2.74036;
    double fSigmaZ_Unreacted = 0.202362;
    double fSigmaAoQ_Unreacted = 0.010049;

    // --- 2D fit cache
    std::map<TString, Fit2DParams> fFitCache;

    // --- Internal helpers
    TreeData LoadTreeData(TTree *tree, bool loadOpa);
    int fOpaCallCount = 0;

    OpaCutResults FitOpaHistogram(TH1F *h, double kRoi, bool storeFit = false);

    double ComputeCrossSectionWithParams(double eff25F, double eff24O,
                                         double k1, double k2, double kOpa);
    TH1F *BuildOpaHistogram(const TreeData &data,
                            double muZ, double muAoQ,
                            double sigmaZ, double sigmaAoQ,
                            double kEllipse, int uniqueId) const;
};