#pragma once

class DataAnalysis
{
public:
    explicit DataAnalysis(std::string dataFilePath, std::string outFileName, bool hasNeutrons);
    ~DataAnalysis();

    void setOffsetsFromTxt(std::string offFile);
    void getData(bool called = true);

    void matchBeta();

private:
    static constexpr double m_neut = 0.939565;                                 // MeV
    static constexpr double dalt = 931.494;                                    // MeV
    static constexpr double Mproj_GeV = (25.0 * 0.93149410242 - 9 * 0.000511); // TODO: read this from the eventFilter tree

    std::string fDataFilePath = "";
    std::vector<double> fNeutronPOffsets{0.0, 0.0};
    std::vector<double> fFragmentPOffsets{0.0, 0.0};
    double fBetaMatchValue = 0.;
    bool fHasNeutrons = false;

    std::string fOutFileName = "";

    TH1F *fErel = nullptr;
    TH1F *fDeltaBeta = nullptr;
};