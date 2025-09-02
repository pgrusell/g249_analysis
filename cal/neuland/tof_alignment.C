#include <ROOT/RDataFrame.hxx>
#include <ROOT/RVec.hxx>
#include <TH2D.h>
#include <TFile.h>
#include <TSystem.h>
#include <TString.h>
#include <TVector3.h>
#include <vector>
#include <iostream>

using ROOT::VecOps::RVec;

// Function to act on every event and hit
RVec<double> ToFCorr(const RVec<double> &t, const RVec<TVector3> &pos)
{
    const double Lref = 1557.0;
    const double c = 29.979;
    const size_t n = std::min(t.size(), pos.size());
    RVec<double> out;
    out.reserve(n);
    for (size_t i = 0; i < n; ++i)
        out.push_back(t[i] - (pos[i].Mag() - Lref) / c);
    return out;
}

void tof_alignment(
    TString absPath = "",
    std::vector<TString> inFiles = {
        "g249_all_det_offline_20250828_210751.root",
        "g249_all_det_offline_20250828_210751_1.root",
        "g249_all_det_offline_20250828_210751_2.root"},
    TString fileOutName = "tof_alignment_losNeuland.root")
{
    // Base path
    TString repopath = gSystem->Getenv("repopath");
    if (absPath.IsNull() || absPath == "")
        absPath = repopath + "/data";
    if (absPath.IsNull() || absPath == "")
    {
        std::cerr << "[WARN] $repopath no definido, uso '.'\n";
        absPath = ".";
    }

    // Files
    TString inDir = absPath;
    std::vector<std::string> files;
    files.reserve(inFiles.size());
    for (const auto &f : inFiles)
        files.emplace_back((inDir + '/' + f).Data());

    // RDF
    ROOT::EnableImplicitMT();
    ROOT::RDataFrame df("evt", files);

    // One point per neuland hit
    auto df2 = df
                   .Define("nHitsNeu", "static_cast<int>(NeulandHits.GetEntriesFast())")
                   .Define("nHitsLos", "static_cast<int>(LosHit.GetEntriesFast())")
                   .Filter("nHitsNeu > 0 && nHitsLos == 1")
                   .Define("paddle", "NeulandHits.fPaddle")
                   .Define("tof", ToFCorr, {"NeulandHits.fT", "NeulandHits.fPosition"});

    // Histo
    auto hTofVsPaddles = df2.Histo2D(
        {"hTofVsPaddles", "ToF vs Paddle;Paddle;ToF (ns)",
         1300, 0.0, 1300.0,
         1000, 46.0, 55.0},
        "paddle", "tof");

    // save
    TString outDir = repopath + "/results/cal";
    TFile out(outDir + "/" + fileOutName, "RECREATE");
    hTofVsPaddles->Write();
    out.Close();
    std::cout << "[OK] Escrito: " << (outDir + "/" + fileOutName) << "\n";
}