#include <ROOT/RDataFrame.hxx>
#include <ROOT/RVec.hxx>
#include <TH2D.h>
#include <TFile.h>
#include <vector>
#include <string>

using ROOT::VecOps::RVec;

void check_eta_correction(TString inFile = "g249_data_incoming_online_20250728_183008.root")
{

    TString repopath = gSystem->Getenv("repopath");

    ROOT::EnableImplicitMT();
    ROOT::RDataFrame df("evt", repopath + "/data/" + inFile);

    // Foot ids
    std::vector<int> detIds = {1, 2, 3, 4, 5, 6, 7, 8};

    // Histogram parameters
    const int nbx = 120;
    const double xmin = 5, xmax = 14; // foot
    const int nby = 600;
    const double ymin = 5, ymax = 14.; // los

    // Vector to store the histograms
    std::vector<ROOT::RDF::RResultPtr<TH2D>> hists;
    hists.reserve(detIds.size());

    for (auto det : detIds)
    {
        // masks
        std::string maskName = Form("mask_det%d", det);
        std::string losChar = Form("los_char_det%d", det);
        std::string footChar = Form("foot_char_det%d", det);

        // filter
        auto df_det = df
                          .Define(maskName,
                                  [](const RVec<UShort_t> &mulStrip)
                                  {
                                      return ROOT::VecOps::Map(mulStrip, [](UShort_t x)
                                                               { return (int)(x > 1); });
                                  },
                                  {"FootHitData.fMulStrip"})
                          .Define(footChar, Form("FootHitData.fZCharge[%s]", maskName.c_str()))
                          .Define(losChar,
                                  [](const RVec<double> &losZ, const RVec<int> &mask)
                                  {
                                      int n = ROOT::VecOps::Sum(mask);
                                      RVec<double> out;
                                      out.reserve(n);
                                      double val = (!losZ.empty()) ? losZ[0] : std::numeric_limits<double>::quiet_NaN();
                                      for (int i = 0; i < n; ++i)
                                          out.push_back(val);
                                      return out;
                                  },
                                  {"LosHit.fZ", maskName})
                          .Filter(Form("ROOT::VecOps::Sum(%s) > 0", maskName.c_str()))
                          .Filter("LosHit.fZ.size() == 1");

        auto h2 = df_det.Histo2D(
            {Form("h2_FootVsLos_z_det%d", det),
             Form("Det %d;Z_{FOOT};Z_{LOS}", det),
             nbx, xmin, xmax, nby, ymin, ymax},
            footChar, losChar);

        hists.push_back(h2);
    }

    // write
    TFile out(repopath + "/results/cal/" + "chargeId_perDet.root", "RECREATE");
    for (auto &h : hists)
    {
        h->Write();
    }
    out.Close();
}
