#include "neuland_alignment.hh"

// OPTION 0:
// - Takes the rootfile as input
// - Calculates ToF vs paddle - this can be already provided
// - Calculates offsets
// - Calculates ToF vs paddle (corrected)
// - Runs the analysis and saves everything

// OPTION 1:
// - Takes as input ToF vs paddle (corrected and uncorrected)
// - Runs the analysis

// OPTION 2:
// Hybrid option (deprecated)

void neuland_alignment(int mode = 0)
{
    /*
    // Histogram with the first parameters
    auto neul = std::make_unique<neulandAlignment>("h_v3_aligned.root", std::vector<TString>{
                                                                            "/nucl_lustre/g249/root_neuland/g249_all_det_offline_20250828_210751.root",
                                                                            "/nucl_lustre/g249/root_neuland/g249_all_det_offline_20250828_210751_1.root",
                                                                            "/nucl_lustre/g249/root_neuland/g249_all_det_offline_20250828_210751_2.root"});

    neul->checkAlignment(10, false);

    // Retrive the offsets from the .txt file
    // neul->calculateOffsets();
    neul->setOffsetsFromTxt();
    */

    // Apply the offsets
    /*
    neul->buildHistogram(std::vector<TString>{
                             "/nucl_lustre/g249/root_neuland/g249_all_det_offline_20250828_210751.root",
                             "/nucl_lustre/g249/root_neuland/g249_all_det_offline_20250828_210751_1.root",
                             "/nucl_lustre/g249/root_neuland/g249_all_det_offline_20250828_210751_2.root"},
                         true);

    neul->checkAlignment();
    */

    std::vector<TString> rootFile = {
        "/nucl_lustre/g249/root_neuland/g249_all_det_offline_20250828_210751.root",
        "/nucl_lustre/g249/root_neuland/g249_all_det_offline_20250828_210751_1.root",
        "/nucl_lustre/g249/root_neuland/g249_all_det_offline_20250828_210751_2.root"};

    switch (mode)
    {
    case 0:
    {

        auto neul = std::make_unique<neulandAlignment>("", rootFile);
        neul->setOffsetsFromTxt();
        neul->buildHistogram(rootFile, true);

        break;
    }

    case 1:
    {

        auto neul = std::make_unique<neulandAlignment>("");
        neul->setCorrectedFromRoot("histogram_builtCorr.root");

        break;
    }

    case 2:
    {
        std::cerr << "[WARN] This option is deprecated\n";
        auto neul = std::make_unique<neulandAlignment>("h_v3_aligned.root", rootFile);
        neul->setOffsetsFromTxt("offsets_v1.txt");
        neul->buildHistogram(rootFile, true);
        break;
    }
    }
}