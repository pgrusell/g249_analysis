#include "neuland_alignment.hh"

void neuland_alignment()
{
    // Get the offsets from the roofile
    auto neul = std::make_unique<neulandAlignment>("tof_alignment_losNeuland_v1.root", std::vector<TString>{
                                                                                           "/nucl_lustre/g249/root_neuland/g249_all_det_offline_20250828_210751.root",
                                                                                           "/nucl_lustre/g249/root_neuland/g249_all_det_offline_20250828_210751_1.root",
                                                                                           "/nucl_lustre/g249/root_neuland/g249_all_det_offline_20250828_210751_2.root"});
    neul->calculateOffsets();

    neul->checkAlignment();

    // Apply the offsets

    neul->buildHistogram(std::vector<TString>{
                             "/nucl_lustre/g249/root_neuland/g249_all_det_offline_20250828_210751.root",
                             "/nucl_lustre/g249/root_neuland/g249_all_det_offline_20250828_210751_1.root",
                             "/nucl_lustre/g249/root_neuland/g249_all_det_offline_20250828_210751_2.root"},
                         true);
}