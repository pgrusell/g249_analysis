#include "footAna.h"

void footAnalysis()
{
    auto footana = std::make_unique<footAna>("/nucl_lustre/pablogrusell/g249/g249_analysis/root_files/etaCorrection/g249_data_incoming_online_20250728_180423.root");
    footana->fillHistograms();
    footana->saveFigures();
}