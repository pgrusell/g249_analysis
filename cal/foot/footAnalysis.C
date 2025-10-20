#include "footAna.h"

void footAnalysis()
{
    // auto footana = std::make_unique<footAna>("/nucl_lustre/pablogrusell/g249/root_files/foothitData_run152*");
    auto footana = std::make_unique<footAna>("/nucl_lustre/pablogrusell/g249/root_files//g249_all_det_offline_0001_20251020_115737.root");
    footana->fillHistograms();
    // footana->saveFigures();
}