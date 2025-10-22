#include "footAna.h"

void footAnalysis()
{
    // auto footana = std::make_unique<footAna>("/nucl_lustre/pablogrusell/g249/root_files//g249_all_det_offline_0001_20251021_152822.root");
    // auto footana = std::make_unique<footAna>("/nucl_lustre/pablogrusell/g249/root_files//g249_all_det_offline_0001_20251021_181004.root");
    auto footana = std::make_unique<footAna>("/nucl_lustre/pablogrusell/g249/root_files//g249_all_det_offline_0001_20251021_195505.root");
    footana->fillHistograms();
    footana->saveFigures();
}