#include "CrossSections.cpp"

void runCrossSections(TString mode = "scratch",
                      TString cacheFile = "fitcache.txt",
                      double kPID = 2.5,
                      double kOpa = 3.5,
                      int nIterations = 1000)
{
    const TString fileFragment = "/nucl_lustre/pablogrusell/g249/g249_analysis/results/final/data_24O_test3.root";
    const TString fileUnreacted = "/nucl_lustre/pablogrusell/g249/g249_analysis/results/final/data_25F_test3.root";

    CrossSections xs;
    xs.SetVerbose(true);
    xs.SetKEllipseFragment(2.5);
    xs.SetKEllipseUnreacted(2.5);
    xs.SetKOpa(3.5);

    // Load trees (needed by all modes)
    xs.SetTrees(fileFragment, fileUnreacted);

    // Try to load cache if it exists
    std::ifstream testCache(cacheFile.Data());
    bool cacheExists = testCache.good();
    testCache.close();

    // Mode scratch
    if (mode == "scratch")
    {
        std::cout << "\n====== MODE: scratch ======\n\n";

        xs.RunFit2DForNucleus("25F", kPID);
        xs.RunFit2DForNucleus("23O", kPID);
        xs.SaveFitCache(cacheFile);

        double sigma = xs.ComputeCrossSection();
        std::cout << "\n>>> Cross section = " << sigma << " mb <<<\n\n";
    }
    // Mode: fromcache
    else if (mode == "fromcache")
    {
        std::cout << "\n====== MODE: fromcache ======\n\n";

        xs.LoadFitCache(cacheFile);

        if (xs.HasCachedFit("25F"))
        {
            auto p = xs.GetCachedFit("25F");
            std::cout << "25F: muZ=" << p.muZ << " sigmaZ=" << p.sigmaZ
                      << " muAoQ=" << p.muAoQ << " sigmaAoQ=" << p.sigmaAoQ << "\n";
        }
        if (xs.HasCachedFit("24O"))
        {
            auto p = xs.GetCachedFit("24O");
            std::cout << "24O: muZ=" << p.muZ << " sigmaZ=" << p.sigmaZ
                      << " muAoQ=" << p.muAoQ << " sigmaAoQ=" << p.sigmaAoQ << "\n";
        }

        double sigma = xs.ComputeCrossSection();
        std::cout << "\n>>> Cross section = " << sigma << " mb <<<\n\n";
    }
    // Mode: plots
    else if (mode == "plots")
    {
        std::cout << "\n====== MODE: plots ======\n\n";

        if (cacheExists)
        {
            xs.LoadFitCache(cacheFile);
            std::cout << "[INFO] Loaded fit cache for plot parameters.\n";
        }

        // Uncomment the plots you want:
        xs.PlotNbOfYieldsVsK();
        xs.PlotOutgoingVsKopa();
    }
    // Mode: systematics
    else if (mode == "systematics")
    {
        std::cout << "\n====== MODE: systematics (N=" << nIterations << ") ======\n\n";

        if (cacheExists)
        {
            xs.LoadFitCache(cacheFile);
            std::cout << "[INFO] Loaded fit cache for systematic parameters.\n";
        }

        xs.SetVerbose(false);

        TH1F *hSyst = xs.ComputeCrossSectionSystematics(
            nIterations,
            /*kEllipseMin=*/2.0, /*kEllipseMax=*/3.0,
            /*kOpaMin=*/2.5, /*kOpaMax=*/3.5,
            /*effP2PMin=*/0.566313, /*effP2PMax=*/0.568788);

        if (hSyst)
        {
            std::cout << "\n>>> Systematic XS: mean = " << hSyst->GetMean()
                      << " mb, RMS = " << hSyst->GetRMS() << " mb <<<\n\n";
        }
    }
    // Mode: analytical
    else if (mode == "analytical")
    {
        std::cout << "\n====== MODE: analytical ======\n\n";

        if (cacheExists)
        {
            xs.LoadFitCache(cacheFile);
            std::cout << "[INFO] Loaded fit cache.\n";
        }

        auto result = xs.ComputeCrossSectionAnalytical(
            /*dNp2p=*/69,
            /*dNproj=*/3700,
            /*dEffP2P=*/0.0003,
            /*dEffFrag=*/0.0022);
    }
    else
    {
        std::cerr << "[ERROR] Unknown mode '" << mode
                  << "'. Use: scratch, fromcache, plots, systematics, analytical\n";
    }
}