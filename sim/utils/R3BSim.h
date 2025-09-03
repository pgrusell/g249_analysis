#pragma once
#include "ISimGeometry.h"
#include "ISimGenerator.h"
#include <FairRunSim.h>
#include <FairParRootFileIo.h>
#include <FairRuntimeDb.h>
#include <TStopwatch.h>
#include <TRandom3.h>
#include <TVirtualMC.h>
#include <FairLogger.h>
#include <TSystem.h>
#include <memory>
#include <string>
#include <optional>

class R3BSim
{
public:
    struct Options
    {
        std::string transport = "TGeant4";
        std::string outFile = "./sim.root";
        std::string parFile = "./par.root";
        std::string materials = "media_r3b.geo";
        int randomSeed = 0;
        int maxNSteps = -15000;
        bool storeTraj = false;
        bool userPList = false;
        std::string logVerbosity = "LOW"; // FairLogger
        std::string logScreen = "INFO";
    };

    explicit R3BSim(Options opts = {}) : fOpts(std::move(opts)) {}

    R3BSim &withGeometry(std::unique_ptr<ISimGeometry> geo)
    {
        fGeo = std::move(geo);
        return *this;
    }
    R3BSim &withGenerator(std::unique_ptr<ISimGenerator> gen)
    {
        fGen = std::move(gen);
        return *this;
    }

    SimPaths resolvePaths() const
    {
        SimPaths p;
        const char *vmc = gSystem->Getenv("VMCWORKDIR");
        p.vmcworkdir = vmc ? vmc : "";
        p.pardir = p.vmcworkdir + "/../params/";
        p.geomdir = p.vmcworkdir + "/geometry/";
        p.configdir = p.vmcworkdir + "/gconfig/";
        return p;
    }

    void run(int nEvents)
    {
        if (!fGeo || !fGen)
        {
            ::Error("R3BSim", "Debes configurar geometría y generador antes de run().");
            return;
        }

        auto paths = resolvePaths();

        // Normalize double bars to bars (/u//land -> /u/land)
        auto norm = [](std::string &s)
        { for (size_t i=1;i<s.size();++i) if(s[i]=='/' && s[i-1]=='/') s.erase(i--,1); };
        norm(paths.pardir);
        norm(paths.geomdir);
        norm(paths.configdir);

        // Export env for R3BRoot
        gSystem->Setenv("GEOMPATH", paths.geomdir.c_str());
        gSystem->Setenv("CONFIG_DIR", paths.configdir.c_str());

        TStopwatch timer;
        timer.Start();

        // --- FairRunSim boilerplate
        auto run = std::make_unique<FairRunSim>();
        run->SetName(fOpts.transport.c_str());
        run->SetOutputFile(fOpts.outFile.c_str());
        auto *rtdb = run->GetRuntimeDb();

        // Use special physics if needed
        if ((fOpts.userPList) && (fOpts.transport.c_str().CompareTo("TGeant4") == 0))
        {
            run->SetUserConfig("g4R3bConfig.C");
            run->SetUserCuts("SetCuts.C");
        }

        run->SetMaterials(fOpts.materials.c_str());

        // Logging
        FairLogger::GetLogger()->SetLogVerbosityLevel(fOpts.logVerbosity.c_str());
        FairLogger::GetLogger()->SetLogScreenLevel(fOpts.logScreen.c_str());

        // Geometry
        fGeo->configure(*run, paths);

        // Generator
        run->SetGenerator(fGen->build());
        run->SetStoreTraj(fOpts.storeTraj);

        // Init
        run->Init();
        TVirtualMC::GetMC()->SetRandom(new TRandom3(fOpts.randomSeed));
        TVirtualMC::GetMC()->SetMaxNStep(fOpts.maxNSteps);

        // Save params
        Bool_t kMerged = kTRUE;
        auto *parOut = new FairParRootFileIo(kMerged);
        parOut->open(fOpts.parFile.c_str());
        rtdb->setOutput(parOut);
        rtdb->saveOutput();
        rtdb->print();

        // Run
        if (nEvents > 0)
        {
            run->Run(nEvents);
        }

        timer.Stop();
        std::cout << "\nMacro finished succesfully.\n"
                  << "Output file: " << fOpts.outFile << "\n"
                  << "Parameter file: " << fOpts.parFile << "\n"
                  << "Real time " << timer.RealTime() << " s, CPU time " << timer.CpuTime() << " s\n\n";
    }

private:
    Options fOpts;
    std::unique_ptr<ISimGeometry> fGeo;
    std::unique_ptr<ISimGenerator> fGen;
};