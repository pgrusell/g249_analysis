#pragma once
#include "ISimGeometry.h"
#include <TVector3.h>
#include <string>

class GladNeulandGeometry : public ISimGeometry
{
public:
    struct Options
    {
        std::string caveGeoFile = "r3b_cave_vacuum.geo";
        std::string gladGeoFile = "glad_v2025.1.geo.root";

        int neulandPlanes = 13;
        TGeoTranslation neulandPos = {0., 0., 1650.}; // mm

        // Magnetic field
        double fieldScale = -0.82;
    };

    GladNeulandGeometry() : fOpt() {}
    explicit GladNeulandGeometry(Options opt) : fOpt(std::move(opt)) {}

    void configure(FairRunSim &run, const SimPaths &) override
    {
        // --- CAVE ---
        auto *cave = new R3BCave("CAVE");
        cave->SetGeometryFileName(fOpt.caveGeoFile.c_str());
        run.AddModule(cave);

        // --- GLAD (magnet geometry) ---
        run.AddModule(new R3BGladMagnet(fOpt.gladGeoFile.c_str()));

        // --- NeuLAND ---
        run.AddModule(new R3BNeuland(fOpt.neulandPlanes, fOpt.neulandPos));

        run.AddTask(
            new R3BNeulandDigitizer(R3BNeulandDigitizer::Options::neulandTamex));
        run.AddTask(new R3BNeulandClusterFinder());
        run.AddTask(new R3BNeulandPrimaryInteractionFinder());
        run.AddTask(new R3BNeulandPrimaryClusterFinder());

        // --- GLAD field map ---
        auto *magField = new R3BGladFieldMap("R3BGladMap");
        magField->SetScale(fOpt.fieldScale);
        run.SetField(magField);

        // save parameters
        auto *rtdb = run.GetRuntimeDb();
        auto *fieldPar = (R3BFieldPar *)rtdb->getContainer("R3BFieldPar");
        if (fieldPar)
        {
            fieldPar->SetParameters(magField);
            fieldPar->setChanged();
        }
    }

private:
    Options fOpt;
};