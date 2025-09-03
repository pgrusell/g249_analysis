#pragma once
#include "ISimGeometry.h"

#include <FairRunSim.h>
#include <R3BCave.h>
#include <R3BGladMagnet.h>
#include <R3BGladFieldMap.h>
#include <R3BFieldPar.h>
#include <R3BNeuland.h>

#include <FairRuntimeDb.h>
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
        TVector3 neulandPos = {0., 0., 1650.}; // mm

        // Magnetic field
        double fieldScale = -0.82;
    };

    explicit GladNeulandGeometry(Options opt = {}) : fOpt(std::move(opt)) {}

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

        // --- GLAD field map ---
        auto *magField = new R3BGladFieldMap("R3BGladMap");
        magField->SetScale(fOpt.fieldScale);
        run.SetField(magField);

        // Empaqueta parámetros del campo en la runtime DB (igual que macro)
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