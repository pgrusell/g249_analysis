// Cling configuration
#pragma cling add_include_path("../utils")
#pragma cling add_include_path("neutronDecay")
#pragma cling add_include_path(".")

#include "../utils/R3BSim.h"
#include "GladNeulandGeometry.h"
#include "../utils/AsciiFileGenerator.h"

void neutronDecaySimulation(Int_t nEvents = 1000,
                            TString asciiFile = "O24_2n_res1.out")
{

    // paths
    TString repopath = getenv("repopath");

    GladNeulandGeometry::Options gopt;
    gopt.caveGeoFile = "r3b_cave_vacuum.geo";
    gopt.gladGeoFile = "glad_v2025.1.geo.root";
    gopt.neulandPlanes = 13;
    gopt.neulandPos = TGeoTranslation(0., 0., 1650.);
    gopt.fieldScale = -0.82;

    R3BSim::Options sopt;
    sopt.transport = "TGeant4";
    sopt.outFile = (repopath + "/results/sim" + "/sim_2ndecay_2res.root").Data();
    sopt.parFile = (repopath + "/results/sim" + "/par_2ndecay_2res.root").Data();
    sopt.materials = "media_r3b.geo";
    sopt.storeTraj = kTRUE;
    sopt.maxNSteps = -15000;
    sopt.userPList = false;

    R3BSim sim(sopt);
    sim.withGeometry(std::make_unique<GladNeulandGeometry>(gopt))
        .withGenerator(std::make_unique<AsciiFileGenerator>((repopath + "/sim/gen/" + asciiFile).Data(), 0., 0., 0.))
        .run(nEvents);
}
