#include "../utils/R3BSim.h"
#include "GladNeulandGeometry.h"
#include "../utils/AsciiFileGenerator.h"

void neutronDecaySimulation(Int_t nEvents = 10000,
                            TString asciiFile = "gen_test.out")
{

    // paths
    TString repopath = getenv("repopath");

    GladNeulandGeometry::Options gopt;
    gopt.caveGeoFile = "r3b_cave_vacuum.geo";
    gopt.gladGeoFile = "glad_v2025.1.geo.root";
    gopt.neulandPlanes = 13;
    gopt.neulandPos = TVector3(0., 0., 1650.);
    gopt.fieldScale = -0.82;

    R3BSim::Options sopt;
    sopt.transport = "TGeant4";
    sopt.outFile = (repopath + "/results/sim/" + "/sim_neutronDecay.root").Data();
    sopt.parFile = (repopath + "/results/sim/" + "/par_neutronDecay.root").Data();
    sopt.materials = "media_r3b.geo";
    sopt.storeTraj = kFALSE;
    sopt.maxNSteps = -15000;
    sopt.userPList = false;

    R3BSim sim(sopt);
    sim.withGeometry(std::make_unique<GladNeulandGeometry>(gopt))
        .withGenerator(std::make_unique<AsciiFileGenerator>((repopath + "/sim/gen/" + asciiFile).Data(), 0., 0., 0.))
        .run(nEvents);
}