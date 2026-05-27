/** --------------------------------------------------------------------
 **
 **  Define the simulation setup for experiment G249
 **
 **  Last Update: 21/05/26
 **  Comments:
 **         - 22/07/25 : Initial setup
 **         - 21/05/26 : Switched primary generator to
 **                      R3BTheoreticalSPOMomdisGenerator. The fragment
 **                      momentum is now sampled from a theoretical SPO
 **                      distribution (TH3 in 25F rest frame), boosted
 **                      with an experimental parent-beta histogram, and
 **                      placed at a vertex drawn from an experimental
 **                      beam-spot (x,y) histogram, uniformly along the
 **                      5 cm target in z.
 **
 **  Inputs expected (TString paths set below):
 **    - fPMom3DFile  : ROOT file containing TH3  "h3_pxpypz" (GeV/c)
 **    - fVertexFile  : ROOT file containing TH2  "hVertexXY" (cm)
 **    - fBetaFile    : ROOT file containing TH1  "h"         (dimensionless)
 **
 **  Execute it as follows:
 **  root -l 'runsim.C(1000)'
 **  where 1000 means the number of events
 **
 **/

void runsim(Int_t nEvents = 100000)
{
    // Timer
    TStopwatch timer;
    timer.Start();

    // Logging
    auto logger = FairLogger::GetLogger();
    logger->SetLogVerbosityLevel("low");
    logger->SetLogScreenLevel("warn");
    logger->SetColoredLog(true);

    // System paths
    const TString workDirectory = getenv("VMCWORKDIR");
    gSystem->Setenv("GEOMPATH", workDirectory + "/geometry");
    gSystem->Setenv("CONFIG_DIR", workDirectory + "/gconfig");

    // Output files
    const TString simufile = "glad.simu.root";
    const TString parafile = "glad.para.root";

    // Input GLAD geometry
    const TString fGladGeo = "glad_v2025.1.geo.root";

    // ---------------------------------------------------------------------
    // Inputs for the custom generator R3BTheoreticalSPOMomdisGenerator
    // ---------------------------------------------------------------------
    // Adjust these three paths to point to your input files.
    const TString fPMom3DFile = "inputs/momentum3d.root";     // TH3 "h3_pxpypz" [GeV/c]
    const TString fVertexFile = "inputs/reaction_point.root"; // TH2 "hVertexXY" [cm]
    const TString fBetaFile = "inputs/beta_dist.root";        // TH1 "h"        [unitless]

    // Fragment to generate (defaults: 24O, fully stripped)
    const Int_t fFragZ = 8;
    const Int_t fFragA = 24;
    const Int_t fFragQ = 8;

    // Target half-thickness in z (cm). 5 cm target -> 2.5 cm.
    const Double_t fTargetHalfZ = 2.5;

    // Basic simulation setup
    auto run = std::make_unique<FairRunSim>();
    run->SetName("TGeant4");
    run->SetStoreTraj(false);
    run->SetMaterials("media_r3b.geo");

    auto config = std::make_unique<FairGenericVMCConfig>();
    run->SetSimulationConfig(std::move(config));
    run->SetSink(std::make_unique<FairRootFileSink>(simufile.Data()));

    // -----   Runtime data base   --------------------------------------------
    auto *rtdb = run->GetRuntimeDb();
    UInt_t runId = 1;
    rtdb->initContainers(runId);

    // ---------------------------------------------------------------------
    // Primary particle generator
    // ---------------------------------------------------------------------
    // Open the three input files and read the histograms. The TFile objects
    // are kept on the heap (not closed) because the generator stores raw
    // TH* pointers and calls GetRandom*() on them throughout the run.
    auto *fPMom = TFile::Open(fPMom3DFile, "READ");
    if (!fPMom || fPMom->IsZombie())
    {
        std::cerr << "ERROR: cannot open momentum file " << fPMom3DFile << std::endl;
        gApplication->Terminate();
    }
    auto *h3_pxpypz = dynamic_cast<TH3 *>(fPMom->Get("h3_pxpypz"));
    if (!h3_pxpypz)
    {
        std::cerr << "ERROR: TH3 'h3_pxpypz' not found in " << fPMom3DFile << std::endl;
        gApplication->Terminate();
    }
    h3_pxpypz->SetDirectory(nullptr); // decouple from the file

    auto *fVertex = TFile::Open(fVertexFile, "READ");
    if (!fVertex || fVertex->IsZombie())
    {
        std::cerr << "ERROR: cannot open vertex file " << fVertexFile << std::endl;
        gApplication->Terminate();
    }
    auto *hVertexXY = dynamic_cast<TH2 *>(fVertex->Get("hVertexXY"));
    if (!hVertexXY)
    {
        std::cerr << "ERROR: TH2 'hVertexXY' not found in " << fVertexFile << std::endl;
        gApplication->Terminate();
    }
    hVertexXY->SetDirectory(nullptr);

    auto *fBeta = TFile::Open(fBetaFile, "READ");
    if (!fBeta || fBeta->IsZombie())
    {
        std::cerr << "ERROR: cannot open beta file " << fBetaFile << std::endl;
        gApplication->Terminate();
    }
    auto *hBeta = dynamic_cast<TH1 *>(fBeta->Get("h"));
    if (!hBeta)
    {
        std::cerr << "ERROR: TH1 'h' not found in " << fBetaFile << std::endl;
        gApplication->Terminate();
    }
    hBeta->SetDirectory(nullptr);

    // We can safely close the input files now: histograms are detached.
    fPMom->Close();
    fVertex->Close();
    fBeta->Close();

    // Build and configure the generator.
    // The name passed to the named constructor is used as the prefix for the
    // truth side branches: "SPO_BetaParent", "SPO_PxRest", "SPO_PyRest",
    // "SPO_PzRest" (all Double_t, GeV/c for the rest-frame momentum).
    auto *spoGen = new R3BTheoreticalSPOMomdisGenerator("SPO");
    spoGen->SetIon(fFragZ, fFragA, fFragQ); // 24O (or change above)
    spoGen->SetRestMomentumHist3D(h3_pxpypz);
    spoGen->SetRestMomentumUnitMeV(); // TH3 axes are in MeV/c
    spoGen->SetParentBetaHist(hBeta);
    spoGen->SetBeamSpotHistXY(hVertexXY);
    spoGen->SetTargetHalfThicknessZ(fTargetHalfZ);
    spoGen->SetTargetCenter(0., 0., 0.);
    spoGen->EnableTruthOutput(kTRUE);              // default; shown for clarity
    spoGen->SetTruthOutputFile("glad.truth.root"); // sits next to glad.simu.root

    auto primGen = std::make_unique<FairPrimaryGenerator>();
    primGen->AddGenerator(spoGen); // ownership transferred
    run->SetGenerator(primGen.release());

    // Register the fragment ion in the PDG database. This MUST happen
    // before run->Init(), otherwise Geant4-VMC builds its particle table
    // without the ion and FairPrimaryGenerator drops every primary with
    // "PDG code <N> not found in database".
    run->AddNewIon(new FairIon(Form("Ion_Z%d_A%d", fFragZ, fFragA),
                               fFragZ, fFragA, fFragQ));

    // Geometry: Cave
    auto cave = std::make_unique<R3BCave>("CAVE");
    cave->SetGeometryFileName("r3b_cave.geo");
    run->AddModule(cave.release());

    // Geometry: Target + vacuum chamber
    R3BTra *tra = new R3BTra("target_area_v2025_5cm.geo.root", {0., 0., 0.});
    tra->SetEnergyCut(1e-6); // 1 keV
    run->AddModule(tra);

    // Geomtry: Fibers
    Bool_t fFi30 = true; // Fi30 detector
    TString fFi30Geo = "fi30_v2022.1.geo.root";

    Bool_t fFi31 = true; // Fi31 detector
    TString fFi31Geo = "fi31_v2022.1.geo.root";

    Bool_t fFi32 = true; // Fi32 detector
    TString fFi32Geo = "fi32_v2022.1.geo.root";

    Bool_t fFi33 = true; // Fi33 detector
    TString fFi33Geo = "fi33_v2022.1.geo.root";

    const double angle = 18. * TMath::DegToRad();

    // run->AddModule(new R3BFiber("Fi30", fFi30Geo, DetectorId::kFI30,
    //                             {-140.760, 0.000, 718.940},
    //                             {"", -90., +18, 90.}));

    // run->AddModule(new R3BFiber("Fi31", fFi31Geo, DetectorId::kFI31,
    //                             {-142.584, 0.000, 805.415},
    //                             {"", -90., +18, 90.}));

    // run->AddModule(new R3BFiber("Fi32", fFi32Geo, DetectorId::kFI32,
    //                             {-132.270, 0.000, 692.780},
    //                             {"", -90., +18, 90.}));

    // run->AddModule(new R3BFiber("Fi33", fFi33Geo, DetectorId::kFI33,
    //                             {-185.566, 0.000, 775.895},
    //                             {"", -90., +18, 90.}));

    run->AddModule(new R3BFiber("Fi30", fFi30Geo, DetectorId::kFI30,
                                {-141.252, 0.000, 722.299},
                                {"", -90., +18, 90.}));

    run->AddModule(new R3BFiber("Fi31", fFi31Geo, DetectorId::kFI31,
                                {-191.109, 0.000, 794.843},
                                {"", -90., +18, 90.}));

    run->AddModule(new R3BFiber("Fi32", fFi32Geo, DetectorId::kFI32,
                                {-132.754, 0.000, 696.144},
                                {"", -90., +18, 90.}));

    run->AddModule(new R3BFiber("Fi33", fFi33Geo, DetectorId::kFI33,
                                {-138.490, 0.000, 794.696},
                                {"", -90., +18, 90.}));

    // Geometry: TofD
    Bool_t fTofD = true; // TofD detector
    TString fTofDGeo = "tofd_v2025.6.geo.root";
    // run->AddModule(new R3BTofD(fTofDGeo,
    //                            {-218.040, 0.000, 956.780},
    //                            {"", -90., +18, 90.}));

    run->AddModule(new R3BTofD(fTofDGeo,
                               {-218.536, 0.000, 960.142},
                               {"", -90., +18, 90.}));

    // Geometry: GLAD
    run->AddModule(new R3BGladMagnet(fGladGeo.Data()));

    // GLAD Filed
    auto *GladField = new R3BGladFieldMap("R3BGladMap");
    GladField->SetFieldfromCurrent(2668.0); // Current in Amperes
    run->SetField(GladField);

    // Init
    run->Init();

    // Save field parameters
    auto *fieldPar = dynamic_cast<R3BFieldPar *>(rtdb->getContainer("R3BFieldPar"));
    fieldPar->SetParameters(GladField);
    fieldPar->setChanged();

    // Output file with parameters
    auto parOut = std::make_unique<FairParRootFileIo>(true);
    parOut->open(parafile.Data());
    rtdb->setOutput(parOut.release());
    rtdb->saveOutput();
    rtdb->print();

    // Simulate
    run->Run(nEvents);

    // Report
    timer.Stop();
    std::cout << "Real time: " << timer.RealTime() << "s, CPU time: " << timer.CpuTime() << "s" << std::endl;
    std::cout << "Macro finished successfully." << std::endl;

    // Clean up the detached histograms (optional; process is about to exit).
    delete h3_pxpypz;
    delete hVertexXY;
    delete hBeta;
}