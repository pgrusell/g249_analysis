/** --------------------------------------------------------------------
 **
 **  Define the simulation setup for experiment G249
 **  Author: <j.l.rodriguez.sanchez@udc.es>
 **
 **  Last Update: 22/07/25
 **  Comments:
 **         - 22/07/25 : Initial setup
 **
 **  Configuration:
 **  (1) Select the right generator "fGenerator"
 **  (2) Select the detectors that you wish for the simulation, for instance,
 *"fCalifa = true"
 **  (3) Look at the file "s455_setup.par" that the positions of your detectors
 *are right
 **
 **  Execute it as follows:
 **  root -l 'runsim.C(1000)'
 **  where 1000 means the number of events
 **
 **/

void runsim(Int_t nEvents = 10)
{
    // Timer
    TStopwatch timer;
    timer.Start();

    // Logging
    auto logger = FairLogger::GetLogger();
    // logger->SetLogVerbosityLevel("low");
    // logger->SetLogScreenLevel("warn");
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

    // Basic simulation setup
    auto run = std::make_unique<FairRunSim>();
    run->SetName("TGeant4");
    run->SetStoreTraj(true);
    run->SetMaterials("media_r3b.geo");

    auto config = std::make_unique<FairGenericVMCConfig>();
    run->SetSimulationConfig(std::move(config));
    run->SetSink(std::make_unique<FairRootFileSink>(simufile.Data()));

    // -----   Runtime data base   --------------------------------------------
    auto *rtdb = run->GetRuntimeDb();
    UInt_t runId = 1;
    rtdb->initContainers(runId);

    // Primary particle generator
    // auto ionGen = std::make_unique<FairIonGenerator>(9, 25, 9, 1, 0., 0., 1.371, 0., 0., 0.);
    auto ionGen = std::make_unique<FairIonGenerator>(8, 24, 8, 1, 0., 0., 1.2516, 0., 0., 0.);
    auto primGen = std::make_unique<FairPrimaryGenerator>();
    primGen->AddGenerator(ionGen.release());
    run->SetGenerator(primGen.release());

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
}