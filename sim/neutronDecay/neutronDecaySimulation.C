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

void neutronDecaySimulation(Int_t nEvents = 10)
{
    // ----------- Configuration area ----------------------------------

    TString OutFile = "sim.root"; // Output file for data
    TString ParFile = "par.root"; // Output file for params

    TString fSetupFile = "Setup.par"; // Input file with the detector positions

    Bool_t fVis = true;             // Store tracks for visualization
    Bool_t fUserPList = false;      // Use of R3B special physics list
    Bool_t fR3BMagnet = false;      // Magnetic field definition
    Bool_t fCalifaDigitizer = true; // Apply hit digitizer task
    Bool_t fCalifaHitFinder = true; // Apply hit finder task

    // MonteCarlo engine: TGeant3, TGeant4, TFluka  --------------------
    TString fMC = "TGeant4";

    // Event generator type: box for particles or ascii&inclroot for p2p-fission
    TString generator1 = "box";
    TString generator2 = "ascii";
    TString generator3 = "inclroot";
    TString generator4 = "roofile";
    TString fGenerator = generator2;

    Bool_t fGlad = true; // Glad Magnet
    TString fGladGeo = "glad_v2025.1.geo.root";

    Bool_t fCalifa = false; // Califa Calorimeter
                            // std::string fCalifaGeo = "califa_full.geo.root";
    // std::string fCalifaGeo = "califa_full_stru.geo.root";
    std::string fCalifaGeo = "califa_v2022.5_stru.geo.root";
    Int_t fCalifaGeoVer = 2024;

    TString fTrackerGeo = "target_area_v2022_5cm.geo.root";

    TString fNeuLandGeo = "neuland_v3_13dp.geo.root";

    // Input event file in the case of ascii/root generator
    TString fEventFile;
    if (fGenerator.CompareTo("ascii") == 0)
        fEventFile = "/nucl_lustre/pablogrusell/g249/g249_analysis/sim/gen/gen_test.out";
    else if (fGenerator.CompareTo("inclroot") == 0)
        fEventFile = "/nucl_lustre/pablogrusell/g249/g249_analysis/sim/gen/gen_test.root";

    Int_t fFieldMap = -1;          // Magentic field map selector
    Double_t fMeasCurrent = 2000.; // Magnetic field current
    Float_t fFieldScale = -0.82;   // Magnetic field scale factor

    // ---------  Detector selection: true - false ---------------------
    // ---- R3B and SOFIA detectors as well as passive elements

    // ---- End of Configuration area   ---------------------------------------

    // ---- Stable part   -----------------------------------------------------
    TString dir = gSystem->Getenv("VMCWORKDIR");

    TString r3b_geomdir = dir + "/geometry/";
    gSystem->Setenv("GEOMPATH", r3b_geomdir.Data());
    r3b_geomdir.ReplaceAll("//", "/");

    TString r3b_confdir = dir + "/gconfig/";
    gSystem->Setenv("CONFIG_DIR", r3b_confdir.Data());
    r3b_confdir.ReplaceAll("//", "/");

    // ----    Debug option   -------------------------------------------------
    gDebug = 0;

    // -----   Timer   --------------------------------------------------------
    TStopwatch timer;
    timer.Start();

    // -----   Create simulation run   ----------------------------------------
    FairRunSim *run = new FairRunSim();
    run->SetName(fMC);                           // Transport engine
    run->SetSink(new FairRootFileSink(OutFile)); // Output file

    // -----   Runtime data base   --------------------------------------------
    FairRuntimeDb *rtdb = run->GetRuntimeDb();

    // -----   Load detector parameters    ------------------------------------

    // ----- Containers

    UInt_t runId = 1;
    rtdb->initContainers(runId);

    // -----   R3B Special Physics List in G4 case
    if ((fUserPList) && (fMC.CompareTo("TGeant4") == 0))
    {
        run->SetUserConfig("g4R3bConfig.C");
        run->SetUserCuts("SetCuts.C");
    }

    // -----   Create media   -------------------------------------------------
    run->SetMaterials("media_r3b.geo"); // Materials

    // -----   Create R3B geometry --------------------------------------------

    // Cave definition
    FairModule *cave = new R3BCave("CAVE");
    cave->SetGeometryFileName("r3b_cave_vacuum.geo");
    run->AddModule(cave);

    fFieldMap = 1;
    fR3BMagnet = true;

    run->AddModule(new R3BGladMagnet(fGladGeo));

    run->AddModule(new R3BNeuland(13, {0., 0., 1650.}));

    // comentar esto si solo se quiere calcular eficiencias para neuland ----
    /*
    run->AddTask(
        new R3BNeulandDigitizer(R3BNeulandDigitizer::Options::neulandTamex));
    run->AddTask(new R3BNeulandClusterFinder());
    run->AddTask(new R3BNeulandPrimaryInteractionFinder());
    run->AddTask(new R3BNeulandPrimaryClusterFinder());
    */
    // ---------------------------------------------------------------------

    R3BGladFieldMap *magField = new R3BGladFieldMap("R3BGladMap");
    magField->SetScale(fFieldScale);

    run->SetField(magField);
    R3BFieldPar *fieldPar = (R3BFieldPar *)rtdb->getContainer("R3BFieldPar");
    fieldPar->SetParameters(magField);
    fieldPar->setChanged();

    // ---- End of field map section

    // -----   Create PrimaryGenerator   --------------------------------------

    // 1 - Create the Main API class for the Generator
    FairPrimaryGenerator *primGen = new FairPrimaryGenerator();

    if (fGenerator.CompareTo("box") == 0)
    {
        // Define the BOX generator
        Int_t pdgId = 2212;     // proton beam
        Double32_t theta1 = 7.; // polar angle distribution
        Double32_t theta2 = 90.;
        Double32_t momentum = 0.70;
        FairBoxGenerator *boxGen = new FairBoxGenerator(pdgId, 1);
        boxGen->SetThetaRange(theta1, theta2);
        boxGen->SetPRange(momentum, 2.0 * momentum);
        boxGen->SetPhiRange(0., 360.);
        boxGen->SetXYZ(0.0, 0.0, 0.0);
        // primGen->AddGenerator(boxGen);

        // 208-Pb fragment
        FairIonGenerator *ionGen =
            new FairIonGenerator(92, 238, 92, 1, 0., 0., 1.09, 0., 0., 0.);
        // primGen->AddGenerator(ionGen);

        auto calgen = new R3BCALIFATestGenerator(22, 1);
        calgen->SetCosTheta();
        calgen->SetThetaRange();
        calgen->SetPhiRange();
        // calgen->SetPRange(momentum, 2.0 * momentum);
        calgen->SetNuclearDecayChain();
        calgen->SetLorentzBoost(0.773);
        calgen->SetDecayChainPoint(0.002, 1);
        // calgen->SetDecayChainPoint(0.0003, 1);
        // calgen->SetDecayChainPoint(0.0005, 1.);
        // calgen->SetDecayChainPoint(0.002, 1.);
        // calgen->SetDecayChainPoint(0.0013325, 0.5);
        primGen->AddGenerator(calgen);
    }

    if (fGenerator.CompareTo("ascii") == 0)
    {
        R3BAsciiGenerator *gen =
            new R3BAsciiGenerator((fEventFile).Data());
        // gen->SetXYZ(targetPar->GetPosX(), targetPar->GetPosY(),
        //             targetPar->GetPosZ());
        gen->SetDxDyDz(0., 0., 0.);
        primGen->AddGenerator(gen);
    }

    if (fGenerator.CompareTo("inclroot") == 0)
    {
        R3BINCLRootGenerator *gen = new R3BINCLRootGenerator(fEventFile.Data());
        //   gen->SetOnlyFission(kTRUE);
        // gen->SetOnlyP2pFission(kTRUE);
        gen->SetXYZ(0., 0., -1.8);
        gen->SetDxDyDz(0.3, 0.3, 2.5);
        primGen->AddGenerator(gen);
    }

    run->SetGenerator(primGen);

    //-------Set visualisation flag to true------------------------------------
    run->SetStoreTraj(fVis);

    FairLogger::GetLogger()->SetLogVerbosityLevel("low");

    // -----   Initialize simulation run   ------------------------------------
    run->Init();

    // -----   Runtime database   ---------------------------------------------
    Bool_t kParameterMerged = kTRUE;
    FairParRootFileIo *parOut = new FairParRootFileIo(kParameterMerged);
    parOut->open(ParFile.Data());
    rtdb->setOutput(parOut);
    rtdb->saveOutput();
    rtdb->print();

    // -----   Start run   ----------------------------------------------------
    if (nEvents > 0)
        run->Run(nEvents);

    // -----   Finish   -------------------------------------------------------
    timer.Stop();
    Double_t rtime = timer.RealTime() / 60.;
    Double_t ctime = timer.CpuTime() / 60.;
    cout << endl
         << endl;
    cout << "Macro finished succesfully." << endl;
    cout << "Output file is " << OutFile << endl;
    cout << "Parameter file is " << ParFile << endl;
    cout << "Real time " << rtime << " min, CPU time " << ctime << " min" << endl
         << endl;

    cout << " Test passed" << endl;
    cout << " All ok " << endl;
    gApplication->Terminate();
}
