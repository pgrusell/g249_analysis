/******************************************************************************
 *   Copyright (C) 2019 GSI Helmholtzzentrum für Schwerionenforschung GmbH    *
 *   Copyright (C) 2019-2025 Members of R3B Collaboration                     *
 *                                                                            *
 *             This software is distributed under the terms of the            *
 *                 GNU General Public Licence (GPL) version 3,                *
 *                    copied verbatim in the file "LICENSE".                  *
 *                                                                            *
 * In applying this license GSI does not waive the privileges and immunities  *
 * granted to it by virtue of its status as an Intergovernmental Organization *
 * or submit itself to any jurisdiction.                                      *
 ******************************************************************************/

// Simulation of CALIFA p2p with the G-249 vacuum chamber + FOOTs geometry.

void califap2psim(Int_t kIter = 1)
{

    int nEvents = 250000;

    TString transport = "TGeant4";
    TString outFile = Form("./sim_%i.root", kIter);
    TString parFile = Form("./par_%i.root", kIter);

    Int_t randomSeed = 0; // 0 for time-dependent random numbers

    // ------------------------------------------------------------------------
    TString dir = gSystem->Getenv("VMCWORKDIR");
    TString pardir = dir + "/../params/";
    pardir.ReplaceAll("//", "/");

    TString r3b_geomdir = dir + "/geometry/";
    gSystem->Setenv("GEOMPATH", r3b_geomdir.Data());
    r3b_geomdir.ReplaceAll("//", "/");

    TString r3b_pardir = pardir + "/geometry/";
    r3b_pardir.ReplaceAll("//", "/");

    TString r3b_confdir = dir + "/gconfig/";
    gSystem->Setenv("CONFIG_DIR", r3b_confdir.Data());
    r3b_confdir.ReplaceAll("//", "/");

    // -----   Timer   --------------------------------------------------------
    TStopwatch timer;
    timer.Start();

    // -----   Create simulation run   ----------------------------------------
    FairRunSim *run = new FairRunSim();
    run->SetName(transport);            // Transport engine
    run->SetOutputFile(outFile.Data()); // Output file
    FairRuntimeDb *rtdb = run->GetRuntimeDb();

    FairParAsciiFileIo *parIo1 = new FairParAsciiFileIo();
    parIo1->open("califaDigi.par", "in");
    rtdb->setFirstInput(parIo1);
    rtdb->print();

    rtdb->getContainer("califaCrystalPars4Sim");

    UInt_t runId = 1;
    rtdb->initContainers(runId);

    // -----   Create media   -------------------------------------------------
    run->SetMaterials("media_r3b.geo"); // Materials

    // -----   Create R3B geometry --------------------------------------------
    // R3B Cave definition
    FairModule *cave = new R3BCave("CAVE");
    cave->SetGeometryFileName("r3b_cave.geo");
    run->AddModule(cave);

    //////////////////// DETECTORS //////////////////////////////

    // ---------------LH2 target + vacuum chamber + FOOT --------
    auto *foots = new R3BTra("target_area_v2025_5cm.geo.root");
    run->AddModule(foots);

    // --------------- CALIFA ------------------------------------
    auto *calsim = new R3BCalifa("califa_v2025.6.geo.root", {0., 0., -2.});
    calsim->SelectGeometryVersion(2025);
    run->AddModule(calsim);

    // -----   Create PrimaryGenerator   --------------------------------------
    TString fEventFile = "/nucl_lustre/g249/sim/p_F25_630.root";
    R3BINCLRootGenerator *gen =
        new R3BINCLRootGenerator((fEventFile).Data());
    gen->SetOnlyP2pSpallation(true);
    gen->SetXYZ(0, 0., 0.);
    gen->SetDxDyDz(0.1, 0.1, 2.5);

    auto primGen = new FairPrimaryGenerator();
    primGen->AddGenerator(gen);

    run->SetGenerator(primGen);
    run->SetStoreTraj(true);

    FairLogger::GetLogger()->SetLogVerbosityLevel("LOW");
    FairLogger::GetLogger()->SetLogScreenLevel("warn");

    // -----   Digitizer: Califa   --------------------------------------------
    auto califaDig = new R3BCalifaDigitizer();
    califaDig->SetNonUniformity(1.);
    califaDig->SetRealConfig(true);

    // califaDig->SetExpGammaEnergyRes(6.); // 5. means 5% at 1 MeV
    // califaDig->SetComponentRes(6.);
    // califaDig->SetDetectionThreshold(0.000010); // in GeV!! 0.000010 means 10 keV

    run->AddTask(califaDig);

    // -----   Clustering: Califa   -------------------------------------------
    // auto Cal2Clus = new R3BCalifaCrystalCal2Cluster();
    // Cal2Clus->SetCrystalThreshold(0.1);       // 100keV
    // Cal2Clus->SetProtonClusterThreshold(18.); // 12MeV
    // Cal2Clus->SetGammaClusterThreshold(0.5);  // 200keV
    // Cal2Clus->SelectGeometryVersion(2025);
    // Cal2Clus->SetRoundWindow(0.4);

    // R3BCalifaCrystalCal2Cluster ----------------------
    auto Cal2Clus = new R3BCalifaCrystalCal2Cluster();
    Cal2Clus->SetCrystalThreshold(500. / 1000.);         // 100keV
    Cal2Clus->SetProtonClusterThreshold(18000. / 1000.); // 12MeV
    Cal2Clus->SetGammaClusterThreshold(500. / 1000.);    // 200keV
    Cal2Clus->SelectGeometryVersion(2025);
    Cal2Clus->SetRoundWindow(0.4);
    Cal2Clus->SetRandomization(kTRUE);
    TString randomization_file = "/nucl_lustre/pablogrusell/g249/R3BParams_g249/califa/ang_dist_par.root";
    randomization_file.ReplaceAll("//", "/");
    Cal2Clus->SetRandomizationFile(randomization_file);

    run->AddTask(Cal2Clus);

    // -----   Initialize simulation run   ------------------------------------
    run->Init();
    TVirtualMC::GetMC()->SetRandom(new TRandom3(randomSeed));

    // ------  Increase nb of step for CALO
    Int_t nSteps = -15000;
    TVirtualMC::GetMC()->SetMaxNStep(nSteps);

    // -----   Runtime database   ---------------------------------------------
    Bool_t kParameterMerged = kTRUE;
    FairParRootFileIo *parOut = new FairParRootFileIo(kParameterMerged);
    parOut->open(parFile.Data());
    rtdb->setOutput(parOut);
    rtdb->saveOutput();
    rtdb->print();

    // -----   Start run   ----------------------------------------------------
    if (nEvents > 0)
    {
        run->Run(nEvents);
    }

    // -----   Finish   -------------------------------------------------------
    timer.Stop();
    Double_t rtime = timer.RealTime();
    Double_t ctime = timer.CpuTime();
    cout << endl
         << endl;
    cout << "Macro finished succesfully." << endl;
    cout << "Output file is " << outFile << endl;
    cout << "Parameter file is " << parFile << endl;
    cout << "Real time " << rtime << " s, CPU time " << ctime << "s" << endl
         << endl;

    cout << " Test passed" << endl;
    cout << " All ok " << endl;
}
