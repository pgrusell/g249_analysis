void eventDisplay()
{
  FairRunAna *fRun = new FairRunAna();
  //fRun->SetSource(new FairFileSource("/nucl_lustre/pablogrusell/g249/g249_analysis/final/macros/erel/sim_results/bg.simu.root"));
  
  fRun->SetSource(new FairFileSource("califa/sim_1.root"));
  
  fRun->SetSink(new FairRootFileSink("vis.root"));

  FairRuntimeDb *rtdb = fRun->GetRuntimeDb();
  FairParRootFileIo *parIo1 = new FairParRootFileIo();
  parIo1->open("califa/par_1.root");
  //parIo1->open("/nucl_lustre/pablogrusell/g249/g249_analysis/final/macros/erel/sim_results/bg.para.root");
  rtdb->setFirstInput(parIo1);
  rtdb->print();

  R3BEventManager *fMan = new R3BEventManager();
  R3BMCTracks *Track = new R3BMCTracks("Monte-Carlo Tracks");

  fMan->AddTask(Track);
  fMan->Init();
  gEve->GetDefaultGLViewer()->SetClearColor(kOrange - 4);
}
