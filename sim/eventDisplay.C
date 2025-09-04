void eventDisplay()
{
  FairRunAna *fRun = new FairRunAna();
  fRun->SetSource(new FairFileSource("/nucl_lustre/pablogrusell/g249/g249_analysis/results/sim/sim_neutronDecay.root"));
  fRun->SetSink(new FairRootFileSink("/nucl_lustre/pablogrusell/g249/g249_analysis/results/sim/par_neutronDecay.root"));

  FairRuntimeDb *rtdb = fRun->GetRuntimeDb();
  FairParRootFileIo *parIo1 = new FairParRootFileIo();
  parIo1->open("/nucl_lustre/pablogrusell/g249/g249_analysis/results/sim/par_neutronDecay.root");
  rtdb->setFirstInput(parIo1);
  rtdb->print();

  R3BEventManager *fMan = new R3BEventManager();
  R3BMCTracks *Track = new R3BMCTracks("Monte-Carlo Tracks");

  fMan->AddTask(Track);
  fMan->Init();
  gEve->GetDefaultGLViewer()->SetClearColor(kOrange - 4);
}
