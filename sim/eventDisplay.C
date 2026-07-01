void eventDisplay()
{
  FairRunAna *fRun = new FairRunAna();
  fRun->SetSource(new FairFileSource("./califa/sim_1.root"));
  fRun->SetSink(new FairRootFileSink("./califa/vis.root"));

  FairRuntimeDb *rtdb = fRun->GetRuntimeDb();
  FairParRootFileIo *parIo1 = new FairParRootFileIo();
  parIo1->open("./califa/par_1.root");
  rtdb->setFirstInput(parIo1);
  rtdb->print();

  R3BEventManager *fMan = new R3BEventManager();
  R3BMCTracks *Track = new R3BMCTracks("Monte-Carlo Tracks");

  fMan->AddTask(Track);
  fMan->Init();
  gEve->GetDefaultGLViewer()->SetClearColor(kOrange - 4);
}
