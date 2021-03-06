void exp1803_reco(Int_t nEvents = 100) {
  //---------------------Files-----------------------------------------------
  TString inFile = "sim_digi.root";
  TString outFile = "reco.root";
  TString parFile = "par.root";
  TString parOutFile = "parOut.root";
  // -----   Timer   --------------------------------------------------------
  TStopwatch timer;
  timer.Start();  
  // -----   Digitization run   ---------------------------------------------
  FairRunAna *run= new FairRunAna();
  run->SetInputFile(inFile);
  run->SetOutputFile(outFile);
  // ------------------------------------------------------------------------
  //-------- Set MC event header --------------------------------------------
  EREventHeader* header = new EREventHeader();
  run->SetEventHeader(header);
  // ------------------------BeamDetTrackFinger--------------------------------
  Int_t verbose = 1; // 1 - only standard log print, 2 - print digi information 
  ERBeamDetTrackFinder* trackFinder = new ERBeamDetTrackFinder(verbose);
  run->AddTask(trackFinder);
  // -----------------------BeamDetTrackPID----------------------------------
  Int_t Z = 2, A = 6, Q = 2;
  TString ionName = "6He";
  ERBeamDetPID* pid = new ERBeamDetPID(verbose);
  pid->SetIonMassNumber(A);
  pid->SetBoxPID(0., 1000., 0., 1000.);
  pid->SetOffsetToF(0.);
  pid->SetProbabilityThreshold(0);
  run->AddTask(pid);  
  // ------- QTelescope TrackFinder -------------------------------------------
  ERQTelescopeTrackFinder* qtelescopeTrackFinder = new ERQTelescopeTrackFinder(verbose);
  qtelescopeTrackFinder->SetHitStation("DoubleSi_SD2_XY_1");
  qtelescopeTrackFinder->SetHitStation("DoubleSi_SD2_XY_3");
  qtelescopeTrackFinder->SetStripEdepRange(0., 100.);          // [GeV]
  qtelescopeTrackFinder->SetTargetPoint(0., 0., 0.);
  // qtelescopeTrackFinder->SetStripEdepRange(0.0097, 100.);   // [GeV]
  // qtelescopeTrackFinder->SetEdepDiffXY(5.);                 // [GeV]
  qtelescopeTrackFinder->SetEdepMaxDiffXY(0.5); 
  run->AddTask(qtelescopeTrackFinder);  
  // -----------Runtime DataBase info ---------------------------------------
  FairRuntimeDb* rtdb = run->GetRuntimeDb();
  FairParRootFileIo*  parInput = new FairParRootFileIo();
  parInput->open(parFile.Data(), "UPDATE");
  rtdb->setFirstInput(parInput);
  // -----   Intialise and run   --------------------------------------------
  FairLogger::GetLogger()->SetLogScreenLevel("DEBUG");
  
  run->Init();
  run->Run(0, nEvents);
  // ------------------------------------------------------------------------;
  rtdb->setOutput(parInput);
  rtdb->saveOutput();
  // -----   Finish   -------------------------------------------------------
  timer.Stop();
  Double_t rtime = timer.RealTime();
  Double_t ctime = timer.CpuTime();
  cout << endl << endl;
  cout << "Macro finished succesfully." << endl;
  cout << "Output file writen:  "    << outFile << endl;
  cout << "Parameter file writen " << parFile << endl;
  cout << "Real time " << rtime << " s, CPU time " << ctime << " s" << endl;
  cout << endl;
  // ------------------------------------------------------------------------

}
