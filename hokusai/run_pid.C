void run_pid
(
  Int_t fRunId = 2898,
  Int_t fSystem = 132
)
{
  TString nameFittedPIDData = Form("pid_prob_summary.sys%d.root",fSystem);
  TString tag = Form("sn%d",fSystem);

  TString parameterFile = "ST.parameters.fullmc.par";
  TString pathToDataConc = "/data/Q20393/recodata/develop_1944_33821f0_conc/";
  //TString pathToDataOutput = "/data/Q19393/recodata/develop_1944_33821f0_pid/";
  TString pathToDataOutput = "/home/ejungwoo/data/pid/";
  TString pathToBeam = "/data/Q20393/Frozen/Aug2019/beam/";
  TString spiritroot = TString(gSystem -> Getenv("VMCWORKDIR"))+"/";
  TString pathToSpiritInput = spiritroot + "/input/";
  TString pathToBeamCut = pathToSpiritInput;
  TString versionConc = "reco.develop.1944.33821f0";
  TString versionOutput; {
    TString name = spiritroot + "VERSION.compiled";
    std::ifstream vfile(name);
    vfile >> versionOutput;
    vfile.close();
  }
  versionOutput = versionOutput + "." + tag;

  TString runName = Form("run%d",fRunId);
  TString nameConc = runName + "." + versionConc;
  TString nameOutput = runName;

  TString par = spiritroot+"parameters/"+parameterFile;
  TString geo = spiritroot+"geometry/geomSpiRIT.man.root"; 
  TString beam = pathToBeam+Form("beam_run%d.ridf.root",fRunId);
  TString beamCut = pathToBeamCut + Form("Sn%d-beamCut.root",fSystem);
  TString conc = pathToDataConc+nameConc+".conc.root"; 
  TString output = pathToDataOutput+nameOutput+"."+versionOutput+".root";
  TString log = pathToDataOutput+nameOutput+"."+versionOutput+".log";
  nameFittedPIDData = pathToSpiritInput + nameFittedPIDData;

  FairLogger *logger = FairLogger::GetLogger();
  logger -> SetLogToScreen(true);

  FairParAsciiFileIo* parReader = new FairParAsciiFileIo();
  parReader -> open(par);

  TChain chain("spirit");
  chain.Add(conc);
 
  FairRunAna* run = new FairRunAna();
  run -> SetGeomFile(geo);
  run -> SetOutputFile(output);
  run -> GetRuntimeDb() -> setSecondInput(parReader);

  auto reader = new STConcReaderTask();
  reader -> SetPersistence(0);
  reader -> SetChain(&chain);

  auto pid = new STPIDProbTask();
  pid -> SetSaveSeparatePIDs(true);
  //pid -> UseG4CutNumCluster(0.4);
  pid -> SetOutputPath(pathToDataOutput);
  pid -> SetPIDData(nameFittedPIDData);
  pid -> SetRIDF(beam);
  pid -> SetBeamA(fSystem); 
  pid -> SetBeamCut(beamCut);

  run -> AddTask(reader);
  run -> AddTask(pid);

  run -> Init();
  run -> Run(0, chain.GetEntries());

  cout << "Conc   : " << conc << endl;
  cout << "PID 0  : " << pid -> GetOutputFileName(0) << endl;
  cout << "PID 1  : " << pid -> GetOutputFileName(1) << endl;
  cout << "PID 2  : " << pid -> GetOutputFileName(2) << endl;
  cout << "PID 3  : " << pid -> GetOutputFileName(3) << endl;
  cout << "PID 4  : " << pid -> GetOutputFileName(4) << endl;
  cout << "PID 5  : " << pid -> GetOutputFileName(5) << endl;
  cout << "PID 5  : " << pid -> GetOutputFileName(6) << endl;
  cout << "PID 5  : " << pid -> GetOutputFileName(7) << endl;

  gApplication -> Terminate();
}
