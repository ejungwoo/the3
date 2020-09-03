void run_analysis_xml(int system, TString anaName="", TString fOutName="", bool iter_unfold=false)
{
  std::srand(std::time(0));

  /***************************************************************
  * Read xml
  ****************************************************************/
  TDOMParser parser;
  parser.SetValidate(false);
  TString xmlFile = Form("analysisInputFiles/analysisConfig/analysisSn%dCM%s.xml", system, anaName.Data());
  parser.ParseFile(xmlFile.Data());
  auto node = parser.GetXMLDocument()->GetRootNode()->GetChildren();
  STAnalysisFactory factory(node);

  /***************************************************************
  *  FairRoot setup
  ****************************************************************/
  auto reader = factory.GetReaderTask();
  int nentries = reader -> GetNEntries();
  TString fPathToData = reader -> GetPathToData();
  TString spiritroot = TString(gSystem -> Getenv("VMCWORKDIR"))+"/";
  TString fVersionOut; std::ifstream(spiritroot+"VERSION.compiled") >> fVersionOut;
  TString fPathToDataOut = Form("/home/ejungwoo/data/ana/%s/Sn%d/",fVersionOut.Data(),system);
  gSystem -> Exec(TString("mkdir -p ")+fPathToDataOut);
  if (fOutName.IsNull())
    fOutName = TString("sys")+system+anaName;

  TString par = spiritroot + "parameters/ST.parameters.par";
  TString geo = spiritroot + "geometry/geomSpiRIT.man.root";
  TString out = fPathToDataOut + fOutName + "." + fVersionOut+".ana.root";
  TString log = fPathToDataOut + fOutName + "." + fVersionOut+".ana.log";

  FairLogger *logger = FairLogger::GetLogger();
  logger -> SetLogToScreen(true);

  FairParAsciiFileIo* parReader = new FairParAsciiFileIo();
  parReader -> open(par);

  FairRunAna* run = new FairRunAna();
  run -> SetGeomFile(geo);
  run -> SetOutputFile(out);
  run -> GetRuntimeDb() -> setSecondInput(parReader);

  /*************************************************************
  * Add tasks
  **************************************************************/

  std::vector<FairTask*> tasks;
  tasks.push_back(reader);
  tasks.push_back(factory.GetFilterEventTask());
//tasks.push_back(factory.GetDivideEventTask());
  tasks.push_back(factory.GetPIDTask());
//tasks.push_back(factory.GetPiProbTask());
  tasks.push_back(factory.GetTransformFrameTask());
  auto eff = factory.GetEfficiencyTask();
  if(eff) eff -> UpdateUnfoldingFile(iter_unfold);
  tasks.push_back(eff);
  tasks.push_back(factory.GetERATTask());
  tasks.push_back(factory.GetReactionPlaneTask());
  tasks.push_back(factory.GetSimpleGraphsTask());

  for(auto task : tasks)
    if(task) run -> AddTask(task);

  auto pst = new STParticleSummaryTask();
  pst -> SetCuts(0.01,0.01,5);
  run -> AddTask(pst);

  run -> Init();
  run -> Run(0, nentries);

  cout /*<< "Log     : "*/ << log << endl;
  cout /*<< "Output  : "*/ << out << endl;
  cout /*<< "Summary : "*/ << pst -> GetSummaryName() << endl;

  //gApplication -> Terminate();
}
