void run_analysis_xml(int system, TString anaName="", TString fOutName="", bool iter_unfold=false)
{
  //Read xml
  TDOMParser parser;
  parser.SetValidate(false);
  TString xmlFile = Form("analysisInputFiles/analysisConfig/analysisSn%dCM%s.xml", system, anaName.Data());
  parser.ParseFile(xmlFile.Data());
  auto node = parser.GetXMLDocument()->GetRootNode()->GetChildren();
  TXMLNode *TaskNode = nullptr;
  TXMLNode *IONode = nullptr;

  for(auto child = node; child; child = child -> GetNextNode()) {
    TString nodename = child -> GetNodeName();
    if(child -> GetNodeType() == TXMLNode::kXMLElementNode)
    {
      if(std::strcmp(child -> GetNodeName(), "TaskList") == 0) TaskNode = child;
      if(std::strcmp(child -> GetNodeName(), "IOInfo") == 0) IONode = child;
    }
  }

  // FairRoot setup
  auto reader = new STConcReaderTask();
  TString fPathToData = reader -> LoadFromXMLNode(IONode).c_str();
  TString spiritroot = TString(gSystem -> Getenv("VMCWORKDIR"))+"/";
  TString fVersionOut;
  std::ifstream(spiritroot+"VERSION.compiled") >> fVersionOut;
  TString fPathToDataOut = Form("/home/ejungwoo/data/ana/%s/Sn%d/",fVersionOut.Data(),system);
  gSystem -> Exec(TString("mkdir -p ")+fPathToDataOut);
  if (fOutName.IsNull())
    fOutName = TString("sys")+system;

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

  // Add tasks
  STAnalysisFactory factory(TaskNode);
  int nentries = reader -> GetNEntries();

  std::vector<FairTask*> tasks;
  tasks.push_back(reader);
  tasks.push_back(factory.GetFilterEventTask());
  tasks.push_back(factory.GetPIDTask());
  tasks.push_back(factory.GetTransformFrameTask());
  auto eff = static_cast<STEfficiencyTask*>(factory.GetEfficiencyTask());
  if(eff) eff -> UpdateUnfoldingFile(iter_unfold);
  tasks.push_back(eff);
  tasks.push_back(factory.GetSimpleGraphsTask());
  tasks.push_back(factory.GetERATTask());

  for(auto task : tasks)
    if(task) run -> AddTask(task);
  auto pst = new STParticleSummaryTask();
  pst -> SetCuts(0.2,0.05,5);
  run -> AddTask(pst);

  run -> Init();
  run -> Run(0, nentries);

  cout /*<< "Log     : "*/ << log << endl;
  cout /*<< "Output  : "*/ << out << endl;
  cout /*<< "Summary : "*/ << pst -> GetSummaryName() << endl;

  //gApplication -> Terminate();
}
