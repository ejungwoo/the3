void run_analysis_xml(const std::string& xmlFile="analysisConfig/analysisSn108CM.xml", TString fOutName="test", bool iter_unfold=false)
{
  /***************************************************************
  * Read xml
  ****************************************************************/
  TDOMParser parser;
  parser.SetValidate(false);
  parser.ParseFile(xmlFile.c_str());
  auto node = parser.GetXMLDocument()->GetRootNode()->GetChildren();
  TXMLNode *TaskNode = nullptr;
  TXMLNode *IONode = nullptr;

  for(auto child = node; child; child = child -> GetNextNode()) {
    TString nodename = child -> GetNodeName();
    /*
    if (nodename == "IOInfo") {
      for(auto child2 = child -> GetChildren(); child2; child2 = child2 -> GetNextNode()) {
        TString nodename2 = child2 -> GetNodeName();
        if (nodename2 == "DataDir") {
          cout << child2 -> GetText() << endl;
        }
        for(auto child3 = child2 -> GetChildren(); child3; child3 = child3 -> GetNextNode()) {
          if (TString(child3 -> GetNodeName())!="text") {
            cout << child3 -> GetNodeName() << " " << child3 -> GetText() << endl;
          }
        }
      }
    }
    */
    if(child -> GetNodeType() == TXMLNode::kXMLElementNode)
    {
      if(std::strcmp(child -> GetNodeName(), "TaskList") == 0) TaskNode = child;
      if(std::strcmp(child -> GetNodeName(), "IOInfo") == 0) IONode = child;
    }
  }

  STAnalysisFactory factory(TaskNode);

  // FairRoot setup
  auto reader = new STConcReaderTask();
  TString spiritroot = TString(gSystem -> Getenv("VMCWORKDIR"))+"/";
  TString fVersionIn = "NewAna.2034.45b9400";
  TString fPathToDataIn = Form("/home/ejungwoo/data/trim/%s/Sn%d/",fVersionIn.Data(),fSystem);
  TString fVersionOut;
  std::ifstream(spiritroot+"VERSION.compiled") >> fVersionOut;
  TString fPathToDataOut = Form("/home/ejungwoo/data/ana/%s/Sn%d/",fVersionOut.Data(),fSystem);
  gSystem -> Exec(TString("mkdir -p ")+fPathToDataOut);

  TString par = spiritroot + "parameters/ST.parameters.par";
  TString geo = spiritroot + "geometry/geomSpiRIT.man.root";
  TString out = fPathToDataOut + fOutName + "." + fVersionOut+"._ana.root";
  TString log = fPathToDataOut + fOutName + "." + fVersionOut+"._ana.log";

  FairLogger *logger = FairLogger::GetLogger();
  logger -> SetLogToScreen(true);

  FairParAsciiFileIo* parReader = new FairParAsciiFileIo();
  parReader -> open(par);

  FairRunAna* run = new FairRunAna();
  run -> SetGeomFile(geo);
  run -> SetOutputFile(out);
  run -> GetRuntimeDb() -> setSecondInput(parReader);

  // Add tasks
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

  run -> Init();
  run -> Run(0, nentries);

  cout << "Log    : " << log << endl;
  cout << "Output : " << out << endl;
}
