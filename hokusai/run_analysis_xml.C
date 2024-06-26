#include <glob.h>

std::vector<std::string> glob(const char *pattern) {
    glob_t g;
    glob(pattern, GLOB_TILDE, nullptr, &g); // one should ensure glob returns 0!
    std::vector<std::string> filelist;
    filelist.reserve(g.gl_pathc);
    for (size_t i = 0; i < g.gl_pathc; ++i) {
        filelist.emplace_back(g.gl_pathv[i]);
    }
    globfree(&g);
    return filelist;
} 

void run_analysis_xml(
    int fSystem = 108,
    int NSplit = 200,
    int splitID = 0,
    int version = 0,
    //TString fOutName = "test_108"
    TString fOutName = ""
    )
{
  //const std::string& xmlFile = Form("analysisConfig/analysisSn%d.xml",fSystem);
  TString xmlFile = Form("analysisConfig/analysisSn%d.xml",fSystem);
  if (version != 0)
    xmlFile = Form("analysisConfig/analysisSn%d_v%d.xml",fSystem,version);

  TString spiritroot = TString(gSystem -> Getenv("VMCWORKDIR"))+"/";
  TString fVersionOut; {
    TString name = spiritroot + "VERSION.compiled";
    std::ifstream vfile(name);
    vfile >> fVersionOut;
    vfile.close();
  }

  TString fVersionIn = "NewAna.2034.45b9400";
  TString fPathToDataIn = Form("/home/ejungwoo/data/trim/%s/Sn%d/",fVersionIn.Data(),fSystem);
  TString fPathToDataOut = Form("/home/ejungwoo/data/ana/%s/Sn%d/",fVersionOut.Data(),fSystem);
  gSystem -> Exec(TString("mkdir -p ")+fPathToDataOut);

  if (fOutName.IsNull())
    fOutName = Form("Sn%d",fSystem);

  //TString fPathToDataIn = "/data/Q20393/production/20191214/data/Sn"; fPathToDataIn = fPathToDataIn + fSystem + "/";
  //TString fVersionIn = "develop.1964.781a3cf";

  /***************************************************************
  * Read xml
  ****************************************************************/
  TDOMParser parser;
  parser.SetValidate(false);
  parser.ParseFile(xmlFile.Data());
  auto node = parser.GetXMLDocument()->GetRootNode()->GetChildren();
  TXMLNode *TaskNode = nullptr;
  //TChain chain("spirit");
  TChain chain("cbmsim");

  for(auto child = node; child; child = child -> GetNextNode())
    if(child -> GetNodeType() == TXMLNode::kXMLElementNode)
    {
      if(std::strcmp(child -> GetNodeName(), "TaskList") == 0) TaskNode = child;
      if(std::strcmp(child -> GetNodeName(), "IOInfo") == 0)
      {
        std::string dataType(static_cast<TXMLAttr*>(child -> GetAttributes() -> At(0)) -> GetValue());
        //for(auto IOInfo = child -> GetChildren(); IOInfo; IOInfo = IOInfo -> GetNextNode())
          //if(std::strcmp(IOInfo -> GetNodeName(), "DataDir") == 0) fPathToDataIn = IOInfo -> GetText(); 

        if(dataType == "Real")
        {
          int start_run, last_run;
          for(auto IOInfo = child -> GetChildren(); IOInfo; IOInfo = IOInfo -> GetNextNode())
          {
            if(std::strcmp(IOInfo -> GetNodeName(), "RunFirst") == 0) start_run = std::atoi(IOInfo -> GetText());
            if(std::strcmp(IOInfo -> GetNodeName(), "RunLast") == 0) last_run = std::atoi(IOInfo -> GetText());
          }
          auto origLevel = gErrorIgnoreLevel;
          gErrorIgnoreLevel = kFatal;
          for(int runNo = start_run; runNo <= last_run; ++runNo)
          {
            //auto fileName = Form("%srun%d_s*.reco.%s.conc.root",fPathToDataIn.Data(),runNo,fVersionIn.Data());
            auto fileName = Form("%srun%d_s*.reco.%s.conc.trimmed.root",fPathToDataIn.Data(),runNo,fVersionIn.Data());
            cout << fileName << endl;
            chain.Add(fileName);
            //TString sRunNo = TString::Itoa(runNo, 10);
            //chain.Add(fPathToDataIn+"run"+sRunNo+"_s*.reco.*trimmed*.root/cbmsim");
            //chain.Add(fPathToDataIn+"run"+sRunNo+"_s*.reco.*conc*.root/spirit");
            //chain.Add(fPathToDataIn+"run"+sRunNo+"_s*.reco.*conc*.root/cbmsim");
            //std::cout << "Reading from file " << fPathToDataIn+"run"+sRunNo+"_s*.reco.*root" << std::endl;
          }
          if(chain.GetEntries() == 0) throw std::runtime_error("No entries is being read from the file!");
          gErrorIgnoreLevel = origLevel;
        }
        /*
        else if(dataType == "Sim")
        {
          TString inputName;
          for(auto IOInfo = child -> GetChildren(); IOInfo; IOInfo = IOInfo -> GetNextNode())
            if(std::strcmp(IOInfo -> GetNodeName(), "InputName") == 0) inputName = IOInfo -> GetText();
          for(const auto &filename : glob(fPathToDataIn+inputName+".conc.root"))
            chain.Add(filename.c_str());
        }
        */
      }
    }
  STAnalysisFactory factory(TaskNode);

  /***************************************************************
  *  FairRoot setup
  ****************************************************************/

  TString par = spiritroot+"parameters/ST.parameters.par";
  TString geo = spiritroot+"geometry/geomSpiRIT.man.root";
  TString out = fPathToDataOut+fOutName+"_" + TString::Itoa(splitID, 10) + "_ana"+"_v"+version+"."+fVersionOut+".root";
  TString log = fPathToDataOut+fOutName+"_" + TString::Itoa(splitID, 10) + "_ana"+"_v"+version+"."+fVersionOut+".log";

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

  int nentries = chain.GetEntries();
  int eventsInSplit = nentries / NSplit + 1;
  int eventsStart = splitID*eventsInSplit;
  if(eventsStart + eventsInSplit > nentries) eventsInSplit = nentries - eventsStart;

  auto reader = new STConcReaderTask();
  reader -> SetChain(&chain);
  reader -> SetEventID(eventsStart);
  run -> AddTask(reader);

  std::vector<FairTask*> tasks;
  tasks.push_back(factory.GetFilterEventTask());
  tasks.push_back(factory.GetPIDTask());
  tasks.push_back(factory.GetTransformFrameTask());
  tasks.push_back(factory.GetEfficiencyTask());

  for(auto task : tasks)
    if(task) run -> AddTask(task);

  run -> Init();
  run -> Run(0, eventsInSplit);

  cout << "Log    : " << log << endl;
  cout << "Output : " << out << endl;
}
