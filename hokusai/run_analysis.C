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

void run_analysis_core(TString par, TString geo, TString out,
                       TString log, TChain& chain, int targetMass,
                       TString meta, TString pidFile, bool iterateMeta, 
                       EfficiencyFactory *effFactory, bool targetELoss=true)
{

  FairLogger *logger = FairLogger::GetLogger();
  logger -> SetLogToScreen(true);

  FairParAsciiFileIo* parReader = new FairParAsciiFileIo();
  parReader -> open(par);

  FairRunAna* run = new FairRunAna();
  run -> SetGeomFile(geo);
  run -> SetOutputFile(out);
  run -> GetRuntimeDb() -> setSecondInput(parReader);

  auto reader = new STConcReaderTask();
  reader -> SetChain(&chain);

  auto eventFilter = new STFilterEventTask();
  //eventFilter -> SetBeamCut("BeamCut.root", "Sn108");
  if (targetMass==112)
    eventFilter -> SetBeamCut("inputFiles/Sn108-beamCut.root", "sigma30");
  else
    eventFilter -> SetBeamCut("inputFiles/Sn132-beamCut.root", "sigma30");
  eventFilter -> SetVertexCut(-18.480, -11.165);
  eventFilter -> SetMultiplicityCut(54, 100);

  //auto filter = new STFilterTask();
  //filter -> SetThetaCut(10);

  //auto PIDCut = new STPIDCutTask();
  //PIDCut -> SetNPitch(1);
  //PIDCut -> SetNYaw(1);
  //PIDCut -> SetCutConditions(15, 15);
  //PIDCut -> SetCutFile("PIDCut_SimSn108.root");

  //auto unfold = new STUnfoldingTask();
  //unfold -> SetMomBins(100, 3700, 50, 25);
  //unfold -> SetThetaBins(0, 1.57, 30, 25);
  //unfold -> LoadMCData("cbmsim", "data/embed_dump/ImQMD_embedCorTriton2/*");
  //unfold -> SetCutConditions(15, 20);
  
  //auto PIDProb = new STPIDMachineLearningTask();
  //PIDProb -> SetChain(&chain);
  //PIDProb -> SetBufferSize(25000);
  //PIDProb -> SetModel("MLForestCut", STAlgorithms::RandomForest);

  auto PIDProb = new STPIDProbTask();
  PIDProb -> SetMetaFile(meta.Data(), iterateMeta);
  PIDProb -> SetPIDFitFile(pidFile.Data());

  auto transform = new STTransformFrameTask();
  transform -> SetDoRotation(true);
  transform -> SetPersistence(true);
  transform -> SetTargetMass(targetMass);
  if(targetELoss)
  {
    transform -> SetTargetThickness(0.8);
    //transform -> SetEnergyLossFile((targetMass == 112)? "../parameters/Sn108Sn112.txt": "../parameters/Sn132Sn124.txt");
    transform -> SetEnergyLossFile((targetMass == 112)? "inputFiles/Sn108Sn112.txt": "inputFiles/Sn132Sn124.txt");
  }

  auto efficiency = new STEfficiencyTask(effFactory);
  for(int pdg : STAnaParticleDB::SupportedPDG)
  {
    auto& settings = efficiency -> AccessSettings(pdg);

    settings.NClusters = 15;
    settings.DPoca = 15;
    //settings.PhiCuts = {{160, 220}, {0, 20}, {320, 360}};
    settings.PhiCuts = {{-180, 140}, {-40, 40}, {140, 180}};

    /*
    settings.NClusters = 10;
    settings.DPoca = 10;
    settings.PhiCuts = {{0,360}};
    */

    settings.ThetaMin = 0; settings.ThetaMax = 90; settings.NThetaBins = 15;
    auto &mBins = settings.NMomBins; auto &mMin = settings.MomMin; auto &mMax = settings.MomMax;
    auto &ptBins = settings.NPtBins; auto &ptMin = settings.PtMin; auto &ptMax = settings.PtMax;
    auto &CMzBins = settings.NCMzBins; auto &CMzMin = settings.CMzMin; auto &CMzMax = settings.CMzMax;
    mBins = 15; ptBins = 20; CMzBins = 20;
    if(pdg == 2212)      { mMin = 100;  mMax = 1500; ptMin = 0; ptMax = 1300; CMzMin = -1000; CMzMax = 1000; }
    if(pdg == 1000010020){ mMin = 200;  mMax = 2200; ptMin = 0; ptMax = 2000; CMzMin = -1300; CMzMax = 1300; }
    if(pdg == 1000010030){ mMin = 500;  mMax = 4500; ptMin = 0; ptMax = 2500; CMzMin = -2500; CMzMax = 2500; }
    if(pdg == 1000020030){ mMin = 500;  mMax = 3200; ptMin = 0; ptMax = 1900; CMzMin = -1900; CMzMax = 1900; }
    if(pdg == 1000020040){ mMin = 1000; mMax = 4200; ptMin = 0; ptMax = 2500; CMzMin = -2500; CMzMax = 2500; }
    if(pdg == 1000020060){ mMin = 2000; mMax = 6000; ptMin = 0; ptMax = 1300; CMzMin = -1000; CMzMax = 1000; }
  }

  run -> AddTask(reader);
  run -> AddTask(eventFilter);
  //run -> AddTask(filter);
  //run -> AddTask(PIDCut);
  run -> AddTask(PIDProb);
  //run -> AddTask(unfold);

  if(!iterateMeta)
  {
    run -> AddTask(transform);
    run -> AddTask(efficiency);
  }

  auto simpleGraphs = new STSimpleGraphsTask();
  run -> AddTask(simpleGraphs);

  run -> Init();
  run -> Run(0, chain.GetEntries());

  cout << "Log    : " << log << endl;
  cout << "Output : " << out << endl;

  gApplication -> Terminate();

}

void run_analysis
(
  int fStartRunNo = 2841,
  int fEndRunNo = 2841,
  int fBeamMass = 132,
  int fTargetMass = 124,

  //TString fPathToOut = "/home/ejungwoo/data/pid2/",
  TString fPathToOut = "/home/ejungwoo/data/pid3/",

  //TString fPathToInput = "/home/ejungwoo/data/conc/system132/",
  //TString fVersionIn = "develop.1944.33821f0",

  TString fPathToInput = "/data/Q20393/production/20191214/data/",//Sn132/",
  TString fVersionIn = "develop.1964.781a3cf",

  TString fMetaFile = "inputFiles/Meta_Sn132KanekoMult50.root", // XXX
  TString fPidFile = "inputFiles/PIDSigma_Sn132KanekoMult50.root", // XXX
  bool fIterateMeta = false
)
{
  TString spiritroot = TString(gSystem -> Getenv("VMCWORKDIR"))+"/";
  if (fPathToInput.IsNull())
    fPathToInput = spiritroot+"macros/data/";
  TString version; {
    TString name = spiritroot + "VERSION.compiled";
    std::ifstream vfile(name);
    vfile >> version;
    vfile.close();
  }
  fPathToInput = fPathToInput + Form("/Sn%d/",fBeamMass);

  TString par = spiritroot+"parameters/ST.parameters.par";
  TString geo = spiritroot+"geometry/geomSpiRIT.man.root"; 

  TString outName = fPathToOut + Form("beam%d_run%d_%d",fBeamMass,fStartRunNo,fEndRunNo);

  TString log = outName + "_ana.log";
  TString out = outName + "_ana.root";

  TChain chain("spirit");
  for (int runNo = fStartRunNo; runNo <= fEndRunNo; ++runNo) {
    TString sRunNo = TString::Itoa(runNo, 10);
    //std::cout << "Reading from file " << in.back() << std::endl;
    chain.Add(fPathToInput+"run"+sRunNo+Form("_s*.reco.%s.conc.root",fVersionIn.Data()));
  }

  auto effFactory = new EfficiencyInCMFactory();
  effFactory -> SetDataBaseForPDG(2212,       "inputFiles/Run2899KanekoNoSC_embedNewProton.root");
  effFactory -> SetDataBaseForPDG(1000010020, "inputFiles/Run2899KanekoNoSC_embedNewDeuteron.root");
  effFactory -> SetDataBaseForPDG(1000010030, "inputFiles/Run2899KanekoNoSC_embedNewTriton.root");
  effFactory -> SetDataBaseForPDG(1000020030, "inputFiles/Run2899KanekoNoSC_embedNewHe3.root");
  effFactory -> SetDataBaseForPDG(1000020040, "inputFiles/Run2899KanekoNoSC_embedNewHe4.root");
  effFactory -> SetDataBaseForPDG(1000020060, "inputFiles/Run2899KanekoNoSC_embedNewHe6.root");


  run_analysis_core(par, geo, out, log, chain, fTargetMass, fMetaFile, fPidFile, fIterateMeta, effFactory);
}
