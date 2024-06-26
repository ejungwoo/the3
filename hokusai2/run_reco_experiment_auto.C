#include <unistd.h>

void readEventList(TString eventListFile, map<Int_t, vector<Int_t> *> &events);

void run_reco_experiment_auto
(
  Int_t fRunNo = 2841,
  Int_t fSplitNo = 0,
  Int_t fNumEventsInSplit = 10,
  Int_t fSys = 132,
  TString fPathToData = "/home/ejungwoo/data/reco/",
  Bool_t fSelectPiEvent = false,
  TString fOut = "",
  TString fLog = "",
  TString fConc = "",
  std::vector<Int_t> fSkipEventArray = {},
  TString fMCFile = "",
  TString fSupplePath = "/data/Q18393/rawdataSupplement",
  TString ggDataPathWithFormat = "/data/Q18393/production/ggNoise/ggNoise_%d.root",
  TString beamDataPathWithFormat = "/data/Q20393/production/20200529/SpiRITROOT/beam/beam_run%d.ridf.root"
)
{
  double fLength = 300;
  /* ======= This part you need initial configuration ========= */
  // Parameter database file - files should be in parameters folder.
  TString systemDB = "systemDB.csv";
  TString runDB = "runDB.csv";

  // Data paths - must have one %d for run number
  // If you don't need either of them, pass it blank.
  //TString ggDataPathWithFormat = "/mnt/spirit/rawdata/misc/Frozen_Information_For_SpiRIT_Analysis/Aug2019/ggNoise/ggNoise_%d.root";
  //TString beamDataPathWithFormat = "/mnt/spirit/rawdata/misc/Frozen_Information_For_SpiRIT_Analysis/Aug2019/BeamData/beam/beam_run%d.ridf.root";

  // Meta data path
  Bool_t fUseMeta = kTRUE;
  //TString fSupplePath = "/mnt/spirit/rawdata/misc/rawdataSupplement";

  // Use relative gain matching - Amplify low gain section
  Bool_t fUseGainMatching = kTRUE;
  /* ========================================================== */

  TString fSpiRITROOTPath = TString(gSystem -> Getenv("VMCWORKDIR"))+"/";
  TString fSystemDB = fSpiRITROOTPath + "parameters/" + systemDB;
  TString fRunDB = fSpiRITROOTPath + "parameters/" + runDB;

  auto fParamSetter = new STParameters(fRunNo, fSystemDB, fRunDB);
  auto fNumEventsInRun = fParamSetter -> GetNumTotalEvents();
  auto fParameterFile = fParamSetter -> GetParameterFile();
  auto fYPedestalOffset = fParamSetter -> GetYPedestal();
  auto fFieldOffsetX = fParamSetter -> GetFieldOffsetX();
  auto fFieldOffsetY = fParamSetter -> GetFieldOffsetY();
  auto fFieldOffsetZ = fParamSetter -> GetFieldOffsetZ();
  auto fSheetChargeDensity = fParamSetter -> GetSheetChargeDensity();
  auto fSystemID = fParamSetter -> GetSystemID();
  auto fTargetZ = fParamSetter -> GetTargetZ();
  auto fBDCOffsetX = fParamSetter -> GetBDCOffsetX();
  auto fBDCOffsetY = fParamSetter -> GetBDCOffsetY();
  auto fGGRunID = fParamSetter -> GetGGRunID();
  auto fRelativeGainRunID = fParamSetter -> GetRelativeGainRunID();

  auto fIsGGDataSet = !ggDataPathWithFormat.IsNull();
  auto fIsBeamDataSet = !beamDataPathWithFormat.IsNull();
  TString fGGData = "";
  TString fBeamData = "";
  if (fIsGGDataSet)   fGGData = Form(ggDataPathWithFormat.Data(), fGGRunID);
  if (fIsBeamDataSet) fBeamData = Form(beamDataPathWithFormat.Data(), fRunNo);
  TString fGainMatchingFile = fSpiRITROOTPath + Form("parameters/RelativeGainRun%d.list", fRelativeGainRunID);

  Int_t start = fSplitNo * fNumEventsInSplit;
  if (start >= fNumEventsInRun) return;
  if (start + fNumEventsInSplit > fNumEventsInRun)
    fNumEventsInSplit = fNumEventsInRun - start;

  TString sRunNo   = TString::Itoa(fRunNo, 10);
  TString sSplitNo = TString::Itoa(fSplitNo, 10);

  if (fPathToData.IsNull())
    fPathToData = fSpiRITROOTPath+"macros/data/";
  else 
    gSystem->mkdir(fPathToData.Data());
  
  TString version; {
    TString name = fSpiRITROOTPath + "VERSION.compiled";
    std::ifstream vfile(name);
    vfile >> version;
    vfile.close();
  }

  TString par = fSpiRITROOTPath+"parameters/"+fParameterFile;
  TString geo = fSpiRITROOTPath+"geometry/geomSpiRIT.man.root";
  TString fRawDataList = TString(gSystem -> Getenv("PWD"))+"/list_run"+sRunNo+".txt";

  if (fOut.IsNull())  fOut  = fPathToData+"run"+sRunNo+"_s"+sSplitNo+".reco."+version+".root";
  if (fLog.IsNull())  fLog  = fPathToData+"run"+sRunNo+"_s"+sSplitNo+".reco."+version+".log";
  if (fConc.IsNull()) fConc = fPathToData+"run"+sRunNo+"_s"+sSplitNo+".reco."+version+".conc.root";

  //TString fOut = fPathToData+"run"+sRunNo+"_s"+sSplitNo+".reco."+version+".root";
  //TString fLog = fPathToData+"run"+sRunNo+"_s"+sSplitNo+"."+version+".log";
  //TString fConc = fPathToData+"run"+sRunNo+"_s"+sSplitNo+".reco."+version+".conc.root";
  
  if (TString(gSystem -> Which(".", fRawDataList)).IsNull() && !fUseMeta)
    gSystem -> Exec("./createList.sh "+sRunNo);

  TString fMetaDataList;
  if (fUseMeta) {
    fRawDataList = Form("%s/run_%04d/dataList.txt", fSupplePath.Data(), fRunNo);
    fMetaDataList = Form("%s/run_%04d/metadataList.txt", fSupplePath.Data(), fRunNo);
  }

  FairLogger *logger = FairLogger::GetLogger();
  logger -> SetLogToScreen(true);

  FairParAsciiFileIo* parReader = new FairParAsciiFileIo();
  parReader -> open(par);

  FairRunAna* run = new FairRunAna();
  run -> SetGeomFile(geo);
  run -> SetOutputFile(fOut);
  run -> GetRuntimeDb() -> setSecondInput(parReader);

  STDecoderTask *decoder = new STDecoderTask();
  decoder -> SetUseSeparatedData(true);
  decoder -> SetPersistence(false);
  // By default, if SetUseGainCalibration(true) is called, reading gain calibration information from parameter file.
  decoder -> SetUseGainCalibration(true);
  /* Manual calibration parameter setters. You need to provide both calibration root file and reference values to match.
   * Use default automatic one if you're not sure what you're doing. */
  // decoder -> SetGainCalibrationData("filename.root");
  // decoder -> SetGainReference(constant, linear, quadratic);
  if (fIsGGDataSet)
    decoder -> SetGGNoiseData(fGGData);
  decoder -> SetDataList(fRawDataList);
  //decoder -> SetEventID(start);
  decoder -> SetTbRange(30, 257); 
  // Low gain calibration. Don't forget you need to uncomment PSA part, too.
  if (fUseGainMatching)
    decoder -> SetGainMatchingData(fGainMatchingFile);
  // Method to select events to reconstruct
  // Format of the input file:
  //        runid eventid
  //        runid eventid
  //        runid eventid
  //        runid eventid

  if (fSelectPiEvent)
  {
    map<Int_t, vector<Int_t> *> events;
    string FileName_PiEvt = Form("./Pick_PiEvt/Sn%d_PiEvt/",fSys);
    FileName_PiEvt = FileName_PiEvt+Form("Sn%d_Run",fSys)+fRunNo+"_PiEvt";
    cout<<"Reading the Event list for the pion events : "<<FileName_PiEvt<<endl;
    readEventList(FileName_PiEvt, events);
    cout <<"Number of events " << fNumEventsInSplit << " starting at " << start <<endl;

    decoder -> SetEventList(*events[fRunNo]);
  }

  decoder -> SetEventID(start);
 
  
  if (fUseMeta) 
  {
    std::ifstream metalistFile(fMetaDataList.Data());
    TString dataFileWithPath;
    for (Int_t iCobo = 0; iCobo < 12; iCobo++) {
      dataFileWithPath.ReadLine(metalistFile);
      dataFileWithPath = Form("%s/run_%04d/%s", fSupplePath.Data(), fRunNo, dataFileWithPath.Data());
      decoder -> SetMetaData(dataFileWithPath, iCobo);
    }
  }

  auto embedTask = new STEmbedTask();
  embedTask -> SetEventID(start);
  embedTask -> SetEmbedFile(fMCFile);

  auto preview = new STEventPreviewTask();
  preview -> SetSkippingEvents(fSkipEventArray);
  preview -> SetPersistence(true);
  //preview -> SetSelectingEvents(*events[fRunNo]);

  auto psa = new STPSAETask();
  psa -> SetPersistence(false);
  // In case pulse cut threshold need to be changed. (Default: 30)
  // psa -> SetThreshold(value);
  // By default, it uses all the layers.
  // psa -> SetLayerCut(-1, 112);
  // Pulse having long tail
  psa -> SetPulserData("pulser_117ns_50tb.dat");
  // Default: Rensheng's peak finding method 1. Finds small peaks
  // psa -> SetPSAPeakFindingOption(0);
  // Low gain calibration. Don't forget you need to uncomment decoder part, too.
  if (fUseGainMatching)
    psa -> SetGainMatchingData(fGainMatchingFile);
  // Statistically matching time jitter in different CoBos
  psa -> SetYOffsets(fSpiRITROOTPath + "parameters/yOffsetCalibration.dat");
  // This is used to match the TPC-Vertex_Y with the BDC_Y.
  psa -> SetYPedestalOffset(fYPedestalOffset); // unit: mm
  
  auto helix = new STHelixTrackingTask();
  helix -> SetPersistence(false);
  helix -> SetClusterPersistence(false);
  // Left, right, top and bottom sides cut
  helix -> SetClusterCutLRTB(420, -420, -64 + fYPedestalOffset, -522 + fYPedestalOffset);
  // High density region cut
  helix -> SetEllipsoidCut(TVector3(0, -260 + fYPedestalOffset, -11.9084), TVector3(120, 100, 220), 5); // current use
  // Changing clustering direction angle and margin. Default: 45 deg with 0 deg margin
  // helix -> SetClusteringAngleAndMargin(35., 3.);

  auto correct = new STCorrectionTask(); //Correct for saturation
  
  // Once created, all the subroutines share one field map. Offset changed in individual instance won't be in effect.
  auto gfBField = STGFBField::GetInstance("samurai_field_map", "A", fFieldOffsetX, fFieldOffsetY, fFieldOffsetZ);   

  auto spaceCharge = new STSpaceChargeCorrectionTask();
  spaceCharge -> SetBField(gfBField -> GetFieldMap());
  if (fSheetChargeDensity != 0) {
    spaceCharge -> SetDriftParameters(-4.355e4, -2.18); // omega tau and mu of the Langevin equation
    const double fLeakageChargeFactor = 10.2;
    spaceCharge -> SetSheetChargeDensity(fSheetChargeDensity, fLeakageChargeFactor*fSheetChargeDensity);
    switch (fSystemID) {
      case 132124: spaceCharge -> SetProjectile(STSpaceCharge::Projectile::Sn132); break;
      case 124112: spaceCharge -> SetProjectile(STSpaceCharge::Projectile::Sn124); break;
      case 112124: spaceCharge -> SetProjectile(STSpaceCharge::Projectile::Sn112); break;
      case 108112: spaceCharge -> SetProjectile(STSpaceCharge::Projectile::Sn108); break;
    }
    spaceCharge -> SetElectronDrift(true);
  } else
    spaceCharge -> SetElectronDrift(false);


  auto genfitPID = new STGenfitPIDTask();
  // In the TPC frame. Here the z position is used when Genfit do the extrapolation.
  genfitPID -> SetTargetPlane(0, 0, fTargetZ); // unit: mm
  genfitPID -> SetPersistence(true);
  genfitPID -> SetBDCFile("");
  // Only for test
  // genfitPID -> SetConstantField();
  genfitPID -> SetListPersistence(true);
  // Removing shorter length tracklet by distance of adjacent clusters.
  // genfitPID -> SetMaxDCluster(60);

  auto genfitVA = new STGenfitVATask();
  genfitVA -> SetPersistence(true);
  // Only for test
  // genfitVA -> SetConstantField();
  genfitVA -> SetListPersistence(true);
  if (fIsBeamDataSet) {
    genfitVA -> SetBeamFile(fBeamData);
    genfitVA -> SetInformationForBDC(fRunNo, fBDCOffsetX, fBDCOffsetY, 0);
  }
  // Uncomment if you want to recalculate the vertex using refit tracks.
  genfitVA -> SetUseRave(true);
  ////////////////////////////////////////////////////////////////////// XXX ////////////////////////////////////////////////////////////////////// 
  genfitVA -> ShiftBDCAfterSC("/home/ejungwoo/spiritroot/parameters/vertex_shift_Sn" + std::to_string(int(fSystemID/1000)) + ".root", 80);
  ////////////////////////////////////////////////////////////////////// XXX ////////////////////////////////////////////////////////////////////// 
  genfitVA -> SetZtoProject(-13.2, 1.7, 3); //(Double_t peakZ, Double_t sigma, Double_t sigmaMultiple), this function will project the BDC on the Target.

  auto embedCorr = new STEmbedCorrelatorTask();
  embedCorr -> SetPersistence(true);

  auto smallOutput = new STSmallOutputTask();
  smallOutput -> SetOutputFile(fConc.Data());
  smallOutput -> SetRun(fRunNo);

  run -> AddTask(decoder);
  if(!fMCFile.IsNull())
    run -> AddTask(embedTask);
  run -> AddTask(preview);
  run -> AddTask(psa);
  run -> AddTask(helix);
  run -> AddTask(correct);
  run -> AddTask(spaceCharge);
  run -> AddTask(genfitPID);
  run -> AddTask(genfitVA);
  if(!fMCFile.IsNull())
    run -> AddTask(embedCorr);
  run -> AddTask(smallOutput);

  run -> Init();
  run -> Run(0,fNumEventsInSplit);
  
  auto outFile = FairRootManager::Instance() -> GetOutFile();
  auto recoHeader = new STRecoHeader("RecoHeader","");
  recoHeader -> SetPar("version", version);
  recoHeader -> SetPar("eventStart", start);
  recoHeader -> SetPar("numEvents", fNumEventsInSplit);
  recoHeader -> SetPar("parameter", fParameterFile);
  if (fIsGGDataSet)
    recoHeader -> SetPar("GGData", fGGData);
  recoHeader -> SetPar("yPedestal", fYPedestalOffset);
  recoHeader -> SetPar("fieldOffsetX", fFieldOffsetX);
  recoHeader -> SetPar("fieldOffsetY", fFieldOffsetY);
  recoHeader -> SetPar("fieldOffsetZ", fFieldOffsetZ);
  recoHeader -> SetPar("sheetChargeDensity", fSheetChargeDensity);
  recoHeader -> SetPar("BDCOffsetX", fBDCOffsetX);
  recoHeader -> SetPar("BDCOffsetY", fBDCOffsetY);
  auto driftVelocityInParameterFile = gSystem -> GetFromPipe("cat " + par + " | grep DriftVelocity | awk '{print $2}'");
  recoHeader -> SetPar("driftVelocity", driftVelocityInParameterFile);
  recoHeader -> Write("RecoHeader");

  cout << "Log    : " << fLog << endl;
  cout << "Input  : " << fRawDataList << endl;
  cout << "Output : " << fOut << endl;

  gApplication -> Terminate();
}

void readEventList(TString eventListFile, map<Int_t, vector<Int_t> *> &events) {
  vector<Int_t> *temp = new vector<Int_t>;

  ifstream eventList(eventListFile.Data());
  Int_t numEvents = 0;
  Int_t oldRunid = 0;
  Int_t runid, eventid;

  while (1) {
    eventList >> runid >> eventid;

    if (eventList.eof()) {
      events.insert(make_pair(oldRunid, temp));

      break;
    }
    
    if (oldRunid == 0)
      oldRunid = runid;

    if (oldRunid != runid) {
      events.insert(make_pair(oldRunid, temp));

      oldRunid = runid;
      temp = new vector<Int_t>;
    }

    temp -> push_back(eventid);
  }
}
