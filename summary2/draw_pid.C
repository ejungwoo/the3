#include "KBGlobal.hh"
#include "init_variables.h"
#include "binning.h"

TString fVName;
TString fVName1;
vector<TCanvas *> fCvsArray;
int cvsXOff = 1300;
//int cvsXOff = 0;
int gIdxSys = 0;
int gIdxParticle = 0;
int countCvs = 0;
bool fUseHandCut = false;
double fTMargin = 0.12;
double fBMargin = 0.20;
//double fLMargin = 0.19;
double fLMargin = 0.25;
double fRMargin = 0.055;
double fRMargin1 = 0.055;

int fSelAna = 0;
int fSelSDType = 0;
double fSDValue = 3;

bool setDrawNZAll = kSet;

TGraph *graphFitPIDMean[fNumSyss][fNumParticles] = {0};
TGraph *graphFitPIDAmp[fNumSyss][fNumParticles] = {0};
TH1F *histPIDMeta[fNumSyss][fNumParticles] = {0};
TF1 *f1PIDMean[fNumSyss][fNumParticles] = {0};
TF1 *f1PIDSigma[fNumSyss][fNumParticles] = {0};
TCutG *cutgPID[fNumSyss][fNumParticles] = {0};
TGraph *graphPIDMean[fNumSyss][fNumParticles] = {0};
TGraph *graphPIDRange[fNumSyss][fNumParticles][2] = {0};

void saveAll();
void savePNG();
void savePDF();
void writeAll();
double calculate_prob(double p_lab, double dedx);
void project(TTree *tree, const char *name, const char *expr, TCut selection, bool verbose=1, TH1 *hist=nullptr);
TCanvas *makeCvs2(const char *name, int w=1050, int h=950);
TCanvas *makeCvs(const char *name, int w=630, int h=550);
void cvsDivide(TCanvas *cvs, int nx, int ny);
void cvsDivideM0(TCanvas *cvs, int nx, int ny);
TLegend *makeLegend(TVirtualPad *cvs, TLegend *legend, TString opt = "", double x_offset=0, double y_offset=0, double width_fixed=0, double height_fixed=0);
void setAtt(TH1 *hist, int iDraw);
void setAtt(TGraph *graph, int iDraw, int option=1);
void setAtt(TF1 *fit, int iDraw);
void setAtt(TMarker *marker, int iDraw, int iColor=-1, int option=1);
const char *makeName(const char *mainName, int iAna, int iLR, int iMult, int iSys, int iCutTheta, int iCutYP, int iPart=-1);
const char *makeTitle(const char *mainName, int iAna, int iLR, int iMult, int iSys, int iCutTheta, int iCutYP, int iPart=5);



void draw_pid
(
    //int selAna = kF132, double selSDValue = 3.0,
    //int selAna = kFSys, double selSDValue = 3.0,
    int selAna = kFNN, double selSDValue = 3.0,

    int selSDType = kSD_xx,
    //int selSDType = kSD_xx_l3,
    //int selSDType = kSD_0_x,

    //int selCutYP0 = kYPFI,
    int selCutYP0 = -2,

    bool checkContent = false,
    bool saveCvsPNG = false,
    bool saveCvsRoot = false
)
{
  if (selAna>=0) fSelAna = selAna;
  if (selSDType>=0) fSelSDType = selSDType;
  if (selSDValue>0) fSDValue = selSDValue;

  gStyle -> SetPaintTextFormat("0.2g");

  
  //TCut cut0 = "eff>0.05&&prob>.7"; fVName = Form("v%s_sd%.1f_p0p7",fAnaNames[selAna],selSDValue);
  //TCut cut0 = "eff>0.05&&prob>.7&&theta_cm*TMath::RadToDeg()<90"; fVName = Form("v%s_sd%.1f_p0p7",fAnaNames[selAna],selSDValue);
  //TCut cut0 = "eff>0.05&&prob>.7&&theta_cm*TMath::RadToDeg()<30"; fVName = Form("v%s_sd%.1f_p0p7",fAnaNames[selAna],selSDValue);

  //TCut cut0 = "eff>0.1&&prob>.1&&theta_cm*TMath::RadToDeg()<90";
  TCut cut0 = "eff>0.05&&prob>.7";
  TCut cut1 = "eff>0.05&&prob>.5";
  //TCut cut0 = "eff>0.05&&prob>.7";
  //TCut cut1 = "eff>0.05&&prob>.1";

  TString probString = "prob";

  ///////////////////////////////////////////////////////////////////////////////////////////////////

  bool drawRawPID = kHold;
  bool setDrawRawPIDProj = kUnSet;
  bool setDrawRawPIDSub = kSet;
  bool setFitdEdx = kSet;

  bool drawCorPID = kHold;
  bool setTestCorPIDJustPDT = kUnSet;
  bool setDrawCorPIDInRaw = kUnSet;
  bool setDrawCorPIDInPtoa = kUnSet;

  bool setDrawPIDLogz = kSet;
  bool setDrawGuideLine = kUnSet;
  bool setDrawGuideLine132 = kUnSet;
  bool setDrawHandCut = fUseHandCut;
  bool setZoomPID0 = kSet;
  bool setZoomPIDZ1 = kUnSet;
  bool setZoomPIDZ2 = kUnSet;

  bool drawCorEKE = kHold;
  bool drawPtoa = kHold;

  bool drawKY = kHold;

  bool drawKT = kDraw;
  bool setDrawKTGrid = kSet;
  bool setFitKE = kUnSet;

  bool drawYP = kDraw;
  bool setDrawYPTogether = kSet;
  bool setDrawYPTheta0 = kSet;
  bool setDrawYPKE = kSet;
  bool setDrawYPPoz = kUnSet;
  bool setDrawYPText = kSet;
  bool setDrawYPGrid = kSet;
  bool setUseCM = kSet;

  bool drawPtoaR21 = kDraw;
  bool drawY0R21 = kDraw;
  bool setDrawR21 = kSet;
  bool setDrawSR = kUnSet;
  bool setDrawDR = kUnSet;
  bool setDrawTemp = kSet;
  bool setDrawR21AMD = kUnSet;
  bool setDrawPN = kUnSet;

  bool drawTLabR21 = kHold;

  bool drawTCMR21 = kDraw;
  bool drawKeoaR21 = kDraw;
  bool setDrawYield = kUnSet;

  bool drawNZR21 = kDraw;
  bool setDrawABNZ = kUnSet;
  bool setDrawABTK = kSet;
  bool setDrawABC = kUnSet;
  bool setDrawABError = kSet;
  //bool setDrawNZInKT = kUnSet;
  bool setDrawNZInKT = kSet;
  bool setFitNZR21TG = kSet;
  bool setDrawRefAB = kSet;
  bool setFitABSep = kUnSet;
  //bool setFitABSep = kSet;
  bool setDrawFit3 = kSet;
  bool setWriteR21ABCvsFile = kSet;
  bool setABErrorToResidual = kSet;
  bool setDrawNZLogy = kSet;
  bool setDrawNZFitLegend = kUnSet;
  bool setIgnoreT = kUnSet;
  /**/ setDrawNZAll = kUnSet;
  bool setDrawSN = kUnSet;

  bool drawDistKeoa = kHold;
  bool drawDistTCM = kHold;

  auto holdAll = [
    &drawRawPID, &drawCorPID,
    &drawCorEKE, &drawKT, &drawPtoa, &drawYP,
    &drawTLabR21, &drawPtoaR21, &drawY0R21, &drawKeoaR21, &drawTCMR21, &drawNZR21,
    &drawDistKeoa
  ]()
  {
    drawRawPID   = kHold;
    drawCorPID   = kHold;
    drawCorEKE   = kHold;
    drawKT       = kHold;
    drawPtoa     = kHold;
    drawYP       = kHold;
    drawTLabR21  = kHold;
    drawPtoaR21  = kHold;
    drawY0R21    = kHold;
    drawTCMR21   = kHold;
    drawKeoaR21  = kHold;
    drawNZR21    = kHold;
    drawDistKeoa = kHold;
  };

  if (setTestCorPIDJustPDT)
    setDrawPIDLogz = kUnSet;

  // sssssssssssssssss //////////////////////////////////////////////////////////////////////////////

  //int selMult = kMultAll;
  int selMult = kMult55;
  //int selMult = kMult45;
  if (selAna==kFNN50)
    selMult = kMult50;

  int selSys = kAll;
  //int selSys = k132;
  //vector<int> selSysIdx = {k132,k108,k112,k124};
  vector<int> selSysIdx = {k132,k108};
  //vector<int> selSysIdx = {k112,k124};

  int selSysComb = kAll;
  //int selSysComb = 0;
  vector<int> selSysCombIdx = {0};
  //vector<int> selSysCombIdx = {0,1,2,3};
  //vector<int> selSysCombIdx = {3};

  int selLR = kLR;
  vector<int> selLRIdx = {kLR,kLeft,kRight};

  int selCutTheta = kThetaAll;
  vector<int> selCutThetaIdx = {kThetaAll, kTheta0, kTheta20, kTheta40, kTheta60, kThetaLT60, kThetaGT60};

  int selCutYP = kYPAll;
  //int selCutYP = kYPUD;
  vector<int> selCutYPIdx = {kYPAll, kY1, kY2, kY3, kY4, kPtoa0, kPtoa50, kPtoa100, kPtoa150, kPtoa200, kPtoa250, kPtoa300, kPtoa350, kYF, kPF3, kYPFI, kPF4, kYPUD, kP0, kP100, kP200, kP300};

  //int selCutYP = kAll;
  //vector<int> selCutYPIdx = {kY2, kY3, kY4};

  if (selCutYP0 >= 0)
    selCutYP = selCutYP0;

  if (drawCorPID && setDrawCorPIDInPtoa)
  {
    setDrawCorPIDInRaw = kSet;

    selSys = kAll;
    selSysIdx.clear();
    for (auto i : {0,1})
      selSysIdx.push_back(i);

    selCutYP = kAll;
    selCutYPIdx.clear();
    for (auto i : {kPtoa0, kPtoa50, kPtoa100, kPtoa150, kPtoa200, kPtoa250, kPtoa300, kPtoa350})
      selCutYPIdx.push_back(i);
  }

  ///////////////////////////////////////////////////////////////////////////////////////////////////

  //double cab132K8[8][3] = {{0.965368,0.28,-0.17},{1.01768,0.31,-0.22},{1.04899,0.30,-0.26},{1.07809,0.25,-0.28},{1.0623,0.17,-0.27},{1.08558,0.09,-0.31},{1.04587,0.02,-0.28},{0.965415,-0.07,-0.22}};

  ///////////////////////////////////////////////////////////////////////////////////////////////////


  int nbinsFrame = 100;


  binning nixPIDProjX(200,800,820); // number of total projections, min x, max y
  binning bnPoz(800,0,3000,"p/Z (MeV/c)");
  binning bndEdx(800,0,2000,"dE/dx");
  //binning bnPIDz(0,0.08,0.0001);
  binning bnPIDz(-1,0.08,0.0001);
  if (setDrawCorPIDInRaw) {
    //bnPIDz.set(1,1,800);
    bnPIDz.set(1,1,50000);
  }

  if (setZoomPIDZ2) {
    //bnPoz.set(100,200,2200);
    //bndEdx.set(200,0,1200);
    bnPoz.set(200,100,800);
    bndEdx.set(200,0,2000);
  }
  else if (setZoomPIDZ1) {
    bnPoz.set(100,0,3000);
    bndEdx.set(100,0,200);
  }
  else if (setZoomPID0) {
    bnPoz.set(200,0,2000);
    bndEdx.set(200,0,800);
  }

  
  
  binning bnY0(80,-1,2,"y_{0}");
  binning bnPoACM(80,0,1000,"p_{c.m.}/A (MeV/c)");
  binning bnPtoa(80,0,800,"p_{T}/A (MeV/c)");
  binning bnKeoa(100,0,400,"KE_{Lab}/A (MeV)");
  binning bnKeoa2(100,0,1500,"KE/A (MeV)");

  binning bnKeoaCM(42,0,210,"KE_{CM}/A (MeV)");
  binning bnRKeoaCM(4,0,120,"KE_{CM}/A (MeV)");
  //binning bnRKeoaCM(5,0,100,"KE_{CM}/A (MeV)");
  //binning bnRKeoaCM(4,0,40,"KE_{CM}/A (MeV)");

  binning bnTCM(60,0,180,"#theta_{CM} (deg)");
  binning bnRTCM(4,0,90,"#theta_{CM} (deg)");

  int dxCvsAB = 150; int dyCvsAB = 200;
  if (bnRKeoaCM.fN==1 && setDrawNZInKT) {
    dxCvsAB = 300; dyCvsAB = 150;
  }

  if (setDrawNZInKT) {}
  else {
    setFitABSep = false;
  }

  binning bnParticleZ(100,0.2,2.8,"Z");
  binning bnParticleN(100,-.8,2.8,"N");
  //binning bnParticleN(100,-1,3,"N");
  //binning bnR21(100,0,2.0,"R_{21}");
  //binning bnR21(100,0.3,1.8,"R_{21}");
  //binning bnR21(100,0.3,2.1,"R_{21}");
  binning bnR21(100,0.3,2.3,"R_{21}");
  //binning bnRY0(5,-0.25,1,"y_{0}"); // selected;
  //binning bnRY0(10,-1,1); // selected;
  binning bnRTLab(4,0,80,"#theta_{Lab.}");

  //binning bnRY0(3,-0.25,1.25,"y_{0}"); // selected;
  //binning bnRPtoa(8,0,400,"p_{T}/A (MeV/c)");
  //binning bnRPtoa(4,0,400,"p_{T}/A (MeV/c)");
  //binning bnRPtoa(10,0,500,"p_{T}/A (MeV/c)");
  //binning bnRY0(3,-.25,1.25,"y_{0}"); // selected;

  //binning bnRY0(4,-.25,1.25,"y_{0}"); // selected;
  //binning bnRY0(4,0,1.00,"y_{0}"); // selected;
  binning bnRY0(4,-.15,1.05,"y_{0}"); // selected;
  binning bnRPtoa(4,0,400,"p_{T}/A (MeV/c)");

  //binning bnRY0(4
  //binning bnRPtoa(4,0,90,"#theta_{CM}");;
  
  if (drawKeoaR21) {
    bnRKeoaCM.set(4,0,120);
    bnRTCM.set(1,0,90);
    //bnR21.set(100,0.3,1.6,"R_{21}");
  }

  if (drawDistKeoa) {
    //selCutThetaIdx.clear();
    //selCutTheta = kAll;
    //for (auto i : {1,2,3,4})
      //selCutThetaIdx.push_back(i);

    holdAll();
    drawDistKeoa = true;
    selSysIdx.clear();
    selSysIdx.push_back(k132);
    //selSysIdx.push_back(k112);
    selSysIdx.push_back(k108);
    bnKeoa.set(100,0,200);
    bnKeoaCM.set(100,0,200);
    bnRKeoaCM.set(100,0,200);
  }

  binning bnRKeoaLab(8,0,400);

  if (drawDistKeoa)
    bnRKeoaLab.fN = 50;

  if (setDrawTemp) {
    //bnRY0.set(6,-.5,1.0);
    //bnRPtoa.set(10,0,400);
  }

  if (setDrawPN) {
    bnRY0.set(10,-1,2);
    bnRPtoa.set(10,0,500);
  }

  if (setDrawSR || setDrawDR) bnRY0.set(21,-.8,1.3);

  binning bnPR(100,0,1.1,"SR");
  binning bnDR(100,1,2.2,"DR");
  binning bnTemp(100,5,15,"T_{app}");
  binning bnPN(100,0,20);

  //binning bnAlpha(100,-.2,.5,"#alpha");
  //binning bnBeta(100,-.5,.2,"#beta");

  binning bnAlpha(100,-.1,.5,"#alpha");
  binning bnBeta(100,-.5,.1,"#beta");

  //binning bnAlpha(100,-.2,.2,"#alpha");
  //binning bnBeta(100,-.2,.2,"#beta");

  double scaleMaxFitdEdx = .05;

  if (drawTLabR21) {
    bool wasDrawYP = drawYP;
    holdAll();
    drawTLabR21 = true;
    drawYP = wasDrawYP;
    setDrawYPTheta0 = kSet;
    setDrawYPText = kSet;
    selCutTheta = kAll;
    selCutThetaIdx.clear(); 
    selCutThetaIdx.push_back(kThetaAll); 
  }

  if (drawNZR21 && setDrawNZAll) {
    holdAll();
    drawNZR21 = true;
    if (setDrawNZInKT) {
      bnRKeoaCM.setN(1);
      bnRTCM.setN(1);
    }
    else {
      bnRY0.setN(1);
      bnRPtoa.setN(1);
    }
  }


  ///////////////////////////////////////////////////////////////////////////////////////////////////

  //fVName = Form("v%s_sd%.1f_p0p1",fAnaNames[selAna],selSDValue);
  fVName = Form("v%s_Y0N%d_PtoaN%d",fAnaNames[selAna],bnRY0.fN,bnRPtoa.fN);
  fVName1 = Form("v%s_Y0N%d_PtoaN1",fAnaNames[selAna],bnRY0.fN);
  if (setDrawNZInKT) {
    fVName = Form("v%s_KeoaN%d_ThetaN%d",fAnaNames[selAna],bnRKeoaCM.fN,bnRTCM.fN);
    fVName1 = Form("v%s_KeoaN%d_ThetaN1",fAnaNames[selAna],bnRKeoaCM.fN);
  }

  cout << endl;
  cout_info << "== " << fVName << " ; " << " " << cut0 << " ; " << probString << endl;

  if (drawRawPID  ) cout_info << "set draw RawPID  " << endl;
  if (drawCorPID  ) cout_info << "set draw CorPID  " << endl;
  if (drawCorEKE  ) cout_info << "set draw CorEKE  " << endl;
  if (drawKT      ) cout_info << "set draw KT      " << endl;
  if (drawPtoa    ) cout_info << "set draw Ptoa    " << endl;
  if (drawYP      ) cout_info << "set draw YP      " << endl;
  if (drawTLabR21 ) cout_info << "set draw ThetaR21" << endl;
  if (drawPtoaR21 ) cout_info << "set draw PtoaR21 " << endl;
  if (drawY0R21   ) cout_info << "set draw Y0R21   " << endl;
  if (drawTCMR21  ) cout_info << "set draw ThetaR21" << endl;
  if (drawKeoaR21 ) cout_info << "set draw KeoaR21 " << endl;
  if (drawNZR21   ) cout_info << "set draw NZR21   " << endl;
  if (drawDistKeoa) cout_info << "set draw DistKeoa" << endl;

  cout_info << "Ana      is = " << setw(10) << (selAna     ==kAll?"everything":     fAnaNames[selAna     ]) << " ;" << setw(3) << selAna      << endl;
  cout_info << "Mult     is = " << setw(10) << (selMult    ==kAll?"everything":    fMultNames[selMult    ]) << " ;" << setw(3) << selMult     << endl;
  cout_info << "LR       is = " << setw(10) << (selLR      ==kAll?"everything":      fLRNames[selLR      ]) << " ;" << setw(3) << selLR       << " ("; for (auto v :       selLRIdx) cout << v << ","; cout << ")" << endl;
  cout_info << "Sys      is = " << setw(10) << (selSys     ==kAll?"everything":     fSysNames[selSys     ]) << " ;" << setw(3) << selSys      << " ("; for (auto v :      selSysIdx) cout << v << ","; cout << ")" << endl;
  cout_info << "CutTheta is = " << setw(10) << (selCutTheta==kAll?"everything":fCutThetaNames[selCutTheta]) << " ;" << setw(3) << selCutTheta << " ("; for (auto v : selCutThetaIdx) cout << v << ","; cout << ")" << endl;
  cout_info << "CutYP    is = " << setw(10) << (selCutYP   ==kAll?"everything":   fCutYPNames[selCutYP   ]) << " ;" << setw(3) << selCutYP    << " ("; for (auto v :    selCutYPIdx) cout << v << ","; cout << ")" << endl;
  cout_info << "SysComb  is = " << setw(10) << (selSysComb ==kAll?"everything": fSysCombNames[selSysComb ]) << " ;" << setw(3) << selSysComb  << " ("; for (auto v :  selSysCombIdx) cout << v << ","; cout << ")" << endl;

  if (checkContent) { 
    TString stopIf0;
    cout << "enter 0 or q to stop: ";
    cin >> stopIf0; cout << endl;
    if (stopIf0=="0"||stopIf0=="q")
      return;
  }

  ///////////////////////////////////////////////////////////////////////////////////////////////////

  gSystem -> mkdir(fVName);
  gSystem -> mkdir(fVName+"/rooto");
  gSystem -> mkdir(fVName+"/figures");
  gSystem -> mkdir(fVName+"/figures_pdf");
  gSystem -> mkdir(fVName+"/others");

  TH2D *histAPrev[4] = {0};
  TH2D *histBPrev[4] = {0};
  TH2D *histCPrev[4] = {0};
  TH2D *histTPrev[4] = {0};

  TH2D *histAPrev1[4] = {0};
  TH2D *histBPrev1[4] = {0};
  TH2D *histCPrev1[4] = {0};
  TH2D *histTPrev1[4] = {0};

  for (auto iComb : selSysCombIdx)
  {
    const char *namePrev = Form("%s/rooto/alpha_beta_%d.root",fVName.Data(),iComb);
    auto filePrev = new TFile(namePrev,"read");
    if (!filePrev->IsZombie()) 
    {
      histAPrev[iComb] = (TH2D *) filePrev -> Get("alpha");
      histBPrev[iComb] = (TH2D *) filePrev -> Get("beta");
      histCPrev[iComb] = (TH2D *) filePrev -> Get("cnorm");
      histTPrev[iComb] = (TH2D *) filePrev -> Get("temperature");

      histAPrev[iComb] -> SetName(Form("icb%d_alphaPrev",iComb));
      histBPrev[iComb] -> SetName(Form("icb%d_betaPrev",iComb));
      histCPrev[iComb] -> SetName(Form("icb%d_cnormPrev",iComb));
      histTPrev[iComb] -> SetName(Form("icb%d_temperaturePrev",iComb));
    }

    const char *namePrev1 = Form("%s/rooto/alpha_beta_%d.root",fVName1.Data(),iComb);
    auto filePrev1 = new TFile(namePrev1,"read");

    if (filePrev1->IsZombie()) 
      cout << namePrev1 << " ?" << endl;
    else
    {
      histAPrev1[iComb] = (TH2D *) filePrev1 -> Get("alpha");
      histBPrev1[iComb] = (TH2D *) filePrev1 -> Get("beta");
      histCPrev1[iComb] = (TH2D *) filePrev1 -> Get("cnorm");
      histTPrev1[iComb] = (TH2D *) filePrev1 -> Get("temperature");

      histAPrev1[iComb] -> SetName(Form("icb%d_alphaPrev1",iComb));
      histBPrev1[iComb] -> SetName(Form("icb%d_betaPrev1",iComb));
      histCPrev1[iComb] -> SetName(Form("icb%d_cnormPrev1",iComb));
      histTPrev1[iComb] -> SetName(Form("icb%d_temperaturePrev1",iComb));

      histAPrev1[iComb] -> SetTitle("alpha");
      histBPrev1[iComb] -> SetTitle("beta");
      histCPrev1[iComb] -> SetTitle("cnorm");
      histTPrev1[iComb] -> SetTitle("common temperature for each KE range");

      makeCvs2(Form("tprev1_icb%d",iComb),500,500);
      histTPrev1[iComb] -> SetMarkerSize(2);
      histTPrev1[iComb] -> Draw("coltext");
    }
  }


  int numEventsInAna[fNumSyss] = {0};

  auto getWeighting = [&numEventsInAna,setDrawRawPIDProj,probString](int iSys, int iParticle, double probValue=-1) {
    const char *finalProbString = probString;
    if (probValue>0)
      finalProbString = Form("%f",probValue);
    else if (setDrawRawPIDProj) {
      gIdxSys = iSys;
      gIdxParticle = iParticle;
      finalProbString = "calculate_prob(p_lab,dedx)";
    }
    TCut weighting = (Form("%s/eff/%d",finalProbString,numEventsInAna[iSys]));
    //TCut weighting = (Form("%s/1/%d",finalProbString,numEventsInAna[iSys]));
    //TCut weighting = (Form("%s/1/%d","1",numEventsInAna[iSys]));
    return weighting;
  };

  ///////////////////////////////////////////////////////////////////////////////////////////////////

  const int nAMDY0 = 20;
  binning bnAMDY0(nAMDY0,-2,2);
  double amdYield[3][2][3][nAMDY0]; // type(default,2sigma_NN, 2sigma_NN+phase-space), system(132,108) particle
  if (setDrawR21AMD)
  {
    std::ifstream amdFile("data2/amd/amd_dndy_table.csv");
    std::string aline;
    double yamd, p132, d132, t132, p108, d108, t108;
    std::getline(amdFile, aline);
    for (auto amdType : {0,1,2}) {
      std::getline(amdFile, aline);
      std::getline(amdFile, aline);
      //cout << amdType << endl;
      for (int iLineY0=0; iLineY0<bnAMDY0.fN; ++iLineY0) {
        std::getline(amdFile, aline);
        stringstream((TString(aline).ReplaceAll(","," ")).Data()) >> yamd >> p132 >> d132 >> t132 >> p108 >> d108 >> t108;
        //cout << iLineY0 << " " << yamd << ": " << p132 << " " << d132 << " " << t132 << ", " << p108 << " " << d108 << " " << t108 << endl;
        amdYield[amdType][0][0][iLineY0] = p132;
        amdYield[amdType][0][1][iLineY0] = d132;
        amdYield[amdType][0][2][iLineY0] = t132;
        amdYield[amdType][1][0][iLineY0] = p108;
        amdYield[amdType][1][1][iLineY0] = d108;
        amdYield[amdType][1][2][iLineY0] = t108;
      }
    }
  }

  ///////////////////////////////////////////////////////////////////////////////////////////////////

  auto f1Sinx = new TF1("satta","sin(x)",0,TMath::Pi());

  TCutG *cutgAll[fNumSyss][fNumLRs][fNumCutThetas][fNumParticles];

  auto GetCutG = [&cutgAll](int iSys, int iLR, int iCutTheta, const char *fileName) {
    auto file = new TFile(fileName);
    for (auto iParticle : {0,1,2,3,4}) {
      auto bound = (TCutG *) file -> Get(Form("bound%d",iParticle));
      bound -> SetName(Form("bound_%d_%d_%d_%d",iSys,iLR,iCutTheta,iParticle));
      bound -> SetVarX("p_lab");
      bound -> SetVarY("dedx");
      cutgAll[iSys][iLR][iCutTheta][iParticle] = bound;
    }
  };

  GetCutG(1,1,0,"vAll/rooto/bound_pidRaw_f7_left_45_100_108_ttaRG0_20_yAll.root");
  GetCutG(1,1,1,"vAll/rooto/bound_pidRaw_f7_left_45_100_108_ttaRG0_20_yAll.root");
  GetCutG(1,1,2,"vAll/rooto/bound_pidRaw_f7_left_45_100_108_ttaRG20_40_yAll.root");
  GetCutG(1,1,3,"vAll/rooto/bound_pidRaw_f7_left_45_100_108_ttaRG40_60_yAll.root");
  GetCutG(1,1,4,"vAll/rooto/bound_pidRaw_f7_left_45_100_108_ttaRG40_60_yAll.root");

  GetCutG(1,2,0,"vAll/rooto/bound_pidRaw_f7_right_45_100_108_ttaRG0_20_yAll.root");
  GetCutG(1,2,1,"vAll/rooto/bound_pidRaw_f7_right_45_100_108_ttaRG0_20_yAll.root");
  GetCutG(1,2,2,"vAll/rooto/bound_pidRaw_f7_right_45_100_108_ttaRG20_40_yAll.root");
  GetCutG(1,2,3,"vAll/rooto/bound_pidRaw_f7_right_45_100_108_ttaRG40_60_yAll.root");
  GetCutG(1,2,4,"vAll/rooto/bound_pidRaw_f7_right_45_100_108_ttaRG40_60_yAll.root");

  GetCutG(0,1,0,"vAll/rooto/bound_pidRaw_f7_left_45_100_132_ttaRG0_20_yAll.root");
  GetCutG(0,1,1,"vAll/rooto/bound_pidRaw_f7_left_45_100_132_ttaRG0_20_yAll.root");
  GetCutG(0,1,2,"vAll/rooto/bound_pidRaw_f7_left_45_100_132_ttaRG20_40_yAll.root");
  GetCutG(0,1,3,"vAll/rooto/bound_pidRaw_f7_left_45_100_132_ttaRG40_60_yAll.root");
  GetCutG(0,1,4,"vAll/rooto/bound_pidRaw_f7_left_45_100_132_ttaRG40_60_yAll.root");

  GetCutG(0,2,0,"vAll/rooto/bound_pidRaw_f7_right_45_100_132_ttaRG0_20_yAll.root");
  GetCutG(0,2,1,"vAll/rooto/bound_pidRaw_f7_right_45_100_132_ttaRG0_20_yAll.root");
  GetCutG(0,2,2,"vAll/rooto/bound_pidRaw_f7_right_45_100_132_ttaRG20_40_yAll.root");
  GetCutG(0,2,3,"vAll/rooto/bound_pidRaw_f7_right_45_100_132_ttaRG40_60_yAll.root");
  GetCutG(0,2,4,"vAll/rooto/bound_pidRaw_f7_right_45_100_132_ttaRG40_60_yAll.root");

  // left and right
  GetCutG(1,0,0,"vAll/rooto/bound_pidRaw_f7_right_45_100_108_ttaRG0_20_yAll.root");
  GetCutG(1,0,1,"vAll/rooto/bound_pidRaw_f7_right_45_100_108_ttaRG0_20_yAll.root");
  GetCutG(1,0,2,"vAll/rooto/bound_pidRaw_f7_right_45_100_108_ttaRG20_40_yAll.root");
  GetCutG(1,0,3,"vAll/rooto/bound_pidRaw_f7_right_45_100_108_ttaRG40_60_yAll.root");
  GetCutG(1,0,4,"vAll/rooto/bound_pidRaw_f7_right_45_100_108_ttaRG40_60_yAll.root");
  GetCutG(0,0,0,"vAll/rooto/bound_pidRaw_f7_right_45_100_132_ttaRG0_20_yAll.root");
  GetCutG(0,0,1,"vAll/rooto/bound_pidRaw_f7_right_45_100_132_ttaRG0_20_yAll.root");
  GetCutG(0,0,2,"vAll/rooto/bound_pidRaw_f7_right_45_100_132_ttaRG20_40_yAll.root");
  GetCutG(0,0,3,"vAll/rooto/bound_pidRaw_f7_right_45_100_132_ttaRG40_60_yAll.root");
  GetCutG(0,0,4,"vAll/rooto/bound_pidRaw_f7_right_45_100_132_ttaRG40_60_yAll.root");

  auto getSelBoundHandCut = [](int iSys, int iParticle) {
    const char *conditionBoundL = "(-0.8<phi_cm&&phi_cm<0.4)";
    const char *conditionBoundR = "(-0.8>phi_cm||phi_cm>0.4)";
    const char *conditionBoundT0  = "(theta_lab*TMath::RadToDeg()>=0&&theta_lab*TMath::RadToDeg()<20)";
    const char *conditionBoundT20 = "(theta_lab*TMath::RadToDeg()>=20&&theta_lab*TMath::RadToDeg()<40)";
    const char *conditionBoundT40 = "(theta_lab*TMath::RadToDeg()>=40&&theta_lab*TMath::RadToDeg()<80)";
    const char *selectionBoundLT0  = Form("bound_%d_%d_%d_%d",iSys,kLeft, kTheta0, iParticle);
    const char *selectionBoundLT20 = Form("bound_%d_%d_%d_%d",iSys,kLeft, kTheta20,iParticle);
    const char *selectionBoundLT40 = Form("bound_%d_%d_%d_%d",iSys,kLeft, kTheta40,iParticle);
    const char *selectionBoundRT0  = Form("bound_%d_%d_%d_%d",iSys,kRight,kTheta0, iParticle);
    const char *selectionBoundRT20 = Form("bound_%d_%d_%d_%d",iSys,kRight,kTheta20,iParticle);
    const char *selectionBoundRT40 = Form("bound_%d_%d_%d_%d",iSys,kRight,kTheta40,iParticle);
    TCut selBoundHandCut = Form("%s*%s*%s+%s*%s*%s+%s*%s*%s+%s*%s*%s+%s*%s*%s+%s*%s*%s"
        ,conditionBoundL,conditionBoundT0 ,selectionBoundLT0
        ,conditionBoundL,conditionBoundT20,selectionBoundLT20
        ,conditionBoundL,conditionBoundT40,selectionBoundLT40
        ,conditionBoundR,conditionBoundT0 ,selectionBoundRT0
        ,conditionBoundR,conditionBoundT20,selectionBoundRT20
        ,conditionBoundR,conditionBoundT40,selectionBoundRT40);
    return selBoundHandCut;
  };

  ///////////////////////////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////////////////////////

  gStyle -> SetOptStat(0);

  for (auto iSys : selSysIdx)
  {
    auto sys = fSysBeams[iSys];
    const char *pidMetaName = Form("data2/Meta_Sn%dKanekoMult50.root",sys);
    const char *pidSigmaName = Form("data2/PIDSigma_Sn%dKanekoMult50.root",sys);
    auto filePIDMeta = new TFile(pidMetaName,"read");
    auto filePIDSigma = new TFile(pidSigmaName,"read");

    for (auto iParticle : fParticleIdx)
    {
      auto pdg = fParticlePDGs[iParticle];
      auto histMeta = (TH1F*) filePIDMeta -> Get(TString::Format("Distribution%d", pdg));
      auto f1Mean = static_cast<TF1*>(filePIDSigma -> Get(TString::Format("BEE%d", pdg)));
      f1Mean -> FixParameter(0, f1Mean -> GetParameter(0));
      f1Mean -> SetRange(0, 4500);
      auto f1Sigma = static_cast<TF1*>(filePIDSigma -> Get(TString::Format("PIDSigma%d", pdg)));
      f1Sigma -> SetRange(0, 4500);
      auto cutg = static_cast<TCutG*>(filePIDSigma -> Get(TString::Format("Limit%d", pdg)));

      histPIDMeta[iSys][iParticle] = histMeta;
      f1PIDMean[iSys][iParticle] = f1Mean;
      f1PIDSigma[iSys][iParticle] = f1Sigma;
      cutgPID[iSys][iParticle] = cutg;

      auto graph = new TGraph();
      auto graph2 = new TGraph();
      auto graph3 = new TGraph();
      for (auto mom=100.; mom<3000.; mom+=2)
      {
        auto dedx = f1PIDMean[iSys][iParticle] -> Eval(mom);
        if (cutg->IsInside(mom,dedx) && dedx<3000) {
          auto dedxSigma = f1PIDSigma[iSys][iParticle] -> Eval(mom);
          graph -> SetPoint(graph->GetN(), mom, dedx);
          graph2 -> SetPoint(graph2->GetN(), mom, dedx-dedxSigma);
          graph3 -> SetPoint(graph3->GetN(), mom, dedx+dedxSigma);
        }
      }

      graphPIDMean[iSys][iParticle] = graph;
      graphPIDMean[iSys][iParticle] -> SetName(Form("guideLine_sys%d_part%d",iSys,iParticle));
      graphPIDMean[iSys][iParticle] -> SetLineColor(kGray+2);
      graphPIDMean[iSys][iParticle] -> SetLineColor(kBlack);
      //graphPIDMean[iSys][iParticle] -> SetLineColor(kRed);
      //if (iSys==0) graphPIDMean[iSys][iParticle] -> SetLineColor(kBlack);
      graphPIDRange[iSys][iParticle][0] = graph2;
      graphPIDRange[iSys][iParticle][1] = graph3;
    }
  }

  for (auto iAna : fAnaIdx)
  {
    if (selAna>=0 && selAna!=iAna) continue;
    const char *anaFName = fAnaFNames[iAna];
    const char *spVersion = fAnaVersion[iAna];

    for (auto iMult : fMultIdx)
    {
      if (selMult>=0 && selMult!=iMult) continue;
      const char *multFName = fMultFNames[iMult];

      TH1D *histKeoaArray[3][fNumSyss][fNumCutThetas][fNumCutYPs][fNumParticles] = {0};
      TH1D *histTCMArray[3][fNumSyss][fNumCutThetas][fNumCutYPs][fNumParticles] = {0};

      for (auto iLR : selLRIdx)
      {
        if (selLR>=0 && selLR!=iLR) continue;
        const char *lrFName = fLRFNames[iLR];

        for (auto iSys : selSysIdx)
        {
          if (selSys>=0 && selSys!=iSys) continue;
          auto sys = fSysBeams[iSys];
          const char *sysTitle = fSysTitles[iSys];

          if (drawRawPID || setDrawRawPIDProj)
          {
            auto nameFileAll = Form("data2/%s/sys%d_%s_%s_%s.%s.ana.all.root",anaFName,sys,anaFName,lrFName,multFName,spVersion);
            cout_info << "File : " << nameFileAll << endl;

            auto treeAll = new TChain("all");
            treeAll -> Add(nameFileAll);

            auto iCutYP = 0;
            //for (auto iCutYP : selCutYPIdx)
            {
              //if (selCutYP>=0 && selCutYP!=iCutYP) continue;
              //TString stringCutYP0 = fCutYPValues[iCutYP].GetTitle();

              for (auto iCutTheta : selCutThetaIdx)
              {
                if (selCutTheta>=0 && selCutTheta!=iCutTheta) continue;
                TCut cutTheta = fCutThetaValues[iCutTheta];

                auto namePIDraw = makeName("pidRaw",iAna,iLR,iMult,iSys,iCutTheta,iCutYP);
                auto titlePIDraw = makeTitle("Raw PID",iAna,iLR,iMult,iSys,iCutTheta,iCutYP);

                auto histPID = (bnPoz*bndEdx).newHist(namePIDraw,titlePIDraw);
                //histPID -> SetMinimum(0.5);
                //histPID -> SetMaximum(800);
                project(treeAll,namePIDraw,"dedx:p_lab",cutTheta);

                TCanvas *cvsPIDraw = nullptr;
                if (drawRawPID) {
                  cvsPIDraw = makeCvs2(namePIDraw,1000,700);
                  //cvsPIDraw = makeCvs2(namePIDraw,950,850);
                  if (setDrawPIDLogz)
                    cvsPIDraw -> SetLogz();
                  histPID -> Draw("colz");
                  if (setDrawGuideLine) {
                    for (auto iParticle : fParticleIdx) {
                      graphPIDMean[iSys][iParticle] -> Draw("samel");
                      if (setDrawGuideLine132)
                        graphPIDMean[0][iParticle] -> Draw("samel");
                    }
                  }
                  if (setDrawHandCut) {
                    for (auto iParticle : fParticleIdx)
                      cutgAll[iSys][iLR][iCutTheta][iParticle] -> Draw("samel");
                  }

                }

                if (setDrawRawPIDProj)
                {
                  //auto binn = binning(histPID);
                  binning binn = bnPoz;
                  int dbin = bnPoz.fN/nixPIDProjX.fN;

                  for (auto iParticle : fParticleIdx) {
                    graphFitPIDMean[iSys][iParticle] = new TGraph();
                    graphFitPIDAmp[iSys][iParticle] = new TGraph();
                  }

                  for (auto iProj=0; iProj<nixPIDProjX.fN; ++iProj)
                  {
                    auto bin1 = (iProj)*dbin+1;
                    auto bin2 = (iProj+1)*dbin;
                    auto poz1 = binn.lowEdge(bin1);
                    auto poz2 = binn.highEdge(bin2);
                    auto pozC = (poz1 + poz2)/2.;

                    auto nameProj = makeName(Form("proj0_dedx_%d",iProj),iAna,iLR,iMult,iSys,iCutTheta,iCutYP);
                    if (setDrawRawPIDProj && setFitdEdx) nameProj = makeName(Form("proj1_dedx_%d",iProj),iAna,iLR,iMult,iSys,iCutTheta,iCutYP);
                    else if (setFitdEdx)                       nameProj = makeName(Form("proj2_dedx_%d",iProj),iAna,iLR,iMult,iSys,iCutTheta,iCutYP);
                    else                                       nameProj = makeName(Form("proj3_dedx_%d",iProj),iAna,iLR,iMult,iSys,iCutTheta,iCutYP);

                    auto histProj = (TH1D *) histPID -> ProjectionY(nameProj,bin1,bin2);
                    histProj -> Rebin(dbin);

                    TCanvas *cvsProj = nullptr;

                    if (nixPIDProjX.isInside(pozC))
                    {
                      TString titleProj = TString("p/Z=")+poz1+"-"+poz2+";dE/dx;";
                      histProj -> SetTitle(titleProj);

                      if (setDrawRawPIDProj) {
                        cvsProj = makeCvs(nameProj,1000,550);
                        cvsProj -> SetLogy();
                        histProj -> Draw();
                      }

                      double scaleAmp = 1.;
                      const char *nameProjFit = makeName(Form("projFit_dedx_%d",iProj),iAna,iLR,iMult,iSys,iCutTheta,iCutYP);
                      auto f1dEdxTotal = new TF1(nameProjFit,"gaus(0)+gaus(3)+gaus(6)+gaus(9)+gaus(12)",bndEdx.fMin,bndEdx.fMax);
                      f1dEdxTotal -> SetNpx(1000);

                      for (auto iParticle : fParticleIdx)
                      {
                        auto dedxRefAmp = histPIDMeta[iSys][iParticle] -> GetBinContent(histPIDMeta[iSys][iParticle] -> GetXaxis() -> FindBin(pozC));
                        auto dedxMean = f1PIDMean[iSys][iParticle] -> Eval(pozC);
                        auto dedxSigma = f1PIDSigma[iSys][iParticle] -> Eval(pozC);

                        double dedxAmp = dedxRefAmp;
                        if (iParticle==0) {
                          dedxAmp = histProj -> GetMaximum();
                          scaleAmp = dedxAmp / dedxRefAmp;
                        } else
                          dedxAmp = dedxRefAmp * scaleAmp;

                        const char *nameProjFitPart = makeName(Form("projFitdEdx0_%d",iProj),iAna,iLR,iMult,iSys,iCutTheta,iCutYP,iParticle);
                        auto f1dEdxParticle = new TF1(nameProjFitPart,"gaus(0)",bndEdx.fMin,bndEdx.fMax);
                        f1dEdxParticle -> SetLineColor(kGray);
                        f1dEdxParticle -> SetParameters(dedxAmp,dedxMean,dedxSigma);
                        f1dEdxParticle -> SetNpx(1000);
                        if (setDrawRawPIDProj) {
                          cvsProj -> cd();
                          f1dEdxParticle -> Draw("samel");
                        }


                        f1dEdxTotal -> SetParameter(0+3*iParticle,dedxAmp);
                        f1dEdxTotal -> SetParameter(1+3*iParticle,dedxMean);
                        f1dEdxTotal -> SetParameter(2+3*iParticle,dedxSigma);

                        if (iParticle==0&&pozC>1500) f1dEdxTotal -> FixParameter(0+3*iParticle,0);
                        if (iParticle==1&&pozC>2200) f1dEdxTotal -> FixParameter(0+3*iParticle,0);
                        if (iParticle==3&&pozC>1300) f1dEdxTotal -> FixParameter(0+3*iParticle,0);
                        if (iParticle==4&&pozC>1800) f1dEdxTotal -> FixParameter(0+3*iParticle,0);

                        if (iParticle==3&&pozC< 300) f1dEdxTotal -> FixParameter(0+3*iParticle,0);

                        if (setFitdEdx) {
                          if (scaleMaxFitdEdx>0) {
                            f1dEdxTotal -> SetParLimits(1+3*iParticle,dedxMean -scaleMaxFitdEdx*dedxMean ,dedxMean +scaleMaxFitdEdx*dedxMean );
                            f1dEdxTotal -> SetParLimits(2+3*iParticle,dedxSigma-scaleMaxFitdEdx*dedxSigma,dedxSigma+scaleMaxFitdEdx*dedxSigma);
                          }
                          else {
                            f1dEdxTotal -> FixParameter(1+3*iParticle,dedxMean);
                            f1dEdxTotal -> FixParameter(2+3*iParticle,dedxSigma);
                          }
                        }
                      }

                      if (setFitdEdx)
                        histProj -> Fit(f1dEdxTotal,"RQ0");

                      for (auto iParticle : fParticleIdx) {
                        const char *nameProjFitPart = makeName(Form("projFitdEdx_%d",iProj),iAna,iLR,iMult,iSys,iCutTheta,iCutYP,iParticle);
                        auto f1dEdxParticle = new TF1(nameProjFit,"gaus(0)",bndEdx.fMin,bndEdx.fMax);
                        f1dEdxParticle -> SetLineColor(kGray+1);
                        auto dedxAmp = f1dEdxTotal -> GetParameter(0+3*iParticle);
                        auto dedxMean = f1dEdxTotal -> GetParameter(1+3*iParticle);
                        auto dedxSigma = f1dEdxTotal -> GetParameter(2+3*iParticle);
                        f1dEdxParticle -> SetParameters(dedxAmp,dedxMean,dedxSigma);
                        f1dEdxParticle -> SetNpx(1000);
                        if (setDrawRawPIDProj) {
                          cvsProj -> cd();
                          f1dEdxParticle -> Draw("samel");
                        }

                        graphFitPIDMean[iSys][iParticle] -> SetPoint(graphFitPIDMean[iSys][iParticle]->GetN(),pozC,dedxMean);
                        if (dedxAmp<0) dedxAmp = 0;
                        auto graphAmp = graphFitPIDAmp[iSys][iParticle];
                        if (!isinf(dedxAmp)) {
                          graphAmp -> SetPoint(graphAmp->GetN(),pozC,dedxAmp);
                        }
                      }

                      if (setDrawRawPIDProj) {
                        cvsProj -> cd();
                        f1dEdxTotal -> Draw("samel");
                      }

                      if (drawRawPID)
                      {
                        cvsPIDraw -> cd();
                        auto line1 = new TLine(poz1,0,poz1,20);
                        line1 -> SetLineColor(kRed);
                        line1 -> Draw("samel");
                        auto line2 = new TLine(poz2,0,poz2,20);
                        line2 -> SetLineColor(kRed);
                        line2 -> Draw("samel");
                      }
                    }
                  }

                  if (setFitdEdx) {
                    for (auto iParticle : fParticleIdx) {
                      cvsPIDraw -> cd();
                      graphFitPIDMean[iSys][iParticle] -> SetLineColor(kGray+1);
                      graphFitPIDMean[iSys][iParticle] -> Draw("samel");
                    }

                    if (setDrawRawPIDProj && setFitdEdx) {
                      makeCvs(Form("refit_%s",namePIDraw) ,1000,550);
                      for (auto iParticle : fParticleIdx) {
                        auto graph = graphFitPIDAmp[iSys][iParticle];
                        setAtt(graph,iParticle);

                        if (iParticle==0) graph -> Draw("apl");
                        else graph -> Draw("samepl");
                      }
                    }
                  }
                }

              }
            }

            delete treeAll;
          }
        } // iSys

        //int numEventsInAna[fNumSyss] = {0};
        TH1D *histY0Array[fNumSyss][fNumCutThetas][fNumCutYPs][fNumParticles] = {0};
        TH1D *histPtoaArray[fNumSyss][fNumCutThetas][fNumCutYPs][fNumParticles] = {0};
        TH1D *histThetaArray[fNumSyss][fNumCutThetas][fNumCutYPs][fNumParticles] = {0};
        TH2D *histXYR21Array[fNumSyss][fNumCutThetas][fNumCutYPs][fNumParticles] = {0};
        TH2D *histXYR21Array2[fNumSyss][fNumCutThetas][fNumCutYPs][fNumParticles] = {0};

        for (auto iSys : selSysIdx)
        {
          if (selSys>=0 && selSys!=iSys) continue;

          auto sys = fSysBeams[iSys];
          const char *sysTitle = fSysTitles[iSys];
          TCut cutSys = fSysCut[iSys];

          //if (numEventsInAna[iSys]==0)
          {
            auto nameFileMult = Form("data2/%s/sys%d_%s_%s_%s.%s.ana.mult.root",anaFName,sys,anaFName,lrFName,multFName,spVersion);

            auto treeMult = new TChain("mult");
            treeMult -> Add(nameFileMult);

            numEventsInAna[iSys] = treeMult -> GetEntries("1");
            cout_info << "Event Multiplicity in " << sys << " = " << numEventsInAna[iSys] << "  (" << nameFileMult << ")" << endl;

            delete treeMult;
          }

          if (drawKT || drawYP || drawKY)
          {
            for (auto iCutYP : selCutYPIdx)
            {
              if (selCutYP>=0 && selCutYP!=iCutYP) continue;
              TString stringCutYP0 = fCutYPValues[iCutYP].GetTitle();

              auto iCutTheta = 0;
              TCut cutTheta = fCutThetaValues[iCutTheta];

              TCanvas *cvsKY = nullptr;
              TCanvas *cvsKT = nullptr;
              TCanvas *cvsYPTogether = nullptr;

              for (auto iParticle : fParticleIdx)
              {
                TCut cutParticlePoz = fParticlePozCut[iParticle];
                TString stringParticleSD = fSDParticleCut[fSelSDType][iParticle].GetTitle();
                stringParticleSD.ReplaceAll("SDVALUE",Form("%.1f",fSDValue));
                TCut cutParticleSD = stringParticleSD.Data();
                TString stringCutYP = stringCutYP0;
                stringCutYP.ReplaceAll("PARTICLEA",Form("%d",fParticleA[iParticle]));
                stringCutYP.ReplaceAll("PARTICLEZ",Form("%d",fParticleZ[iParticle]));
                stringCutYP.ReplaceAll("PARTICLEM",Form("%f",fParticleMass[iParticle]));
                TCut cutParticleYP = stringCutYP.Data();

                const char *nameParticle = fParticleNames[iParticle];
                const char *titleParticle = fParticleTitles[iParticle];
                auto nameFileParticle = Form("data2/%s/sys%d_%s_%s_%s.%s.ana.%s.root",anaFName,sys,anaFName,lrFName,multFName,spVersion,nameParticle);

                auto treeParticle = new TChain(nameParticle);
                treeParticle -> Add(nameFileParticle);

                if (drawKY)
                {
                  auto nameKYCorPart = makeName("kyCor",iAna,iLR,iMult,iSys,iCutTheta,iCutYP,iParticle);
                  auto titleKYCorPart = makeTitle("Cor.",iAna,iLR,iMult,iSys,iCutTheta,iCutYP,iParticle);

                  auto partz = fParticleZ[iParticle];
                  auto partm = fParticleMass[iParticle];
                  auto parta = fParticleA[iParticle];

                  TCut selection = cut0 * getWeighting(iSys,iParticle) * cutSys * cutParticleYP * cutTheta * cutParticlePoz * cutParticleSD;

                  auto histKYCorPart = (bnY0*bnKeoaCM).newHist(nameKYCorPart,titleKYCorPart);
                  const char *expression = Form("(sqrt((p_cm)*(p_cm)+%f*%f)-%f)/%d:fy_cm/(by_cm/2)",partm,partm,partm,parta);
                  //auto histKYCorPart = (bnY0*bnPoACM).newHist(nameKYCorPart,titleKYCorPart);
                  //const char *expression = Form("p_cm:fy_cm/(by_cm/2)");

                  project(treeParticle,nameKYCorPart,expression,selection);

                  if (cvsKY==nullptr) {
                    cvsKY = makeCvs2(nameKYCorPart,1545,873);
                    cvsDivide(cvsKY,3,2);
                  }

                  cvsKY -> cd(iParticle+1);
                  histKYCorPart -> Draw("colz");

                  /*
                  if (setDrawKTGrid)
                  {
                    for (auto ibinKeoa=0; ibinKeoa<=bnRKeoaCM.fN; ++ibinKeoa) {
                      auto binKeoa = ibinKeoa + 1;
                      auto line = new TLine(bnRKeoaCM.lowEdge(binKeoa),bnRTCM.lowEdge(),bnRKeoaCM.lowEdge(binKeoa),bnRTCM.highEdge());
                      line -> Draw("samel");
                    }

                    for (auto ibinTCM=0; ibinTCM<=bnRTCM.fN; ++ibinTCM) {
                      auto binTCM = ibinTCM + 1;
                      auto line = new TLine(bnRKeoaCM.lowEdge(),bnRTCM.lowEdge(binTCM),bnRKeoaCM.highEdge(),bnRTCM.lowEdge(binTCM));
                      line -> Draw("samel");
                    }

                    for (auto ibinKeoa=0; ibinKeoa<bnRKeoaCM.fN; ++ibinKeoa) {
                      for (auto ibinTCM=0; ibinTCM<bnRTCM.fN; ++ibinTCM) {
                        auto binKeoa = ibinKeoa + 1;
                        auto binTCM = ibinTCM + 1;
                        auto x = bnRKeoaCM.getCenter(binKeoa);
                        auto y = bnRTCM.getCenter(binTCM);
                        auto text = new TLatex(x,y,Form("%d,%d",ibinKeoa,ibinTCM));
                        text -> Draw();
                        text -> SetTextAlign(22);
                        text -> SetTextSize(0.03);
                      }
                    }
                  }
                  */

                  /*
                  if (setFitKE)
                  {
                    int nbinsRTCM = bnRTCM.fN*(bnTCM.fMax/bnRTCM.fMax);
                    int nbinsPerR = bnTCM.fN/nbinsRTCM;
                    bnRTCM.reset();
                    TH1D * histMain = nullptr;
                    double yMaxKT = 0;
                    auto nameKTPY0 = makeName("KTPYF",iAna,iLR,iMult,iSys,iCutTheta,iCutYP,iParticle);
                    auto cvsKTPY = makeCvs(nameKTPY0);
                    auto legendKTPY = new TLegend();
                    while (bnRTCM.next())
                    {
                      int idx = bnRTCM.fIdx;
                      int bin1 = 1 + (idx-1)*nbinsPerR;
                      int bin2 = idx*nbinsPerR;
                      auto nameKTPY = makeName(Form("KTPY%d",idx),iAna,iLR,iMult,iSys,iCutTheta,iCutYP,iParticle);
                      auto histKTPY = (TH1D *) histKTCorPart -> ProjectionX(nameKTPY,bin1,bin2);
                      setAtt(histKTPY,idx-1);
                      if (idx==1) { histMain = histKTPY; histKTPY -> Draw("hist"); }
                      else histKTPY -> Draw("samehist");
                      if (yMaxKT < histKTPY -> GetMaximum())
                        yMaxKT = histKTPY -> GetMaximum();
                      legendKTPY -> AddEntry(histKTPY,Form("%s=%.0f",bnRTCM.fTitle,bnRTCM.getValue()));
                    }
                    histMain -> SetMaximum(yMaxKT*1.1);
                    makeLegend(cvsKTPY,legendKTPY) -> Draw();
                  }
                  */
                }

                if (drawKT)
                {
                  auto nameKTCorPart = makeName("ktCor",iAna,iLR,iMult,iSys,iCutTheta,iCutYP,iParticle);
                  auto titleKTCorPart = makeTitle("Cor.",iAna,iLR,iMult,iSys,iCutTheta,iCutYP,iParticle);

                  auto histKTCorPart = (bnKeoaCM*bnTCM).newHist(nameKTCorPart,titleKTCorPart);
                  TCut selection = cut0 * getWeighting(iSys,iParticle) * cutSys * cutParticleYP * cutTheta * cutParticlePoz * cutParticleSD;
                  auto partz = fParticleZ[iParticle];
                  auto partm = fParticleMass[iParticle];
                  auto parta = fParticleA[iParticle];
                  const char *expression = Form("theta_cm*TMath::RadToDeg():(sqrt((p_cm)*(p_cm)+%f*%f)-%f)/%d",partm,partm,partm,parta);
                  //const char *expression = Form("theta_cm*TMath::RadToDeg():(sqrt((p_cm*%d)*(p_cm*%d)+%f*%f)-%f)/%d",partz,partz,partm,partm,partm,1);
                  project(treeParticle,nameKTCorPart,expression,selection);

                  if (cvsKT==nullptr) {
                    cvsKT = makeCvs2(nameKTCorPart,1545,873);
                    cvsDivide(cvsKT,3,2);
                  }

                  cvsKT -> cd(iParticle+1);
                  histKTCorPart -> Draw("colz");

                  if (setDrawKTGrid)
                  {
                    for (auto ibinKeoa=0; ibinKeoa<=bnRKeoaCM.fN; ++ibinKeoa) {
                      auto binKeoa = ibinKeoa + 1;
                      auto line = new TLine(bnRKeoaCM.lowEdge(binKeoa),bnRTCM.lowEdge(),bnRKeoaCM.lowEdge(binKeoa),bnRTCM.highEdge());
                      line -> Draw("samel");
                    }

                    for (auto ibinTCM=0; ibinTCM<=bnRTCM.fN; ++ibinTCM) {
                      auto binTCM = ibinTCM + 1;
                      auto line = new TLine(bnRKeoaCM.lowEdge(),bnRTCM.lowEdge(binTCM),bnRKeoaCM.highEdge(),bnRTCM.lowEdge(binTCM));
                      line -> Draw("samel");
                    }

                    for (auto ibinKeoa=0; ibinKeoa<bnRKeoaCM.fN; ++ibinKeoa) {
                      for (auto ibinTCM=0; ibinTCM<bnRTCM.fN; ++ibinTCM) {
                        auto binKeoa = ibinKeoa + 1;
                        auto binTCM = ibinTCM + 1;
                        auto x = bnRKeoaCM.getCenter(binKeoa);
                        auto y = bnRTCM.getCenter(binTCM);
                        auto text = new TLatex(x,y,Form("%d,%d",ibinKeoa,ibinTCM));
                        text -> Draw();
                        text -> SetTextAlign(22);
                        text -> SetTextSize(0.03);
                      }
                    }
                  }

                  if (setFitKE)
                  {
                    int nbinsRTCM = bnRTCM.fN*(bnTCM.fMax/bnRTCM.fMax);
                    int nbinsPerR = bnTCM.fN/nbinsRTCM;
                    bnRTCM.reset();
                    TH1D * histMain = nullptr;
                    double yMaxKT = 0;
                    auto nameKTPY0 = makeName("KTPYF",iAna,iLR,iMult,iSys,iCutTheta,iCutYP,iParticle);
                    auto cvsKTPY = makeCvs(nameKTPY0);
                    auto legendKTPY = new TLegend();
                    while (bnRTCM.next())
                    {
                      int idx = bnRTCM.fIdx;
                      int bin1 = 1 + (idx-1)*nbinsPerR;
                      int bin2 = idx*nbinsPerR;
                      auto nameKTPY = makeName(Form("KTPY%d",idx),iAna,iLR,iMult,iSys,iCutTheta,iCutYP,iParticle);
                      auto histKTPY = (TH1D *) histKTCorPart -> ProjectionX(nameKTPY,bin1,bin2);
                      setAtt(histKTPY,idx-1);
                      if (idx==1) { histMain = histKTPY; histKTPY -> Draw("hist"); }
                      else histKTPY -> Draw("samehist");
                      if (yMaxKT < histKTPY -> GetMaximum())
                        yMaxKT = histKTPY -> GetMaximum();
                      legendKTPY -> AddEntry(histKTPY,Form("%s=%.0f",bnRTCM.fTitle,bnRTCM.getValue()));
                    }
                    histMain -> SetMaximum(yMaxKT*1.1);
                    makeLegend(cvsKTPY,legendKTPY) -> Draw();
                  }
                }

                if (drawYP)
                {
                  auto nameYPPart = makeName("yp",iAna,iLR,iMult,iSys,iCutTheta,iCutYP,iParticle);
                  auto titleYPPart = makeTitle("",iAna,iLR,iMult,iSys,iCutTheta,iCutYP,iParticle);

                  auto histYPPart  = (bnY0*bnPtoa).newHist(nameYPPart,titleYPPart);
                  TCut selection = cut0 * getWeighting(iSys,iParticle) * cutSys * cutParticleYP * cutTheta * cutParticlePoz * cutParticleSD;
                  project(treeParticle,nameYPPart,Form("pt_cm/%d:fy_cm/(by_cm/2)",fParticleA[iParticle]),selection);

                  TCanvas *cvsYP = nullptr;
                  if (setDrawYPTogether) {
                    if (cvsYPTogether==nullptr) {
                      auto nameYPPartT = makeName("yp",iAna,iLR,iMult,iSys,iCutTheta,iCutYP,-1);
                      cvsYPTogether = makeCvs2(nameYPPartT,1545,873);
                      cvsDivide(cvsYPTogether,3,2);
                    }
                    cvsYP = (TCanvas *) cvsYPTogether -> cd(iParticle+1);
                  }
                  else
                    cvsYP = makeCvs2(nameYPPart);
                  histYPPart -> Draw("colz");
                  //histYPPart -> Draw("contz");

                  if (setDrawYPGrid) {
                    for (auto ibinY0=0; ibinY0<=bnRY0.fN; ++ibinY0) {
                      auto binY0 = ibinY0 + 1;
                      auto line = new TLine(bnRY0.lowEdge(binY0),bnRPtoa.lowEdge(),bnRY0.lowEdge(binY0),bnRPtoa.highEdge());
                      line -> Draw("samel");
                    }

                    for (auto ibinPtoa=0; ibinPtoa<=bnRPtoa.fN; ++ibinPtoa) {
                      auto binPtoa = ibinPtoa + 1;
                      auto line = new TLine(bnRY0.lowEdge(),bnRPtoa.lowEdge(binPtoa),bnRY0.highEdge(),bnRPtoa.lowEdge(binPtoa));
                      line -> Draw("samel");
                    }

                    for (auto ibinY0=0; ibinY0<bnRY0.fN; ++ibinY0) {
                      for (auto ibinPtoa=0; ibinPtoa<bnRPtoa.fN; ++ibinPtoa) {
                        auto binY0 = ibinY0 + 1;
                        auto binPtoa = ibinPtoa + 1;
                        auto x = bnRY0.getCenter(binY0);
                        auto y = bnRPtoa.getCenter(binPtoa);
                        auto text = new TLatex(x,y,Form("%d,%d",ibinY0,ibinPtoa));
                        text -> Draw();
                        text -> SetTextAlign(22);
                        text -> SetTextSize(0.03);
                      }
                    }
                  }

                  if (setDrawYPTheta0) {
                    if (setUseCM) {
                      //for (auto iCutTheta0 : {2,4,6})
                      //for (auto iCutTheta0 : {2,4,6,8})
                      for (auto iCutTheta0 : {4})
                      {
                        int theta0;
                        if (iCutTheta0==0) theta0 = 20;
                        if (iCutTheta0==1) theta0 = 40;
                        if (iCutTheta0==2) theta0 = 60;
                        if (iCutTheta0==3) theta0 = 80;
                        if (iCutTheta0==4) theta0 = 45;
                        auto nameYPPart2 = makeName(Form("yptta%d",theta0),iAna,iLR,iMult,iSys,iCutTheta,iCutYP,iParticle);

                        auto nameYPPart3 = makeName(Form("fityptta%d",theta0),iAna,iLR,iMult,iSys,iCutTheta,iCutYP,iParticle);
                        TF1 *fitPYTT0 = new TF1(nameYPPart3,"pol3",bnY0.fMin,bnY0.fMax);
                             if (iParticle==0&&iCutTheta0==0) fitPYTT0 -> SetParameters(0.0915153, 126.913, -0.424865, 4.6369);
                        else if (iParticle==0&&iCutTheta0==1) fitPYTT0 -> SetParameters(-1.07311, 303.558, -24.6869, 38.1692);
                        else if (iParticle==0&&iCutTheta0==2) fitPYTT0 -> SetParameters(-3.77352, 654.414, -164.357, 309.893);
                        else if (iParticle==0&&iCutTheta0==3) fitPYTT0 -> SetParameters(38.8146, 1209.69, 5337.41, -6494.49);
                        else if (iParticle==1&&iCutTheta0==0) fitPYTT0 -> SetParameters(-0.94412, 130.742, -5.10702, 5.9206);
                        else if (iParticle==1&&iCutTheta0==1) fitPYTT0 -> SetParameters(0.476904, 292.275, -5.44058, 27.1222);
                        else if (iParticle==1&&iCutTheta0==2) fitPYTT0 -> SetParameters(0.171048, 619.522, -78.4747, 235.322);
                        else if (iParticle==1&&iCutTheta0==3) fitPYTT0 -> SetParameters(28.9653, 1476.94, 3282.18, -3035.38);
                        else if (iParticle==2&&iCutTheta0==0) fitPYTT0 -> SetParameters(-1.73135, 134.42, -9.37208, 7.51902);
                        else if (iParticle==2&&iCutTheta0==1) fitPYTT0 -> SetParameters(1.13857, 286.933, 7.29332, 18.2867);
                        else if (iParticle==2&&iCutTheta0==2) fitPYTT0 -> SetParameters(2.31779, 593.928, 4.92044, 152.886);
                        else if (iParticle==2&&iCutTheta0==3) fitPYTT0 -> SetParameters(22.2969, 1691.53, 1088.54, 2441.7);
                        else if (iParticle==3&&iCutTheta0==0) fitPYTT0 -> SetParameters(2.27599, 114.306, 19.944, -5.33683);
                        else if (iParticle==3&&iCutTheta0==1) fitPYTT0 -> SetParameters(-0.491297, 296.275, -8.51022, 27.2079);
                        else if (iParticle==3&&iCutTheta0==2) fitPYTT0 -> SetParameters(3.3637, 587.639, 17.6572, 146.537);
                        else if (iParticle==3&&iCutTheta0==3) fitPYTT0 -> SetParameters(11.3119, 1998.41, -922.841, 6156.03);
                        else if (iParticle==4&&iCutTheta0==0) fitPYTT0 -> SetParameters(0.146599, 127.271, -3.01256, 5.65602);
                        else if (iParticle==4&&iCutTheta0==1) fitPYTT0 -> SetParameters(0.515944, 291.465, -5.51404, 25.7714);
                        else if (iParticle==4&&iCutTheta0==2) fitPYTT0 -> SetParameters(-1.0333, 625.622, -103.089, 249.629);
                        else if (iParticle==4&&iCutTheta0==3) fitPYTT0 -> SetParameters(16.5863, 1869.43, -633.24, 5342.01);

                        else if (iParticle==0&&iCutTheta0==4) fitPYTT0 -> SetParameters(-2.67251, 372.468, -51.1756, 66.962);
                        else if (iParticle==1&&iCutTheta0==4) fitPYTT0 -> SetParameters(0.074858, 352.396, -15.6081, 46.5267);
                        else if (iParticle==2&&iCutTheta0==4) fitPYTT0 -> SetParameters(0.688515, 350.105, -14.442, 46.9903);
                        else if (iParticle==3&&iCutTheta0==4) fitPYTT0 -> SetParameters(1.11233, 350.236, -14.3401, 46.3099);
                        else if (iParticle==4&&iCutTheta0==4) fitPYTT0 -> SetParameters(-2.17183, 367.625, -50.6081, 66.114);

                        else {
                          auto histYPThetaPart = (bnY0*bnPtoa).newHist(nameYPPart2,titleYPPart);
                          TCut cutTheta0 = Form("theta_cm>= %d*TMath::DegToRad()&&theta_cm<%d*TMath::DegToRad()",theta0-1,theta0+1);
                          selection = cut0 * getWeighting(iSys,iParticle) * cutSys * cutParticleYP * cutTheta * cutTheta0 * cutParticlePoz * cutParticleSD;
                          project(treeParticle,nameYPPart2,Form("pt_cm/%d:fy_cm/(by_cm/2)",fParticleA[iParticle]),selection,0);
                          //if (iCutTheta0==2) histYPThetaPart -> Draw("colz");
                          histYPThetaPart -> Draw("samecol");

                          histYPThetaPart -> Fit(fitPYTT0,"RQ0");
                          cout << "if (iParticle==" << iParticle << "&&iCutTheta0==" << iCutTheta0 << ") fitPYTT0 -> SetParameters("
                            << fitPYTT0 -> GetParameter(0)
                            << ", " << fitPYTT0 -> GetParameter(1)
                            << ", " << fitPYTT0 -> GetParameter(2)
                            << ", " << fitPYTT0 -> GetParameter(3) << ");" << endl;
                        }
                        fitPYTT0 -> SetLineStyle(2);
                        fitPYTT0 -> SetLineColor(kSpring-6);
                        fitPYTT0 -> SetRange(0,fitPYTT0 -> GetX(bnPtoa.fMax));
                        fitPYTT0 -> Draw("samel");

                        if (setDrawYPText) {
                          auto xAt = bnY0.fMax;
                          auto yAt = fitPYTT0 -> Eval(xAt);

                          if (yAt<bnPtoa.fMax) {
                            TLatex *text = new TLatex(xAt,yAt,Form("#theta_{CM}=%d#circ",theta0));
                            text -> SetTextColor(kSpring-6);
                            text -> SetTextFont(132);
                            text -> SetTextAlign(31);
                            text -> SetTextSize(.03);
                            text -> Draw("samel");
                          }
                          else {
                            yAt = bnPtoa.fMax;
                            xAt = fitPYTT0 -> GetX(yAt);
                            TLatex *text = new TLatex(xAt,yAt,Form("#theta_{Lab}=%d#circ",theta0));
                            text -> SetTextColor(kSpring-6);
                            text -> SetTextFont(132);
                            text -> SetTextAlign(23);
                            text -> SetTextSize(.03);
                            text -> Draw("samel");
                          }
                        }
                      }
                    }
                    else {
                      //left
                      if (0)
                      for (auto iCutTheta0 : {2,4,6})
                      {
                        int theta0 = 10*iCutTheta0;
                        auto nameYPPart3 = makeName(Form("fitypttaL%d",theta0),iAna,iLR,iMult,iSys,iCutTheta,iCutYP,iParticle);
                        TF1 *fitPYTT0 = new TF1(nameYPPart3,"pol3",bnY0.fMin,bnY0.fMax);
                        if (iParticle==0&&iCutTheta0==2) fitPYTT0 -> SetParameters(153.764, 159.678, 15.8344, 10.4356);
                        if (iParticle==0&&iCutTheta0==4) fitPYTT0 -> SetParameters(359.115, 419.475, 133.317, 74.8092);
                        if (iParticle==0&&iCutTheta0==6) fitPYTT0 -> SetParameters(957.397, 1485.75, 805.723, 249.215);
                        if (iParticle==1&&iCutTheta0==2) fitPYTT0 -> SetParameters(153.125, 159.625, 17.2045, 7.16013);
                        if (iParticle==1&&iCutTheta0==4) fitPYTT0 -> SetParameters(358.199, 418.164, 119.323, 53.0897);
                        if (iParticle==1&&iCutTheta0==6) fitPYTT0 -> SetParameters(925.684, 1350.42, 613.342, 161.085);
                        if (iParticle==2&&iCutTheta0==2) fitPYTT0 -> SetParameters(153.191, 159.935, 16.7452, 6.05917);
                        if (iParticle==2&&iCutTheta0==4) fitPYTT0 -> SetParameters(357.397, 415.437, 123.999, 64.0913);
                        if (iParticle==2&&iCutTheta0==6) fitPYTT0 -> SetParameters(1002.93, 1812.56, 1456.12, 642.88);
                        if (iParticle==3&&iCutTheta0==2) fitPYTT0 -> SetParameters(153.429, 159.803, 15.2702, 8.4324);
                        if (iParticle==3&&iCutTheta0==4) fitPYTT0 -> SetParameters(358.102, 413.971, 113.696, 54.2012);
                        if (iParticle==3&&iCutTheta0==6) fitPYTT0 -> SetParameters(940.464, 1413.94, 654.907, 122.295);
                        if (iParticle==4&&iCutTheta0==2) fitPYTT0 -> SetParameters(152.517, 158.431, 15.5089, 7.31443);
                        if (iParticle==4&&iCutTheta0==4) fitPYTT0 -> SetParameters(356.386, 413.661, 109.878, 42.6489);
                        if (iParticle==4&&iCutTheta0==6) fitPYTT0 -> SetParameters(847.918, 971.275, -41.59, -227.475);
                        fitPYTT0 -> SetLineStyle(2);
                        fitPYTT0 -> SetLineColor(kSpring-6);

                        double xMin = 0;
                        nameYPPart3 = makeName(Form("fitypket"),iAna,iLR,iMult,iSys,iCutTheta,iCutYP,iParticle);
                        auto fitPoz = new TF1(nameYPPart3,"+[3]*sqrt(1-((x-[0])/[2])*((x-[0])/[2]))+[1]",-1,2);
                        if (iParticle==0) fitPoz -> SetParameters( -1, 0, .25, 100);
                        else if (iParticle==1) fitPoz -> SetParameters( -1, 0, .25, 100);
                        else if (iParticle==2) fitPoz -> SetParameters( -1, 0, 0.35, 138);
                        else if (iParticle==3) fitPoz -> SetParameters( -1, 0, 0.72, 260);
                        else if (iParticle==4) fitPoz -> SetParameters( -1, 0, 0.55, 200);
                        fitPoz -> SetRange(-1,fitPoz->GetParameter(0)+fitPoz->GetParameter(2));
                        double xxx2 = fitPoz->GetParameter(2);
                        for (double xxx=-0.9; xxx<xxx2; xxx+=0.01)
                        {
                          //cout << xxx << " " << fitPYTT0 -> Eval(xxx) << " " <<  fitPoz -> Eval(xxx) << endl;
                          if (fitPYTT0 -> Eval(xxx) > fitPoz -> Eval(xxx)) {
                            xMin = xxx;
                            break;
                          }
                        }

                        if (setDrawYPText) {
                          auto xAt = bnY0.fMax;
                          auto yAt = fitPYTT0 -> Eval(xAt);
                          if (yAt<bnPtoa.fMax) {
                            fitPYTT0 -> SetRange(xMin,xAt);
                            TLatex *text = new TLatex(xAt,yAt,Form("#theta_{L}=%d#circ",theta0));
                            text -> SetTextColor(kSpring-6);
                            text -> SetTextFont(132);
                            text -> SetTextAlign(31);
                            text -> SetTextSize(.03);
                            text -> Draw("samel");
                          }
                          else {
                            yAt = 795;
                            xAt = fitPYTT0 -> GetX(yAt);
                            fitPYTT0 -> SetRange(xMin,xAt);
                            TLatex *text = new TLatex(xAt,yAt,Form("#theta_{L}=%d#circ",theta0));
                            text -> SetTextColor(kSpring-6);
                            text -> SetTextFont(132);
                            text -> SetTextAlign(33);
                            text -> SetTextSize(.03);
                            text -> Draw("samel");
                          }
                        }
                        fitPYTT0 -> Draw("samel");
                      }

                      //right
                      if (1)
                      for (auto iCutTheta0 : {2,4,6})
                      {
                        int theta0 = 10*iCutTheta0;
                        auto nameYPPart3 = makeName(Form("fitypttaR%d",theta0),iAna,iLR,iMult,iSys,iCutTheta,iCutYP,iParticle);
                        TF1 *fitPYTT0 = new TF1(nameYPPart3,"pol3",bnY0.fMin,bnY0.fMax);
                        if (iParticle==0&&iCutTheta0==2) fitPYTT0 -> SetParameters(118.019, 121.496, 10.2241, 5.87606);
                        if (iParticle==0&&iCutTheta0==4) fitPYTT0 -> SetParameters(297.32, 331.255, 79.7603, 47.3711);
                        if (iParticle==0&&iCutTheta0==6) fitPYTT0 -> SetParameters(722.575, 1061.13, 585.208, 234.343);
                        if (iParticle==1&&iCutTheta0==2) fitPYTT0 -> SetParameters(117.573, 121.143, 10.8934, 4.71274);
                        if (iParticle==1&&iCutTheta0==4) fitPYTT0 -> SetParameters(296.376, 331.322, 74.7611, 34.9029);
                        if (iParticle==1&&iCutTheta0==6) fitPYTT0 -> SetParameters(708.596, 999.862, 483.362, 175.898);
                        if (iParticle==2&&iCutTheta0==2) fitPYTT0 -> SetParameters(117.642, 121.275, 10.9714, 3.74235);
                        if (iParticle==2&&iCutTheta0==4) fitPYTT0 -> SetParameters(296.231, 331.227, 74.1494, 33.1961);
                        if (iParticle==2&&iCutTheta0==6) fitPYTT0 -> SetParameters(699.746, 957.293, 409.222, 131.047);
                        if (iParticle==3&&iCutTheta0==2) fitPYTT0 -> SetParameters(117.588, 121.761, 10.9632, 4.1624);
                        if (iParticle==3&&iCutTheta0==4) fitPYTT0 -> SetParameters(296.624, 330.273, 71.8979, 37.5958);
                        if (iParticle==3&&iCutTheta0==6) fitPYTT0 -> SetParameters(708.005, 984.894, 411.327, 85.049);
                        if (iParticle==4&&iCutTheta0==2) fitPYTT0 -> SetParameters(117.062, 120.505, 10.0153, 4.21134);
                        if (iParticle==4&&iCutTheta0==4) fitPYTT0 -> SetParameters(294.529, 327.531, 68.039, 29.1065);
                        if (iParticle==4&&iCutTheta0==6) fitPYTT0 -> SetParameters(657.24, 738.002, 10.851, -113.578);
                        fitPYTT0 -> SetLineStyle(2);
                        fitPYTT0 -> SetLineColor(kSpring-6);

                        double xMin = 0;
                        nameYPPart3 = makeName(Form("fitypket"),iAna,iLR,iMult,iSys,iCutTheta,iCutYP,iParticle);
                        auto fitPoz = new TF1(nameYPPart3,"+[3]*sqrt(1-((x-[0])/[2])*((x-[0])/[2]))+[1]",-1,2);
                        if (iParticle==0) fitPoz -> SetParameters( -1, 0, .25, 100);
                        else if (iParticle==1) fitPoz -> SetParameters( -1, 0, .25, 100);
                        else if (iParticle==2) fitPoz -> SetParameters( -1, 0, 0.35, 138);
                        else if (iParticle==3) fitPoz -> SetParameters( -1, 0, 0.72, 260);
                        else if (iParticle==4) fitPoz -> SetParameters( -1, 0, 0.55, 200);
                        fitPoz -> SetRange(-1,fitPoz->GetParameter(0)+fitPoz->GetParameter(2));
                        double xxx2 = fitPoz->GetParameter(2);
                        for (double xxx=-0.9; xxx<xxx2; xxx+=0.01)
                        {
                          //cout << xxx << " " << fitPYTT0 -> Eval(xxx) << " " <<  fitPoz -> Eval(xxx) << endl;
                          if (fitPYTT0 -> Eval(xxx) > fitPoz -> Eval(xxx)) {
                            xMin = xxx;
                            break;
                          }
                        }

                        if (setDrawYPText) {
                          auto xAt = bnY0.fMax;
                          auto yAt = fitPYTT0 -> Eval(xAt);
                          if (yAt<bnPtoa.fMax) {
                            fitPYTT0 -> SetRange(xMin,xAt);
                            TLatex *text = new TLatex(xAt,yAt,Form("#theta_{R}=%d#circ",theta0));
                            text -> SetTextColor(kSpring-6);
                            text -> SetTextFont(132);
                            text -> SetTextAlign(31);
                            text -> SetTextSize(.03);
                            text -> Draw("samel");
                          }
                          else {
                            yAt = 795;
                            xAt = fitPYTT0 -> GetX(yAt);
                            fitPYTT0 -> SetRange(xMin,xAt);
                            TLatex *text = new TLatex(xAt,yAt,Form("#theta_{R}=%d#circ",theta0));
                            text -> SetTextColor(kSpring-6);
                            text -> SetTextFont(132);
                            text -> SetTextAlign(13);
                            text -> SetTextSize(.03);
                            text -> Draw("samel");
                          }
                        }

                        fitPYTT0 -> Draw("samel");
                      }
                    }
                  }

                  if (setDrawYPPoz)
                  {
                    for (auto iCutPoz : {2,5,10,15})
                    {
                      double iPoz;
                      iPoz = iCutPoz *2;
                      auto partz = fParticleZ[iParticle];
                      auto partm = fParticleMass[iParticle];
                      auto parta = fParticleA[iParticle];
                      const char *exprPoz = "p_lab";
                      auto dPoz = (iPoz/2.)*2;
                      TCut cutPoz = Form("%s>%f&&%s<%f",exprPoz,iPoz*100-dPoz,exprPoz,iPoz*100+dPoz);

                      auto nameYPPart2 = makeName(Form("yppoz%d",iCutPoz),iAna,iLR,iMult,iSys,iCutTheta,iCutYP,iParticle);

                      auto histYPPozPart = (bnY0*bnPtoa).newHist(nameYPPart2,titleYPPart);
                      selection = cut0 * getWeighting(iSys,iParticle) * cutSys * cutParticleYP * cutTheta * cutPoz * cutParticlePoz * cutParticleSD;
                      project(treeParticle,nameYPPart2,Form("pt_cm/%d:fy_cm/(by_cm/2)",fParticleA[iParticle]),selection);

                      auto nameYPPart3 = makeName(Form("fitypke"),iAna,iLR,iMult,iSys,iCutTheta,iCutYP,iParticle);
                      auto fitPoz = new TF1(nameYPPart3,"+[3]*sqrt(1-((x-[0])/[2])*((x-[0])/[2]))+[1]",-1,2);
                           if (iParticle==0 && iCutPoz==2)  fitPoz -> SetParameters( -0.997, 0, 1.1, 400);
                      else if (iParticle==0 && iCutPoz==5)  fitPoz -> SetParameters( -1, -50, 2.5, 950);
                      else if (iParticle==1 && iCutPoz==2)  fitPoz -> SetParameters( -1, 0, 0.55, 200);
                      else if (iParticle==1 && iCutPoz==5)  fitPoz -> SetParameters( -0.97, 0, 1.32, 500);
                      else if (iParticle==1 && iCutPoz==10) fitPoz -> SetParameters( -0.997, -50, 2.48, 950);
                      else if (iParticle==2 && iCutPoz==2)  fitPoz -> SetParameters( -0.997, 0, 0.35, 138);
                      else if (iParticle==2 && iCutPoz==5)  fitPoz -> SetParameters( -1, 0, 0.9, 340);
                      else if (iParticle==2 && iCutPoz==10) fitPoz -> SetParameters( -1, -50, 1.75, 710);
                      else if (iParticle==2 && iCutPoz==15) fitPoz -> SetParameters( -0.997, -50, 2.48, 950);
                      else if (iParticle==3 && iCutPoz==2)  fitPoz -> SetParameters( -1, 0, 0.72, 260);
                      else if (iParticle==3 && iCutPoz==5)  fitPoz -> SetParameters( -1, -50, 1.75, 710);
                      else if (iParticle==4 && iCutPoz==2)  fitPoz -> SetParameters( -1, 0, 0.55, 200);
                      else if (iParticle==4 && iCutPoz==5)  fitPoz -> SetParameters( -1, 0, 1.36, 500);
                      else if (iParticle==4 && iCutPoz==10) fitPoz -> SetParameters( -1, 0, 2.5, 850);
                      else
                        continue;
                      histYPPozPart -> Fit(fitPoz,"Q0");
                      fitPoz -> SetRange(-1,fitPoz->GetParameter(0)+fitPoz->GetParameter(2));
                      //fitPoz -> SetNpx(1000);
                      fitPoz -> SetLineStyle(2);
                      fitPoz -> SetLineColor(kRed+1);
                      fitPoz -> Draw("samel");


                      if (setDrawYPText) {
                        auto xAt = fitPoz->GetX(0);
                        //TLatex *text = new TLatex(xAt,10,Form("p_{Lab}/Z=%d",int(iPoz*100)));
                        TLatex *text = new TLatex(xAt,10,Form("p/Z=%d",int(iPoz*100)));
                        text -> SetTextFont(132);
                        text -> SetTextAngle(90);
                        text -> SetTextAlign(13);
                        text -> SetTextColor(kRed+1);
                        text -> SetTextSize(.03);
                        text -> Draw("samel");
                      }
                    }
                  }

                  if (setDrawYPKE)
                  {
                    if (setUseCM)
                    {
                      for (auto iCutKE : {0,1,2})
                      //for (auto iCutKE : {2})
                      {
                        double iKE;
                        if (1) {//iParticle==0||iParticle==1||iParticle==2) {
                          if (iCutKE==0) iKE = .5;
                          else if (iCutKE==1) iKE = 2;
                          else if (iCutKE==2) iKE = .1;
                        }
                        else {
                          if (iCutKE==0) iKE = .5;
                          else if (iCutKE==1) iKE = 2;
                          else if (iCutKE==2) iKE = 5;
                        }
                        auto partz = fParticleZ[iParticle];
                        auto partm = fParticleMass[iParticle];
                        auto parta = fParticleA[iParticle];
                        const char *exprKE = Form("(sqrt((p_cm)*(p_cm)+%f*%f)-%f)/%d",partm,partm,partm,parta);
                        //const char *exprKE = Form("(sqrt((p_cm)*(p_cm)+%f*%f)-%d)/%d",partm,partm,0,parta);
                        auto dKE = (iKE/2.)*5;
                        //TCut cutKE = Form("%s>%f&&%s<%f",exprKE,iKE*100-dKE,exprKE,iKE*100+dKE);
                        TCut cutKE = Form("abs(%s-%f)<%f",exprKE,iKE*100,dKE);

                        auto nameYPPart2 = makeName(Form("ypke%d",iCutKE),iAna,iLR,iMult,iSys,iCutTheta,iCutYP,iParticle);

                        auto nameYPPart3 = makeName(Form("fitypke"),iAna,iLR,iMult,iSys,iCutTheta,iCutYP,iParticle);
                        auto fitKE = new TF1(nameYPPart3,"+[1]*sqrt(1-((x)/[0])*((x)/[0]))",-1,2);
                        fitKE -> SetParameters(0.4,130);

                        if (iCutKE==0) fitKE -> SetParameters(0.893751,304.401);
                        else if (iCutKE==1) fitKE -> SetParameters(1.68403,640);
                        else if (iCutKE==2) fitKE -> SetParameters(0.406264,132.714);
                        /*
                        if (iParticle==0&&iCutKE==0) fitKE -> SetParameters(0.893751,304.401);
                        else if (iParticle==0&&iCutKE==1) fitKE -> SetParameters(1.68403,640);
                        else if (iParticle==0&&iCutKE==2) fitKE -> SetParameters(0.406264,132.714);
                        else if (iParticle==1&&iCutKE==0) fitKE -> SetParameters(0.93125,299.581);
                        else if (iParticle==1&&iCutKE==1) fitKE -> SetParameters(1.68144,635.28);
                        else if (iParticle==1&&iCutKE==2) fitKE -> SetParameters(0.418753,130.503);
                        else if (iParticle==2&&iCutKE==0) fitKE -> SetParameters(0.931254,303.145);
                        else if (iParticle==2&&iCutKE==1) fitKE -> SetParameters(1.79376,623.697);
                        else if (iParticle==2&&iCutKE==2) fitKE -> SetParameters(0.406257,134.492);

                        else if (iParticle==3&&iCutKE==2) fitKE -> SetParameters(0.443752,155.708);
                        else if (iParticle==3&&iCutKE==0) fitKE -> SetParameters(0.931251,320.457);
                        else if (iParticle==3&&iCutKE==1) fitKE -> SetParameters(1.68403,640);//fitKE -> SetParameters(1.53125,540);

                        else if (iParticle==4&&iCutKE==2) fitKE -> SetParameters(0.456258,152.717);
                        else if (iParticle==4&&iCutKE==0) fitKE -> SetParameters(0.931259,316.565);
                        else if (iParticle==4&&iCutKE==1) fitKE -> SetParameters(1.68403,640);//fitKE -> SetParameters(1.5,540);
                        */
                        else
                        {
                          auto histYPKEPart = (bnY0*bnPtoa).newHist(nameYPPart2,titleYPPart);
                          selection = cut0 * getWeighting(iSys,iParticle) * cutSys * cutParticleYP * cutTheta * cutKE * cutParticlePoz * cutParticleSD;
                          project(treeParticle,nameYPPart2,Form("pt_cm/%d:fy_cm/(by_cm/2)",fParticleA[iParticle]),selection);
                          auto graph = new TGraphErrors();
                          bnY0.reset();
                          while (bnY0.next()) {
                            auto histPj = (TH1D *) histYPKEPart->ProjectionY(Form("%s_p%d",histYPKEPart->GetName(),bnY0.fIdx),bnY0.fIdx,bnY0.fIdx);
                            auto yRep = histPj -> GetMean();
                            auto yRepError = histPj -> GetStdDev();
                            auto xRep = bnY0.getValue();
                            if (yRep>0) {
                              auto gpidx = graph -> GetN();
                              graph -> SetPoint(gpidx,xRep,yRep);
                              graph -> SetPointError(gpidx,.5*bnY0.getW(),3*bnPtoa.getW());
                              //graph -> SetPointError(gpidx,.5*bnY0.getW(),yRepError);
                            }
                          }
                          //if (iCutKE==0) histYPKEPart  -> Draw("col");
                          //histYPKEPart  -> Draw("col");
                          graph -> Draw("samep");
                          fitKE -> SetParLimits(0,fitKE->GetParameter(0)-.1,fitKE->GetParameter(0)+.1);
                          fitKE -> SetParLimits(1,fitKE->GetParameter(1)-50,fitKE->GetParameter(1)+50);
                          fitKE -> SetNpx(100);
                          //graph -> Fit(fitKE,"QB0");
                          cout << "if (iParticle=="<<iParticle<<"&&iCuKE=="<<iCutKE<<") fitKE -> SetParameters("<<fitKE->GetParameter(0) << "," << fitKE->GetParameter(1) << ");"<<endl;
                        }

                        fitKE -> SetRange(-fitKE->GetParameter(0),fitKE->GetParameter(0));
                        //fitKE -> SetNpx(1000);
                        fitKE -> SetLineStyle(2);
                        fitKE -> SetLineColor(kBlack);
                        fitKE -> Draw("samel");

                        if (setDrawYPText) {
                          auto xAt = fitKE -> GetParameter(0);
                          TLatex *text = new TLatex(xAt,10,Form("KE/A=%d",int(iKE*100)));
                          text -> SetTextColor(kBlack);
                          text -> SetTextFont(132);
                          text -> SetTextAngle(90);
                          text -> SetTextAlign(13);
                          if (iParticle==4) text -> SetTextAlign(11);
                          text -> SetTextSize(.03);
                          text -> Draw("samel");
                        }
                      }
                    }
                    else
                      for (auto iCutKE : {0,1,2})
                      {
                        double iKE;
                        if (iParticle==0) {
                          if (iCutKE==0) iKE = .5;
                          else if (iCutKE==1) iKE = 2;
                          else if (iCutKE==2) iKE = 5;
                        }
                        else {
                          if (iCutKE==0) iKE = 2;
                          else if (iCutKE==1) iKE = 5;
                          else if (iCutKE==2) iKE = 10;
                        }
                        auto partz = fParticleZ[iParticle];
                        auto partm = fParticleMass[iParticle];
                        auto parta = fParticleA[iParticle];
                        //const char *exprKE = Form("(sqrt((p_lab*%d)*(p_lab*%d)+%f*%f)-%f)",partz,partz,partm,partm,partm);
                        const char *exprKE = Form("(sqrt((p_lab*%d)*(p_lab*%d)+%f*%f)-%f)/%d",partz,partz,partm,partm,partm,parta);
                        //const char *exprKE = Form("(sqrt((p_cm*%d)*(p_cm*%d)+%f*%f)-%f)",partz,partz,partm,partm,partm);
                        auto dKE = (iKE/2.)*5;
                        TCut cutKE = Form("%s>%f&&%s<%f",exprKE,iKE*100-dKE,exprKE,iKE*100+dKE);

                        auto nameYPPart2 = makeName(Form("ypke%d",iCutKE),iAna,iLR,iMult,iSys,iCutTheta,iCutYP,iParticle);

                        auto histYPKEPart = (bnY0*bnPtoa).newHist(nameYPPart2,titleYPPart);
                        //if (iCutKE==0) histYPKEPart  -> Draw("col");
                        //else histYPKEPart  -> Draw("colsame");
                        selection = cut0 * getWeighting(iSys,iParticle) * cutSys * cutParticleYP * cutTheta * cutKE * cutParticlePoz * cutParticleSD;
                        project(treeParticle,nameYPPart2,Form("pt_cm/%d:fy_cm/(by_cm/2)",fParticleA[iParticle]),selection);

                        auto nameYPPart3 = makeName(Form("fitypke"),iAna,iLR,iMult,iSys,iCutTheta,iCutYP,iParticle);
                        auto fitKE = new TF1(nameYPPart3,"+[3]*sqrt(1-((x-[0])/[2])*((x-[0])/[2]))+[1]",-1,2);
                        if (iParticle==0 && iCutKE==2) fitKE -> SetParameters(-1.3958, 0, 3, 1000);
                        else if (iParticle==0 && iCutKE==1) fitKE -> SetParameters(-1.0967, -50, 1.8, 700);
                        else if (iParticle==0 && iCutKE==0) fitKE -> SetParameters(-1.00013, 0.397535, 0.864614, 299.648);
                        else if (iParticle==1 && iCutKE==2) fitKE -> SetParameters(-0.997, -100, 2.6, 1100);
                        else if (iParticle==1 && iCutKE==1) fitKE -> SetParameters(-1.78481, -98.1845, 2.77286, 863.169);
                        else if (iParticle==1 && iCutKE==0) fitKE -> SetParameters(-1, 0, 1.2, 450);
                        else if (iParticle==2 && iCutKE==2) fitKE -> SetParameters(-1, -100, 2.2, 930);
                        else if (iParticle==2 && iCutKE==1) fitKE -> SetParameters(-1.0074, -102.414, 1.65051, 668.946);
                        else if (iParticle==2 && iCutKE==0) fitKE -> SetParameters(-1, 0, 1, 350);
                        else if (iParticle==3 && iCutKE==2) fitKE -> SetParameters(-1, -100, 2.2, 930);
                        else if (iParticle==3 && iCutKE==1) fitKE -> SetParameters(-1.0088, -97.9979, 1.62158, 671.03);
                        else if (iParticle==3 && iCutKE==0) fitKE -> SetParameters(-1, 0, 1, 350);
                        else if (iParticle==4 && iCutKE==2) fitKE -> SetParameters(-1, 0, 1.9, 690);
                        else if (iParticle==4 && iCutKE==1) fitKE -> SetParameters(-1, -50, 1.35, 540);
                        else if (iParticle==4 && iCutKE==0) fitKE -> SetParameters(-1.00195, 0.436524, 0.879741, 300.285);

                        histYPKEPart -> Fit(fitKE,"Q0");
                        fitKE -> SetRange(-1,fitKE->GetParameter(0)+fitKE->GetParameter(2));
                        //fitKE -> SetNpx(1000);
                        fitKE -> SetLineStyle(2);
                        fitKE -> SetLineColor(kBlack);
                        fitKE -> Draw("samel");

                        if (setDrawYPText) {
                          auto xAt = fitKE->GetX(0);
                          TLatex *text = new TLatex(xAt,10,Form("KE=%d",int(iKE*100)));
                          text -> SetTextColor(kBlack);
                          text -> SetTextFont(132);
                          text -> SetTextAngle(90);
                          text -> SetTextAlign(13);
                          if (iParticle==4) text -> SetTextAlign(11);
                          text -> SetTextSize(.03);
                          text -> Draw("samel");
                        }
                      }
                  }

                }

                delete treeParticle;
              }
            }
          }


          if (drawCorEKE || drawCorPID || drawTLabR21 || drawPtoaR21 || drawDistKeoa || drawDistTCM || drawTCMR21 || drawKeoaR21 || drawY0R21 || drawNZR21)
          {
            for (auto iCutYP : selCutYPIdx)
            {
              if (selCutYP>=0 && selCutYP!=iCutYP) continue;
              TString stringCutYP0 = fCutYPValues[iCutYP].GetTitle();

              for (auto iCutTheta : selCutThetaIdx)
              {
                if (selCutTheta>=0 && selCutTheta!=iCutTheta) continue;
                TCut cutTheta = fCutThetaValues[iCutTheta];

                TH2D *histPIDCor = nullptr;
                TH2D *histEKECor = nullptr;
                vector<TH2D *> histEKECorArray;// = nullptr;

                for (auto iParticle : fParticleIdx)
                {
                  TCut cutParticlePoz = fParticlePozCut[iParticle];
                  TString stringParticleSD = fSDParticleCut[fSelSDType][iParticle].GetTitle();
                  stringParticleSD.ReplaceAll("SDVALUE",Form("%.1f",fSDValue));
                  TCut cutParticleSD = stringParticleSD.Data();
                  TString stringCutYP = stringCutYP0;
                  stringCutYP.ReplaceAll("PARTICLEA",Form("%d",fParticleA[iParticle]));
                  stringCutYP.ReplaceAll("PARTICLEZ",Form("%d",fParticleZ[iParticle]));
                  stringCutYP.ReplaceAll("PARTICLEM",Form("%f",fParticleMass[iParticle]));
                  TCut cutParticleYP = stringCutYP.Data();

                  const char *nameParticle = fParticleNames[iParticle];
                  const char *titleParticle = fParticleTitles[iParticle];

                  auto nameFileParticle = Form("data2/%s/sys%d_%s_%s_%s.%s.ana.%s.root",anaFName,sys,anaFName,lrFName,multFName,spVersion,nameParticle);

                  auto treeParticle = new TChain(nameParticle);
                  treeParticle -> Add(nameFileParticle);

                  if (drawCorPID)
                  {
                    auto namePIDCorPart = makeName("pidCor",iAna,iLR,iMult,iSys,iCutTheta,iCutYP,iParticle);
                    auto titlePIDCorPart = makeTitle("Cor.",iAna,iLR,iMult,iSys,iCutTheta,iCutYP,iParticle);

                    auto histPIDCorPart = (bnPoz*bndEdx).newHist(namePIDCorPart,titlePIDCorPart);
                    histPIDCorPart -> SetMinimum(0.5);
                    histPIDCorPart -> SetMaximum(800);
                    TCut selection = cut0 * cutSys * cutParticleYP * cutTheta * getWeighting(iSys,iParticle) * cutParticlePoz * cutParticleSD;

                    if (setDrawCorPIDInRaw)
                      selection = cut0 * cutSys * cutParticleYP * cutTheta * cutParticlePoz * cutParticleSD;

                    project(treeParticle,namePIDCorPart,"dedx:p_lab",selection);

                    if (setTestCorPIDJustPDT) {
                      selection = cut0 * cutSys * cutParticleYP * cutTheta * cutParticlePoz * cutParticleSD * Form("1/%f",histPIDCorPart->Integral());
                      project(treeParticle,namePIDCorPart,"dedx:p_lab",selection);
                    }

                    if (histPIDCor==nullptr) {
                      auto namePIDCor = makeName("pidCor",iAna,iLR,iMult,iSys,iCutTheta,iCutYP);
                      histPIDCor = (TH2D *) histPIDCorPart -> Clone(namePIDCor);
                      histPIDCor -> SetTitle(Form("%s;p/Z (MeV/c);dE/dx;",titlePIDCorPart));
                    }
                    else {
                      if (setTestCorPIDJustPDT) {
                        if (iParticle==0||iParticle==1||iParticle==2)
                          histPIDCor -> Add(histPIDCorPart);
                      }
                      else
                        histPIDCor -> Add(histPIDCorPart);
                    }
                  }

                  if (drawCorEKE)
                  {
                    auto nameEKECorPart = makeName("ekeCor",iAna,iLR,iMult,iSys,iCutTheta,iCutYP,iParticle);
                    auto titleEKECorPart = makeTitle("Cor.",iAna,iLR,iMult,iSys,iCutTheta,iCutYP,iParticle);

                    auto histEKECorPart = (bnKeoa2*bndEdx).newHist(nameEKECorPart,titleEKECorPart);
                    TCut selection = cut0 * getWeighting(iSys,iParticle) * cutSys * cutParticleYP * cutTheta * cutParticlePoz * cutParticleSD;
                    auto partz = fParticleZ[iParticle];
                    auto partm = fParticleMass[iParticle];
                    auto parta = fParticleA[iParticle];
                    const char *expression = Form("dedx:(sqrt((p_lab*%d)*(p_lab*%d)+%f*%f)-%f)",partz,partz,partm,partm,partm);
                    project(treeParticle,nameEKECorPart,expression,selection);

                    if (histEKECor==nullptr) {
                      auto nameEKECor = makeName("ekeCor",iAna,iLR,iMult,iSys,iCutTheta,iCutYP);
                      histEKECor = (TH2D *) histEKECorPart -> Clone(nameEKECor);
                      histEKECor -> SetTitle(Form("%s;KE_{Lab} (MeV);dE/dx;",titleEKECorPart));
                    }
                    else
                      histEKECor -> Add(histEKECorPart);

                    histEKECorArray.push_back(histEKECorPart);
                  }

                  if (drawTLabR21)
                  {
                    auto nameThetaPart = makeName("theta",iAna,iLR,iMult,iSys,iCutTheta,iCutYP,iParticle);
                    auto titleThetaPart = makeTitle("Cor.",iAna,iLR,iMult,iSys,iCutTheta,iCutYP,iParticle);

                    auto histTheta = bnRTLab.newHist(nameThetaPart,titleThetaPart);
                    TCut selection = cut0 * getWeighting(iSys,iParticle) * cutSys * cutParticleYP * cutTheta  * cutParticlePoz * cutParticleSD;
                    if (fUseHandCut) selection = cut0 * getWeighting(iSys,iParticle,1) * getSelBoundHandCut(iSys,iParticle) * cutSys * cutParticleYP * cutTheta * cutParticlePoz * cutParticleSD;

                    project(treeParticle,nameThetaPart,("theta_cm*TMath::RadToDeg()"),selection);

                    histThetaArray[iSys][iCutTheta][iCutYP][iParticle] = histTheta;
                  }

                  if (drawPtoaR21)
                  {
                    auto namePtoaPart = makeName("ptoa",iAna,iLR,iMult,iSys,iCutTheta,iCutYP,iParticle);
                    auto titlePtoaPart = makeTitle("Cor.",iAna,iLR,iMult,iSys,iCutTheta,iCutYP,iParticle);

                    auto histPtoa = bnRPtoa.newHist(namePtoaPart,titlePtoaPart);
                    TCut selection = cut0 * getWeighting(iSys,iParticle) * cutSys * cutParticleYP * cutTheta  * cutParticlePoz * cutParticleSD;
                    if (fUseHandCut) selection = cut0 * getWeighting(iSys,iParticle,1) * getSelBoundHandCut(iSys,iParticle) * cutSys * cutParticleYP * cutTheta * cutParticlePoz * cutParticleSD;

                    project(treeParticle,namePtoaPart,Form("pt_cm/%d",fParticleA[iParticle]),selection);

                    histPtoaArray[iSys][iCutTheta][iCutYP][iParticle] = histPtoa;
                  }

                  if (drawDistKeoa || drawKeoaR21)
                  {
                    auto nameKeoaPart = makeName("keoa",iAna,iLR,iMult,iSys,iCutTheta,iCutYP,iParticle);
                    auto titleKeoaPart = makeTitle("Cor.",iAna,iLR,iMult,iSys,iCutTheta,iCutYP,iParticle);

                    auto histKeoa = bnRKeoaCM.newHist(nameKeoaPart,titleKeoaPart);
                    TCut selection = cut0 * getWeighting(iSys,iParticle) * cutSys * cutParticleYP * cutTheta * cutParticlePoz * cutParticleSD * TCut(Form("theta_cm*TMath::RadToDeg()<%f",bnRTCM.fMax));
                    if (fUseHandCut) selection = cut0 * getWeighting(iSys,iParticle,1) * getSelBoundHandCut(iSys,iParticle) * cutSys * cutParticleYP * cutTheta * cutParticlePoz * cutParticleSD;
                    auto partz = fParticleZ[iParticle];
                    auto partm = fParticleMass[iParticle];
                    auto parta = fParticleA[iParticle];
                    //const char *expression = Form("(sqrt((p_lab*%d)*(p_lab*%d)+%f*%f)-%f)/%d",partz,partz,partm,partm,partm,parta);
                    //const char *expression = Form("(sqrt((p_cm*%d)*(p_cm*%d)+%f*%f)-%f)/%d",partz,partz,partm,partm,partm,parta);
                    const char *expression = Form("(sqrt((p_cm)*(p_cm)+%f*%f)-%f)/%d",partm,partm,partm,parta);

                    project(treeParticle,nameKeoaPart,expression,selection);

                    histKeoaArray[iLR][iSys][iCutTheta][iCutYP][iParticle] = histKeoa;
                  }

                  if (drawDistTCM||drawTCMR21)
                  {
                    auto nameTCMPart = makeName("tcm2",iAna,iLR,iMult,iSys,iCutTheta,iCutYP,iParticle);
                    auto titleTCMPart = makeTitle("Cor.",iAna,iLR,iMult,iSys,iCutTheta,iCutYP,iParticle);

                    TH1D *histTCM = nullptr;
                    if (drawDistTCM && !drawTCMR21)
                      histTCM = bnTCM.newHist(nameTCMPart,titleTCMPart);
                    else
                      histTCM = bnRTCM.newHist(nameTCMPart,titleTCMPart);
                    TCut selection = cut0 * getWeighting(iSys,iParticle) * cutSys * cutParticleYP * cutTheta * cutParticlePoz * cutParticleSD;
                    if (fUseHandCut) selection = cut0 * getWeighting(iSys,iParticle,1) * getSelBoundHandCut(iSys,iParticle) * cutSys * cutParticleYP * cutTheta * cutParticlePoz * cutParticleSD;
                    const char *expression = "theta_cm*TMath::RadToDeg()";

                    project(treeParticle,nameTCMPart,expression,selection);

                    histTCMArray[iLR][iSys][iCutTheta][iCutYP][iParticle] = histTCM;
                  }
                  
                  if (drawY0R21)
                  {
                    auto nameY0R21Part = makeName("y0",iAna,iLR,iMult,iSys,iCutTheta,iCutYP,iParticle);
                    auto titleY0R21Part = makeTitle("Cor.",iAna,iLR,iMult,iSys,iCutTheta,iCutYP,iParticle);

                    auto histY0R21 = new TH1D(nameY0R21Part,Form("%s;y_{0};",titleY0R21Part),bnRY0.fN,bnRY0.fMin,bnRY0.fMax);
                    TCut selection = cut0 * getWeighting(iSys,iParticle) * cutSys * cutParticleYP * cutTheta  * cutParticlePoz * cutParticleSD;
                    if (fUseHandCut) selection = cut0 * getWeighting(iSys,iParticle,1) * getSelBoundHandCut(iSys,iParticle) * cutSys * cutParticleYP * cutTheta * cutParticlePoz * cutParticleSD;

                    project(treeParticle,nameY0R21Part,"fy_cm/(by_cm/2)",selection);

                    histY0Array[iSys][iCutTheta][iCutYP][iParticle] = histY0R21;
                  }

                  if (drawNZR21)
                  {
                    auto bnX = bnRY0;
                    auto bnY = bnRPtoa;

                    auto partz = fParticleZ[iParticle];
                    auto partm = fParticleMass[iParticle];
                    auto parta = fParticleA[iParticle];
                    const char *expressionNZ = Form("pt_cm/%d:fy_cm/(by_cm/2)",fParticleA[iParticle]);

                    if (setDrawNZInKT) {
                      expressionNZ = Form("theta_cm*TMath::RadToDeg():(sqrt((p_cm)*(p_cm)+%f*%f)-%f)/%d",partm,partm,partm,parta);
                      bnX = bnRKeoaCM;
                      bnY = bnRTCM;
                    }

                    auto nameYPR21Part = makeName("ypR21",iAna,iLR,iMult,iSys,iCutTheta,iCutYP,iParticle);
                    auto titleYPR21Part = makeTitle("Cor.",iAna,iLR,iMult,iSys,iCutTheta,iCutYP,iParticle);
                    TCut selection  = cut0 * getWeighting(iSys,iParticle) * cutSys * cutParticleYP * cutTheta * cutParticlePoz * cutParticleSD;
                    auto histXYR21Part = (bnX*bnY).newHist(nameYPR21Part,titleYPR21Part);
                    project(treeParticle,nameYPR21Part,expressionNZ,selection,histXYR21Part);

                    auto nameYPR21Part2 = makeName("ypR21_2",iAna,iLR,iMult,iSys,iCutTheta,iCutYP,iParticle);
                    TCut selection2 = cut1 * getWeighting(iSys,iParticle,1) * cutSys * cutParticleYP * cutTheta * cutParticlePoz * cutParticleSD;
                    auto histXYR21Part2 = (bnX*bnY).newHist(nameYPR21Part2,titleYPR21Part);
                    project(treeParticle,nameYPR21Part2,expressionNZ,selection2,histXYR21Part2);

                    bnX.reset();
                    bnY.reset();
                    while (bnX.next()) {
                      while (bnY.next()) {
                        auto value1 = histXYR21Part -> GetBinContent(bnX.fIdx,bnY.fIdx);
                        auto value2 = histXYR21Part2 -> GetBinContent(bnX.fIdx,bnY.fIdx);
                        auto error2 = abs(value1 - value2);
                        auto error = error2;
                      }
                    }

                    histXYR21Array[iSys][iCutTheta][iCutYP][iParticle] = histXYR21Part;
                    histXYR21Array2[iSys][iCutTheta][iCutYP][iParticle] = histXYR21Part2;
                  }

                  delete treeParticle;
                }

                if (drawCorPID)
                {
                  auto namePIDCor = makeName("pidCor",iAna,iLR,iMult,iSys,iCutTheta,iCutYP);
                  auto cvsPIDCor = makeCvs2(namePIDCor,1000,700);
                  if (setDrawPIDLogz)
                    cvsPIDCor -> SetLogz();
                  if (bnPIDz.fN>0) histPIDCor -> GetZaxis() -> SetRangeUser(bnPIDz.fMin,bnPIDz.fMax);
                  histPIDCor -> Draw("colz");
                  if (setDrawGuideLine) {
                    for (auto iParticle : fParticleIdx) {
                      graphPIDMean[iSys][iParticle] -> Draw("samel");
                      if (setDrawGuideLine132)
                        graphPIDMean[0][iParticle] -> Draw("samel");
                    }
                  }
                }

                if (drawCorEKE)
                {
                  auto nameEKECor = makeName("ekeCor",iAna,iLR,iMult,iSys,iCutTheta,iCutYP);
                  auto cvsEKECor = makeCvs2(nameEKECor,1500,1000);
                  cvsDivide(cvsEKECor,3,2);

                  for (auto i=0; i<5; ++i) {
                    auto cvs = cvsEKECor -> cd(i+1);
                    histEKECorArray.at(i) -> Draw("colz");
                    cvs -> SetLogz();
                  }

                  auto cvs = cvsEKECor -> cd(6);
                  cvs-> SetLogz();
                  histEKECor -> Draw("colz");
                }

                if (drawPtoa)
                {
                  const char *namePtoa = makeName("ptoaAll",iAna,iLR,iMult,iSys,iCutTheta,iCutYP);
                  auto cvsR21 = makeCvs(namePtoa);

                  double histMax = 0;
                  for (auto iParticle : fParticleIdx) {
                    auto max0 = histPtoaArray[iSys][iCutTheta][iCutYP][iParticle] -> GetMaximum();
                    if (histMax < max0)
                      histMax = max0;
                  }
                  for (auto iParticle : fParticleIdx) {
                    auto hist = histPtoaArray[iSys][iCutTheta][iCutYP][iParticle];
                    hist -> SetMaximum(histMax*1.05);
                    hist -> SetMinimum(0);
                    setAtt(hist,iParticle);
                    if (iParticle==0) hist -> Draw("plhist");
                    else hist -> Draw("sameplhist");
                  }
                }

              }
            }
          }
        }


        if (drawTLabR21)
        {
          for (auto iCutYP : selCutYPIdx)
          {
            if (selCutYP>=0 && selCutYP!=iCutYP) continue;

            auto iCutTheta = kThetaAll;
            {

              for (auto iComb : selSysCombIdx)
              {
                if (selSysComb>=0 && selSysComb!=iComb) continue;

                auto iSys2 = fSysCombIdx[iComb][0];
                auto iSys1 = fSysCombIdx[iComb][1];
                auto iSysComb = 100 + iSys2*10 + iSys1;

                const char *nameR21 = makeName("r21_theta",iAna,iLR,iMult,iSysComb,iCutTheta,iCutYP);
                auto titleR21 = makeTitle("Cor.",iAna,iLR,iMult,iSysComb,iCutTheta,iCutYP);

                auto cvs = makeCvs(nameR21);

                auto hist = (bnRTLab*bnR21).newHist(nameR21,titleR21);
                hist -> Draw();

                auto legend = new TLegend();
                for (auto iParticle : fParticleIdx)
                {
                  const char *nameParticle = fParticleNames[iParticle];

                  auto hist1 = histThetaArray[iSys2][iCutTheta][iCutYP][iParticle];
                  auto hist2 = histThetaArray[iSys1][iCutTheta][iCutYP][iParticle];

                  const char *nameR21Part = Form("%s_%s",nameR21,nameParticle);

                  auto hist0 = (TH1D *) hist1 -> Clone(nameR21Part);
                  hist0 -> Divide(hist2);
                  auto graph = new TGraph();
                  for (auto bin=1; bin<=bnRTLab.fN; ++bin)
                  {
                    if (hist1->GetBinContent(bin)<0.04||hist2->GetBinContent(bin)<0.04) {
                      hist0 -> SetBinContent(bin,0);
                    }
                    else
                      graph -> SetPoint(graph -> GetN(), hist0 -> GetXaxis() -> GetBinCenter(bin), hist0 -> GetBinContent(bin));
                  }

                  graph -> SetMarkerSize(1.3);
                  setAtt(graph,iParticle);

                  graph -> Draw("samepl");
                  legend -> AddEntry(graph,fParticleNames[iParticle],"pl");
                }
                makeLegend(cvs,legend) -> Draw();
              }

              /*
              if (setDrawTemp)
              {
                const char *nameTemp = makeName("temp_theta",iAna,iLR,iMult,kSysAll,iCutTheta,iCutYP);
                auto titleTemp = makeTitle("",iAna,iLR,iMult,kSysAll,iCutTheta,iCutYP);
                auto cvsTemp = makeCvs(nameTemp);
                auto hist = (bnRTLab*bnTemp).newHist(nameTemp,titleTemp);
                hist -> Draw();

                auto legend = new TLegend();
                for (auto iSys : selSysIdx)
                {
                  auto histd = histThetaArray[iSys][iCutTheta][iCutYP][kD];
                  auto hist4 = histThetaArray[iSys][iCutTheta][iCutYP][kHe4];
                  auto histt = histThetaArray[iSys][iCutTheta][iCutYP][kT];
                  auto hist3 = histThetaArray[iSys][iCutTheta][iCutYP][kHe3];

                  const char *nameHistTemp = Form("%s_sys%d_t",nameTemp,iSys);
                  auto histTemp = (TH1D *) histd -> Clone(nameHistTemp);
                  histTemp -> Multiply(hist4);
                  histTemp -> Divide(histt);
                  histTemp -> Divide(hist3);
                  setAtt(histTemp,iSys);

                  binning bnTempFinal(histTemp);
                  while (bnTempFinal.next())
                  {
                    auto value = histTemp -> GetBinContent(bnTempFinal.fIdx);
                    value = 14.3/TMath::Log(1.59*value);
                    histTemp -> SetBinContent(bnTempFinal.fIdx,value);
                  }
                  histTemp -> Draw("samehistpl");
                  legend -> AddEntry(histTemp ,Form("%d",fSysBeams[iSys]),"pl");
                }
                makeLegend(cvsTemp,legend,"lt") -> Draw();

              }
              */

            }
          }
        }

        if (drawPtoaR21)
        {
          for (auto iCutYP : selCutYPIdx)
          {
            if (selCutYP>=0 && selCutYP!=iCutYP) continue;

            for (auto iCutTheta : selCutThetaIdx)
            {
              if (selCutTheta>=0 && selCutTheta!=iCutTheta) continue;

              for (auto iComb : selSysCombIdx)
              {
                if (selSysComb>=0 && selSysComb!=iComb) continue;

                auto iSys2 = fSysCombIdx[iComb][0];
                auto iSys1 = fSysCombIdx[iComb][1];
                auto iSysComb = 100 + iSys2*10 + iSys1;

                const char *nameR21 = makeName("r21_ptoa",iAna,iLR,iMult,iSysComb,iCutTheta,iCutYP);
                auto titleR21 = makeTitle("Cor.",iAna,iLR,iMult,iSysComb,iCutTheta,iCutYP);

                auto cvs = makeCvs(nameR21);

                auto hist = (bnRPtoa*bnR21).newHist(nameR21,titleR21);
                hist -> Draw();

                auto legend = new TLegend();
                for (auto iParticle : fParticleIdx)
                {
                  const char *nameParticle = fParticleNames[iParticle];

                  auto hist1 = histPtoaArray[iSys2][iCutTheta][iCutYP][iParticle];
                  auto hist2 = histPtoaArray[iSys1][iCutTheta][iCutYP][iParticle];

                  const char *nameR21Part = Form("%s_%s",nameR21,nameParticle);

                  auto hist0 = (TH1D *) hist1 -> Clone(nameR21Part);
                  hist0 -> Divide(hist2);
                  auto graph = new TGraph();
                  //for (auto bin=1; bin<=8; ++bin)
                  for (auto bin=1; bin<=bnRPtoa.fN; ++bin)
                  {
                    if (hist1->GetBinContent(bin)<0.04||hist2->GetBinContent(bin)<0.04) {
                      hist0 -> SetBinContent(bin,0);
                    }
                    else
                      graph -> SetPoint(graph -> GetN(), hist0 -> GetXaxis() -> GetBinCenter(bin), hist0 -> GetBinContent(bin));
                  }

                  graph -> SetMarkerSize(1.3);
                  setAtt(graph,iParticle);

                  graph -> Draw("samepl");
                  legend -> AddEntry(graph,fParticleNames[iParticle],"pl");
                }
                makeLegend(cvs,legend) -> Draw();
              }

              if (setDrawTemp)
              {
                const char *nameTemp = makeName("temp_ptoa",iAna,iLR,iMult,kSysAll,iCutTheta,iCutYP);
                auto titleTemp = makeTitle("",iAna,iLR,iMult,kSysAll,iCutTheta,iCutYP);
                auto cvsTemp = makeCvs(nameTemp);
                auto hist = (bnRPtoa*bnTemp).newHist(nameTemp,titleTemp);
                hist -> Draw();

                auto legend = new TLegend();
                for (auto iSys : selSysIdx)
                {
                  auto histd = histPtoaArray[iSys][iCutTheta][iCutYP][kD];
                  auto hist4 = histPtoaArray[iSys][iCutTheta][iCutYP][kHe4];
                  auto histt = histPtoaArray[iSys][iCutTheta][iCutYP][kT];
                  auto hist3 = histPtoaArray[iSys][iCutTheta][iCutYP][kHe3];

                  const char *nameHistTemp = Form("%s_sys%d_t",nameTemp,iSys);
                  auto histTemp = (TH1D *) histd -> Clone(nameHistTemp);
                  histTemp -> Multiply(hist4);
                  histTemp -> Divide(histt);
                  histTemp -> Divide(hist3);
                  setAtt(histTemp,iSys);

                  //binning bnTempFinal(histTemp);
                  binning bnTempFinal = bnTemp;//(histTemp);
                  while (bnTempFinal.next())
                  {
                    auto value = histTemp -> GetBinContent(bnTempFinal.fIdx);
                    value = 14.3/TMath::Log(1.59*value);
                    histTemp -> SetBinContent(bnTempFinal.fIdx,value);
                  }
                  histTemp -> Draw("samehistpl");
                  legend -> AddEntry(histTemp ,Form("%d",fSysBeams[iSys]),"pl");
                }
                makeLegend(cvsTemp,legend,"lt") -> Draw();

              }

              if (setDrawPN)
              {
                /*
                const char *namePN = makeName("psuedoN_ptoa",iAna,iLR,iMult,kSysAll,iCutTheta,iCutYP);
                auto titlePN = makeTitle("",iAna,iLR,iMult,kSysAll,iCutTheta,iCutYP);
                auto cvsPN = makeCvs(namePN);
                auto hist = (bnRPtoa*bnPN).newHist(namePN,titlePN);
                hist -> Draw();

                auto legend = new TLegend();
                for (auto iSys : selSysIdx)
                {
                  auto histp = histPtoaArray[iSys][iCutTheta][iCutYP][kP];
                  auto histt = histPtoaArray[iSys][iCutTheta][iCutYP][kT];
                  auto hist3 = histPtoaArray[iSys][iCutTheta][iCutYP][kHe3];

                  const char *nameHistPN = Form("%s_sys%d_pn",namePN,iSys);
                  auto histPN = (TH1D *) histt -> Clone(nameHistPN);
                  histPN -> Multiply(histp);
                  histPN -> Divide(hist3);
                  setAtt(histPN,iSys);

                  binning bnPNFinal(histPN);
                  while (bnPNFinal.next())
                  {
                    auto value = histPN -> GetBinContent(bnPNFinal.fIdx);
                    //value = value / bnPNFinal.fW; 
                    histPN -> SetBinContent(bnPNFinal.fIdx,value);
                  }

                  if (iSys==0) histPN -> Draw("histpl");
                  else
                    histPN -> Draw("samehistpl");
                  legend -> AddEntry(histPN ,Form("%d %.2f",fSysBeams[iSys],histPN->Integral()),"pl");
                }
                makeLegend(cvsPN,legend,"lt") -> Draw();
                */

                for (auto iSys : selSysIdx)
                {
                  const char *namePN = makeName("psuedoN_ptoa",iAna,iLR,iMult,iSys,iCutTheta,iCutYP);
                  auto titlePN = makeTitle("",iAna,iLR,iMult,iSys,iCutTheta,iCutYP);
                  auto cvsPN = makeCvs(namePN);
                  auto hist = (bnRPtoa*bnPN).newHist(namePN,titlePN);
                  hist -> Draw();

                  auto legend = new TLegend();

                  auto histp = histPtoaArray[iSys][iCutTheta][iCutYP][kP];
                  auto histt = histPtoaArray[iSys][iCutTheta][iCutYP][kT];
                  auto hist3 = histPtoaArray[iSys][iCutTheta][iCutYP][kHe3];
                  setAtt(histp,kP);
                  setAtt(histt,kT);
                  setAtt(hist3,kHe3);

                  const char *nameHistPN = Form("%s_sys%d_pn",namePN,iSys);
                  auto histPN = (TH1D *) histt -> Clone(nameHistPN);
                  histPN -> Multiply(histp);
                  histPN -> Divide(hist3);
                  setAtt(histPN,5);

                  binning bnPNFinal = bnRPtoa;//(histPN);
                  while (bnPNFinal.next())
                  {
                    //auto value = histPN -> GetBinContent(bnPNFinal.fIdx);
                    //value = value / bnPNFinal.fW; 
                    //histPN -> SetBinContent(bnPNFinal.fIdx,value);
                  }

                  histPN -> Draw("samepl");
                  histp -> Draw("samepl");
                  histt -> Draw("samepl");
                  hist3 -> Draw("samepl");

                  legend -> AddEntry(histPN ,Form("n %d %.2f",fSysBeams[iSys],histPN->Integral()),"pl");
                  legend -> AddEntry(histp ,Form("p %d %.2f",fSysBeams[iSys],histp->Integral()),"pl");
                  legend -> AddEntry(histt ,Form("t %d %.2f",fSysBeams[iSys],histt->Integral()),"pl");
                  legend -> AddEntry(hist3 ,Form("he3 %d %.2f",fSysBeams[iSys],hist3->Integral()),"pl");

                  makeLegend(cvsPN,legend,"lt") -> Draw();
                }
              }

            }
          }
        }

        if (drawY0R21)
        {
          for (auto iCutYP : selCutYPIdx)
          {
            if (selCutYP>=0 && selCutYP!=iCutYP) continue;

            for (auto iCutTheta : selCutThetaIdx)
            {
              if (selCutTheta>=0 && selCutTheta!=iCutTheta) continue;

              for (auto iComb : selSysCombIdx)
              {
                if (selSysComb>=0 && selSysComb!=iComb) continue;

                auto iSys2 = fSysCombIdx[iComb][0];
                auto iSys1 = fSysCombIdx[iComb][1];
                auto iSysComb = 100 + iSys2*10 + iSys1;

                if (setDrawR21)
                {
                  const char *nameR21 = makeName("r21_y0All",iAna,iLR,iMult,iSysComb,iCutTheta,iCutYP);
                  auto titleR21 = makeTitle("Cor.",iAna,iLR,iMult,iSysComb,iCutTheta,iCutYP);

                  auto cvs = makeCvs(nameR21);
                  auto hist = (bnRY0*bnR21).newHist(nameR21,titleR21);
                  hist -> Draw();

                  auto legend = new TLegend();
                  for (auto iParticle : fParticleIdx)
                  {
                    const char *nameParticle = fParticleNames[iParticle];

                    auto hist1 = histY0Array[iSys2][iCutTheta][iCutYP][iParticle];
                    auto hist2 = histY0Array[iSys1][iCutTheta][iCutYP][iParticle];

                    const char *nameR21Part = Form("%s_%s",nameR21,nameParticle);

                    auto hist0 = (TH1D *) hist1 -> Clone(nameR21Part);
                    hist0 -> Divide(hist2);
                    auto graph = new TGraph();
                    for (auto bin=1; bin<=bnRY0.fN; ++bin) {
                      //if (hist1->GetBinContent(bin)<0.04||hist2->GetBinContent(bin)<0.04) hist0 -> SetBinContent(bin,0);
                      //else
                        if (hist1->GetBinContent(bin)==0) {}
                      else graph -> SetPoint(graph -> GetN(), hist0 -> GetXaxis() -> GetBinCenter(bin), hist0 -> GetBinContent(bin));
                    }

                    graph -> SetMarkerSize(1.3);
                    setAtt(graph,iParticle);

                    graph -> Draw("samepl");
                    legend -> AddEntry(graph,fParticleNames[iParticle],"pl");
                  }
                  makeLegend(cvs,legend) -> Draw();
                }

                if (setDrawSR || setDrawDR)
                {
                  for (auto iPR : fPRIdx)
                  {
                    auto iParticle1 = fPRParticleIdx[iPR][0];
                    auto iParticle2 = fPRParticleIdx[iPR][1];

                    const char *namePR = makeName(Form("ratioSR_%s",fPRNames[iPR]),iAna,iLR,iMult,iSysComb,iCutTheta,iCutYP);
                    auto titlePR = makeTitle("",iAna,iLR,iMult,iSysComb,iCutTheta,iCutYP);

                    TCanvas *cvsPR;
                    if (setDrawSR) {
                      cvsPR = makeCvs(namePR);
                      auto hist = (bnRY0*bnPR).newHist(namePR,titlePR);
                      hist -> Draw();
                    }

                    auto legend = new TLegend();
                    TH1D *histSys[fNumSyss] = {0};
                    for (auto iSys : {iSys2,iSys1})
                    {
                      auto hist_p1_s1 = histY0Array[iSys][iCutTheta][iCutYP][iParticle1];
                      auto hist_p2_s1 = histY0Array[iSys][iCutTheta][iCutYP][iParticle2];

                      const char *namePRSys = Form("%s_sys%d",namePR,iSys);
                      auto hist_s1 = (TH1D *) hist_p1_s1 -> Clone(namePRSys);
                      hist_s1 -> Divide(hist_p2_s1);
                      histSys[iSys] = hist_s1;

                      if (setDrawSR) {
                        setAtt(hist_s1,iSys);
                        hist_s1 -> Draw("samehistpl");
                        legend -> AddEntry(hist_s1,Form("%s (%d)",fPRTitles[iPR],fSysBeams[iSys]),"pl");
                      }
                    }
                    if (setDrawSR)
                      makeLegend(cvsPR,legend,"",0,0,0.4) -> Draw();

                    if (setDrawDR) {
                      const char *nameDR = makeName(Form("ratioDR_%s",fPRNames[iPR]),iAna,iLR,iMult,iSysComb,iCutTheta,iCutYP);
                      auto titleDR = makeTitle("",iAna,iLR,iMult,iSysComb,iCutTheta,iCutYP);

                      auto cvsDR = makeCvs(nameDR);//,500,510);
                      const char *ttlDR = Form("%s;y_{0};DR(%s/%s)[%d/%d]",titlePR,fParticleNames[iParticle1],fParticleNames[iParticle2],fSysBeams[iSys2],fSysBeams[iSys1]);
                      if (iPR==kDOP) {
                        auto hist = new TH2D(nameDR,ttlDR,nbinsFrame,bnRY0.fMin,bnRY0.fMax,bnDR.fN,1.1,1.5);
                        hist -> Draw();
                      }
                      else if (iPR==kTOP) {
                        auto hist = new TH2D(nameDR,ttlDR,nbinsFrame,bnRY0.fMin,bnRY0.fMax,bnDR.fN,1.3,2.1);
                        hist -> Draw();
                      }
                      else {
                        auto hist = (bnRY0*bnDR).newHist(nameDR,ttlDR);
                        hist -> Draw();
                      }

                      auto histDR = (TH1D *) histSys[iSys2] -> Clone(nameDR);
                      histDR -> Divide(histSys[iSys1]);
                      histDR -> SetMarkerStyle(24);
                      histDR -> SetMarkerColor(kBlack);
                      histDR -> SetLineColor(kBlack);
                      histDR -> Draw("samehistpl");

                      auto legend = new TLegend();
                      legend -> AddEntry(histDR ,Form("%s (%d/%d)",fPRTitles[iPR],fSysBeams[iSys2],fSysBeams[iSys1]),"pl");
                      makeLegend(cvsDR,legend,"",0,0,0.4) -> Draw();
                    }
                  }
                }

              }

              if (setDrawTemp)
              {
                const char *nameTemp = makeName("temp_y0",iAna,iLR,iMult,kSysAll,iCutTheta,iCutYP);
                auto titleTemp = makeTitle("",iAna,iLR,iMult,kSysAll,iCutTheta,iCutYP);
                auto cvsTemp = makeCvs(nameTemp);
                auto hist = (bnRY0*bnTemp).newHist(nameTemp,titleTemp);
                hist -> Draw();

                auto legend = new TLegend();
                for (auto iSys : selSysIdx)
                {
                  auto histd = histY0Array[iSys][iCutTheta][iCutYP][kD];
                  auto hist4 = histY0Array[iSys][iCutTheta][iCutYP][kHe4];
                  auto histt = histY0Array[iSys][iCutTheta][iCutYP][kT];
                  auto hist3 = histY0Array[iSys][iCutTheta][iCutYP][kHe3];

                  const char *nameHistTemp = Form("%s_sys%d_t",nameTemp,iSys);
                  auto histTemp = (TH1D *) histd -> Clone(nameHistTemp);
                  histTemp -> Multiply(hist4);
                  histTemp -> Divide(histt);
                  histTemp -> Divide(hist3);
                  setAtt(histTemp,iSys);

                  //binning bnTempFinal(histTemp);
                  binning bnTempFinal = bnTemp;//(histTemp);
                  while (bnTempFinal.next())
                  {
                    auto value = histTemp -> GetBinContent(bnTempFinal.fIdx);
                    value = 14.3/TMath::Log(1.59*value);
                    histTemp -> SetBinContent(bnTempFinal.fIdx,value);
                  }
                  histTemp -> Draw("samehistpl");
                  legend -> AddEntry(histTemp ,Form("%d",fSysBeams[iSys]),"pl");
                }
                makeLegend(cvsTemp,legend,"lt") -> Draw();
              }

              if (setDrawPN)
              {
                const char *namePN = makeName("psuedoN_y0",iAna,iLR,iMult,kSysAll,iCutTheta,iCutYP);
                auto titlePN = makeTitle("",iAna,iLR,iMult,kSysAll,iCutTheta,iCutYP);
                auto cvsPN = makeCvs(namePN);
                auto hist = (bnRY0*bnPN).newHist(namePN,titlePN);
                hist -> Draw();

                auto legend = new TLegend();
                for (auto iSys : selSysIdx)
                {
                  auto histp = histY0Array[iSys][iCutTheta][iCutYP][kP];
                  auto histt = histY0Array[iSys][iCutTheta][iCutYP][kT];
                  auto hist3 = histY0Array[iSys][iCutTheta][iCutYP][kHe3];

                  const char *nameHistPN = Form("%s_sys%d_pn",namePN,iSys);
                  auto histPN = (TH1D *) histt -> Clone(nameHistPN);
                  histPN -> Multiply(histp);
                  histPN -> Divide(hist3);
                  setAtt(histPN,iSys);

                  //binning bnPNFinal(histPN);
                  binning bnPNFinal = bnRY0;//(histPN);
                  while (bnPNFinal.next())
                  {
                    auto value = histPN -> GetBinContent(bnPNFinal.fIdx);
                    //value = value / bnPNFinal.fW; 
                    histPN -> SetBinContent(bnPNFinal.fIdx,value);
                  }

                  if (iSys==0) histPN -> Draw("histpl");
                  else
                    histPN -> Draw("samehistpl");
                  legend -> AddEntry(histPN ,Form("%d %.2f",fSysBeams[iSys],histPN->Integral()),"pl");
                }
                makeLegend(cvsPN,legend,"lt") -> Draw();
              }
            }
          }
        }


        if (drawKeoaR21)
        {
          for (auto iCutYP : selCutYPIdx)
          {
            if (selCutYP>=0 && selCutYP!=iCutYP) continue;

            for (auto iCutTheta : selCutThetaIdx)
            {
              if (selCutTheta>=0 && selCutTheta!=iCutTheta) continue;

              for (auto iComb : selSysCombIdx)
              {
                if (selSysComb>=0 && selSysComb!=iComb) continue;

                auto iSys2 = fSysCombIdx[iComb][0];
                auto iSys1 = fSysCombIdx[iComb][1];
                auto iSysComb = 100 + iSys2*10 + iSys1;

                const char *nameR21 = makeName("r21_keoa",iAna,iLR,iMult,iSysComb,iCutTheta,iCutYP);
                auto titleR21 = makeTitle("Cor.",iAna,iLR,iMult,iSysComb,iCutTheta,iCutYP);

                TCanvas *cvs;
                if (setDrawYield) {
                  cvs = makeCvs(nameR21,1000,500);
                  cvs -> Divide(2,1);
                }
                else
                  cvs = makeCvs(nameR21);

                auto hist = (bnRKeoaCM*bnR21).newHist(nameR21,titleR21);
                if (setDrawYield) cvs -> cd(1);
                hist -> Draw();

                auto legend = new TLegend();
                double ymax = 0;
                for (auto iParticle : fParticleIdx)
                {
                  auto hist2 = histKeoaArray[iLR][iSys2][iCutTheta][iCutYP][iParticle];
                  auto hist1 = histKeoaArray[iLR][iSys1][iCutTheta][iCutYP][iParticle];
                  if (ymax < hist2 -> GetMaximum()) ymax = hist2 -> GetMaximum(); 
                  if (ymax < hist1 -> GetMaximum()) ymax = hist1 -> GetMaximum(); 
                }
                for (auto iParticle : fParticleIdx)
                {
                  const char *nameParticle = fParticleNames[iParticle];

                  auto hist2 = histKeoaArray[iLR][iSys2][iCutTheta][iCutYP][iParticle];
                  auto hist1 = histKeoaArray[iLR][iSys1][iCutTheta][iCutYP][iParticle];

                  const char *nameR21Part = Form("%s_%s",nameR21,nameParticle);

                  auto hist0 = (TH1D *) hist2 -> Clone(nameR21Part);
                  hist0 -> Divide(hist1);

                  setAtt(hist2,iParticle);
                  setAtt(hist1,iParticle);
                  setAtt(hist0,iParticle);

                  if (setDrawYield) {
                    cvs -> cd(2);
                    if (iParticle==0) {
                      hist2 -> Draw("histpl");
                      hist2 -> SetMaximum(ymax*1.05);
                      hist2 -> SetMinimum(0);
                    }
                    else
                      hist2 -> Draw("histplsame");
                    hist1 -> SetLineStyle(2);
                    hist1 -> Draw("histplsame");

                    cvs -> cd(1);
                  }
                  //hist0 -> Draw("same histpl");
                  hist0 -> Draw("same histp");
                  legend -> AddEntry(hist0,fParticleNames[iParticle],"p");

                  if (histAPrev1[iComb]!=nullptr)
                  {
                    auto graph = new TGraph();
                    bnRKeoaCM.reset();
                    while (bnRKeoaCM.next())
                    {
                      auto idx = bnRKeoaCM.fIdx;
                      auto n0 = fParticleN[iParticle];
                      auto z0 = fParticleZ[iParticle];
                      auto c0 = histCPrev1[iComb] -> GetBinContent(idx,1);
                      auto alpha0 = histAPrev1[iComb] -> GetBinContent(idx,1);
                      auto beta0 = histBPrev1[iComb] -> GetBinContent(idx,1);
                      graph -> SetPoint(graph->GetN(),bnRKeoaCM.getValue(),c0*exp(n0*alpha0+z0*beta0));
                    }
                    setAtt(graph,iParticle);
                    graph -> SetLineStyle(2);
                    graph -> Draw("samel");
                  }
                }
                if (histAPrev1[iComb]!=nullptr) {
                  auto graph00 = new TGraph;
                  graph00 -> SetLineStyle(2);
                  legend -> AddEntry(graph00,"(fit)","l");
                }
                makeLegend(cvs,legend) -> Draw();

                if (histAPrev1[iComb]!=nullptr)
                {
                  const char *nameABComb = makeName("abcomb_keoa",iAna,iLR,iMult,iSysComb,iCutTheta,iCutYP);
                  auto titleABComb = makeTitle("",iAna,iLR,iMult,iSysComb,iCutTheta,iCutYP);
                  auto cvsABComb = makeCvs(nameABComb);
                  {
                  auto hist = (bnRKeoaCM*binning(100,-.5,1,"#alpha+#beta")).newHist(nameABComb,titleABComb);
                  hist -> Draw();
                  }

                  auto line0 = new TLine(bnRKeoaCM.lowEdge(),0,bnRKeoaCM.highEdge(),0);
                  auto line05 = new TLine(bnRKeoaCM.lowEdge(),0.5,bnRKeoaCM.highEdge(),0.5);
                  line0  -> SetLineColor(kGray);
                  line05 -> SetLineColor(kGray);
                  line0  -> SetLineStyle(2);
                  line05 -> SetLineStyle(2);
                  line0  -> Draw("samel");
                  line05 -> Draw("samel");

                  const char *nameABComb2 = makeName("abcomb2_keoa",iAna,iLR,iMult,iSysComb,iCutTheta,iCutYP);
                  auto titleABComb2 = makeTitle("",iAna,iLR,iMult,iSysComb,iCutTheta,iCutYP);
                  auto cvsABComb2 = makeCvs(nameABComb2);
                  {
                  auto hist = (bnRKeoaCM*binning(100,-.5,1,"#alpha-#beta")).newHist(nameABComb2,titleABComb2);
                  hist -> Draw();
                  }

                  line0 = new TLine(bnRKeoaCM.lowEdge(),0,bnRKeoaCM.highEdge(),0);
                  line05 = new TLine(bnRKeoaCM.lowEdge(),0.5,bnRKeoaCM.highEdge(),0.5);
                  line0  -> SetLineColor(kGray);
                  line05 -> SetLineColor(kGray);
                  line0  -> SetLineStyle(2);
                  line05 -> SetLineStyle(2);
                  line0  -> Draw("samel");
                  line05 -> Draw("samel");

                  binning bnX = bnRKeoaCM;

                  auto legendABComb2 = new TLegend();

                  {
                    auto graph_ake = new TGraph();
                    auto graph_bke = new TGraph();
                    auto graph_apbke = new TGraph();
                    auto graph_ambke = new TGraph();
                    bnX.reset();
                    while (bnX.next()) {
                      auto alpha0 = histAPrev1[iComb] -> GetBinContent(bnX.bi(),1);
                      auto beta0 = histBPrev1[iComb] -> GetBinContent(bnX.bi(),1);
                      graph_ake -> SetPoint(graph_ake->GetN(),bnX.getValue(),alpha0);
                      graph_bke -> SetPoint(graph_bke->GetN(),bnX.getValue(),beta0);
                      graph_apbke -> SetPoint(graph_apbke->GetN(),bnX.getValue(),alpha0+beta0);
                      graph_ambke -> SetPoint(graph_ambke->GetN(),bnX.getValue(),alpha0-beta0);
                    }

                    cvsABComb -> cd();
                    setAtt(graph_ake,1);   graph_ake   -> Draw("samelp");
                    setAtt(graph_bke,2);   graph_bke   -> Draw("samelp");
                    setAtt(graph_apbke,0,4); graph_apbke -> Draw("samelp");
                    auto legendABComb = new TLegend();
                    legendABComb -> AddEntry(graph_ake,"#alpha","p");
                    legendABComb -> AddEntry(graph_bke,"#beta","p");
                    legendABComb -> AddEntry(graph_apbke,"#alpha+#beta","p");
                    makeLegend(cvsABComb,legendABComb,"",0,0,0.4,0.3) -> Draw();

                    cvsABComb2 -> cd();
                    setAtt(graph_ake,1);   graph_ake   -> Draw("samelp");
                    setAtt(graph_bke,2);   graph_bke   -> Draw("samelp");
                    setAtt(graph_ambke,0,4); graph_ambke -> Draw("samelp");
                    auto legendABComb2 = new TLegend();
                    legendABComb2 -> AddEntry(graph_ake,"#alpha","p");
                    legendABComb2 -> AddEntry(graph_bke,"#beta","p");
                    legendABComb2 -> AddEntry(graph_ambke,"#alpha-#beta","p");
                    makeLegend(cvsABComb2,legendABComb2,"",0,0,0.4,0.3) -> Draw();
                  }
                }

              }

              if (setDrawTemp)
              {
                const char *nameTemp = makeName("temp_keoa",iAna,iLR,iMult,kSysAll,iCutTheta,iCutYP);
                auto titleTemp = makeTitle("",iAna,iLR,iMult,kSysAll,iCutTheta,iCutYP);
                auto cvsTemp = makeCvs(nameTemp);
                auto hist = (bnRKeoaCM*bnTemp).newHist(nameTemp,titleTemp);
                hist -> Draw();

                auto legend = new TLegend();
                for (auto iSys : selSysIdx)
                {
                  auto histd = histKeoaArray[iLR][iSys][iCutTheta][iCutYP][kD];
                  auto hist4 = histKeoaArray[iLR][iSys][iCutTheta][iCutYP][kHe4];
                  auto histt = histKeoaArray[iLR][iSys][iCutTheta][iCutYP][kT];
                  auto hist3 = histKeoaArray[iLR][iSys][iCutTheta][iCutYP][kHe3];

                  const char *nameHistTemp = Form("%s_sys%d_t",nameTemp,iSys);
                  auto histTemp = (TH1D *) histd -> Clone(nameHistTemp);
                  histTemp -> Multiply(hist4);
                  histTemp -> Divide(histt);
                  histTemp -> Divide(hist3);
                  setAtt(histTemp,iSys);

                  //binning bnTempFinal(histTemp);
                  binning bnTempFinal = bnTemp;//(histTemp);
                  while (bnTempFinal.next())
                  {
                    auto value = histTemp -> GetBinContent(bnTempFinal.fIdx);
                    value = 14.3/TMath::Log(1.59*value);
                    histTemp -> SetBinContent(bnTempFinal.fIdx,value);
                  }
                  histTemp -> Draw("samehistpl");
                  legend -> AddEntry(histTemp ,Form("%d",fSysBeams[iSys]),"pl");
                }
                makeLegend(cvsTemp,legend,"lt") -> Draw();

              }
            }
          }
        }

        if (drawTCMR21)
        {
          for (auto iCutYP : selCutYPIdx)
          {
            if (selCutYP>=0 && selCutYP!=iCutYP) continue;

            for (auto iCutTheta : selCutThetaIdx)
            {
              if (selCutTheta>=0 && selCutTheta!=iCutTheta) continue;

              for (auto iComb : selSysCombIdx)
              {
                if (selSysComb>=0 && selSysComb!=iComb) continue;

                auto iSys2 = fSysCombIdx[iComb][0];
                auto iSys1 = fSysCombIdx[iComb][1];
                auto iSysComb = 100 + iSys2*10 + iSys1;

                const char *nameR21 = makeName("r21_tcm",iAna,iLR,iMult,iSysComb,iCutTheta,iCutYP);
                auto titleR21 = makeTitle("Cor.",iAna,iLR,iMult,iSysComb,iCutTheta,iCutYP);

                auto cvs = makeCvs(nameR21,680,550);

                auto hist = (bnRTCM*bnR21).newHist(nameR21,titleR21);
                hist -> Draw();

                auto legend = new TLegend();
                for (auto iParticle : fParticleIdx)
                {
                  const char *nameParticle = fParticleNames[iParticle];

                  auto hist1 = histTCMArray[iLR][iSys2][iCutTheta][iCutYP][iParticle];
                  auto hist2 = histTCMArray[iLR][iSys1][iCutTheta][iCutYP][iParticle];

                  const char *nameR21Part = Form("%s_%s",nameR21,nameParticle);

                  auto hist0 = (TH1D *) hist1 -> Clone(nameR21Part);
                  hist0 -> Divide(hist2);

                  setAtt(hist0,iParticle);

                  hist0 -> Draw("same histpl");
                  legend -> AddEntry(hist0,fParticleNames[iParticle],"l");
                }
                makeLegend(cvs,legend) -> Draw();
              }
            }
          }
        }

        if (drawNZR21)
        {
          auto bnX = bnRY0;
          auto bnY = bnRPtoa;
          if (setDrawNZInKT) {
            bnX = bnRKeoaCM;
            bnY = bnRTCM;
          }

          for (auto iCutYP : selCutYPIdx)
          {
            if (selCutYP>=0 && selCutYP!=iCutYP) continue;

            for (auto iCutTheta : selCutThetaIdx)
            {
              if (selCutTheta>=0 && selCutTheta!=iCutTheta) continue;

              for (auto iComb : selSysCombIdx)
              {
                if (selSysComb>=0 && selSysComb!=iComb) continue;

                auto iSys2 = fSysCombIdx[iComb][0];
                auto iSys1 = fSysCombIdx[iComb][1];
                auto iSysComb = 100 + iSys2*10 + iSys1;

                const char *nameAlpha = makeName("alpha",iAna,iLR,iMult,iSysComb,iCutTheta,iCutYP);
                const char *nameMBeta = makeName("mbeta",iAna,iLR,iMult,iSysComb,iCutTheta,iCutYP);
                if (setDrawNZAll) {
                  nameAlpha = makeName("nzall_alpha",iAna,iLR,iMult,iSysComb,iCutTheta,iCutYP);
                  nameMBeta = makeName("nzall_mbeta",iAna,iLR,iMult,iSysComb,iCutTheta,iCutYP);
                }

                TCanvas *cvsAlpha;
                TCanvas *cvsMBeta;
                if (setDrawNZAll) {
                  cvsAlpha = makeCvs(nameAlpha);
                  cvsMBeta = makeCvs(nameMBeta);
                }
                else {
                  cvsAlpha = makeCvs(nameAlpha,100+dxCvsAB*bnX.fN,dyCvsAB*bnY.fN);
                  cvsMBeta = makeCvs(nameMBeta,100+dxCvsAB*bnX.fN,dyCvsAB*bnY.fN);
                  cvsDivideM0(cvsAlpha, bnX.fN, bnY.fN);
                  cvsDivideM0(cvsMBeta, bnX.fN, bnY.fN);
                }

                double temperatureArray[10][10][3] = {{0}};
                double r21Array[10][10][fNumParticles] = {{0}};
                double r21ErrorArray[10][10][fNumParticles] = {{0}};
                double alphaArray[10][10] = {{0}};
                double mbetaArray[10][10] = {{0}};
                double muNArray[10][10] = {{0}};
                double muZArray[10][10] = {{0}};
                double cnormArray[10][10] = {{0}};
                double alphaErrorArray[10][10] = {{0}};
                double mbetaErrorArray[10][10] = {{0}};

                TH2D *histYPR21[2][fNumParticles] = {{0}}; // 2/1, 1, 2
                TH2D *histYield[2][fNumParticles] = {{0}}; // 2/1, 1, 2

                /*
                const char *nameSys2 = makeName("sys2",iAna,iLR,iMult,iSysComb,iCutTheta,iCutYP);
                const char *nameSys1 = makeName("sys1",iAna,iLR,iMult,iSysComb,iCutTheta,iCutYP);
                auto cvsSys2 = makeCvs(nameSys2);
                auto cvsSys1 = makeCvs(nameSys1);
                */

                for (auto iParticle : fParticleIdx) {
                  const char *hname = makeName("hnzr21",iAna,iLR,iMult,iSysComb,iCutTheta,iCutYP,iParticle);
                  auto hist2 = (TH2D *) histXYR21Array[iSys2][iCutTheta][iCutYP][iParticle];
                  auto hist1 = (TH2D *) histXYR21Array[iSys1][iCutTheta][iCutYP][iParticle];
                  auto hist21 = (TH2D *) hist2 -> Clone(hname);
                  hist21 -> Divide(hist1);
                  //cvsSys2 -> cd(); hist2 -> Draw("colz");
                  //cvsSys1 -> cd(); hist1 -> Draw("colz");

                  const char *ename = makeName("enzr21",iAna,iLR,iMult,iSysComb,iCutTheta,iCutYP,iParticle);
                  auto eist2 = (TH2D *) histXYR21Array2[iSys2][iCutTheta][iCutYP][iParticle];
                  auto eist1 = (TH2D *) histXYR21Array2[iSys1][iCutTheta][iCutYP][iParticle];
                  auto eist21 = (TH2D *) eist2 -> Clone(ename);
                  eist21 -> Divide(eist1);

                  histYPR21[0][iParticle] = hist21;
                  histYPR21[1][iParticle] = eist21;

                  histYield[0][iParticle] = hist2;
                  histYield[1][iParticle] = hist1;
                }

                for (auto ibinY0=0; ibinY0<bnX.fN; ++ibinY0)
                {
                  for (auto ibinPtoa=0; ibinPtoa<bnY.fN; ++ibinPtoa)
                  {
                    auto binY0 = ibinY0 + 1;
                    auto binPtoa = ibinPtoa + 1;

                    bool setDrawLegend = false;
                    //bool setDrawLegend = ((!setDrawNZAll) && ibinY0==bnX.fN-1&&ibinPtoa==bnY.fN-1);

                    auto graphAll = new TGraph2DErrors();
                    auto graph_z1_inN = new TGraphErrors(); graph_z1_inN -> SetMarkerStyle(10);
                    auto graph_z2_inN = new TGraphErrors(); graph_z2_inN -> SetMarkerStyle(10);
                    auto graph_n0_inZ = new TGraphErrors(); graph_n0_inZ -> SetMarkerStyle(10);
                    auto graph_n1_inZ = new TGraphErrors(); graph_n1_inZ -> SetMarkerStyle(10);
                    auto graph_n2_inZ = new TGraphErrors(); graph_n2_inZ -> SetMarkerStyle(10);
                    vector<TMarker *> markersAlpha;
                    vector<TMarker *> markersBeta;

                    for (auto iParticle : fParticleIdx)
                    {
                      auto r21Value = histYPR21[0][iParticle] -> GetBinContent(binY0, binPtoa);
                      auto r21Value_2 = histYPR21[1][iParticle] -> GetBinContent(binY0, binPtoa);
                      auto r21Error = abs(r21Value - r21Value_2);
                      r21Array[binY0][binPtoa][iParticle] = r21Value;
                      r21ErrorArray[binY0][binPtoa][iParticle] = r21Error;

                      //auto value2 = histXYR21Array[iSys2][iCutTheta][iCutYP][iParticle] -> GetBinContent(binY0, binPtoa);
                      //auto value1 = histXYR21Array[iSys1][iCutTheta][iCutYP][iParticle] -> GetBinContent(binY0, binPtoa);
                      //auto r21Value = value2/value1;
                      //r21Array[binY0][binPtoa][iParticle] = r21Value;
                      //auto r21Error = sqrt((1/sqrt(value1))*(1/sqrt(value1)) * (1/sqrt(value2))*(1/sqrt(value2)));
                      //auto r21Error = 0;//0.001*sqrt((1/sqrt(value1))*(1/sqrt(value1)) * (1/sqrt(value2))*(1/sqrt(value2)));
                      //r21ErrorArray[binY0][binPtoa][iParticle] = r21Error;

                      auto zValue = fParticleZ[iParticle];
                      auto nValue = fParticleN[iParticle];

                      if (zValue==1) { graph_z1_inN -> SetPoint(graph_z1_inN->GetN(), nValue, r21Value); }
                      if (zValue==2) { graph_z2_inN -> SetPoint(graph_z2_inN->GetN(), nValue, r21Value); }
                      if (nValue==0) { graph_n0_inZ -> SetPoint(graph_n0_inZ->GetN(), zValue, r21Value); }
                      if (nValue==1) { graph_n1_inZ -> SetPoint(graph_n1_inZ->GetN(), zValue, r21Value); }
                      if (nValue==2) { graph_n2_inZ -> SetPoint(graph_n2_inZ->GetN(), zValue, r21Value); }

                      if (zValue==1) { graph_z1_inN -> SetPointError(graph_z1_inN->GetN()-1, 0, r21Error); }
                      if (zValue==2) { graph_z2_inN -> SetPointError(graph_z2_inN->GetN()-1, 0, r21Error); }
                      if (nValue==0) { graph_n0_inZ -> SetPointError(graph_n0_inZ->GetN()-1, 0, r21Error); }
                      if (nValue==1) { graph_n1_inZ -> SetPointError(graph_n1_inZ->GetN()-1, 0, r21Error); }
                      if (nValue==2) { graph_n2_inZ -> SetPointError(graph_n2_inZ->GetN()-1, 0, r21Error); }

                      auto mAlpha = new TMarker(zValue, r21Value, 20); setAtt(mAlpha,iParticle);
                      auto mBeta = new TMarker(nValue, r21Value, 20); setAtt(mBeta,iParticle);
                      markersAlpha.push_back(mAlpha);
                      markersBeta.push_back(mBeta);

                      if (setIgnoreT&&iParticle!=2)
                        graphAll -> SetPoint(graphAll->GetN(), nValue, zValue, r21Value);
                      else {
                        graphAll -> SetPoint(graphAll->GetN(), nValue, zValue, r21Value);
                        graphAll -> SetPointError(graphAll->GetN()-1, 0, 0, r21Value);
                      }
                    }

                    for (auto jSys : {0,1})
                    {
                      auto yieldd = histYield[jSys][kD]   -> GetBinContent(binY0, binPtoa);
                      auto yield4 = histYield[jSys][kHe4] -> GetBinContent(binY0, binPtoa);
                      auto yieldt = histYield[jSys][kT]   -> GetBinContent(binY0, binPtoa);
                      auto yield3 = histYield[jSys][kHe3] -> GetBinContent(binY0, binPtoa);
                      temperatureArray[binY0][binPtoa][jSys] = 14.3/TMath::Log(1.59*(yieldd * yield4 / yieldt / yield3));
                    }
                    temperatureArray[binY0][binPtoa][2] = (temperatureArray[binY0][binPtoa][0] + temperatureArray[binY0][binPtoa][1]) / 2.;

                    auto idxCvs = (bnY.fN-ibinPtoa-1)*bnX.fN+ibinY0+1;

                    auto cvsAlpha0 = cvsAlpha -> cd();
                    if (!setDrawNZAll)
                      cvsAlpha0 = cvsAlpha -> cd(idxCvs);
                    if (setDrawNZLogy) cvsAlpha0 -> SetLogy();
                    //auto frameAlpha = new TH2D(Form("%s_%d",cvsAlpha->GetName(),idxCvs),";N;R21",100,-1,3,100,bnR21.fMin,bnR21.fMax);
                    auto frameAlpha = (bnParticleN*bnR21).newHist(Form("%s_%d",cvsAlpha->GetName(),idxCvs));
                    if (!setDrawNZAll) {
                      frameAlpha -> GetXaxis() -> SetTitleSize(0.10);
                      frameAlpha -> GetYaxis() -> SetTitleSize(0.10);
                      frameAlpha -> GetXaxis() -> SetLabelSize(0.10);
                      frameAlpha -> GetYaxis() -> SetLabelSize(0.10);
                      frameAlpha -> GetXaxis() -> SetNdivisions(506);
                      frameAlpha -> GetYaxis() -> SetNdivisions(506);
                      frameAlpha -> GetXaxis() -> CenterTitle();
                      frameAlpha -> GetYaxis() -> CenterTitle();
                    }
                    frameAlpha -> GetXaxis() -> SetNdivisions(4,0,0,true);///*n=*/fAxisNDivisions/*,optim=true*/);
                    frameAlpha -> Draw();
                    TLegend *legendAlpha = new TLegend();

                    double alphaValue = 0;
                    double betaValue = 0;
                    double alphaError = 0;
                    double betaError = 0;
                    double cnormValue = 0;
                    double muNValue = 0;
                    double muZValue = 0;

                    TF2 *fitIsocaling = new TF2("fitIsocaling","[2]*exp([0]*x+[1]*y)",1,2,0,2);
                    fitIsocaling -> SetParameters(1,.3,-.3);
                    if (setFitNZR21TG) {
                      auto fitResultPtr = graphAll -> Fit(fitIsocaling,"","RQ0");

                      alphaValue = fitIsocaling -> GetParameter(0);
                      betaValue = fitIsocaling -> GetParameter(1);
                      cnormValue = fitIsocaling -> GetParameter(2);

                      double tempMean = (temperatureArray[binY0][binPtoa][0] + temperatureArray[binY0][binPtoa][1])/2.;
                      muNValue = alphaValue*tempMean;
                      muZValue = betaValue*tempMean;

                      alphaArray[binY0][binPtoa] = alphaValue;
                      mbetaArray[binY0][binPtoa] = -betaValue;
                      cnormArray[binY0][binPtoa] = cnormValue;
                      //muNArray[binY0][binPtoa] = muNValue;
                      //muZArray[binY0][binPtoa] = muZValue;

                      if (setABErrorToResidual)
                      {
                        double aberror = 0.;
                        /*
                        for (auto iParticle : fParticleIdx)
                        {
                          if (setIgnoreT&&iParticle==2) continue;
                          auto nValue = fParticleN[iParticle];
                          auto zValue = fParticleZ[iParticle];
                          auto r21Value = r21Array[binY0][binPtoa][iParticle];
                          auto residual = r21Value - fitIsocaling -> Eval(nValue, zValue);
                          aberror += (residual/r21Value)*(residual/r21Value);
                        }
                        if (setIgnoreT)
                          aberror = sqrt(aberror/4);
                        else
                          aberror = sqrt(aberror/5);

                        alphaErrorArray[binY0][binPtoa] = aberror;
                        mbetaErrorArray[binY0][binPtoa] = aberror;
                          */

                        double aError = 0.;
                        double bError = 0.;
                        for (auto iParticle : fParticleIdx)
                        {
                          auto nValue = fParticleN[iParticle];
                          auto zValue = fParticleZ[iParticle];
                          auto r21Value = r21Array[binY0][binPtoa][iParticle];
                          auto r21Error = r21ErrorArray[binY0][binPtoa][iParticle];

                          auto r21Eval = fitIsocaling -> Eval(nValue, zValue);

                          //auto aError1 = (r21Value - r21Eval) / r21Value;
                          //auto bError1 = (r21Value - r21Eval) / r21Value;
                          //auto aError1 = (r21Value - r21Eval) / r21Error / r21Value / nValue;
                          //auto bError1 = (r21Value - r21Eval) / r21Error / r21Value / zValue;
                          //aError += aError1*aError1;
                          //bError += bError1*bError1;

                          auto aError1 = (r21Value - r21Eval) / r21Value / nValue;
                          auto bError1 = (r21Value - r21Eval) / r21Value / zValue;
                          if (nValue!=0) aError += aError1*aError1;
                          if (zValue!=0) bError += bError1*bError1;
                        }
                        aError = sqrt(aError / 4);
                        bError = sqrt(bError / 5);

                        alphaErrorArray[binY0][binPtoa] = aError;
                        mbetaErrorArray[binY0][binPtoa] = bError;

                        //XXX
                        //alphaErrorArray[binY0][binPtoa] = 0;
                        //mbetaErrorArray[binY0][binPtoa] = 0;
                      }
                      else {
                        alphaErrorArray[binY0][binPtoa] = fitIsocaling -> GetParError(0);
                        mbetaErrorArray[binY0][binPtoa] = fitIsocaling -> GetParError(1);
                        /*
                        double rmsR21 = 0;
                        for (auto iParticle : fParticleIdx)
                          rmsR21 =+ r21ErrorArray[binY0][binPtoa][iParticle]*r21ErrorArray[binY0][binPtoa][iParticle];
                        rmsR21 = sqrt(rmsR21);
                        alphaErrorArray[binY0][binPtoa] = rmsR21;
                        mbetaErrorArray[binY0][binPtoa] = rmsR21;
                        */
                      }

                      alphaError = alphaErrorArray[binY0][binPtoa];
                      betaError = mbetaErrorArray[binY0][binPtoa];

                      TF1 *fit1 = new TF1("fit1","[2]*exp([0]*x+[1])",-.5,2.5);
                      fit1 -> SetParameters(alphaValue,betaValue*1,cnormValue);
                      fit1 -> SetLineColor(kGray+1);
                      fit1 -> SetLineStyle(2);

                      TF1 *fit2 = new TF1("fit2","[2]*exp([0]*x+[1])",0.5,2.5);
                      //TF1 *fit2 = new TF1("fit2","[2]*exp([0]*x+[1])",-.5,2.5);
                      fit2 -> SetParameters(alphaValue,betaValue*2,cnormValue);
                      fit2 -> SetLineColor(kGray+1);
                      //fit2 -> SetLineStyle(2);

                      if (!setDrawLegend) {
                        fit1 -> Draw("samel");
                        fit2 -> Draw("samel");
                      }
                      if (!setDrawLegend) {
                        graph_z1_inN -> Draw("same e");
                        graph_z2_inN -> Draw("same e");
                        for (auto m0 : markersBeta)
                          m0 -> Draw("p");
                      }

                      auto mm = new TMarker(0,0,fDrawMStyle2[ibinPtoa]);
                      //mm -> SetMarkerColor(fDrawColor[ibinY0]);
                      setAtt(mm,ibinPtoa,ibinY0,2);
                      legendAlpha -> AddEntry(mm, Form("#alpha = %.3f",alphaValue),"p");
                    }

                    if (setDrawLegend)
                    {
                      auto lengendIsotope = new TLegend();
                      int iParticle = 0;
                      for (auto m0 : markersAlpha) {
                        lengendIsotope -> AddEntry(m0,fParticleTitles[iParticle],"p");
                        iParticle++;
                      }
                      makeLegend(cvsAlpha0,lengendIsotope,"lt",.05,0,.98,.98) -> Draw();
                    }
                    else {
                      //makeLegend(cvsAlpha0,legendAlpha,"lt",0.2,0,.8,.20) -> Draw();
                      auto ma1 = new TMarker(-.4,bnR21.fMax-0.25*bnR21.getFullWidth(),fDrawMStyle2[ibinPtoa]);
                      setAtt(ma1,ibinPtoa,ibinY0,2);
                      ma1 -> Draw();
                      //auto aa1 = new TLatex(-.2,bnR21.fMax-0.15*bnR21.getFullWidth(),Form("#alpha = %.2f #pm %.3f",alphaValue,alphaError));
                      auto aa1 = new TLatex(-.1,bnR21.fMax-0.25*bnR21.getFullWidth(),Form("#alpha = %.2f",alphaValue));
                      aa1 -> SetTextAlign(12);
                      aa1 -> SetTextSize(0.1);
                      aa1 -> SetTextFont(42);
                      //aa1 -> SetTextColor(kRed+2);
                      aa1 -> Draw();
                      auto tt1 = new TLatex(-.6,bnR21.fMin+0.1*bnR21.getFullWidth(),Form("T_{%d} = %.1f",fSysCombBeam[iComb][0],temperatureArray[binY0][binPtoa][0]));
                      tt1 -> SetTextAlign(11);
                      tt1 -> SetTextSize(0.1);
                      tt1 -> SetTextFont(42);
                      tt1 -> Draw();
                      int dTOT2PC = floor(0.5+100*abs(temperatureArray[binY0][binPtoa][0] - temperatureArray[binY0][binPtoa][1]) / temperatureArray[binY0][binPtoa][0]);
                      auto tt2 = new TLatex(-.6,bnR21.fMin+0.03*bnR21.getFullWidth(),Form("T_{%d} = %.1f (%d%s)",fSysCombBeam[iComb][1],temperatureArray[binY0][binPtoa][1],dTOT2PC,"%"));
                      tt2 -> SetTextAlign(11);
                      tt2 -> SetTextSize(0.1);
                      tt2 -> SetTextFont(42);
                      tt2 -> Draw();

                      auto lineR211 = new TLine(bnParticleN.fMin,1,bnParticleN.fMax,1);
                      lineR211 -> SetLineColor(kGray);
                      lineR211 -> SetLineStyle(2);
                      lineR211 -> Draw("samel");
                    }

                    auto cvsMBeta0 = cvsMBeta -> cd();
                    if (!setDrawNZAll)
                      cvsMBeta0 = cvsMBeta -> cd(idxCvs);
                    if (setDrawNZLogy) cvsMBeta0 -> SetLogy();
                    //auto frameBeta = new TH2D(Form("%s_%d",cvsMBeta->GetName(),idxCvs),";Z;R21",100,0,3,100,bnR21.fMin,bnR21.fMax);
                    auto frameBeta = (bnParticleZ*bnR21).newHist(Form("%s_%d",cvsMBeta->GetName(),idxCvs));
                    if (!setDrawNZAll) {
                      frameBeta -> GetXaxis() -> SetTitleSize(0.10);
                      frameBeta -> GetYaxis() -> SetTitleSize(0.10);
                      frameBeta -> GetXaxis() -> SetLabelSize(0.10);
                      frameBeta -> GetYaxis() -> SetLabelSize(0.10);
                      frameBeta -> GetXaxis() -> SetNdivisions(506);
                      frameBeta -> GetYaxis() -> SetNdivisions(506);
                      frameBeta -> GetXaxis() -> CenterTitle();
                      frameBeta -> GetYaxis() -> CenterTitle();
                    }
                    frameBeta -> GetXaxis() -> SetNdivisions(4,0,0,true);///*n=*/fAxisNDivisions/*,optim=true*/);
                    frameBeta -> Draw();
                    TLegend *legendBeta = new TLegend();

                    if (setFitNZR21TG) {
                      //TF1 *fit0 = new TF1("fit0","[2]*exp([0]+[1]*x)",.5,2.5);
                      TF1 *fit0 = new TF1("fit0","[2]*exp([0]+[1]*x)",.5,1.5);
                      fit0 -> SetParameters(alphaValue*0,betaValue,cnormValue);
                      fit0 -> SetLineColor(kGray+1);

                      TF1 *fit1 = new TF1("fit1","[2]*exp([0]+[1]*x)",.5,2.5);
                      fit1 -> SetParameters(alphaValue*1,betaValue,cnormValue);
                      fit1 -> SetLineColor(kGray+1);
                      fit1 -> SetLineStyle(2);

                      TF1 *fit2 = new TF1("fit2","[2]*exp([0]+[1]*x)",.5,2.5);
                      fit2 -> SetParameters(alphaValue*2,betaValue,cnormValue);
                      fit2 -> SetLineColor(kGray+1);
                      //fit2 -> SetLineStyle(10);

                      if (!setDrawLegend) {
                        fit0 -> Draw("samel");
                        fit1 -> Draw("samel");
                        fit2 -> Draw("samel");
                      }
                      if (!setDrawLegend) {
                        graph_n0_inZ -> Draw("same e");
                        graph_n1_inZ -> Draw("same e");
                        graph_n2_inZ -> Draw("same e");
                        for (auto m0 : markersAlpha)
                          m0 -> Draw("p");
                      }

                      auto mm = new TMarker(0,0,fDrawMStyle2[ibinPtoa]);
                      //mm -> SetMarkerColor(fDrawColor[ibinY0]);
                      //setAtt(mm,ibinY0,-1,2);
                      setAtt(mm,ibinPtoa,ibinY0,2);
                      //mm -> SetMarkerSize(1.8);
                      legendBeta -> AddEntry(mm, Form("-#beta = %.3f",-betaValue),"p");
                    }
                    if (setDrawLegend)
                    {
                      auto lengendIsotope = new TLegend();
                      int iParticle = 0;
                      for (auto m0 : markersBeta) {
                        lengendIsotope -> AddEntry(m0,fParticleTitles[iParticle],"p");
                        iParticle++;
                      }
                      makeLegend(cvsMBeta0,lengendIsotope,"lt",.05,0,.98,.98) -> Draw();
                    }
                    else {
                      //makeLegend(cvsMBeta0,legendBeta,"lt",0.2,0,.8,.16) -> Draw();
                      auto mb1 = new TMarker(.5,bnR21.fMax-0.25*bnR21.getFullWidth(),fDrawMStyle2[ibinPtoa]);
                      setAtt(mb1,ibinPtoa,ibinY0,2);
                      mb1 -> Draw();
                      //auto aa1 = new TLatex(-.2,bnR21.fMax-0.15*bnR21.getFullWidth(),Form("#alpha = %.2f #pm %.3f",alphaValue,alphaError));
                      auto bb1 = new TLatex(.7,bnR21.fMax-0.25*bnR21.getFullWidth(),Form("#beta = %.2f",betaValue));
                      //auto bb1 = new TLatex(.2,bnR21.fMax-0.15*bnR21.getFullWidth(),Form("#beta = %.2f #pm %.3f",alphaValue,alphaError));
                      bb1 -> SetTextAlign(12);
                      bb1 -> SetTextSize(0.1);
                      bb1 -> SetTextFont(42);
                      //bb1 -> SetTextColor(kBlue+2);
                      bb1 -> Draw();
                      auto tt1 = new TLatex(0.4,bnR21.fMin+0.1*bnR21.getFullWidth(),Form("T_{%d} = %.1f",fSysCombBeam[iComb][0],temperatureArray[binY0][binPtoa][0]));
                      tt1 -> SetTextAlign(11);
                      tt1 -> SetTextSize(0.1);
                      tt1 -> SetTextFont(42);
                      tt1 -> Draw();
                      int dTOT2PC = floor(0.5+100*abs(temperatureArray[binY0][binPtoa][0] - temperatureArray[binY0][binPtoa][1]) / temperatureArray[binY0][binPtoa][0]);
                      auto tt2 = new TLatex(0.4,bnR21.fMin+0.03*bnR21.getFullWidth(),Form("T_{%d} = %.1f (%d%s)",fSysCombBeam[iComb][1],temperatureArray[binY0][binPtoa][1],dTOT2PC,"%"));
                      tt2 -> SetTextAlign(11);
                      tt2 -> SetTextSize(0.1);
                      tt2 -> SetTextFont(42);
                      tt2 -> Draw();

                      auto lineR211 = new TLine(bnParticleZ.fMin,1,bnParticleZ.fMax,1);
                      lineR211 -> SetLineColor(kGray);
                      lineR211 -> SetLineStyle(2);
                      lineR211 -> Draw("samel");
                    }
                  }
                }

                auto nameAlphaYP = makeName("alaph_yp",iAna,iLR,iMult,iSysComb,iCutTheta,iCutYP);
                auto histAlpha = (bnX*bnY).newHist(nameAlphaYP,Form("#alpha %s",fSysCombTitles[iComb]));

                auto nameMBetaYP = makeName("mbeta_yp",iAna,iLR,iMult,iSysComb,iCutTheta,iCutYP);
                auto histMBeta = (bnX*bnY).newHist(nameMBetaYP,Form("-#beta %s",fSysCombTitles[iComb]));

                auto nameAlphaTYP = makeName("alphaT_yp",iAna,iLR,iMult,iSysComb,iCutTheta,iCutYP);
                auto histAlphaT = (bnX*bnY).newHist(nameAlphaTYP,Form("#alpha#timesT %s",fSysCombTitles[iComb]));

                auto nameMBetaTYP = makeName("mbetaT_yp",iAna,iLR,iMult,iSysComb,iCutTheta,iCutYP);
                auto histMBetaT = (bnX*bnY).newHist(nameMBetaTYP,Form("-#beta#timesT %s",fSysCombTitles[iComb]));

                auto nameCNormYP = makeName("cnorm_yp",iAna,iLR,iMult,iSysComb,iCutTheta,iCutYP);
                auto histCNorm = (bnX*bnY).newHist(nameCNormYP,Form("c %s",fSysCombTitles[iComb]));

                auto nameErrorYP = makeName("error_yp",iAna,iLR,iMult,iSysComb,iCutTheta,iCutYP);
                auto histdT = (bnX*bnY).newHist(nameErrorYP,Form("100#times(T_{2}-T_{1})/T_{2} %s",fSysCombTitles[iComb]));

                auto nameKYP = makeName("k_yp",iAna,iLR,iMult,iSysComb,iCutTheta,iCutYP);
                auto histK = (bnX*bnY).newHist(nameKYP,Form("#alpha+#beta %s",fSysCombTitles[iComb]));

                histAlpha -> GetZaxis() -> SetRangeUser(bnAlpha.fMin,bnAlpha.fMax);
                //histMBeta -> GetZaxis() -> SetRangeUser(-bnBeta.fMax,-bnBeta.fMin);
                histMBeta -> GetZaxis() -> SetRangeUser(0.05,-bnBeta.fMin);
                //histMBeta -> GetZaxis() -> SetRangeUser(-0.05,0.05);
                histCNorm -> GetZaxis() -> SetRangeUser(0.9,1.2);
                //histdT -> GetZaxis() -> SetRangeUser(0,0.2);
                //histK -> GetZaxis() -> SetRangeUser(-0.22,0.1);


                /////////////////////////////////////////////////////////////////////////////////////////////
                /////////////////////////////////////////////////////////////////////////////////////////////
                /////////////////////////////////////////////////////////////////////////////////////////////
                /////////////////////////////////////////////////////////////////////////////////////////////
                /////////////////////////////////////////////////////////////////////////////////////////////
                /////////////////////////////////////////////////////////////////////////////////////////////

                for (bool iMM : {false,true})
                {
                  auto nameR21AB = makeName("r21_ab",iAna,iLR,iMult,iSysComb,iCutTheta,iCutYP);
                  binning bnA = bnAlpha;
                  binning bnB = bnBeta;

                  if (iMM)
                  {
                    nameR21AB = makeName("r21_mm",iAna,iLR,iMult,iSysComb,iCutTheta,iCutYP);
                    if (iComb!=3) {
                      bnA.set(100,-1,5,"#Delta#mu_{n} = #alpha #times <T>_{KE}");
                      bnB.set(100,-5,1,"#Delta#mu_{p} = #beta #times <T>_{KE}");
                    }
                    else {
                      bnA.set(100,-2,2,"#Delta#mu_{n} = #alpha #times <T>");
                      bnB.set(100,-2,2,"#Delta#mu_{p} = #beta #times <T>");
                    }
                  }

                  auto cvsABMM = makeCvs(nameR21AB,630,640);
                  auto titleR21AB = makeTitle("",iAna,iLR,iMult,iSysComb,iCutTheta,iCutYP);
                  auto frameR21AB = (bnA*bnB).newHist(nameR21AB,titleR21AB);
                  frameR21AB -> GetYaxis() -> SetTitleOffset(1.4);
                  frameR21AB -> Draw();

                  bool drawMarkers = true;

                  auto legendAB = new TLegend();
                  vector <TMarker*> markers;
                  TGraphErrors *graphABError[10] = {0};
                  TGraph       *graphABNoError[10] = {0};

                  bnY.end();
                  while (bnY.prev())
                  {
                    bnX.end();
                    while (bnX.prev())
                    {
                      if (graphABError[bnX.bi()]==nullptr) {
                        graphABError[bnX.bi()] = new TGraphErrors();
                        graphABNoError[bnX.bi()] = new TGraph();
                      }
                      auto alpha = alphaArray[bnX.bi()][bnY.bi()];
                      auto mbeta = mbetaArray[bnX.bi()][bnY.bi()];
                      auto beta = -mbetaArray[bnX.bi()][bnY.bi()];
                      auto cnorm = cnormArray[bnX.bi()][bnY.bi()];
                      auto alphaError = alphaErrorArray[bnX.bi()][bnY.bi()];
                      auto mbetaError = mbetaErrorArray[bnX.bi()][bnY.bi()];

                      double tempMean = 0;// = (temperatureArray[bnX.bi()][bnY.bi()][0] + temperatureArray[bnX.bi()][bnY.bi()][1])/2.;
                      if (histTPrev1[iComb]==nullptr) {
                        cout_error << "run once more" << endl;
                        tempMean = (temperatureArray[bnX.bi()][bnY.bi()][0] + temperatureArray[bnX.bi()][bnY.bi()][1])/2.;
                      }
                      else
                        tempMean = histTPrev1[iComb] -> GetBinContent(bnX.bi(),1);

                      auto muNValue = alpha * tempMean;
                      auto muZValue = beta * tempMean;

                      if (iMM)
                      {
                        alpha = muNValue;
                        beta = muZValue;
                        mbeta = -muZValue;
                        alphaError = alphaErrorArray[bnX.bi()][bnY.bi()] * tempMean;
                        mbetaError = mbetaErrorArray[bnX.bi()][bnY.bi()] * tempMean;
                      }
                      auto betaError = mbetaError;
                      TMarker *markerAB = nullptr;

                      if (bnA.isInside(alpha) && bnB.isInside(beta)) {
                        markerAB = new TMarker(alpha,beta,fDrawMStyle[bnY.ii()]);
                        setAtt(markerAB,bnY.ii(),bnX.ii(),2);
                        markers.push_back(markerAB);
                      }

                      graphABError[bnX.bi()] -> SetPoint(graphABError[bnX.bi()]->GetN(),alpha,beta);
                      graphABError[bnX.bi()] -> SetPointError(graphABError[bnX.bi()]->GetN()-1,alphaError,betaError);
                      graphABNoError[bnX.bi()] -> SetPoint(graphABNoError[bnX.bi()]->GetN(),alpha,beta);

                      if (!iMM) {
                        histAlpha -> SetBinContent(bnX.bi(),bnY.bi(),alpha);
                        histMBeta -> SetBinContent(bnX.bi(),bnY.bi(),mbeta);
                        histCNorm -> SetBinContent(bnX.bi(),bnY.bi(),cnorm);
                        int dTOT2PC = (100*abs(temperatureArray[bnX.bi()][bnY.bi()][0] - temperatureArray[bnX.bi()][bnY.bi()][1]) / temperatureArray[bnX.bi()][bnY.bi()][0]);
                        histdT -> SetBinContent(bnX.bi(),bnY.bi(),dTOT2PC);
                        histK -> SetBinContent(bnX.bi(),bnY.bi(),alpha-mbeta);
                        histAlphaT -> SetBinContent(bnX.bi(),bnY.bi(),muNValue);
                        histMBetaT -> SetBinContent(bnX.bi(),bnY.bi(),-muZValue);
                      }
                    }
                  }

                  bnX.end();
                  while (bnX.prev())
                  {
                    TF1 *fitSep = nullptr;
                    if (setDrawNZInKT) fitSep = new TF1(Form("absep_linf_%s",nameR21AB),"-x+[0]",bnA.fMin,bnA.fMax);
                    else fitSep = new TF1(Form("absep_linf_%s",nameR21AB),"[1]*x+[0]",bnA.fMin,bnA.fMax);
                    graphABError[bnX.bi()] -> Fit(fitSep,"RQ0N");

                    if (setFitABSep && iComb!=3)
                    {
                      auto graphSep = new TGraphErrors(); 
                      graphSep -> SetPoint(0,bnA.fMin,fitSep->Eval(bnA.fMin));
                      graphSep -> SetPoint(1,bnA.fMax,fitSep->Eval(bnA.fMax));
                      auto errorSep = fitSep -> GetParError(0);
                      /*
                      if (!setDrawNZInKT) {
                        errorSep = 0.01;
                        if (bnX.bi()==1)
                          errorSep = 0;
                      }
                      */
                      graphSep -> SetPointError(0,0,errorSep);
                      graphSep -> SetPointError(1,0,errorSep);
                      setAtt(graphSep,bnX.ii(),3);

                      if (!setDrawFit3) {
                        graphSep -> SetPointError(0,0,bnA.getFullWidth()*0.01);
                        graphSep -> SetPointError(1,0,bnA.getFullWidth()*0.01);
                        cvsABMM -> cd();
                      }
                      graphSep -> Draw("same3");

                    }

                    setAtt(graphABError[bnX.bi()],bnX.ii(),2);
                    graphABError[bnX.bi()] -> SetMarkerStyle(10);
                    cvsABMM -> cd();
                    if (setDrawNZInKT) graphABError[bnX.bi()] -> Draw("sameezl");
                    else graphABError[bnX.bi()] -> Draw("sameezl");
                  }

                  TF1 *f1RefAB = nullptr;
                  f1RefAB = new TF1(Form("ref_%s",nameR21AB),"-x",bnA.fMin,bnA.fMax);
                  f1RefAB -> SetLineColor(kGray);
                  if (setDrawRefAB)
                  {
                    cvsABMM -> cd();

                    auto refxy0 = new TGraph;
                    refxy0 -> SetPoint(0,bnA.fMin,-bnA.fMin);
                    refxy0 -> SetPoint(1,bnA.fMax,-bnA.fMax);
                    refxy0 -> SetLineColor(kGray);
                    refxy0 -> SetLineStyle(2);
                    refxy0 -> Draw("samel");

                    auto refx0 = new TGraph();
                    refx0 -> SetPoint(0,0,bnB.fMin);
                    refx0 -> SetPoint(1,0,bnB.fMax);
                    refx0 -> SetLineColor(kGray);
                    refx0 -> SetLineStyle(2);
                    refx0 -> Draw("samel");

                    auto refy0 = new TGraph();
                    refy0 -> SetPoint(0,bnA.fMin,0);
                    refy0 -> SetPoint(1,bnA.fMax,0);
                    refy0 -> SetLineColor(kGray);
                    refy0 -> SetLineStyle(2);
                    refy0 -> Draw("samel");
                  }

                  cvsABMM -> cd();
                  for (auto mab : markers) mab -> Draw("p");

                  bnY.reset();
                  while (bnY.next()) {
                    auto markerAB = new TMarker(0,0,0);
                    setAtt(markerAB,bnY.ii(),0,2);
                    if (setDrawNZInKT) legendAB -> AddEntry(markerAB,Form("%s = %.0f ~ %.0f",bnY.getTitle(),bnY.cLow(),bnY.cHigh()),"p");
                    else legendAB -> AddEntry(markerAB,Form("%s = %.0f ~ %.0f",bnY.getTitle(),bnY.cLow(),bnY.cHigh()),"p");
                  }

                  bnX.reset();
                  while (bnX.next()) {
                    auto lineAB = new TGraph();
                    setAtt(lineAB,bnX.ii());
                    lineAB -> SetLineWidth(2);
                    if (setDrawNZInKT) legendAB -> AddEntry(lineAB,Form("%s = %.0f ~ %.0f",bnX.getTitle(),bnX.cLow(),bnX.cHigh()),"l");
                    else legendAB -> AddEntry(lineAB,Form("y_{0} = %.2f ~ %.2f",bnX.cLow(),bnX.cHigh()),"l");
                  }

                  cvsABMM -> cd();
                  makeLegend(cvsABMM,legendAB,"rt",0,0,.5,.35) -> Draw();
                }

                /////////////////////////////////////////////////////////////////////////////////////////////
                /////////////////////////////////////////////////////////////////////////////////////////////
                /////////////////////////////////////////////////////////////////////////////////////////////
                /////////////////////////////////////////////////////////////////////////////////////////////
                /////////////////////////////////////////////////////////////////////////////////////////////
                /////////////////////////////////////////////////////////////////////////////////////////////

                if (!setDrawNZAll)
                {
                  auto nameABYP = makeName("ab_yp",iAna,iLR,iMult,iSysComb,iCutTheta,iCutYP);
                  TCanvas *cvsAB;
                  if (setDrawABNZ||setDrawABTK)  {
                    cvsAB = makeCvs2(nameABYP,1000,1000);
                    cvsDivide(cvsAB,2,2);
                  }
                  else if (setDrawABC) {
                    cvsAB = makeCvs2(nameABYP,600,1100);
                    cvsDivide(cvsAB,1,3);
                  }
                  else {
                    cvsAB = makeCvs2(nameABYP,600,1100);
                    cvsDivide(cvsAB,1,2);
                  }
                  histAlpha -> SetMarkerSize(2);
                  histMBeta -> SetMarkerSize(2);
                  cvsAB -> cd(1); histAlpha -> Draw("textcolz");
                  cvsAB -> cd(2); histMBeta -> Draw("textcolz");

                  if (setDrawABNZ)  {
                    histAlphaT -> SetMarkerSize(2);
                    cvsAB -> cd(3);
                    histAlphaT -> Draw("textcolz");

                    histMBetaT -> SetMarkerSize(2);
                    cvsAB -> cd(4);
                    histMBetaT -> Draw("textcolz");
                  }
                  if (setDrawABTK)  {
                    histdT -> SetMarkerSize(2);
                    cvsAB -> cd(3);
                    histdT -> Draw("textcolz");

                    histK -> SetMarkerSize(2);
                    cvsAB -> cd(4);
                    histK -> Draw("textcolz");
                  }
                  if (setDrawABC) {
                    histCNorm -> SetMarkerSize(2);
                    cvsAB -> cd(3);
                    histCNorm -> Draw("textcolz");
                  }
                }

                TCanvas *cvsR21SN;
                if (setDrawNZAll && setDrawSN)
                {
                  // s(n) = R21(N,Z) exp(-beta Z)
                  const char *nameR21SN = makeName("r21sn",iAna,iLR,iMult,iSysComb,iCutTheta,iCutYP);
                  cvsR21SN = makeCvs(nameR21SN);
                  auto histSNN = (bnParticleN*bnR21).newHist(nameR21SN);
                  histSNN -> Draw();
                  auto graphSN = new TGraph();
                  graphSN -> SetMarkerSize(1.5);
                  graphSN -> SetMarkerStyle(20);
                  graphSN -> SetMarkerColor(kBlack);
                  for (auto iParticle : fParticleIdx) {
                    auto nValue = fParticleN[iParticle];
                    auto zValue = fParticleZ[iParticle];
                    auto r21Value = r21Array[1][1][iParticle] * exp(mbetaArray[1][1]*zValue);
                    graphSN -> SetPoint(graphSN->GetN(),nValue,r21Value);
                  }
                  graphSN -> Draw("samp");
                }

                if (setWriteR21ABCvsFile)
                {
                  const char *nameAlphaBetaRoot = Form("%s/rooto/alpha_beta_%d.root",fVName.Data(),iComb);
                  cout << nameAlphaBetaRoot << endl;
                  auto file = new TFile(nameAlphaBetaRoot,"recreate");
                  auto histA = (bnX*bnY).newHist("alpha");
                  auto histB = (bnX*bnY).newHist("beta");
                  auto histC = (bnX*bnY).newHist("cnorm");
                  auto histT = (bnX*bnY).newHist("temperature");
                  bnX.reset();
                  while (bnX.next()) {
                    bnY.reset();
                    while (bnY.next()) {
                      histA -> SetBinContent(bnX.fIdx, bnY.fIdx, alphaArray[bnX.fIdx][bnY.fIdx]);
                      histB -> SetBinContent(bnX.fIdx, bnY.fIdx, -mbetaArray[bnX.fIdx][bnY.fIdx]);
                      histC -> SetBinContent(bnX.fIdx, bnY.fIdx, cnormArray[bnX.fIdx][bnY.fIdx]);
                      histT -> SetBinContent(bnX.fIdx, bnY.fIdx, temperatureArray[bnX.fIdx][bnY.fIdx][2]);
                    }
                  }
                  histA -> Write();
                  histB -> Write();
                  histC -> Write();
                  histT -> Write();


                  auto nameR21YP = makeName("r21_yp",iAna,iLR,iMult,iSysComb,iCutTheta,iCutYP);
                  const char *nameAlphaBetaOut = Form("%s/others/alpha_beta_%s.txt",fVName.Data(),nameR21YP);
                  cout << nameAlphaBetaOut << endl;
                  std::ofstream fileAB(nameAlphaBetaOut);

                  fileAB << "cnorm" << endl;
                  fileAB << ","; bnX.reset(); while (bnX.next()) { fileAB << ", " << bnX.fIdx; } fileAB << endl;
                  fileAB << ","; bnX.reset(); while (bnX.next()) { fileAB << ", " << bnY0.fTitle << "=" << Form("%.2f",bnX.cLow()) << "~" << Form("%.2f",bnX.cHigh()); } fileAB << endl;
                  bnY.end(); while (bnY.prev()) {
                    fileAB << bnY.fIdx << ",  " << bnY.fTitle << "=" << int(bnY.cLow()) << "~" << int(bnY.cHigh());
                    bnX.reset(); while (bnX.next()) { fileAB << ", " << cnormArray[bnX.fIdx][bnY.fIdx]; } fileAB << endl;
                  } fileAB << endl;

                  fileAB << "temperature" << endl;
                  fileAB << ","; bnX.reset(); while (bnX.next()) { fileAB << ", " << bnX.fIdx; } fileAB << endl;
                  fileAB << ","; bnX.reset(); while (bnX.next()) { fileAB << ", " << bnY0.fTitle << "=" << Form("%.2f",bnX.cLow()) << "~" << Form("%.2f",bnX.cHigh()); } fileAB << endl;
                  bnY.end(); while (bnY.prev()) {
                    fileAB << bnY.fIdx << ",  " << bnY.fTitle << "=" << int(bnY.cLow()) << "~" << int(bnY.cHigh());
                    bnX.reset(); while (bnX.next()) { fileAB << ", " << temperatureArray[bnX.fIdx][bnY.fIdx][2]; } fileAB << endl;
                  } fileAB << endl;

                  fileAB << "alpha" << endl;
                  fileAB << ","; bnX.reset(); while (bnX.next()) { fileAB << ", " << bnX.fIdx; } fileAB << endl;
                  fileAB << ","; bnX.reset(); while (bnX.next()) { fileAB << ", " << bnY0.fTitle << "=" << Form("%.2f",bnX.cLow()) << "~" << Form("%.2f",bnX.cHigh()); } fileAB << endl;
                  bnY.end(); while (bnY.prev()) {
                    fileAB << bnY.fIdx << ",  " << bnY.fTitle << "=" << int(bnY.cLow()) << "~" << int(bnY.cHigh());
                    bnX.reset(); while (bnX.next()) { fileAB << ", " << alphaArray[bnX.fIdx][bnY.fIdx]; } fileAB << endl;
                  } fileAB << endl;

                  fileAB << "beta" << endl;
                  fileAB << ","; bnX.reset(); while (bnX.next()) { fileAB << ", " << bnX.fIdx; } fileAB << endl;
                  fileAB << ","; bnX.reset(); while (bnX.next()) { fileAB << ", " << bnY0.fTitle << "=" << Form("%.2f",bnX.cLow()) << "~" << Form("%.2f",bnX.cHigh()); } fileAB << endl;
                  bnY.end(); while (bnY.prev()) {
                    fileAB << bnY.fIdx << ",  " << bnY.fTitle << "=" << int(bnY.cLow()) << "~" << int(bnY.cHigh());
                    bnX.reset(); while (bnX.next()) { fileAB << ", " << -mbetaArray[bnX.fIdx][bnY.fIdx]; } fileAB << endl;
                  } fileAB << endl;

                  for (auto iParticle : fParticleIdx)
                  {
                    fileAB << "R21_" << fParticleNames[iParticle] << endl;
                    fileAB << ","; bnX.reset(); while (bnX.next()) { fileAB << ", " << bnX.fIdx; } fileAB << endl;
                    fileAB << ","; bnX.reset(); while (bnX.next()) { fileAB << ", " << bnY0.fTitle << "=" << Form("%.2f",bnX.cLow()) << "~" << Form("%.2f",bnX.cHigh()); } fileAB << endl;
                    bnY.end(); while (bnY.prev()) {
                      fileAB << bnY.fIdx << ",  " << bnY.fTitle << "=" << int(bnY.cLow()) << "~" << int(bnY.cHigh());
                      bnX.reset(); while (bnX.next()) { fileAB << ", " << r21Array[bnX.fIdx][bnY.fIdx][iParticle]; } fileAB << endl;
                    } fileAB << endl;
                  }

                  fileAB.close();
                }
              }
            }
          }
        }


      } // iLR

      if (drawDistTCM)
      {
        for (auto iParticle : fParticleIdx)
        {
          const char *nameParticle = fParticleNames[iParticle];

          for (auto iCutYP : selCutYPIdx)
          {
            if (selCutYP>=0 && selCutYP!=iCutYP) continue;

            {
              {
                auto nameDist = makeName("dist_tcm",iAna,0,iMult,4,0,iCutYP,iParticle);
                auto titleDist = makeTitle("",iAna,0,iMult,4,0,iCutYP,iParticle);
                auto cvs = makeCvs(nameDist,620,550);
                auto hist = (bnTCM*binning(100,0.001,1,"dN/d(TCM)/#Delta#Omega")).newHist(nameDist);
                cvs -> SetLogy();

                hist -> Draw();

                TLegend *legend = new TLegend();
                for (auto iLR : selLRIdx)
                {
                  if (selLR>=0 && selLR!=iLR) continue;
                  for (auto iCutTheta : selCutThetaIdx)
                  {
                    if (selCutTheta>=0 && selCutTheta!=iCutTheta) continue;

                    for (auto iSys : selSysIdx)
                    {
                      auto hist1 = histTCMArray[iLR][iSys][iCutTheta][iCutYP][iParticle];

                      double dPhi = 120;
                      double theta1 = TMath::DegToRad()*fCutThetaRanges[iCutTheta][0];
                      double theta2 = TMath::DegToRad()*fCutThetaRanges[iCutTheta][1];
                      auto solidAngleInPi = 2 * (dPhi/360) * f1Sinx -> Integral(theta1, theta2);
                      auto solidAngle = solidAngleInPi * TMath::Pi();
                      auto nbins = hist1 -> GetXaxis() -> GetNbins();
                      auto binWidth = hist1 -> GetXaxis() -> GetBinWidth(1);
                      for (auto bin=1; bin<=nbins; ++bin)
                        hist1 -> SetBinContent(bin,hist1->GetBinContent(bin)/solidAngle/binWidth);

                      setAtt(hist1,iSys);
                      if (iLR==kLeft) hist1 -> SetLineWidth(2);

                      hist1 -> Draw("samehist");
                      legend -> AddEntry(hist1,Form("%d %s",fSysBeams[iSys],fLRNames[iLR]),"l");
                    }
                  }
                }
                makeLegend(cvs,legend,"",0,0,0,0.30) -> Draw();
              }
            }
          }
        }
      }

      if (drawDistKeoa)
      {
        for (auto iParticle : fParticleIdx)
        {
          const char *nameParticle = fParticleNames[iParticle];

          for (auto iCutYP : selCutYPIdx)
          {
            if (selCutYP>=0 && selCutYP!=iCutYP) continue;

            //for (auto iComb : selSysCombIdx)
            {
              //if (selSysComb>=0 && selSysComb!=iComb) continue;

              //auto iSys2 = fSysCombIdx[iComb][0];
              //auto iSys1 = fSysCombIdx[iComb][1];
              //auto iSysComb = 100 + iSys2*10 + iSys1;

              //for (auto iSys : {iSys2, iSys1})
              {
                auto nameR21 = makeName("dist_keoa",iAna,0,iMult,4,0,iCutYP,iParticle);
                auto titleR21 = makeTitle("",iAna,0,iMult,4,0,iCutYP,iParticle);
                //auto nameR21 = makeName("dist_keoa",iAna,0,iMult,iSysComb,0,iCutYP,iParticle);
                //auto titleR21 = makeTitle("",iAna,0,iMult,iSysComb,0,iCutYP,iParticle);
                //auto nameR21 = makeName("dist_keoa",iAna,0,iMult,iSys,0,iCutYP,iParticle);
                //auto titleR21 = makeTitle("",iAna,0,iMult,iSys,0,iCutYP,iParticle);
                auto cvs = makeCvs(nameR21,620,550);
                //auto hist = new TH2D(nameR21,Form("%s;KE_{Lab}/A (MeV);dN/d(KE_{Lab}/A)/#Delta#Omega",titleR21),100,0,400.fMax,100,0.001,1);
                auto hist = new TH2D(nameR21,Form("%s;KE_{CM}/A (MeV);dN/d(KE_{CM}/A)/#Delta#Omega",titleR21),100,0,bnKeoaCM.fMax,100,0.001,1);
                cvs -> SetLogy();

                hist -> Draw();

                TLegend *legend = new TLegend();
                for (auto iLR : selLRIdx)
                {
                  if (selLR>=0 && selLR!=iLR) continue;
                  for (auto iCutTheta : selCutThetaIdx)
                  {
                    if (selCutTheta>=0 && selCutTheta!=iCutTheta) continue;

                    for (auto iSys : selSysIdx)
                    //for (auto iSys : {iSys2, iSys1})
                    {
                      auto hist1 = histKeoaArray[iLR][iSys][iCutTheta][iCutYP][iParticle];

                      double dPhi = 120;
                      double theta1 = TMath::DegToRad()*fCutThetaRanges[iCutTheta][0];
                      double theta2 = TMath::DegToRad()*fCutThetaRanges[iCutTheta][1];
                      auto solidAngleInPi = 2 * (dPhi/360) * f1Sinx -> Integral(theta1, theta2);
                      auto solidAngle = solidAngleInPi * TMath::Pi();
                      auto nbins = hist1 -> GetXaxis() -> GetNbins();
                      auto binWidth = hist1 -> GetXaxis() -> GetBinWidth(1);
                      for (auto bin=1; bin<=nbins; ++bin)
                        hist1 -> SetBinContent(bin,hist1->GetBinContent(bin)/solidAngle/binWidth);

                      setAtt(hist1,iSys);
                      if (iLR==kLeft) hist1 -> SetLineWidth(2);

                      hist1 -> Draw("samehist");
                      if (iCutTheta==kTheta0) legend -> AddEntry(hist1,Form("%d %s",fSysBeams[iSys],fLRNames[iLR]),"l");
                    }
                  }
                }
                makeLegend(cvs,legend,"",0,0,0,0.30) -> Draw();
              }
            }
          }
        }
      }

    } // iMult
  } // iAna

  if (saveCvsPNG)
    saveAll();

  if (saveCvsRoot)
    writeAll();
}

TLegend *makeLegend(TVirtualPad *cvs, TLegend *legend, TString fLegendDrawStyle="rt", double x_offset, double y_offset, double width_fixed, double height_fixed) // jumpto_makel
{
  double fWidthPerLengthLegend = 0.012;
  double fHeightPerEntryLegend = 0.05;
  double fWidthDefaultLegend = 0.15;
  int fFillStyleLegend = 0;
  int fBorderSizeLegend = 0;

  auto lmargin_cvs = fLMargin;
  auto rmargin_cvs = fRMargin;
  auto bmargin_cvs = fBMargin;
  auto tmargin_cvs = fTMargin;
  if (cvs != nullptr) {
    lmargin_cvs = cvs -> GetLeftMargin();
    rmargin_cvs = cvs -> GetRightMargin();
    bmargin_cvs = cvs -> GetBottomMargin();
    tmargin_cvs = cvs -> GetTopMargin();
  }
  else if (gPad != nullptr) {
    lmargin_cvs = ((TCanvas *) gPad) -> GetLeftMargin();
    rmargin_cvs = ((TCanvas *) gPad) -> GetRightMargin();
    bmargin_cvs = ((TCanvas *) gPad) -> GetBottomMargin();
    tmargin_cvs = ((TCanvas *) gPad) -> GetTopMargin();
  }

  auto x1_box =     lmargin_cvs;
  auto x2_box = 1.- rmargin_cvs;
  auto y1_box =     bmargin_cvs;
  auto y2_box = 1.- tmargin_cvs;
  auto unit_x = x2_box - x1_box;
  auto unit_y = y2_box - y1_box;

  auto y1_stat = y2_box;
  auto y2_stat = y2_box;
  auto x1_stat = x1_box;
  auto x2_stat = x2_box;
  auto dx_stat = x2_stat - x1_stat;
  auto dy_stat = y2_stat - y1_stat;

  if (fLegendDrawStyle=="xx1")
  {
    auto legend0 = legend;
    legend = new TLegend();
    auto count = 0;
    TIter next_entry(legend0->GetListOfPrimitives());
    while (auto label=(TLegendEntry*)next_entry()) {
      auto name = label -> GetLabel();
      auto obj = label -> GetObject();
      auto opt = label -> GetOption();
      if (count<5)
        legend -> AddEntry(obj,name,opt);
      count++;
    }
  }
  else if (fLegendDrawStyle=="xx2")
  {
    auto legend0 = legend;
    legend = new TLegend();
    auto count = 0;
    TIter next_entry(legend0->GetListOfPrimitives());
    while (auto label=(TLegendEntry*)next_entry()) {
      auto name = label -> GetLabel();
      auto obj = label -> GetObject();
      auto opt = label -> GetOption();
      if (count>=5)
        legend -> AddEntry(obj,name,opt);
      count++;
    }
  }

  auto length_text_max = 0;
  TIter next_entry(legend->GetListOfPrimitives());
  while (auto label=(TLegendEntry*)next_entry()) {
    auto lenth_text = TString(label->GetLabel()).Length();
    if (length_text_max<lenth_text)
      length_text_max = lenth_text;
  }

  auto y2_legend = y1_stat;
  auto y1_legend = y2_legend - legend->GetNRows() * fHeightPerEntryLegend;
  if (height_fixed>0) y1_legend = y2_legend - height_fixed*unit_y;

  auto x2_legend = x2_box;
  auto x1_legend = x2_legend - (length_text_max * fWidthPerLengthLegend) - fWidthDefaultLegend*unit_x;
  if (width_fixed>0) x1_legend = x2_legend - width_fixed*unit_x;

  if (fLegendDrawStyle=="rr")
  {
    x1_legend = x2_box;
    auto www = (length_text_max * fWidthPerLengthLegend) + fWidthDefaultLegend*unit_x;
    x2_legend = x1_legend + www;
    //cvs -> SetRightMargin(www);
  }
  if (fLegendDrawStyle=="lt")
  {
    x1_legend = x1_box;
    x2_legend = x1_legend + (length_text_max * fWidthPerLengthLegend) + fWidthDefaultLegend*unit_x;
    if (width_fixed>0) x2_legend = x1_legend + width_fixed*unit_x;
    if (x_offset==0) x_offset = 0.02;
  }
  else if (fLegendDrawStyle=="lb")
  {
    x1_legend = x1_box;
    x2_legend = x1_legend + (length_text_max * fWidthPerLengthLegend) + fWidthDefaultLegend*unit_x;
    if (width_fixed>0) x2_legend = x1_legend + width_fixed*unit_x;

    y1_legend = y1_box;
    y2_legend = y1_legend + legend->GetNRows() * fHeightPerEntryLegend;
    if (height_fixed>0) y1_legend = y2_legend - height_fixed*unit_y;
    if (x_offset==0) x_offset = 0.02;
    if (y_offset==0) y_offset = 0.04;
  }
  else if (fLegendDrawStyle=="fb")
  {
    x1_legend = x1_box;
    x2_legend = x2_box;

    y1_legend = y1_box;
    y2_legend = y1_legend + legend->GetNRows() * fHeightPerEntryLegend*1.2;
    if (height_fixed>0) y1_legend = y2_legend - height_fixed*unit_y;
    if (x_offset==0) x_offset = 0.02;
    if (y_offset==0) y_offset = 0.04;
  }

  else if (fLegendDrawStyle=="rb")
  {
    x2_legend = x2_box;
    x1_legend = x2_legend - (length_text_max * fWidthPerLengthLegend) - fWidthDefaultLegend*unit_x;
    if (width_fixed>0) x1_legend = x2_legend - width_fixed*unit_x;

    y1_legend = y1_box;
    y2_legend = y1_legend + legend->GetNRows() * fHeightPerEntryLegend;
    if (height_fixed>0) y1_legend = y2_legend - height_fixed*unit_y;
    if (y_offset==0) y_offset = 0.04;
  }
  else if (fLegendDrawStyle=="xx1")
  {
    x1_legend = x1_box;
    x2_legend = x1_legend + (length_text_max * fWidthPerLengthLegend) + fWidthDefaultLegend*unit_x;
    if (width_fixed>0) x2_legend = x1_legend + width_fixed*unit_x;

    y1_legend = y1_box;
    y2_legend = y1_legend + legend->GetNRows() * fHeightPerEntryLegend;
    if (height_fixed>0) y1_legend = y2_legend - height_fixed*unit_y;

    x_offset = 0.02;
    y_offset = 0.04;
  }
  else if (fLegendDrawStyle=="xx2")
  {
    x2_legend = x2_box;
    x1_legend = x2_legend - (length_text_max * fWidthPerLengthLegend) - fWidthDefaultLegend*unit_x;
    if (width_fixed>0) x1_legend = x2_legend - width_fixed*unit_x;

    y1_legend = y1_box;
    y2_legend = y1_legend + legend->GetNRows() * fHeightPerEntryLegend;
    if (height_fixed>0) y1_legend = y2_legend - height_fixed*unit_y;

    x_offset = 0;
    y_offset = 0.04;
  }

  x1_legend += x_offset*unit_x;
  x2_legend += x_offset*unit_x;
  y1_legend += y_offset*unit_y;
  y2_legend += y_offset*unit_y;

  auto dx_legend = x2_legend - x1_legend;

  if (dx_stat > 0) {
    if (abs(dx_legend-dx_stat)/dx_stat < 0.4)
      dx_legend = dx_stat;

    x1_legend = x2_legend - dx_legend;
  }

  legend -> SetX1(x1_legend);
  legend -> SetX2(x2_legend);
  legend -> SetY1(y1_legend);
  legend -> SetY2(y2_legend);

  legend -> SetFillStyle(fFillStyleLegend);
  legend -> SetBorderSize(fBorderSizeLegend);

  return legend;
}

double calculate_prob(double p_lab, double dedx)
{
  double probParticle = 0;
  double probSum = 0;
  for (auto iParticle : fParticleIdx) {
    double amp = graphFitPIDAmp[gIdxSys][gIdxParticle] -> Eval(p_lab);
    double mean =  f1PIDMean[gIdxSys][gIdxParticle] -> Eval(p_lab);
    double sigma = f1PIDSigma[gIdxSys][gIdxParticle] -> Eval(p_lab);
    auto prob = amp*TMath::Gaus(dedx, mean, sigma, true);
    probSum += prob;
    if (iParticle==gIdxParticle)
      probParticle = prob;
  }

  return (probParticle/probSum);
}

void project(TTree *tree, const char *name, const char *expr, TCut selection, bool verbose, TH1 *hist)
{
  tree -> Project(name,expr,selection);
  if (verbose) {
    if (hist!=nullptr)
      cout_info << tree -> GetName() << " -> " << name << " << " << expr << " @ " << selection << hist -> Integral() << endl;
    else
      cout_info << tree -> GetName() << " -> " << name << " << " << expr << " @ " << selection << endl;
  }
}

TCanvas *makeCvs2(const char *name, int w, int h) {
  const char *name0 = name;
  auto cvs = new TCanvas(name0,name0,cvsXOff+20*(countCvs+1), 20*(countCvs+1), w, h);
  cvs -> SetRightMargin(0.155);
  cvs -> SetLeftMargin(0.10);
  countCvs++;
  fCvsArray.push_back(cvs);
  return cvs;
}

TCanvas *makeCvs(const char *name, int w, int h) {
  const char *name0 = name;
  auto cvs = new TCanvas(name0,name0,cvsXOff+20*(countCvs+1), 20*(countCvs+1), w, h);
  cvs -> SetLeftMargin(0.10);
  countCvs++;
  fCvsArray.push_back(cvs);
  return cvs;
}

void savePDF() {
  for (auto cvs : fCvsArray) {
    cvs -> cd();
    cvs -> SaveAs(fVName+"/figures_pdf/"+cvs->GetName()+".pdf"); 
  }
}

void saveAll() {
  savePNG();
  savePDF();
}

void savePNG() {
  for (auto cvs : fCvsArray) {
    cvs -> cd();
    cvs -> SaveAs(fVName+"/figures/"+cvs->GetName()+".png"); 
  }
}

void writeAll()
{
  for (auto cvs : fCvsArray) {
    cvs -> cd();
    cvs -> SaveAs(fVName+"/rooto/"+cvs->GetName()+".root"); 
  }
}


void cvsDivideM0(TCanvas *cvs, int nx, int ny)
{
  // TODO
  auto sMargin = 0;
  cvs->Divide(nx,ny,0,0);
  for (auto ix=1; ix<=nx; ++ix) {
    for (auto iy=1; iy<=ny; ++iy) {
      auto i = ix + nx*(iy-1);
      if (iy==1&&iy==ny) {
        if (ix==1)       cvs->cd(i)->SetMargin(fLMargin,       0,fBMargin,0);
        else if (ix==nx) cvs->cd(i)->SetMargin(       0,fRMargin,fBMargin,0);
        else             cvs->cd(i)->SetMargin(       0,       0,fBMargin,0);
      }
      else if (iy==1) {
        if (ix==1&&ix==nx) cvs->cd(i)->SetMargin(fLMargin,fRMargin,sMargin,fTMargin);
        else if (ix==1)    cvs->cd(i)->SetMargin(fLMargin,       0,sMargin,fTMargin);
        else if (ix==nx)   cvs->cd(i)->SetMargin(       0,fRMargin,sMargin,fTMargin);
        else               cvs->cd(i)->SetMargin(       0,       0,sMargin,fTMargin);
      }
      else if (iy==ny) {
        if (ix==1&&ix==nx) cvs->cd(i)->SetMargin(fLMargin,fRMargin,fBMargin,sMargin);
        else if (ix==1)    cvs->cd(i)->SetMargin(fLMargin,       0,fBMargin,sMargin);
        else if (ix==nx)   cvs->cd(i)->SetMargin(       0,fRMargin,fBMargin,sMargin);
        else               cvs->cd(i)->SetMargin(       0,       0,fBMargin,sMargin);
      }
      else {
        if (ix==1&&ix==nx) cvs->cd(i)->SetMargin(fLMargin,fRMargin,sMargin,sMargin);
        else if (ix==1)    cvs->cd(i)->SetMargin(fLMargin,       0,sMargin,sMargin);
        else if (ix==nx)   cvs->cd(i)->SetMargin(       0,fRMargin,sMargin,sMargin);
        else               cvs->cd(i)->SetMargin(       0,       0,sMargin,sMargin);
     }
    }
  }
}

void cvsDivide(TCanvas *cvs, int nx, int ny)
{
  auto lMargin = cvs -> GetLeftMargin();
  auto rMargin = cvs -> GetRightMargin();
  auto bMargin = cvs -> GetBottomMargin();
  auto tMargin = cvs -> GetTopMargin();

  cvs -> Divide(nx,ny,0.01,0.001);

  for (auto ix=1; ix<=nx; ++ix) {
    for (auto iy=1; iy<=ny; ++iy) {
      auto i = ix + nx*(iy-1);
      cvs -> cd(i) -> SetMargin(lMargin,rMargin,bMargin,tMargin);
    }
  }
}

const char *makeName(const char *mainName, int iAna, int iLR, int iMult, int iSys, int iCutTheta, int iCutYP, int iPart)
{
  TString sdName = fSDNames[fSelSDType];
  sdName.ReplaceAll("SDVALUE",Form("%.1f",fSDValue));
  sdName.ReplaceAll(".","p");

  //mainName = Form("%s_%s_",mainName,sdName.Data());
  //mainName = Form("%s_%s_",mainName,sdName.Data());

  if (fUseHandCut)
    mainName = Form("%s_HPID_",mainName);

  const char *systemName;
  if (iSys>=100) {
    iSys = iSys - 100;
    int iSys2 = int(iSys/10);
    int iSys1 = iSys - 10*iSys2;
    systemName = Form("%so%s",fSysNames[iSys2],fSysNames[iSys1]);
  }
  else
    systemName = Form("%s",fSysNames[iSys]);

  if (iPart<0) {
    const char *name = Form("%s_%s_%s_%s_%s_%s_%s_%s",mainName,fAnaNames[iAna],sdName.Data(),fLRNames[iLR],fMultNames[iMult],systemName,fCutThetaNames[iCutTheta],fCutYPNames[iCutYP]);
    return name;
  }

  const char *name = Form("%s_%s_%s_%s_%s_%s_%s_%s_%s",mainName,fAnaNames[iAna],sdName.Data(),fLRNames[iLR],fMultNames[iMult],systemName,fCutThetaNames[iCutTheta],fCutYPNames[iCutYP],fParticleNames[iPart]);
  return name;
}

const char *makeTitle(const char *mainName, int iAna, int iLR, int iMult, int iSys, int iCutTheta, int iCutYP, int iPart)
{
  if (fUseHandCut)
    mainName = Form("%s, Hand-Cut-PID",mainName);

  const char *systemTitle;
  if (iSys>=100) {
    iSys = iSys - 100;
    int iSys2 = int(iSys/10);
    int iSys1 = iSys - 10*iSys2;
    //systemTitle = Form("%d/%d",fSysBeams[iSys2],fSysBeams[iSys1]);
    systemTitle = Form("%s/%s",fSysTitles[iSys2],fSysTitles[iSys1]);
  }
  else
    systemTitle = Form("%s",fSysTitles[iSys]);


  const char *partTitle = "";
  if (iPart!=5)
    partTitle = Form(", %s",fParticleNames[iPart]);

  const char *y0Title = "";
  if (iCutYP!=0)
    y0Title = Form(", %s",fCutYPTitles[iCutYP]);

  TString sdTitle = fSDTitles[fSelSDType];
  sdTitle.ReplaceAll("SDVALUE",Form("%.1f",fSDValue));
  sdTitle.ReplaceAll("ANATTL2",fAnaTitles2[fSelAna]);

  //const char *title = Form("%s, %s, %s, %s, %s, %s%s%s",mainName,fAnaTitles[iAna],fLRTitles[iLR],fMultTitles[iMult],systemTitle,fCutThetaTitles[iCutTheta],y0Title,partTitle);
  //const char *title = Form("%s",fCutThetaTitles[iCutTheta]); if (fUseHandCut) title = Form("Hand-Cut-PID, %s",mainName);
  //const char *title = Form("%s, %s",fMultTitles[iMult],systemTitle); if (fUseHandCut) title = Form("Hand-Cut-PID, %s",mainName);
  //const char *title = Form("%s",fCutYPTitles[iCutYP]); if (fUseHandCut) title = Form("Hand-Cut-PID, %s",mainName);
  //const char *title = Form("%s, %s, %s, %s",mainName,fMultTitles[iMult],systemTitle,fCutThetaTitles[iCutTheta]);
  //const char *title = Form("%s, %s, %s",mainName,fMultTitles[iMult],systemTitle);
  //const char *title = Form("%s, %s%s",fMultTitles[iMult],systemTitle,partTitle);
  const char *title = Form("%s, %s, %s, %s%s%s",sdTitle.Data(), fLRTitles[iLR], fMultTitles[iMult],systemTitle,y0Title,partTitle);

  return title;

}

void setAtt(TH1 *hist, int iDraw) {
  hist -> SetLineColor(fDrawColor[iDraw]);
  hist -> SetMarkerColor(fDrawColor[iDraw]);
  hist -> SetMarkerStyle(fDrawMStyle[iDraw]);
  hist -> SetMarkerSize(fDrawMSize[iDraw]);
}

void setAtt(TGraph *graph, int iDraw, int option) {
  if (option==3) {
    graph -> SetLineColor(fDrawColor[iDraw]);
    graph -> SetMarkerColor(fDrawColor[iDraw]);
    graph -> SetFillColor(fDrawColor3[iDraw]);
    graph -> SetMarkerStyle(fDrawMStyle2[iDraw]);
    graph -> SetMarkerSize(fDrawMSize[iDraw]);
  }
  else if (option==2) {
    graph -> SetLineColor(fDrawColor2[iDraw]);
    graph -> SetMarkerColor(fDrawColor2[iDraw]);
    graph -> SetFillColor(fDrawColor3[iDraw]);
    graph -> SetMarkerStyle(fDrawMStyle2[iDraw]);
    graph -> SetMarkerSize(fDrawMSize2[iDraw]);
  }
  else if (option==4) {
    graph -> SetLineColor(fDrawColor[iDraw]);
    graph -> SetMarkerColor(fDrawColor[iDraw]);
    graph -> SetFillColor(fDrawColor[iDraw]);
    graph -> SetMarkerStyle(fDrawMStyle3[iDraw]);
    graph -> SetMarkerSize(fDrawMSize[iDraw]);
  }
  else { //3?
    graph -> SetLineColor(fDrawColor[iDraw]);
    graph -> SetMarkerColor(fDrawColor[iDraw]);
    graph -> SetFillColor(fDrawColor3[iDraw]);
    graph -> SetMarkerStyle(fDrawMStyle[iDraw]);
    graph -> SetMarkerSize(fDrawMSize3[iDraw]);
  }
}

void setAtt(TMarker *marker, int iDraw, int iColor, int option) {
  if (iColor<0) iColor = iDraw;

  if (option==2) {
    marker -> SetMarkerColor(fDrawColor[iColor]);
    marker -> SetMarkerStyle(fDrawMStyle2[iDraw]);
    marker -> SetMarkerSize(fDrawMSize2[iDraw]);
    if (setDrawNZAll) marker -> SetMarkerSize(2*fDrawMSize2[iDraw]);
  }
  else {
    marker -> SetMarkerColor(fDrawColor[iColor]);
    marker -> SetMarkerStyle(fDrawMStyle[iDraw]);
    marker -> SetMarkerSize(fDrawMSize[iDraw]);
    if (setDrawNZAll) marker -> SetMarkerSize(2*fDrawMSize[iDraw]);
  }
}

void setAtt(TF1 *fit, int iDraw) {
  fit -> SetLineColor(fDrawColor[iDraw]);
}
