#include "KBGlobal.hh"
#include "init_variables.h"
#include "binning.h"

TString fVName;
vector<TCanvas *> fCvsArray;
int cvsXOff = 1000;
int gIdxSys = 0;
int gIdxParticle = 0;
int countCvs = 0;
bool fUseHandCut = false;
double fTMargin = 0.12;
double fBMargin = 0.20;
double fLMargin = 0.19;
double fRMargin = 0.055;
double fRMargin1 = 0.055;

int fSelAna = 0;
int fSelSDType = 0;
double fSDValue = 3;

TGraph *graphFitPIDMean[fNumSyss][fNumParticles] = {0};
TGraph *graphFitPIDAmp[fNumSyss][fNumParticles] = {0};
TH1F *histPIDMeta[fNumSyss][fNumParticles] = {0};
TF1 *f1PIDMean[fNumSyss][fNumParticles] = {0};
TF1 *f1PIDSigma[fNumSyss][fNumParticles] = {0};
TCutG *cutgPID[fNumSyss][fNumParticles] = {0};
TGraph *graphPIDMean[fNumSyss][fNumParticles] = {0};
TGraph *graphPIDRange[fNumSyss][fNumParticles][2] = {0};

void saveAll();
void writeAll();
double calculate_prob(double p_lab, double dedx);
void project(TTree *tree, const char *name, const char *expr, TCut selection, bool verbose=1);
TCanvas *makeCvs2(const char *name, int w=1050, int h=950);
TCanvas *makeCvs(const char *name, int w=630, int h=550);
void cvsDivide(TCanvas *cvs, int nx, int ny);
void cvsDivideM0(TCanvas *cvs, int nx, int ny);
TLegend *makeLegend(TVirtualPad *cvs, TLegend *legend, TString opt = "", double x_offset=0, double y_offset=0, double width_fixed=0, double height_fixed=0);
void setDrawAtt(TH1 *hist, int iDraw);
void setDrawAtt(TGraph *graph, int iDraw);
void setDrawAtt(TF1 *fit, int iDraw);
void setDrawAtt(TMarker *marker, int iDraw, int iColor);
const char *makeName(const char *mainName, int iAna, int iLR, int iMult, int iSys, int iCutTheta, int iCutYP, int iPart=-1);
const char *makeTitle(const char *mainName, int iAna, int iLR, int iMult, int iSys, int iCutTheta, int iCutYP, int iPart=5);



void draw_pid
(
    //int selAna = kF132, double selSDValue = 0.8,
    //int selAna = kF132, double selSDValue = 1.0,
    int selAna = kFSys, double selSDValue = 3.0,

    int selSDType = kSD_xx,
    //int selSDType = kSD_xx_l3,
    //int selSDType = kSD_0_x,

    //int selCutYP0 = kYPAll,
    //int selCutYP0 = -2,
    int selCutYP0 = kYPF3,
    bool checkContent = false,
    bool saveCvsPNG = false,
    bool saveCvsRoot = false
)
{
  if (selAna>=0) fSelAna = selAna;
  if (selSDType>=0) fSelSDType = selSDType;
  if (selSDValue>0) fSDValue = selSDValue;

  TCut cut0 = "prob>.7"; fVName = Form("v%s_sd%.1f_p0p7",fAnaNames[selAna],selSDValue);
  //TCut cut0 = "prob>.7"; fVName = "vExtreme";
  //TCut cut0 = "prob>.6"; fVName = "vProb6";
  //TCut cut0 = ""; fVName = "vAll";

  TString probString = "prob";

  ///////////////////////////////////////////////////////////////////////////////////////////////////

  bool drawRawPID = kHold;
  bool setDrawRawPIDProj = kUnSet;
  bool setDrawRawPIDSub = kSet;
  bool setFitdEdx = kSet;

  bool drawCorPID = kHold;
  bool setDrawCorPIDInRaw = kSet;
  bool setDrawCorPIDInPtoa = kUnSet;

  bool setDrawGuideLine = kSet;
  bool setDrawGuideLine132 = kSet;
  bool setDrawHandCut = fUseHandCut;
  bool setZoomPID0 = kSet;
  bool setZoomPIDZ1 = kUnSet;
  bool setZoomPIDZ2 = kUnSet;

  bool drawCorEKE = kHold;
  bool drawKT = kHold;
  bool drawPtoa = kHold;

  bool drawYP = kHold;
  bool setDrawYPTogether = kSet;
  bool setDrawYPTheta0 = kUnSet;
  bool setDrawYPKE = kUnSet;
  bool setDrawYPPoz = kUnSet;
  bool setDrawYPText = kUnSet;
  bool setDrawYPGrid = kSet;

  bool drawPtoaR21 = kHold;
  bool drawY0R21 = kHold;
  bool setDrawR21 = kSet;
  bool setDrawSR = kUnSet;
  bool setDrawDR = kUnSet;
  bool setDrawTemp = kSet;
  bool setDrawR21AMD = kSet;

  bool drawNZR21 = kDraw;
  bool setFitNZR21TG = kSet;
  bool setDrawInvB = kUnSet;
  bool setFitABSep = kUnSet;
  bool setFitAB = kUnSet;
  bool setWriteR21ABCvsFile = kSet;
  bool setDrawABError = kSet;
  bool setABErrorToResidual = kSet;
  bool setDrawNZLogy = kUnSet;
  bool setDrawNZAll = kSet;
  bool setDrawSN = kSet;

  bool drawKeoaR21 = kHold;

  bool drawDistKeoa = kHold;

  auto holdAll = [&drawRawPID, &drawCorPID, &drawCorEKE, &drawKT, &drawPtoa, &drawYP, &drawPtoaR21, &drawY0R21, &drawKeoaR21, &drawNZR21, &drawDistKeoa]() {
    drawRawPID   = kHold;
    drawCorPID   = kHold;
    drawCorEKE   = kHold;
    drawKT       = kHold;
    drawPtoa     = kHold;
    drawYP       = kHold;
    drawPtoaR21  = kHold;
    drawY0R21    = kHold;
    drawKeoaR21  = kHold;
    drawNZR21    = kHold;
    drawDistKeoa = kHold;
  };

  // sssssssssssssssss //////////////////////////////////////////////////////////////////////////////

  int selMult = kMult55;

  int selSys = kAll;
  vector<int> selSysIdx = {k132,k108,k112,k124};
  //vector<int> selSysIdx = {k132,k108,k112};

  int selLR = kLR;
  vector<int> selLRIdx = {kLR,kLeft,kRight};

  int selCutTheta = kThetaAll;
  vector<int> selCutThetaIdx = {kThetaAll, kTheta0, kTheta20, kTheta40, kTheta60, kThetaLT60, kThetaGT60};

  int selCutYP = kYPAll;
  //int selCutYP = kYPF4;
  vector<int> selCutYPIdx = {kYPAll, kY1, kY2, kY3, kY4, kPtoa0, kPtoa50, kPtoa100, kPtoa150, kPtoa200, kPtoa250, kPtoa300, kPtoa350, kYF, kPF3, kYPF3, kPF4, kYPF4};

  //int selCutYP = kAll;
  //vector<int> selCutYPIdx = {kY2, kY3, kY4};

  if (selCutYP0 >= 0)
    selCutYP = selCutYP0;

  //int selSysComb = kAll;
  int selSysComb = 0;
  vector<int> selSysCombIdx = {0,1,2,3};

  if (drawDistKeoa) {
    selCutThetaIdx.clear();
    selCutTheta = kAll;
    for (auto i : {1,2,3,4})
      selCutThetaIdx.push_back(i);

    holdAll();
    drawDistKeoa = true;
    selSysIdx.clear();
    selSysIdx.push_back(k132);
    selSysIdx.push_back(k112);
    selSysIdx.push_back(k108);
  }

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

  int nbinsFrame = 100;

  binning nixPIDProjX(200,800,820); // number of total projections, min x, max y
  binning bnPoz(800,0,3000,"p/Z (MeV/c)");
  binning bndEdx(800,0,2000,"dE/dx");
  binning bnPIDz(0,0.08,0.0001);
  if (setDrawCorPIDInRaw) {
    //bnPIDz.set(1,1,800);
    bnPIDz.set(1,1,50000);
  }

  if (setZoomPIDZ2) {
    bnPoz.set(300,200,2200);
    bndEdx.set(300,0,800);
  }
  else if (setZoomPIDZ1) {
    bnPoz.set(400,0,3000);
    bndEdx.set(400,0,200);
  }
  else if (setZoomPID0) {
    bnPoz.set(200,0,2000);
    bndEdx.set(200,0,800);
  }

  
  
  binning bnY0(100,-1,2,"y_{0}");
  binning bnPtoa(100,0,800,"p_{T}/A (MeV/c)");
  binning bnKeoa(100,0,400,"KE_{Lab}/A (MeV)");
  binning bnKeoa2(100,0,1500,"KE/A (MeV)");
  binning bnTheta(100,0,90,"#theta_{Lab} (deg)");

  binning bnParticleN(100,-1,3,"N");
  binning bnR21(100,0,2.0,"R_{21}");
  //binning bnRY0(5,-0.25,1); // selected;
  binning bnRY0(3,-0.25,1.25,"y_{0}"); // selected;
  binning bnRPtoa(8,0,400);
  //binning bnRPtoa(10,0,500);

  if (drawNZR21 && setDrawNZAll) {
    holdAll();
    drawNZR21  = true;
    bnRY0.set(1,-0.25,0.75); // selected;
    bnRPtoa.set(1,0,300);
  }

  binning bnRKeoa(8,0,400);

  if (drawDistKeoa)
    bnRKeoa.fN = 50;

  if (setDrawTemp) {
    //bnRY0.set(6,-.5,1.0);
    //bnRPtoa.set(10,0,400);
  }

  if (setDrawSR || setDrawDR) bnRY0.set(21,-.8,1.3);

  binning bnPR(100,0,1.1,"SR");
  binning bnDR(100,1,2.2,"DR");
  binning bnTemp(100,5,15,"T_{app}");

  binning bnAlpha(100,0,.4,"#alpha");
  binning bnBeta(100,-.4,0,"#beta");

  if (setDrawInvB) {
    binning bnBeta(100,0,0.4);
  }

  double scaleMaxFitdEdx = .05;

  ///////////////////////////////////////////////////////////////////////////////////////////////////

  cout << endl;
  cout_info << "== " << fVName << " ; " << " " << cut0 << " ; " << probString << endl;

  if (drawRawPID  ) cout_info << "set draw RawPID  " << endl;
  if (drawCorPID  ) cout_info << "set draw CorPID  " << endl;
  if (drawCorEKE  ) cout_info << "set draw CorEKE  " << endl;
  if (drawKT      ) cout_info << "set draw KT      " << endl;
  if (drawPtoa    ) cout_info << "set draw Ptoa    " << endl;
  if (drawYP      ) cout_info << "set draw YP      " << endl;
  if (drawPtoaR21 ) cout_info << "set draw PtoaR21 " << endl;
  if (drawY0R21   ) cout_info << "set draw Y0R21"    << endl;
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
  gSystem -> mkdir(fVName+"/others");

  auto getProbString = [setDrawRawPIDProj,probString](int iSys, int iParticle) {
    const char *finalProbString = probString;
    if (setDrawRawPIDProj) {
      gIdxSys = iSys;
      gIdxParticle = iParticle;
      finalProbString = "calculate_prob(p_lab,dedx)";
    }
    return finalProbString;
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
      //graphPIDMean[iSys][iParticle] -> SetLineColor(kGray+2);
      graphPIDMean[iSys][iParticle] -> SetLineColor(kRed);
      if (iSys==0) graphPIDMean[iSys][iParticle] -> SetLineColor(kBlack);
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
                /*
                if (setZoomPIDZ1) {
                  namePIDraw = makeName("pidRawZoomed",iAna,iLR,iMult,iSys,iCutTheta,iCutYP);
                  titlePIDraw = makeTitle("Raw PID (Zoomed)",iAna,iLR,iMult,iSys,iCutTheta,iCutYP);
                }
                */

                auto histPID = (bnPoz*bndEdx).newHist(namePIDraw,titlePIDraw);
                //histPID -> SetMinimum(0.5);
                //histPID -> SetMaximum(800);
                project(treeAll,namePIDraw,"dedx:p_lab",cutTheta);

                TCanvas *cvsPIDraw = nullptr;
                if (drawRawPID) {
                  cvsPIDraw = makeCvs2(namePIDraw,1000,700);
                  //cvsPIDraw = makeCvs2(namePIDraw,950,850);
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
                  auto binn = binning(histPID);
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
                        setDrawAtt(graph,iParticle);

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

        int numEventsInAna[fNumSyss] = {0};
        TH1D *histY0Array[fNumSyss][fNumCutThetas][fNumCutYPs][fNumParticles] = {0};
        TH1D *histPtoaArray[fNumSyss][fNumCutThetas][fNumCutYPs][fNumParticles] = {0};
        TH2D *histYPR21Array[fNumSyss][fNumCutThetas][fNumCutYPs][fNumParticles] = {0};

        for (auto iSys : selSysIdx)
        {
          if (selSys>=0 && selSys!=iSys) continue;

          auto sys = fSysBeams[iSys];
          const char *sysTitle = fSysTitles[iSys];

          if (numEventsInAna[iSys]==0) {
            auto nameFileMult = Form("data2/%s/sys%d_%s_%s_%s.%s.ana.mult.root",anaFName,sys,anaFName,lrFName,multFName,spVersion);

            auto treeMult = new TChain("mult");
            treeMult -> Add(nameFileMult);

            numEventsInAna[iSys] = treeMult -> GetEntries("1");
            cout_info << "Event Multiplicity in " << sys << " = " << numEventsInAna[iSys] << "  (" << nameFileMult << ")" << endl;

            delete treeMult;
          }

          if (drawKT || drawYP)
          {
            for (auto iCutYP : selCutYPIdx)
            {
              if (selCutYP>=0 && selCutYP!=iCutYP) continue;
              TString stringCutYP0 = fCutYPValues[iCutYP].GetTitle();

              auto iCutTheta = 0;
              TCut cutTheta = fCutThetaValues[iCutTheta];

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
                TCut cutParticleYP = stringCutYP.Data();

                const char *nameParticle = fParticleNames[iParticle];
                const char *titleParticle = fParticleTitles[iParticle];
                auto nameFileParticle = Form("data2/%s/sys%d_%s_%s_%s.%s.ana.%s.root",anaFName,sys,anaFName,lrFName,multFName,spVersion,nameParticle);

                auto treeParticle = new TChain(nameParticle);
                treeParticle -> Add(nameFileParticle);

                if (drawKT)
                {
                  auto nameKTCorPart = makeName("ktCor",iAna,iLR,iMult,iSys,iCutTheta,iCutYP,iParticle);
                  auto titleKTCorPart = makeTitle("Cor.",iAna,iLR,iMult,iSys,iCutTheta,iCutYP,iParticle);

                  auto histKTCorPart = (bnKeoa*bndEdx).newHist(nameKTCorPart,titleKTCorPart);
                  TCut selection = cut0 * TCut(Form("%s/eff/%d",getProbString(iSys,iParticle),numEventsInAna[iSys])) * cutParticleYP * cutTheta * cutParticlePoz * cutParticleSD;
                  auto partz = fParticleZ[iParticle];
                  auto partm = fParticleMass[iParticle];
                  auto parta = fParticleA[iParticle];
                  const char *expression = Form("theta_lab*TMath::RadToDeg():(sqrt((p_lab*%d)*(p_lab*%d)+%f*%f)-%f)/%d",partz,partz,partm,partm,partm,parta);
                  project(treeParticle,nameKTCorPart,expression,selection);

                  if (cvsKT==nullptr) {
                    cvsKT = makeCvs2(nameKTCorPart,1000,700);
                    cvsDivide(cvsKT,3,2);
                  }

                  cvsKT -> cd(iParticle+1);
                  histKTCorPart -> Draw("colz");
                }

                if (drawYP)
                {
                  auto nameYPPart = makeName("yp",iAna,iLR,iMult,iSys,iCutTheta,iCutYP,iParticle);
                  auto titleYPPart = makeTitle("",iAna,iLR,iMult,iSys,iCutTheta,iCutYP,iParticle);

                  auto histYPPart  = (bnY0*bnPtoa).newHist(nameYPPart,titleYPPart);
                  TCut selection = cut0 * TCut(Form("%s/eff/%d",getProbString(iSys,iParticle),numEventsInAna[iSys])) * cutParticleYP * cutTheta * cutParticlePoz * cutParticleSD;
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
                    double par0 = 3, par1 = -20;
                    if (0) {
                      for (auto iCutTheta0 : {2,4,6})
                      {
                        int theta0 = 10*iCutTheta0;
                        auto nameYPPart2 = makeName(Form("yptta%d",theta0),iAna,iLR,iMult,iSys,iCutTheta,iCutYP,iParticle);

                        auto histYPThetaPart = (bnY0*bnPtoa).newHist(nameYPPart2,titleYPPart);
                        TCut cutTheta0 = Form("theta_lab>= %d*TMath::DegToRad()&&theta_lab<%d*TMath::DegToRad()",theta0-1,theta0+1); //fCutTheta0Values[iCutTheta0];
                        selection = cut0 * TCut(Form("%s/eff/%d",getProbString(iSys,iParticle),numEventsInAna[iSys])) * cutParticleYP * cutTheta * cutTheta0 * cutParticlePoz * cutParticleSD;
                        project(treeParticle,nameYPPart2,Form("pt_cm/%d:fy_cm/(by_cm/2)",fParticleA[iParticle]),selection,0);
                        if (iCutTheta0==2) histYPThetaPart -> Draw("colz");
                        else histYPThetaPart -> Draw("samecol");

                        auto nameYPPart3 = makeName(Form("fityptta%d",theta0),iAna,iLR,iMult,iSys,iCutTheta,iCutYP,iParticle);
                        TF1 *fitPYTT0 = new TF1(nameYPPart3,"pol3",bnY0.fMin,bnY0.fMax);
                        fitPYTT0 -> SetParameter(par0, par1);
                        histYPThetaPart -> Fit(fitPYTT0,"RQ0");
                        cout << "if (iParticle==" << iParticle << "&&iCutTheta0==" << iCutTheta0 << ") fitPYTT0 -> SetParameters(" << fitPYTT0 -> GetParameter(0)
                          << ", " << fitPYTT0 -> GetParameter(1)
                          << ", " << fitPYTT0 -> GetParameter(2) 
                          << ", " << fitPYTT0 -> GetParameter(3) << ");" << endl;
                        par0 = fitPYTT0 -> GetParameter(0);
                        par1 = fitPYTT0 -> GetParameter(1);
                        fitPYTT0 -> SetLineStyle(2);
                        fitPYTT0 -> SetLineColor(kSpring-6);
                        fitPYTT0 -> Draw("samel");

                        if (setDrawYPText) {
                          auto xAt = bnY0.fMax;
                          auto yAt = fitPYTT0 -> Eval(xAt);

                          if (yAt<bnPtoa.fMax) {
                            TLatex *text = new TLatex(xAt,yAt,Form("#theta_{Lab}=%d#circ",theta0));
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
                      selection = cut0 * TCut(Form("%s/eff/%d",getProbString(iSys,iParticle),numEventsInAna[iSys])) * cutParticleYP * cutTheta * cutPoz * cutParticlePoz * cutParticleSD;
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
                      const char *exprKE = Form("(sqrt((p_lab*%d)*(p_lab*%d)+%f*%f)-%f)",partz,partz,partm,partm,partm);
                      auto dKE = (iKE/2.)*5;
                      TCut cutKE = Form("%s>%f&&%s<%f",exprKE,iKE*100-dKE,exprKE,iKE*100+dKE);

                      auto nameYPPart2 = makeName(Form("ypke%d",iCutKE),iAna,iLR,iMult,iSys,iCutTheta,iCutYP,iParticle);

                      auto histYPKEPart = (bnY0*bnPtoa).newHist(nameYPPart2,titleYPPart);
                      selection = cut0 * TCut(Form("%s/eff/%d",getProbString(iSys,iParticle),numEventsInAna[iSys])) * cutParticleYP * cutTheta * cutKE * cutParticlePoz * cutParticleSD;
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


          if (drawKT || drawCorEKE || drawCorPID || drawPtoaR21 || drawDistKeoa || drawKeoaR21 || drawY0R21 || drawNZR21)
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
                    TCut selection = cut0 * cutParticleYP * cutTheta * TCut(Form("%s/eff/%d",getProbString(iSys,iParticle),numEventsInAna[iSys])) * cutParticlePoz * cutParticleSD;

                    if (setDrawCorPIDInRaw)
                      selection = cut0 * cutParticleYP * cutTheta * cutParticlePoz * cutParticleSD;

                    project(treeParticle,namePIDCorPart,"dedx:p_lab",selection);

                    if (histPIDCor==nullptr) {
                      auto namePIDCor = makeName("pidCor",iAna,iLR,iMult,iSys,iCutTheta,iCutYP);
                      histPIDCor = (TH2D *) histPIDCorPart -> Clone(namePIDCor);
                      histPIDCor -> SetTitle(Form("%s;p/Z (MeV/c);dE/dx;",titlePIDCorPart));
                    }
                    else
                      histPIDCor -> Add(histPIDCorPart);
                  }

                  if (drawCorEKE)
                  {
                    auto nameEKECorPart = makeName("ekeCor",iAna,iLR,iMult,iSys,iCutTheta,iCutYP,iParticle);
                    auto titleEKECorPart = makeTitle("Cor.",iAna,iLR,iMult,iSys,iCutTheta,iCutYP,iParticle);

                    auto histEKECorPart = (bnKeoa2*bndEdx).newHist(nameEKECorPart,titleEKECorPart);
                    TCut selection = cut0 * TCut(Form("%s/eff/%d",getProbString(iSys,iParticle),numEventsInAna[iSys])) * cutParticleYP * cutTheta * cutParticlePoz * cutParticleSD;
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

                  if (drawPtoaR21)
                  {
                    auto namePtoaPart = makeName("ptoa",iAna,iLR,iMult,iSys,iCutTheta,iCutYP,iParticle);
                    auto titlePtoaPart = makeTitle("Cor.",iAna,iLR,iMult,iSys,iCutTheta,iCutYP,iParticle);

                    auto histPtoa = bnRPtoa.newHist(namePtoaPart,titlePtoaPart);
                    TCut selection = cut0 * TCut(Form("%s/eff/%d",getProbString(iSys,iParticle),numEventsInAna[iSys])) * cutParticleYP * cutTheta  * cutParticlePoz * cutParticleSD;
                    if (fUseHandCut) selection = cut0 * TCut(Form("%s/eff/%d","1",numEventsInAna[iSys])) * getSelBoundHandCut(iSys,iParticle) * cutParticleYP * cutTheta * cutParticlePoz * cutParticleSD;

                    project(treeParticle,namePtoaPart,Form("pt_cm/%d",fParticleA[iParticle]),selection);

                    histPtoaArray[iSys][iCutTheta][iCutYP][iParticle] = histPtoa;
                  }

                  if (drawDistKeoa || drawKeoaR21)
                  {
                    auto nameKeoaPart = makeName("keoa",iAna,iLR,iMult,iSys,iCutTheta,iCutYP,iParticle);
                    auto titleKeoaPart = makeTitle("Cor.",iAna,iLR,iMult,iSys,iCutTheta,iCutYP,iParticle);

                    auto histKeoa = bnRKeoa.newHist(nameKeoaPart,titleKeoaPart);
                    TCut selection = cut0 * TCut(Form("%s/eff/%d",getProbString(iSys,iParticle),numEventsInAna[iSys])) * cutParticleYP * cutTheta * cutParticlePoz * cutParticleSD;
                    if (fUseHandCut) selection = cut0 * TCut(Form("%s/eff/%d","1",numEventsInAna[iSys])) * getSelBoundHandCut(iSys,iParticle) * cutParticleYP * cutTheta * cutParticlePoz * cutParticleSD;
                    auto partz = fParticleZ[iParticle];
                    auto partm = fParticleMass[iParticle];
                    auto parta = fParticleA[iParticle];
                    const char *expression = Form("(sqrt((p_lab*%d)*(p_lab*%d)+%f*%f)-%f)/%d",partz,partz,partm,partm,partm,parta);

                    project(treeParticle,nameKeoaPart,expression,selection);

                    histKeoaArray[iLR][iSys][iCutTheta][iCutYP][iParticle] = histKeoa;
                  }
                  
                  if (drawY0R21)
                  {
                    auto nameY0R21Part = makeName("y0",iAna,iLR,iMult,iSys,iCutTheta,iCutYP,iParticle);
                    auto titleY0R21Part = makeTitle("Cor.",iAna,iLR,iMult,iSys,iCutTheta,iCutYP,iParticle);

                    auto histY0R21 = new TH1D(nameY0R21Part,Form("%s;y_{0};",titleY0R21Part),bnRY0.fN,bnRY0.fMin,bnRY0.fMax);
                    TCut selection = cut0 * TCut(Form("%s/eff/%d",getProbString(iSys,iParticle),numEventsInAna[iSys])) * cutParticleYP * cutTheta  * cutParticlePoz * cutParticleSD;
                    if (fUseHandCut) selection = cut0 * TCut(Form("%s/eff/%d","1",numEventsInAna[iSys])) * getSelBoundHandCut(iSys,iParticle) * cutParticleYP * cutTheta * cutParticlePoz * cutParticleSD;

                    project(treeParticle,nameY0R21Part,"fy_cm/(by_cm/2)",selection);

                    histY0Array[iSys][iCutTheta][iCutYP][iParticle] = histY0R21;
                  }

                  if (drawNZR21)
                  {
                    auto nameYPR21Part = makeName("ypR21",iAna,iLR,iMult,iSys,iCutTheta,iCutYP,iParticle);
                    auto titleYPR21Part = makeTitle("Cor.",iAna,iLR,iMult,iSys,iCutTheta,iCutYP,iParticle);

                    auto histYPR21Part = (bnRY0*bnRPtoa).newHist(nameYPR21Part,titleYPR21Part);
                    TCut selection = cut0 * TCut(Form("%s/eff/%d",getProbString(iSys,iParticle),numEventsInAna[iSys])) * cutParticleYP * cutTheta * cutParticlePoz * cutParticleSD;
                    if (fUseHandCut) selection = cut0 * TCut(Form("%s/eff/%d","1",numEventsInAna[iSys])) * getSelBoundHandCut(iSys,iParticle) * cutParticleYP * cutTheta * cutParticlePoz * cutParticleSD;
                    project(treeParticle,nameYPR21Part,Form("pt_cm/%d:fy_cm/(by_cm/2)",fParticleA[iParticle]),selection);

                    histYPR21Array[iSys][iCutTheta][iCutYP][iParticle] = histYPR21Part;
                  }

                  delete treeParticle;
                }

                if (drawCorPID)
                {
                  auto namePIDCor = makeName("pidCor",iAna,iLR,iMult,iSys,iCutTheta,iCutYP);
                  auto cvsPIDCor = makeCvs2(namePIDCor,1000,700);
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
                    setDrawAtt(hist,iParticle);
                    if (iParticle==0) hist -> Draw("plhist");
                    else hist -> Draw("sameplhist");
                  }
                }

              }
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
                  setDrawAtt(graph,iParticle);

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
                  setDrawAtt(histTemp,iSys);

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
                  auto hist = (bnRY0,bnR21).newHist(nameR21,titleR21);
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
                    setDrawAtt(graph,iParticle);

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
                      auto hist = (bnRY0,bnPR).newHist(namePR,titlePR);
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
                        setDrawAtt(hist_s1,iSys);
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
                  setDrawAtt(histTemp,iSys);

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

                auto cvs = makeCvs(nameR21,680,550);
                //cvs -> SetGrid(1,1);

                auto hist = (bnKeoa*bnR21).newHist(nameR21,titleR21);
                hist -> Draw();

                auto legend = new TLegend();
                for (auto iParticle : fParticleIdx)
                {
                  const char *nameParticle = fParticleNames[iParticle];

                  auto hist1 = histKeoaArray[iLR][iSys2][iCutTheta][iCutYP][iParticle];
                  auto hist2 = histKeoaArray[iLR][iSys1][iCutTheta][iCutYP][iParticle];

                  const char *nameR21Part = Form("%s_%s",nameR21,nameParticle);

                  auto hist0 = (TH1D *) hist1 -> Clone(nameR21Part);
                  hist0 -> Divide(hist2);

                  setDrawAtt(hist0,iParticle);

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

                TCanvas *cvsAlpha;
                TCanvas *cvsMBeta;
                if (setDrawNZAll) {
                  cvsAlpha = makeCvs(nameAlpha);
                  cvsMBeta = makeCvs(nameMBeta);
                }
                else {
                  cvsAlpha = makeCvs(nameAlpha,100+240*bnRY0.fN,100*bnRPtoa.fN);
                  cvsMBeta = makeCvs(nameMBeta,100+240*bnRY0.fN,100*bnRPtoa.fN);
                  cvsDivideM0(cvsAlpha, bnRY0.fN, bnRPtoa.fN);
                  cvsDivideM0(cvsMBeta, bnRY0.fN, bnRPtoa.fN);
                }

                double r21Array[10][10][fNumParticles] = {{0}};
                double r21ErrorArray[10][10][fNumParticles] = {{0}};
                double alphaArray[10][10] = {{0}};
                double mbetaArray[10][10] = {{0}};
                double alphaErrorArray[10][10] = {{0}};
                double mbetaErrorArray[10][10] = {{0}};

                for (auto ibinY0=0; ibinY0<bnRY0.fN; ++ibinY0)
                {
                  for (auto ibinPtoa=0; ibinPtoa<bnRPtoa.fN; ++ibinPtoa)
                  {
                    auto binY0 = ibinY0 + 1;
                    auto binPtoa = ibinPtoa + 1;

                    auto graphAll = new TGraph2D();
                    auto graph_z1_inN = new TGraph(); graph_z1_inN -> SetMarkerSize(1.5); graph_z1_inN -> SetMarkerStyle(20); graph_z1_inN -> SetMarkerColor(kBlack);
                    auto graph_z2_inN = new TGraph(); graph_z2_inN -> SetMarkerSize(1.5); graph_z2_inN -> SetMarkerStyle(25); graph_z2_inN -> SetMarkerColor(kRed  );
                    auto graph_n0_inZ = new TGraph(); graph_n0_inZ -> SetMarkerSize(1.5); graph_n0_inZ -> SetMarkerStyle(20); graph_n0_inZ -> SetMarkerColor(kBlack);
                    auto graph_n1_inZ = new TGraph(); graph_n1_inZ -> SetMarkerSize(1.5); graph_n1_inZ -> SetMarkerStyle(25); graph_n1_inZ -> SetMarkerColor(kRed  );
                    auto graph_n2_inZ = new TGraph(); graph_n2_inZ -> SetMarkerSize(1.5); graph_n2_inZ -> SetMarkerStyle(22); graph_n2_inZ -> SetMarkerColor(kBlue );

                    for (auto iParticle : fParticleIdx)
                    {
                      auto value2 = histYPR21Array[iSys2][iCutTheta][iCutYP][iParticle] -> GetBinContent(binY0, binPtoa);
                      auto value1 = histYPR21Array[iSys1][iCutTheta][iCutYP][iParticle] -> GetBinContent(binY0, binPtoa);
                      auto r21Value = value2/value1;
                      r21Array[binY0][binPtoa][iParticle] = r21Value;

                      //auto r21Error = sqrt((1/sqrt(value1))*(1/sqrt(value1)) * (1/sqrt(value2))*(1/sqrt(value2)));
                      auto r21Error = 0.001*sqrt((1/sqrt(value1))*(1/sqrt(value1)) * (1/sqrt(value2))*(1/sqrt(value2)));
                      r21ErrorArray[binY0][binPtoa][iParticle] = r21Error;

                      auto zValue = fParticleZ[iParticle];
                      auto nValue = fParticleN[iParticle];

                      if (zValue==1) graph_z1_inN -> SetPoint(graph_z1_inN->GetN(), nValue, r21Value);
                      if (zValue==2) graph_z2_inN -> SetPoint(graph_z2_inN->GetN(), nValue, r21Value);
                      if (nValue==0) graph_n0_inZ -> SetPoint(graph_n0_inZ->GetN(), zValue, r21Value);
                      if (nValue==1) graph_n1_inZ -> SetPoint(graph_n1_inZ->GetN(), zValue, r21Value);
                      if (nValue==2) graph_n2_inZ -> SetPoint(graph_n2_inZ->GetN(), zValue, r21Value);
                      graphAll -> SetPoint(graphAll->GetN(), nValue, zValue, r21Value);
                    }

                    auto idxCvs = (bnRPtoa.fN-ibinPtoa-1)*bnRY0.fN+ibinY0+1;

                    auto cvsAlpha0 = cvsAlpha -> cd();
                    if (!setDrawNZAll)
                      cvsAlpha0 = cvsAlpha -> cd(idxCvs);
                    if (setDrawNZLogy) cvsAlpha0 -> SetLogy();
                    auto frameAlpha = new TH2D(Form("%s_%d",cvsAlpha->GetName(),idxCvs),";N;R21",100,-1,3,100,0.5,2.0);
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
                    frameAlpha -> Draw();
                    if (!setDrawNZAll && ibinY0==0&&ibinPtoa==bnRPtoa.fN-1) {}
                    else {
                      graph_z1_inN -> Draw("samep");
                      graph_z2_inN -> Draw("samep");
                    }
                    TLegend *legendAlpha = new TLegend();

                    double alphaValue = 0;
                    double betaValue = 0;

                    TF2 *fitIsocaling = new TF2("fitIsocaling","[2]*exp([0]*x+[1]*y)",1,2,0,2);
                    if (setFitNZR21TG) {
                      auto fitResultPtr = graphAll -> Fit(fitIsocaling);

                      alphaValue = fitIsocaling -> GetParameter(0);
                      betaValue = fitIsocaling -> GetParameter(1);
                      alphaArray[binY0][binPtoa] = alphaValue;
                      mbetaArray[binY0][binPtoa] = -betaValue;

                      if (setABErrorToResidual)
                      {
                        double rmsResidualR21 = 0.;
                        for (auto iParticle : fParticleIdx)
                        {
                          auto nValue = fParticleN[iParticle];
                          auto zValue = fParticleZ[iParticle];
                          auto r21Value = r21Array[binY0][binPtoa][iParticle];
                          auto residual = r21Value - fitIsocaling -> Eval(nValue, zValue);
                          rmsResidualR21 += (residual/r21Value)*(residual/r21Value);
                        }
                        rmsResidualR21 = sqrt(rmsResidualR21/fNumParticles);
                        alphaErrorArray[binY0][binPtoa] = rmsResidualR21;
                        mbetaErrorArray[binY0][binPtoa] = rmsResidualR21;
                      }
                      else {
                        double rmsR21 = 0;
                        for (auto iParticle : fParticleIdx)
                          rmsR21 =+ r21ErrorArray[binY0][binPtoa][iParticle]*r21ErrorArray[binY0][binPtoa][iParticle];
                        rmsR21 = sqrt(rmsR21);
                        alphaErrorArray[binY0][binPtoa] = rmsR21;
                        mbetaErrorArray[binY0][binPtoa] = rmsR21;
                      }

                      TF1 *fit1 = new TF1("fit1","exp([0]*x+[1])",-.5,2.5);
                      fit1 -> SetParameters(alphaValue,betaValue*1);
                      fit1 -> SetLineColor(kGray+1);

                      TF1 *fit2 = new TF1("fit1","exp([0]*x+[1])",-.5,2.5);
                      fit2 -> SetParameters(alphaValue,betaValue*2);
                      fit2 -> SetLineColor(kGray+1);

                      if (!setDrawNZAll && ibinY0==0&&ibinPtoa==bnRPtoa.fN-1) {}
                      else {
                        fit1 -> Draw("samel");
                        fit2 -> Draw("samel");
                      }

                      auto mm = new TMarker(0,0,fDrawMStyle[ibinPtoa]);
                      mm -> SetMarkerColor(fDrawColor[ibinY0]);
                      mm -> SetMarkerSize(1.8);
                      legendAlpha -> AddEntry(mm, Form("#alpha = %.3f",alphaValue),"p");
                    }

                    if (!setDrawNZAll && ibinY0==0&&ibinPtoa==bnRPtoa.fN-1)
                    {
                      auto legendGraph = new TLegend();
                      legendGraph -> AddEntry(graph_z1_inN,"Z=1","p");
                      legendGraph -> AddEntry(graph_z2_inN,"Z=2","p");
                      makeLegend(cvsAlpha0,legendGraph,"lt",.2,0,.6,.8) -> Draw();
                    }
                    else
                      makeLegend(cvsAlpha0,legendAlpha,"lt",0,0,.5,.25) -> Draw();

                    auto cvsMBeta0 = cvsMBeta -> cd();
                    if (!setDrawNZAll)
                      cvsMBeta0 = cvsMBeta -> cd(idxCvs);
                    if (setDrawNZLogy) cvsMBeta0 -> SetLogy();
                    auto frameBeta = new TH2D(Form("%s_%d",cvsMBeta->GetName(),idxCvs),";Z;R21",100,0,3,100,0.5,2.0);
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
                    frameBeta -> Draw();
                    if (!setDrawNZAll && ibinY0==0&&ibinPtoa==bnRPtoa.fN-1) {}
                    else {
                      graph_n0_inZ -> Draw("samep");
                      graph_n1_inZ -> Draw("samep");
                      graph_n2_inZ -> Draw("samep");
                    }
                    TLegend *legendBeta = new TLegend();

                    if (setFitNZR21TG) {
                      TF1 *fit1 = new TF1("fit1","exp([0]+[1]*x)",.5,2.5);
                      fit1 -> SetParameters(alphaValue*1,betaValue);
                      fit1 -> SetLineColor(kGray+1);

                      TF1 *fit2 = new TF1("fit1","exp([0]+[1]*x)",.5,2.5);
                      fit2 -> SetParameters(alphaValue*2,betaValue);
                      fit2 -> SetLineColor(kGray+1);

                      if (!setDrawNZAll && ibinY0==0&&ibinPtoa==bnRPtoa.fN-1) {}
                      else {
                        fit1 -> Draw("samel");
                        fit2 -> Draw("samel");
                      }

                      auto mm = new TMarker(0,0,fDrawMStyle[ibinPtoa]);
                      mm -> SetMarkerColor(fDrawColor[ibinY0]);
                      mm -> SetMarkerSize(1.8);
                      legendBeta -> AddEntry(mm, Form("-#beta = %.3f",-betaValue),"p");
                    }
                    if (!setDrawNZAll && ibinY0==0&&ibinPtoa==bnRPtoa.fN-1)
                    {
                      auto legendGraph = new TLegend();
                      legendGraph -> AddEntry(graph_n0_inZ,"N=0","p");
                      legendGraph -> AddEntry(graph_n1_inZ,"N=1","p");
                      legendGraph -> AddEntry(graph_n2_inZ,"N=2","p");
                      makeLegend(cvsMBeta0,legendGraph,"lt",.2,0,.6,.8) -> Draw();
                    }
                    else
                      makeLegend(cvsMBeta0,legendBeta,"rt",0,0,0.5,.25) -> Draw();
                  }
                }

                auto nameAlphaYP = makeName("alaph_yp",iAna,iLR,iMult,iSysComb,iCutTheta,iCutYP);
                auto nameMBetaYP = makeName("mbeta_yp",iAna,iLR,iMult,iSysComb,iCutTheta,iCutYP);
                auto titleAlphaYP = makeTitle("alaph_yp",iAna,iLR,iMult,iSysComb,iCutTheta,iCutYP);
                auto titleMBetaYP = makeTitle("mbeta_yp",iAna,iLR,iMult,iSysComb,iCutTheta,iCutYP);
                auto histAlpha = (bnRY0*bnRPtoa).newHist(nameAlphaYP,Form("#alpha %s",fSysCombTitles[iComb]));
                auto histMBeta = (bnRY0*bnRPtoa).newHist(nameMBetaYP,Form("-#beta %s",fSysCombTitles[iComb]));
                histAlpha -> GetZaxis() -> SetRangeUser(0,0.4);
                histMBeta -> GetZaxis() -> SetRangeUser(0,0.4);

                if (!setDrawNZAll) {
                  auto nameR21AB = makeName("r21_ab",iAna,iLR,iMult,iSysComb,iCutTheta,iCutYP);
                  auto cvsR21AB = makeCvs(nameR21AB,630,640);
                  auto titleR21AB = makeTitle("",iAna,iLR,iMult,iSysComb,iCutTheta,iCutYP);
                  auto frameR21AB = (bnAlpha*bnBeta).newHist(nameR21AB,titleR21AB);
                  if (setDrawInvB) frameR21AB -> SetTitle(Form("%s;#alpha;|#beta|",titleR21AB));
                  frameR21AB -> GetYaxis() -> SetTitleOffset(1.4);
                  frameR21AB -> Draw();

                  TGraph *graphAlpha[10];
                  TGraph *graphMBeta[10];

                  bool drawMarkers = true;

                  auto legendAB = new TLegend();
                  vector <TMarker*> markers;
                  TGraphErrors *graphABForFit = new TGraphErrors();
                  TGraphErrors *graphABY0[10] = {0};
                  for (auto ibinPtoa=bnRPtoa.fN-1; ibinPtoa>=0; --ibinPtoa) {
                    auto binPtoa = ibinPtoa + 1;
                    TGraphErrors *graphAB = new TGraphErrors();
                    graphAlpha[ibinPtoa] = new TGraph();
                    graphMBeta[ibinPtoa] = new TGraph();
                    for (auto ibinY0=bnRY0.fN-1; ibinY0>=0; --ibinY0) {
                      if (graphABY0[ibinY0]==nullptr) {
                        graphABY0[ibinY0] = new TGraphErrors();
                      }
                      auto binY0 = ibinY0 + 1;
                      auto alpha = alphaArray[binY0][binPtoa];
                      auto mbeta = mbetaArray[binY0][binPtoa];
                      auto beta = -mbeta;
                      auto alphaError = alphaErrorArray[binY0][binPtoa];
                      auto mbetaError = mbetaErrorArray[binY0][binPtoa];
                      auto betaError = mbetaError;
                      //auto r21Error = r21ErrorArray[binY0][binPtoa][0] /100;
                      TMarker *markerAB = nullptr;
                      if (drawMarkers) {
                        if (setDrawInvB) {
                          if (bnAlpha.isInside(alpha)&&bnBeta.isInside(beta)) {
                            markerAB = new TMarker(alpha,mbeta,fDrawMStyle[ibinPtoa]);
                            setDrawAtt(markerAB,ibinPtoa,ibinY0);
                            markerAB -> SetMarkerSize(1.6);
                            markers.push_back(markerAB);
                          }
                        }
                        else {
                          if (bnAlpha.isInside(alpha)&&bnBeta.isInside(beta)) {
                            markerAB = new TMarker(alpha,beta,fDrawMStyle[ibinPtoa]);
                            setDrawAtt(markerAB,ibinPtoa,ibinY0);
                            markerAB -> SetMarkerSize(1.6);
                            markers.push_back(markerAB);
                          }
                        }
                      }
                      if (setDrawInvB) {
                        graphAB -> SetPoint(graphAB->GetN(),alpha,mbeta);
                        if (setDrawABError) graphAB -> SetPointError(graphAB->GetN()-1,alphaError,mbetaError);
                        graphABY0[ibinY0] -> SetPoint(graphABY0[ibinY0]->GetN(),alpha,mbeta);
                        if (setDrawABError) graphABY0[ibinY0] -> SetPointError(graphABY0[ibinY0]->GetN()-1,alphaError,mbetaError);
                        graphABForFit -> SetPoint(graphABForFit->GetN(),alpha,mbeta);
                        if (setDrawABError) graphABForFit -> SetPointError(graphABForFit->GetN()-1,alphaError,mbetaError);
                      }
                      else {
                        graphAB -> SetPoint(graphAB->GetN(),alpha,beta);
                        if (setDrawABError) graphAB -> SetPointError(graphAB->GetN()-1,alphaError,betaError);
                        graphABY0[ibinY0] -> SetPoint(graphABY0[ibinY0]->GetN(),alpha,beta);
                        if (setDrawABError) graphABY0[ibinY0] -> SetPointError(graphABY0[ibinY0]->GetN()-1,alphaError,betaError);
                        graphABForFit -> SetPoint(graphABForFit->GetN(),alpha,beta);
                        if (setDrawABError) graphABForFit -> SetPointError(graphABForFit->GetN()-1,alphaError,betaError);
                      }
                      histAlpha -> SetBinContent(ibinY0+1,ibinPtoa+1,alpha);
                      histMBeta -> SetBinContent(ibinY0+1,ibinPtoa+1,mbeta);
                      graphAlpha[ibinPtoa] -> SetPoint(graphAlpha[ibinPtoa]->GetN(), bnRY0.getCenter(binY0), alpha);
                      graphMBeta[ibinPtoa] -> SetPoint(graphMBeta[ibinPtoa]->GetN(), bnRY0.getCenter(binY0), mbeta);
                    }
                    for (auto graph : {graphAlpha[ibinPtoa],graphMBeta[ibinPtoa]}) {
                      graph -> SetLineColor(ibinPtoa);
                      graph -> SetMarkerColor(ibinPtoa);
                      graph -> SetMarkerStyle(20);
                    }
                    if (!drawMarkers) {
                      setDrawAtt(graphAB,ibinPtoa);
                      graphAB -> SetMarkerStyle(25);
                      graphAB -> SetMarkerSize(1.6);
                      graphAB -> Draw("samep");
                      legendAB -> AddEntry(graphAB,Form("%d) p_{T}/A = %.0f ~ %.0f MeV/c",ibinPtoa, bnRPtoa.lowEdge(binPtoa),bnRPtoa.highEdge(binPtoa)));
                    }
                    else {
                      graphAB -> SetMarkerStyle(1);
                      graphAB -> SetMarkerSize(.5);
                      graphAB -> SetLineColor(kGray);
                      graphAB -> Draw("ez");
                      //graphAB -> SetFillColor(kGray);
                      //graphAB -> Draw("a2");
                    }
                  }
                  if (setFitABSep) {
                    for (auto ibinY0=bnRY0.fN-1; ibinY0>=0; --ibinY0)
                    {
                      //graphABY0[ibinY0] -> SetMarkerStyle(20);
                      //graphABY0[ibinY0] -> Draw("samep");
                      TF1 *fitSep = nullptr;
                      if (setDrawInvB) fitSep = new TF1(Form("absep_invf_%s",nameR21AB),"x+[0]",0,.4);
                      else             fitSep = new TF1(Form("absep_linf_%s",nameR21AB),"-x+[0]",0,.4);
                      graphABY0[ibinY0] -> Fit(fitSep);
                      setDrawAtt(fitSep,ibinY0);
                      fitSep -> SetLineStyle(2);
                      fitSep -> Draw("samel");
                    }
                  }
                  if (setFitAB) {
                    TF1 *fitLinear = nullptr;
                    if (setDrawInvB) fitLinear = new TF1(Form("ab_invf_%s",nameR21AB),"x+[0]",0,.4);
                    else             fitLinear = new TF1(Form("ab_linf_%s",nameR21AB),"-x+[0]",0,.4);
                    graphABForFit -> SetMarkerStyle(20);
                    //graphABForFit -> SetMarkerColor(kSpring-6);
                    graphABForFit -> SetMarkerColor(kBlack);
                    //graphABForFit -> Draw("samep");
                    graphABForFit -> Fit(fitLinear,"","RQ0");
                    //fitLinear -> SetLineColor(kRed-4);
                    //fitLinear -> SetLineColor(kSpring-6);
                    fitLinear -> SetLineColor(kSpring-1);
                    fitLinear -> SetLineStyle(2);
                    fitLinear -> Draw("samel");
                  }
                  if (drawMarkers) {
                    for (auto mab : markers)
                      mab -> Draw("p");
                    for (auto ibinPtoa=0; ibinPtoa<bnRPtoa.fN; ++ibinPtoa) {
                      auto binPtoa = ibinPtoa + 1;
                      auto markerAB = new TMarker(0,0,0);
                      setDrawAtt(markerAB,ibinPtoa,0);
                      legendAB -> AddEntry(markerAB,Form("p_{T}/A = %.0f ~ %.0f MeV/c",bnRPtoa.lowEdge(binPtoa),bnRPtoa.highEdge(binPtoa)),"p");
                    }
                    for (auto ibinY0=0; ibinY0<bnRY0.fN; ++ibinY0) {
                      auto binY0 = ibinY0 + 1;
                      auto lineAB = new TGraph();
                      setDrawAtt(lineAB,ibinY0);
                      lineAB -> SetLineWidth(2);
                      legendAB -> AddEntry(lineAB,Form("y_{0} = %.2f ~ %.2f",bnRY0.lowEdge(binY0),bnRY0.highEdge(binY0)),"l");
                    }
                  }

                  if (setDrawInvB) makeLegend(cvsR21AB,legendAB,"lt",0,0,.35,.4) -> Draw();
                  else             makeLegend(cvsR21AB,legendAB,"rt",0,0,.35,.4) -> Draw();

                  TF1 *f1RefAB = nullptr;
                  if (setDrawInvB) f1RefAB = new TF1(Form("ref_%s",nameR21AB),"x",0,.4);
                  else f1RefAB = new TF1(Form("ref_%s",nameR21AB),"-x",0,.4);
                  f1RefAB -> SetLineStyle(2);
                  f1RefAB -> SetLineColor(kGray);
                  f1RefAB -> Draw("samel");
                }

                if (!setDrawNZAll) {
                  auto nameABYP = makeName("ab_yp",iAna,iLR,iMult,iSysComb,iCutTheta,iCutYP);
                  auto cvsAB = makeCvs2(nameABYP,600,1100);
                  cvsDivide(cvsAB,1,2);
                  cvsAB -> cd(1); histAlpha -> Draw("colz");
                  cvsAB -> cd(2); histMBeta -> Draw("colz");
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
                  auto nameR21YP = makeName("r21_yp",iAna,iLR,iMult,iSysComb,iCutTheta,iCutYP);
                  const char *nameAlphaBetaOut = Form("%s/others/alpha_beta_%s.txt",fVName.Data(),nameR21YP);
                  cout << nameAlphaBetaOut << endl;
                  std::ofstream fileAB(nameAlphaBetaOut);

                  fileAB << "alpha" << endl;
                  fileAB << ","; for (auto ibinY0=0; ibinY0<bnRY0.fN; ++ibinY0) { auto binY0 = ibinY0 + 1; fileAB << ", " << binY0; } fileAB << endl;
                  fileAB << ","; for (auto ibinY0=0; ibinY0<bnRY0.fN; ++ibinY0) { auto binY0 = ibinY0 + 1; fileAB << ", y0=" << Form("%.2f",bnRY0.lowEdge(binY0)) << "~" << Form("%.2f",bnRY0.highEdge(binY0)); } fileAB << endl;
                  for (auto ibinPtoa=bnRPtoa.fN-1; ibinPtoa>=0; --ibinPtoa) {
                    auto binPtoa = ibinPtoa + 1;
                    fileAB << binPtoa << ",  pt/a=" << int(bnRPtoa.lowEdge(binPtoa)) << "~" << int(bnRPtoa.highEdge(binPtoa));
                    for (auto ibinY0=0; ibinY0<bnRY0.fN; ++ibinY0) { auto binY0 = ibinY0 + 1; fileAB << ", " << alphaArray[binY0][binPtoa]; } fileAB << endl;
                  } fileAB << endl;

                  fileAB << "-beta" << endl;
                  fileAB << ","; for (auto ibinY0=0; ibinY0<bnRY0.fN; ++ibinY0) { auto binY0 = ibinY0 + 1; fileAB << ", " << binY0; } fileAB << endl;
                  fileAB << ","; for (auto ibinY0=0; ibinY0<bnRY0.fN; ++ibinY0) { auto binY0 = ibinY0 + 1; fileAB << ", y0=" << Form("%.2f",bnRY0.lowEdge(binY0)) << "~" << Form("%.2f",bnRY0.highEdge(binY0)); } fileAB << endl;
                  for (auto ibinPtoa=bnRPtoa.fN-1; ibinPtoa>=0; --ibinPtoa) {
                    auto binPtoa = ibinPtoa + 1;
                    fileAB << binPtoa << ",  pt/a=" << int(bnRPtoa.lowEdge(binPtoa)) << "~" << int(bnRPtoa.highEdge(binPtoa));
                    for (auto ibinY0=0; ibinY0<bnRY0.fN; ++ibinY0) { auto binY0 = ibinY0 + 1; fileAB << ", " << mbetaArray[binY0][binPtoa]; } fileAB << endl;
                  } fileAB << endl;

                  for (auto iParticle : fParticleIdx)
                  {
                    fileAB << "R21_" << fParticleNames[iParticle] << endl;
                    fileAB << ","; for (auto ibinY0=0; ibinY0<bnRY0.fN; ++ibinY0) { auto binY0 = ibinY0 + 1; fileAB << ", " << binY0; } fileAB << endl;
                    fileAB << ","; for (auto ibinY0=0; ibinY0<bnRY0.fN; ++ibinY0) { auto binY0 = ibinY0 + 1; fileAB << ", y0=" << Form("%.2f",bnRY0.lowEdge(binY0)) << "~" << Form("%.2f",bnRY0.highEdge(binY0)); } fileAB << endl;
                    for (auto ibinPtoa=bnRPtoa.fN-1; ibinPtoa>=0; --ibinPtoa) {
                      auto binPtoa = ibinPtoa + 1;
                      fileAB << binPtoa << ",  pt/a=" << int(bnRPtoa.lowEdge(binPtoa)) << "~" << int(bnRPtoa.highEdge(binPtoa));
                      for (auto ibinY0=0; ibinY0<bnRY0.fN; ++ibinY0) { auto binY0 = ibinY0 + 1; fileAB << ", " << r21Array[binY0][binPtoa][iParticle]; } fileAB << endl;
                    } fileAB << endl;
                  }

                  fileAB.close();
                }
              }
            }
          }
        }


      } // iLR

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
                auto hist = new TH2D(nameR21,Form("%s;KE_{Lab}/A (MeV);dN/d(KE_{Lab}/A)/#Delta#Omega",titleR21),100,0,400,100,0.001,1);
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

                      setDrawAtt(hist1,iSys);
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

void project(TTree *tree, const char *name, const char *expr, TCut selection, bool verbose)
{
  if (verbose) cout_info << tree -> GetName() << " -> " << name << " << " << expr << " @ " << selection << endl;
  tree -> Project(name,expr,selection);
}

TCanvas *makeCvs2(const char *name, int w, int h) {
  const char *name0 = name;
  auto cvs = new TCanvas(name0,name0,cvsXOff+20*(countCvs+1), 20*(countCvs+1), w, h);
  cvs -> SetRightMargin(0.155);
  cvs -> SetLeftMargin(0.11);
  countCvs++;
  fCvsArray.push_back(cvs);
  return cvs;
}

TCanvas *makeCvs(const char *name, int w, int h) {
  const char *name0 = name;
  auto cvs = new TCanvas(name0,name0,cvsXOff+20*(countCvs+1), 20*(countCvs+1), w, h);
  cvs -> SetLeftMargin(0.11);
  countCvs++;
  fCvsArray.push_back(cvs);
  return cvs;
}

void saveAll() {
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

void setDrawAtt(TH1 *hist, int iDraw) {
  hist -> SetLineColor(fDrawColor[iDraw]);
  hist -> SetMarkerColor(fDrawColor[iDraw]);
  hist -> SetMarkerStyle(fDrawMStyle[iDraw]);
  hist -> SetMarkerSize(fDrawMSize[iDraw]);
}
void setDrawAtt(TGraph *graph, int iDraw) {
  graph -> SetLineColor(fDrawColor[iDraw]);
  graph -> SetMarkerColor(fDrawColor[iDraw]);
  graph -> SetMarkerStyle(fDrawMStyle[iDraw]);
  graph -> SetMarkerSize(fDrawMSize[iDraw]);
}
void setDrawAtt(TMarker *marker, int iDraw, int iColor) {
  if (iColor<0) iColor = iDraw;
  marker -> SetMarkerColor(fDrawColor[iColor]);
  marker -> SetMarkerStyle(fDrawMStyle[iDraw]);
  marker -> SetMarkerSize(fDrawMSize[iDraw]);
}
void setDrawAtt(TF1 *fit, int iDraw) {
  fit -> SetLineColor(fDrawColor[iDraw]);
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
    systemTitle = Form("%d/%d",fSysBeams[iSys2],fSysBeams[iSys1]);
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
