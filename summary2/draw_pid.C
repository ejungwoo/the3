#include "KBGlobal.hh"
#include "init_variables.h"
#include "binning.h"

TString fVName;
vector<TCanvas *> fCvsArray;
int cvsXOff = 1300;
int gIdxSys = 0;
int gIdxParticle = 0;
int countCvs = 0;
bool fUseHandCut = false;
double fTMargin = 0.12;
double fBMargin = 0.20;
double fLMargin = 0.19;
double fRMargin = 0.055;
double fRMargin1 = 0.055;

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
void setParticleAttributes(TH1 *hist, int iParticle);
void setParticleAttributes(TGraph *graph, int iParticle);
void setSysAttributes(TH1 *hist, int iSys);
void setSysAttributes(TGraph *graph, int iSys);
const char *makeName(const char *mainName, int iAna, int iLR, int iMult, int iSys, int iCutTheta, int iCutYP, int iPart=-1);
const char *makeTitle(const char *mainName, int iAna, int iLR, int iMult, int iSys, int iCutTheta, int iCutYP, int iPart=5);



void draw_pid()
{
  bool saveCvsPNG = false;
  bool saveCvsRoot = false;

  TString probString;
  TCut cut0 = "prob>.6"; fVName = "vProb6"; probString = "prob";
  //TCut cut0 = ""; fVName = "vAll"; probString = "prob";

  ///////////////////////////////////////////////////////////////////////////////////////////////////

  bool drawRawPID = kHold;
  bool setZoomPID = kUnSet;
  bool setFitdEdx = kSet;
  bool setDrawGuideLine = kSet;
  bool setDrawRawPIDSub = kSet;
  bool setDrawHandCut = fUseHandCut;
  bool setDrawRawPIDProjection = kUnSet;

  bool drawCorPID = kHold;
  bool setDrawCorPIDInRaw = kSet;

  bool drawCorEKE = kHold;
  bool drawKT = kHold;
  bool drawPtoa = kHold;
  bool drawYP = kHold;

  bool setDrawYPTheta0 = kSet;
  bool setDrawYPKE = kSet;
  bool setDrawYPPoz = kSet;
  bool setDrawYPText = kSet;
  bool setDrawYPGrid = kUnSet;

  bool drawPtoaR21 = kDraw;
  bool drawY0R21 = kHold;
  bool setDrawR21 = kSet;
  bool setDrawSR = kUnSet;
  bool setDrawDR = kUnSet;
  bool setDrawTemp = kUnSet;

  bool drawKeoaR21 = kHold;
  bool drawNZR21 = kHold;
  bool setFitNZR21TG = kSet;

  bool drawDistKeoa = kHold;

  // sssssssssssssssss //////////////////////////////////////////////////////////////////////////////

  int selAna = kf7;
  //int selAna = kx0;

  int selLR = kLR;
  int selMult = kMult55;
  //int selMult = kAll;
  int selSys = kAll;

  int selCutTheta = kThetaAll;
  //vector<int> selCutThetaIdx = {0};
  vector<int> selCutThetaIdx = {0,1,2,3,4};
  //vector<int> selCutThetaIdx = {1,2,3,4};

  int selCutYP = kypAll;
  //int selCutYP = kyF;
  //int selCutYP = kpF;
  //int selCutYP = kAll;
  //vector<int> selCutYPIdx[] = {kptoa0, kptoa50, kptoa100, kptoa150, kptoa200, kptoa250, kptoa300};
  vector<int> selCutYPIdx = {kypAll, ky02, ky0, ky04, ky0610, kptoa0, kptoa50, kptoa100, kptoa150, kptoa200, kptoa250, kptoa300, kptoa350, kyF, kpF};

  vector<int> selLRIdx = {0};
  //vector<int> selLRIdx = {1};
  //vector<int> selLRIdx = {1,2};
  //vector<int> selLRIdx = {0,1,2};

  vector<int> selSysIdx = {0,1,2,3};
  //vector<int> selSysIdx = {0,1};
  vector<int> selSysCombIndx = {0};

  int selSysR21 = selSys;
  int selCombR21 = kAll;

  gSystem -> mkdir(fVName);
  gSystem -> mkdir(fVName+"/rooto");
  gSystem -> mkdir(fVName+"/figures");
  gSystem -> mkdir(fVName+"/others");

  if (drawDistKeoa) {
    selCutThetaIdx.clear();
    selCutTheta = kAll;
    for (auto i : {1,2,3,4})
      selCutThetaIdx.push_back(i);
  }

  ///////////////////////////////////////////////////////////////////////////////////////////////////

  int nbinsFrame = 100;

  binning nixPIDProjX(200,800,820); // number of total projections, min x, max y
  binning bnPoz(800,0,3000);  if (setZoomPID) bnPoz.set(800,400,1400);
  binning bndEdx(800,0,2000); if (setZoomPID) bndEdx.set(800,0,600);
  
  binning bnY0(100,-1,2);
  binning bnPtoa(100,0,800);
  binning bnKeoa(100,0,400);
  binning bnKeoa2(100,0,1500);
  binning bnTheta(100,0,90);

  binning bnR21(100,0,2.0);
  binning bnRY0(8,-1,1);
  binning bnSY0(8,-0.5,1); // selected

  binning bnRPt(8,0,400);
  binning bnRPtoa(8,0,400);
  binning bnRKeoa(8,0,400);

  if (drawDistKeoa)
    bnRKeoa.fN = 50;

  if (setDrawTemp) {
    bnSY0.set(6,-.5,1.0);
    bnRPtoa.set(10,0,400);
  }

  if (setDrawSR || setDrawDR)
    bnSY0.set(21,-.8,1.3);

  binning bnPR(100,0,1.1);
  binning bnDR(100,1,2.2);
  binning bnTemp(100,5,15);

  double scaleMaxFitdEdx = .05;
  int idxAnaY01 = 4;
  int idxAnaPtoa1 = 1;
  int numDivY0 = 5;
  int numDivPtoa = 6;

  ///////////////////////////////////////////////////////////////////////////////////////////////////

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
  cout_info << "SysComb  is = " << setw(10) <<                    "everything"                              << " ;" << setw(3) << ""          << " ("; for (auto v :    selSysCombIndx) cout << v << ","; cout << ")" << endl;

  TString stopIf0;
  cout << "enter 0 to stop: "; cin >> stopIf0; cout << endl;
  if (stopIf0=="0")
    return;

  ///////////////////////////////////////////////////////////////////////////////////////////////////

  auto getProbString = [setDrawRawPIDProjection,probString](int iSys, int iParticle) {
    const char *finalProbString = probString;
    if (setDrawRawPIDProjection) {
      gIdxSys = iSys;
      gIdxParticle = iParticle;
      finalProbString = "calculate_prob(p_lab,dedx)";
    }
    return finalProbString;
  };

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

          if (drawRawPID || setDrawRawPIDProjection)
          {
            auto nameFileAll = Form("data2/%s/sys%d_%s_%s_%s.%s.ana.all.root",anaFName,sys,anaFName,lrFName,multFName,spVersion);
            cout_info << "File : " << nameFileAll << endl;

            auto treeAll = new TChain("all");
            treeAll -> Add(nameFileAll);

            auto iCutYP = 0;
            //for (auto iCutYP : fCutY0Idx)
            {
              //if (selCutYP>=0 && selCutYP!=iCutYP) continue;
              TCut cutYP = fCutYPValues[iCutYP];

              for (auto iCutTheta : selCutThetaIdx)
              {
                if (selCutTheta>=0 && selCutTheta!=iCutTheta) continue;
                TCut cutTheta = fCutThetaValues[iCutTheta];

                auto namePIDRaw = makeName("pidRaw",iAna,iLR,iMult,iSys,iCutTheta,iCutYP);
                auto titlePIDRaw = makeTitle("Raw PID",iAna,iLR,iMult,iSys,iCutTheta,iCutYP);
                if (setZoomPID) {
                  namePIDRaw = makeName("pidRawZoomed",iAna,iLR,iMult,iSys,iCutTheta,iCutYP);
                  titlePIDRaw = makeTitle("Raw PID (Zoomed)",iAna,iLR,iMult,iSys,iCutTheta,iCutYP);
                }

                auto histPID = new TH2D(namePIDRaw,Form("%s;p/Z (MeV/c);dE/dx;",titlePIDRaw),bnPoz.fN,bnPoz.fMin,bnPoz.fMax,bndEdx.fN,bndEdx.fMin,bndEdx.fMax);
                //histPID -> SetMinimum(0.5);
                //histPID -> SetMaximum(800);
                project(treeAll,namePIDRaw,"dedx:p_lab",cutYP*cutTheta);

                TCanvas *cvsPIDRaw = nullptr;
                if (drawRawPID) {
                  cvsPIDRaw = makeCvs2(namePIDRaw,1000,700);
                  //cvsPIDRaw = makeCvs2(namePIDRaw,950,850);
                  cvsPIDRaw -> SetLogz();
                  histPID -> Draw("colz");
                  if (setDrawGuideLine)
                    for (auto iParticle : fParticleIdx)
                      graphPIDMean[iSys][iParticle] -> Draw("samel");
                  if (setDrawHandCut) {
                    for (auto iParticle : fParticleIdx)
                      cutgAll[iSys][iLR][iCutTheta][iParticle] -> Draw("samel");
                  }

                }

                if (setDrawRawPIDProjection)
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
                    if (setDrawRawPIDProjection && setFitdEdx) nameProj = makeName(Form("proj1_dedx_%d",iProj),iAna,iLR,iMult,iSys,iCutTheta,iCutYP);
                    else if (setFitdEdx)                       nameProj = makeName(Form("proj2_dedx_%d",iProj),iAna,iLR,iMult,iSys,iCutTheta,iCutYP);
                    else                                       nameProj = makeName(Form("proj3_dedx_%d",iProj),iAna,iLR,iMult,iSys,iCutTheta,iCutYP);

                    auto histProj = (TH1D *) histPID -> ProjectionY(nameProj,bin1,bin2);
                    histProj -> Rebin(dbin);

                    TCanvas *cvsProj = nullptr;

                    if (nixPIDProjX.isInside(pozC))
                    {
                      TString titleProj = TString("p/Z=")+poz1+"-"+poz2+";dE/dx;";
                      histProj -> SetTitle(titleProj);

                      if (setDrawRawPIDProjection) {
                        cvsProj = makeCvs(nameProj,1000,550);
                        cvsProj -> SetLogy();
                        histProj -> Draw();
                      }

                      double scaleAmp = 1.;
                      const char *nameProjFit = makeName(Form("projFit_dedx_%d",iProj),iAna,iLR,iMult,iSys,iCutTheta,iCutYP);
                      auto f1dEdxTotal = new TF1(nameProjFit,"gaus(0)+gaus(3)+gaus(6)+gaus(9)+gaus(12)",0,bndEdx.fMax);
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
                        auto f1dEdxParticle = new TF1(nameProjFitPart,"gaus(0)",0,bndEdx.fMax);
                        f1dEdxParticle -> SetLineColor(kGray);
                        f1dEdxParticle -> SetParameters(dedxAmp,dedxMean,dedxSigma);
                        f1dEdxParticle -> SetNpx(1000);
                        if (setDrawRawPIDProjection) {
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
                        auto f1dEdxParticle = new TF1(nameProjFit,"gaus(0)",0,bndEdx.fMax);
                        f1dEdxParticle -> SetLineColor(kGray+1);
                        auto dedxAmp = f1dEdxTotal -> GetParameter(0+3*iParticle);
                        auto dedxMean = f1dEdxTotal -> GetParameter(1+3*iParticle);
                        auto dedxSigma = f1dEdxTotal -> GetParameter(2+3*iParticle);
                        f1dEdxParticle -> SetParameters(dedxAmp,dedxMean,dedxSigma);
                        f1dEdxParticle -> SetNpx(1000);
                        if (setDrawRawPIDProjection) {
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

                      if (setDrawRawPIDProjection) {
                        cvsProj -> cd();
                        f1dEdxTotal -> Draw("samel");
                      }

                      if (drawRawPID)
                      {
                        cvsPIDRaw -> cd();
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
                      cvsPIDRaw -> cd();
                      graphFitPIDMean[iSys][iParticle] -> SetLineColor(kGray+1);
                      graphFitPIDMean[iSys][iParticle] -> Draw("samel");
                    }

                    if (setDrawRawPIDProjection && setFitdEdx) {
                      makeCvs(Form("refit_%s",namePIDRaw) ,1000,550);
                      for (auto iParticle : fParticleIdx) {
                        auto graph = graphFitPIDAmp[iSys][iParticle];
                        setParticleAttributes(graph,iParticle);

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
          if (selSysR21>=0 && selSysR21!=iSys) continue;

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
            for (auto iCutYP : fCutY0Idx)
            {
              if (selCutYP>=0 && selCutYP!=iCutYP) continue;
              TCut cutYP = fCutYPValues[iCutYP];

              auto iCutTheta = 0;
              TCut cutTheta = fCutThetaValues[iCutTheta];

              TCanvas *cvsKT = nullptr;

              for (auto iParticle : fParticleIdx)
              {
                const char *nameParticle = fParticleNames[iParticle];
                const char *titleParticle = fParticleTitles[iParticle];
                auto nameFileParticle = Form("data2/%s/sys%d_%s_%s_%s.%s.ana.%s.root",anaFName,sys,anaFName,lrFName,multFName,spVersion,nameParticle);

                auto treeParticle = new TChain(nameParticle);
                treeParticle -> Add(nameFileParticle);

                if (drawKT)
                {
                  auto nameKTCorPart = makeName("ktCor",iAna,iLR,iMult,iSys,iCutTheta,iCutYP,iParticle);
                  auto titleKTCorPart = makeTitle("Cor.",iAna,iLR,iMult,iSys,iCutTheta,iCutYP,iParticle);

                  auto histKTCorPart = new TH2D(nameKTCorPart,Form("%s;KE_{Lab}/A (MeV);#theta_{lab};",titleKTCorPart),bnKeoa.fN,bnKeoa.fMin,bnKeoa.fMax,bnTheta.fN,bnTheta.fMin,bnTheta.fMax);
                  TCut selection = cut0 * TCut(Form("%s/eff/%d",getProbString(iSys,iParticle),numEventsInAna[iSys])) * cutYP * cutTheta * fParticlePozCut[iParticle];
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

                  auto histYPPart = new TH2D(nameYPPart,Form("%s;y_{0};p_{T}/A;",titleYPPart),bnY0.fN,bnY0.fMin,bnY0.fMax,bnPtoa.fN,bnPtoa.fMin,bnPtoa.fMax);
                  TCut selection = cut0 * TCut(Form("%s/eff/%d",getProbString(iSys,iParticle),numEventsInAna[iSys])) * cutYP * cutTheta * fParticlePozCut[iParticle];
                  project(treeParticle,nameYPPart,Form("pt_cm/%d:fy_cm/(by_cm/2)",fParticleA[iParticle]),selection);

                  auto cvs = makeCvs2(nameYPPart);
                  histYPPart -> Draw("colz");
                  TLine *line0 = new TLine(0,0,0,bnPtoa.fMax);
                  line0 -> SetLineStyle(9);
                  //line0 -> Draw("samel");

                  if (setDrawYPGrid) {
                    for (auto ibinY0=0; ibinY0<=numDivY0; ++ibinY0) {
                      auto binY0 = ibinY0+idxAnaY01;
                      auto line = new TLine(bnRY0.lowEdge(binY0),0,bnRY0.lowEdge(binY0),300);
                      line -> Draw("samel");
                    }

                    for (auto ibinPtoa=0; ibinPtoa<=numDivPtoa; ++ibinPtoa) {
                      auto binPtoa = ibinPtoa+idxAnaPtoa1;
                      auto line = new TLine(-.25,bnRPt.lowEdge(binPtoa),1,bnRPt.lowEdge(binPtoa));
                      line -> Draw("samel");
                    }

                    for (auto ibinY0=0; ibinY0<numDivY0; ++ibinY0) {
                      for (auto ibinPtoa=0; ibinPtoa<numDivPtoa; ++ibinPtoa) {
                        auto binY0 = ibinY0+idxAnaY01;
                        auto binPtoa = ibinPtoa+idxAnaPtoa1;
                        auto x = bnRY0.getCenter(binY0);
                        auto y = bnRPt.getCenter(binPtoa);
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

                        auto histYPPart = new TH2D(nameYPPart2,Form("%s;y_{0};p_{T}/A;",titleYPPart),bnY0.fN,bnY0.fMin,bnY0.fMax,bnPtoa.fN,bnPtoa.fMin,bnPtoa.fMax);
                        TCut cutTheta0 = Form("theta_lab>= %d*TMath::DegToRad()&&theta_lab<%d*TMath::DegToRad()",theta0-1,theta0+1); //fCutTheta0Values[iCutTheta0];
                        selection = cut0 * TCut(Form("%s/eff/%d",getProbString(iSys,iParticle),numEventsInAna[iSys])) * cutYP * cutTheta * cutTheta0 * fParticlePozCut[iParticle];
                        project(treeParticle,nameYPPart2,Form("pt_cm/%d:fy_cm/(by_cm/2)",fParticleA[iParticle]),selection,0);
                        if (iCutTheta0==2) histYPPart -> Draw("colz");
                        else histYPPart -> Draw("samecol");

                        auto nameYPPart3 = makeName(Form("fityptta%d",theta0),iAna,iLR,iMult,iSys,iCutTheta,iCutYP,iParticle);
                        TF1 *fitPYTT0 = new TF1(nameYPPart3,"pol3",bnY0.fMin,bnY0.fMax);
                        fitPYTT0 -> SetParameter(par0, par1);
                        histYPPart -> Fit(fitPYTT0,"RQ0");
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

                      auto histYPPart = new TH2D(nameYPPart2,Form("%s;y_{0};p_{T}/A;",titleYPPart),bnY0.fN,bnY0.fMin,bnY0.fMax,bnPtoa.fN,bnPtoa.fMin,bnPtoa.fMax);
                      selection = cut0 * TCut(Form("%s/eff/%d",getProbString(iSys,iParticle),numEventsInAna[iSys])) * cutYP * cutTheta * cutPoz * fParticlePozCut[iParticle];
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
                      histYPPart -> Fit(fitPoz,"Q0");
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

                      auto histYPPart = new TH2D(nameYPPart2,Form("%s;y_{0};p_{T}/A;",titleYPPart),bnY0.fN,bnY0.fMin,bnY0.fMax,bnPtoa.fN,bnPtoa.fMin,bnPtoa.fMax);
                      selection = cut0 * TCut(Form("%s/eff/%d",getProbString(iSys,iParticle),numEventsInAna[iSys])) * cutYP * cutTheta * cutKE * fParticlePozCut[iParticle];
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

                      histYPPart -> Fit(fitKE,"Q0");
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
              TCut cutYP = fCutYPValues[iCutYP];

              for (auto iCutTheta : selCutThetaIdx)
              {
                if (selCutTheta>=0 && selCutTheta!=iCutTheta) continue;
                TCut cutTheta = fCutThetaValues[iCutTheta];

                TH2D *histPIDCor = nullptr;
                TH2D *histEKECor = nullptr;
                vector<TH2D *> histEKECorArray;// = nullptr;

                for (auto iParticle : fParticleIdx)
                {
                  const char *nameParticle = fParticleNames[iParticle];
                  const char *titleParticle = fParticleTitles[iParticle];

                  auto nameFileParticle = Form("data2/%s/sys%d_%s_%s_%s.%s.ana.%s.root",anaFName,sys,anaFName,lrFName,multFName,spVersion,nameParticle);

                  auto treeParticle = new TChain(nameParticle);
                  treeParticle -> Add(nameFileParticle);

                  if (drawCorPID)
                  {
                    auto namePIDCorPart = makeName("pidCor",iAna,iLR,iMult,iSys,iCutTheta,iCutYP,iParticle);
                    auto titlePIDCorPart = makeTitle("Cor.",iAna,iLR,iMult,iSys,iCutTheta,iCutYP,iParticle);

                    auto histPIDCorPart = new TH2D(namePIDCorPart,Form("%s;p/Z (MeV/c);dE/dx;",titlePIDCorPart),bnPoz.fN,0,bnPoz.fMax,bndEdx.fN,0,bndEdx.fMax);
                    histPIDCorPart -> SetMinimum(0.5);
                    histPIDCorPart -> SetMaximum(800);
                    TCut selection = cut0 * cutYP * cutTheta * TCut(Form("%s/eff/%d",getProbString(iSys,iParticle),numEventsInAna[iSys])) * fParticlePozCut[iParticle];
                    if (setDrawCorPIDInRaw)
                      selection = cut0 * cutYP * cutTheta * TCut(Form("%s",getProbString(iSys,iParticle))) * fParticlePozCut[iParticle];
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

                    auto histEKECorPart = new TH2D(nameEKECorPart,Form("%s;KE_{Lab} (MeV);dE/dx;",titleEKECorPart),bnKeoa2.fN,bnKeoa2.fMin,bnKeoa2.fMax,bndEdx.fN,0,100);
                    TCut selection = cut0 * TCut(Form("%s/eff/%d",getProbString(iSys,iParticle),numEventsInAna[iSys])) * cutYP * cutTheta * fParticlePozCut[iParticle];
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

                    auto histPtoa = new TH1D(namePtoaPart,Form("%s;p_{T}/A (MeV/c);",titlePtoaPart),bnRPtoa.fN,bnRPtoa.fMin,bnRPtoa.fMax);
                    TCut selection = cut0 * TCut(Form("%s/eff/%d",getProbString(iSys,iParticle),numEventsInAna[iSys])) * cutYP * cutTheta  * fParticlePozCut[iParticle];
                    if (fUseHandCut) selection = cut0 * TCut(Form("%s/eff/%d","1",numEventsInAna[iSys])) * getSelBoundHandCut(iSys,iParticle) * cutYP * cutTheta * fParticlePozCut[iParticle];
                    project(treeParticle,namePtoaPart,Form("pt_cm/%d",fParticleA[iParticle]),selection);

                    histPtoaArray[iSys][iCutTheta][iCutYP][iParticle] = histPtoa;
                  }

                  if (drawDistKeoa || drawKeoaR21)
                  {
                    auto nameKeoaPart = makeName("keoa",iAna,iLR,iMult,iSys,iCutTheta,iCutYP,iParticle);
                    auto titleKeoaPart = makeTitle("Cor.",iAna,iLR,iMult,iSys,iCutTheta,iCutYP,iParticle);

                    auto histKeoa = new TH1D(nameKeoaPart,Form("%s;KE_{Lab}/A (MeV);",titleKeoaPart),bnRKeoa.fN,bnRKeoa.fMin,bnRKeoa.fMax);
                    TCut selection = cut0 * TCut(Form("%s/eff/%d",getProbString(iSys,iParticle),numEventsInAna[iSys])) * cutYP * cutTheta * fParticlePozCut[iParticle];
                    if (fUseHandCut) selection = cut0 * TCut(Form("%s/eff/%d","1",numEventsInAna[iSys])) * getSelBoundHandCut(iSys,iParticle) * cutYP * cutTheta * fParticlePozCut[iParticle];
                    auto partz = fParticleZ[iParticle];
                    auto partm = fParticleMass[iParticle];
                    auto parta = fParticleA[iParticle];
                    const char *expression = Form("(sqrt((p_lab*%d)*(p_lab*%d)+%f*%f)-%f)/%d",partz,partz,partm,partm,partm,parta);

                    TString selectionValue = selection.GetTitle();
                    selectionValue.ReplaceAll("PARTICLEA",Form("%d",fParticleA[iParticle]));
                    selection = selectionValue;

                    project(treeParticle,nameKeoaPart,expression,selection);

                    histKeoaArray[iLR][iSys][iCutTheta][iCutYP][iParticle] = histKeoa;
                  }
                  
                  if (drawY0R21)
                  {
                    auto nameY0R21Part = makeName("y0",iAna,iLR,iMult,iSys,iCutTheta,iCutYP,iParticle);
                    auto titleY0R21Part = makeTitle("Cor.",iAna,iLR,iMult,iSys,iCutTheta,iCutYP,iParticle);

                    auto histY0R21 = new TH1D(nameY0R21Part,Form("%s;y_{0};",titleY0R21Part),bnSY0.fN,bnSY0.fMin,bnSY0.fMax);
                    TCut selection = cut0 * TCut(Form("%s/eff/%d",getProbString(iSys,iParticle),numEventsInAna[iSys])) * cutYP * cutTheta  * fParticlePozCut[iParticle];
                    if (fUseHandCut) selection = cut0 * TCut(Form("%s/eff/%d","1",numEventsInAna[iSys])) * getSelBoundHandCut(iSys,iParticle) * cutYP * cutTheta * fParticlePozCut[iParticle];
                    TString selectionValue = selection.GetTitle();
                    selectionValue.ReplaceAll("PARTICLEA",Form("%d",fParticleA[iParticle]));
                    selection = selectionValue;

                    project(treeParticle,nameY0R21Part,"fy_cm/(by_cm/2)",selection);

                    histY0Array[iSys][iCutTheta][iCutYP][iParticle] = histY0R21;
                  }

                  if (drawNZR21)
                  {
                    auto nameYPR21Part = makeName("ypR21",iAna,iLR,iMult,iSys,iCutTheta,iCutYP,iParticle);
                    auto titleYPR21Part = makeTitle("Cor.",iAna,iLR,iMult,iSys,iCutTheta,iCutYP,iParticle);

                    auto histYPR21Part = new TH2D(nameYPR21Part,Form("%s;y_{0};p_{T}/A;",titleYPR21Part),bnRY0.fN,bnRY0.fMin,bnRY0.fMax, bnRPt.fN,bnRPt.fMin,bnRPt.fMax);
                    TCut selection = cut0 * TCut(Form("%s/eff/%d",getProbString(iSys,iParticle),numEventsInAna[iSys])) * cutYP * cutTheta * fParticlePozCut[iParticle];
                    if (fUseHandCut) selection = cut0 * TCut(Form("%s/eff/%d","1",numEventsInAna[iSys])) * getSelBoundHandCut(iSys,iParticle) * cutYP * cutTheta * fParticlePozCut[iParticle];
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
                  if (!setDrawCorPIDInRaw) {
                    histPIDCor -> SetMaximum(0.08);
                    histPIDCor -> SetMinimum(0.0001);
                  }
                  histPIDCor -> Draw("colz");
                  if (setDrawGuideLine)
                    for (auto iParticle : fParticleIdx)
                      graphPIDMean[iSys][iParticle] -> Draw("samel");
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
                    setParticleAttributes(hist,iParticle);
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
          for (auto iCutYP : fCutY0Idx)
          {
            if (selCutYP>=0 && selCutYP!=iCutYP) continue;

            for (auto iCutTheta : selCutThetaIdx)
            {
              if (selCutTheta>=0 && selCutTheta!=iCutTheta) continue;

              for (auto iComb : selSysCombIndx)
              {
                if (selCombR21>=0 && selCombR21!=iComb) continue;

                auto iSys2 = fSysCombIdx[iComb][0];
                auto iSys1 = fSysCombIdx[iComb][1];
                auto iSysComb = 100 + iSys2*10 + iSys1;

                const char *nameR21 = makeName("r21_ptoa",iAna,iLR,iMult,iSysComb,iCutTheta,iCutYP);
                auto titleR21 = makeTitle("Cor.",iAna,iLR,iMult,iSysComb,iCutTheta,iCutYP);

                //auto cvs = makeCvs(nameR21,680,550);
                auto cvs = makeCvs(nameR21);

                auto hist = new TH2D(nameR21,Form("%s;p_{T}/A (MeV/c);R_{21}",titleR21),nbinsFrame,bnRPtoa.fMin,bnRPtoa.fMax,bnR21.fN,bnR21.fMin,bnR21.fMax);
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
                  setParticleAttributes(graph,iParticle);

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
                auto hist = new TH2D(nameTemp,Form("%s;p_{T}/A;T",titleTemp),nbinsFrame,bnRPtoa.fMin,bnRPtoa.fMax,bnTemp.fN,bnTemp.fMin,bnTemp.fMax);
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
                  setSysAttributes(histTemp,iSys);

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

              for (auto iComb : selSysCombIndx)
              {
                if (selCombR21>=0 && selCombR21!=iComb) continue;

                auto iSys2 = fSysCombIdx[iComb][0];
                auto iSys1 = fSysCombIdx[iComb][1];
                auto iSysComb = 100 + iSys2*10 + iSys1;

                if (setDrawR21) {
                  const char *nameR21 = makeName("r21_y0All",iAna,iLR,iMult,iSysComb,iCutTheta,iCutYP);
                  auto titleR21 = makeTitle("Cor.",iAna,iLR,iMult,iSysComb,iCutTheta,iCutYP);

                  auto cvs = makeCvs(nameR21);
                  auto hist = new TH2D(nameR21,Form("%s;y_{0};R_{21}",titleR21),nbinsFrame,bnSY0.fMin,bnSY0.fMax,bnR21.fN,bnR21.fMin,bnR21.fMax);
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
                    for (auto bin=1; bin<=bnSY0.fN; ++bin) {
                      //if (hist1->GetBinContent(bin)<0.04||hist2->GetBinContent(bin)<0.04) hist0 -> SetBinContent(bin,0);
                      //else
                        if (hist1->GetBinContent(bin)==0) {}
                      else graph -> SetPoint(graph -> GetN(), hist0 -> GetXaxis() -> GetBinCenter(bin), hist0 -> GetBinContent(bin));
                    }

                    graph -> SetMarkerSize(1.3);
                    setParticleAttributes(graph,iParticle);

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
                      auto hist = new TH2D(namePR,Form("%s;y_{0};SR(%s/%s)",titlePR,fParticleNames[iParticle1],fParticleNames[iParticle2]),nbinsFrame,bnSY0.fMin,bnSY0.fMax,bnPR.fN,bnPR.fMin,bnPR.fMax);
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
                        setSysAttributes(hist_s1,iSys);
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
                        auto hist = new TH2D(nameDR,ttlDR,nbinsFrame,bnSY0.fMin,bnSY0.fMax,bnDR.fN,1.1,1.5);
                        hist -> Draw();
                      }
                      else if (iPR==kTOP) {
                        auto hist = new TH2D(nameDR,ttlDR,nbinsFrame,bnSY0.fMin,bnSY0.fMax,bnDR.fN,1.3,2.1);
                        hist -> Draw();
                      }
                      else {
                        auto hist = new TH2D(nameDR,ttlDR,nbinsFrame,bnSY0.fMin,bnSY0.fMax,bnDR.fN,bnDR.fMin,bnDR.fMax);
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
                auto hist = new TH2D(nameTemp,Form("%s;y_{0};T",titleTemp),nbinsFrame,bnSY0.fMin,bnSY0.fMax,bnTemp.fN,bnTemp.fMin,bnTemp.fMax);
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
                  setSysAttributes(histTemp,iSys);

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
          for (auto iCutYP : fCutY0Idx)
          {
            if (selCutYP>=0 && selCutYP!=iCutYP) continue;

            for (auto iCutTheta : selCutThetaIdx)
            {
              if (selCutTheta>=0 && selCutTheta!=iCutTheta) continue;

              for (auto iComb : selSysCombIndx)
              {
                if (selCombR21>=0 && selCombR21!=iComb) continue;

                auto iSys2 = fSysCombIdx[iComb][0];
                auto iSys1 = fSysCombIdx[iComb][1];
                auto iSysComb = 100 + iSys2*10 + iSys1;

                const char *nameR21 = makeName("r21_keoa",iAna,iLR,iMult,iSysComb,iCutTheta,iCutYP);
                auto titleR21 = makeTitle("Cor.",iAna,iLR,iMult,iSysComb,iCutTheta,iCutYP);

                auto cvs = makeCvs(nameR21,680,550);
                //cvs -> SetGrid(1,1);

                auto hist = new TH2D(nameR21,Form("%s;KE_{Lab}/A (MeV);R_{21}",titleR21),100,0,400,bnR21.fN,bnR21.fMin,bnR21.fMax);
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

                  setParticleAttributes(hist0,iParticle);

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
          for (auto iCutYP : fCutY0Idx)
          {
            if (selCutYP>=0 && selCutYP!=iCutYP) continue;

            for (auto iCutTheta : selCutThetaIdx)
            {
              if (selCutTheta>=0 && selCutTheta!=iCutTheta) continue;

              for (auto iComb : selSysCombIndx)
              {
                if (selCombR21>=0 && selCombR21!=iComb) continue;

                auto iSys2 = fSysCombIdx[iComb][0];
                auto iSys1 = fSysCombIdx[iComb][1];
                auto iSysComb = 100 + iSys2*10 + iSys1;

                const char *nameAlpha = makeName("alpha",iAna,iLR,iMult,iSysComb,iCutTheta,iCutYP);
                const char *nameMBeta = makeName("mbeta",iAna,iLR,iMult,iSysComb,iCutTheta,iCutYP);

                auto cvsAlpha = makeCvs(nameAlpha,1200,780);
                auto cvsMBeta = makeCvs(nameMBeta,1200,780);
                cvsDivideM0(cvsAlpha, numDivY0, numDivPtoa);
                cvsDivideM0(cvsMBeta, numDivY0, numDivPtoa);

                double r21Array[10][10][fNumParticles] = {{0}};
                double r21ErrorArray[10][10][fNumParticles] = {{0}};
                double alphaArray[10][10] = {{0}};
                double mbetaArray[10][10] = {{0}};
                double alphaErrorArray[10][10] = {{0}};
                double mbetaErrorArray[10][10] = {{0}};

                //auto binnx = binning(histYPR21Array[iSys2][iCutTheta][iCutYP][0],1);
                //auto binny = binning(histYPR21Array[iSys2][iCutTheta][iCutYP][0],2);

                for (auto ibinY0=0; ibinY0<numDivY0; ++ibinY0)
                {
                  for (auto ibinPtoa=0; ibinPtoa<numDivPtoa; ++ibinPtoa)
                  {
                    auto binY0 = ibinY0+idxAnaY01;
                    auto binPtoa = ibinPtoa+idxAnaPtoa1;

                    cout << ibinY0 << " " << ibinPtoa << " | " << binY0 << " " << binPtoa << endl;

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
                      auto r21Error = 1./sqrt(value1);
                      r21Array[binY0][binPtoa][iParticle] = r21Value;
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

                    auto idxCvs = (numDivPtoa-ibinPtoa-1)*numDivY0+ibinY0+1;

                    auto cvs = cvsAlpha -> cd(idxCvs);
                    cvs -> SetLogy();
                    auto frameAlpha = new TH2D(Form("%s_%d",cvsAlpha->GetName(),idxCvs),";N;R21",100,-1,3,100,0.5,2.0);
                    frameAlpha -> Draw();
                    frameAlpha -> GetXaxis() -> SetTitleSize(0.10);
                    frameAlpha -> GetYaxis() -> SetTitleSize(0.10);
                    frameAlpha -> GetXaxis() -> SetLabelSize(0.10);
                    frameAlpha -> GetYaxis() -> SetLabelSize(0.10);
                    frameAlpha -> GetXaxis() -> SetNdivisions(506);
                    frameAlpha -> GetYaxis() -> SetNdivisions(506);
                    frameAlpha -> GetXaxis() -> CenterTitle();
                    frameAlpha -> GetYaxis() -> CenterTitle();
                    //if (ibinY0==1&&ibinPtoa==8) {}
                    if (ibinY0==0&&ibinPtoa==numDivPtoa-1) {}
                    else {
                      graph_z1_inN -> Draw("samep");
                      graph_z2_inN -> Draw("samep");
                    }
                    TLegend *legendAlpha = new TLegend();
                    //legendAlpha -> AddEntry(graph_z2_inN, Form("(%d,%d)",binY0,binPtoa),"");
                    legendAlpha -> AddEntry(graph_z2_inN, Form("(%d,%d)",ibinY0,ibinPtoa),"");

                    double alphaValue = 0;
                    double betaValue = 0;

                    double alphaError = 0;
                    double betaError = 0;

                    TF2 *fitAll = new TF2("fitAll","[2]*exp([0]*x+[1]*y)",1,2,0,2);
                    if (setFitNZR21TG) {
                      graphAll -> Fit(fitAll);

                      alphaValue = fitAll -> GetParameter(0);
                      betaValue = fitAll -> GetParameter(1);
                      alphaArray[binY0][binPtoa] = alphaValue;
                      mbetaArray[binY0][binPtoa] = -betaValue;

                      alphaError = fitAll -> GetParError(0);
                      betaError = fitAll -> GetParError(1);
                      alphaErrorArray[binY0][binPtoa] = alphaError;
                      mbetaErrorArray[binY0][binPtoa] = betaError;

                      TF1 *fit1 = new TF1("fit1","exp([0]*x+[1])",-.5,2.5);
                      fit1 -> SetParameters(alphaValue,betaValue*1);
                      fit1 -> SetLineColor(kGray+1);

                      TF1 *fit2 = new TF1("fit1","exp([0]*x+[1])",-.5,2.5);
                      fit2 -> SetParameters(alphaValue,betaValue*2);
                      fit2 -> SetLineColor(kGray+1);

                      if (ibinY0==0&&ibinPtoa==numDivPtoa-1) {}
                      else {
                        fit1 -> Draw("samel");
                        fit2 -> Draw("samel");
                      }

                      legendAlpha -> AddEntry(fit2, Form("#alpha = %.3f",alphaValue),"");
                    }

                    if (ibinY0==0&&ibinPtoa==numDivPtoa-1)
                    {
                      auto legendGraph = new TLegend();
                      legendGraph -> AddEntry(graph_z1_inN,"Z=1","p");
                      legendGraph -> AddEntry(graph_z2_inN,"Z=2","p");
                      makeLegend(cvs,legendGraph,"lt",.2,0,.6,.8) -> Draw();
                    }
                    else
                      makeLegend(cvs,legendAlpha,"lt",-.1,0,.5,.42) -> Draw();

                    cvs = cvsMBeta -> cd(idxCvs);
                    cvs -> SetLogy();
                    auto frameBeta = new TH2D(Form("%s_%d",cvsMBeta->GetName(),idxCvs),";Z;R21",100,0,3,100,0.5,2.0);
                    frameBeta -> Draw();
                    frameBeta -> GetXaxis() -> SetTitleSize(0.10);
                    frameBeta -> GetYaxis() -> SetTitleSize(0.10);
                    frameBeta -> GetXaxis() -> SetLabelSize(0.10);
                    frameBeta -> GetYaxis() -> SetLabelSize(0.10);
                    frameBeta -> GetXaxis() -> SetNdivisions(506);
                    frameBeta -> GetYaxis() -> SetNdivisions(506);
                    frameBeta -> GetXaxis() -> CenterTitle();
                    frameBeta -> GetYaxis() -> CenterTitle();
                    if (ibinY0==0&&ibinPtoa==numDivPtoa-1) {}
                    else {
                      graph_n0_inZ -> Draw("samep");
                      graph_n1_inZ -> Draw("samep");
                      graph_n2_inZ -> Draw("samep");
                    }
                    TLegend *legendBeta = new TLegend();
                    legendBeta -> AddEntry(graph_n2_inZ , Form("(%d,%d)",ibinY0,ibinPtoa),"");

                    if (setFitNZR21TG) {
                      TF1 *fit1 = new TF1("fit1","exp([0]+[1]*x)",.5,2.5);
                      fit1 -> SetParameters(alphaValue*1,betaValue);
                      fit1 -> SetLineColor(kGray+1);

                      TF1 *fit2 = new TF1("fit1","exp([0]+[1]*x)",.5,2.5);
                      fit2 -> SetParameters(alphaValue*2,betaValue);
                      fit2 -> SetLineColor(kGray+1);

                      if (ibinY0==0&&ibinPtoa==numDivPtoa-1) {}
                      else {
                        fit1 -> Draw("samel");
                        fit2 -> Draw("samel");
                      }

                      legendBeta -> AddEntry(fit2, Form("-#beta = %.3f",-betaValue),"");
                    }
                    if (ibinY0==0&&ibinPtoa==numDivPtoa-1)
                    {
                      auto legendGraph = new TLegend();
                      legendGraph -> AddEntry(graph_n0_inZ,"N=0","p");
                      legendGraph -> AddEntry(graph_n1_inZ,"N=1","p");
                      legendGraph -> AddEntry(graph_n2_inZ,"N=2","p");
                      makeLegend(cvs,legendGraph,"lt",.2,0,.6,.8) -> Draw();
                    }
                    else
                      makeLegend(cvs,legendBeta,"rt",0,0,.5,.42) -> Draw();
                  }
                }

                auto nameR21YP = makeName("r21_yp",iAna,iLR,iMult,iSysComb,iCutTheta,iCutYP);
                const char *nameAlphaBetaOut = Form("%s/others/alpha_beta_%s.txt",fVName.Data(),nameR21YP);
                cout << nameAlphaBetaOut << endl;
                std::ofstream fileAB(nameAlphaBetaOut);

                auto nameR21AB = makeName("r21_ab",iAna,iLR,iMult,iSysComb,iCutTheta,iCutYP);
                auto cvsR21AB = makeCvs(nameR21AB);
                auto titleR21AB = makeTitle("",iAna,iLR,iMult,iSysComb,iCutTheta,iCutYP);
                (new TH2D(nameR21AB,Form("%s;#alpha;-#beta",titleR21AB),100,0,.4,100,0,.4)) -> Draw();
                auto nameAlphaYP = makeName("alaph_yp",iAna,iLR,iMult,iSysComb,iCutTheta,iCutYP);
                auto nameMBetaYP = makeName("mbeta_yp",iAna,iLR,iMult,iSysComb,iCutTheta,iCutYP);
                auto titleAlphaYP = makeTitle("alaph_yp",iAna,iLR,iMult,iSysComb,iCutTheta,iCutYP);
                auto titleMBetaYP = makeTitle("mbeta_yp",iAna,iLR,iMult,iSysComb,iCutTheta,iCutYP);
                auto histAlpha = new TH2D(nameAlphaYP,Form("#alpha %s;y0;pt/a",fSysCombTitles[iComb]),numDivY0,bnRY0.lowEdge(idxAnaY01),bnRY0.highEdge(idxAnaY01+numDivY0-1),numDivPtoa,bnRPt.lowEdge(idxAnaPtoa1),bnRPt.highEdge(idxAnaPtoa1+numDivY0-1));
                auto histMBeta = new TH2D(nameMBetaYP,Form("-#beta %s;y0;pt/a",fSysCombTitles[iComb]),numDivY0,bnRY0.lowEdge(idxAnaY01),bnRY0.highEdge(idxAnaY01+numDivY0-1),numDivPtoa,bnRPt.lowEdge(idxAnaPtoa1),bnRPt.highEdge(idxAnaPtoa1+numDivY0-1));
                //cout << nameAlphaYP << numDivY0 << " " << bnRY0.lowEdge(idxAnaY01) << " " << bnRY0.highEdge(idxAnaY01+numDivY0-1) << endl;
                histAlpha -> GetZaxis() -> SetRangeUser(0,0.4);
                histMBeta -> GetZaxis() -> SetRangeUser(0,0.4);

                TGraph *graphAlpha[10];
                TGraph *graphMBeta[10];

                auto legendAB = new TLegend();
                for (auto ibinPtoa=numDivPtoa-1; ibinPtoa>=0; --ibinPtoa) {
                  auto binPtoa = ibinPtoa+idxAnaPtoa1;
                  auto graphAB = new TGraphErrors();
                  graphAlpha[ibinPtoa] = new TGraph();
                  graphMBeta[ibinPtoa] = new TGraph();
                  for (auto ibinY0=numDivY0-1; ibinY0>=0; --ibinY0) {
                    auto binY0 = ibinY0+idxAnaY01;
                    auto alpha = alphaArray[binY0][binPtoa];
                    auto mbeta = mbetaArray[binY0][binPtoa];
                    auto alphaError = alphaErrorArray[binY0][binPtoa];
                    auto mbetaError = mbetaErrorArray[binY0][binPtoa];
                    auto r21Error = r21ErrorArray[binY0][binPtoa][0] /100;
                    graphAB -> SetPoint(graphAB->GetN(),alpha,mbeta);
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
                  graphAB -> SetLineColor(ibinPtoa);
                  graphAB -> SetMarkerColor(ibinPtoa);
                  graphAB -> SetMarkerStyle(25);
                  graphAB -> SetMarkerSize(1.6);
                  graphAB -> Draw("samep");
                  legendAB -> AddEntry(graphAB,Form("%d) p_{T}/A = %.0f ~ %.0f MeV/c",ibinPtoa, bnRPt.lowEdge(binPtoa),bnRPt.highEdge(binPtoa)));
                }
                makeLegend(cvsR21AB,legendAB,"lt",0,0,0,.3) -> Draw();
                auto f1Linear = new TF1(Form("invf_%s",nameR21AB),"x",0,.4);
                f1Linear -> SetLineStyle(2);
                f1Linear -> SetLineColor(kGray);
                f1Linear -> Draw("samel");

                auto nameABYP = makeName("ab_yp",iAna,iLR,iMult,iSysComb,iCutTheta,iCutYP);
                auto cvsAB = makeCvs2(nameABYP,600,1100);
                cvsDivide(cvsAB,1,2);
                cvsAB -> cd(1); histAlpha -> Draw("colz");
                cvsAB -> cd(2); histMBeta -> Draw("colz");


                fileAB << "alpha" << endl;
                fileAB << ","; for (auto binY0=3; binY0<=8; ++binY0) fileAB << ", " << binY0; fileAB << endl;
                fileAB << ","; for (auto binY0=3; binY0<=8; ++binY0) fileAB << ", y0=" << Form("%.2f",bnRY0.lowEdge(binY0)) << "~" << Form("%.2f",bnRY0.highEdge(binY0)); fileAB << endl;
                for (auto binPtoa=8; binPtoa>=1; --binPtoa) {
                  fileAB << binPtoa << ",  pt/a=" << int(bnRPt.lowEdge(binPtoa)) << "~" << int(bnRPt.highEdge(binPtoa));
                  for (auto binY0=3; binY0<=8; ++binY0) fileAB << ", " << alphaArray[binY0][binPtoa]; fileAB << endl;
                } fileAB << endl;

                fileAB << "-beta" << endl;
                fileAB << ","; for (auto binY0=3; binY0<=8; ++binY0) fileAB << ", " << binY0; fileAB << endl;
                fileAB << ","; for (auto binY0=3; binY0<=8; ++binY0) fileAB << ", y0=" << Form("%.2f",bnRY0.lowEdge(binY0)) << "~" << Form("%.2f",bnRY0.highEdge(binY0)); fileAB << endl;
                for (auto binPtoa=8; binPtoa>=1; --binPtoa) {
                  fileAB << binPtoa << ",  pt/a=" << int(bnRPt.lowEdge(binPtoa)) << "~" << int(bnRPt.highEdge(binPtoa));
                  for (auto binY0=3; binY0<=8; ++binY0) fileAB << ", " << mbetaArray[binY0][binPtoa]; fileAB << endl;
                } fileAB << endl;

                for (auto iParticle : fParticleIdx)
                {
                  fileAB << "R21_" << fParticleNames[iParticle] << endl;
                  fileAB << ","; for (auto binY0=3; binY0<=8; ++binY0) fileAB << ", " << binY0; fileAB << endl;
                  fileAB << ","; for (auto binY0=3; binY0<=8; ++binY0) fileAB << ", y0=" << Form("%.2f",bnRY0.lowEdge(binY0)) << "~" << Form("%.2f",bnRY0.highEdge(binY0)); fileAB << endl;
                  for (auto binPtoa=8; binPtoa>=1; --binPtoa) {
                    fileAB << binPtoa << ",  pt/a=" << int(bnRPt.lowEdge(binPtoa)) << "~" << int(bnRPt.highEdge(binPtoa));
                    for (auto binY0=3; binY0<=8; ++binY0) fileAB << ", " << r21Array[binY0][binPtoa][iParticle]; fileAB << endl;
                  } fileAB << endl;
                }

                fileAB.close();
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

          for (auto iCutYP : fCutY0Idx)
          {
            if (selCutYP>=0 && selCutYP!=iCutYP) continue;

            //for (auto iComb : selSysCombIndx)
            {
              //if (selCombR21>=0 && selCombR21!=iComb) continue;

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

                      setSysAttributes(hist1,iSys);
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

void setParticleAttributes(TH1 *hist, int iParticle) {
  hist -> SetLineColor(fParticleColor[iParticle]);
  hist -> SetMarkerColor(fParticleColor[iParticle]);
  hist -> SetMarkerStyle(fParticleMStyle[iParticle]);
  hist -> SetMarkerSize(fParticleMSize[iParticle]);
}
void setParticleAttributes(TGraph *graph, int iParticle) {
  graph -> SetLineColor(fParticleColor[iParticle]);
  graph -> SetMarkerColor(fParticleColor[iParticle]);
  graph -> SetMarkerStyle(fParticleMStyle[iParticle]);
  graph -> SetMarkerSize(fParticleMSize[iParticle]);
}
void setSysAttributes(TH1 *hist, int iSys) {
  hist -> SetLineColor(fSysColor[iSys]);
  hist -> SetMarkerColor(fSysColor[iSys]);
  hist -> SetMarkerStyle(fSysMStyle[iSys]);
  hist -> SetMarkerSize(fSysMSize[iSys]);
}
void setSysAttributes(TGraph *graph, int iSys) {
  graph -> SetLineColor(fSysColor[iSys]);
  graph -> SetMarkerColor(fSysColor[iSys]);
  graph -> SetMarkerStyle(fSysMStyle[iSys]);
  graph -> SetMarkerSize(fSysMSize[iSys]);
}

const char *makeName(const char *mainName, int iAna, int iLR, int iMult, int iSys, int iCutTheta, int iCutYP, int iPart)
{
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
    const char *name = Form("%s_%s_%s_%s_%s_%s_%s",mainName,fAnaNames[iAna],fLRNames[iLR],fMultNames[iMult],systemName,fCutThetaNames[iCutTheta],fCutYPNames[iCutYP]);
    return name;
  }

  const char *name = Form("%s_%s_%s_%s_%s_%s_%s_%s",mainName,fAnaNames[iAna],fLRNames[iLR],fMultNames[iMult],systemName,fCutThetaNames[iCutTheta],fCutYPNames[iCutYP],fParticleNames[iPart]);
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

  //const char *title = Form("%s, %s, %s, %s, %s, %s%s%s",mainName,fAnaTitles[iAna],fLRTitles[iLR],fMultTitles[iMult],systemTitle,fCutThetaTitles[iCutTheta],y0Title,partTitle);
  //const char *title = Form("%s",fCutThetaTitles[iCutTheta]); if (fUseHandCut) title = Form("Hand-Cut-PID, %s",mainName);
  //const char *title = Form("%s, %s",fMultTitles[iMult],systemTitle); if (fUseHandCut) title = Form("Hand-Cut-PID, %s",mainName);
  //const char *title = Form("%s",fCutYPTitles[iCutYP]); if (fUseHandCut) title = Form("Hand-Cut-PID, %s",mainName);
  //const char *title = Form("%s, %s, %s, %s",mainName,fMultTitles[iMult],systemTitle,fCutThetaTitles[iCutTheta]);
  //const char *title = Form("%s, %s, %s",mainName,fMultTitles[iMult],systemTitle);
  //const char *title = Form("%s, %s%s",fMultTitles[iMult],systemTitle,partTitle);
  const char *title = Form("%s, %s, %s%s%s",fLRTitles[iLR], fMultTitles[iMult],systemTitle,y0Title,partTitle);

  return title;
}
