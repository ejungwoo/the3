#include "KBGlobal.hh"
#include "init_variables.h"
#include "binning.h"

int cvsXOff = 1300;
int giSys = 0;
int giParticle = 0;
vector<TCanvas *> cvsArray;

TString vName;

double fTMargin = 0.12;
double fBMargin = 0.20;
double fLMargin = 0.19;
double fRMargin = 0.055;
double fRMargin1 = 0.055;

TGraph *graphFitPIDMean[fNumSystems][fNumParticles] = {0};
TGraph *graphFitPIDAmp[fNumSystems][fNumParticles] = {0};
TH1F *histPIDMeta[fNumSystems][fNumParticles] = {0};
TF1 *f1PIDMean[fNumSystems][fNumParticles] = {0};
TF1 *f1PIDSigma[fNumSystems][fNumParticles] = {0};
TCutG *cutgPID[fNumSystems][fNumParticles] = {0};
TGraph *graphPIDMean[fNumSystems][fNumParticles] = {0};
TGraph *graphPIDRange[fNumSystems][fNumParticles][2] = {0};

void CvsDivide(TCanvas *cvs, int nx, int ny);
void CvsDivideM0(TCanvas *cvs, int nx, int ny);
TLegend *makeLegend(TVirtualPad *cvs, TLegend *legend, TString opt = "", double x_offset=0, double y_offset=0, double width_fixed=0, double height_fixed=0);

int fParIndex1[2] = {0,2};
int fParIndex2[2] = {1,2};
struct GlobalChi2 {
  const ROOT::Math::IMultiGenFunction *fChi21;
  const ROOT::Math::IMultiGenFunction *fChi22;
  double operator() (const double *parIn) const {
    double par1[2] = {parIn[0], parIn[2]};
    double par2[2] = {parIn[1], parIn[2]};
    return (*fChi21)(par1) + (*fChi22)(par2);
  }
  GlobalChi2(ROOT::Math::IMultiGenFunction &chi21, ROOT::Math::IMultiGenFunction &chi22) : fChi21(&chi21), fChi22(&chi22) {}
};
double fPol1Function(double *x, double *par) {
  double value = TMath::Exp(par[0] + par[1]*x[0]);
  return value;
}


double calculate_prob(double p_lab, double dedx)
{
  double probParticle = 0;
  double probSum = 0;
  for (auto iParticle : fParticleIdx) {
    double amp = graphFitPIDAmp[giSys][giParticle] -> Eval(p_lab);
    double mean =  f1PIDMean[giSys][giParticle] -> Eval(p_lab);
    double sigma = f1PIDSigma[giSys][giParticle] -> Eval(p_lab);
    auto prob = amp*TMath::Gaus(dedx, mean, sigma, true);
    probSum += prob;
    if (iParticle==giParticle)
      probParticle = prob;
  }

  return (probParticle/probSum);
}

void project(TTree *tree, const char *name, const char *expr, TCut selection)
{
  cout_info << tree -> GetName() << " -> " << name << " << " << expr << " @ " << selection << endl;
  tree -> Project(name,expr,selection);
}

int countCvs = 0;

//TCanvas *makeCvs2(const char *name, int w=680, int h=550) {
TCanvas *makeCvs2(const char *name, int w=1050, int h=950) {
  const char *name0 = Form("cvs_%s",name);
  auto cvs = new TCanvas(name0,name0,cvsXOff+20*(countCvs+1), 20*(countCvs+1), w, h);
  cvs -> SetRightMargin(0.155);
  countCvs++;
  cvsArray.push_back(cvs);
  return cvs;
}

//TCanvas *makeCvs(const char *name, int w=680, int h=550) {
TCanvas *makeCvs(const char *name, int w=1050, int h=950) {
  const char *name0 = Form("cvs_%s",name);
  auto cvs = new TCanvas(name0,name0,cvsXOff+20*(countCvs+1), 20*(countCvs+1), w, h);
  //cvs -> SetRightMargin(0.155);
  countCvs++;
  cvsArray.push_back(cvs);
  return cvs;
}

void saveAll() {
  for (auto cvs : cvsArray) {
    cvs -> cd();
    cvs -> SaveAs(vName+"/figures/"+cvs->GetName()+".png"); 
  }
}

void writeAll()
{
  for (auto cvs : cvsArray) {
    cvs -> cd();
    cvs -> SaveAs(vName+"/rooto/"+cvs->GetName()+".root"); 
  }
}

bool fUseHandCut = false;

const char *makeName(const char *mainName, int iAna, int iLR, int iMult, int iSys, int iCutTTA, int iCutYP, int iPart=-1)
{
  if (fUseHandCut)
    mainName = Form("%s_HPID_",mainName);

  const char *systemName;
  if (iSys>=100) {
    iSys = iSys - 100;
    int iSys2 = int(iSys/10);
    int iSys1 = iSys - 10*iSys2;
    systemName = Form("%so%s",fSystemNames[iSys2],fSystemNames[iSys1]);
  }
  else
    systemName = Form("%s",fSystemNames[iSys]);

  if (iPart<0) {
    const char *name = Form("%s_%s_%s_%s_%s_%s_%s",mainName,fAnaNames[iAna],fLRNames[iLR],fMultNames[iMult],systemName,fCutTTANames[iCutTTA],fCutYPNames[iCutYP]);
    //const char *name = Form("%s_%s_%s_%s_%s_%s_%s",mainName,fAnaNames[iAna],fLRNames[iLR],fMultNames[iMult],systemName,fCutTTANames[iCutTTA],fCutYPTitles[iCutYP]);
    return name;
  }

  const char *name = Form("%s_%s_%s_%s_%s_%s_%s_%s",mainName,fAnaNames[iAna],fLRNames[iLR],fMultNames[iMult],systemName,fCutTTANames[iCutTTA],fCutYPNames[iCutYP],fParticleNames[iPart]);
  return name;
}

const char *makeTitle(const char *mainName, int iAna, int iLR, int iMult, int iSys, int iCutTTA, int iCutYP, int iPart=5)
{
  if (fUseHandCut)
    mainName = Form("%s, Hand-Cut-PID",mainName);

  const char *systemTitle;
  if (iSys>=100) {
    iSys = iSys - 100;
    int iSys2 = int(iSys/10);
    int iSys1 = iSys - 10*iSys2;
    systemTitle = Form("%d/%d",fSystems[iSys2],fSystems[iSys1]);
  }
  else
    systemTitle = Form("%d",fSystems[iSys]);

  const char *partTitle = "";
  if (iPart!=5)
    partTitle = Form(", %s",fParticleNames[iPart]);

  const char *y0Title = "";
  if (iCutYP!=0)
    //y0Title = Form(", %s",fCutYPNames[iCutYP]);
    y0Title = Form(", %s",fCutYPTitles[iCutYP]);

  const char *title = Form("%s, %s, %s, %s, %s, %s%s%s",mainName,fAnaTitles[iAna],fLRTitles[iLR],fMultTitles[iMult],systemTitle,fCutTTATitles[iCutTTA],y0Title,partTitle);
  //const char *title = Form("%s",fCutTTATitles[iCutTTA]); if (fUseHandCut) title = Form("Hand-Cut-PID, %s",mainName);
  //const char *title = Form("%s, %s",fMultTitles[iMult],systemTitle); if (fUseHandCut) title = Form("Hand-Cut-PID, %s",mainName);
  //const char *title = Form("%s",fCutYPTitles[iCutYP]); if (fUseHandCut) title = Form("Hand-Cut-PID, %s",mainName);

  return title;
}

void draw_pid()
{
  //int selAna = kf7;
  int selAna = kx0;
  int selLR = kall;
  //int selMult = kn55;
  int selMult = knAll;
  int selSys = kall;
  //int selSys = k132;

  const int selCutTTAIdx[] = {0};
  //const int selCutTTAIdx[] = {1,2,3,4};
  int selCutTTA = kall;

  //int selCutYP = ky02;
  int selCutYP = kya;
  //int selCutYP = kyF;
  //int selCutYP = kyH;
  int selCutPtoaIdx[] = {kptoa0, kptoa50, kptoa100, kptoa150, kptoa200, kptoa250, kptoa300, kptoa350};

  int selLRIdx[] = {0};
  //int selLRIdx[] = {1,2};
  //int selLRIdx[] = {1};

  vector<int> selSystemIdxR21 = {0,1,2,3};
  vector<int> selSysCombIndx = {0,1,2,3};
  //vector<int> selSystemIdxR21 = {0,1};
  //vector<int> selSysCombIndx = {0};

  bool saveCvsPNG = false;
  bool saveCvsRoot = false;

  TString probString;

  //TCut cut0 = ""; vName = "vAll"; probString = "prob";
  TCut cut0 = "prob>.6"; vName = "vProb6"; probString = "1";
  //TCut cut0 = "prob>.5"; vName = "vProb5"; probString = "1";
  //TCut cut0 = "prob>.6&&fy_cm/(by_cm/2)>-.25&&fy_cm/(by_cm/2)<1"; vName = "vACuts"; probString = "1";
  //TCut cut0 = "fy_cm/(by_cm/2)>0.7"; vName = "y0GTp7"; probString = "1";

  gSystem -> mkdir(vName);
  gSystem -> mkdir(vName+"/rooto");
  gSystem -> mkdir(vName+"/figures");
  gSystem -> mkdir(vName+"/others");

  ///////////////////////////////////////////////////////////////////////////////////////////////////

  bool drawRawPID = false;
  bool drawGuideLine = false;
  bool drawHandCut = fUseHandCut;

  bool anaRawPIDProjection = false;
  bool drawRawPIDProjection = false;
  bool drawRawPIDSub = false;
  int numProjections = 300;
  double pidProjRange1 = 500;
  double pidProjRange2 = 520;
  bool fitdEdx = false;
  double scaleMaxFitdEdx = .05;

  bool drawCorPID = false;
  bool drawCorPID22222 = true;
  bool drawCorPIDOXParticle = false;

  bool drawCorEKE = false;
  bool drawYP = false;
  bool drawYPTTA0 = true;
  bool drawYPKE = true;
  bool drawYPPoz = true;
  bool drawYPText = true;
  bool drawYPGrid = false;
  bool drawKT = false;
  bool drawPtoa = false;

  bool drawPtoaR21 = false;
  bool drawY0R21 = false;
  bool drawKeoaR21 = false;
  bool drawNZR21 = true;
  bool fitNZR21 = true;
  bool fitNZR21TG = true;

  int nbinsPID = 800;
  double dEdxMin = 0;
  double dEdxMax = 600;
  double pozMin = 0;
  double pozMax = 3000;
  //double dEdxMin = 0;
  //double dEdxMax = 600;
  //double pozMin = 400;
  //double pozMax = 1400;
  int nbinsY0 = 200;
  double y0Min = -1;
  double y0Max = 2;
  int nbinsPt = 200;
  double maxPt = 800;
  int nbinsKeoa = 100;
  double keoaMax = 400.;
  int nbinsTheta = 100;
  double thetaMax = 90;
  int nbinsPtoa = 8;
  double ptoaMax = 400;
  int nbinsKeoaR21 = 4;
  int nbinsNZYP = 8;
  int nbinsY0R21 = 8;

  ///////////////////////////////////////////////////////////////////////////////////////////////////

  //const char *finalProbString = "prob";
  const char *finalProbString = "calculate_prob(p_lab,dedx)";
  if (!anaRawPIDProjection)
    finalProbString = probString;

  int selSysR21 = selSys;
  int selCombR21 = kall;
  //int selCombR21 = 3;
  /*
  if (drawPtoaR21) {
    selCombR21 = 0;
    selSystemIdxR21.push_back(fSysCombIdx[selCombR21][0]);
    selSystemIdxR21.push_back(fSysCombIdx[selCombR21][1]);
    selSysR21 = kall;
  }
  else
    for (auto sys : {0,1})
      selSystemIdxR21.push_back(sys);
      */

  if (drawYPTTA0)
  {
    nbinsY0 = 200;
    nbinsPt = 200;
  }

  ///////////////////////////////////////////////////////////////////////////////////////////////////

  TCutG *cutgAll[fNumSystems][fNumLRs][fNumCutTTAs][fNumParticles];

  auto GetCutG = [&cutgAll](int iSys, int iLR, int iCutTTA, const char *fileName) {
    auto file = new TFile(fileName);
    for (auto iParticle : {0,1,2,3,4}) {
      auto bound = (TCutG *) file -> Get(Form("bound%d",iParticle));
      bound -> SetName(Form("bound_%d_%d_%d_%d",iSys,iLR,iCutTTA,iParticle));
      bound -> SetVarX("p_lab");
      bound -> SetVarY("dedx");
      cutgAll[iSys][iLR][iCutTTA][iParticle] = bound;
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


  ///////////////////////////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////////////////////////

  gStyle -> SetOptStat(0);

  //for (auto iSys : fSystemIdx)
  for (auto iSys : selSystemIdxR21)
  {
    auto sys = fSystems[iSys];
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
      //graphPIDMean[iSys][iParticle] -> SetLineStyle(2);
      //graphPIDMean[iSys][iParticle] -> SetLineColor(kGray+2);
      graphPIDMean[iSys][iParticle] -> SetLineColor(kRed-4);
      graphPIDRange[iSys][iParticle][0] = graph2;
      graphPIDRange[iSys][iParticle][1] = graph3;
    }
  }

  for (auto iAna : fAnaIdx)
  {
    if (selAna>=0 && selAna!=iAna) continue;
    const char *anaFName = fAnaFNames[iAna];
    const char *spVersion = fAnaVersion[iAna];

    for (auto iLR : selLRIdx)
    {
      if (selLR>=0 && selLR!=iLR) continue;
      const char *lrFName = fLRFNames[iLR];

      for (auto iMult : fMultIdx)
      {
        if (selMult>=0 && selMult!=iMult) continue;
        const char *multFName = fMultFNames[iMult];

        int numEventsInAna[fNumSystems] = {0};
        TH1D *histY0Array[fNumSystems][fNumCutTTAs][fNumCutYPs][fNumParticles] = {0};
        TH1D *histPtoaArray[fNumSystems][fNumCutTTAs][fNumCutYPs][fNumParticles] = {0};
        TH1D *histKeoaArray[fNumSystems][fNumCutTTAs][fNumCutYPs][fNumParticles] = {0};
        TH2D *histYPR21Array[fNumSystems][fNumCutTTAs][fNumCutYPs][fNumParticles] = {0};

        //for (auto iSys : fSystemIdx)
        for (auto iSys : selSystemIdxR21)
        {
          if (selSys>=0 && selSys!=iSys) continue;
          auto sys = fSystems[iSys];
          const char *sysTitle = fSystemTitles[iSys];

          if (drawRawPID || anaRawPIDProjection)
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

              for (auto iCutTTA : selCutTTAIdx)
              {
                if (selCutTTA>=0 && selCutTTA!=iCutTTA) continue;
                TCut cutTTA = fCutTTAValues[iCutTTA];

                auto namePIDRaw = makeName("pidRaw",iAna,iLR,iMult,iSys,iCutTTA,iCutYP);
                auto titlePIDRaw = makeTitle("Raw",iAna,iLR,iMult,iSys,iCutTTA,iCutYP);

                auto histPID = new TH2D(namePIDRaw,Form("%s;p/Z (MeV/c);dE/dx;",titlePIDRaw),nbinsPID,pozMin,pozMax,nbinsPID,dEdxMin,dEdxMax);
                //histPID -> SetMinimum(0.5);
                //histPID -> SetMaximum(800);
                project(treeAll,namePIDRaw,"dedx:p_lab",cutYP*cutTTA);

                TCanvas *cvsPIDRaw = nullptr;
                if (drawRawPID) {
                  cvsPIDRaw = makeCvs2(namePIDRaw,1000,700);
                  //cvsPIDRaw = makeCvs2(namePIDRaw,950,850);
                  cvsPIDRaw -> SetLogz();
                  histPID -> Draw("colz");
                  if (drawGuideLine)
                    for (auto iParticle : fParticleIdx)
                      graphPIDMean[iSys][iParticle] -> Draw("samel");
                  if (drawHandCut) {
                    for (auto iParticle : fParticleIdx)
                      cutgAll[iSys][iLR][iCutTTA][iParticle] -> Draw("samel");
                  }

                }

                if (anaRawPIDProjection)
                {
                  auto binn = binning(histPID);
                  int dbin = nbinsPID/numProjections;

                  for (auto iParticle : fParticleIdx) {
                    graphFitPIDMean[iSys][iParticle] = new TGraph();
                    graphFitPIDAmp[iSys][iParticle] = new TGraph();
                  }

                  for (auto iProj=0; iProj<numProjections; ++iProj)
                  {
                    auto bin1 = (iProj)*dbin+1;
                    auto bin2 = (iProj+1)*dbin;
                    auto poz1 = binn.lowEdge(bin1);
                    auto poz2 = binn.highEdge(bin2);
                    auto pozC = (poz1 + poz2)/2.;

                    auto nameProj = makeName(Form("proj0_dedx_%d",iProj),iAna,iLR,iMult,iSys,iCutTTA,iCutYP);
                    if (drawRawPIDSub && fitdEdx) nameProj = makeName(Form("proj1_dedx_%d",iProj),iAna,iLR,iMult,iSys,iCutTTA,iCutYP);
                    else if (fitdEdx)             nameProj = makeName(Form("proj2_dedx_%d",iProj),iAna,iLR,iMult,iSys,iCutTTA,iCutYP);
                    else                          nameProj = makeName(Form("proj3_dedx_%d",iProj),iAna,iLR,iMult,iSys,iCutTTA,iCutYP);

                    auto histProj = (TH1D *) histPID -> ProjectionY(nameProj,bin1,bin2);
                    histProj -> Rebin(dbin);

                    TCanvas *cvsProj = nullptr;

                    if (pozC>pidProjRange1 && pozC<pidProjRange2)
                    {
                      //cout_info << countCvs << " " << poz1 << " " << poz2 << endl;
                      TString titleProj = TString("p/Z=")+poz1+"-"+poz2+";dE/dx;";
                      histProj -> SetTitle(titleProj);

                      if (drawRawPIDProjection) {
                        cvsProj = makeCvs(nameProj,1000,550);
                        cvsProj -> SetLogy();
                        histProj -> Draw();
                      }

                      double scaleAmp = 1.;
                      const char *nameProjFit = makeName(Form("projFit_dedx_%d",iProj),iAna,iLR,iMult,iSys,iCutTTA,iCutYP);
                      auto f1dEdxTotal = new TF1(nameProjFit,"gaus(0)+gaus(3)+gaus(6)+gaus(9)+gaus(12)",0,dEdxMax);
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

                        const char *nameProjFitPart = makeName(Form("projFitdEdx0_%d",iProj),iAna,iLR,iMult,iSys,iCutTTA,iCutYP,iParticle);
                        auto f1dEdxParticle = new TF1(nameProjFitPart,"gaus(0)",0,dEdxMax);
                        f1dEdxParticle -> SetLineColor(kGray);
                        f1dEdxParticle -> SetParameters(dedxAmp,dedxMean,dedxSigma);
                        f1dEdxParticle -> SetNpx(1000);
                        if (drawRawPIDProjection && drawRawPIDSub) {
                          cvsProj -> cd();
                          //f1dEdxParticle -> Draw("samel");
                        }


                        f1dEdxTotal -> SetParameter(0+3*iParticle,dedxAmp);
                        f1dEdxTotal -> SetParameter(1+3*iParticle,dedxMean);
                        f1dEdxTotal -> SetParameter(2+3*iParticle,dedxSigma);

                        if (iParticle==0&&pozC>1500) f1dEdxTotal -> FixParameter(0+3*iParticle,0);
                        if (iParticle==1&&pozC>2200) f1dEdxTotal -> FixParameter(0+3*iParticle,0);
                        if (iParticle==3&&pozC>1300) f1dEdxTotal -> FixParameter(0+3*iParticle,0);
                        if (iParticle==4&&pozC>1800) f1dEdxTotal -> FixParameter(0+3*iParticle,0);

                        if (iParticle==3&&pozC< 300) f1dEdxTotal -> FixParameter(0+3*iParticle,0);

                        if (fitdEdx) {
                          if (scaleMaxFitdEdx>0) {
                            //f1dEdxTotal -> SetParLimits(0+3*iParticle,dedxAmp  -scaleMaxFitdEdx*dedxAmp  ,dedxAmp  +scaleMaxFitdEdx*dedxAmp  );
                            f1dEdxTotal -> SetParLimits(1+3*iParticle,dedxMean -scaleMaxFitdEdx*dedxMean ,dedxMean +scaleMaxFitdEdx*dedxMean );
                            f1dEdxTotal -> SetParLimits(2+3*iParticle,dedxSigma-scaleMaxFitdEdx*dedxSigma,dedxSigma+scaleMaxFitdEdx*dedxSigma);
                          }
                          else {
                            f1dEdxTotal -> FixParameter(1+3*iParticle,dedxMean);
                            f1dEdxTotal -> FixParameter(2+3*iParticle,dedxSigma);
                          }
                        }
                      }

                      if (fitdEdx)
                        histProj -> Fit(f1dEdxTotal,"RQ0");

                      for (auto iParticle : fParticleIdx) {
                        const char *nameProjFitPart = makeName(Form("projFitdEdx_%d",iProj),iAna,iLR,iMult,iSys,iCutTTA,iCutYP,iParticle);
                        auto f1dEdxParticle = new TF1(nameProjFit,"gaus(0)",0,dEdxMax);
                        f1dEdxParticle -> SetLineColor(kGray+1);
                        //f1dEdxParticle -> SetLineColor(kBlue);
                        auto dedxAmp = f1dEdxTotal -> GetParameter(0+3*iParticle);
                        auto dedxMean = f1dEdxTotal -> GetParameter(1+3*iParticle);
                        auto dedxSigma = f1dEdxTotal -> GetParameter(2+3*iParticle);
                        f1dEdxParticle -> SetParameters(dedxAmp,dedxMean,dedxSigma);
                        f1dEdxParticle -> SetNpx(1000);
                        if (drawRawPIDProjection && drawRawPIDSub) {
                          cvsProj -> cd();
                          f1dEdxParticle -> Draw("samel");
                        }

                        graphFitPIDMean[iSys][iParticle] -> SetPoint(graphFitPIDMean[iSys][iParticle]->GetN(),pozC,dedxMean);
                        if (dedxAmp<0) dedxAmp = 0;
                        auto graphAmp = graphFitPIDAmp[iSys][iParticle];
                        //if (iParticle==2) cout_info << sys << " " << graphAmp->GetN() << " " << pozC << " " << dedxAmp << endl;
                        if (!isinf(dedxAmp)) {
                          graphAmp -> SetPoint(graphAmp->GetN(),pozC,dedxAmp);
                        }
                      }

                      if (drawRawPIDProjection && drawRawPIDSub) {
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

                  if (fitdEdx) {
                    for (auto iParticle : fParticleIdx) {
                      cvsPIDRaw -> cd();
                      graphFitPIDMean[iSys][iParticle] -> SetLineColor(kGray+1);
                      graphFitPIDMean[iSys][iParticle] -> Draw("samel");
                    }

                    makeCvs(Form("refit_%s",namePIDRaw) ,1000,550);
                    for (auto iParticle : fParticleIdx) {
                      auto graph = graphFitPIDAmp[iSys][iParticle];
                      if (iParticle==0) { graph -> SetMarkerStyle(24); graph -> SetMarkerColor(kRed     ); graph -> SetLineColor(kRed     ); }
                      if (iParticle==1) { graph -> SetMarkerStyle(25); graph -> SetMarkerColor(kBlue    ); graph -> SetLineColor(kBlue    ); }
                      if (iParticle==2) { graph -> SetMarkerStyle(26); graph -> SetMarkerColor(kSpring-6); graph -> SetLineColor(kSpring-6); }
                      if (iParticle==3) { graph -> SetMarkerStyle(30); graph -> SetMarkerColor(kOrange-3); graph -> SetLineColor(kOrange-3); }
                      if (iParticle==4) { graph -> SetMarkerStyle(28); graph -> SetMarkerColor(kViolet-5); graph -> SetLineColor(kViolet-5); }

                      if (iParticle==0)
                        graph -> Draw("apl");
                      else
                        graph -> Draw("samepl");
                    }
                  }
                }

              }
            }
          }
        }

        for (auto iSys : selSystemIdxR21)
        {
          if (selSysR21>=0 && selSysR21!=iSys) continue;
          auto sys = fSystems[iSys];
          const char *sysTitle = fSystemTitles[iSys];

          auto nameFileMult = Form("data2/%s/sys%d_%s_%s_%s.%s.ana.mult.root",anaFName,sys,anaFName,lrFName,multFName,spVersion);
          //cout_info << "File : " << nameFileMult << endl;

          auto treeMult = new TChain("mult");
          treeMult -> Add(nameFileMult);

          numEventsInAna[iSys] = treeMult -> GetEntries("1");
          cout_info << "Multiplicity in " << sys << " : " << numEventsInAna[iSys] << endl;

          if (drawKT || drawYP)
          {
            for (auto iCutYP : fCutY0Idx)
            {
              if (selCutYP>=0 && selCutYP!=iCutYP) continue;
              TCut cutYP = fCutYPValues[iCutYP];

              auto iCutTTA = 0;
              TCut cutTTA = fCutTTAValues[iCutTTA];

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
                  auto nameKTCorPart = makeName("ktCor",iAna,iLR,iMult,iSys,iCutTTA,iCutYP,iParticle);
                  auto titleKTCorPart = makeTitle("Cor.",iAna,iLR,iMult,iSys,iCutTTA,iCutYP,iParticle);

                  auto histKTCorPart = new TH2D(nameKTCorPart,Form("%s;KE_{Lab}/A (MeV);#theta_{lab};",titleKTCorPart),nbinsKeoa,0,keoaMax,nbinsTheta,0,thetaMax);
                  giSys = iSys;
                  giParticle = iParticle;
                  TCut selection = cut0 * TCut(Form("%s/eff/%d",finalProbString,numEventsInAna[iSys])) * cutYP * cutTTA * fParticlePozCut[iParticle];
                  auto partz = fParticleZ[iParticle];
                  auto partm = fParticleMass[iParticle];
                  auto parta = fParticleA[iParticle];
                  const char *expression = Form("theta_lab*TMath::RadToDeg():(sqrt((p_lab*%d)*(p_lab*%d)+%f*%f)-%f)/%d",partz,partz,partm,partm,partm,parta);
                  project(treeParticle,nameKTCorPart,expression,selection);

                  if (cvsKT==nullptr) {
                    cvsKT = makeCvs2(nameKTCorPart,1000,700);
                    CvsDivide(cvsKT,3,2);
                  }

                  cvsKT -> cd(iParticle+1);
                  histKTCorPart -> Draw("colz");
                }

                if (drawYP)
                {
                  auto nameYPPart = makeName("yp",iAna,iLR,iMult,iSys,iCutTTA,iCutYP,iParticle);
                  auto titleYPPart = makeTitle("Cor.",iAna,iLR,iMult,iSys,iCutTTA,iCutYP,iParticle);

                  //auto histYPPart = new TH2D(nameYPPart,Form("%s;y_{0};p_{T}/A;",titleYPPart),nbinsY0,y0Min,y0Max,nbinsPt,0,maxPt);
                  auto histYPPart = new TH2D(nameYPPart,Form("%s;y_{0};p_{T}/A;",titleYPPart),100,y0Min,y0Max,100,0,maxPt);
                  giSys = iSys;
                  giParticle = iParticle;
                  TCut selection = cut0 * TCut(Form("%s/eff/%d",finalProbString,numEventsInAna[iSys])) * cutYP * cutTTA * fParticlePozCut[iParticle];
                  project(treeParticle,nameYPPart,Form("pt_cm/%d:fy_cm/(by_cm/2)",fParticleA[iParticle]),selection);

                  auto cvs = makeCvs2(nameYPPart);
                  histYPPart -> Draw("colz");
                  TLine *line0 = new TLine(0,0,0,maxPt);
                  line0 -> SetLineStyle(9);
                  //line0 -> Draw("samel");

                  if (drawYPGrid) {
                    auto hist0 = new TH2D(Form("%s_frame",nameYPPart),"",nbinsNZYP,-1,1,nbinsNZYP,0,400);
                    auto binnx = binning(hist0,1);
                    auto binny = binning(hist0,2);
                    //for (auto ibinX=1; ibinX<=7; ++ibinX) {
                    for (auto ibinX=2; ibinX<=7; ++ibinX) {
                      auto binX = ibinX+2;
                      auto line = new TLine(binnx.lowEdge(binX),0,binnx.lowEdge(binX),300);
                      line -> Draw("samel");
                    }
                    //for (auto ibinY=1; ibinY<=9; ++ibinY) {
                    for (auto ibinY=1; ibinY<=7; ++ibinY) {
                      auto binY = ibinY;
                      auto line = new TLine(-.25,binny.lowEdge(binY),1,binny.lowEdge(binY));
                      line -> Draw("samel");
                    }

                    //for (auto ibinX=1; ibinX<=6; ++ibinX) {
                    for (auto ibinX=2; ibinX<=6; ++ibinX) {
                      //for (auto ibinY=1; ibinY<=8; ++ibinY) {
                      for (auto ibinY=1; ibinY<=6; ++ibinY) {
                        auto binX = ibinX+2;
                        auto binY = ibinY;
                        auto x = binnx.getCenter(binX);
                        auto y = binny.getCenter(binY);
                        auto text = new TLatex(x,y,Form("%d,%d",binX,binY));
                        text -> Draw();
                        text -> SetTextAlign(22);
                        text -> SetTextSize(0.03);
                      }
                    }
                  }

                  if (drawYPTTA0) {
                    double par0 = 3, par1 = -20;
                    //for (auto iCutTTA0 : {2,4,6})
                    for (auto iCutTTA0 : {4,6})
                    {
                      int theta0 = 10*iCutTTA0;
                      auto nameYPPart2 = makeName(Form("yptta%d",theta0),iAna,iLR,iMult,iSys,iCutTTA,iCutYP,iParticle);

                      auto histYPPart = new TH2D(nameYPPart2,Form("%s;y_{0};p_{T}/A;",titleYPPart),nbinsY0,y0Min,y0Max,nbinsPt,0,maxPt);
                      giSys = iSys;
                      giParticle = iParticle;
                      selection = cut0 * TCut(Form("%s/eff/%d",finalProbString,numEventsInAna[iSys])) * cutYP * cutTTA * fCutTTA0Values[iCutTTA0] * fParticlePozCut[iParticle];
                      project(treeParticle,nameYPPart2,Form("pt_cm/%d:fy_cm/(by_cm/2)",fParticleA[iParticle]),selection);

                      auto nameYPPart3 = makeName(Form("fityptta%d",theta0),iAna,iLR,iMult,iSys,iCutTTA,iCutYP,iParticle);
                      TF1 *fitPYTT0 = new TF1(nameYPPart3,"[0]*(x+1)*(x-[1])",y0Min,y0Max);
                      fitPYTT0 -> SetParameter(par0, par1);
                      fitPYTT0 -> SetParLimits(0,0,1500);
                      fitPYTT0 -> SetParLimits(1,-30,0);
                      histYPPart -> Fit(fitPYTT0,"RQ0");
                      par0 = fitPYTT0 -> GetParameter(0);
                      par1 = fitPYTT0 -> GetParameter(1);
                      fitPYTT0 -> SetLineStyle(2);
                      //fitPYTT0 -> SetLineColor(kRed-3);
                      fitPYTT0 -> SetLineColor(kSpring-6);
                      fitPYTT0 -> Draw("samel");

                      if (drawYPText) {
                        auto xAt = y0Max;
                        auto yAt = fitPYTT0 -> Eval(xAt);

                        if (yAt<maxPt) {
                          TLatex *text = new TLatex(xAt,yAt,Form("#theta_{Lab}=%d#circ",theta0));
                          text -> SetTextColor(kSpring-6);
                          text -> SetTextFont(132);
                          text -> SetTextAlign(31);
                          text -> SetTextSize(.03);
                          text -> Draw("samel");
                        }
                        else {
                          yAt = maxPt;
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

                  if (drawYPPoz)
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

                      auto nameYPPart2 = makeName(Form("yppoz%d",iCutPoz),iAna,iLR,iMult,iSys,iCutTTA,iCutYP,iParticle);

                      //auto histYPPart = new TH2D(nameYPPart2,Form("%s;y_{0};p_{T}/A;",titleYPPart),nbinsY0,y0Min,y0Max,nbinsPoz,0,maxPoz);
                      auto histYPPart = new TH2D(nameYPPart2,Form("%s;y_{0};p_{T}/A;",titleYPPart),nbinsY0,y0Min,y0Max,nbinsPt,0,maxPt);
                      giSys = iSys;
                      giParticle = iParticle;
                      selection = cut0 * TCut(Form("%s/eff/%d",finalProbString,numEventsInAna[iSys])) * cutYP * cutTTA * cutPoz * fParticlePozCut[iParticle];
                      project(treeParticle,nameYPPart2,Form("pt_cm/%d:fy_cm/(by_cm/2)",fParticleA[iParticle]),selection);

                      //if (iCutPoz==2) histYPPart -> Draw("colz");
                      //else histYPPart -> Draw("same col");

                      auto nameYPPart3 = makeName(Form("fitypke"),iAna,iLR,iMult,iSys,iCutTTA,iCutYP,iParticle);
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


                      if (drawYPText) {
                        auto xAt = fitPoz->GetX(0);
                        //TLatex *text = new TLatex(xAt,10,Form("p_{Lab}/Z=%d",int(iPoz*100)));
                        TLatex *text = new TLatex(xAt,10,Form("p/Z=%d",int(iPoz*100)));
                        text -> SetTextFont(132);
                        text -> SetTextAngle(90);
                        text -> SetTextAlign(13);
                        //text -> SetTextAlign(22);
                        text -> SetTextColor(kRed+1);
                        text -> SetTextSize(.03);
                        text -> Draw("samel");
                      }
                    }
                  }

                  if (drawYPKE)
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

                      auto nameYPPart2 = makeName(Form("ypke%d",iCutKE),iAna,iLR,iMult,iSys,iCutTTA,iCutYP,iParticle);

                      auto histYPPart = new TH2D(nameYPPart2,Form("%s;y_{0};p_{T}/A;",titleYPPart),nbinsY0,y0Min,y0Max,nbinsPt,0,maxPt);
                      giSys = iSys;
                      giParticle = iParticle;
                      selection = cut0 * TCut(Form("%s/eff/%d",finalProbString,numEventsInAna[iSys])) * cutYP * cutTTA * cutKE * fParticlePozCut[iParticle];
                      project(treeParticle,nameYPPart2,Form("pt_cm/%d:fy_cm/(by_cm/2)",fParticleA[iParticle]),selection);

                      //if (iCutKE==0) histYPPart -> Draw("colz");
                      //else histYPPart -> Draw("same col");

                      auto nameYPPart3 = makeName(Form("fitypke"),iAna,iLR,iMult,iSys,iCutTTA,iCutYP,iParticle);
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

                      if (drawYPText) {
                        auto xAt = fitKE->GetX(0);
                        TLatex *text = new TLatex(xAt,10,Form("KE=%d",int(iKE*100)));
                        text -> SetTextColor(kBlack);
                        text -> SetTextFont(132);
                        text -> SetTextAngle(90);
                        text -> SetTextAlign(13);
                        if (iParticle==4) text -> SetTextAlign(11);
                        //text -> SetTextAlign(22);
                        text -> SetTextSize(.03);
                        text -> Draw("samel");
                      }
                    }
                  }

                }
              }
            }
          }

          if (drawY0R21)
          {
            for (auto iCutYP : selCutPtoaIdx)
            {
              //if (selCutYP>=0 && selCutYP!=iCutYP) continue;
              TCut cutYP = fCutYPValues[iCutYP];

              for (auto iCutTTA : selCutTTAIdx)
              {
                TCut cutTTA = fCutTTAValues[iCutTTA];

                for (auto iParticle : fParticleIdx)
                {
                  const char *nameParticle = fParticleNames[iParticle];
                  const char *titleParticle = fParticleTitles[iParticle];

                  auto nameFileParticle = Form("data2/%s/sys%d_%s_%s_%s.%s.ana.%s.root",anaFName,sys,anaFName,lrFName,multFName,spVersion,nameParticle);
                  //cout_info << "File : " << nameFileParticle << endl;

                  auto treeParticle = new TChain(nameParticle);
                  treeParticle -> Add(nameFileParticle);

                  if (drawY0R21)
                  {
                    auto nameY0Part = makeName("y0",iAna,iLR,iMult,iSys,iCutTTA,iCutYP,iParticle);
                    auto titleY0Part = makeTitle("Cor.",iAna,iLR,iMult,iSys,iCutTTA,iCutYP,iParticle);

                    //auto histY0 = new TH1D(nameY0Part,Form("%s;y_{0};",titleY0Part),nbinsY0R21,0,1);
                    auto histY0 = new TH1D(nameY0Part,Form("%s;y_{0};",titleY0Part),6,-.5,1);
                    giSys = iSys;
                    giParticle = iParticle;
                    TCut selection = cut0 * TCut(Form("%s/eff/%d",finalProbString,numEventsInAna[iSys])) * cutYP * cutTTA  * fParticlePozCut[iParticle];

                    if (fUseHandCut) {
                      const char *nameBound = Form("bound_%d_%d_%d_%d",iSys,iLR,iCutTTA,iParticle);
                      auto selectionBound = nameBound;

                      if (iLR==klr)
                      {
                        const char *conditionBoundL = "(-0.8<phi_cm&&phi_cm<0.4)";
                        const char *conditionBoundR = "(-0.8>phi_cm||phi_cm>0.4)";
                        const char *selectionBoundL = Form("bound_%d_%d_%d_%d",iSys,kleft,iCutTTA,iParticle);
                        const char *selectionBoundR = Form("bound_%d_%d_%d_%d",iSys,kright,iCutTTA,iParticle);

                        if (iCutTTA==kttaAll)
                        {
                          const char *conditionBoundT0  = "(theta_lab*TMath::RadToDeg()>=0&&theta_lab*TMath::RadToDeg()<20)";
                          const char *conditionBoundT20 = "(theta_lab*TMath::RadToDeg()>=20&&theta_lab*TMath::RadToDeg()<40)";
                          const char *conditionBoundT40 = "(theta_lab*TMath::RadToDeg()>=40&&theta_lab*TMath::RadToDeg()<80)";

                          const char *selectionBoundLT0  = Form("bound_%d_%d_%d_%d",iSys,kleft, ktta0, iParticle);
                          const char *selectionBoundLT20 = Form("bound_%d_%d_%d_%d",iSys,kleft, ktta20,iParticle);
                          const char *selectionBoundLT40 = Form("bound_%d_%d_%d_%d",iSys,kleft, ktta40,iParticle);
                          const char *selectionBoundRT0  = Form("bound_%d_%d_%d_%d",iSys,kright,ktta0, iParticle);
                          const char *selectionBoundRT20 = Form("bound_%d_%d_%d_%d",iSys,kright,ktta20,iParticle);
                          const char *selectionBoundRT40 = Form("bound_%d_%d_%d_%d",iSys,kright,ktta40,iParticle);

                          selectionBound = Form("%s*%s*%s+%s*%s*%s+%s*%s*%s+%s*%s*%s+%s*%s*%s+%s*%s*%s"
                              ,conditionBoundL,conditionBoundT0 ,selectionBoundLT0
                              ,conditionBoundL,conditionBoundT20,selectionBoundLT20
                              ,conditionBoundL,conditionBoundT40,selectionBoundLT40
                              ,conditionBoundR,conditionBoundT0 ,selectionBoundRT0
                              ,conditionBoundR,conditionBoundT20,selectionBoundRT20
                              ,conditionBoundR,conditionBoundT40,selectionBoundRT40);
                        }
                        else {
                          selectionBound = Form("%s*%s+%s*%s",conditionBoundL,selectionBoundL,conditionBoundR,selectionBoundR);
                        }
                      }
                      else if (iCutTTA==kttaAll) {
                        const char *conditionBoundT0  = "(theta_lab*TMath::RadToDeg()>=0&&theta_lab*TMath::RadToDeg()<20)";
                        const char *conditionBoundT20 = "(theta_lab*TMath::RadToDeg()>=20&&theta_lab*TMath::RadToDeg()<40)";
                        const char *conditionBoundT40 = "(theta_lab*TMath::RadToDeg()>=40&&theta_lab*TMath::RadToDeg()<80)";
                        const char *selectionBoundT0  = Form("bound_%d_%d_%d_%d",iSys,iLR,ktta0,iParticle);
                        const char *selectionBoundT20 = Form("bound_%d_%d_%d_%d",iSys,iLR,ktta20,iParticle);
                        const char *selectionBoundT40 = Form("bound_%d_%d_%d_%d",iSys,iLR,ktta40,iParticle);
                        selectionBound = Form("%s*%s+%s*%s+%s*%s",conditionBoundT0,selectionBoundT0,conditionBoundT20,selectionBoundT20,conditionBoundT40,selectionBoundT40);
                      }
                      selection = cut0 * TCut(Form("%s/eff/%d","1",numEventsInAna[iSys])) * selectionBound * cutYP * cutTTA * fParticlePozCut[iParticle];
                    }

                    TString selectionValue = selection.GetTitle();
                    selectionValue.ReplaceAll("PARTICLEA",Form("%d",fParticleA[iParticle]));
                    selection = selectionValue;

                    project(treeParticle,nameY0Part,"fy_cm/(by_cm/2)",selection);

                    histY0Array[iSys][iCutTTA][iCutYP][iParticle] = histY0;
                  }
                }
              }
            }

          }

          if (drawKT || drawCorEKE || drawCorPID || drawPtoaR21 || drawKeoaR21 || drawNZR21)
          {
            for (auto iCutYP : fCutY0Idx)
            {
              if (selCutYP>=0 && selCutYP!=iCutYP) continue;
              TCut cutYP = fCutYPValues[iCutYP];

              for (auto iCutTTA : selCutTTAIdx)
              {
                TCut cutTTA = fCutTTAValues[iCutTTA];

                TH2D *histPIDCor = nullptr;
                TH2D *histEKECor = nullptr;
                vector<TH2D *> histEKECorArray;// = nullptr;

                for (auto iParticle : fParticleIdx)
                {
                  const char *nameParticle = fParticleNames[iParticle];
                  const char *titleParticle = fParticleTitles[iParticle];

                  auto nameFileParticle = Form("data2/%s/sys%d_%s_%s_%s.%s.ana.%s.root",anaFName,sys,anaFName,lrFName,multFName,spVersion,nameParticle);
                  //cout_info << "File : " << nameFileParticle << endl;

                  auto treeParticle = new TChain(nameParticle);
                  treeParticle -> Add(nameFileParticle);

                  if (drawCorPID)
                  {
                    auto namePIDCorPart = makeName("pidCor",iAna,iLR,iMult,iSys,iCutTTA,iCutYP,iParticle);
                    auto titlePIDCorPart = makeTitle("Cor.",iAna,iLR,iMult,iSys,iCutTTA,iCutYP,iParticle);

                    auto histPIDCorPart = new TH2D(namePIDCorPart,Form("%s;p/Z (MeV/c);dE/dx;",titlePIDCorPart),nbinsPID,0,pozMax,nbinsPID,0,dEdxMax);
                    histPIDCorPart -> SetMinimum(0.5);
                    histPIDCorPart -> SetMaximum(800);
                    giSys = iSys;
                    giParticle = iParticle;
                    TCut selection = cut0 * cutYP * cutTTA * TCut(Form("%s/eff/%d",finalProbString,numEventsInAna[iSys])) * fParticlePozCut[iParticle];
                    if (drawCorPID22222)
                      selection = cut0 * cutYP * cutTTA * TCut(Form("%s",finalProbString)) * fParticlePozCut[iParticle];
                    project(treeParticle,namePIDCorPart,"dedx:p_lab",selection);

                    if (histPIDCor==nullptr) {
                      auto namePIDCor = makeName("pidCor",iAna,iLR,iMult,iSys,iCutTTA,iCutYP);
                      histPIDCor = (TH2D *) histPIDCorPart -> Clone(namePIDCor);
                      histPIDCor -> SetTitle(Form("%s;p/Z (MeV/c);dE/dx;",titlePIDCorPart));
                    }
                    else
                      histPIDCor -> Add(histPIDCorPart);
                  }

                  if (drawCorEKE)
                  {
                    auto nameEKECorPart = makeName("ekeCor",iAna,iLR,iMult,iSys,iCutTTA,iCutYP,iParticle);
                    auto titleEKECorPart = makeTitle("Cor.",iAna,iLR,iMult,iSys,iCutTTA,iCutYP,iParticle);

                    //auto histEKECorPart = new TH2D(nameEKECorPart,Form("%s;KE_{Lab} (MeV);dE/dx;",titleEKECorPart),nbinsKeoa,0,4*keoaMax,nbinsPID,0,.5*dEdxMax);
                    auto histEKECorPart = new TH2D(nameEKECorPart,Form("%s;KE_{Lab} (MeV);dE/dx;",titleEKECorPart),nbinsKeoa,0,1500,nbinsPID,0,100);
                    giSys = iSys;
                    giParticle = iParticle;
                    TCut selection = cut0 * TCut(Form("%s/eff/%d",finalProbString,numEventsInAna[iSys])) * cutYP * cutTTA * fParticlePozCut[iParticle];
                    auto partz = fParticleZ[iParticle];
                    auto partm = fParticleMass[iParticle];
                    auto parta = fParticleA[iParticle];
                    //const char *expression = Form("dedx:(sqrt((p_lab*%d)*(p_lab*%d)+%f*%f)-%f)/%d",partz,partz,partm,partm,partm,parta);
                    const char *expression = Form("dedx:(sqrt((p_lab*%d)*(p_lab*%d)+%f*%f)-%f)",partz,partz,partm,partm,partm);
                    project(treeParticle,nameEKECorPart,expression,selection);

                    if (histEKECor==nullptr) {
                      auto nameEKECor = makeName("ekeCor",iAna,iLR,iMult,iSys,iCutTTA,iCutYP);
                      histEKECor = (TH2D *) histEKECorPart -> Clone(nameEKECor);
                      histEKECor -> SetTitle(Form("%s;KE_{Lab} (MeV);dE/dx;",titleEKECorPart));
                    }
                    else
                      histEKECor -> Add(histEKECorPart);

                    histEKECorArray.push_back(histEKECorPart);
                  }

                  if (drawPtoaR21)
                  {
                    auto namePtoaPart = makeName("ptoa",iAna,iLR,iMult,iSys,iCutTTA,iCutYP,iParticle);
                    auto titlePtoaPart = makeTitle("Cor.",iAna,iLR,iMult,iSys,iCutTTA,iCutYP,iParticle);

                    //auto histPtoa = new TH1D(namePtoaPart,Form("%s;p_{T}/A (MeV/c);",titlePtoaPart),nbinsPtoa,0,ptoaMax);
                    auto histPtoa = new TH1D(namePtoaPart,Form("%s;p_{T}/A (MeV/c);",titlePtoaPart),nbinsPtoa,0,ptoaMax);
                    giSys = iSys;
                    giParticle = iParticle;
                    TCut selection = cut0 * TCut(Form("%s/eff/%d",finalProbString,numEventsInAna[iSys])) * cutYP * cutTTA  * fParticlePozCut[iParticle];

                    if (fUseHandCut) {
                      const char *nameBound = Form("bound_%d_%d_%d_%d",iSys,iLR,iCutTTA,iParticle);
                      auto selectionBound = nameBound;

                      if (iLR==klr)
                      {
                        const char *conditionBoundL = "(-0.8<phi_cm&&phi_cm<0.4)";
                        const char *conditionBoundR = "(-0.8>phi_cm||phi_cm>0.4)";
                        const char *selectionBoundL = Form("bound_%d_%d_%d_%d",iSys,kleft,iCutTTA,iParticle);
                        const char *selectionBoundR = Form("bound_%d_%d_%d_%d",iSys,kright,iCutTTA,iParticle);

                        //cout << iCutTTA << " " << kttaAll << endl;
                        if (iCutTTA==kttaAll)
                        {
                          const char *conditionBoundT0  = "(theta_lab*TMath::RadToDeg()>=0&&theta_lab*TMath::RadToDeg()<20)";
                          const char *conditionBoundT20 = "(theta_lab*TMath::RadToDeg()>=20&&theta_lab*TMath::RadToDeg()<40)";
                          const char *conditionBoundT40 = "(theta_lab*TMath::RadToDeg()>=40&&theta_lab*TMath::RadToDeg()<80)";

                          const char *selectionBoundLT0  = Form("bound_%d_%d_%d_%d",iSys,kleft, ktta0, iParticle);
                          const char *selectionBoundLT20 = Form("bound_%d_%d_%d_%d",iSys,kleft, ktta20,iParticle);
                          const char *selectionBoundLT40 = Form("bound_%d_%d_%d_%d",iSys,kleft, ktta40,iParticle);
                          const char *selectionBoundRT0  = Form("bound_%d_%d_%d_%d",iSys,kright,ktta0, iParticle);
                          const char *selectionBoundRT20 = Form("bound_%d_%d_%d_%d",iSys,kright,ktta20,iParticle);
                          const char *selectionBoundRT40 = Form("bound_%d_%d_%d_%d",iSys,kright,ktta40,iParticle);

                          selectionBound = Form("%s*%s*%s+%s*%s*%s+%s*%s*%s+%s*%s*%s+%s*%s*%s+%s*%s*%s"
                              ,conditionBoundL,conditionBoundT0 ,selectionBoundLT0
                              ,conditionBoundL,conditionBoundT20,selectionBoundLT20
                              ,conditionBoundL,conditionBoundT40,selectionBoundLT40
                              ,conditionBoundR,conditionBoundT0 ,selectionBoundRT0
                              ,conditionBoundR,conditionBoundT20,selectionBoundRT20
                              ,conditionBoundR,conditionBoundT40,selectionBoundRT40);
                        }
                        else {
                          selectionBound = Form("%s*%s+%s*%s",conditionBoundL,selectionBoundL,conditionBoundR,selectionBoundR);
                        }
                      }
                      else if (iCutTTA==kttaAll) {
                        const char *conditionBoundT0  = "(theta_lab*TMath::RadToDeg()>=0&&theta_lab*TMath::RadToDeg()<20)";
                        const char *conditionBoundT20 = "(theta_lab*TMath::RadToDeg()>=20&&theta_lab*TMath::RadToDeg()<40)";
                        const char *conditionBoundT40 = "(theta_lab*TMath::RadToDeg()>=40&&theta_lab*TMath::RadToDeg()<80)";
                        const char *selectionBoundT0  = Form("bound_%d_%d_%d_%d",iSys,iLR,ktta0,iParticle);
                        const char *selectionBoundT20 = Form("bound_%d_%d_%d_%d",iSys,iLR,ktta20,iParticle);
                        const char *selectionBoundT40 = Form("bound_%d_%d_%d_%d",iSys,iLR,ktta40,iParticle);
                        selectionBound = Form("%s*%s+%s*%s+%s*%s",conditionBoundT0,selectionBoundT0,conditionBoundT20,selectionBoundT20,conditionBoundT40,selectionBoundT40);
                      }
                      selection = cut0 * TCut(Form("%s/eff/%d","1",numEventsInAna[iSys])) * selectionBound * cutYP * cutTTA * fParticlePozCut[iParticle];
                    }

                    project(treeParticle,namePtoaPart,Form("pt_cm/%d",fParticleA[iParticle]),selection);

                    histPtoaArray[iSys][iCutTTA][iCutYP][iParticle] = histPtoa;
                  }

                  if (drawKeoaR21)
                  {
                    auto nameKeoaPart = makeName("keoa",iAna,iLR,iMult,iSys,iCutTTA,iCutYP,iParticle);
                    auto titleKeoaPart = makeTitle("Cor.",iAna,iLR,iMult,iSys,iCutTTA,iCutYP,iParticle);

                    auto histKeoa = new TH1D(nameKeoaPart,Form("%s;KE_{Lab}/A (MeV);",titleKeoaPart),nbinsKeoaR21,0,keoaMax);
                    giSys = iSys;
                    giParticle = iParticle;
                    TCut selection = cut0 * TCut(Form("%s/eff/%d",finalProbString,numEventsInAna[iSys])) * cutYP * cutTTA * fParticlePozCut[iParticle];

                    if (fUseHandCut) {
                      const char *nameBound = Form("bound_%d_%d_%d_%d",iSys,iLR,iCutTTA,iParticle);
                      auto selectionBound = nameBound;

                      if (iLR==klr)
                      {
                        const char *conditionBoundL = "(-0.8<phi_cm&&phi_cm<0.4)";
                        const char *conditionBoundR = "(-0.8>phi_cm||phi_cm>0.4)";
                        const char *selectionBoundL = Form("bound_%d_%d_%d_%d",iSys,kleft,iCutTTA,iParticle);
                        const char *selectionBoundR = Form("bound_%d_%d_%d_%d",iSys,kright,iCutTTA,iParticle);

                        //cout << iCutTTA << " " << kttaAll << endl;
                        if (iCutTTA==kttaAll)
                        {
                          const char *conditionBoundT0  = "(theta_lab*TMath::RadToDeg()>=0&&theta_lab*TMath::RadToDeg()<20)";
                          const char *conditionBoundT20 = "(theta_lab*TMath::RadToDeg()>=20&&theta_lab*TMath::RadToDeg()<40)";
                          const char *conditionBoundT40 = "(theta_lab*TMath::RadToDeg()>=40&&theta_lab*TMath::RadToDeg()<80)";

                          const char *selectionBoundLT0  = Form("bound_%d_%d_%d_%d",iSys,kleft, ktta0, iParticle);
                          const char *selectionBoundLT20 = Form("bound_%d_%d_%d_%d",iSys,kleft, ktta20,iParticle);
                          const char *selectionBoundLT40 = Form("bound_%d_%d_%d_%d",iSys,kleft, ktta40,iParticle);
                          const char *selectionBoundRT0  = Form("bound_%d_%d_%d_%d",iSys,kright,ktta0, iParticle);
                          const char *selectionBoundRT20 = Form("bound_%d_%d_%d_%d",iSys,kright,ktta20,iParticle);
                          const char *selectionBoundRT40 = Form("bound_%d_%d_%d_%d",iSys,kright,ktta40,iParticle);

                          selectionBound = Form("%s*%s*%s+%s*%s*%s+%s*%s*%s+%s*%s*%s+%s*%s*%s+%s*%s*%s"
                              ,conditionBoundL,conditionBoundT0 ,selectionBoundLT0
                              ,conditionBoundL,conditionBoundT20,selectionBoundLT20
                              ,conditionBoundL,conditionBoundT40,selectionBoundLT40
                              ,conditionBoundR,conditionBoundT0 ,selectionBoundRT0
                              ,conditionBoundR,conditionBoundT20,selectionBoundRT20
                              ,conditionBoundR,conditionBoundT40,selectionBoundRT40);
                        }
                        else {
                          selectionBound = Form("%s*%s+%s*%s",conditionBoundL,selectionBoundL,conditionBoundR,selectionBoundR);
                        }
                      }
                      else if (iCutTTA==kttaAll) {
                        const char *conditionBoundT0  = "(theta_lab*TMath::RadToDeg()>=0&&theta_lab*TMath::RadToDeg()<20)";
                        const char *conditionBoundT20 = "(theta_lab*TMath::RadToDeg()>=20&&theta_lab*TMath::RadToDeg()<40)";
                        const char *conditionBoundT40 = "(theta_lab*TMath::RadToDeg()>=40&&theta_lab*TMath::RadToDeg()<80)";
                        const char *selectionBoundT0  = Form("bound_%d_%d_%d_%d",iSys,iLR,ktta0,iParticle);
                        const char *selectionBoundT20 = Form("bound_%d_%d_%d_%d",iSys,iLR,ktta20,iParticle);
                        const char *selectionBoundT40 = Form("bound_%d_%d_%d_%d",iSys,iLR,ktta40,iParticle);
                        selectionBound = Form("%s*%s+%s*%s+%s*%s",conditionBoundT0,selectionBoundT0,conditionBoundT20,selectionBoundT20,conditionBoundT40,selectionBoundT40);
                      }
                      selection = cut0 * TCut(Form("%s/eff/%d","1",numEventsInAna[iSys])) * selectionBound * cutYP * cutTTA * fParticlePozCut[iParticle];
                    }

                    auto partz = fParticleZ[iParticle];
                    auto partm = fParticleMass[iParticle];
                    auto parta = fParticleA[iParticle];
                    const char *expression = Form("(sqrt((p_lab*%d)*(p_lab*%d)+%f*%f)-%f)/%d",partz,partz,partm,partm,partm,parta);
                    project(treeParticle,nameKeoaPart,expression,selection);

                    histKeoaArray[iSys][iCutTTA][iCutYP][iParticle] = histKeoa;
                  }

                  if (drawNZR21)
                  {
                    auto nameYPR21Part = makeName("ypR21",iAna,iLR,iMult,iSys,iCutTTA,iCutYP,iParticle);
                    auto titleYPR21Part = makeTitle("Cor.",iAna,iLR,iMult,iSys,iCutTTA,iCutYP,iParticle);

                    auto histYPR21Part = new TH2D(nameYPR21Part,Form("%s;y_{0};p_{T}/A;",titleYPR21Part),nbinsNZYP,-1,1,nbinsNZYP,0,400);
                    giSys = iSys;
                    giParticle = iParticle;
                    TCut selection = cut0 * TCut(Form("%s/eff/%d",finalProbString,numEventsInAna[iSys])) * cutYP * cutTTA * fParticlePozCut[iParticle];
                    if (fUseHandCut) {
                      const char *nameBound = Form("bound_%d_%d_%d_%d",iSys,iLR,iCutTTA,iParticle);
                      auto selectionBound = nameBound;

                      if (iLR==klr)
                      {
                        const char *conditionBoundL = "(-0.8<phi_cm&&phi_cm<0.4)";
                        const char *conditionBoundR = "(-0.8>phi_cm||phi_cm>0.4)";
                        const char *selectionBoundL = Form("bound_%d_%d_%d_%d",iSys,kleft,iCutTTA,iParticle);
                        const char *selectionBoundR = Form("bound_%d_%d_%d_%d",iSys,kright,iCutTTA,iParticle);

                        //cout << iCutTTA << " " << kttaAll << endl;
                        if (iCutTTA==kttaAll)
                        {
                          const char *conditionBoundT0  = "(theta_lab*TMath::RadToDeg()>=0&&theta_lab*TMath::RadToDeg()<20)";
                          const char *conditionBoundT20 = "(theta_lab*TMath::RadToDeg()>=20&&theta_lab*TMath::RadToDeg()<40)";
                          const char *conditionBoundT40 = "(theta_lab*TMath::RadToDeg()>=40&&theta_lab*TMath::RadToDeg()<80)";

                          const char *selectionBoundLT0  = Form("bound_%d_%d_%d_%d",iSys,kleft, ktta0, iParticle);
                          const char *selectionBoundLT20 = Form("bound_%d_%d_%d_%d",iSys,kleft, ktta20,iParticle);
                          const char *selectionBoundLT40 = Form("bound_%d_%d_%d_%d",iSys,kleft, ktta40,iParticle);
                          const char *selectionBoundRT0  = Form("bound_%d_%d_%d_%d",iSys,kright,ktta0, iParticle);
                          const char *selectionBoundRT20 = Form("bound_%d_%d_%d_%d",iSys,kright,ktta20,iParticle);
                          const char *selectionBoundRT40 = Form("bound_%d_%d_%d_%d",iSys,kright,ktta40,iParticle);

                          selectionBound = Form("%s*%s*%s+%s*%s*%s+%s*%s*%s+%s*%s*%s+%s*%s*%s+%s*%s*%s"
                              ,conditionBoundL,conditionBoundT0 ,selectionBoundLT0
                              ,conditionBoundL,conditionBoundT20,selectionBoundLT20
                              ,conditionBoundL,conditionBoundT40,selectionBoundLT40
                              ,conditionBoundR,conditionBoundT0 ,selectionBoundRT0
                              ,conditionBoundR,conditionBoundT20,selectionBoundRT20
                              ,conditionBoundR,conditionBoundT40,selectionBoundRT40);
                        }
                        else {
                          selectionBound = Form("%s*%s+%s*%s",conditionBoundL,selectionBoundL,conditionBoundR,selectionBoundR);
                        }
                      }
                      else if (iCutTTA==kttaAll) {
                        const char *conditionBoundT0  = "(theta_lab*TMath::RadToDeg()>=0&&theta_lab*TMath::RadToDeg()<20)";
                        const char *conditionBoundT20 = "(theta_lab*TMath::RadToDeg()>=20&&theta_lab*TMath::RadToDeg()<40)";
                        const char *conditionBoundT40 = "(theta_lab*TMath::RadToDeg()>=40&&theta_lab*TMath::RadToDeg()<80)";
                        const char *selectionBoundT0  = Form("bound_%d_%d_%d_%d",iSys,iLR,ktta0,iParticle);
                        const char *selectionBoundT20 = Form("bound_%d_%d_%d_%d",iSys,iLR,ktta20,iParticle);
                        const char *selectionBoundT40 = Form("bound_%d_%d_%d_%d",iSys,iLR,ktta40,iParticle);
                        selectionBound = Form("%s*%s+%s*%s+%s*%s",conditionBoundT0,selectionBoundT0,conditionBoundT20,selectionBoundT20,conditionBoundT40,selectionBoundT40);
                      }
                      selection = cut0 * TCut(Form("%s/eff/%d","1",numEventsInAna[iSys])) * selectionBound * cutYP * cutTTA * fParticlePozCut[iParticle];
                    }
                    project(treeParticle,nameYPR21Part,Form("pt_cm/%d:fy_cm/(by_cm/2)",fParticleA[iParticle]),selection);

                    histYPR21Array[iSys][iCutTTA][iCutYP][iParticle] = histYPR21Part;
                  }
                }

                if (drawCorPID)
                {
                  auto namePIDCor = makeName("pidCor",iAna,iLR,iMult,iSys,iCutTTA,iCutYP);
                  auto cvsPIDCor = makeCvs2(namePIDCor,1000,700);
                  cvsPIDCor -> SetLogz();
                  if (!drawCorPID22222) {
                    histPIDCor -> SetMaximum(0.08);
                    histPIDCor -> SetMinimum(0.0001);
                  }
                  histPIDCor -> Draw("colz");
                  if (drawGuideLine)
                    for (auto iParticle : fParticleIdx)
                      graphPIDMean[iSys][iParticle] -> Draw("samel");
                }

                if (drawCorEKE)
                {
                  auto nameEKECor = makeName("ekeCor",iAna,iLR,iMult,iSys,iCutTTA,iCutYP);
                  auto cvsEKECor = makeCvs2(nameEKECor,1500,1000);
                  CvsDivide(cvsEKECor,3,2);

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
                  const char *namePtoa = makeName("ptoaAll",iAna,iLR,iMult,iSys,iCutTTA,iCutYP);
                  auto cvsR21 = makeCvs(namePtoa);

                  double histMax = 0;
                  for (auto iParticle : fParticleIdx) {
                    auto max0 = histPtoaArray[iSys][iCutTTA][iCutYP][iParticle] -> GetMaximum();
                    if (histMax < max0)
                      histMax = max0;
                  }
                  for (auto iParticle : fParticleIdx) {
                    auto hist = histPtoaArray[iSys][iCutTTA][iCutYP][iParticle];
                    hist -> SetMaximum(histMax*1.05);
                    hist -> SetMinimum(0);
                    hist -> SetLineColor(iParticle+1);
                    hist -> SetMarkerColor(iParticle+1);
                    hist -> SetMarkerStyle(20+iParticle);

                    if (iParticle==0) { hist -> SetMarkerStyle(24); hist -> SetMarkerColor(kRed     ); hist -> SetLineColor(kRed     ); }
                    if (iParticle==1) { hist -> SetMarkerStyle(25); hist -> SetMarkerColor(kBlue    ); hist -> SetLineColor(kBlue    ); }
                    if (iParticle==2) { hist -> SetMarkerStyle(26); hist -> SetMarkerColor(kSpring-6); hist -> SetLineColor(kSpring-6); }
                    if (iParticle==3) { hist -> SetMarkerStyle(30); hist -> SetMarkerColor(kOrange-3); hist -> SetLineColor(kOrange-3); }
                    if (iParticle==4) { hist -> SetMarkerStyle(28); hist -> SetMarkerColor(kViolet-5); hist -> SetLineColor(kViolet-5); }

                    if (iParticle==0)
                      hist -> Draw("plhist");
                    else
                      hist -> Draw("sameplhist");
                  }
                }

              }
            }
          }
        }


        if (drawY0R21)
        {
          for (auto iCutYP : selCutPtoaIdx)
          {
            //if (selCutYP>=0 && selCutYP!=iCutYP) continue;

            for (auto iCutTTA : selCutTTAIdx)
            {
              if (selCutTTA>=0 && selCutTTA!=iCutTTA) continue;

              for (auto iComb : selSysCombIndx)
              {
                if (selCombR21>=0 && selCombR21!=iComb) continue;

                auto iSys2 = fSysCombIdx[iComb][0];
                auto iSys1 = fSysCombIdx[iComb][1];
                auto iSysComb = 100 + iSys2*10 + iSys1;

                const char *nameR21 = makeName("r21_y0",iAna,iLR,iMult,iSysComb,iCutTTA,iCutYP);
                auto titleR21 = makeTitle("Cor.",iAna,iLR,iMult,iSysComb,iCutTTA,iCutYP);

                auto cvs = makeCvs(nameR21);

                auto hist = new TH2D(nameR21,Form("%s;y_{0};R_{21}",titleR21),100,-.5,1,100,0,2.0);
                hist -> Draw();

                auto legend = new TLegend();
                for (auto iParticle : fParticleIdx)
                {
                  const char *nameParticle = fParticleNames[iParticle];

                  auto hist1 = histY0Array[iSys2][iCutTTA][iCutYP][iParticle];
                  auto hist2 = histY0Array[iSys1][iCutTTA][iCutYP][iParticle];
                  /*
                     if (iParticle==0) {
                     hist1 -> Draw("hist");
                     hist2 -> Draw("histsamel");
                     }
                     else {
                     hist1 -> Draw("histsamel");
                     hist2 -> Draw("histsamel");
                     }
                   */

                  const char *nameR21Part = Form("%s_%s",nameR21,nameParticle);

                  auto hist0 = (TH1D *) hist1 -> Clone(nameR21Part);
                  hist0 -> Divide(hist2);
                  auto graph = new TGraph();
                  for (auto bin=1; bin<=6; ++bin)
                  {
                    graph -> SetPoint(graph -> GetN(), hist0 -> GetXaxis() -> GetBinCenter(bin), hist0 -> GetBinContent(bin));
                  }

                  if (iParticle==0) { graph -> SetMarkerStyle(24); graph -> SetMarkerColor(kRed     ); graph -> SetLineColor(kRed     ); }
                  if (iParticle==1) { graph -> SetMarkerStyle(25); graph -> SetMarkerColor(kBlue    ); graph -> SetLineColor(kBlue    ); }
                  if (iParticle==2) { graph -> SetMarkerStyle(26); graph -> SetMarkerColor(kSpring-6); graph -> SetLineColor(kSpring-6); }
                  if (iParticle==3) { graph -> SetMarkerStyle(30); graph -> SetMarkerColor(kOrange-3); graph -> SetLineColor(kOrange-3); }
                  if (iParticle==4) { graph -> SetMarkerStyle(28); graph -> SetMarkerColor(kViolet-5); graph -> SetLineColor(kViolet-5); }

                  graph -> Draw("samepl");
                  legend -> AddEntry(graph,fParticleNames[iParticle],"l");
                }
                makeLegend(cvs,legend) -> Draw();
              }
            }
          }
        }

        if (drawPtoaR21)
        {
          for (auto iCutYP : fCutY0Idx)
          {
            if (selCutYP>=0 && selCutYP!=iCutYP) continue;

            for (auto iCutTTA : selCutTTAIdx)
            {
              if (selCutTTA>=0 && selCutTTA!=iCutTTA) continue;

              for (auto iComb : selSysCombIndx)
              {
                if (selCombR21>=0 && selCombR21!=iComb) continue;

                auto iSys2 = fSysCombIdx[iComb][0];
                auto iSys1 = fSysCombIdx[iComb][1];
                auto iSysComb = 100 + iSys2*10 + iSys1;

                const char *nameR21 = makeName("r21_ptoa",iAna,iLR,iMult,iSysComb,iCutTTA,iCutYP);
                auto titleR21 = makeTitle("Cor.",iAna,iLR,iMult,iSysComb,iCutTTA,iCutYP);

                //auto cvs = makeCvs(nameR21,680,550);
                auto cvs = makeCvs(nameR21);

                auto hist = new TH2D(nameR21,Form("%s;p_{T}/A (MeV/c);R_{21}",titleR21),100,0,ptoaMax,100,0,2.0);
                hist -> Draw();

                auto legend = new TLegend();
                for (auto iParticle : fParticleIdx)
                {
                  const char *nameParticle = fParticleNames[iParticle];

                  auto hist1 = histPtoaArray[iSys2][iCutTTA][iCutYP][iParticle];
                  auto hist2 = histPtoaArray[iSys1][iCutTTA][iCutYP][iParticle];

                  const char *nameR21Part = Form("%s_%s",nameR21,nameParticle);

                  auto hist0 = (TH1D *) hist1 -> Clone(nameR21Part);
                  hist0 -> Divide(hist2);
                  auto graph = new TGraph();
                  //for (auto bin=1; bin<=8; ++bin)
                  for (auto bin=1; bin<=nbinsPtoa; ++bin)
                  {
                    if (hist1->GetBinContent(bin)<0.04||hist2->GetBinContent(bin)<0.04) {
                      hist0 -> SetBinContent(bin,0);
                    }
                    else
                      graph -> SetPoint(graph -> GetN(), hist0 -> GetXaxis() -> GetBinCenter(bin), hist0 -> GetBinContent(bin));
                  }

                  if (iParticle==0) { graph -> SetMarkerStyle(24); graph -> SetMarkerColor(kRed     ); graph -> SetLineColor(kRed     ); }
                  if (iParticle==1) { graph -> SetMarkerStyle(25); graph -> SetMarkerColor(kBlue    ); graph -> SetLineColor(kBlue    ); }
                  if (iParticle==2) { graph -> SetMarkerStyle(26); graph -> SetMarkerColor(kSpring-6); graph -> SetLineColor(kSpring-6); }
                  if (iParticle==3) { graph -> SetMarkerStyle(30); graph -> SetMarkerColor(kOrange-3); graph -> SetLineColor(kOrange-3); }
                  if (iParticle==4) { graph -> SetMarkerStyle(28); graph -> SetMarkerColor(kViolet-5); graph -> SetLineColor(kViolet-5); }

                  graph -> Draw("samepl");
                  legend -> AddEntry(graph,fParticleNames[iParticle],"l");
                }
                makeLegend(cvs,legend) -> Draw();
              }
            }
          }
        }


        if (drawKeoaR21)
        {
          for (auto iCutYP : fCutY0Idx)
          {
            if (selCutYP>=0 && selCutYP!=iCutYP) continue;

            for (auto iCutTTA : selCutTTAIdx)
            {
              if (selCutTTA>=0 && selCutTTA!=iCutTTA) continue;

              for (auto iComb : selSysCombIndx)
              {
                if (selCombR21>=0 && selCombR21!=iComb) continue;

                auto iSys2 = fSysCombIdx[iComb][0];
                auto iSys1 = fSysCombIdx[iComb][1];
                auto iSysComb = 100 + iSys2*10 + iSys1;

                const char *nameR21 = makeName("r21_keoa",iAna,iLR,iMult,iSysComb,iCutTTA,iCutYP);
                auto titleR21 = makeTitle("Cor.",iAna,iLR,iMult,iSysComb,iCutTTA,iCutYP);

                auto cvs = makeCvs(nameR21,680,550);
                //cvs -> SetGrid(1,1);

                auto hist = new TH2D(nameR21,Form("%s;KE_{Lab}/A (MeV);R_{21}",titleR21),100,0,400,100,0,2.0);
                hist -> Draw();

                auto legend = new TLegend();
                for (auto iParticle : fParticleIdx)
                {
                  const char *nameParticle = fParticleNames[iParticle];

                  auto hist1 = histKeoaArray[iSys2][iCutTTA][iCutYP][iParticle];
                  auto hist2 = histKeoaArray[iSys1][iCutTTA][iCutYP][iParticle];

                  const char *nameR21Part = Form("%s_%s",nameR21,nameParticle);

                  auto hist0 = (TH1D *) hist1 -> Clone(nameR21Part);
                  hist0 -> Divide(hist2);

                  hist0 -> SetMarkerSize(1.3);

                  if (iParticle==0) { hist0 -> SetMarkerStyle(24); hist0 -> SetMarkerColor(kBlack   ); hist0 -> SetLineColor(kBlack   ); }
                  if (iParticle==1) { hist0 -> SetMarkerStyle(25); hist0 -> SetMarkerColor(kRed     ); hist0 -> SetLineColor(kRed     ); }
                  if (iParticle==2) { hist0 -> SetMarkerStyle(26); hist0 -> SetMarkerColor(kBlue    ); hist0 -> SetLineColor(kBlue    ); }
                  if (iParticle==3) { hist0 -> SetMarkerStyle(30); hist0 -> SetMarkerColor(kSpring-6); hist0 -> SetLineColor(kSpring-6); }
                  if (iParticle==4) { hist0 -> SetMarkerStyle(28); hist0 -> SetMarkerColor(kOrange-3); hist0 -> SetLineColor(kOrange-3); }

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

            for (auto iCutTTA : selCutTTAIdx)
            {
              if (selCutTTA>=0 && selCutTTA!=iCutTTA) continue;

              for (auto iComb : selSysCombIndx)
              {
                if (selCombR21>=0 && selCombR21!=iComb) continue;

                auto iSys2 = fSysCombIdx[iComb][0];
                auto iSys1 = fSysCombIdx[iComb][1];
                auto iSysComb = 100 + iSys2*10 + iSys1;

                const char *nameAlpha = makeName("alpha",iAna,iLR,iMult,iSysComb,iCutTTA,iCutYP);
                const char *nameMBeta = makeName("mbeta",iAna,iLR,iMult,iSysComb,iCutTTA,iCutYP);

                auto cvsAlpha = makeCvs(nameAlpha,2500,2000);
                auto cvsMBeta = makeCvs(nameMBeta,2500,2000);
                //CvsDivideM0(cvsAlpha,6,8);
                //CvsDivideM0(cvsMBeta,6,8);
                CvsDivideM0(cvsAlpha,5,6);
                CvsDivideM0(cvsMBeta,5,6);

                double r21Array[10][10][fNumParticles];
                double r21ErrorArray[10][10][fNumParticles];
                double alphaArray[10][10];
                double mbetaArray[10][10];
                double alphaErrorArray[10][10];
                double mbetaErrorArray[10][10];

                auto binnx = binning(histYPR21Array[iSys2][iCutTTA][iCutYP][0],1);
                auto binny = binning(histYPR21Array[iSys2][iCutTTA][iCutYP][0],2);

                //for (auto ibinX=1; ibinX<=6; ++ibinX)
                for (auto ibinX=2; ibinX<=6; ++ibinX)
                {
                  //for (auto ibinY=1; ibinY<=8; ++ibinY)
                  for (auto ibinY=1; ibinY<=6; ++ibinY)
                  {
                    auto binX = ibinX+2;
                    auto binY = ibinY;

                    //auto graphAll = new TGraph2DErrors();
                    auto graphAll = new TGraph2D();
                    auto graph_z1_inN = new TGraph(); graph_z1_inN -> SetMarkerSize(1.5); graph_z1_inN -> SetMarkerStyle(20); graph_z1_inN -> SetMarkerColor(kBlack);
                    auto graph_z2_inN = new TGraph(); graph_z2_inN -> SetMarkerSize(1.5); graph_z2_inN -> SetMarkerStyle(25); graph_z2_inN -> SetMarkerColor(kRed  );
                    auto graph_n0_inZ = new TGraph(); graph_n0_inZ -> SetMarkerSize(1.5); graph_n0_inZ -> SetMarkerStyle(20); graph_n0_inZ -> SetMarkerColor(kBlack);
                    auto graph_n1_inZ = new TGraph(); graph_n1_inZ -> SetMarkerSize(1.5); graph_n1_inZ -> SetMarkerStyle(25); graph_n1_inZ -> SetMarkerColor(kRed  );
                    auto graph_n2_inZ = new TGraph(); graph_n2_inZ -> SetMarkerSize(1.5); graph_n2_inZ -> SetMarkerStyle(22); graph_n2_inZ -> SetMarkerColor(kBlue );

                    for (auto iParticle : fParticleIdx)
                    {
                      auto value2 = histYPR21Array[iSys2][iCutTTA][iCutYP][iParticle] -> GetBinContent(binX, binY);
                      auto value1 = histYPR21Array[iSys1][iCutTTA][iCutYP][iParticle] -> GetBinContent(binX, binY);
                      auto r21Value = value2/value1;
                      auto r21Error = 1./sqrt(value1);
                      r21Array[binX][binY][iParticle] = r21Value;
                      r21ErrorArray[binX][binY][iParticle] = r21Error;

                      auto zValue = fParticleZ[iParticle];
                      auto nValue = fParticleN[iParticle];

                      if (zValue==1) graph_z1_inN -> SetPoint(graph_z1_inN->GetN(), nValue, r21Value);
                      if (zValue==2) graph_z2_inN -> SetPoint(graph_z2_inN->GetN(), nValue, r21Value);
                      if (nValue==0) graph_n0_inZ -> SetPoint(graph_n0_inZ->GetN(), zValue, r21Value);
                      if (nValue==1) graph_n1_inZ -> SetPoint(graph_n1_inZ->GetN(), zValue, r21Value);
                      if (nValue==2) graph_n2_inZ -> SetPoint(graph_n2_inZ->GetN(), zValue, r21Value);
                      graphAll -> SetPoint(graphAll->GetN(), nValue, zValue, r21Value);

                      //if (zValue==1) graph_z1_inN -> SetPointError(graph_z1_inN->GetN()-1, r21Error, r21Error);
                      //if (zValue==2) graph_z2_inN -> SetPointError(graph_z2_inN->GetN()-1, r21Error, r21Error);
                      //if (nValue==0) graph_n0_inZ -> SetPointError(graph_n0_inZ->GetN()-1, r21Error, r21Error);
                      //if (nValue==1) graph_n1_inZ -> SetPointError(graph_n1_inZ->GetN()-1, r21Error, r21Error);
                      //if (nValue==2) graph_n2_inZ -> SetPointError(graph_n2_inZ->GetN()-1, r21Error, r21Error);
                      //graphAll -> SetPointError(graphAll->GetN()-1, r21Error, r21Error,0);
                      //graphAll -> SetPointError(graphAll->GetN()-1, .5, .5, r21Error);
                      //graphAll -> SetPointError(graphAll->GetN()-1, .5, .5, .5);
                    }

                    //auto idxCvs = (8-ibinY)*6+ibinX;
                    auto idxCvs = (6-ibinY)*5+(ibinX-1);

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
                    //if (ibinX==1&&ibinY==8) {}
                    if (ibinX==2&&ibinY==6) {}
                    else {
                      graph_z1_inN -> Draw("samep");
                      graph_z2_inN -> Draw("samep");
                    }
                    TLegend *legendAlpha = new TLegend();
                    legendAlpha -> AddEntry(graph_z2_inN, Form("(%d,%d)",binX,binY),"");
                    //legendAlpha -> AddEntry(graph_z2_inN, Form("(%d,%d),(%.2f,%d)",binX,binY,binnx.getCenter(ibinX),int(binny.getCenter(ibinY))),"");

                    double alphaValue = 0;
                    double betaValue = 0;

                    double alphaError = 0;
                    double betaError = 0;

                    //TF2 *fitAll = new TF2("fitAll","exp([0]*x+[1]*y)",1,2,0,2);
                    TF2 *fitAll = new TF2("fitAll","[2]*exp([0]*x+[1]*y)",1,2,0,2);
                    if (fitNZR21TG) {
                      graphAll -> Fit(fitAll);

                      alphaValue = fitAll -> GetParameter(0);
                      betaValue = fitAll -> GetParameter(1);
                      alphaArray[binX][binY] = alphaValue;
                      mbetaArray[binX][binY] = -betaValue;

                      alphaError = fitAll -> GetParError(0);
                      betaError = fitAll -> GetParError(1);
                      alphaErrorArray[binX][binY] = alphaError;
                      mbetaErrorArray[binX][binY] = betaError;

                      TF1 *fit1 = new TF1("fit1","exp([0]*x+[1])",-.5,2.5);
                      fit1 -> SetParameters(alphaValue,betaValue*1);
                      fit1 -> SetLineColor(kGray+1);

                      TF1 *fit2 = new TF1("fit1","exp([0]*x+[1])",-.5,2.5);
                      fit2 -> SetParameters(alphaValue,betaValue*2);
                      fit2 -> SetLineColor(kGray+1);

                      //if (ibinX==1&&ibinY==8) {}
                      if (ibinX==2&&ibinY==6) {}
                      else {
                        fit1 -> Draw("samel");
                        fit2 -> Draw("samel");
                      }

                      legendAlpha -> AddEntry(fit2, Form("#alpha = %.3f",alphaValue),"");
                    }
                    else if (fitNZR21) {
                      auto graph1 = graph_z1_inN;
                      auto graph2 = graph_z2_inN;
                      TF1 *fit1 = new TF1("fit1",fPol1Function,0,100,2);
                      TF1 *fit2 = new TF1("fit2",fPol1Function,0,100,2);
                      double fitRange1 = -.5;
                      double fitRange2 = 2.5;
                      ROOT::Math::WrappedMultiTF1 wfit1(*fit1,1);
                      ROOT::Math::WrappedMultiTF1 wfit2(*fit2,1);
                      ROOT::Fit::DataOptions option;
                      ROOT::Fit::DataRange range1(fitRange1,fitRange2);
                      ROOT::Fit::DataRange range2(fitRange1,fitRange2);
                      ROOT::Fit::BinData data1(option, range1);
                      ROOT::Fit::BinData data2(option, range2);
                      ROOT::Fit::FillData(data1, graph1);
                      ROOT::Fit::FillData(data2, graph2);
                      ROOT::Fit::Chi2Function chi21(data1, wfit1);
                      ROOT::Fit::Chi2Function chi22(data2, wfit2);
                      GlobalChi2 globalChi2(chi21, chi22);
                      ROOT::Fit::Fitter fitter;
                      vector<ROOT::Fit::ParameterSettings> parSetting = {
                        ROOT::Fit::ParameterSettings("intercepti1", 0, 0.001, -100, 100),
                        ROOT::Fit::ParameterSettings("intercepti2", 0, 0.001, -100, 100),
                        ROOT::Fit::ParameterSettings("slope"      , 0, 0.001, -100, 100)};
                      fitter.Config().SetParamsSettings(parSetting);
                      fitter.Config().MinimizerOptions().SetPrintLevel(0);
                      fitter.Config().SetMinimizer("Minuit","Minimize");
                      fitter.FitFCN(3, globalChi2, 0, data1.Size()+data2.Size(), true);
                      ROOT::Fit::FitResult result = fitter.Result();
                      fit1 -> SetFitResult(result, fParIndex1);
                      fit1 -> SetRange(range1().first, range1().second);
                      fit1 -> SetLineColor(kGray+1);
                      graph1 -> GetListOfFunctions() -> Add(fit1);
                      fit2 -> SetFitResult(result, fParIndex2);
                      fit2 -> SetRange(range2().first, range2().second);
                      fit2 -> SetLineColor(kGray+1);
                      graph2 -> GetListOfFunctions() -> Add(fit2);
                      //if (ibinX==1&&ibinY==8) {}
                      if (ibinX==2&&ibinY==6) {}
                      else {
                        fit1 -> Draw("samel");
                        fit2 -> Draw("samel");
                      }
                      int countPoints = 0;
                      double rms = 0;
                      double x0, y0;
                      for (auto i=0; i<graph1->GetN(); i++) {
                        graph1 -> GetPoint(i,x0,y0);
                        auto dy = y0 - fit1 -> Eval(x0);
                        rms = rms + dy*dy;
                        countPoints++;
                      }
                      for (auto i=0; i<graph2->GetN(); i++) {
                        graph2 -> GetPoint(i,x0,y0);
                        auto dy = y0 - fit2 -> Eval(x0);
                        rms = rms + dy*dy;
                        countPoints++;
                      }
                      rms = rms/countPoints; 
                      rms = sqrt(rms);
                      legendAlpha -> AddEntry(fit2, Form("#alpha = %.3f",fit2->GetParameter(1)),"");
                      legendAlpha -> AddEntry(fit2, Form("rms = %.3f",rms),"");
                      alphaArray[binX][binY] = fit2->GetParameter(1);
                    }

                    //if (ibinX==1&&ibinY==8)
                    if (ibinX==2&&ibinY==6)
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
                    //if (ibinX==1&&ibinY==8) {}
                    if (ibinX==2&&ibinY==6) {}
                    else {
                      graph_n0_inZ -> Draw("samep");
                      graph_n1_inZ -> Draw("samep");
                      graph_n2_inZ -> Draw("samep");
                    }
                    TLegend *legendBeta = new TLegend();
                    legendBeta -> AddEntry(graph_n2_inZ , Form("(%d,%d)",binX,binY),"");
                    //legendBeta -> AddEntry(graph_n2_inZ, Form("(%d,%d),(%.2f,%d)",binX,binY,binnx.getCenter(ibinX),int(binny.getCenter(ibinY))),"");

                    if (fitNZR21TG) {
                      TF1 *fit1 = new TF1("fit1","exp([0]+[1]*x)",.5,2.5);
                      fit1 -> SetParameters(alphaValue*1,betaValue);
                      fit1 -> SetLineColor(kGray+1);

                      TF1 *fit2 = new TF1("fit1","exp([0]+[1]*x)",.5,2.5);
                      fit2 -> SetParameters(alphaValue*2,betaValue);
                      fit2 -> SetLineColor(kGray+1);

                      //if (ibinX==1&&ibinY==8) {}
                      if (ibinX==2&&ibinY==6) {}
                      else {
                        fit1 -> Draw("samel");
                        fit2 -> Draw("samel");
                      }

                      legendBeta -> AddEntry(fit2, Form("-#beta = %.3f",-betaValue),"");
                    }
                    else if (fitNZR21) {
                      auto graph3 = graph_n1_inZ;
                      auto graph4 = graph_n2_inZ;
                      TF1 *fit3 = new TF1("fit3",fPol1Function,0,100,2);
                      TF1 *fit4 = new TF1("fit4",fPol1Function,0,100,2);
                      double fitRange1 = .5;
                      double fitRange2 = 2.5;
                      ROOT::Math::WrappedMultiTF1 wfit1(*fit3,1);
                      ROOT::Math::WrappedMultiTF1 wfit2(*fit4,1);
                      ROOT::Fit::DataOptions option;
                      ROOT::Fit::DataRange range1(fitRange1,fitRange2);
                      ROOT::Fit::DataRange range2(fitRange1,fitRange2);
                      ROOT::Fit::BinData data1(option, range1);
                      ROOT::Fit::BinData data2(option, range2);
                      ROOT::Fit::FillData(data1, graph3);
                      ROOT::Fit::FillData(data2, graph4);
                      ROOT::Fit::Chi2Function chi21(data1, wfit1);
                      ROOT::Fit::Chi2Function chi22(data2, wfit2);
                      GlobalChi2 globalChi2(chi21, chi22);
                      ROOT::Fit::Fitter fitter;
                      vector<ROOT::Fit::ParameterSettings> parSetting = {
                        ROOT::Fit::ParameterSettings("intercepti1", 0, 0.001, -100, 100),
                        ROOT::Fit::ParameterSettings("intercepti2", 0, 0.001, -100, 100),
                        ROOT::Fit::ParameterSettings("slope"      , 0, 0.001, -100, 100)};
                      fitter.Config().SetParamsSettings(parSetting);
                      fitter.Config().MinimizerOptions().SetPrintLevel(0);
                      fitter.Config().SetMinimizer("Minuit","Minimize");
                      fitter.FitFCN(3, globalChi2, 0, data1.Size()+data2.Size(), true);
                      ROOT::Fit::FitResult result = fitter.Result();
                      fit3 -> SetFitResult(result, fParIndex1);
                      fit3 -> SetRange(range1().first, range1().second);
                      fit3 -> SetLineColor(kGray+1);
                      graph3 -> GetListOfFunctions() -> Add(fit3);
                      fit4 -> SetFitResult(result, fParIndex2);
                      fit4 -> SetRange(range2().first, range2().second);
                      fit4 -> SetLineColor(kGray+1);
                      graph4 -> GetListOfFunctions() -> Add(fit4);
                      //if (ibinX==1&&ibinY==8) {}
                      if (ibinX==2&&ibinY==6) {}
                      else {
                        fit3 -> Draw("samel");
                        fit4 -> Draw("samel");
                      }
                      int countPoints = 0;
                      double rms = 0;
                      double x0, y0;
                      for (auto i=0; i<graph3->GetN(); i++) {
                        graph3 -> GetPoint(i,x0,y0);
                        auto dy = y0 - fit3 -> Eval(x0);
                        rms = rms + dy*dy;
                        countPoints++;
                      }
                      for (auto i=0; i<graph4->GetN(); i++) {
                        graph4 -> GetPoint(i,x0,y0);
                        auto dy = y0 - fit4 -> Eval(x0);
                        rms = rms + dy*dy;
                        countPoints++;
                      }
                      rms = rms/countPoints; 
                      rms = sqrt(rms);
                      legendBeta -> AddEntry(fit4, Form("-#beta = %.3f",-(fit4->GetParameter(1))), "");
                      legendBeta -> AddEntry(fit4, Form("rms = %.3f",rms),"");
                      mbetaArray[binX][binY] = -(fit4->GetParameter(1));
                    }
                    //if (ibinX==1&&ibinY==8)
                    if (ibinX==2&&ibinY==6)
                    {
                      auto legendGraph = new TLegend();
                      legendGraph -> AddEntry(graph_n0_inZ,"N=0","p");
                      legendGraph -> AddEntry(graph_n1_inZ,"N=1","p");
                      legendGraph -> AddEntry(graph_n2_inZ,"N=2","p");
                      makeLegend(cvs,legendGraph,"lt",.2,0,.6,.8) -> Draw();
                    }
                    else
                      //makeLegend(cvs,legendBeta,"rt",0.11,0,.5,.42) -> Draw();
                      makeLegend(cvs,legendBeta,"rt",0,0,.5,.42) -> Draw();
                  }
                }

                auto nameR21YP = makeName("r21_yp",iAna,iLR,iMult,iSysComb,iCutTTA,iCutYP);
                const char *nameAlphaBetaOut = Form("%s/others/alpha_beta_%s.txt",vName.Data(),nameR21YP);
                cout << nameAlphaBetaOut << endl;
                std::ofstream fileAB(nameAlphaBetaOut);

                auto nameR21AB = makeName("r21_ab",iAna,iLR,iMult,iSysComb,iCutTTA,iCutYP);
                auto cvsR21AB = makeCvs(nameR21AB);
                auto titleR21AB = makeTitle("",iAna,iLR,iMult,iSysComb,iCutTTA,iCutYP);
                //(new TH2D(nameR21AB,Form("%s;#alpha;-#beta",titleR21AB),100,0,.5,100,0,.5)) -> Draw();
                (new TH2D(nameR21AB,Form("%s;#alpha;-#beta",titleR21AB),100,0,.4,100,0,.4)) -> Draw();
                auto nameAlphaYP = makeName("alaph_yp",iAna,iLR,iMult,iSysComb,iCutTTA,iCutYP);
                auto nameMBetaYP = makeName("mbeta_yp",iAna,iLR,iMult,iSysComb,iCutTTA,iCutYP);
                auto titleAlphaYP = makeTitle("alaph_yp",iAna,iLR,iMult,iSysComb,iCutTTA,iCutYP);
                auto titleMBetaYP = makeTitle("mbeta_yp",iAna,iLR,iMult,iSysComb,iCutTTA,iCutYP);
                //auto histAlpha = new TH2D(nameAlphaYP,"#alpha;y0;pt/a",nbinsNZYP-2,-.5,1,nbinsNZYP,0,400);
                //auto histMBeta = new TH2D(nameMBetaYP,"-#beta;y0;pt/a",nbinsNZYP-2,-.5,1,nbinsNZYP,0,400);
                auto histAlpha = new TH2D(nameAlphaYP,Form("#alpha %s;y0;pt/a",fSysCombTitles[iComb]),nbinsNZYP-3,-.25,1,nbinsNZYP-2,0,300);
                auto histMBeta = new TH2D(nameMBetaYP,Form("-#beta %s;y0;pt/a",fSysCombTitles[iComb]),nbinsNZYP-3,-.25,1,nbinsNZYP-2,0,300);
                histAlpha -> GetZaxis() -> SetRangeUser(0,0.4);
                histMBeta -> GetZaxis() -> SetRangeUser(0,0.4);

                TGraph *graphAlpha[8];
                TGraph *graphMBeta[8];

                auto legendAB = new TLegend();
                //for (auto ibinY=8; ibinY>=1; --ibinY) {
                for (auto ibinY=6; ibinY>=1; --ibinY) {
                  auto graphAB = new TGraphErrors();
                  //auto graphAB = new TGraph();
                  graphAlpha[ibinY] = new TGraph();
                  graphMBeta[ibinY] = new TGraph();
                  //for (auto ibinX=6; ibinX>=1; --ibinX) {
                  for (auto ibinX=6; ibinX>=2; --ibinX) {
                    auto binX = ibinX+2;
                    auto binY = ibinY;
                    auto alpha = alphaArray[binX][binY];
                    auto mbeta = mbetaArray[binX][binY];
                    auto alphaError = alphaErrorArray[binX][binY];
                    auto mbetaError = mbetaErrorArray[binX][binY];
                    auto r21Error = r21ErrorArray[binX][binY][0] /100;
                    graphAB -> SetPoint(graphAB->GetN(),alpha,mbeta);
                    //graphAB -> SetPointError(graphAB->GetN()-1,r21Error,r21Error);
                    histAlpha -> SetBinContent(ibinX-1,ibinY,alpha);
                    histMBeta -> SetBinContent(ibinX-1,ibinY,mbeta);
                    graphAlpha[ibinY] -> SetPoint(graphAlpha[ibinY]->GetN(), binnx.getCenter(binX), alpha);
                    graphMBeta[ibinY] -> SetPoint(graphMBeta[ibinY]->GetN(), binnx.getCenter(binX), mbeta);
                  }
                  for (auto graph : {graphAlpha[ibinY],graphMBeta[ibinY]}) {
                    graph -> SetLineColor(ibinY);
                    graph -> SetMarkerColor(ibinY);
                    graph -> SetMarkerStyle(20);
                  }
                  graphAB -> SetLineColor(ibinY);
                  graphAB -> SetMarkerColor(ibinY);
                  graphAB -> SetMarkerStyle(25);
                  graphAB -> SetMarkerSize(1.6);
                  graphAB -> Draw("samep");
                  legendAB -> AddEntry(graphAB,Form("%d) p_{T}/A = %.0f ~ %.0f MeV/c",ibinY, binny.lowEdge(ibinY),binny.highEdge(ibinY)));
                }
                makeLegend(cvsR21AB,legendAB,"lt",0,0,0,.3) -> Draw();
                auto f1Linear = new TF1(Form("invf_%s",nameR21AB),"x",0,.4);
                f1Linear -> SetLineStyle(2);
                f1Linear -> SetLineColor(kGray);
                f1Linear -> Draw("samel");

                auto nameABYP = makeName("ab_yp",iAna,iLR,iMult,iSysComb,iCutTTA,iCutYP);
                //auto cvsAB = makeCvs2(nameABYP,1200,550);
                auto cvsAB = makeCvs2(nameABYP,600,1100);
                CvsDivide(cvsAB,1,2);
                cvsAB -> cd(1); histAlpha -> Draw("colz");
                cvsAB -> cd(2); histMBeta -> Draw("colz");

                /*
                auto nameAlpha1 = makeName("alpha1",iAna,iLR,iMult,iSysComb,iCutTTA,iCutYP);
                auto nameMBeta1 = makeName("mbeta1",iAna,iLR,iMult,iSysComb,iCutTTA,iCutYP);
                titleR21AB = makeTitle("",iAna,iLR,iMult,iSysComb,iCutTTA,iCutYP);
                nameAlphaYP = makeName("alaph1_yp",iAna,iLR,iMult,iSysComb,iCutTTA,iCutYP);
                nameMBetaYP = makeName("mbeta1_yp",iAna,iLR,iMult,iSysComb,iCutTTA,iCutYP);
                auto histAlpha1 = new TH2D(nameAlphaYP,Form("%s;y0;#alpha",titleR21AB),100,-.5,1,100,0,0.5);
                auto histMBeta1 = new TH2D(nameMBetaYP,Form("%s;y0;-#beta",titleR21AB),100,-.5,1,100,0,0.5);

                auto cvsAlpha1 = makeCvs(nameAlpha1);
                histAlpha1 -> Draw();
                for (auto ibinY=8; ibinY>=1; --ibinY)
                  graphAlpha[ibinY] -> Draw("samelp");
                legendAB -> Draw();

                auto cvsMBeta1 = makeCvs(nameMBeta1);
                histMBeta1 -> Draw();
                for (auto ibinY=8; ibinY>=1; --ibinY)
                  graphMBeta[ibinY] -> Draw("samelp");
                legendAB -> Draw();
                */



                fileAB << "alpha" << endl;
                fileAB << ","; for (auto binX=3; binX<=8; ++binX) fileAB << ", " << binX; fileAB << endl;
                fileAB << ","; for (auto binX=3; binX<=8; ++binX) fileAB << ", y0=" << Form("%.2f",binnx.lowEdge(binX)) << "~" << Form("%.2f",binnx.highEdge(binX)); fileAB << endl;
                for (auto binY=8; binY>=1; --binY) {
                  fileAB << binY << ",  pt/a=" << int(binny.lowEdge(binY)) << "~" << int(binny.highEdge(binY));
                  for (auto binX=3; binX<=8; ++binX) fileAB << ", " << alphaArray[binX][binY]; fileAB << endl;
                } fileAB << endl;

                fileAB << "-beta" << endl;
                fileAB << ","; for (auto binX=3; binX<=8; ++binX) fileAB << ", " << binX; fileAB << endl;
                fileAB << ","; for (auto binX=3; binX<=8; ++binX) fileAB << ", y0=" << Form("%.2f",binnx.lowEdge(binX)) << "~" << Form("%.2f",binnx.highEdge(binX)); fileAB << endl;
                for (auto binY=8; binY>=1; --binY) {
                  fileAB << binY << ",  pt/a=" << int(binny.lowEdge(binY)) << "~" << int(binny.highEdge(binY));
                  for (auto binX=3; binX<=8; ++binX) fileAB << ", " << mbetaArray[binX][binY]; fileAB << endl;
                } fileAB << endl;

                for (auto iParticle : fParticleIdx)
                {
                  fileAB << "R21_" << fParticleNames[iParticle] << endl;
                  fileAB << ","; for (auto binX=3; binX<=8; ++binX) fileAB << ", " << binX; fileAB << endl;
                  fileAB << ","; for (auto binX=3; binX<=8; ++binX) fileAB << ", y0=" << Form("%.2f",binnx.lowEdge(binX)) << "~" << Form("%.2f",binnx.highEdge(binX)); fileAB << endl;
                  for (auto binY=8; binY>=1; --binY) {
                    fileAB << binY << ",  pt/a=" << int(binny.lowEdge(binY)) << "~" << int(binny.highEdge(binY));
                    for (auto binX=3; binX<=8; ++binX) fileAB << ", " << r21Array[binX][binY][iParticle]; fileAB << endl;
                  } fileAB << endl;
                }

                fileAB.close();
              }
            }
          }
        }



      }
    }
  }

  if (saveCvsPNG)
    saveAll();

  if (saveCvsRoot)
    writeAll();
}

TLegend *makeLegend(TVirtualPad *cvs, TLegend *legend, TString fLegendDrawStyle="rt", double x_offset, double y_offset, double width_fixed, double height_fixed) // jumpto_makel
{
  double fWidthPerLengthLegend = 0.008;
  double fHeightPerEntryLegend = 0.05;
  double fWidthDefaultLegend = 0.10;
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

void CvsDivideM0(TCanvas *cvs, int nx, int ny)
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
        //cvs->cd(i)->SetMargin(fLMargin,fRMargin,fBMargin,fTMargin);
     }
    }
  }
  //return cvs;
}

void CvsDivide(TCanvas *cvs, int nx, int ny)
{
  auto lMargin = cvs -> GetLeftMargin();
  auto rMargin = cvs -> GetRightMargin();
  auto bMargin = cvs -> GetBottomMargin();
  auto tMargin = cvs -> GetTopMargin();

  cvs -> Divide(nx,ny,0.01,0.001);

  for (auto ix=1; ix<=nx; ++ix)
    for (auto iy=1; iy<=ny; ++iy) {
      auto i = ix + nx*(iy-1);
      cvs -> cd(i) -> SetMargin(lMargin,rMargin,bMargin,tMargin);
    }
}
