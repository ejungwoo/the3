#include "KBGlobal.hh"
#include "init_variables.h"
#include "binning.h"

int cvsXOff = 1300;
int giSys = 0;
int giParticle = 0;
vector<TCanvas *> cvsArray;

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

void Divide0(TCanvas *cvs, int nx, int ny);
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

TCanvas *makeCvs2(const char *name, int w=680, int h=550) {
  const char *name0 = Form("cvs_%s",name);
  auto cvs = new TCanvas(name0,name0,cvsXOff+20*(countCvs+1), 20*(countCvs+1), w, h);
  cvs -> SetRightMargin(0.155);
  countCvs++;
  cvsArray.push_back(cvs);
  return cvs;
}

TCanvas *makeCvs(const char *name, int w=680, int h=550) {
  const char *name0 = Form("cvs_%s",name);
  auto cvs = new TCanvas(name0,name0,cvsXOff+20*(countCvs+1), 20*(countCvs+1), w, h);
  //cvs -> SetRightMargin(0.155);
  countCvs++;
  cvsArray.push_back(cvs);
  return cvs;
}

void writeAll()
{
  for (auto cvs : cvsArray)
    cvs -> SaveAs(TString("data_cvs/")+cvs->GetName()+".root");
}

bool fUseHandCut = false;

const char *makeName(const char *mainName, int iAna, int iLR, int iMult, int iSys, int iCutTTA, int iY0, int iPart=-1)
{
  if (fUseHandCut)
    mainName = Form("%s_HPID_",mainName);

  const char *systemName;
  if (iSys>=100) {
    iSys = iSys - 100;
    int iSys2 = int(iSys/10);
    int iSys1 = iSys - iSys2;
    systemName = Form("%so%s",fSystemNames[iSys2],fSystemNames[iSys1]);
  }
  else
    systemName = Form("%s",fSystemNames[iSys]);

  if (iPart<0) {
    const char *name = Form("%s_%s_%s_%s_%s_%s_%s",mainName,fAnaNames[iAna],fLRNames[iLR],fMultNames[iMult],systemName,fCutTTANames[iCutTTA],fCutY0Names[iY0]);
    return name;
  }

  const char *name = Form("%s_%s_%s_%s_%s_%s_%s_%s",mainName,fAnaNames[iAna],fLRNames[iLR],fMultNames[iMult],systemName,fCutTTANames[iCutTTA],fCutY0Names[iY0],fParticleNames[iPart]);
  return name;
}

const char *makeTitle(const char *mainName, int iAna, int iLR, int iMult, int iSys, int iCutTTA, int iY0, int iPart=5)
{
  if (fUseHandCut)
    mainName = Form("%s, Hand-Cut-PID",mainName);

  const char *systemTitle;
  if (iSys>=100) {
    iSys = iSys - 100;
    int iSys2 = int(iSys/10);
    int iSys1 = iSys - iSys2;
    systemTitle = Form("%d/%d",fSystems[iSys2],fSystems[iSys1]);
  }
  else
    systemTitle = Form("%d",fSystems[iSys]);

  const char *partTitle = "";
  if (iPart<0)
    partTitle = Form(", %s",fCutY0Titles[iY0]);

  const char *y0Title = "";
  if (iY0!=0)
    y0Title = Form(", %s",fCutTTATitles[iCutTTA]);

  //const char *title = Form("%s, %s, %s, %s, %s, %s%s%s",mainName,fAnaTitles[iAna],fLRTitles[iLR],fMultTitles[iMult],systemTitle,fCutTTATitles[iCutTTA],y0Title,partTitle);

  const char *title = Form("%s",fCutTTATitles[iCutTTA]);
  if (fUseHandCut)
    title = Form("Hand-Cut-PID, %s",mainName);

  return title;
}

void draw_pid()
{
  int selAna = kf7;
  int selLR = kall;
  int selMult = kn55;
  int selSys = kall;

  const int selCutTTAIdx[] = {0};
  //const int selCutTTAIdx[] = {1,2,3,4};
  int selCutTTA = kall;

  //int selCutY0 = ky02;
  int selCutY0 = kya;

  int selLRIdx[] = {0};
  //int selLRIdx[] = {1};

  bool saveFigures = false;
  TString pathToFigures = "figures3/";
  TString figureFormat = ".png";

  TCut cut0 = "";

  ///////////////////////////////////////////////////////////////////////////////////////////////////

  bool drawRawPID = false;
  bool drawGuideLine = true;
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
  bool drawCorPIDOXParticle = false;

  bool drawCorEKE = false;
  bool drawYP = true;
  bool drawYPGrid = true;
  bool drawYPTTA0 = true;
  bool drawKT = false;
  bool drawPtoa = false;

  bool drawPtoaR21 = false;
  bool drawKeoaR21 = false;
  bool drawNZR21 = true;
  bool fitNZR21 = true;

  int nbinsPID = 800;
  double maxdEdx = 2000;
  double maxPoz = 3000;
  int nbinsY0 = 200;
  double minY0 = -1;
  double maxY0 = 2;
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

  ///////////////////////////////////////////////////////////////////////////////////////////////////

  //const char *wProbString = "prob";
  const char *wProbString = "calculate_prob(p_lab,dedx)";
  if (!anaRawPIDProjection)
    wProbString = "prob";

  vector<int> selSystemIdxR21;
  int selSysR21 = selSys;
  int selCombR21 = 0;
  if (drawPtoaR21) {
    selCombR21 = 0;
    selSystemIdxR21.push_back(fSysCombIdx[selCombR21][0]);
    selSystemIdxR21.push_back(fSysCombIdx[selCombR21][1]);
    selSysR21 = kall;
  }
  else
    for (auto sys : {0,1})
      selSystemIdxR21.push_back(sys);

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

  GetCutG(1,1,0,"data_bound/bound_pidRaw_f7_left_45_100_108_ttaRG0_20_yAll.root");
  GetCutG(1,1,1,"data_bound/bound_pidRaw_f7_left_45_100_108_ttaRG0_20_yAll.root");
  GetCutG(1,1,2,"data_bound/bound_pidRaw_f7_left_45_100_108_ttaRG20_40_yAll.root");
  GetCutG(1,1,3,"data_bound/bound_pidRaw_f7_left_45_100_108_ttaRG40_60_yAll.root");
  GetCutG(1,1,4,"data_bound/bound_pidRaw_f7_left_45_100_108_ttaRG40_60_yAll.root");

  GetCutG(1,2,0,"data_bound/bound_pidRaw_f7_right_45_100_108_ttaRG0_20_yAll.root");
  GetCutG(1,2,1,"data_bound/bound_pidRaw_f7_right_45_100_108_ttaRG0_20_yAll.root");
  GetCutG(1,2,2,"data_bound/bound_pidRaw_f7_right_45_100_108_ttaRG20_40_yAll.root");
  GetCutG(1,2,3,"data_bound/bound_pidRaw_f7_right_45_100_108_ttaRG40_60_yAll.root");
  GetCutG(1,2,4,"data_bound/bound_pidRaw_f7_right_45_100_108_ttaRG40_60_yAll.root");

  GetCutG(0,1,0,"data_bound/bound_pidRaw_f7_left_45_100_132_ttaRG0_20_yAll.root");
  GetCutG(0,1,1,"data_bound/bound_pidRaw_f7_left_45_100_132_ttaRG0_20_yAll.root");
  GetCutG(0,1,2,"data_bound/bound_pidRaw_f7_left_45_100_132_ttaRG20_40_yAll.root");
  GetCutG(0,1,3,"data_bound/bound_pidRaw_f7_left_45_100_132_ttaRG40_60_yAll.root");
  GetCutG(0,1,4,"data_bound/bound_pidRaw_f7_left_45_100_132_ttaRG40_60_yAll.root");

  GetCutG(0,2,0,"data_bound/bound_pidRaw_f7_right_45_100_132_ttaRG0_20_yAll.root");
  GetCutG(0,2,1,"data_bound/bound_pidRaw_f7_right_45_100_132_ttaRG0_20_yAll.root");
  GetCutG(0,2,2,"data_bound/bound_pidRaw_f7_right_45_100_132_ttaRG20_40_yAll.root");
  GetCutG(0,2,3,"data_bound/bound_pidRaw_f7_right_45_100_132_ttaRG40_60_yAll.root");
  GetCutG(0,2,4,"data_bound/bound_pidRaw_f7_right_45_100_132_ttaRG40_60_yAll.root");

  // left and right
  GetCutG(1,0,0,"data_bound/bound_pidRaw_f7_right_45_100_108_ttaRG0_20_yAll.root");
  GetCutG(1,0,1,"data_bound/bound_pidRaw_f7_right_45_100_108_ttaRG0_20_yAll.root");
  GetCutG(1,0,2,"data_bound/bound_pidRaw_f7_right_45_100_108_ttaRG20_40_yAll.root");
  GetCutG(1,0,3,"data_bound/bound_pidRaw_f7_right_45_100_108_ttaRG40_60_yAll.root");
  GetCutG(1,0,4,"data_bound/bound_pidRaw_f7_right_45_100_108_ttaRG40_60_yAll.root");
  GetCutG(0,0,0,"data_bound/bound_pidRaw_f7_right_45_100_132_ttaRG0_20_yAll.root");
  GetCutG(0,0,1,"data_bound/bound_pidRaw_f7_right_45_100_132_ttaRG0_20_yAll.root");
  GetCutG(0,0,2,"data_bound/bound_pidRaw_f7_right_45_100_132_ttaRG20_40_yAll.root");
  GetCutG(0,0,3,"data_bound/bound_pidRaw_f7_right_45_100_132_ttaRG40_60_yAll.root");
  GetCutG(0,0,4,"data_bound/bound_pidRaw_f7_right_45_100_132_ttaRG40_60_yAll.root");


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
        TH1D *histPtoaArray[fNumSystems][fNumCutTTAs][fNumParticles] = {0};
        TH1D *histKeoaArray[fNumSystems][fNumCutTTAs][fNumParticles] = {0};
        TH2D *histYPR21Array[fNumSystems][fNumCutTTAs][fNumParticles] = {0};

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

            auto iCutY0 = 0;
            //for (auto iCutY0 : fCutY0Idx)
            {
              //if (selCutY0>=0 && selCutY0!=iCutY0) continue;
              TCut cutY0 = fCutY0Values[iCutY0];

              for (auto iCutTTA : selCutTTAIdx)
              {
                if (selCutTTA>=0 && selCutTTA!=iCutTTA) continue;
                TCut cutTTA = fCutTTAValues[iCutTTA];

                auto namePIDRaw = makeName("pidRaw",iAna,iLR,iMult,iSys,iCutTTA,iCutY0);
                auto titlePIDRaw = makeTitle("Raw",iAna,iLR,iMult,iSys,iCutTTA,iCutY0);

                auto histPID = new TH2D(namePIDRaw,Form("%s;p/Z (MeV/c);dE/dx;",titlePIDRaw),nbinsPID,0,maxPoz,nbinsPID,0,maxdEdx);
                //histPID -> SetMinimum(0.5);
                //histPID -> SetMaximum(800);
                project(treeAll,namePIDRaw,"dedx:p_lab",cutY0*cutTTA);

                TCanvas *cvsPIDRaw = nullptr;
                if (drawRawPID) {
                  cvsPIDRaw = makeCvs2(namePIDRaw,1000,700);
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

                    auto nameProj = makeName(Form("proj0_dedx_%d",iProj),iAna,iLR,iMult,iSys,iCutTTA,iCutY0);
                    if (drawRawPIDSub && fitdEdx) nameProj = makeName(Form("proj1_dedx_%d",iProj),iAna,iLR,iMult,iSys,iCutTTA,iCutY0);
                    else if (fitdEdx)             nameProj = makeName(Form("proj2_dedx_%d",iProj),iAna,iLR,iMult,iSys,iCutTTA,iCutY0);
                    else                          nameProj = makeName(Form("proj3_dedx_%d",iProj),iAna,iLR,iMult,iSys,iCutTTA,iCutY0);

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
                      const char *nameProjFit = makeName(Form("projFit_dedx_%d",iProj),iAna,iLR,iMult,iSys,iCutTTA,iCutY0);
                      auto f1dEdxTotal = new TF1(nameProjFit,"gaus(0)+gaus(3)+gaus(6)+gaus(9)+gaus(12)",0,maxdEdx);
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

                        const char *nameProjFitPart = makeName(Form("projFitdEdx0_%d",iProj),iAna,iLR,iMult,iSys,iCutTTA,iCutY0,iParticle);
                        auto f1dEdxParticle = new TF1(nameProjFitPart,"gaus(0)",0,maxdEdx);
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
                        const char *nameProjFitPart = makeName(Form("projFitdEdx_%d",iProj),iAna,iLR,iMult,iSys,iCutTTA,iCutY0,iParticle);
                        auto f1dEdxParticle = new TF1(nameProjFit,"gaus(0)",0,maxdEdx);
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

                      if (drawRawPIDProjection)
                      {
                        cvsProj -> cd();
                        cvsProj -> SaveAs(pathToFigures+nameProj+figureFormat); 
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
            for (auto iCutY0 : fCutY0Idx)
            {
              if (selCutY0>=0 && selCutY0!=iCutY0) continue;
              TCut cutY0 = fCutY0Values[iCutY0];

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
                  auto nameKTCorPart = makeName("ktCor",iAna,iLR,iMult,iSys,iCutTTA,iCutY0,iParticle);
                  auto titleKTCorPart = makeTitle("Cor.",iAna,iLR,iMult,iSys,iCutTTA,iCutY0,iParticle);

                  auto histKTCorPart = new TH2D(nameKTCorPart,Form("%s;KE_{Lab}/A (MeV);#theta_{lab};",titleKTCorPart),nbinsKeoa,0,keoaMax,nbinsTheta,0,thetaMax);
                  giSys = iSys;
                  giParticle = iParticle;
                  TCut selection = cut0 * TCut(Form("%s/eff/%d",wProbString,numEventsInAna[iSys])) * cutY0 * cutTTA;
                  auto partz = fParticleZ[iParticle];
                  auto partm = fParticleMass[iParticle];
                  auto parta = fParticleA[iParticle];
                  const char *expression = Form("theta_lab*TMath::RadToDeg():(sqrt((p_lab*%d)*(p_lab*%d)+%f*%f)-%f)/%d",partz,partz,partm,partm,partm,parta);
                  project(treeParticle,nameKTCorPart,expression,selection);

                  if (cvsKT==nullptr) {
                    cvsKT = makeCvs2(nameKTCorPart,1000,700);
                    cvsKT -> Divide(3,2);
                  }

                  cvsKT -> cd(iParticle+1);
                  histKTCorPart -> Draw("colz");
                }

                if (drawYP)
                {
                  auto nameYPPart = makeName("yp",iAna,iLR,iMult,iSys,iCutTTA,iCutY0,iParticle);
                  auto titleYPPart = makeTitle("Cor.",iAna,iLR,iMult,iSys,iCutTTA,iCutY0,iParticle);

                  auto histYPPart = new TH2D(nameYPPart,Form("%s;y_{0};p_{T}/A;",titleYPPart),nbinsY0,minY0,maxY0,nbinsPt,0,maxPt);
                  giSys = iSys;
                  giParticle = iParticle;
                  TCut selection = cut0 * TCut(Form("%s/eff/%d",wProbString,numEventsInAna[iSys])) * cutY0 * cutTTA;
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
                    for (auto ibinX=1; ibinX<=7; ++ibinX) {
                      auto binX = ibinX+2;
                      auto line = new TLine(binnx.lowEdge(binX),0,binnx.lowEdge(binX),400);
                      line -> Draw("samel");
                    }
                    //for (auto ibinY=1; ibinY<=8; ++ibinY) {
                    for (auto ibinY=1; ibinY<=9; ++ibinY) {
                      auto binY = ibinY;
                      auto line = new TLine(-.5,binny.lowEdge(binY),1,binny.lowEdge(binY));
                      line -> Draw("samel");
                    }

                    for (auto ibinX=1; ibinX<=6; ++ibinX) {
                      for (auto ibinY=1; ibinY<=8; ++ibinY) {
                        auto binX = ibinX+2;
                        auto binY = ibinY;
                        auto x = binnx.getCenter(binX);
                        auto y = binny.getCenter(binY);
                        auto text = new TText(x,y,Form("%d,%d",binX,binY));
                        text -> Draw();
                        text -> SetTextAlign(22);
                        text -> SetTextSize(0.03);
                      }
                    }
                  }

                  if (drawYPTTA0) {
                    double par0 = 3, par1 = -20;
                    //for (auto iTTA0 : {5,6,7,8})
                    for (auto iTTA0 : {4})
                    //for (auto iTTA0 : {1,2,3,4,5,6,7,8})
                    {
                      int theta0 = 10*iTTA0;
                      auto nameYPPart2 = makeName(Form("yptta%d",theta0),iAna,iLR,iMult,iSys,iCutTTA,iCutY0,iParticle);

                      auto histYPPart = new TH2D(nameYPPart2,Form("%s;y_{0};p_{T}/A;",titleYPPart),nbinsY0,minY0,maxY0,nbinsPt,0,maxPt);
                      giSys = iSys;
                      giParticle = iParticle;
                      selection = cut0 * TCut(Form("%s/eff/%d",wProbString,numEventsInAna[iSys])) * cutY0 * cutTTA * fCutTTA0Values[iTTA0];
                      project(treeParticle,nameYPPart2,Form("pt_cm/%d:fy_cm/(by_cm/2)",fParticleA[iParticle]),selection);

                      auto nameYPPart3 = makeName(Form("fityptta%d",theta0),iAna,iLR,iMult,iSys,iCutTTA,iCutY0,iParticle);
                      TF1 *fitPYTT0 = new TF1(nameYPPart3,"[0]*(x+1)*(x-[1])",minY0,maxY0);
                      fitPYTT0 -> SetParameter(par0, par1);
                      fitPYTT0 -> SetParLimits(0,0,1500);
                      fitPYTT0 -> SetParLimits(1,-30,0);
                      histYPPart -> Fit(fitPYTT0,"RQ0");
                      par0 = fitPYTT0 -> GetParameter(0);
                      par1 = fitPYTT0 -> GetParameter(1);
                      fitPYTT0 -> Draw("samel");
                    }
                  }

                }
              }
            }
          }

          if (drawKT || drawCorEKE || drawCorPID || drawPtoaR21 || drawKeoaR21 || drawNZR21)
          {
            for (auto iCutY0 : fCutY0Idx)
            {
              if (selCutY0>=0 && selCutY0!=iCutY0) continue;
              TCut cutY0 = fCutY0Values[iCutY0];

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
                    auto namePIDCorPart = makeName("pidCor",iAna,iLR,iMult,iSys,iCutTTA,iCutY0,iParticle);
                    auto titlePIDCorPart = makeTitle("Cor.",iAna,iLR,iMult,iSys,iCutTTA,iCutY0,iParticle);

                    auto histPIDCorPart = new TH2D(namePIDCorPart,Form("%s;p/Z (MeV/c);dE/dx;",titlePIDCorPart),nbinsPID,0,maxPoz,nbinsPID,0,maxdEdx);
                    histPIDCorPart -> SetMinimum(0.5);
                    histPIDCorPart -> SetMaximum(800);
                    giSys = iSys;
                    giParticle = iParticle;
                    TCut selection = cut0 * cutY0 * cutTTA * TCut(Form("%s/eff/%d",wProbString,numEventsInAna[iSys]));
                    project(treeParticle,namePIDCorPart,"dedx:p_lab",selection);

                    if (histPIDCor==nullptr) {
                      auto namePIDCor = makeName("pidCor",iAna,iLR,iMult,iSys,iCutTTA,iCutY0);
                      histPIDCor = (TH2D *) histPIDCorPart -> Clone(namePIDCor);
                      histPIDCor -> SetTitle(Form("%s;p/Z (MeV/c);dE/dx;",titlePIDCorPart));
                    }
                    else
                      histPIDCor -> Add(histPIDCorPart);
                  }

                  if (drawCorEKE)
                  {
                    auto nameEKECorPart = makeName("ekeCor",iAna,iLR,iMult,iSys,iCutTTA,iCutY0,iParticle);
                    auto titleEKECorPart = makeTitle("Cor.",iAna,iLR,iMult,iSys,iCutTTA,iCutY0,iParticle);

                    auto histEKECorPart = new TH2D(nameEKECorPart,Form("%s;KE_{Lab} (MeV);dE/dx;",titleEKECorPart),nbinsKeoa,0,4*keoaMax,nbinsPID,0,.5*maxdEdx);
                    giSys = iSys;
                    giParticle = iParticle;
                    TCut selection = cut0 * TCut(Form("%s/eff/%d",wProbString,numEventsInAna[iSys])) * cutY0 * cutTTA;
                    auto partz = fParticleZ[iParticle];
                    auto partm = fParticleMass[iParticle];
                    auto parta = fParticleA[iParticle];
                    //const char *expression = Form("dedx:(sqrt((p_lab*%d)*(p_lab*%d)+%f*%f)-%f)/%d",partz,partz,partm,partm,partm,parta);
                    const char *expression = Form("dedx:(sqrt((p_lab*%d)*(p_lab*%d)+%f*%f)-%f)",partz,partz,partm,partm,partm);
                    project(treeParticle,nameEKECorPart,expression,selection);

                    if (histEKECor==nullptr) {
                      auto nameEKECor = makeName("ekeCor",iAna,iLR,iMult,iSys,iCutTTA,iCutY0);
                      histEKECor = (TH2D *) histEKECorPart -> Clone(nameEKECor);
                      histEKECor -> SetTitle(Form("%s;KE_{Lab} (MeV);dE/dx;",titleEKECorPart));
                    }
                    else
                      histEKECor -> Add(histEKECorPart);

                    histEKECorArray.push_back(histEKECorPart);
                  }

                  if (drawPtoaR21)
                  {
                    auto namePtoaPart = makeName("ptoa",iAna,iLR,iMult,iSys,iCutTTA,iCutY0,iParticle);
                    auto titlePtoaPart = makeTitle("Cor.",iAna,iLR,iMult,iSys,iCutTTA,iCutY0,iParticle);

                    auto histPtoa = new TH1D(namePtoaPart,Form("%s;p_{T}/A (MeV/c);",titlePtoaPart),nbinsPtoa,0,ptoaMax);
                    giSys = iSys;
                    giParticle = iParticle;
                    TCut selection = cut0 * TCut(Form("%s/eff/%d",wProbString,numEventsInAna[iSys])) * cutY0 * cutTTA ;

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
                      selection = cut0 * TCut(Form("%s/eff/%d","1",numEventsInAna[iSys])) * selectionBound * cutY0 * cutTTA;
                    }

                    project(treeParticle,namePtoaPart,Form("pt_cm/%d",fParticleA[iParticle]),selection);

                    histPtoaArray[iSys][iCutTTA][iParticle] = histPtoa;
                  }

                  if (drawKeoaR21)
                  {
                    auto nameKeoaPart = makeName("keoa",iAna,iLR,iMult,iSys,iCutTTA,iCutY0,iParticle);
                    auto titleKeoaPart = makeTitle("Cor.",iAna,iLR,iMult,iSys,iCutTTA,iCutY0,iParticle);

                    auto histKeoa = new TH1D(nameKeoaPart,Form("%s;KE_{Lab}/A (MeV);",titleKeoaPart),nbinsKeoaR21,0,keoaMax);
                    giSys = iSys;
                    giParticle = iParticle;
                    TCut selection = cut0 * TCut(Form("%s/eff/%d",wProbString,numEventsInAna[iSys])) * cutY0 * cutTTA;

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
                      selection = cut0 * TCut(Form("%s/eff/%d","1",numEventsInAna[iSys])) * selectionBound * cutY0 * cutTTA;
                    }

                    auto partz = fParticleZ[iParticle];
                    auto partm = fParticleMass[iParticle];
                    auto parta = fParticleA[iParticle];
                    const char *expression = Form("(sqrt((p_lab*%d)*(p_lab*%d)+%f*%f)-%f)/%d",partz,partz,partm,partm,partm,parta);
                    project(treeParticle,nameKeoaPart,expression,selection);

                    histKeoaArray[iSys][iCutTTA][iParticle] = histKeoa;
                  }

                  if (drawNZR21)
                  {
                    auto nameYPR21Part = makeName("ypR21",iAna,iLR,iMult,iSys,iCutTTA,iCutY0,iParticle);
                    auto titleYPR21Part = makeTitle("Cor.",iAna,iLR,iMult,iSys,iCutTTA,iCutY0,iParticle);

                    auto histYPR21Part = new TH2D(nameYPR21Part,Form("%s;y_{0};p_{T}/A;",titleYPR21Part),nbinsNZYP,-1,1,nbinsNZYP,0,400);
                    giSys = iSys;
                    giParticle = iParticle;
                    TCut selection = cut0 * TCut(Form("%s/eff/%d",wProbString,numEventsInAna[iSys])) * cutY0 * cutTTA;
                    project(treeParticle,nameYPR21Part,Form("pt_cm/%d:fy_cm/(by_cm/2)",fParticleA[iParticle]),selection);

                    histYPR21Array[iSys][iCutTTA][iParticle] = histYPR21Part;
                  }
                }

                if (drawCorPID)
                {
                  auto namePIDCor = makeName("pidCor",iAna,iLR,iMult,iSys,iCutTTA,iCutY0);
                  auto cvsPIDCor = makeCvs2(namePIDCor,1000,700);
                  cvsPIDCor -> SetLogz();
                  histPIDCor -> SetMaximum(0.08);
                  histPIDCor -> SetMinimum(0.0001);
                  histPIDCor -> Draw("colz");
                  if (drawGuideLine)
                    for (auto iParticle : fParticleIdx)
                      graphPIDMean[iSys][iParticle] -> Draw("samel");
                }

                if (drawCorEKE)
                {
                  auto nameEKECor = makeName("ekeCor",iAna,iLR,iMult,iSys,iCutTTA,iCutY0);
                  auto cvsEKECor = makeCvs2(nameEKECor,1500,1000);
                  cvsEKECor -> Divide(3,2);

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
                  const char *namePtoa = makeName("ptoaAll",iAna,iLR,iMult,iSys,iCutTTA,iCutY0);
                  auto cvsR21 = makeCvs(namePtoa);

                  double histMax = 0;
                  for (auto iParticle : fParticleIdx) {
                    auto max0 = histPtoaArray[iSys][iCutTTA][iParticle] -> GetMaximum();
                    if (histMax < max0)
                      histMax = max0;
                  }
                  for (auto iParticle : fParticleIdx) {
                    auto hist = histPtoaArray[iSys][iCutTTA][iParticle];
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


        if (drawPtoaR21)
        {
          for (auto iCutY0 : fCutY0Idx)
          {
            if (selCutY0>=0 && selCutY0!=iCutY0) continue;

            for (auto iCutTTA : selCutTTAIdx)
            {
              if (selCutTTA>=0 && selCutTTA!=iCutTTA) continue;

              for (auto iComb : fSysCombIndx)
              {
                if (selCombR21>=0 && selCombR21!=iComb) continue;

                auto iSys2 = fSysCombIdx[selCombR21][0];
                auto iSys1 = fSysCombIdx[selCombR21][1];
                auto iSysComb = 100 + iSys2*10 + iSys1;
                
                const char *nameR21 = makeName("r21_ptoa",iAna,iLR,iMult,iSysComb,iCutTTA,iCutY0);
                auto titleR21 = makeTitle("Cor.",iAna,iLR,iMult,iSysComb,iCutTTA,iCutY0);

                auto cvs = makeCvs(nameR21,680,550);

                auto hist = new TH2D(nameR21,Form("%s;p_{T}/A (MeV/c);R_{21}",titleR21),100,0,400,100,0,2.0);
                hist -> Draw();

                auto legend = new TLegend();
                for (auto iParticle : fParticleIdx)
                {
                  const char *nameParticle = fParticleNames[iParticle];

                  auto hist1 = histPtoaArray[iSys2][iCutTTA][iParticle];
                  auto hist2 = histPtoaArray[iSys1][iCutTTA][iParticle];

                  const char *nameR21Part = Form("%s_%s",nameR21,nameParticle);

                  auto hist0 = (TH1D *) hist1 -> Clone(nameR21Part);
                  hist0 -> Divide(hist2);
                  auto graph = new TGraph();
                  for (auto bin=1; bin<=8; ++bin)
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
          for (auto iCutY0 : fCutY0Idx)
          {
            if (selCutY0>=0 && selCutY0!=iCutY0) continue;

            for (auto iCutTTA : selCutTTAIdx)
            {
              if (selCutTTA>=0 && selCutTTA!=iCutTTA) continue;

              for (auto iComb : fSysCombIndx)
              {
                if (selCombR21>=0 && selCombR21!=iComb) continue;

                auto iSys2 = fSysCombIdx[selCombR21][0];
                auto iSys1 = fSysCombIdx[selCombR21][1];
                auto iSysComb = 100 + iSys2*10 + iSys1;

                const char *nameR21 = makeName("r21_keoa",iAna,iLR,iMult,iSysComb,iCutTTA,iCutY0);
                auto titleR21 = makeTitle("Cor.",iAna,iLR,iMult,iSysComb,iCutTTA,iCutY0);

                auto cvs = makeCvs(nameR21,680,550);
                //cvs -> SetGrid(1,1);

                auto hist = new TH2D(nameR21,Form("%s;KE_{Lab}/A (MeV);R_{21}",titleR21),100,0,400,100,0,2.0);
                hist -> Draw();

                auto legend = new TLegend();
                for (auto iParticle : fParticleIdx)
                {
                  const char *nameParticle = fParticleNames[iParticle];

                  auto hist1 = histKeoaArray[iSys2][iCutTTA][iParticle];
                  auto hist2 = histKeoaArray[iSys1][iCutTTA][iParticle];

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
          for (auto iCutY0 : fCutY0Idx)
          {
            if (selCutY0>=0 && selCutY0!=iCutY0) continue;

            for (auto iCutTTA : selCutTTAIdx)
            {
              if (selCutTTA>=0 && selCutTTA!=iCutTTA) continue;

              for (auto iComb : fSysCombIndx)
              {
                if (selCombR21>=0 && selCombR21!=iComb) continue;

                auto iSys2 = fSysCombIdx[selCombR21][0];
                auto iSys1 = fSysCombIdx[selCombR21][1];
                auto iSysComb = 100 + iSys2*10 + iSys1;

                const char *nameAlpha = makeName("alpha",iAna,iLR,iMult,iSysComb,iCutTTA,iCutY0);
                const char *nameBeta = makeName("beta",iAna,iLR,iMult,iSysComb,iCutTTA,iCutY0);

                auto cvsAlpha = makeCvs(nameAlpha,2500,2000);
                auto cvsBeta  = makeCvs(nameBeta,2500,2000);
                //cvsAlpha -> Divide(6,8);
                //cvsBeta  -> Divide(6,8);
                Divide0(cvsAlpha,6,8);// -> Divide(6,8);
                Divide0(cvsBeta,6,8);//  -> Divide(6,8);

                double alphaArray[10][10];
                double betaArray[10][10];

                auto binnx = binning(histYPR21Array[iSys2][iCutTTA][0],1);
                auto binny = binning(histYPR21Array[iSys2][iCutTTA][0],2);

                for (auto ibinX=1; ibinX<=6; ++ibinX)
                {
                  for (auto ibinY=1; ibinY<=8; ++ibinY)
                  {
                    auto binX = ibinX+2;
                    auto binY = ibinY;

                    auto graph_z1_inN = new TGraph(); graph_z1_inN -> SetMarkerSize(1.5); graph_z1_inN -> SetMarkerStyle(20); graph_z1_inN -> SetMarkerColor(kBlack);
                    auto graph_z2_inN = new TGraph(); graph_z2_inN -> SetMarkerSize(1.5); graph_z2_inN -> SetMarkerStyle(25); graph_z2_inN -> SetMarkerColor(kRed  );
                    auto graph_n0_inZ = new TGraph(); graph_n0_inZ -> SetMarkerSize(1.5); graph_n0_inZ -> SetMarkerStyle(20); graph_n0_inZ -> SetMarkerColor(kBlack);
                    auto graph_n1_inZ = new TGraph(); graph_n1_inZ -> SetMarkerSize(1.5); graph_n1_inZ -> SetMarkerStyle(25); graph_n1_inZ -> SetMarkerColor(kRed  );
                    auto graph_n2_inZ = new TGraph(); graph_n2_inZ -> SetMarkerSize(1.5); graph_n2_inZ -> SetMarkerStyle(22); graph_n2_inZ -> SetMarkerColor(kBlue );

                    for (auto iParticle : fParticleIdx)
                    {
                      auto value2 = histYPR21Array[iSys2][iCutTTA][iParticle] -> GetBinContent(binX, binY);
                      auto value1 = histYPR21Array[iSys1][iCutTTA][iParticle] -> GetBinContent(binX, binY);
                      auto r21Value = value2/value1;

                      auto zValue = fParticleZ[iParticle];
                      auto nValue = fParticleN[iParticle];

                      if (zValue==1) graph_z1_inN -> SetPoint(graph_z1_inN -> GetN(), nValue, r21Value);
                      if (zValue==2) graph_z2_inN -> SetPoint(graph_z2_inN -> GetN(), nValue, r21Value);
                      if (nValue==0) graph_n0_inZ -> SetPoint(graph_n0_inZ -> GetN(), zValue, r21Value);
                      if (nValue==1) graph_n1_inZ -> SetPoint(graph_n1_inZ -> GetN(), zValue, r21Value);
                      if (nValue==2) graph_n2_inZ -> SetPoint(graph_n2_inZ -> GetN(), zValue, r21Value);
                    }

                    auto idxCvs = (8-ibinY)*6+ibinX;

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
                    graph_z1_inN -> Draw("samep");
                    graph_z2_inN -> Draw("samep");
                    TLegend *legendAlpha = new TLegend();
                    legendAlpha -> AddEntry(graph_z2_inN, Form("(%d,%d)",binX,binY),"");

                    if (fitNZR21) {
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
                      fit1 -> Draw("samel");
                      fit2 -> Draw("samel");
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
                      alphaArray[ibinX][ibinY] = fit2->GetParameter(1);
                    }
                    makeLegend(cvs,legendAlpha,"lt",-.1,0,.5,.42) -> Draw();

                    cvs = cvsBeta -> cd(idxCvs);
                    cvs -> SetLogy();
                    auto frameBeta = new TH2D(Form("%s_%d",cvsBeta->GetName(),idxCvs),";Z;R21",100,0,3,100,0.5,2.0);
                    frameBeta -> Draw();
                    frameBeta -> GetXaxis() -> SetTitleSize(0.10);
                    frameBeta -> GetYaxis() -> SetTitleSize(0.10);
                    frameBeta -> GetXaxis() -> SetLabelSize(0.10);
                    frameBeta -> GetYaxis() -> SetLabelSize(0.10);
                    frameBeta -> GetXaxis() -> SetNdivisions(506);
                    frameBeta -> GetYaxis() -> SetNdivisions(506);
                    frameBeta -> GetXaxis() -> CenterTitle();
                    frameBeta -> GetYaxis() -> CenterTitle();
                    graph_n0_inZ -> Draw("samep");
                    graph_n1_inZ -> Draw("samep");
                    graph_n2_inZ -> Draw("samep");
                    TLegend *legendBeta = new TLegend();
                    legendBeta -> AddEntry(graph_n2_inZ , Form("(%d,%d)",binX,binY),"");

                    if (fitNZR21) {
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
                      fit3 -> Draw("samel");
                      fit4 -> Draw("samel");
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
                      legendBeta -> AddEntry(fit4, Form("#beta = %.3f",fit4->GetParameter(1)), "");
                      legendBeta -> AddEntry(fit4, Form("rms = %.3f",rms),"");
                      betaArray[ibinX][ibinY] = fit4->GetParameter(1);
                    }
                    makeLegend(cvs,legendBeta,"rt",0.11,0,.5,.42) -> Draw();

                  }
                }

                auto nameR21YP = makeName("r21_yp",iAna,iLR,iMult,iSysComb,iCutTTA,iCutY0);
                const char *nameAlphaBetaOut = Form("data_ab/alpha_beta_%s.txt",nameR21YP);
                cout << nameAlphaBetaOut << endl;
                std::ofstream fileR21YP(nameAlphaBetaOut);

                auto nameR21AB = makeName("r21_ab",iAna,iLR,iMult,iSysComb,iCutTTA,iCutY0);
                auto cvsAB = makeCvs(nameR21AB);
                (new TH2D(nameR21AB,";#alpha;#beta",100,0,.6,100,-.5,.1)) -> Draw();

                auto legendAB = new TLegend();
                for (auto ibinY=8; ibinY>=1; --ibinY) {
                  auto graphAB = new TGraph();
                  for (auto ibinX=6; ibinX>=1; --ibinX) {
                    auto alpha = alphaArray[ibinX][ibinY];
                    auto beta = betaArray[ibinX][ibinY];
                    graphAB -> SetPoint(graphAB->GetN(),alpha,beta);
                  }
                  graphAB -> SetLineColor(ibinY);
                  graphAB -> SetMarkerColor(ibinY);
                  graphAB -> SetMarkerStyle(24);
                  graphAB -> Draw("samep");
                  legendAB -> AddEntry(graphAB,Form("%d) p_{T}/A = %.0f ~ %.0f MeV/c",ibinY, binny.lowEdge(ibinY),binny.highEdge(ibinY)));
                }
                makeLegend(cvsAB,legendAB) -> Draw();
                auto invf = new TF1(Form("invf_%s",nameR21AB),"-x",0,.5);
                invf -> SetLineStyle(2);
                invf -> SetLineColor(kGray);
                invf -> Draw("samel");

                fileR21YP << "alpha" << endl;
                cout << "alpha" << endl;
                for (auto ibinX=6; ibinX>=1; --ibinX) {
                  for (auto ibinY=8; ibinY>=1; --ibinY) {
                    fileR21YP << alphaArray[ibinX][ibinY] << ", ";
                    cout << alphaArray[ibinX][ibinY] << ", ";
                  }
                  fileR21YP << endl;
                  cout << endl;
                }

                fileR21YP << "beta" << endl;
                cout << "beta" << endl;
                for (auto ibinX=6; ibinX>=1; --ibinX) {
                  for (auto ibinY=8; ibinY>=1; --ibinY) {
                    fileR21YP << betaArray[ibinX][ibinY] << ", ";
                    cout << betaArray[ibinX][ibinY] << ", ";
                  }
                  fileR21YP << endl;
                  cout << endl;
                }

              }
            }
          }
        }



      }
    }
  }

  if (saveFigures)
    for (auto cvs : cvsArray) {
      cvs -> cd();
      cvs -> SaveAs(pathToFigures+cvs->GetName()+figureFormat); 
    }
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

void Divide0(TCanvas *cvs, int nx, int ny)
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
