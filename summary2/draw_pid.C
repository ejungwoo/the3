#include "KBGlobal.hh"
#include "init_variables.h"
#include "/Users/ejungwoo/config/ejungwoo.h"

int cvsXOff = 1300;
int giSys = 0;
int giParticle = 0;
vector<TCanvas *> cvsArray;

TGraph *graphFitPIDMean[fNumSystems][fNumParticles] = {0};
TGraph *graphFitPIDAmp[fNumSystems][fNumParticles] = {0};
TH1F *histPIDMeta[fNumSystems][fNumParticles] = {0};
TF1 *f1PIDMean[fNumSystems][fNumParticles] = {0};
TF1 *f1PIDSigma[fNumSystems][fNumParticles] = {0};
TCutG *cutgPID[fNumSystems][fNumParticles] = {0};
TGraph *graphPIDMean[fNumSystems][fNumParticles] = {0};
TGraph *graphPIDRange[fNumSystems][fNumParticles][2] = {0};

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
  cout_info << tree -> GetName() << " " << name << " " << expr << " " << selection << endl;
  tree -> Project(name,expr,selection);
}

int countCvs = 0;

TCanvas *newCvs(const char *name, int w=680, int h=550) {
  //const char *name0 = Form("cvs%d_%s",countCvs,name);
  const char *name0 = Form("cvs_%s",name);
  auto cvs = new TCanvas(name0,name0,cvsXOff+20*(countCvs+1), 20*(countCvs+1), w, h);
  cvs -> SetRightMargin(0.155);
  countCvs++;
  cvsArray.push_back(cvs);
  return cvs;
}

const char *makeName(const char *mainName, int idxAna, int idxLR, int idxMult, int idxSys, int idxTTA, int iParticle=-1)
{
  if (idxSys>=100) {
    idxSys = idxSys - 100;
    int idxSys1 = int(idxSys/10);
    int idxSys2 = idxSys - idxSys1;

    const char *name = Form("%s_%s_%s_%s_%so%s_%s",mainName,fAnaNames[idxAna],fLRNames[idxLR],fMultNames[idxMult],fSystemNames[idxSys1],fSystemNames[idxSys2],fCutTTANames[idxTTA]);
  }

  if (iParticle<0) {
    const char *name = Form("%s_%s_%s_%s_%s_%s",mainName,fAnaNames[idxAna],fLRNames[idxLR],fMultNames[idxMult],fSystemNames[idxSys],fCutTTANames[idxTTA]);
    return name;
  }
  const char *name = Form("%s_%s_%s_%s_%s_%s_%s",mainName,fAnaNames[idxAna],fLRNames[idxLR],fMultNames[idxMult],fSystemNames[idxSys],fCutTTANames[idxTTA],fParticleNames[iParticle]);
  return name;
}

const char *makeTitle(const char *mainName, int idxAna, int idxLR, int idxMult, int idxSys, int idxTTA, int iParticle=5)
{
  if (idxSys>=100) {
    idxSys = idxSys - 100;
    int idxSys1 = int(idxSys/10);
    int idxSys2 = idxSys - idxSys1;

    const char *title = Form("%s, %s, %s, %s, %s/%s, %s, %s",mainName,fAnaTitles[idxAna],fLRTitles[idxLR],fMultTitles[idxMult],fSystemTitles[idxSys1],fSystemTitles[idxSys2],fCutTTATitles[idxTTA],fParticleTitles[iParticle]);
    return title;
  }

  const char *title = Form("%s, %s, %s, %s, %s, %s, %s",mainName,fAnaTitles[idxAna],fLRTitles[idxLR],fMultTitles[idxMult],fSystemTitles[idxSys],fCutTTATitles[idxTTA],fParticleTitles[iParticle]);
  return title;
}

void draw_pid()
{
  int selAna = kf7;
  int selLR = klr;
  int selMult = kn55;
  int selSys = kall;

  //const int cutTTAIdx[] = {0,1,2,3,4};
  //const int cutTTAIdx[] = {3,4};
  const int cutTTAIdx[] = {5,6};
  int selTTA = kall;

  bool saveFigures = true;
  TString pathToFigures = "figures2/";
  TString figureFormat = ".png";

  //TCut cut0 = "fy_cm/(by_cm/2)>0 && fy_cm/(by_cm/2)<.4";
  TCut cut0 = "fy_cm/(by_cm/2)>0";
  //TCut cut0 = "1";

  ///////////////////////////////////////////////////////////////////////////////////////////////////

  bool drawRawPID = false;
  bool drawGuideLine = false;
  bool drawHandCut = false;
  int nbinsPID = 500;
  double maxdEdx = 2000;
  double maxPoz = 3000;

  int nbinsY0 = 200;
  double minY0 = -1;
  double maxY0 = 2;
  int nbinsPt = 200;
  double maxPt = 800;

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

  bool drawPY = false;
  bool drawPYTTA0 = true;

  bool drawKT = false;
  int nbinsKeoa = 100;
  double keoaMax = 400.;
  int nbinsTheta = 100;
  double thetaMax = 90;

  bool drawPtoaR21 = true;
  bool useHandCut = false;
  bool drawPtoa = false;
  int nbinsPtoa = 8;
  double ptoaMax = 400;

  bool drawKeoaR21 = false;
  int nbinsKeoaR21 = 4;

  ///////////////////////////////////////////////////////////////////////////////////////////////////

  //const char *wProbString = "prob";
  const char *wProbString = "calculate_prob(p_lab,dedx)";
  if (!anaRawPIDProjection)
    wProbString = "prob";

  vector<int> systemIdxR21;
  int selSysR21 = selSys;
  int selCombR21 = 0;
  if (drawPtoaR21) {
    selCombR21 = 0;
    systemIdxR21.push_back(fSysCombIdx[selCombR21][0]);
    systemIdxR21.push_back(fSysCombIdx[selCombR21][1]);
    selSysR21 = kall;
  }
  else
    for (auto sys : {0,1})
    //for (auto sys : {0})
      systemIdxR21.push_back(sys);

  if (drawPYTTA0)
  {
    nbinsY0 = 200;
    nbinsPt = 200;
  }

  ///////////////////////////////////////////////////////////////////////////////////////////////////

  auto sys132filegraph0 = new TFile("c0.root");
  auto sys132filegraph1 = new TFile("c1.root");
  auto sys132filegraph2 = new TFile("c2.root");
  auto sys132filegraph3 = new TFile("c3.root");

  auto sys132graph0 = (TGraph *) sys132filegraph0 -> Get("c0");
  auto sys132graph1 = (TGraph *) sys132filegraph1 -> Get("c1");
  auto sys132graph2 = (TGraph *) sys132filegraph2 -> Get("c2");
  auto sys132graph3 = (TGraph *) sys132filegraph3 -> Get("c3");

  auto sys132graph_cut1 = ejungwoo::make_cutg("cutg0_132",sys132graph0, sys132graph1); sys132graph_cut1 -> SetLineColor(kRed);
  auto sys132graph_cut2 = ejungwoo::make_cutg("cutg1_132",sys132graph1, sys132graph2); sys132graph_cut2 -> SetLineColor(kRed);
  auto sys132graph_cut3 = ejungwoo::make_cutg("cutg2_132",sys132graph2, sys132graph3); sys132graph_cut3 -> SetLineColor(kRed);

  auto sys108filegraph0 = new TFile("d0.root");
  auto sys108filegraph1 = new TFile("d1.root");
  auto sys108filegraph2 = new TFile("d2.root");
  auto sys108filegraph3 = new TFile("d3.root");

  auto sys108graph0 = (TGraph *) sys108filegraph0 -> Get("d0");
  auto sys108graph1 = (TGraph *) sys108filegraph1 -> Get("d1");
  auto sys108graph2 = (TGraph *) sys108filegraph2 -> Get("d2");
  auto sys108graph3 = (TGraph *) sys108filegraph3 -> Get("d3");

  auto sys108graph_cut1 = ejungwoo::make_cutg("cutg0_108",sys108graph0, sys108graph1); sys108graph_cut1 -> SetLineColor(kRed);
  auto sys108graph_cut2 = ejungwoo::make_cutg("cutg1_108",sys108graph1, sys108graph2); sys108graph_cut2 -> SetLineColor(kRed);
  auto sys108graph_cut3 = ejungwoo::make_cutg("cutg2_108",sys108graph2, sys108graph3); sys108graph_cut3 -> SetLineColor(kRed);

  for (auto graph_cut : {sys132graph_cut1,sys132graph_cut2,sys132graph_cut3,sys108graph_cut1,sys108graph_cut2,sys108graph_cut3}) {
    graph_cut -> SetVarX("p_lab");
    graph_cut -> SetVarY("dedx");
  }

  ///////////////////////////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////////////////////////

  gStyle -> SetOptStat(0);

  //for (auto iSys : fSystemIdx)
  for (auto iSys : systemIdxR21)
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
      graphPIDRange[iSys][iParticle][0] = graph2;
      graphPIDRange[iSys][iParticle][1] = graph3;
    }
  }

  for (auto iAna : fAnaIdx)
  {
    if (selAna>=0 && selAna!=iAna) continue;
    const char *anaFName = fAnaFNames[iAna];
    const char *spVersion = fAnaVersion[iAna];

    for (auto iLR : fLRIdx)
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

        //for (auto iSys : fSystemIdx)
        for (auto iSys : systemIdxR21)
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

            for (auto iCutTTA : cutTTAIdx)
            {
              if (selTTA>=0 && selTTA!=iCutTTA) continue;
              const char *ttaName = fCutTTANames[iCutTTA];
              const char *ttaTitle = fCutTTATitles[iCutTTA];

              TCut cutTTA = fCutTTAValues[iCutTTA];

              auto namePIDRaw = makeName("pidRaw",iAna,iLR,iMult,iSys,iCutTTA);
              auto titlePIDRaw = makeTitle("Raw",iAna,iLR,iMult,iSys,iCutTTA);

              auto histPID = new TH2D(namePIDRaw,Form("%s;p/Z (MeV/c);dE/dx;",titlePIDRaw),nbinsPID,0,maxPoz,nbinsPID,0,maxdEdx);
              //histPID -> SetMinimum(0.5);
              //histPID -> SetMaximum(800);
              //treeAll -> Project(namePIDRaw,"dedx:p_lab",cutTTA);
              project(treeAll,namePIDRaw,"dedx:p_lab",cutTTA);

              TCanvas *cvsPIDRaw = nullptr;
              if (drawRawPID) {
                cvsPIDRaw = newCvs(namePIDRaw,1000,700);
                cvsPIDRaw -> SetLogz();
                histPID -> Draw("colz");
                if (drawGuideLine)
                  for (auto iParticle : fParticleIdx)
                    graphPIDMean[iSys][iParticle] -> Draw("samel");
                if (drawHandCut) {
                  if (sys==108) {
                    sys108graph_cut1 -> Draw("samel");
                    sys108graph_cut2 -> Draw("samel");
                    sys108graph_cut3 -> Draw("samel");
                  }
                  if (sys==132) {
                    sys132graph_cut1 -> Draw("samel");
                    sys132graph_cut2 -> Draw("samel");
                    sys132graph_cut3 -> Draw("samel");
                  }
                }

              }

              if (anaRawPIDProjection)
              {
                auto binn = ebinning(histPID);
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

                  auto nameProj = makeName(Form("proj0_dedx_%d",iProj),iAna,iLR,iMult,iSys,iCutTTA);
                  if (drawRawPIDSub && fitdEdx) nameProj = makeName(Form("proj1_dedx_%d",iProj),iAna,iLR,iMult,iSys,iCutTTA);
                  else if (fitdEdx)             nameProj = makeName(Form("proj2_dedx_%d",iProj),iAna,iLR,iMult,iSys,iCutTTA);
                  else                          nameProj = makeName(Form("proj3_dedx_%d",iProj),iAna,iLR,iMult,iSys,iCutTTA);

                  auto histProj = (TH1D *) histPID -> ProjectionY(nameProj,bin1,bin2);
                  histProj -> Rebin(dbin);

                  TCanvas *cvsProj = nullptr;

                  if (pozC>pidProjRange1 && pozC<pidProjRange2)
                  {
                    //cout_info << countCvs << " " << poz1 << " " << poz2 << endl;
                    TString titleProj = TString("p/Z=")+poz1+"-"+poz2+";dE/dx;";
                    histProj -> SetTitle(titleProj);

                    if (drawRawPIDProjection) {
                      cvsProj = newCvs(nameProj,1000,550);
                      cvsProj -> SetLogy();
                      histProj -> Draw();
                    }

                    double scaleAmp = 1.;
                    const char *nameProjFit = makeName(Form("projFit_dedx_%d",iProj),iAna,iLR,iMult,iSys,iCutTTA);
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

                      const char *nameProjFitPart = makeName(Form("projFitdEdx0_%d",iProj),iAna,iLR,iMult,iSys,iCutTTA,iParticle);
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
                      const char *nameProjFitPart = makeName(Form("projFitdEdx_%d",iProj),iAna,iLR,iMult,iSys,iCutTTA,iParticle);
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

                  newCvs(Form("refit_%s",namePIDRaw) ,1000,550);
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

        for (auto iSys : systemIdxR21)
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

          if (drawKT || drawPY)
          {
            auto iCutTTA = 0;
            TCut cutTTA = fCutTTAValues[iCutTTA];

            TCanvas *cvsKT = nullptr;

            for (auto iParticle : fParticleIdx)
            {
              const char *nameParticle = fParticleNames[iParticle];
              const char *titleParticle = fParticleTitles[iParticle];
              auto nameFileParticle = Form("data2/%s/sys%d_%s_%s_%s.%s.ana.%s.root",anaFName,sys,anaFName,lrFName,multFName,spVersion,nameParticle);
              //cout_info << "File : " << nameFileParticle << endl;

              auto treeParticle = new TChain(nameParticle);
              treeParticle -> Add(nameFileParticle);


              //if (selTTA>=0 && selTTA!=iCutTTA) continue;
              const char *ttaName = fCutTTANames[iCutTTA];
              const char *ttaTitle = fCutTTATitles[iCutTTA];

              if (drawKT)
              {
                auto nameKTCorPart = makeName("ktCor",iAna,iLR,iMult,iSys,iCutTTA,iParticle);
                auto titleKTCorPart = makeTitle("Cor.",iAna,iLR,iMult,iSys,iCutTTA,iParticle);

                auto histKTCorPart = new TH2D(nameKTCorPart,Form("%s;KE_{Lab}/A (MeV);#theta_{lab};",titleKTCorPart),nbinsKeoa,0,keoaMax,nbinsTheta,0,thetaMax);
                giSys = iSys;
                giParticle = iParticle;
                TCut selection = cut0 * TCut(Form("%s/eff/%d",wProbString,numEventsInAna[iSys])) * cutTTA;
                auto partz = fParticleZ[iParticle];
                auto partm = fParticleMass[iParticle];
                auto parta = fParticleA[iParticle];
                const char *expression = Form("theta_lab*TMath::RadToDeg():(sqrt((p_lab*%d)*(p_lab*%d)+%f*%f)-%f)/%d",partz,partz,partm,partm,partm,parta);
                project(treeParticle,nameKTCorPart,expression,selection);

                if (cvsKT==nullptr) {
                  cvsKT = newCvs(nameKTCorPart,1000,700);
                  cvsKT -> Divide(3,2);
                }

                cvsKT -> cd(iParticle+1);
                histKTCorPart -> Draw("colz");
              }

              if (drawPY)
              {
                auto namePYPart = makeName("py",iAna,iLR,iMult,iSys,iCutTTA,iParticle);
                auto titlePYPart = makeTitle("Cor.",iAna,iLR,iMult,iSys,iCutTTA,iParticle);

                auto histPYPart = new TH2D(namePYPart,Form("%s;y_{0};p_{T}/A;",titlePYPart),nbinsY0,minY0,maxY0,nbinsPt,0,maxPt);
                giSys = iSys;
                giParticle = iParticle;
                TCut selection = cut0 * TCut(Form("%s/eff/%d",wProbString,numEventsInAna[iSys])) * cutTTA;
                //treeParticle -> Project(namePYPart,Form("pt_cm/%d:fy_cm/(by_cm/2)",fParticleA[iParticle]),selection);
                project(treeParticle,namePYPart,Form("pt_cm/%d:fy_cm/(by_cm/2)",fParticleA[iParticle]),selection);

                auto cvs = newCvs(namePYPart);
                histPYPart -> Draw("colz");
                TLine *line0 = new TLine(0,0,0,maxPt);
                line0 -> SetLineStyle(9);
                line0 -> Draw("samel");

                if (drawPYTTA0) {
                  double par0 = 3, par1 = -20;
                  //for (auto iTTA0 : {5,6,7,8})
                  for (auto iTTA0 : {1,2,3,4,5,6,7,8})
                  {
                    int theta0 = 10*iTTA0;
                    auto namePYPart2 = makeName(Form("pytta%d",theta0),iAna,iLR,iMult,iSys,iCutTTA,iParticle);

                    auto histPYPart = new TH2D(namePYPart2,Form("%s;y_{0};p_{T}/A;",titlePYPart),nbinsY0,minY0,maxY0,nbinsPt,0,maxPt);
                    giSys = iSys;
                    giParticle = iParticle;
                    selection = cut0 * TCut(Form("%s/eff/%d",wProbString,numEventsInAna[iSys])) * (cutTTA + fCutTTA0Values[iTTA0]);
                    //treeParticle -> Project(namePYPart2,Form("pt_cm/%d:fy_cm/(by_cm/2)",fParticleA[iParticle]),selection);
                    project(treeParticle,namePYPart2,Form("pt_cm/%d:fy_cm/(by_cm/2)",fParticleA[iParticle]),selection);
                    //histPYPart -> Draw("samecol");

                    auto namePYPart3 = makeName(Form("fitpytta%d",theta0),iAna,iLR,iMult,iSys,iCutTTA,iParticle);
                    TF1 *fitPYTT0 = new TF1(namePYPart3,"[0]*(x+1)*(x-[1])",minY0,maxY0);
                    fitPYTT0 -> SetParameter(par0, par1);
                    fitPYTT0 -> SetParLimits(0,0,1500);
                    fitPYTT0 -> SetParLimits(1,-30,0);
                    histPYPart -> Fit(fitPYTT0,"RQ0");
                    par0 = fitPYTT0 -> GetParameter(0);
                    par1 = fitPYTT0 -> GetParameter(1);
                    fitPYTT0 -> Draw("samel");
                  }
                }

              }
            }
          }

          if (drawKT || drawCorEKE || drawCorPID || drawPtoaR21 || drawKeoaR21)
          {
            for (auto iCutTTA : cutTTAIdx)
            //for (auto iCutTTA : {1,2,3,4})
            //for (auto iCutTTA : {0})
            {
              //if (selTTA>=0 && selTTA!=iCutTTA) continue;
              const char *ttaName = fCutTTANames[iCutTTA];
              const char *ttaTitle = fCutTTATitles[iCutTTA];

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
                  auto namePIDCorPart = makeName("pidCor",iAna,iLR,iMult,iSys,iCutTTA,iParticle);
                  auto titlePIDCorPart = makeTitle("Cor.",iAna,iLR,iMult,iSys,iCutTTA,iParticle);

                  auto histPIDCorPart = new TH2D(namePIDCorPart,Form("%s;p/Z (MeV/c);dE/dx;",titlePIDCorPart),nbinsPID,0,maxPoz,nbinsPID,0,maxdEdx);
                  histPIDCorPart -> SetMinimum(0.5);
                  histPIDCorPart -> SetMaximum(800);
                  giSys = iSys;
                  giParticle = iParticle;
                  TCut selection = cut0 * cutTTA * TCut(Form("%s/eff/%d",wProbString,numEventsInAna[iSys]));
                  //treeParticle -> Project(namePIDCorPart,"dedx:p_lab",selection);
                  project(treeParticle,namePIDCorPart,"dedx:p_lab",selection);

                  if (histPIDCor==nullptr) {
                    auto namePIDCor = makeName("pidCor",iAna,iLR,iMult,iSys,iCutTTA);
                    histPIDCor = (TH2D *) histPIDCorPart -> Clone(namePIDCor);
                    histPIDCor -> SetTitle(Form("%s;p/Z (MeV/c);dE/dx;",titlePIDCorPart));
                  }
                  else
                    histPIDCor -> Add(histPIDCorPart);
                }

                if (drawCorEKE)
                {
                  auto nameEKECorPart = makeName("ekeCor",iAna,iLR,iMult,iSys,iCutTTA,iParticle);
                  auto titleEKECorPart = makeTitle("Cor.",iAna,iLR,iMult,iSys,iCutTTA,iParticle);

                  auto histEKECorPart = new TH2D(nameEKECorPart,Form("%s;KE_{Lab} (MeV);dE/dx;",titleEKECorPart),nbinsKeoa,0,4*keoaMax,nbinsPID,0,.5*maxdEdx);
                  giSys = iSys;
                  giParticle = iParticle;
                  TCut selection = cut0 * TCut(Form("%s/eff/%d",wProbString,numEventsInAna[iSys])) * cutTTA;
                  auto partz = fParticleZ[iParticle];
                  auto partm = fParticleMass[iParticle];
                  auto parta = fParticleA[iParticle];
                  //const char *expression = Form("dedx:(sqrt((p_lab*%d)*(p_lab*%d)+%f*%f)-%f)/%d",partz,partz,partm,partm,partm,parta);
                  const char *expression = Form("dedx:(sqrt((p_lab*%d)*(p_lab*%d)+%f*%f)-%f)",partz,partz,partm,partm,partm);
                  project(treeParticle,nameEKECorPart,expression,selection);

                  if (histEKECor==nullptr) {
                    auto nameEKECor = makeName("ekeCor",iAna,iLR,iMult,iSys,iCutTTA);
                    histEKECor = (TH2D *) histEKECorPart -> Clone(nameEKECor);
                    histEKECor -> SetTitle(Form("%s;KE_{Lab} (MeV);dE/dx;",titleEKECorPart));
                  }
                  else
                    histEKECor -> Add(histEKECorPart);

                  histEKECorArray.push_back(histEKECorPart);
                }

                if (drawPtoaR21)
                {
                  auto namePtoaPart = makeName("ptoa",iAna,iLR,iMult,iSys,iCutTTA,iParticle);
                  auto titlePtoaPart = makeTitle("Cor.",iAna,iLR,iMult,iSys,iCutTTA,iParticle);

                  auto histPtoa = new TH1D(namePtoaPart,Form("%s;p_{T}/A (MeV/c);",titlePtoaPart),nbinsPtoa,0,ptoaMax);
                  giSys = iSys;
                  giParticle = iParticle;
                  TCut selection = cut0 * TCut(Form("%s/eff/%d",wProbString,numEventsInAna[iSys])) * cutTTA ;
                  //TCut selection = cut0 * TCut(Form("%s/eff/%d","1",numEventsInAna[iSys])) * cutTTA ;
                  //treeParticle -> Project(namePtoaPart,Form("pt_cm/%d",fParticleA[iParticle]),selection);
                  if (useHandCut) {
                    if (sys==132) {
                      if (iParticle==0) selection = cut0 * TCut("cutg0_132") * TCut(Form("%s/eff/%d","1",numEventsInAna[iSys])) * cutTTA;
                      if (iParticle==1) selection = cut0 * TCut("cutg1_132") * TCut(Form("%s/eff/%d","1",numEventsInAna[iSys])) * cutTTA;
                      if (iParticle==2) selection = cut0 * TCut("cutg2_132") * TCut(Form("%s/eff/%d","1",numEventsInAna[iSys])) * cutTTA;
                      if (iParticle==3) selection = "0";
                      if (iParticle==4) selection = "0";
                    }
                    if (sys==108) {
                      if (iParticle==0) selection = cut0 * TCut("cutg0_108") * TCut(Form("%s/eff/%d","1",numEventsInAna[iSys])) * cutTTA;
                      if (iParticle==1) selection = cut0 * TCut("cutg1_108") * TCut(Form("%s/eff/%d","1",numEventsInAna[iSys])) * cutTTA;
                      if (iParticle==2) selection = cut0 * TCut("cutg2_108") * TCut(Form("%s/eff/%d","1",numEventsInAna[iSys])) * cutTTA;
                      if (iParticle==3) selection = "0";
                      if (iParticle==4) selection = "0";
                    }
                  }
                  project(treeParticle,namePtoaPart,Form("pt_cm/%d",fParticleA[iParticle]),selection);

                  histPtoaArray[iSys][iCutTTA][iParticle] = histPtoa;
                }

                if (drawKeoaR21)
                {
                  auto nameKeoaPart = makeName("ptoa",iAna,iLR,iMult,iSys,iCutTTA,iParticle);
                  auto titleKeoaPart = makeTitle("Cor.",iAna,iLR,iMult,iSys,iCutTTA,iParticle);

                  auto histKeoa = new TH1D(nameKeoaPart,Form("%s;KE_{Lab}/A (MeV);",titleKeoaPart),nbinsKeoaR21,0,keoaMax);
                  giSys = iSys;
                  giParticle = iParticle;
                  TCut selection = cut0 * TCut(Form("%s/eff/%d",wProbString,numEventsInAna[iSys])) * cutTTA ;
                  //treeParticle -> Project(nameKeoaPart,Form("pt_cm/%d",fParticleA[iParticle]),selection);
                  auto partz = fParticleZ[iParticle];
                  auto partm = fParticleMass[iParticle];
                  auto parta = fParticleA[iParticle];
                  const char *expression = Form("(sqrt((p_lab*%d)*(p_lab*%d)+%f*%f)-%f)/%d",partz,partz,partm,partm,partm,parta);
                  project(treeParticle,nameKeoaPart,expression,selection);

                  histKeoaArray[iSys][iCutTTA][iParticle] = histKeoa;
                }

              }

              if (drawCorPID)
              {
                auto namePIDCor = makeName("pidCor",iAna,iLR,iMult,iSys,iCutTTA);
                auto cvsPIDCor = newCvs(namePIDCor,1000,700);
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
                auto nameEKECor = makeName("ekeCor",iAna,iLR,iMult,iSys,iCutTTA);
                auto cvsEKECor = newCvs(nameEKECor,1500,1000);
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
                const char *namePtoa = makeName("ptoa",iAna,iLR,iMult,iSys,iCutTTA);
                auto cvsR21 = newCvs(namePtoa);

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


        if (drawPtoaR21)
        {
          for (auto iCutTTA : cutTTAIdx)
          {
            if (selTTA>=0 && selTTA!=iCutTTA) continue;
            const char *ttaName = fCutTTANames[iCutTTA];
            const char *ttaTitle = fCutTTATitles[iCutTTA];

            for (auto iComb : fSysCombIndx)
            {
              if (selCombR21>=0 && selCombR21!=iComb) continue;

              auto iSys1 = fSysCombIdx[selCombR21][0];
              auto iSys2 = fSysCombIdx[selCombR21][1];
              auto iSys = 100 + iSys1*10 + iSys2;
              
              const char *nameR21 = makeName("r21",iAna,iLR,iMult,iSys,iCutTTA);
              auto titleR21 = makeTitle("Cor.",iAna,iLR,iMult,iSys,iCutTTA);

              auto cvs = newCvs(nameR21,680,550);

              auto hist = new TH2D(nameR21,Form("%s;p_{T}/A (MeV/c);",titleR21),100,0,400,100,0,1.8);
              hist -> Draw();
              for (auto iParticle : fParticleIdx)
              {
                const char *nameParticle = fParticleNames[iParticle];

                auto hist1 = histPtoaArray[iSys1][iCutTTA][iParticle];
                auto hist2 = histPtoaArray[iSys2][iCutTTA][iParticle];

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
              }
            }
          }
        }


        if (drawKeoaR21)
        {
          for (auto iCutTTA : cutTTAIdx)
          {
            if (selTTA>=0 && selTTA!=iCutTTA) continue;
            const char *ttaName = fCutTTANames[iCutTTA];
            const char *ttaTitle = fCutTTATitles[iCutTTA];

            for (auto iComb : fSysCombIndx)
            {
              if (selCombR21>=0 && selCombR21!=iComb) continue;

              auto iSys1 = fSysCombIdx[selCombR21][0];
              auto iSys2 = fSysCombIdx[selCombR21][1];
              auto iSys = iSys1*10 + iSys2;

              
              const char *nameR21 = makeName("r21_keoa",iAna,iLR,iMult,iSys,iCutTTA);
              auto titleR21 = makeTitle("Cor.",iAna,iLR,iMult,iSys,iCutTTA);

              auto cvs = newCvs(nameR21,680,550);
              //cvs -> SetGrid(1,1);

              auto hist = new TH2D(nameR21,Form("%s;KE_{Lab}/A (MeV);",titleR21),100,0,400,100,0,2.0);
              hist -> Draw();
              for (auto iParticle : fParticleIdx)
              {
                const char *nameParticle = fParticleNames[iParticle];

                auto hist1 = histKeoaArray[iSys1][iCutTTA][iParticle];
                auto hist2 = histKeoaArray[iSys2][iCutTTA][iParticle];

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
