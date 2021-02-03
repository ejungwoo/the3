#include "KBGlobal.hh"
#include "init_variables.h"

int cvsXOff = 0; //1300;

void draw_pid()
{
  gStyle -> SetOptStat(0);

  //int selAna = kx6;
  int selAna = kf6;
  //int selAna = kall;
  int selLR = klr;
  int selMult = knAll;
  int selSys = k132;
  int selTTA = kttaAll;

  bool saveFigures = false;

  bool drawRawPID = true;
  int nbinsPID = 500;
  double dedxMax = 2000;
  double pozMax = 2000;

  bool drawRawPIDProjection = false;
  double pidProjRange1 = 1000;
  double pidProjRange2 = 1400;
  bool fitdEdx = true;
  double scaleMaxFitdEdx = .05;

  bool drawCorPID = true;

  TString pathToFigures = "figures/";





  TH1F *histPIDMeta[fNumSystems][fNumParticles] = {0};
  TF1 *f1PIDMean[fNumSystems][fNumParticles] = {0};
  TF1 *f1PIDSigma[fNumSystems][fNumParticles] = {0};
  TCutG *cutgPID[fNumSystems][fNumParticles] = {0};
  TGraph *graphPIDMean[fNumSystems][fNumParticles] = {0};
  TGraph *graphPIDRange[fNumSystems][fNumParticles][2] = {0};

  for (auto iSys : fSystemIdx)
  {
    auto sys = fSystems[iSys];
    TString pidMetaName = Form("data2/Meta_Sn%dKanekoMult50.root",sys);
    TString pidSigmaName = Form("data2/PIDSigma_Sn%dKanekoMult50.root",sys);
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
    const char *anaName = fAnaNames[iAna];
    const char *anaOName = fAnaONames[iAna];
    const char *anaTitle = fAnaTitles[iAna];
    const char *spVersion = fAnaVersion[iAna];

    for (auto iLR : fLRIdx)
    {
      if (selLR>=0 && selLR!=iLR) continue;
      const char *lrName = fLRNames[iLR].Data();
      const char *lrFName = fLRFNames[iLR].Data();
      const char *lrTitle = fLRTitles[iLR].Data();

      for (auto iMult : fMultIdx)
      {
        if (selMult>=0 && selMult!=iMult) continue;
        const char *multName = fMultNames[iMult].Data();
        const char *multFName = fMultFNames[iMult].Data();
        const char *multTitle = fMultTitles[iMult].Data();

        int numEventsInAna[fNumSystems] = {0};

        for (auto iSys : fSystemIdx)
        {
          if (selSys>=0 && selSys!=iSys) continue;
          auto sys = fSystems[iSys];
          const char *sysTitle = fSystemTitles[iSys];

          if (drawRawPID || drawRawPIDProjection)
          {
            auto nameFileAll = Form("data2/%s/sys%d_%s_%s_%s.%s.ana.all.root",anaName,sys,anaName,lrFName,multFName,spVersion);
            cout_info << "File : " << nameFileAll << endl;

            auto treeAll = new TChain("all");
            treeAll -> Add(nameFileAll);

            for (auto iTTA : fCutTTAIdx)
            {
              if (selTTA>=0 && selTTA!=iTTA) continue;
              const char *ttaName = fCutTTANames[iTTA];
              const char *ttaTitle = fCutTTATitles[iTTA];

              TCut selection = fCutTTAValues[iTTA];

              const char *namePIDRaw = Form("pidRaw_%s_%s_%s_%d_%s",anaName,lrName,multName,sys,ttaName);
              TString outPIDRaw = Form("pidRaw_%s_%s_%s_%d_%s",anaOName,lrName,multName,sys,ttaName);
              TString title = Form("[%s]  (Raw)  %s, %s, %s, %s",anaTitle,lrTitle,multTitle,sysTitle,ttaTitle);

              auto histPID = new TH2D(namePIDRaw,title+";p/Z (MeV);dE/dx;",nbinsPID,0,pozMax,nbinsPID,0,dedxMax);
              histPID -> SetMinimum(0.5);
              histPID -> SetMaximum(800);
              treeAll -> Project(namePIDRaw,"dedx:p_lab",selection);

              TCanvas *cvsPIDRaw = nullptr;
              if (drawRawPID) {
                cvsPIDRaw = new TCanvas(namePIDRaw,namePIDRaw,cvsXOff,100,1000,700);
                cvsPIDRaw -> SetLogz();
                histPID -> Draw("colz");
                for (auto iParticle : fParticleIdx)
                  graphPIDMean[iSys][iParticle] -> Draw("samel");

                if (saveFigures)
                  cvsPIDRaw -> SaveAs(pathToFigures+outPIDRaw+".png"); 
              }

              if (drawRawPIDProjection)
              {
                auto binn = ebinning(histPID);
                //int numProjections = 200;
                int numProjections = nbinsPID;
                int dbin = nbinsPID/numProjections;
                int countProj = 0;

                TGraph *graphFitPIDMean[5] = {0};
                for (auto iParticle : fParticleIdx)
                  graphFitPIDMean[iParticle] = new TGraph();

                for (auto idxProjection=0; idxProjection<numProjections; ++idxProjection)
                {
                  auto bin1 = (idxProjection)*dbin+1;
                  auto bin2 = (idxProjection+1)*dbin;
                  auto poz1 = binn.lowEdge(bin1);
                  auto poz2 = binn.highEdge(bin2);
                  auto pozC = (poz1 + poz2)/2.;

                  const char *nameProj = Form("%s_dedx_%d",namePIDRaw,idxProjection);
                  auto histProj = (TH1D *) histPID -> ProjectionY(nameProj,bin1,bin2);

                  if (pozC>pidProjRange1 && pozC<pidProjRange2)
                  {
                    cout_info << poz1 << " " << poz2 << endl;
                    TString titleProj = TString("p/Z=")+poz1+"~"+poz2+";dE/dx;";
                    histProj -> SetTitle(titleProj);

                    auto cvsProj = new TCanvas(nameProj,nameProj,cvsXOff+20*(countProj+1),20*(countProj+1),1000,550);
                    cvsProj -> SetLogy();
                    histProj -> Draw();

                    double scaleAmp = 1.;
                    const char *nameProjParticle = Form("%s_dedx_%d_dist",namePIDRaw,idxProjection);
                    auto f1dEdxTotal = new TF1(nameProjParticle,"gaus(0)+gaus(3)+gaus(6)+gaus(9)+gaus(12)",0,dedxMax);
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

                      nameProjParticle = Form("%s_dedx_%d_dist_ori_%s",namePIDRaw,idxProjection,fParticleNames[iParticle].Data());
                      auto f1dEdxParticle = new TF1(nameProjParticle,"gaus(0)",0,dedxMax);
                      f1dEdxParticle -> SetLineColor(kGray);
                      f1dEdxParticle -> SetParameters(dedxAmp,dedxMean,dedxSigma);
                      f1dEdxParticle -> SetNpx(1000);
                      f1dEdxParticle -> Draw("samel");


                      f1dEdxTotal -> SetParameter(0+3*iParticle,dedxAmp);
                      f1dEdxTotal -> SetParameter(1+3*iParticle,dedxMean);
                      f1dEdxTotal -> SetParameter(2+3*iParticle,dedxSigma);

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
                      nameProjParticle = Form("%s_dedx_%d_dist_fit_%s",namePIDRaw,idxProjection,fParticleNames[iParticle].Data());
                      auto f1dEdxParticle = new TF1(nameProjParticle,"gaus(0)",0,dedxMax);
                      f1dEdxParticle -> SetLineColor(kGray+1);
                      auto dedxAmp = f1dEdxTotal -> GetParameter(0+3*iParticle);
                      auto dedxMean = f1dEdxTotal -> GetParameter(1+3*iParticle);
                      auto dedxSigma = f1dEdxTotal -> GetParameter(2+3*iParticle);
                      f1dEdxParticle -> SetParameters(dedxAmp,dedxMean,dedxSigma);
                      f1dEdxParticle -> SetNpx(1000);
                      f1dEdxParticle -> Draw("samel");

                      graphFitPIDMean[iParticle] -> SetPoint(graphFitPIDMean[iParticle]->GetN(),pozC,dedxMean);
                    }

                    f1dEdxTotal -> Draw("samel");

                    if (drawRawPID)
                    {
                      cvsPIDRaw -> cd();
                      auto line1 = new TLine(poz1,0,poz1,20);
                      line1 -> SetLineColor(kRed);
                      line1 -> Draw("samel");
                      auto line2 = new TLine(poz2,0,poz2,20);
                      line2 -> SetLineColor(kRed);
                      line2 -> Draw("samel");

                      if (fitdEdx)
                        for (auto iParticle : fParticleIdx) {
                          //graphFitPIDMean[iParticle] -> SetMarkerStyle(25);
                          //graphFitPIDMean[iParticle] -> Draw("samep");
                          graphFitPIDMean[iParticle] -> SetLineColor(kGray+1);
                          graphFitPIDMean[iParticle] -> Draw("samel");
                        }

                    }
                    countProj++;
                  }

                }
              }

            }
          }

          auto nameFileMult = Form("data2/%s/sys%d_%s_%s_%s.%s.ana.mult.root",anaName,sys,anaName,lrFName,multFName,spVersion);
          cout_info << "File : " << nameFileMult << endl;

          auto treeMult = new TChain("mult");
          treeMult -> Add(nameFileMult);

          numEventsInAna[iSys] = treeMult -> GetEntries("1");
          cout_info << "Multiplicity in " << sys << " : " << numEventsInAna[iSys] << endl;

          if (drawCorPID)
          {
            for (auto iTTA : fCutTTAIdx)
            {
              if (selTTA>=0 && selTTA!=iTTA) continue;
              const char *ttaName = fCutTTANames[iTTA];
              const char *ttaTitle = fCutTTATitles[iTTA];

              TCut selection = fCutTTAValues[iTTA];

              TH2D *histPIDCor = nullptr;
              const char *namePIDCor = Form("pidCor_%s_%s_%s_%d_%s",anaName,lrName,multName,sys,ttaName);
              TString outPIDCor = Form("pidCor_%s_%s_%s_%d_%s",anaOName,lrName,multName,sys,ttaName);
              TString titlePIDCor = Form("[%s]  (#timesP/E/N)  %s, %s, %s, %s",anaTitle,lrTitle,multTitle,sysTitle,ttaTitle);
              titlePIDCor += ";p/Z (MeV);dE/dx;";

              for (auto iParticle : fParticleIdx)
              {
                const char *nameParticle = fParticleNames[iParticle];
                const char *titleParticle = fParticleTitles[iParticle];
                auto nameFileParticle = Form("data2/%s/sys%d_%s_%s_%s.%s.ana.%s.root",anaName,sys,anaName,lrFName,multFName,spVersion,nameParticle);
                cout_info << "File : " << nameFileParticle << endl;

                auto treeParticle = new TChain(nameParticle);
                treeParticle -> Add(nameFileParticle);

                const char *namePIDCorPart = Form("pidCor_%s_%s_%s_%d_%s_%s",anaName,lrName,multName,sys,ttaName,nameParticle);
                TString title = Form("[%s]  %s, %s, %s, %s, %s",anaTitle,lrTitle,multTitle,sysTitle,ttaTitle,titleParticle);

                auto histPIDCorPart = new TH2D(namePIDCorPart,title+";p/Z (MeV);dE/dx;",nbinsPID,0,pozMax,nbinsPID,0,dedxMax);
                histPIDCorPart -> SetMinimum(0.5);
                histPIDCorPart -> SetMaximum(800);
                TCut selection = Form("prob/eff/%d",numEventsInAna[iSys]);
                treeParticle -> Project(namePIDCorPart,"dedx:p_lab",selection);

                if (iParticle==0) {
                  histPIDCor = (TH2D *) histPIDCorPart -> Clone(namePIDCor);
                  histPIDCor -> SetTitle(titlePIDCor);
                }
                else
                  histPIDCor -> Add(histPIDCorPart);
              }

              TCanvas *cvsPIDCor = nullptr;
              if (drawRawPID) {
                cvsPIDCor = new TCanvas(namePIDCor,namePIDCor,cvsXOff,100,1000,700);
                cvsPIDCor -> SetLogz();
                histPIDCor -> Draw("colz");
                for (auto iParticle : fParticleIdx)
                  graphPIDMean[iSys][iParticle] -> Draw("samel");

                if (saveFigures)
                  cvsPIDCor -> SaveAs(pathToFigures+outPIDCor+".png"); 
              }
            }
          }


        }
      }
    }
  }
}
