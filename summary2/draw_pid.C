#include "init_variables.h"

void draw_pid()
{
  gStyle -> SetOptStat(0);

  int selAna = kall;
  int selLR = klr;
  int selMult = knAll;
  int selSys = k132;
  int selTTA = kttaAll;

  bool saveFigures = false;

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
    auto filePIDMeta = new TFile(pidSigmaName,"read");
    auto filePIDSigma = new TFile(pidSigmaName,"read");

    TString name = Form("system%d",iSys);

    for (auto iParticle : fParticleIdx)
    {
      auto pdg = fParticlePDGs[iParticle];
      auto histMeta = = (TH1F*) filePIDMeta -> Get(TString::Format("Distribution%d", pdg));
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

        for (auto iSys : fSystemIdx)
        {
          if (selSys>=0 && selSys!=iSys) continue;
          auto sys = fSystems[iSys];
          const char *sysTitle = fSystemTitles[iSys];

          auto nameFileAll = Form("data2/%s/sys%d_%s_%s_%s.%s.ana.all.root",anaName,sys,anaName,lrFName,multFName,spVersion);
          cout << "file: " << nameFileAll << endl;

          auto treeAll = new TChain("all");
          treeAll -> Add(nameFileAll);

          for (auto iTTA : fCutTTAIdx)
          {
            if (selTTA>=0 && selTTA!=iTTA) continue;
            const char *ttaName = fCutTTANames[iTTA];
            const char *ttaTitle = fCutTTATitles[iTTA];

            TCut cut = "(nc+nl>15)&&(dpoca)<15";
            cut = cut + fCutTTAValues[iTTA];

            const char *name = Form("%s_%s_%s_%d_%s",anaName,lrName,multName,sys,ttaName);
            TString outName = Form("pid_%s_%s_%s_%d_%s",anaOName,lrName,multName,sys,ttaName);
            TString title = Form("[%s]  %s, %s, %s, %s",anaTitle,lrTitle,multTitle,sysTitle,ttaTitle);

            cout << "+drawing " << name << endl;

            auto cvs = new TCanvas(name,name,1300,100,1000,700);
            cvs -> SetLogz();

            //auto hist = new TH2D(name,title+";p/Z (MeV);dE/dx;",1000,0,3000,1000,0,1000);
            auto hist = new TH2D(name,title+";p/Z (MeV);dE/dx;",500,0,3000,500,0,1000);
            treeAll -> Draw(Form("dedx:p_lab>>%s",name),cut,"colz");

            for (auto iParticle : fParticleIdx) {
              graphPIDMean[iSys][iParticle] -> Draw("samel");
              //graphPIDRange[iSys][iParticle][0] -> SetLineStyle(2);
              //graphPIDRange[iSys][iParticle][1] -> SetLineStyle(2);
              //graphPIDRange[iSys][iParticle][0] -> Draw("samel");
              //graphPIDRange[iSys][iParticle][1] -> Draw("samel");
            }

            if (saveFigures)
              cvs -> SaveAs(pathToFigures+outName+".png"); 
          }
        }
      }
    }
  }
}
