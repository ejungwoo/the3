#include "/Users/ejungwoo/config/ejungwoo.h"
using ejungwoo::variable;
using ejungwoo::binning;
using ejungwoo::titles;
using ejungwoo::setpar;
using ejungwoo::drawall;
using ejungwoo::drawsaveall;

TObjArray graph_is(int idx_particle, TH1 *hist132, TH1 *hist108, double num_tracks_cut);
TH1D *hist_is(int idx_particle, TH1 *hist132, TH1 *hist108, double num_tracks_cut);

void draw_all(
  TString fRunOption = "full"
  //TString fRunOption = "midy"
  //TString fRunOption = "y1"
  //TString fRunOption = "y2"
  //TString fRunOption = "box2"
)
{
  ejungwoo::gcvspos(1000);
  ejungwoo::gstat(0);
  ejungwoo::gsave(0);
  //ejungwoo::gautosavetp(1);
  ejungwoo::glegendleft();

  //double drawFitSaveThetaCM = 30;
  //double drawFitSaveThetaCM = 60;
  //double drawFitSaveThetaCM = 120;
  //double drawFitSaveThetaCM = 150;
  double drawFitSaveThetaCM = 0;

  bool ana_p = 1;
  int printTree = 0;
  int drawPID = 0;
  int drawPYDist0 = 0;
  int drawPYDistPBoundary = 1;
  int drawPYDistRangeBox = 1;
  int drawPYTTBoundary = 1;
  int drawPYDist1 = 0;
  int drawPtoA = 0;
  int drawAB = 0;
  int drawR21NZ = 2;
  int drawTemperature = 1;

  double probLL = 0.2;
  double effLL = 0.05;
  double sdHL = 2.0;
  double ptoaRange[] = {0,400};
  double y0Range[] = {0,0.4};
  int nBinsY0 = 100;
  //int rebinY0 = 100;
  int rebinY0 = 10;

  double thetaCMRange[] = {0,180};
  double thetaLabRange[] = {30,60};



  bool ana_e = 0;
  int drawETheta = 1;
  int drawEDistTheta = 1;
  int drawER21 = 1;

  int numTheta = 4;
  double thetaFull = 90;
  double thetaMax = 80;
  double dTheta = thetaMax / numTheta;



  auto setParameters = [&ptoaRange, &y0Range](double p0, double p1, double y0, double y1, int nx=100) {
    ptoaRange[0] = p0;
    ptoaRange[1] = p1;
    y0Range[0] = y0;
    y0Range[1] = y1;
  };

  if (fRunOption=="full") setParameters( 0, 400,  -1, 2.0);
  if (fRunOption=="midy") setParameters( 0, 400, -.2,  .2);
  if (fRunOption=="y1")   setParameters( 0, 400,   0,  .4);
  if (fRunOption=="y2")   setParameters( 0, 400,  .4,  .8);
  if (fRunOption=="box1") setParameters( 0, 400,  .4,  .8);
  if (fRunOption=="box2") setParameters( 100, 200,  0,  .4, 1);

  if (drawFitSaveThetaCM>0) {
    thetaCMRange[0] = drawFitSaveThetaCM-1;
    thetaCMRange[1] = drawFitSaveThetaCM+1;
    printTree = 0;
    drawPID = 0;
    drawPYDist0 = 0;
    drawPYDist1 = 2;
  }

  TString fileTags[] = {"all_45_49", "all_50_54"};
  int mult0 = 45;
  int mult1 = 54;
  double dPhi = 120;

  //TString fileTags[] = {"all_55_100"};
  //int mult0 = 55;
  //int mult1 = 99;
  //double dPhi = 120;


  auto funcTemperature = new TF1("temperature_r_","14.3 / (TMath::Log(1.59*x[0]))",0,10000);
  ejungwoo::gversion(Form("e%d",mult0));

  const int nSyss = 4;
  const int iSyss[] = {0,1,2,3};
  const int sysBeam[] = {132,108,124,112};
  const int sysTarget[] = {124,112,112,124};
  const TString sysNames0[] = {"sys132","sys108","sys124","sys112"};
  const TString sysNames1[] = {"132+124","108+112","124+112","112+124"};
  const TString sysEnergyLossFile[] = {"data/Sn132Sn124.txt","data/Sn108Sn112.txt","data/Sn112Sn124.txt","data/Sn112Sn124.txt"};
  int numEvents[] = {0,0,0,0};

  const int nPIDs = 5;
  const int iPIDs[] = {0,1,2,3,4};
  const int particleZ[] = {1,1,1,2,2};
  const int particleN[] = {0,1,2,1,2};
  const int particleA[] = {1,2,3,3,4};
  const double particleMass[] = {938.272, 1875.61, 2808.92, 2808.39, 3727.38};
  const TString particleNames[] = {"p","d","t","he3","he4"};

  const int nTests = 6;
  const int iTests[] = {0,1,2,3,4,5};
  const double testPValues[] = {500,1000,1500,2000,2500,3000};
  const int particleiTestCut[] = {1, 2, 4, 4, 5};
  const double testThetas[] = {0.1,0.2,0.5,1,3,5,10,15,20,25,30,35,40,45,50,55,60,65,70,75,80,85,87,89,89.9};

  TGraph *graphELossInTarget[4];
  for (auto iSys : iSyss) {
    graphELossInTarget[iSys] = new TGraph(sysEnergyLossFile[iSys], "%lg %*lg %lg");
  }

  auto setSystem = [&numEvents](int iSys) {
    setpar("nevents",numEvents[iSys]);
  }; setSystem(0);

  auto setParticle = [&particleZ, &particleA, &particleMass](int iPID) {
    setpar("a",particleA[iPID]);
    setpar("z",particleZ[iPID]);
    setpar("m",particleMass[iPID]);
  }; setParticle(0);

  auto setPIDAtt = [](TH1 *hist, int iPID, int rebin=0, double xMax=0, bool makedndx=0, TString title="") {
    if (rebin>1) hist -> Rebin(rebin);
    if (xMax > 0) hist -> GetXaxis() -> SetRangeUser(0,xMax);
    if (makedndx) hist = ejungwoo::dndx(hist);
    if (!title.IsNull()) hist -> SetTitle(title);
    double msize; int mstyle, mcolor, lcolor;
    if      (iPID==0) { mstyle=20, msize=1, mcolor=kBlack, lcolor=mcolor; }
    else if (iPID==1) { mstyle=21, msize=1, mcolor=kRed+1, lcolor=mcolor; }
    else if (iPID==2) { mstyle=22, msize=1.15, mcolor=kBlue+1, lcolor=mcolor; }
    else if (iPID==3) { mstyle=23, msize=1.15, mcolor=kGreen+3, lcolor=mcolor; }
    //else if (iPID==4) { mstyle=30, msize=1.2, mcolor=kViolet-5, lcolor=mcolor; }
    else if (iPID==4) { mstyle=30, msize=1.2, mcolor=kBlack, lcolor=mcolor; }
    hist -> SetMarkerStyle(mstyle);
    hist -> SetMarkerColor(mcolor);
    hist -> SetMarkerSize(msize);
    hist -> SetLineColor(lcolor);
    return hist;
  };

  auto setCombAtt = [](TH1 *hist, int iComb, int rebin=0, double xMax=0, bool makedndx=0, TString title="") {
    if (rebin>1) hist -> Rebin(rebin);
    if (xMax > 0) hist -> GetXaxis() -> SetRangeUser(0,xMax);
    if (makedndx) hist = ejungwoo::dndx(hist);
    if (!title.IsNull()) hist -> SetTitle(title);
    double msize; int mstyle, mcolor, lcolor;
    if      (iComb==0) { mstyle=20, msize=1, mcolor=kBlack, lcolor=mcolor; }
    else if (iComb==1) { mstyle=21, msize=1, mcolor=kRed+1, lcolor=mcolor; }
    else if (iComb==2) { mstyle=22, msize=1.15, mcolor=kBlue+1, lcolor=mcolor; }
    else if (iComb==3) { mstyle=23, msize=1.15, mcolor=kGreen+3, lcolor=mcolor; }
    hist -> SetMarkerStyle(mstyle);
    hist -> SetMarkerColor(mcolor);
    hist -> SetMarkerSize(msize);
    hist -> SetLineColor(lcolor);
    return hist;
  };

  auto setSysAtt = [](TH1 *hist, int iSys, int rebin=0, double xMax=0, bool makedndx=0, TString title="") {
    if (rebin>1) hist -> Rebin(rebin);
    if (xMax > 0) hist -> GetXaxis() -> SetRangeUser(0,xMax);
    if (makedndx) hist = ejungwoo::dndx(hist);
    if (!title.IsNull()) hist -> SetTitle(title);
    double msize; int mstyle, mcolor, lcolor;
    if      (iSys==0) { mstyle=24, msize=1, mcolor=kRed, lcolor=mcolor; }
    else if (iSys==1) { mstyle=25, msize=1, mcolor=kAzure+8, lcolor=mcolor; }
    else if (iSys==2) { mstyle=26, msize=1, mcolor=kSpring-6, lcolor=mcolor; }
    else if (iSys==3) { mstyle=30, msize=1.2, mcolor=kViolet-5, lcolor=mcolor; }
    hist -> SetMarkerStyle(mstyle);
    hist -> SetMarkerColor(mcolor);
    hist -> SetMarkerSize(msize);
    hist -> SetLineColor(lcolor);
    return hist;
  };

  setpar("cut_mult","1");
  TChain *chain[nSyss][nPIDs];
  TChain *mult[nSyss];
  for (auto iSys : iSyss) {
    for (auto iPID : iPIDs) {
      chain[iSys][iPID] = new TChain(particleNames[iPID]);
      for (auto tag : fileTags) {
        auto fileName = Form("data/%s.%s.NewAna.2100.de29163.ana.particle.root",sysNames0[iSys].Data(),tag.Data());
        chain[iSys][iPID] -> Add(fileName);
        if (printTree) {
          chain[iSys][iPID] -> Print("toponly");
          if (printTree==2) return;
          printTree = 0;
        }
      }
    }

    mult[iSys] = new TChain("mult");
    for (auto tag : fileTags) {
      auto fileName = Form("data/%s.%s.NewAna.2100.de29163.ana.particle.root",sysNames0[iSys].Data(),tag.Data());
      mult[iSys] -> Add(fileName);
    }
    numEvents[iSys] = mult[iSys] -> GetEntries(ejungwoo::getpar("cut_mult"));
  }

  variable vptoa("ptoa","pt_cm/$$(a)");
  variable vy0("y0","fy_cm/(by_cm/2)");

  setpar("cut_PENS",Form("prob/eff/$$(nevents)*(prob>%.2f&&eff>%.2f&&abs(sd)<%.2f)",probLL,effLL,sdHL));
  setpar("cut_ptoa",Form("($$(ptoa)>%.2f&&$$(ptoa)<%.2f)",ptoaRange[0],ptoaRange[1]));
  setpar("cut_y0",Form("($$(y0)>%.2f&&$$(y0)<%.2f)",y0Range[0],y0Range[1]));
  setpar("cut_tl",Form("(theta_lab>%.4f&&theta_lab<%.4f)",TMath::DegToRad()*thetaLabRange[0],TMath::DegToRad()*thetaLabRange[1]));
  setpar("cut_tc",Form("(theta_cm>%.4f&&theta_cm<%.4f)",TMath::DegToRad()*thetaCMRange[0],TMath::DegToRad()*thetaCMRange[1]));

  setpar("cut_allp","$$(cut_PENS)*$$(cut_ptoa)*$$(cut_y0)");

  double x2 = 400;

  //auto vpydist0 = new variable("pydist0","pt_cm/$$(a):fy_cm/(by_cm/2)","$$(cut_PENS)",titles("","y0","pt/a"),binning(100,-1,2),binning(nBinsY0,0,x2));
  auto vpydist0 = new variable("pydist0","pt_cm/$$(a):fy_cm/(by_cm/2)","$$(cut_allp)",titles("","y0","pt/a"),binning(100,-1,2),binning(nBinsY0,0,x2));
  auto vpydist1 = new variable("pydist1","pt_cm/$$(a):fy_cm/(by_cm/2)","$$(cut_allp)" ,titles("","y0","pt/a"),binning(100,-1,2),binning(nBinsY0,0,x2));
  //auto vpydist1 = new variable("pydist1","pt_cm/$$(a):fy_cm/(by_cm/2)","$$(cut_PENS)" ,titles("","y0","pt/a"),binning(100,-1,2),binning(nBinsY0,0,x2));
  auto vpid = new variable("pid", "dedx:p_lab", "$$(cut_PENS)", titles(" ","p_{Lab}/Z (MeV/c^{2})","dE/dx"), binning(400,0,2500),binning(400,0,600)); 

  double nucleonMass = 931.5;
  double targetThickness = 0.8;
  double tempBeamEnergyTargetPlane = 280;

  while (ana_p)
  {
    TGraph *graphPYDistAtP[nSyss][nPIDs][nTests];
    for (auto iSys : iSyss)
    {
      double energyLossInTarget = graphELossInTarget[iSys] -> Eval(tempBeamEnergyTargetPlane)*targetThickness*1000/2.;
      double energyPerN = (tempBeamEnergyTargetPlane*sysBeam[iSys] - energyLossInTarget)/sysBeam[iSys];
      double EBeam = energyPerN*sysBeam[iSys] + sysBeam[iSys]*nucleonMass;
      double PBeam = sqrt(EBeam*EBeam - sysBeam[iSys]*sysBeam[iSys]*nucleonMass*nucleonMass);
      TLorentzVector LV(0,0,PBeam,EBeam);
      double beta = PBeam/(LV.Gamma()*sysBeam[iSys]*nucleonMass + sysTarget[iSys]*nucleonMass);
      double by = LV.Rapidity();
      TVector3 vBeam(0,0,-beta);

      for (auto iPID : iPIDs) {
        setParticle(iPID);
        for (auto iTest : iTests) {
          auto pValue = testPValues[iTest];
          graphPYDistAtP[iSys][iPID][iTest] = new TGraph();
          auto graph = graphPYDistAtP[iSys][iPID][iTest];
          for (auto theta : testThetas) {
            TLorentzVector pCM(pValue*TMath::Sin(TMath::DegToRad()*theta),0,pValue*TMath::Cos(TMath::DegToRad()*theta),particleMass[iPID]);
            pCM.Boost(vBeam);
            auto fy = pCM.Rapidity();
            if (fy!=fy) continue;
            graph -> SetPoint(graph->GetN(), fy/(by/2), pCM.Pt()/particleA[iPID]);
          }
        }
      }
    }

    TH2D *histPY[nSyss][nPIDs] = {0};
    TH2D *histAB[nSyss][nPIDs] = {0};
    TH1D *histMP[nSyss][nPIDs] = {0};

    for (auto iSys : iSyss)
    {
      setSystem(iSys);

      //================================================================================================================================================
      if (drawPID)
        for (auto var : {vpid}) {
          for (auto iPID : iPIDs) {
            setParticle(iPID);
            var -> drawaddhist(sysNames0[iSys]+"_"+var->getName(), chain[iSys][iPID],"logz") -> expandRange();
          }
        } if (drawPID==2) continue;

      //================================================================================================================================================
      if (drawPYDist0)
        for (auto var : {vpydist0}) {
          TString groupName = sysNames0[iSys] + "_" + var -> getName();
          for (auto iPID : iPIDs) {
            setParticle(iPID);
            auto hist = (TH2D *) var -> draw(chain[iSys][iPID]);
            ejungwoo::addnext(groupName, hist);

            if (drawPYDistPBoundary)
              for (auto iTest : iTests) {
                if (iTest > particleiTestCut[iPID]) break;
                auto graph = graphPYDistAtP[iSys][iPID][iTest];
                graph -> SetLineColor(kRed-7);
                graph -> SetLineStyle(2);
                if (iTest%2==1) {
                  graph -> SetLineStyle(1);
                }
                ejungwoo::addsame(groupName, graph, "lcsame colorx addx");
              }

            if (drawPYDistRangeBox) {
              auto graphRangeBox = new TGraph();
              graphRangeBox -> SetLineColor(kGray+2);
              graphRangeBox -> SetPoint(graphRangeBox->GetN(), y0Range[0], ptoaRange[0]);
              graphRangeBox -> SetPoint(graphRangeBox->GetN(), y0Range[0], ptoaRange[1]);
              graphRangeBox -> SetPoint(graphRangeBox->GetN(), y0Range[1], ptoaRange[1]);
              graphRangeBox -> SetPoint(graphRangeBox->GetN(), y0Range[1], ptoaRange[0]);
              graphRangeBox -> SetPoint(graphRangeBox->GetN(), y0Range[0], ptoaRange[0]);
              ejungwoo::addsame(groupName,graphRangeBox,"colorx samel addx");
            }

            if (drawPYTTBoundary) {
              for (auto tt : {30,60,90,120,150}){
                if (tt==90) {
                  auto pyboundary = new TGraph();
                  pyboundary -> SetPoint(0,0,0);
                  pyboundary -> SetPoint(1,0,1000);
                  pyboundary -> SetLineColor(kBlue-3);
                  pyboundary -> SetLineStyle(2);
                  ejungwoo::addsame(groupName,pyboundary,"colorx samel addx");
                }
                else {
                  auto file = new TFile(Form("data/pyBoundary.ttcm%d.sys%d.%s.root",int(tt),sysBeam[iSys],particleNames[iPID].Data()),"read");
                  auto pyboundary = (TF1 *) file -> Get("pyfit");
                  pyboundary -> SetRange(-1,2);
                  pyboundary -> SetLineColor(kBlue-3);
                  pyboundary -> SetLineStyle(2);
                  pyboundary -> SetLineWidth(1);
                  ejungwoo::addsame(groupName,pyboundary,"colorx samel addx");
                }
              }
            }
          }
        } if (drawPYDist0==2) continue;

      //================================================================================================================================================
      for (auto var : {vpydist1}) {
        TString groupName = sysNames0[iSys] + "_" + var -> getName();
        for (auto iPID : iPIDs) {
          setParticle(iPID);
          histPY[iSys][iPID] = (TH2D *) var -> draw(chain[iSys][iPID]);
          if (drawPYDist1)
            ejungwoo::addnext(groupName, histPY[iSys][iPID]);
          if (drawFitSaveThetaCM) {
            TGraph *graph;
            if (abs(drawFitSaveThetaCM-90)<45)
              graph = (TGraph *) ejungwoo::histana2(histPY[iSys][iPID], 2, 20, "m").graphana();
            else
              graph = (TGraph *) ejungwoo::histana2(histPY[iSys][iPID], 1, 20, "m").graphana();
            graph -> SetMarkerColor(kBlack);
            graph -> SetMarkerStyle(24);
            graph -> SetMarkerSize(2);
            ejungwoo::addsame(groupName, graph, "samep colorx");

            auto fit = ejungwoo::fitgraph(graph, "[0]*x*x+[1]*x", "RQ0");
            fit -> SetName("pyfit");
            ejungwoo::addsame(groupName, fit, "samel colorx");

            auto file = new TFile(Form("data/pyBoundary.ttcm%d.sys%d.%s.root",int(drawFitSaveThetaCM),sysBeam[iSys],particleNames[iPID].Data()),"recreate");
            fit -> Write();
          }
        }
      } if (drawPYDist1==2) continue;

      TString groupName = sysNames0[iSys] + "_ab";
      for (auto iPID : iPIDs) {
        histAB[iSys][iPID] = (TH2D *) histPY[iSys][iPID] -> Clone();
        histAB[iSys][iPID] -> Rebin2D(10,10);
        if (drawAB)
          ejungwoo::addnext(groupName, histAB[iSys][iPID]);
      }

      //================================================================================================================================================
      for (auto iPID : iPIDs) {
        TString groupName = sysNames0[iSys] + "_ptoa";
        setParticle(iPID);
        histMP[iSys][iPID] = (TH1D *) histPY[iSys][iPID] -> ProjectionY();
        auto hist = setPIDAtt(histMP[iSys][iPID],iPID,rebinY0,ptoaRange[1],true,groupName+"_"+particleNames[iPID]+";pt/a;dN/d(pt/a)");
        if (drawPtoA)
          ejungwoo::addsame(groupName, hist, "histpl colorx", particleNames[iPID]) -> expandRange();
      } if (drawPtoA==2) continue;

    } // iSys

    if (drawR21NZ) {
      for (auto icomb : {0}) {
        int iSys1 = 0;
        int iSys2 = 1;

        TString groupName = TString("combiZ") + icomb + "_R21";
        //TString groupNameZ = TString("combiZ") + icomb + "_R21";
        //TString groupNameN = TString("combiN") + icomb + "_R21";
        TString groupName1 = TString("combiN") + icomb + "_R21";

        //const int countGroup = 2;
        //const int indexNext = 1;
        const int countGroup = 1;
        const int indexNext = 0;
        int offsetZ = -3;
        int indexGroup = 0;

        auto vR21N = new variable("histN","","",";Z;R21",100,-1,3,100,0,2);
        auto vR21Z = new variable("histZ","","",";N;R21",100,-1,3,100,0,2);
        auto vR21A = new variable("histA","","",";N;R21",100,-4,3,100,0,2);

        //for (auto binx : {4,5,6}) { for (auto biny : {2,3,4}) {
        for (auto binx : {4,5}) { for (auto biny : {2,3}) {
        //for (auto binx : {5}) { for (auto biny : {3}) {
            auto graphN0 = ejungwoo::new_g(); graphN0 -> SetMarkerStyle(ejungwoo::markeri(0)); //graphN0 -> SetMarkerSize(1.5);
            auto graphN1 = ejungwoo::new_g(); graphN1 -> SetMarkerStyle(ejungwoo::markeri(1)); //graphN1 -> SetMarkerSize(1.5);
            auto graphN2 = ejungwoo::new_g(); graphN2 -> SetMarkerStyle(ejungwoo::markeri(2)); //graphN2 -> SetMarkerSize(1.5);
            auto graphZ1 = ejungwoo::new_g(); graphZ1 -> SetMarkerStyle(ejungwoo::markeri(3)); //graphZ1 -> SetMarkerSize(1.5);
            auto graphZ2 = ejungwoo::new_g(); graphZ2 -> SetMarkerStyle(ejungwoo::markeri(4)); //graphZ2 -> SetMarkerSize(1.5);

            auto histR21N = vR21N -> draw(); histR21N -> SetTitle(Form("%d,%d;Z;R21",binx,biny));
            auto histR21Z = vR21Z -> draw(); histR21Z -> SetTitle(Form("%d,%d;N;R21",binx,biny));
            auto histR21A = vR21A -> draw(); histR21A -> SetTitle(Form("%d,%d;N;R21",binx,biny));

            ejungwoo::add(groupName1, countGroup*indexGroup+indexNext, histR21A, "pl colorx") -> legendlt();
            //ejungwoo::add(groupName1, countGroup*indexGroup+indexNext, histR21N, "pl colorx") -> legendlt();
            //ejungwoo::add(groupName1, countGroup*indexGroup+0, histR21Z, "pl colorx") -> legendlt();

            for (auto iPID : iPIDs) {
              auto val2 = histAB[iSys1][iPID] -> GetBinContent(binx,biny);
              auto val1 = histAB[iSys2][iPID] -> GetBinContent(binx,biny);
              double r21 = val2 / val1;

              if (particleN[iPID]==0) graphN0 -> SetPoint(graphN0 -> GetN(), particleZ[iPID], r21);
              if (particleN[iPID]==1) graphN1 -> SetPoint(graphN1 -> GetN(), particleZ[iPID], r21);
              if (particleN[iPID]==2) graphN2 -> SetPoint(graphN2 -> GetN(), particleZ[iPID], r21);
              if (particleZ[iPID]==1) graphZ1 -> SetPoint(graphZ1 -> GetN(), particleN[iPID]+offsetZ , r21);
              if (particleZ[iPID]==2) graphZ2 -> SetPoint(graphZ2 -> GetN(), particleN[iPID]+offsetZ , r21);

              double xN = particleN[iPID]; double yN = r21;
              double xZ = particleZ[iPID]; double yZ = r21;
              if (iPID==0) { yN = yN + 0.2; xZ = xZ - 0.4; xN = xN + offsetZ; }
              if (iPID==1) { yN = yN + 0.2; xZ = xZ - 0.4; xN = xN + offsetZ; }
              if (iPID==2) { yN = yN + 0.2; xZ = xZ - 0.4; xN = xN + offsetZ; }
              if (iPID==3) { yN = yN - 0.2; xZ = xZ + 0.4; xN = xN + offsetZ; }
              if (iPID==4) { yN = yN - 0.2; xZ = xZ + 0.4; xN = xN + offsetZ; }

              auto ttN = new TText(xN, yN, particleNames[iPID]); ttN -> SetTextAlign(22); ttN -> SetTextColor(kBlue);
              auto ttZ = new TText(xZ, yZ, particleNames[iPID]); ttZ -> SetTextAlign(22); ttZ -> SetTextColor(kBlue);

              ejungwoo::add(groupName1, countGroup*indexGroup+0, ttN, "addx same");
              ejungwoo::add(groupName1, countGroup*indexGroup+indexNext, ttZ, "addx same");
            }
            auto fitN1 = ejungwoo::fitgraph(graphN1,"pol1");
            auto beta = fitN1 -> GetParameter(1);
            auto fitZ1 = ejungwoo::fitgraph(graphZ1,"pol1");
            fitZ1 -> SetLineColor(kGreen+3);
            auto alpha = fitZ1 -> GetParameter(1);

            ejungwoo::add(groupName1, countGroup*indexGroup+indexNext, fitN1, "pl colorx", Form("#beta=%.2f",beta));
            ejungwoo::add(groupName1, countGroup*indexGroup+indexNext, graphN0, "pl colorx", "N=0");
            ejungwoo::add(groupName1, countGroup*indexGroup+indexNext, graphN1, "pl colorx", "N=1");
            ejungwoo::add(groupName1, countGroup*indexGroup+indexNext, graphN2, "pl colorx", "N=2");

            ejungwoo::add(groupName1, countGroup*indexGroup+0, fitZ1, "pl colorx", Form("#alpha=%.2f",alpha));
            ejungwoo::add(groupName1, countGroup*indexGroup+0, graphZ1, "pl colorx", "Z=1");
            ejungwoo::add(groupName1, countGroup*indexGroup+0, graphZ2, "pl colorx", "Z=2");

            indexGroup++;
          }
        }
      }
    } if (drawR21NZ) break;

    //================================================================================================================================================
    if (drawTemperature && histMP[0][0]!=nullptr) {
      TString groupNameR = "rheh";
      TString groupNameTemp = "temp";
      TString groupNameTRatio = "tratio";

      const char *xTitle = "p_{T}/A (MeV/c)";
      const char *xVar = "p_{T}/A";
      auto yrheh = Form("R_{He-H} = #scale[0.8]{#frac{dM_{d}/d#Omegad(%s) #times dM_{he4}/d#Omegad(%s)}{dM_{t}/d#Omegad(%s) #times dM_{he3}/d#Omegad(%s)}}",xVar,xVar,xVar,xVar);
      variable frheh("rheh", "", "", titles(" ",xTitle,yrheh), binning(8,0,400), binning(100,0,5.));
      variable ftemp("temp", "", "", titles(" ",xTitle,"T = #frac{14.3}{log[1.59*R_{He-H}]}"), binning(8,0,400), binning(100,0,25));
      variable ftratio("tratio", "", "", titles(" ",xTitle,"T_{2} / T_{1}"), binning(8,0,400), binning(100,0,2));

      frheh.drawadd(groupNameR) -> legendlt();
      ftemp.drawadd(groupNameTemp) -> legendlt();
      ftratio.drawadd(groupNameTRatio);

      TH1D *histT[4] = {0};
      for (auto iSys : iSyss) {
        auto histR = (TH1D *) histMP[iSys][1] -> Clone();
        histR -> Multiply(histMP[iSys][4]);
        histR -> Divide(histMP[iSys][2]);
        histR -> Divide(histMP[iSys][3]);

        auto binnx = binning(histR);
        histT[iSys] = (TH1D *) histMP[iSys][1] -> Clone(TString("temperature")+iSys);
        histT[iSys] -> Reset();
        histT[iSys] -> SetMinimum(0);
        binnx.resetb();
        while (binnx.nextb()) {
          if (binnx.value>ptoaRange[1]) continue;
          auto rvalue = histR -> GetBinContent(binnx.idx);
          auto temperature = funcTemperature -> Eval(rvalue); 14.3 / (TMath::Log(1.59*rvalue));
          histT[iSys] -> Fill(binnx.value,temperature);
        }

        histR = (TH1D *) setSysAtt((TH1 *) histR, iSys);
        histT[iSys] = (TH1D *) setSysAtt((TH1 *) histT[iSys], iSys);

        ejungwoo::addsame(groupNameR,histR,"plhist gridx gridy",TString("R_{He-H}  sys")+sysBeam[iSys]);
        ejungwoo::addsame(groupNameTemp,histT[iSys],"plhist gridx gridy",TString("temp.  sys")+sysBeam[iSys]);
      }

      for (auto icomb : {0,1}) {
        int iSys1 = 0; int iSys2 = 1;
        if (icomb==1) { iSys1 = 3; iSys2 = 2; }

        if (histT[iSys1]!=nullptr||histT[iSys2]!=nullptr) {
          auto histTRatio = (TH1D *) histT[iSys1] -> Clone();
          histTRatio -> Divide(histT[iSys2]);
          histTRatio = (TH1D *) setCombAtt((TH1 *) histTRatio, icomb);
          ejungwoo::addsame(groupNameTRatio,histTRatio,"plhist gridx gridy",Form("T%d / T%d",sysBeam[iSys1],sysBeam[iSys2]));
        }
      }
    }

    break;
  }


  while (ana_e)
  {
    setpar("cut_all","$$(cut_PENS)");
    auto binnTheta = binning(40,0,80);
    auto binnKE = binning(40,0,400);
    auto vetheta = new variable("vetheta", "theta_lab*TMath::RadToDeg():(sqrt((p_lab*$$(z))*(p_lab*$$(z))+$$(m)*$$(m))-$$(m))/$$(a)", "$$(cut_all)", titles("", "KE_{Lab}/A (MeV)", "#theta_{Lab} (Deg.)"), binnKE,binnTheta);

    setpar("cut_alle","$$(cut_PENS)*$$(cut_tl)");
    auto vedistlab = new variable("vedistlab", "(sqrt((p_lab*$$(z))*(p_lab*$$(z))+$$(m)*$$(m))-$$(m))/$$(a)", "$$(cut_alle)", titles("", "KE_{Lab}/A (MeV)", ""), binning(40,0,400));
    auto f1Sinx = new TF1("satta","sin(x)",0,TMath::Pi());

    TH2D *histKET[nSyss][nPIDs];

    if (drawEDistTheta || drawETheta || drawER21)
    {
      for (auto iSys : iSyss)
      {
        setSystem(iSys);

        auto groupName = sysNames0[iSys]+"_ketheta";
        for (auto iPID : iPIDs) {
          setParticle(iPID);

          TString mainTitle = Form("mult=%d~%d,  sys%d,  %s",  mult0,  mult1,  sysBeam[iSys],  particleNames[iPID].Data());
          vetheta -> setmaint(mainTitle);
          histKET[iSys][iPID] = (TH2D *) vetheta -> draw(chain[iSys][iPID]);
          auto hist = (TH1 *) histKET[iSys][iPID];
          if (drawETheta) {
            ejungwoo::addnext(groupName, hist, "rangex colz");
            {
              auto graph = new TGraph();
              graph -> SetLineWidth(2);
              graph -> SetLineColor(kRed);
              graph -> SetPoint(0,0,0);
              graph -> SetPoint(1,400,0);
              graph -> SetPoint(2,400,20);
              graph -> SetPoint(3,0,20);
              graph -> SetPoint(4,0,0);
              ejungwoo::addsame(groupName, graph, "addx colorx l ");
            }
            {
              auto graph = new TGraph();
              graph -> SetLineWidth(2);
              graph -> SetLineColor(kRed);
              graph -> SetPoint(0,0,20);
              graph -> SetPoint(1,400,20);
              graph -> SetPoint(2,400,40);
              graph -> SetPoint(3,0,40);
              graph -> SetPoint(4,0,20);
              ejungwoo::addsame(groupName, graph, "addx colorx l ");
            }
            {
              auto graph = new TGraph();
              graph -> SetLineWidth(2);
              graph -> SetLineColor(kRed);
              graph -> SetPoint(0,0,40);
              graph -> SetPoint(1,400,40);
              graph -> SetPoint(2,400,60);
              graph -> SetPoint(3,0,60);
              graph -> SetPoint(4,0,40);
              ejungwoo::addsame(groupName, graph, "addx colorx l ");
            }
            {
              auto graph = new TGraph();
              graph -> SetLineWidth(2);
              graph -> SetLineColor(kRed);
              graph -> SetPoint(0,0,60);
              graph -> SetPoint(1,400,60);
              graph -> SetPoint(2,400,80);
              graph -> SetPoint(3,0,80);
              graph -> SetPoint(4,0,60);
              ejungwoo::addsame(groupName, graph, "addx colorx l ");
            }
          }
        }
      }
    } if (drawETheta==2) break;

    TH1D *histKEN[nSyss][nPIDs][numTheta];
    TH1D *histKE[nSyss][nPIDs][numTheta];

    if (drawEDistTheta || drawER21)
    {
      for (auto iSys : iSyss)
      {
        setSystem(iSys);

        for (auto iPID : iPIDs) {
          setParticle(iPID);

          auto numBins = binnTheta.n/numTheta;
          for (auto iTheta=0; iTheta<numTheta; ++iTheta) {
            double theta1 = TMath::DegToRad()*dTheta*(iTheta);
            double theta2 = TMath::DegToRad()*dTheta*(iTheta+1);
            auto solidAngleInPi = 2 * (dPhi/360) * f1Sinx -> Integral(theta1, theta2);
            auto solidAngle = solidAngleInPi * TMath::Pi();

            auto hist2 = histKET[iSys][iPID];
            auto hist0 = (TH1 *) hist2 -> ProjectionX(Form("ke_%d_%s_tt%d_",sysBeam[iSys],particleNames[iPID].Data(),iTheta),1+iTheta*numBins,(iTheta+1)*numBins);

            bool drawInLogY = false;

            if (1) {
              TString mainTitleN = Form("mult=%d~%d,  sys%d,  %s;KE_{Lab}/A (MeV);#frac{dN}{#Delta#Omega d(KE/A)}",  mult0,  mult1,  sysBeam[iSys],  particleNames[iPID].Data());
              auto hist1 = ejungwoo::mult_y(hist0,1./solidAngle);
              hist1 = setPIDAtt(hist1, iTheta, 0, 0, true,  mainTitleN);
              histKEN[iSys][iPID][iTheta] = (TH1D *) hist1;
              auto groupNameN = "ken_"+sysNames0[iSys];
              if (drawEDistTheta) {
                if (drawInLogY) {
                  hist1 -> SetMaximum(0.2);
                  hist1 -> SetMinimum(0.001);
                  if (iTheta==0) ejungwoo::addnext(groupNameN, hist1, "rangex logy colorx hist", particleNames[iPID] + Form(", #theta = %d~%d",int(dTheta*iTheta),int(dTheta*(iTheta+1))));
                  else           ejungwoo::addsame(groupNameN, hist1, "rangex logy colorx hist", particleNames[iPID] + Form(", #theta = %d~%d",int(dTheta*iTheta),int(dTheta*(iTheta+1))));
                } else {
                  hist1 -> SetMinimum(0);
                  hist1 -> SetMaximum(0.2);
                  if (iTheta==0) ejungwoo::addnext(groupNameN, hist1, "rangex colorx hist gridx gridy", particleNames[iPID] + Form(", #theta = %d~%d",int(dTheta*iTheta),int(dTheta*(iTheta+1))));
                  else           ejungwoo::addsame(groupNameN, hist1, "rangex colorx hist gridx gridy", particleNames[iPID] + Form(", #theta = %d~%d",int(dTheta*iTheta),int(dTheta*(iTheta+1))));
                }
              }
            }

            int rebinValue = 10;
            if (1) {
              TString mainTitleN = Form("mult=%d~%d,  sys%d,  %s;KE_{Lab}/A (MeV);#frac{dN}{d(KE/A)}",  mult0,  mult1,  sysBeam[iSys],  particleNames[iPID].Data());
              auto hist1 = setPIDAtt(hist0, iTheta, 0, 0, true,  mainTitleN);
              histKE[iSys][iPID][iTheta] = (TH1D *) hist1;
              //////////////////////////////////////////////////////////////////////////////
              histKE[iSys][iPID][iTheta] -> Rebin(rebinValue);
              //////////////////////////////////////////////////////////////////////////////
              auto groupNameN = "ke_"+sysNames0[iSys];
              if (drawEDistTheta) {
                if (drawInLogY) {
                  hist1 -> SetMaximum(0.2);
                  hist1 -> SetMinimum(0.001);
                  if (iTheta==0) ejungwoo::addnext(groupNameN, hist1, "rangex logy colorx hist", particleNames[iPID] + Form(", #theta = %d~%d",int(dTheta*iTheta),int(dTheta*(iTheta+1))));
                  else           ejungwoo::addsame(groupNameN, hist1, "rangex logy colorx hist", particleNames[iPID] + Form(", #theta = %d~%d",int(dTheta*iTheta),int(dTheta*(iTheta+1))));
                } else {
                  hist1 -> SetMinimum(0);
                  hist1 -> SetMaximum(0.05*rebinValue);
                  if (iTheta==0) ejungwoo::addnext(groupNameN, hist1, "rangex colorx gridx gridy", particleNames[iPID] + Form(", #theta = %d~%d",int(dTheta*iTheta),int(dTheta*(iTheta+1))));
                  else           ejungwoo::addsame(groupNameN, hist1, "rangex colorx gridx gridy", particleNames[iPID] + Form(", #theta = %d~%d",int(dTheta*iTheta),int(dTheta*(iTheta+1))));
                }
              }
            }

          }
        }
      }
    }

    if (drawER21) {
      auto fptoar21 = variable("fptoar21", "", "", titles(" ","","R21"), binnKE, binning(100,0.,2.));
      for (auto iSysComb : {0,1,2,3})
      {
        int iSys1, iSys2;
        if (iSysComb==0) { iSys1=0; iSys2=1; }
        else if (iSysComb==1) { iSys1=3; iSys2=2; }
        else if (iSysComb==2) { iSys1=0; iSys2=2; }
        else if (iSysComb==3) { iSys1=3; iSys2=1; }

        TString systga1 = sysNames1[iSys1];
        TString systga2 = sysNames1[iSys2];

        for (auto iTheta=0; iTheta<numTheta; ++iTheta)
        //for (auto iTheta=0; iTheta<1; ++iTheta)
        {
          double theta1 = TMath::DegToRad()*dTheta*(iTheta);
          double theta2 = TMath::DegToRad()*dTheta*(iTheta+1);

          TString ename = Form("r21_%d_%d",sysBeam[iSys1],sysBeam[iSys2]);
          TString mainTitle = Form("mult=%d~%d, %d/%d, #theta=%d~%d", mult0, mult1, sysBeam[iSys1], sysBeam[iSys2], int(dTheta*iTheta), int(dTheta*(iTheta+1)));
          TString xyTitle = Form(";KE_{Lab}/A (MeV);R_{21}( #frac{%s}{%s} )", sysNames1[iSys1].Data(), sysNames1[iSys2].Data());
          auto histFrame = fptoar21.draw();
          histFrame -> SetTitle(mainTitle+xyTitle);
          ejungwoo::addnext(ename,histFrame,"gridx gridy");
          for (int iPID : iPIDs) {
            double num_tracks_per_event_cut = 0.0001;
            cout << iSys1 << " " << iSys2 << " " << iTheta << "" << iPID << endl;
            auto hist = hist_is(iPID, histKE[iSys1][iPID][iTheta], histKE[iSys2][iPID][iTheta], num_tracks_per_event_cut);
            ejungwoo::addsame(ename, hist, "colorx EX0 pl", particleNames[iPID]);
            ejungwoo::addsame(ename, hist, "colorx histl addx");
          }
        }
      }
    } if (drawER21==2) break;


    break;
  }



  drawsaveall("cvsl","png");
}

TObjArray graph_is(int idx_particle, TH1 *hist132, TH1 *hist108, double num_tracks_cut)
{
  auto graph0 = new TGraphErrors();
  auto graph_v132 = new TGraphErrors();
  auto graph_v108 = new TGraphErrors();
  binning binn(hist132);

  for (auto bin=1; bin<=binn.n; ++bin)
  {
    auto binc = binn.getc(bin);
    auto v132 = hist132 -> GetBinContent(bin) * binn.w;
    auto v108 = hist108 -> GetBinContent(bin) * binn.w;
    auto idxvv = graph_v132 -> GetN();
    graph_v132 -> SetPoint(idxvv, binc, v132);
    graph_v108 -> SetPoint(idxvv, binc, v108);
    if (v132 < num_tracks_cut || v108 < num_tracks_cut)
      continue;
    auto idxis = graph0 -> GetN();
    graph0 -> SetPoint(idxis, binc, v132/v108);
  }

  TObjArray array;
  for (auto graph : {graph0, graph_v132, graph_v108}) {
    graph -> SetMarkerStyle(ejungwoo::markeri(idx_particle));
    graph -> SetLineColor(ejungwoo::colori(idx_particle));
    //graph -> SetLineWidth(2);
    array.Add(graph);
  }

  return array;
}

TH1D *hist_is(int idx_particle, TH1 *hist132, TH1 *hist108, double num_tracks_cut)
{
  binning binn(hist132);
  auto hist2 = (TH1D *) hist132 -> Clone();
  auto hist1 = (TH1D *) hist108 -> Clone();
  for (auto bin=1; bin<=binn.n; ++bin) {
    auto content = hist1 -> GetBinContent(bin);
    if (content==0) hist1 -> SetBinContent(bin,-1);
  }
  hist2 -> Divide(hist1);
  hist2 -> SetMarkerStyle(ejungwoo::markeri(idx_particle));
  hist2 -> SetMarkerColor(ejungwoo::colori(idx_particle));
  hist2 -> SetLineColor(ejungwoo::colori(idx_particle));

  return hist2;
}
