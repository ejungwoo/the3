#include "init_variables.h"

void draw_all()
{
  gStyle -> SetOptStat(0);

  auto cvsOvv = new ecanvas("cvsOvv");
  auto cvsPID = new ecanvas("pid_p05");
  auto cvsPoz = new ecanvas("cvsPoz");

  //cvsOvv -> off();
  //cvsPID -> off();
  //cvsPoz -> off();

  cvsOvv -> setPadSizeSingle(1500,1000);
  auto h1_prob      = cvsOvv -> histNext("prob",      "[pname], [sys], [tag], [multttl] (raw);PID probability",         "prob",            "[cut_raw]", 100, 0, 1,      "stats=0,legend=0,logy=1");
  auto h1_eff       = cvsOvv -> histNext("eff",       "[pname], [sys], [tag], [multttl] (raw);Track efficiency",        "eff",             "[cut_raw]", 100, 0, 1,      "stats=0,legend=0,logy=1");
  auto h1_sd        = cvsOvv -> histNext("sd",        "[pname], [sys], [tag], [multttl] (raw);PID sigma-dist(#sigma)",  "sd",              "[cut_raw]", 100, -5, 5,     "stats=0,legend=0,draw=hist");

  auto h1_y0        = cvsOvv -> histNext("y0",        "[pname], [sys], [tag], [multttl];Rapidity y_{0}",    "fy_cm/(by_cm/2)",             "[cut_epc]", 100, -1, 2,     "stats=0,legend=0,draw=hist");
  auto h1_ptoa      = cvsOvv -> histNext("ptoa",      "[pname], [sys], [tag], [multttl];p_{T}/A (MeV/c)",   "pt_cm/[a]",                   "[cut_epc]", 100, 0, 1000,   "stats=0,legend=0,draw=hist");
  auto h1_ke_cm     = cvsOvv -> histNext("ke_cm",     "[pname], [sys], [tag], [multttl];KE_{CM} (MeV)",     "ke_cm",                       "[cut_epc]", 100, 0, 500,    "stats=0,legend=0,draw=hist");
  auto h1_poz       = cvsOvv -> histNext("poz",       "[pname], [sys], [tag], [multttl];p_{Lab}/Z (MeV/c)", "p_lab",                       "[cut_epc]", 400, 0, 2500,   "stats=0,legend=0,draw=hist");
  auto h1_dedx      = cvsOvv -> histNext("dedx",      "[pname], [sys], [tag], [multttl];dE/dx",             "dedx",                        "[cut_epc]", 400, 0, 1200,   "stats=0,legend=0,draw=hist");
  auto h1_theta_cm  = cvsOvv -> histNext("theta_cm",  "[pname], [sys], [tag], [multttl];#theta_{CM}",       "TMath::RadToDeg()*theta_cm",  "[cut_epc]", 100, 0, 180,    "stats=0,legend=0,draw=hist");
  auto h1_phi_cm    = cvsOvv -> histNext("phi_cm",    "[pname], [sys], [tag], [multttl];#phi_{CM}",         "TMath::RadToDeg()*phi_cm",    "[cut_epc]", 100, -180, 180, "stats=0,legend=0,draw=hist");
  auto h1_theta_lab = cvsOvv -> histNext("theta_lab", "[pname], [sys], [tag], [multttl];#theta_{Lab}",      "TMath::RadToDeg()*theta_lab", "[cut_epc]", 100, 0, 180,    "stats=0,legend=0,draw=hist");
  auto h1_phi_lab   = cvsOvv -> histNext("phi_lab",   "[pname], [sys], [tag], [multttl];#phi_{Lab}",        "TMath::RadToDeg()*phi_lab",   "[cut_epc]", 100, -180, 180, "stats=0,legend=0,draw=hist");

  auto h2_pid       = new ehist("pid",  h1_poz, h1_dedx, "","stats=0,legend=0,logz=1,zmax=0.05,zmin=0.00005,colz");
  auto h2_pid1      = new ehist("pid1", h1_poz, h1_dedx, "[cut_tta1]","stats=0,legend=0,logz=1,colz");
  auto h2_pid2      = new ehist("pid2", h1_poz, h1_dedx, "[cut_tta2]","stats=0,legend=0,logz=1,colz");
  auto h2_pid3      = new ehist("pid3", h1_poz, h1_dedx, "[cut_tta3]","stats=0,legend=0,logz=1,colz");
  auto h2_pid4      = new ehist("pid4", h1_poz, h1_dedx, "[cut_tta4]","stats=0,legend=0,logz=1,colz");
  auto h2_ypt       = new ehist("ypt",  h1_y0,  h1_ptoa, "","stats=0,legend=0,colz");

  cvsOvv -> addNext(h2_pid);
  cvsOvv -> addNext(h2_pid1);
  cvsOvv -> addNext(h2_pid2);
  cvsOvv -> addNext(h2_pid3);
  cvsOvv -> addNext(h2_pid4);
  cvsOvv -> addNext(h2_ypt);

  cvsOvv -> setPar("cut_tta1", "theta_lab> 0*TMath::DegToRad() && theta_lab<20*TMath::DegToRad()");
  cvsOvv -> setPar("cut_tta2", "theta_lab>20*TMath::DegToRad() && theta_lab<40*TMath::DegToRad()");
  cvsOvv -> setPar("cut_tta3", "theta_lab>40*TMath::DegToRad() && theta_lab<60*TMath::DegToRad()");
  cvsOvv -> setPar("cut_tta4", "theta_lab>60*TMath::DegToRad() && theta_lab<80*TMath::DegToRad()");

  cvsPID -> setPadSizeSingle(1500,1000);

  //for (auto iParticle : fParticleIdx) {
  /*
  for (auto iParticle : {0}) {
    const char *name = Form("pid_guideline_%d",iParticle);
    auto fileName = Form("data/%s.root",name);
    auto file = new TFile(fileName);
    auto graph = (TGraph *) file -> Get(name);
    cvsPID -> addToAll(graph,"legend=0,draw=samel");
  }
  */

  TString cut_all = "prob/eff/[nevents]*(prob>.1&&eff>.05&&abs(sd)<3)";

  cvsPoz -> fixPar("pname", "p");
  auto h1_poz2 = cvsPoz -> histSame("poz2", "[pname], [sys], [tag], [multttl];p_{Lab}/Z (MeV/c)", "p_lab",   "[cut_raw]", 100, 0, 1500,   "stats=0,legend=1,gridx=1,gridy=1,draw=hist");
  auto h1_poz3 = cvsPoz -> histSame("poz3", "[pname], [sys], [tag], [multttl];p_{Lab}/Z (MeV/c)", "p_lab",   "[cut_raw]", 100, 0, 1500,   "stats=0,legend=1,gridx=1,gridy=1,draw=hist");
  //h1_poz2 -> addTask("tfitGaus",.5,"","legend=1,gausInfo=1,gausLegend=1,draw=samel");
  //h1_poz3 -> addTask("tfitGaus",.5,"","legend=1,gausInfo=1,gausLegend=1,draw=samel");

  //for (auto iTag : fTagIdx2)
  for (auto iTag : {0})
  {
    const char *fileTag = fTags[iTag];

    //for (auto iMult : fMultOptions)
    for (auto iMult : {0})
    {
      Long64_t numEvents[4] = {0};
      //for (auto iSys : fSystemIdx) {
      for (auto iSys : {0,1}) {
        auto system = fSystems[iSys];
        auto fileNameMult = Form("data2/sys%d_%s_%d_%d.%s.ana.%s.root", system, fileTag, fMultLL[iMult], fMultHL[iMult], fSpiritVersion, "mult");
        auto fileMult = new TFile(fileNameMult,"read");
        auto treeMult = (TTree *) fileMult -> Get("mult");
        numEvents[iSys] = treeMult -> GetEntries("1");
        cout_info << "Number of events in " << fileMult->GetName() << " : " << numEvents[iSys] << endl;
      }

      //for (auto iSys : fSystemIdx)
      for (auto iSys : {0})
      {
        auto system = fSystems[iSys];
        auto anaTagSystem = fTagName[iTag] + fMultName[iMult] + system;

        for (auto cvs : {cvsOvv,cvsPID, cvsPoz}) {
          cvs -> setPad(0);
          cvs -> setTag(anaTagSystem);
          cvs -> setPar("cut_raw","");
          cvs -> setPar("cut_epc",cut_all);
          cvs -> setPar("nevents", numEvents[iSys]);
          cvs -> setPar("sys", TString::Itoa(system,10));
          cvs -> setPar("pname", "all");
          cvs -> setPar("tag", fTagName2[iTag]);
          cvs -> setPar("multttl", fMultName3[iMult]);
        }

        //auto fileNameData = Form("data/sys%d_%s_%d_%d.%s.ana.%s.root", system, fileTag, fMultLL[iMult], fMultHL[iMult], fSpiritVersion, fParticleNames[iParticle]);
        //auto fileData = new TFile(fileNameData,"read");
        //auto fileData = (TTree *) fileData -> Get(fParticleNames[iParticle].Data());

        //cvsPID -> draw(treeP05);
        //cvsPoz -> draw(treeP05);

        auto cvsPID2 = new ecanvas(Form("pid%d_epc",system));
        cvsPID2 -> setPadSizeSingle(1500,1000);
        cvsPID2 -> setTag(anaTagSystem);

        for (auto iParticle : fParticleIdx)
        {
          auto particleName = fParticleNames[iParticle];
          auto anaTagParticle = fTagName[iTag] + fMultName[iMult] + system + particleName;

          auto fileNameParticle = Form("data2/sys%d_%s_%d_%d.%s.ana.%s.root", system, fileTag, fMultLL[iMult], fMultHL[iMult], fSpiritVersion, fParticleNames[iParticle].Data());
          auto fileParticle = new TFile(fileNameParticle,"read");
          auto treeParticle = (TTree *) fileParticle -> Get(fParticleNames[iParticle].Data());

          cvsOvv -> setPar("a", fParticleA[iParticle]);
          cvsOvv -> setPar("z", fParticleZ[iParticle]);
          cvsOvv -> setPar("pname", particleName);

          cvsOvv -> setPad(0);
          cvsOvv -> setTag(anaTagParticle);

          cvsOvv -> draw(treeParticle);
          return;

          if (iParticle==0) {
                 if (iTag==0) { h1_poz2 -> make(treeParticle); h1_poz2 -> getHist() -> SetLineColor(kRed);}
            else if (iTag==1) { h1_poz3 -> make(treeParticle); }
          }
        }

        /*
        for (auto iParticle : fParticleIdx) {
          const char *name = Form("pid_guideline_%d",iParticle);
          auto fileNameGuide = Form("data/%s.root",name);
          auto file = new TFile(fileNameGuide);
          auto graph = (TGraph *) file -> Get(name);
          cvsPID2 -> addToAll(graph,"legend=0,draw=samel");
        }
        */

        //cvsPID2 -> draw();
        //cvsPID2 -> save("png","figures");
      }
    }
  }

  //cvsPoz -> draw();
  //cvsPoz -> save("png","figures");
}

//#ifdef ALIJLKBJHSWOIEITOQOPSKLFJSDFQW
//  setpar("cut_PENS",Form("prob/eff/$$(nevents)*(prob>%.2f&&eff>%.2f&&abs(sd)<%.2f)",probLL,effLL,sdHL));
//  setpar("cut_ptoa",Form("($$(ptoa)>%.2f&&$$(ptoa)<%.2f)",ptoaRange[0],ptoaRange[1]));
//  setpar("cut_y0",Form("($$(y0)>%.2f&&$$(y0)<%.2f)",y0Range[0],y0Range[1]));
//  setpar("cut_tl",Form("(theta_lab>%.4f&&theta_lab<%.4f)",TMath::DegToRad()*thetaLabRange[0],TMath::DegToRad()*thetaLabRange[1]));
//  setpar("cut_tc",Form("(theta_cm>%.4f&&theta_cm<%.4f)",TMath::DegToRad()*thetaCMRange[0],TMath::DegToRad()*thetaCMRange[1]));
//  setpar("cut_phic",Form("(phi_cm>%.4f&&phi_cm<%.4f)",TMath::DegToRad()*phiCMRange[0],TMath::DegToRad()*phiCMRange[1]));
//  setpar("cut_dedx",Form("(dedx>%.4f&&dedx<%.4f)",dedx1,dedx2));
//#endif
