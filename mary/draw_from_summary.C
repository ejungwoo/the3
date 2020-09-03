#include "/Users/ejungwoo/config/ejungwoo.h"
#include "init_variables.h"
#include "functions.h"

//TString fVersions[] = {"kleft"};
TString fVersions[] = {"tleftmid"};

const int fIPIDs[] = {0,1,2,3,4};
//const int fIPIDs[] = {0,1};
//const int fIPIDs[] = {0};
const int fISystems[] = {0,1};
//const int fISystems[] = {0};

TObjArray draw_is(int idx_particle, TH1 *hist132, TH1 *hist108, double num_tracks_cut=0); 
TObjArray draw_ci(vector<TH1 *>hists, double ds, int isys);
TGraphErrors *draw_ratio(TGraphErrors *graph1, TGraphErrors *graph2, double max);

void draw_from_summary()
{
  ejungwoo::gcvspos(1300);
  ejungwoo::gstat(0);
  ejungwoo::gsave(0);
  //ejungwoo::gfast();
  ejungwoo::gsetvmark(0);
  //ejungwoo::gdummytp();

  bool anaY0 = 1;
  bool anaAll = 0;
  bool anaIS = 1;
  bool anaCI = 0;
  bool anaCorr = 0;
  bool anaAngle = 0;
  bool anaMult = 0;
  bool anaPid = 0;

  if (0) { anaIS = 1; anaCI = 1; anaCorr = 1; anaAngle = 1; anaMult = 1; anaPid = 1; anaAll = 1; ejungwoo::gdummytp(); }

  double num_tracks_per_event_cut = 0.01;
  int nbinsIS = 20;
  double dKECI = 5.;

  setpar("cut_x",                              "1.*($$(mult_cut))*($$(rap_cut))*($$(pes_cut))*($$(poz_cut))*($$(angle_cut))");
  setpar("cut_p",         "$$(prob)/$$(num_events)*($$(mult_cut))*($$(rap_cut))*($$(pes_cut))*($$(poz_cut))*($$(angle_cut))");
  setpar("cut_e",       "1./$$(eff)/$$(num_events)*($$(mult_cut))*($$(rap_cut))*($$(pes_cut))*($$(poz_cut))*($$(angle_cut))");
  setpar("cut_pe","$$(prob)/$$(eff)/$$(num_events)*($$(mult_cut))*($$(rap_cut))*($$(pes_cut))*($$(poz_cut))*($$(angle_cut))");

  //vector<variable*> varListAll = { &fvar_foby_cm, };
  vector<variable*> varListAll = {
    //&(fvar_prob+"y"),
    &fvar_eff, &fvar_sd,
    //&fvar_ptoa_cm,
    // &fvar_pt_cm, &fvar_fy_cm, &fvar_p_cm, &fvar_ke_cm,
    //&fvar_p_lab, &fvar_dedx,
    //&fvar_dpoca, &fvar_nr, &fvar_nl,
    &fvar_tta_cm, &fvar_phi_cm,
    //&fvar_tta_lab, &fvar_phi_lab,
    //&fvar_p_kab+&fvar_dedx,
    //&fvar_ny_cm,
    &fvar_foby_cm,
    &fvar_fy_cm,
    //fvar_foby_cm.add(fvar_eff),
    //fvar_poz_lab.add(fvar_eff),
  };

  auto fvar_pid = fvar_poz_lab+fvar_dedx+"z";

  /******************************************************************************************/

  auto filem = new TFile("data/particle_cutg_s20_p2.pidcut.root");
  TCutG *cutgs[6];
  TF1 *means[6];
  for (auto ipid : fIPIDs) {
    cutgs[ipid] = (TCutG *) filem -> Get(Form("CUTG%d",fParticlePDGs[ipid]));
    cutgs[ipid] -> SetLineColor(kRed);
    means[ipid] = (TF1 *) filem -> Get(Form("BEE%d",fParticlePDGs[ipid]));
    means[ipid] -> SetLineColor(kGray+2);
    means[ipid] -> SetLineStyle(1);
  }

  /******************************************************************************************/

  auto file = new TFile("data/dndy_kaneko.root");
  TH1D *histKP[2] = {0};
  TH1D *histKD[2] = {0};
  TH1D *histKT[2] = {0};
  histKP[0] = (TH1D *) file -> Get("h1dndy_132Sn_Proton");   histKP[0] -> SetMarkerStyle(24); histKP[0] -> SetMarkerColor(ejungwoo::colori(0));
  histKD[0] = (TH1D *) file -> Get("h1dndy_132Sn_Deuteron"); histKD[0] -> SetMarkerStyle(24); histKD[0] -> SetMarkerColor(ejungwoo::colori(1));
  histKT[0] = (TH1D *) file -> Get("h1dndy_132Sn_Triton");   histKT[0] -> SetMarkerStyle(24); histKT[0] -> SetMarkerColor(ejungwoo::colori(2));
  histKP[1] = (TH1D *) file -> Get("h1dndy_108Sn_Proton");   histKP[1] -> SetMarkerStyle(24); histKP[1] -> SetMarkerColor(ejungwoo::colori(0));
  histKD[1] = (TH1D *) file -> Get("h1dndy_108Sn_Deuteron"); histKD[1] -> SetMarkerStyle(24); histKD[1] -> SetMarkerColor(ejungwoo::colori(1));
  histKT[1] = (TH1D *) file -> Get("h1dndy_108Sn_Triton");   histKT[1] -> SetMarkerStyle(24); histKT[1] -> SetMarkerColor(ejungwoo::colori(2));

  /******************************************************************************************/

  for (auto aversion : fVersions)
  {
    auto condition_array = setversion(aversion);
    //vector<variable*> varListDndv   = { &fvar_foby_cm };
    vector<variable*> varListAngle  = { fvar_phi_lab.add(fvar_tta_lab), fvar_phi_cm.add(fvar_tta_cm), fvar_phi_lab.add(fvar_phi_cm), fvar_tta_lab.add(fvar_tta_cm), };
    vector<variable*> varListCorr   = { fvar_foby_cm.add(fvar_ptoa_cm), fvar_foby_cm.add(fvar_poz_lab), };
    vector<variable*> varListIS     = { &fvar_foby_cm };
    vector<variable*> varListCI     = { &fvar_keoa_cm };
    vector<variable*> varListN1     = { &fvar_n, &fvar_np, &fvar_nd, &fvar_nt, &fvar_nhe3, &fvar_nhe4, };

    TTree *tree_pid[4][5] = {0};
    TTree *tree_mult[4] = {0};
    
    for (auto isys : {0,1}) {
      setpar_syspid(isys);
      auto sys = fSystems[isys];
      TString fileName = Form("data/sys%d.%s.ana.particle.root",sys,ejungwoo::versioncc());
      auto file = new TFile(fileName);
      tree_mult[isys] = (TTree *) file -> Get("mult");
      for (int ipid : fIPIDAll)
        tree_pid[isys][ipid] = (TTree *) file -> Get(fParticleTreeNames[ipid]);
      fNumEvents[isys] = tree_mult[isys] -> GetEntries(ejungwoo::getpar("mult_cut2"));
      cout << fileName << " " << sys << " " << fNumEvents[isys] << endl;
    }

    for (auto isys : fISystems)
    {
      setpar_syspid(isys);
      //auto condition_array = setversion(version);
      auto legend_cond = new TLegend();
      for (auto condition : condition_array) {
        legend_cond -> AddEntry((TObject *)nullptr,ejungwoo::replace_parameters(condition),"");
        ejungwoo::make_l(legend_cond,-0.65);
      }

      auto sys = fSystems[isys];

      TTree *tree = tree_mult[isys];
      if (anaMult) {
        cout << "===================================== anaMult " << endl;
        setpar_syspid(isys);
        for (auto var : varListN1) { var -> drawaddnext(tree,TString("no")+sys,"hist"); }
      }

      if (anaAll) {
        cout << "===================================== anaAll " << endl;
        for (int ipid : fIPIDs) {
          setpar_syspid(isys,ipid);
          auto tree = tree_pid[isys][ipid];
          //TString ename = TString("all__")+sys;
          for (auto var : varListAll)
          {
            TString ename = TString("all__")+fParticleNames[ipid]+sys;

            var -> setCut("");
            auto integral = var -> drawaddnext(tree,ename,"hist","raw") -> fLastDrawing -> GetHistogram() -> Integral("width");

            var -> setCut("$$(cut_x)");
            auto integral2 = var -> drawaddnext(tree,ename,"hist","cut") -> fLastDrawing -> GetHistogram() -> Integral("width");

            var -> setCut("$$(cut_pe)");
            auto hist3 = ejungwoo::norm_integral(var -> draw(tree),integral2);
            ejungwoo::addsame(ename,hist3,"hist","cpe (scale-cut-i)");
          }

        }
      }

      if (anaY0) {
        bool hasAddedKanekoResult = true;
        cout << "===================================== anaY0 " << endl;
        for (auto var : {&fvar_foby_cm}) {
          for (int ipid : fIPIDs) {
            setpar_syspid(isys,ipid);
            auto tree = tree_pid[isys][ipid];

            TString ename = TString("dndv__")+var -> getName()+"__"+sys;
            auto hist1 = ejungwoo::dndx(var -> draw(tree));
            hist1 -> SetMarkerStyle(20);
            hist1 -> SetLineColor(ejungwoo::colori(ipid));
            hist1 -> SetLineWidth(2);
            hist1 -> SetMarkerColor(ejungwoo::colori(ipid));
            hist1 -> SetTitle(Form("%d;y0;dN/dy0",sys));
            hist1 -> SetMaximum(19);
            ejungwoo::addsame(ename, hist1, "hist rangex colorx gridx gridy", fParticleNames[ipid] + " " + sys);
            ejungwoo::addsame(ename, legend_cond, "addl");
            if (!hasAddedKanekoResult) {
              hasAddedKanekoResult = true;
              ejungwoo::addsame(ename, histKP[isys], "p colorx addx", fParticleNames[ipid] + " " + sys + " kaneko");
              ejungwoo::addsame(ename, histKD[isys], "p colorx addx", fParticleNames[ipid] + " " + sys + " kaneko");
              ejungwoo::addsame(ename, histKT[isys], "p colorx addx", fParticleNames[ipid] + " " + sys + " kaneko");
            }
          }
        }
      }

      if (anaAngle) {
        cout << "===================================== anaAngle " << endl;
        for (int ipid : fIPIDs) {
          auto tree = tree_pid[isys][ipid];
          setpar_syspid(isys,ipid);

          for (auto var : varListAngle)
            var -> drawadd(tree,var -> getName()+sys,ipid);
        }
      }

      if (anaCorr) {
        cout << "===================================== anaCorr " << endl;
        for (int ipid : fIPIDs) {
          auto tree = tree_pid[isys][ipid];
          setpar_syspid(isys,ipid);

          for (auto var : varListCorr)
            var -> drawadd(tree,var -> getName()+sys,ipid);
        }
      }

      if (anaPid) {
        cout << "===================================== anaPid " << endl;
        for (int ipid : fIPIDs) {
          auto tree = tree_pid[isys][ipid];
          setpar_syspid(isys,ipid);

          auto var = &fvar_pid;

          //var -> setCut(""); var -> drawaddhist(tree, TString("pid_raw_all")+sys, "logz");
          //var -> setCut(""); var -> drawaddnext(tree, TString("pid_raw_each")+sys, "logz");

          //var -> setCut("$$(cut_x)"); var -> drawaddhist(tree, TString("pid_cut_all")+sys, "logz");
          var -> setCut("$$(cut_x)"); var -> drawaddnext(tree, TString("pid_cut_all")+sys, "logz");

          //var -> setCut("$$(cut_pe"); var -> drawaddhist(tree, TString("pid_pecut_all")+sys, "logz");
          var -> setCut("$$(cut_pe)"); var -> drawaddnext(tree, TString("pid_pecut_all")+sys, "logz");
        }
      }
    }

    if (anaIS)
    {
      cout << "===================================== anaIS " << endl;
      for (auto ii : {0})
      {
        int isys1, isys2;
        if (ii==0) { isys1=0; isys2=1; }
        else       { isys1=3; isys2=2; }

        auto systgA1 = Form("%d+%d",fSystems[isys1],fTargetA[isys1]);
        auto systgA2 = Form("%d+%d",fSystems[isys2],fTargetA[isys2]);

        for (auto var : varListIS) {
          setpar("pname","all");
          TString cname_nn = Form("NN%d_",ii)+var->getName()+"__"+aversion;
          TString cname_is = Form("IS%d_",ii)+var->getName()+"__"+aversion;

          var->setn(nbinsIS);
          var->setmaint(Form("$$(pname) multiplicity %d vs %d",fSystems[isys1],fSystems[isys2]));
          auto name00 = var->getHistName();
          auto title00 = titles(Form("isoscaling plot %d/%d",fSystems[isys1],fSystems[isys2]),var->getTitle().x,Form("p(%s)/p(%s)",systgA1,systgA2));
          ejungwoo::add(cname_is,0,new_h(name00,title00,var->getBinn(),binning(100,0,2)),"addx");
          double ymax1 = 0, ymax2 = 0;
          for (int ipid : fIPIDs) {
            setpar_syspid(isys1,ipid); auto hist1 = var -> draw(tree_pid[isys1][ipid]);
            setpar_syspid(isys2,ipid); auto hist2 = var -> draw(tree_pid[isys2][ipid]);
            //ejungwoo::add(cname_nn,ipid,hist1,"hist",Form("%d",fSystems[isys1]));
            //ejungwoo::add(cname_nn,ipid,hist2,"hist",Form("%d",fSystems[isys2]));
            auto graphs = draw_is(ipid, hist1, hist2, num_tracks_per_event_cut);
            ejungwoo::add(cname_is, 0, graphs.At(0), "pl", fParticleNames[ipid]);
            //ejungwoo::add(cname_is, 1, graphs.At(1), "pl", fParticleNames[ipid]);
            //ejungwoo::add(cname_is, 2, graphs.At(2), "pl", fParticleNames[ipid]);
            double ymax1_ = ejungwoo::y2_g((TGraph *) graphs.At(1)); if (ymax1 < ymax1_) ymax1 = ymax1_;
            double ymax2_ = ejungwoo::y2_g((TGraph *) graphs.At(2)); if (ymax2 < ymax2_) ymax2 = ymax2_;
          }

          /*
          auto line1 = new TLine(var->getBinn().min,num_tracks_per_event_cut,var->getBinn().max,num_tracks_per_event_cut); line1 -> SetLineStyle(2);
          auto line2 = new TLine(var->getBinn().min,num_tracks_per_event_cut,var->getBinn().max,num_tracks_per_event_cut); line2 -> SetLineStyle(2);
          ejungwoo::add(cname_is,1,line1);
          ejungwoo::add(cname_is,2,line1);

          titles ttlIS132(systgA1,var->getTitle().x,fttly_ntne);
          titles ttlIS108(systgA2,var->getTitle().x,fttly_ntne);
          if (0) {
            ejungwoo::add(cname_is,1,new_h(var->getHistName()+Form("_frame%d",fSystems[isys1]),ttlIS132,var->getBinn(),binning(100,ymax1*0.0001,ymax1*5)),"addx rangex logy");
            ejungwoo::add(cname_is,2,new_h(var->getHistName()+Form("_frame%d",fSystems[isys2]),ttlIS108,var->getBinn(),binning(100,ymax2*0.0001,ymax2*5)),"addx rangex logy");
          } else {
            ejungwoo::add(cname_is,1,new_h(var->getHistName()+Form("_frame%d",fSystems[isys1]),ttlIS132,var->getBinn(),binning(100,0,ymax1*1.05)),"addx rangex");
            ejungwoo::add(cname_is,2,new_h(var->getHistName()+Form("_frame%d",fSystems[isys2]),ttlIS108,var->getBinn(),binning(100,0,ymax2*1.05)),"addx rangex");
          }
          */
        }
      }
    }

    if (anaCI)
    {
      cout << "===================================== anaCI " << endl;
      titles ttlCIM("","KE_{CM} (MeV)","#frac{dM_{n,CI}}{dKE_{CM} #Delta#Omega_{CM}}");
      titles ttlCIR("","KE_{CM} (MeV)","R = #frac{dM_{n,CI}}{dKE_{CM} #Delta#Omega_{CM}} / #frac{dM_{p,CI}}{dKE_{CM} #Delta#Omega_{CM}}");
      binning binnCIx(40,0,200);
      binning binnCIy(100,0.00005,1);
      binning binnCIy2(100,0.001,1);
      binning binnCIr(100,0,1.2);
      binning binnCIdr(100,0,2);

      for (auto var : varListCI) {
        var->setw(dKECI);
        var->setName(TString("ci_")+var->getName()); 
        var->setmaint("");
        TString cname_ci = var->getName()+"__"+aversion;
        ejungwoo::add(cname_ci, 0, new_h(var->getHistName()+"_cin_allframe", ttlCIM, binnCIx, binnCIy2), "addx rangex logy ");
        ejungwoo::add(cname_ci, 1, new_h(var->getHistName()+"_cip_allframe", ttlCIM, binnCIx, binnCIy2), "addx rangex logy ");
        ejungwoo::add(cname_ci, 2, new_h(var->getHistName()+"_cir_allframe", ttlCIR, binnCIx, binnCIr), "addx rangex ");

        TGraphErrors *graphs_rnp[4];
        for (auto isys : {0,1,2,3}) {
          auto sys = fSystems[isys];
          auto tga = fTargetA[isys];
          var->setmaint("($$(sys)) particle multiplicity");
          titles ttlCIMN("($$(sys)) CI neutrons","KE_{CM} (MeV)",Form("#frac{dM(%d)_{n,CI}}{dKE_{CM} #Delta#Omega_{CM}}",sys));
          titles ttlCIMP("($$(sys)) CI protons","KE_{CM} (MeV)",Form("#frac{dM(%d)_{p,CI}}{dKE_{CM} #Delta#Omega_{CM}}",sys));
          titles ttlCIRs("($$(sys)) CI n/p","KE_{CM} (MeV)",Form("R(%d) = #frac{dM(%d)_{n,CI}}{dKE_{CM} #Delta#Omega_{CM}} / #frac{dM(%d)_{p,CI}}{dKE_{CM} #Delta#Omega_{CM}}",sys,sys,sys));
          vector<TH1*> hists;
          TString cname_ci_sys = var->getName()+"_"+sys+"__"+aversion;
          for (auto ipid : fIPIDAll) {
          /****************************************************************************************************************/
            setpar_syspid(isys,ipid);
            auto hist = var -> draw(tree_pid[isys][ipid]);
            hists.push_back(hist);
            ejungwoo::add(cname_ci_sys, 1, hist, "histl logy grid", fParticleNames[ipid]);
          /****************************************************************************************************************/
          }
          setpar("pname","CI");
          auto graphs = draw_ci(hists, fSolidAngle, isys);
          /****************************************************************************************************************/
          ejungwoo::add(cname_ci_sys, 0, new_h(var->getHistName()+"_cin_frame", ttlCIMN, binnCIx, binnCIy), "addx rangex logy gridy");
          ejungwoo::add(cname_ci_sys, 0, graphs.At(0), "colorx p", "CI neutrons");
          /****************************************************************************************************************/
          ejungwoo::add(cname_ci_sys, 2, new_h(var->getHistName()+"_cip_frame", ttlCIMP, binnCIx, binnCIy), "addx rangex logy gridy");
          ejungwoo::add(cname_ci_sys, 2, graphs.At(1), "colorx p", "CI protons");
          /****************************************************************************************************************/
          ejungwoo::add(cname_ci_sys, 3, new_h(var->getHistName()+"_cir_frame", ttlCIRs, binnCIx, binnCIr), "addx rangex gridy");
          ejungwoo::add(cname_ci_sys, 3, graphs.At(2), "colorx p", "CI n/p");
          /****************************************************************************************************************/

          ejungwoo::add(cname_ci, 0, graphs.At(0), "colorx l", Form("(%d+%d) CI neutrons",sys,tga));
          ejungwoo::add(cname_ci, 1, graphs.At(1), "colorx l", Form("(%d+%d) CI protons",sys,tga));
          ejungwoo::add(cname_ci, 2, graphs.At(2), "colorx l", Form("(%d+%d) CI n/p",sys,tga));

          graphs_rnp[isys] = (TGraphErrors *) graphs.At(2);
        }

        titles ttlCID("Double ratio","KE_{CM} (MeV)","DR(n/p)_{#frac{132}{108}} = #frac{R(132)}{R(108)}");
        ejungwoo::add(var->getName()+"_r132o108__"+aversion, new_h(var->getHistName()+"_cidr_allframe", ttlCID, binnCIx, binnCIdr), "addx rangex gridx");
        ejungwoo::add(var->getName()+"_r132o108__"+aversion, draw_ratio(graphs_rnp[0], graphs_rnp[1], binnCIx.max), "p addx");
      }
    }
  }

  /******************************************************************************************/

  ejungwoo::drawsaveall("cvsl","png");
}
