TObjArray draw_is(int idx_particle, TH1 *hist132, TH1 *hist108, double num_tracks_cut=0); 

void draw_from_summary()
{
  ejungwoo::gstat(0);

  bool just_is = 1;
  bool ana_mult = 1;
  bool ana_val = 1;

  double num_tracks_per_event_cut = 0.02;

  vector<TString> versionlist = { "NewAna.2070.0231288", };

  if (just_is) { ana_mult = 0; ana_val = 0; }

  TString fParticleTreeNames[6] = {"p","d","t","he3","he4","he6"};
  TString fParticleNames[6] = {"ptn","dtn","ttn","he3","he4","he6"};
  double fParticleMass[6] = {938.272, 1871.06, 2809.41, 2809.41, 3728.4, 5606.55};
  int fParticleA[6] = {1,2,3,3,4,6};
  int fParticleZ[6] = {1,1,1,2,2,2};

  ejungwoo::binning bmult  (100,0,100);
  ejungwoo::binning bmult2 (25,0,25);
  ejungwoo::binning bmult3 (40,0,40);
  ejungwoo::binning bkeoa  (40,0,400);
  ejungwoo::binning bptoa  (200,0,1200);
  ejungwoo::binning bny    (200,-2.,2.0);
  ejungwoo::binning btta   (200,0,180);
  ejungwoo::binning bphi   (200,-180,180);
  ejungwoo::binning bpozlab(200,0,2500);
  ejungwoo::binning bdedx  (200,0,1500);

  ejungwoo::variable("prob","prob");
  ejungwoo::variable("eff","eff");

  ejungwoo::variable var_keoa_cm ("keoa_cm"   ,"ke_cm/$$(a)"                ,"$$(cuti)" ,"$$(pname)_$$(sys);KE_{CM}/A (MeV);N_{tracks} / N_{events}"    ,bkeoa);
  ejungwoo::variable var_ptoa_cm ("ptoa_cm"   ,"pt_cm/$$(a)"                ,"$$(cuti)" ,"$$(pname)_$$(sys);p_{T}/A (MeV/c);N_{tracks} / N_{events}"    ,bptoa);
  ejungwoo::variable var_ny_cm   ("ny_cm"     ,"ny_cm"                      ,"$$(cuti)" ,"$$(pname)_$$(sys);y_{CM}/y_{beam,CM};N_{tracks} / N_{events}" ,bny);
  ejungwoo::variable var_tta_cm  ("tta_cm"    ,"theta_cm*TMath::RadToDeg()" ,"$$(cuta)" ,"$$(pname)_$$(sys);#theta_{CM} (Deg.)"                         ,btta);
  ejungwoo::variable var_phi_cm  ("phi_cm"    ,"phi_cm*TMath::RadToDeg()"   ,"$$(cuta)" ,"$$(pname)_$$(sys);#phi_{CM} (Deg.)"                           ,bphi);
  ejungwoo::variable var_poz_lab ("poz_lab"   ,"p_lab/$$(z)"                ,"$$(cuta)" ,"$$(pname)_$$(sys);p_{Lab}/Z (MeV/c^{2})"                      ,bpozlab);
  ejungwoo::variable var_dedx    ("dedx"      ,"dedx"                       ,"$$(cuta)" ,"$$(pname)_$$(sys);dE/dx"                                      ,bdedx);

  ejungwoo::variable var_n       ("n"         ,"np+nd+nt+nhe3+nhe4"         ,"$$(cutn)" ,"all_$$(sys);N_{prob > 0.5}" ,bmult);
  ejungwoo::variable var_np      ("np"        ,"np"                         ,"$$(cutn)" ,  "p_$$(sys);N_{prob > 0.5}" ,bmult3);
  ejungwoo::variable var_nd      ("nd"        ,"nd"                         ,"$$(cutn)" ,  "d_$$(sys);N_{prob > 0.5}" ,bmult3);
  ejungwoo::variable var_nt      ("nt"        ,"nt"                         ,"$$(cutn)" ,  "t_$$(sys);N_{prob > 0.5}" ,bmult3);
  ejungwoo::variable var_nhe3    ("nhe3"      ,"nhe3"                       ,"$$(cutn)" ,"he3_$$(sys);N_{prob > 0.5}" ,bmult3);
  ejungwoo::variable var_nhe4    ("nhe4"      ,"nhe4"                       ,"$$(cutn)" ,"he4_$$(sys);N_{prob > 0.5}" ,bmult3);

  ejungwoo::variable var_ng      ("n_good"    ,"n_good"    ,"$$(cutn)"      ,"all_$$(sys);N_{prob > $$(prob_cut), eff > $$(eff_cut), abs(sd) < $$(sd_cut)}" ,bmult3);
  ejungwoo::variable var_ngp     ("np_good"   ,"np_good"   ,"$$(cutn)"      ,  "p_$$(sys);N_{prob > $$(prob_cut), eff > $$(eff_cut), abs(sd) < $$(sd_cut)}" ,bmult2);
  ejungwoo::variable var_ngd     ("nd_good"   ,"nd_good"   ,"$$(cutn)"      ,  "d_$$(sys);N_{prob > $$(prob_cut), eff > $$(eff_cut), abs(sd) < $$(sd_cut)}" ,bmult2);
  ejungwoo::variable var_ngt     ("nt_good"   ,"nt_good"   ,"$$(cutn)"      ,  "t_$$(sys);N_{prob > $$(prob_cut), eff > $$(eff_cut), abs(sd) < $$(sd_cut)}" ,bmult2);
  ejungwoo::variable var_nghe3   ("nhe3_good" ,"nhe3_good" ,"$$(cutn)"      ,"he3_$$(sys);N_{prob > $$(prob_cut), eff > $$(eff_cut), abs(sd) < $$(sd_cut)}" ,bmult2);
  ejungwoo::variable var_nghe4   ("nhe4_good" ,"nhe4_good" ,"$$(cutn)"      ,"he4_$$(sys);N_{prob > $$(prob_cut), eff > $$(eff_cut), abs(sd) < $$(sd_cut)}" ,bmult2);

  //vector<ejungwoo::variable> varListIS = { var_ny_cm, var_keoa_cm, var_ptoa_cm };
  vector<ejungwoo::variable> varListIS = { var_ny_cm };
  vector<ejungwoo::variable> varListAll = { var_phi_cm+var_tta_cm, var_ny_cm+var_ptoa_cm, var_poz_lab+var_dedx, };
  vector<ejungwoo::variable> varListN = { var_n+var_ng, var_np+var_ngp, var_nd+var_ngd, var_nt+var_ngt, var_nhe3+var_nghe3, var_nhe4+var_nghe4, };
  vector<ejungwoo::variable> varListN1 = { var_np, var_nd, var_nt, var_nhe3, var_nhe4, };
  vector<ejungwoo::variable> varListN2 = { var_ngp, var_ngd, var_ngt, var_nghe3, var_nghe4, };

  //ejungwoo::setpar("cuti","1.*$$(prob)/$$(eff)/$$(num_events)*(($$(tta_cm)<110)&&($$(phi_cm)<-100||$$(phi_cm)>100))");
  //ejungwoo::setpar("cuti","1.*$$(prob)/$$(eff)/$$(num_events)*($$(phi_cm)<-100||$$(phi_cm)>100)");
  //ejungwoo::setpar("cuti","1.*$$(prob)/$$(eff)/$$(num_events)*($$(tta_cm)<110)");
  //ejungwoo::setpar("cuti","1.*$$(prob)/$$(eff)/$$(num_events)");
  ejungwoo::setpar("cuti","1.*$$(prob)/$$(eff)/$$(num_events)");
  ejungwoo::setpar("cuta","$$(cuti)");
  ejungwoo::setpar("cutn","");

  for (auto version : versionlist)
  {
    TString vshort = TString("v")+ejungwoo::tok(version,".",1);
    ejungwoo::gversion(version);

    TH1D *histKE[4][6][2] = {0};
    int num_events[2] = {0};

    for (auto sys : {132,108})
    {
      int isys = (sys==132?0:1);
      ejungwoo::setpar("sys",sys);

      TString fileName = Form("data__%s/sys%d.%s.ana.particle.root",ejungwoo::versioncc(),sys,ejungwoo::versioncc());
      auto file = new TFile(fileName);
      auto tree_mult = (TTree *) file -> Get("mult");
      TTree *tree_pid[5] = {0};
      for (int ipid : {0,1,2,3,4})
        tree_pid[ipid] = (TTree *) file -> Get(fParticleTreeNames[ipid]);

      ejungwoo::setpar("prob_cut", ((TParameter<double> *) file -> Get("prob_cut")) -> GetVal());
      ejungwoo::setpar("eff_cut", ((TParameter<double> *) file -> Get("eff_cut")) -> GetVal());
      ejungwoo::setpar("sd_cut", ((TParameter<double> *) file -> Get("sd_cut")) -> GetVal());

      TTree *tree = tree_mult;
      num_events[isys] = tree_mult -> GetEntries();
      ejungwoo::setpar("num_events",num_events[isys]);
      cout << sys << " " << num_events[isys] << endl;

      if (ana_mult) {
        ejungwoo::gheader(vshort+"_"+sys+"_");

        TPad *cvs = nullptr;
        if (varListN.size()>0)
          cvs = ejungwoo::div(ejungwoo::cc(TString("mult"),1500,1000),3,2);

        int ivar = 0;
        for (auto var : varListN) {
          auto hist = ejungwoo::tp(var,tree);
          TString cname = var.name;
          ejungwoo::addto(cname,(cvs->cd(ivar+1)));
          ejungwoo::addto(cname, hist);
          ivar++;
        }
        for (auto var : varListN1) ejungwoo::addto(TString("n"),ejungwoo::tp(var,tree),"hist");
        for (auto var : varListN2) ejungwoo::addto(TString("ng"),ejungwoo::tp(var,tree),"hist");
      }

      if (ana_val)
      {
        int ivar = 0;
        TPad *cvs[10] = {0};
        for (auto var : varListAll) {
          cvs[ivar] = ejungwoo::div(ejungwoo::cc((var.name),1500,1000),3,2);
          ivar++;
        }

        for (int ipid : {0,1,2,3,4})
        {
          auto tree = tree_pid[ipid];
          ejungwoo::setpar("pname",fParticleNames[ipid]);
          ejungwoo::gheader(fParticleNames[ipid]+"_"+vshort+"_"+sys+"_");
          ejungwoo::setpar("a",fParticleA[ipid]);
          ejungwoo::setpar("z",fParticleZ[ipid]);

          int ivar = 0;
          for (auto var : varListAll) {
            auto hist = ejungwoo::tp(var,tree);
            
            TString cname = var.name+"_"+ipid;
            ejungwoo::addto(cname,(cvs[ivar]->cd(ipid+1)));
            ejungwoo::addto(cname, hist);
            ivar++;
          }
        }
      }

      if (just_is || ana_val)
      {
        for (int ipid : {0,1,2,3,4}) {
          auto tree = tree_pid[ipid];
          ejungwoo::setpar("pname",fParticleNames[ipid]);
          ejungwoo::gheader(fParticleNames[ipid]+"_"+vshort+"_"+sys+"_");
          ejungwoo::setpar("a",fParticleA[ipid]);
          ejungwoo::setpar("z",fParticleZ[ipid]);

          int ivar = 0;
          for (auto var : varListIS) {
            var.binn.n=50;
            histKE[ivar][ipid][isys] = (TH1D *) ejungwoo::tp(var,tree);
            ivar++;
          }
        }
      }
    }

    if (just_is || ana_val)
    {
      ejungwoo::gheader("");

      int ivar = 0;
      for (auto var : varListIS)
      {
        ejungwoo::setpar("pname","all");
        TString cname_is = TString("IS_")+var.name+"__"+vshort;
        cout << cname_is << endl;
        auto cvs = ejungwoo::div(ejungwoo::cc(cname_is+"_div",1500,1000),3,2);

        var.title.y = "p(132+124)/p(108+112)";
        ejungwoo::addto(cname_is,var.new_h(),"addx");

        for (int ipid : {0,1,2,3,4})
        {
          ejungwoo::setpar("pname",fParticleNames[ipid]);
          ejungwoo::addto(cname_is+"_"+ipid,(cvs->cd(ipid+1)));
          ejungwoo::addto(cname_is+"_"+ipid,histKE[ivar][ipid][0],"hist","132");
          ejungwoo::addto(cname_is+"_"+ipid,histKE[ivar][ipid][1],"hist","108");

          auto graphs = draw_is(ipid, histKE[ivar][ipid][0], histKE[ivar][ipid][1],num_tracks_per_event_cut);
          auto graphis = (TGraphErrors *) graphs.At(0);
          ejungwoo::addto(cname_is, graphis, "pl", fParticleNames[ipid]);
        }
        ivar++;
      }
    }
  }

  ejungwoo::drawall();
}

TObjArray draw_is(int idx_particle, TH1 *hist132, TH1 *hist108, double num_tracks_cut)
{
  auto graph_is = new TGraphErrors();
  auto graph_v132 = new TGraphErrors();
  auto graph_v108 = new TGraphErrors();
  ejungwoo::binning binn(hist132);

  for (auto bin=1; bin<=binn.n; ++bin)
  {
    auto binc = binn.getc(bin);
    auto v132 = hist132 -> GetBinContent(bin);
    auto v108 = hist108 -> GetBinContent(bin);

    auto idxvv = graph_v132 -> GetN();
    graph_v132 -> SetPoint(idxvv, binc, v132);
    graph_v108 -> SetPoint(idxvv, binc, v108);

    if (v132 < num_tracks_cut || v108 < num_tracks_cut)
      continue;

    auto idxis = graph_is -> GetN();
    graph_is -> SetPoint(idxis, binc, v132/v108);
  }

  TObjArray array;
  for (auto graph : {graph_is, graph_v132, graph_v108}) {
    graph -> SetMarkerStyle(ejungwoo::markeri(idx_particle));
    //graph -> SetMarkerSize(1.2);
    graph -> SetLineColor(ejungwoo::colori(idx_particle));
    graph -> SetLineWidth(2);
    array.Add(graph);
  }

  return array;
}
