//TGraphErrors *draw_is(int idx_particle, TH1 *hist132, int num_events_132, TH1 *hist108, int num_events_108, int num_tracks_cut=10000); 
TGraphErrors *draw_is(int idx_particle, TH1 *hist132, TH1 *hist108, double num_tracks_cut=0); 

void draw_from_summary()
{
  ejungwoo::gstat(0);

  bool ana_mult = true;
  bool ana_val = true;

  TString fParticleNames[6] = {"p","d","t","he3","he4","he6"};
  double fParticleMass[6] = {938.272, 1871.06, 2809.41, 2809.41, 3728.4, 5606.55};
  int fParticleA[6] = {1,2,3,3,4,6};
  int fParticleZ[6] = {1,1,1,2,2,2};

  ejungwoo::binning bmult  (100,0,100);
  ejungwoo::binning bmult2 (25,0,25);
  ejungwoo::binning bmult3 (40,0,40);
  ejungwoo::binning bkeoa  (40,0,400);
  ejungwoo::binning bptoa  (200,0,1200);
  ejungwoo::binning bny    (200,-1.5,3.0);
  ejungwoo::binning btta   (200,0,180);
  ejungwoo::binning bphi   (200,-180,180);
  ejungwoo::binning bpozlab(200,0,2500);
  ejungwoo::binning bdedx  (200,0,1500);

  ejungwoo::variable var_eff     ("eff","eff");
  ejungwoo::variable var_keoa_cm ("keoa_cm" ,"ke_cm/$$(a)"                ,"$$(gcut)" ,"$$(pname)_$$(sys);KE_{CM}/A (MeV);N_{tracks} / N_{events}"    ,bkeoa);
  ejungwoo::variable var_ptoa_cm ("ptoa_cm" ,"pt_cm/$$(a)"                ,"$$(gcut)" ,"$$(pname)_$$(sys);p_{T}/A (MeV/c);N_{tracks} / N_{events}"    ,bptoa);
  ejungwoo::variable var_ny_cm   ("ny_cm"   ,"ny_cm"                      ,"$$(gcut)" ,"$$(pname)_$$(sys);y_{CM}/y_{beam,CM};N_{tracks} / N_{events}" ,bny);
  ejungwoo::variable var_tta_cm  ("tta_cm"  ,"theta_cm*TMath::RadToDeg()" ,""         ,"$$(pname)_$$(sys);#theta_{CM} (Deg.)"                         ,btta);
  ejungwoo::variable var_phi_cm  ("phi_cm"  ,"phi_cm*TMath::RadToDeg()"   ,""         ,"$$(pname)_$$(sys);#phi_{CM} (Deg.)"                           ,bphi);
  ejungwoo::variable var_poz_lab ("poz_lab" ,"p_lab/$$(z)"                ,""         ,"$$(pname)_$$(sys);p_{Lab}/Z (MeV/c^{2})"                      ,bpozlab);
  ejungwoo::variable var_dedx    ("dedx"    ,"dedx"                       ,""         ,"$$(pname)_$$(sys);dE/dx"                                      ,bdedx);

  ejungwoo::variable var_n    ("n"    ,"n"    ,"" ,"all_$$(sys);N_{prob > 0.5}",bmult );
  ejungwoo::variable var_np   ("np"   ,"np"   ,"" ,  "p_$$(sys);N_{prob > 0.5}",bmult3);
  ejungwoo::variable var_nd   ("nd"   ,"nd"   ,"" ,  "d_$$(sys);N_{prob > 0.5}",bmult3);
  ejungwoo::variable var_nt   ("nt"   ,"nt"   ,"" ,  "t_$$(sys);N_{prob > 0.5}",bmult3);
  ejungwoo::variable var_nhe3 ("nhe3" ,"nhe3" ,"" ,"he3_$$(sys);N_{prob > 0.5}",bmult3);
  ejungwoo::variable var_nhe4 ("nhe4" ,"nhe4" ,"" ,"he4_$$(sys);N_{prob > 0.5}",bmult3);
  ejungwoo::variable var_ng    ("n_good"    ,"n_good"    ,"" ,"all_$$(sys);N_{prob > 0.8, eff > 0.1, abs(sd) < 3}",bmult2);
  ejungwoo::variable var_ngp   ("np_good"   ,"np_good"   ,"" ,  "p_$$(sys);N_{prob > 0.8, eff > 0.1, abs(sd) < 3}",bmult2);
  ejungwoo::variable var_ngd   ("nd_good"   ,"nd_good"   ,"" ,  "d_$$(sys);N_{prob > 0.8, eff > 0.1, abs(sd) < 3}",bmult2);
  ejungwoo::variable var_ngt   ("nt_good"   ,"nt_good"   ,"" ,  "t_$$(sys);N_{prob > 0.8, eff > 0.1, abs(sd) < 3}",bmult2);
  ejungwoo::variable var_nghe3 ("nhe3_good" ,"nhe3_good" ,"" ,"he3_$$(sys);N_{prob > 0.8, eff > 0.1, abs(sd) < 3}",bmult2);
  ejungwoo::variable var_nghe4 ("nhe4_good" ,"nhe4_good" ,"" ,"he4_$$(sys);N_{prob > 0.8, eff > 0.1, abs(sd) < 3}",bmult2);


  vector<ejungwoo::variable> varlist = { var_ny_cm+var_ptoa_cm, var_poz_lab+var_dedx, };
  vector<ejungwoo::variable> varlist2 = { var_n+var_ng, var_np+var_ngp, var_nd+var_ngd, var_nt+var_ngt, var_nhe3+var_nghe3, var_nhe4+var_nghe4, };
  vector<ejungwoo::variable> varlist3 = { var_np, var_nd, var_nt, var_nhe3, var_nhe4, };
  vector<ejungwoo::variable> varlist4 = { var_ngp, var_ngd, var_ngt, var_nghe3, var_nghe4, };

  ejungwoo::gsetcut(TCut("($$(tta_cm)<110)/$$(eff)/$$(num_events)"));
  //ejungwoo::gsetcut(TCut("($$(tta_cm)<110)"));

  //for (auto version : {"NewAna.2034.45b9400" ,"NewAna.2057.254c9ba"})
  //for (auto version : {"NewAna.2057.254c9ba"})
  for (auto version : {"NewAna.2034.45b9400"})
  {
    TString vshort = TString("v")+ejungwoo::tok(version,".",1);
    ejungwoo::gversion(version);

    TH1D *histKE[6][2] = {0};
    int num_events[2] = {0};

    for (auto sys : {132,108})
    {
      int isys = (sys==132?0:1);
      ejungwoo::setpar("sys",sys);

      //TString fileName = Form("data_summary/mult_all_%d.%s.root",sys,ejungwoo::versioncc());
      TString fileName = Form("data_summary2/summary_all_%d.%s.root",sys,ejungwoo::versioncc());
      auto file = new TFile(fileName);
      auto tree = (TTree *) file -> Get("mult");
      num_events[isys] = tree -> GetEntries();
      ejungwoo::setpar("num_events",num_events[isys]);
      cout << sys << " " << num_events[isys] << endl;

      if (ana_mult) {
        ejungwoo::gheader(vshort+"_"+sys+"_");

        TCanvas *cvs = ejungwoo::div(ejungwoo::cc((vshort+"_"+sys),1500,1000),3,2);

        int ivar = 0;
        for (auto var : varlist2) {
          auto hist = ejungwoo::tp(var,tree);
          TString cname = var.name+"_"+vshort+"_"+sys;
          ejungwoo::addto(cname,(cvs->cd(ivar+1)));
          ejungwoo::addto(cname, hist);
          ivar++;
        }
        for (auto var : varlist3) ejungwoo::addto(TString("n")+vshort+"_"+sys,ejungwoo::tp(var,tree),"hist");
        for (auto var : varlist4) ejungwoo::addto(TString("ng")+vshort+"_"+sys,ejungwoo::tp(var,tree),"hist");
      }

      if (ana_val)
      {
        TString fileName = Form("data_summary2/summary_all_%d.%s.root",sys,ejungwoo::versioncc());
        auto file = new TFile(fileName);

        int ivar = 0;
        TCanvas *cvs[10] = {0};
        for (auto var : varlist) {
          cvs[ivar] = ejungwoo::div(ejungwoo::cc((var.name+"_"+vshort+"_"+sys),1500,1000),3,2);
          ivar++;
        }

        for (int ipid : {0,1,2,3,4})
        {
          ejungwoo::setpar("pname",fParticleNames[ipid]);
          ejungwoo::gheader(fParticleNames[ipid]+"_"+vshort+"_"+sys+"_");
          ejungwoo::setpar("a",fParticleA[ipid]);
          ejungwoo::setpar("z",fParticleZ[ipid]);
          auto tree = (TTree *) file -> Get(fParticleNames[ipid]);

          for (auto var : {var_keoa_cm})  {
          //histKE[ipid][isys] = (TH1D *) (ejungwoo::drawv(var,tree)).hist;
            histKE[ipid][isys] = (TH1D *) ejungwoo::tp(var,tree);
          }

          int ivar = 0;
          for (auto var : varlist) {
            auto hist = ejungwoo::tp(var,tree);
            
            TString cname = var.name+"_"+vshort+"_"+sys+"_"+ipid;
            ejungwoo::addto(cname,(cvs[ivar]->cd(ipid+1)));
            ejungwoo::addto(cname, hist);
            ivar++;
          }
        }
      }
    }

    if (ana_val)
    {
      TString cname_ke = TString("KE_")+vshort;
      auto cvs = ejungwoo::div(ejungwoo::cc(cname_ke,1500,1000),3,2);
      TString cname_is = TString("IS_")+vshort;
      //ejungwoo::addto(cname_is,graphis, "pl", fParticleNames[ipid]);
      for (int ipid : {0,1,2,3,4})
      {
        ejungwoo::setpar("pname",fParticleNames[ipid]);
        ejungwoo::addto(cname_ke+"_"+ipid,(cvs->cd(ipid+1)));
        ejungwoo::addto(cname_ke+"_"+ipid,histKE[ipid][0],"hist","132");
        ejungwoo::addto(cname_ke+"_"+ipid,histKE[ipid][1],"hist","108");

        //auto graphis = draw_is(ipid, histKE[ipid][0], num_events[0], histKE[ipid][1], num_events[1]);
        auto graphis = draw_is(ipid, histKE[ipid][0], histKE[ipid][1],0.05);
        ejungwoo::addto(cname_is,graphis, "pl", fParticleNames[ipid]);
      }
    }
  }

  ejungwoo::drawall();
}

//TGraphErrors *draw_is(int idx_particle, TH1 *hist132, int num_events_132, TH1 *hist108, int num_events_108, int num_tracks_cut)
TGraphErrors *draw_is(int idx_particle, TH1 *hist132, TH1 *hist108, double num_tracks_cut)
{
  auto graph_is = new TGraphErrors();
  ejungwoo::binning binn(hist132);

  for (auto bin=1; bin<=binn.n; ++bin)
  {
    auto v132 = hist132 -> GetBinContent(bin);
    auto v108 = hist108 -> GetBinContent(bin);

    auto num_tracks_132 = /*num_events_132*/v132;
    auto num_tracks_108 = /*num_events_108*/v108;

    //cout << num_tracks_132 << " " << num_tracks_108 << endl;

    if (num_tracks_132 < num_tracks_cut || num_tracks_108 < num_tracks_cut)
      continue;

    Double_t ratio = v132/v108;
    auto idx_point = graph_is -> GetN();
    auto bincenter = hist132 -> GetXaxis() -> GetBinCenter(bin);
    graph_is -> SetPoint(idx_point, bincenter, ratio);
    //graph_is -> SetPointError(idx_point, binn.w/2, 0);
  }

  graph_is -> SetMarkerStyle(ejungwoo::markeri(idx_particle));
  graph_is -> SetMarkerSize(1.2);
  graph_is -> SetLineColor(ejungwoo::colori(idx_particle));
  graph_is -> SetLineWidth(2);

  return graph_is;
}
