TGraphErrors *draw_is(int idx_particle, TH1 *hist132, int num_events_132, TH1 *hist108, int num_events_108, int num_tracks_cut=10000); 

void draw_from_summary()
{
  bool ana_mult = false;
  bool ana_val = true;

  ejungwoo::gstat(0);

  TString fParticleNames[6] = {"p","d","t","he3","he4","he6"};
  double fParticleMass[6] = {938.272, 1871.06, 2809.41, 2809.41, 3728.4, 5606.55};
  int fParticleA[6] = {1,2,3,3,4,6};
  int fParticleZ[6] = {1,1,1,2,2,2};


  ejungwoo::binning bmult(100,0,100);
  ejungwoo::binning bmult2(12,0,12);
  ejungwoo::binning bkeoa(40,0,400);
  ejungwoo::binning bptoa(200,0,400);
  ejungwoo::binning bny(200,-1.5,3.0);
  ejungwoo::binning btta(200,0,180);
  ejungwoo::binning bphi(200,-180,180);
  
  ejungwoo::variable var_eff("eff","eff");
  ejungwoo::variable var_keoa_cm("keoa_cm","ke_cm/$$(a)","$$(gcut)","KE_{CM}/A (MeV)",bkeoa);
  ejungwoo::variable var_ptoa_cm("ptoa_cm","pt_cm/$$(a)","$$(gcut)","p_{T}/A (MeV/c)",bptoa);
  ejungwoo::variable var_ny_cm("ny_cm","ny_cm","$$(gcut)","y_{CM}/y_{beam,CM}",bny);
  ejungwoo::variable var_tta_cm("tta_cm","theta_cm*TMath::RadToDeg()","$$(gcut)","#theta_{CM} (Deg.)",btta);
  ejungwoo::variable var_phi_cm("phi_cm","phi_cm*TMath::RadToDeg()","$$(gcut)","#phi_{CM} (Deg.)",bphi);

  vector<ejungwoo::variable> varlist = {
    var_keoa_cm,
    //var_phi_cm+var_tta_cm,
    //var_ny_cm+var_ptoa_cm,
  };

  vector<ejungwoo::variable> varlist2 = {
    ejungwoo::variable("n_good:n","",bmult,bmult)
  };

  vector<ejungwoo::variable> varlist3 = {
    ejungwoo::variable("p",  "np_good",  "","mult",bmult2),
    ejungwoo::variable("d",  "nd_good",  "","mult",bmult2),
    ejungwoo::variable("t",  "nt_good",  "","mult",bmult2),
    ejungwoo::variable("he3","nhe3_good","","mult",bmult2),
    ejungwoo::variable("he4","nhe4_good","","mult",bmult2),
  };

  //ejungwoo::gsetcut(TCut("($$(tta_cm)<110)/$$(eff)"));
  ejungwoo::gsetcut(TCut("($$(tta_cm)<110)"));


  //for (auto version : {"NewAna.2057.254c9ba"})
  for (auto version : {"NewAna.2034.45b9400" ,"NewAna.2057.254c9ba"})
  //for (auto version : {"NewAna.2034.45b9400"})
  {
    TString vshort = TString("v")+ejungwoo::tok(version,".",1);
    ejungwoo::gversion(version);

    TH1D *histKE[6][2] = {0};
    int num_events[2] = {0};

    //for (auto sys : {108,132})
    for (auto sys : {132,108})
    {
      int isys = (sys==132?0:1);
      ejungwoo::setpar("sys",sys);

      {
        TString fileName = Form("data_summary/mult_all_%d.%s.root",sys,ejungwoo::versioncc());
        auto file = new TFile(fileName);
        auto tree = (TTree *) file -> Get("mult");
        num_events[isys] = tree -> GetEntries();
        cout << sys << " " << num_events[isys] << endl;

        if (ana_mult) {
          for (auto var : varlist2) 
            ejungwoo::drawv(var,tree);
          for (auto var : varlist3)
            ejungwoo::addto(TString("np")+sys+"_"+vshort,ejungwoo::tp(var,tree),"c");
        }
      }

      if (ana_val)
      {
        TString fileName = Form("data_summary/summary_all_%d.%s.root",sys,ejungwoo::versioncc());
        auto file = new TFile(fileName);
        for (int ipid : {0,1,2,3,4})
        {
          ejungwoo::gheader(fParticleNames[ipid]+"_"+sys+"_");
          ejungwoo::setpar("a",fParticleA[ipid]);
          auto tree = (TTree *) file -> Get(fParticleNames[ipid]);
          for (auto var : varlist)  {
            //auto ch = ejungwoo::drawv(var,tree);
            if (var.name=="keoa_cm")
              histKE[ipid][isys] = (TH1D *) ejungwoo::tp(var,tree);
          }
        }
      }
    }

    TString cvsname = TString("IS_")+vshort;
    for (int ipid : {0,1,2,3,4}) {
      auto graphis = draw_is(ipid, histKE[ipid][0], num_events[0], histKE[ipid][1], num_events[1]);
      ejungwoo::addto(cvsname,graphis, "pl", fParticleNames[ipid]);
    }
  }

  ejungwoo::drawall();
}

TGraphErrors *draw_is(int idx_particle, TH1 *hist132, int num_events_132, TH1 *hist108, int num_events_108, int num_tracks_cut) {
  auto graph_is = new TGraphErrors();
  ejungwoo::binning binn(hist132);

  for (auto bin=1; bin<=binn.n; ++bin)
  {
    auto v132 = hist132 -> GetBinContent(bin);
    auto v108 = hist108 -> GetBinContent(bin);

    auto num_tracks_132 = num_events_132*v132;
    auto num_tracks_108 = num_events_108*v108;

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
