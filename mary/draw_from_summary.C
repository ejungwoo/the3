#include "init_variables.h"

using ejungwoo::variable;
using ejungwoo::binning;
using ejungwoo::titles;
using ejungwoo::setpar;

TObjArray draw_is(int idx_particle, TH1 *hist132, TH1 *hist108, double num_tracks_cut=0); 

void draw_from_summary()
{
  ejungwoo::gcvspos(1300);
  ejungwoo::gstat(0);
  ejungwoo::gsave(0);
  //ejungwoo::gdummytp();
  ejungwoo::gfast();

  bool ana_is   = 0;
  bool ana_mult = 0;
  bool ana_val  = 1;
  bool ana_pid  = 0;
  bool ana_all  = 0;
  double num_tracks_per_event_cut = 0.02;

  //setpar("angle_cut","($$(phi_cm)<-100||$$(phi_cm)>100)");
  //setpar("angle_cut","($$(tta_cm)<110)");
  setpar("angle_cut","1");

  //setpar("poz_cut","($$(p_lab)>$$(pozll))");
  setpar("poz_cut","1");

  //setpar("pes_cut","($$(asd)<3)");
  setpar("pes_cut","1");

  vector<TString> spversions = { "NewAna.2070.0231288", };
  int ipids[] = {0,1,2,3,4};
  int isyss[] = {0,1,2,3};
  int systems[] = {132   , 108   , 112   , 124   };
  double yAAs[] = {0.3822, 0.3647, 0.3538, 0.3902};
  double yNNs[] = {0.3696, 0.3697, 0.3705, 0.3706};
  Long64_t num_events[4] = {0,};

  TString particleTreeNames[6] = {"p","d","t","he3","he4","he6"};
  TString particleNames[6] = {"p","d","t","he3","he4","he6"};
  double particleMass[6] = {938.272, 1871.06, 2809.41, 2809.41, 3728.4, 5606.55};
  double particlePozCut[6] = {100,300,300,100,200,200};
  int particlePDGs[6] = {2212, 1000010030, 1000010020, 1000020040, 1000020030,};
  int particleA[6] = {1,2,3,3,4,6};
  int particleZ[6] = {1,1,1,2,2,2};
  int numProtons[6]  = {1,1,1,2,2,2};
  int numNeutrons[6] = {0,1,2,1,2,4};

  auto setpar_syspid = [
    &systems, &yAAs, &yNNs, &num_events, 
    &particleNames, &particleMass, &particlePozCut, &particlePDGs, 
    &particleA, &particleZ, &numProtons, &numNeutrons
  ] (int isys, int ipid)
  {
    setpar("sys",systems[isys]);
    setpar("yaa",yAAs[isys]);
    setpar("ynn",yNNs[isys]);
    setpar("num_events",num_events[isys]);
    setpar("pname",particleNames[ipid]);
    setpar("pozll",particlePozCut[ipid]);
    setpar("pdg",particlePDGs[ipid]);
    setpar("m",particleMass[ipid]);
    setpar("a",particleA[ipid]);
    setpar("z",particleZ[ipid]);
    setpar("nump",numProtons[ipid]);
    setpar("numn",numNeutrons[ipid]);
  };

  ejungwoo::cutg("/Users/ejungwoo/spirit/the3/mary/data__NewAna.2070.0231288/p_ypt_test_cutg.NewAna.2070.0231288.root","p_ypt_test_cutg","$$(ny_cm)","pt_cm");
  ejungwoo::cutg("/Users/ejungwoo/spirit/the3/mary/data__NewAna.2070.0231288/d_ypt_test_cutg.NewAna.2070.0231288.root","d_ypt_test_cutg","$$(ny_cm)","pt_cm/2");
  ejungwoo::cutg("/Users/ejungwoo/spirit/the3/mary/data__NewAna.2070.0231288/he3_ypt_test_cutg.NewAna.2070.0231288.root","he3_ypt_test_cutg","$$(ny_cm)","pt_cm/3");
  ejungwoo::cutg("/Users/ejungwoo/spirit/the3/mary/data__NewAna.2070.0231288/ptest1.NewAna.2070.0231288.root","ptest1","$$(ny_cm)","pt_cm/1");
  ejungwoo::cutg("/Users/ejungwoo/spirit/the3/mary/data__NewAna.2070.0231288/ttest1.NewAna.2070.0231288.root","ttest1","$$(ny_cm)","pt_cm/3");

  setpar("cut_n",                "1./$$(num_events)*$$(pes_cut)*$$(poz_cut)");
  setpar("cut_p",          "$$(prob)/$$(num_events)*$$(pes_cut)*$$(poz_cut)*$$(angle_cut)");
  setpar("cut_e",        "1./$$(eff)/$$(num_events)*$$(pes_cut)*$$(poz_cut)*$$(angle_cut)");
  setpar("cut_pe", "$$(prob)/$$(eff)/$$(num_events)*$$(pes_cut)*$$(poz_cut)*$$(angle_cut)");

  setpar("cut_prob6", "$$(prob)>.6");
  setpar("cut_ptest1", "(ptest1)*$$(cut_pe)");
  setpar("cut_ttest1", "(ttest1)*$$(cut_pe)");
  //setpar("cut0", "$$(cut_ptest1)");
  setpar("cut_compare", "$$(cut_ptest1)");
  //setpar("cut0", "$$(cut_compare)");
  setpar("cut0", "$$(cut_p)");

  vector<variable> varListAll = {
    fvar_prob, fvar_eff, fvar_sd,
    fvar_pt_cm, fvar_fy_cm, fvar_p_cm, fvar_ke_cm, fvar_p_lab, fvar_dedx,
    fvar_p_lab+fvar_dedx,
    //fvar_dpoca, fvar_nr, fvar_nl,
    //fvar_tta_cm, fvar_phi_cm, fvar_tta_lab, fvar_phi_lab,
  };
  vector<variable> varListMain = { fvar_phi_lab+fvar_tta_lab, fvar_phi_cm+fvar_tta_cm, fvar_ny_cm+fvar_ptoa_cm, fvar_ny_cm+fvar_eff, fvar_tta_cm+fvar_eff};
  vector<variable> varListIS   = { fvar_ny_cm, fvar_keoa_cm, fvar_ptoa_cm };
  vector<variable> varListN1   = { fvar_np, fvar_nd, fvar_nt, fvar_nhe3, fvar_nhe4, };
  vector<variable> varListN2   = { fvar_ngp, fvar_ngd, fvar_ngt, fvar_nghe3, fvar_nghe4, };
  vector<variable> varListN    = { fvar_n+fvar_ng, fvar_np+fvar_ngp, fvar_nd+fvar_ngd, fvar_nt+fvar_ngt, fvar_nhe3+fvar_nghe3, fvar_nhe4+fvar_nghe4, };
  /*
  for (auto varList : {varListAll, varListIS, varListMain, varListN, varListN1, varListN2}) {
    for (auto var : varList)
      var.cut = ejungwoo::getpar("cut_pe");
  }
  */
  auto fvar_pid = fvar_p_lab+fvar_dedx+"z";
  fvar_pid.cut = ejungwoo::getpar("cut_prob6");

  /******************************************************************************************/

  auto filem = new TFile("particle_cutg_s20_p2.pidcut.root");
  TCutG *cutgs[6];
  TF1 *means[6];
  for (auto ipid : ipids) {
    cutgs[ipid] = (TCutG *) filem -> Get(Form("CUTG%d",particlePDGs[ipid]));
    cutgs[ipid] -> SetLineColor(kRed);
    means[ipid] = (TF1 *) filem -> Get(Form("BEE%d",particlePDGs[ipid]));
    means[ipid] -> SetLineColor(kGray+2);
    means[ipid] -> SetLineStyle(1);
    means[ipid] -> Print();
  }

  /******************************************************************************************/

  for (auto version : spversions)
  {
    ejungwoo::gversion(version);
    TString vshort = TString("v")+ejungwoo::tok(version,".",1);
    setpar("vshort", vshort);

    TTree *tree_pid[4][5] = {0};
    TTree *tree_mult[4] = {0};
    for (auto isys : {0,1,2,3}) {
      auto sys = systems[isys];
      TString fileName = Form("data__%s/sys%d.%s.ana.particle.root",ejungwoo::versioncc(),sys,ejungwoo::versioncc());
      auto file = new TFile(fileName);
      tree_mult[isys] = (TTree *) file -> Get("mult");
      for (int ipid : ipids) {
        tree_pid[isys][ipid] = (TTree *) file -> Get(particleTreeNames[ipid]);
      }
      num_events[isys] = tree_mult[isys] -> GetEntries();
      cout << fileName << " " << sys << " " << num_events[isys] << endl;
    }

    for (auto isys : isyss)
    {
      auto sys = systems[isys];

      TTree *tree = tree_mult[isys];
      if (ana_mult) {
        for (auto var : varListN) {
          var.drawaddnext(tree,ejungwoo::fHeader+"nn"+sys);
        }
        for (auto var : varListN1) var.drawadd(tree,TString("no")+sys,0,"hist");
        for (auto var : varListN2) var.drawadd(tree,TString("ng")+sys,0,"hist");
      }

      if (ana_all) {
        for (int ipid : ipids) {
          setpar_syspid(isys,ipid);
          auto tree = tree_pid[isys][ipid];
          for (auto var : varListAll) {
            TString ename = TString("all__")+particleNames[ipid]+sys;
            var.cut = "";
            auto hist1 = ejungwoo::norm_max(ejungwoo::tp(var,tree));
            ejungwoo::addnext(ename, hist1, "hist");
            var.name = var.name+2;
            var.cut = "$$(cut_compare)";
            auto hist2 = ejungwoo::norm_max(ejungwoo::tp(var,tree));
            ejungwoo::addsame(ename, hist2, "hist");
          }
        }
      }

      for (int ipid : ipids)
      {
        auto tree = tree_pid[isys][ipid];
        setpar_syspid(isys,ipid);

        if (ana_val)
          for (auto var : varListMain)
            var.drawadd(tree,var.name+sys,ipid);

        if (ana_pid) {
          fvar_pid.cut = ejungwoo::getpar("cut_prob6");
          auto hist_pid = ejungwoo::tp(fvar_pid,tree);
          ejungwoo::set_vor0(hist_pid,ipid+1);
          ejungwoo::addhist(TString("pid")+sys, hist_pid, "");
          ejungwoo::add(TString("pid")+sys, cutgs[ipid],"addx colorx samel");
          ejungwoo::add(TString("pid")+sys, means[ipid],"addx colorx samel");

          fvar_pid.cut = ejungwoo::getpar("cut_pe");
          hist_pid = ejungwoo::tp(fvar_pid,tree);
          ejungwoo::addhist(TString("pid_pe_")+sys, hist_pid, "logz");
          ejungwoo::add(TString("pid_pe_")+sys, cutgs[ipid],"addx colorx samel");
          ejungwoo::add(TString("pid_pe_")+sys, means[ipid],"addx colorx samel");

          fvar_pid.cut = ejungwoo::getpar("cut_compare");
          hist_pid = ejungwoo::tp(fvar_pid,tree);
          ejungwoo::addhist(TString("pid_compare_")+sys, hist_pid, "logz");
          ejungwoo::add(TString("pid_compare_")+sys, cutgs[ipid],"addx colorx samel");
          ejungwoo::add(TString("pid_compare_")+sys, means[ipid],"addx colorx samel");

          fvar_pid.cut = ejungwoo::getpar("cut_p");
          hist_pid = ejungwoo::tp(fvar_pid,tree);
          ejungwoo::addhist(TString("pid_p_")+sys, hist_pid, "logz");
          ejungwoo::add(TString("pid_p_")+sys, cutgs[ipid],"addx colorx samel");
          ejungwoo::add(TString("pid_p_")+sys, means[ipid],"addx colorx samel");
        }
      }
    }

    if (ana_is)
    {
      ejungwoo::gheader("");
      int ivar = 0;
      for (auto var : varListIS) {
        setpar("pname","all");
        TString cname_nn = TString("NN_")+var.name+"__"+vshort;
        TString cname_is = TString("IS_")+var.name+"__"+vshort;

        var.binn.n=50;
        var.title.y = "p(132+124)/p(108+112)";
        ejungwoo::add(cname_is,new_h(var.hist_name,var.title,var.binn,binning(100,0,3)),"addx");
        for (int ipid : ipids) {
          setpar_syspid(0,ipid); auto hist1 = ejungwoo::tp(var,tree_pid[0][ipid]); ejungwoo::add(cname_nn,ipid,hist1,"hist","132");
          setpar_syspid(1,ipid); auto hist2 = ejungwoo::tp(var,tree_pid[1][ipid]); ejungwoo::add(cname_nn,ipid,hist2,"hist","108");
          auto graphs = draw_is(ipid, hist1, hist2, num_tracks_per_event_cut);
          auto graphis = (TGraphErrors *) graphs.At(0);
          ejungwoo::add(cname_is, graphis, "pl", particleNames[ipid]);
        }
        ivar++;
      }
    }
  }

  /******************************************************************************************/

  ejungwoo::gheader("");

  ejungwoo::drawsaveall("cvsl","png");
}

TObjArray draw_is(int idx_particle, TH1 *hist132, TH1 *hist108, double num_tracks_cut)
{
  auto graph_is = new TGraphErrors();
  auto graph_v132 = new TGraphErrors();
  auto graph_v108 = new TGraphErrors();
  binning binn(hist132);

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
    graph -> SetLineColor(ejungwoo::colori(idx_particle));
    graph -> SetLineWidth(2);
    array.Add(graph);
  }

  return array;
}
