#include "init_variables.h"

using ejungwoo::variable;
using ejungwoo::binning;
using ejungwoo::titles;
using ejungwoo::setpar;

const int iversions[] = {-1};
const int ipids[] = {0,1,2,3,4};
const int isyss[] = {0};
const int numProtons[6]  = {1,1,1,2,2,2};
const int numNeutrons[6] = {0,1,2,1,2,4};

TString spversion;
TString vshort;
TString vout;
double phiLL;
double phiHL;
double ttaLL;
double ttaHL;
double dsdinpi;
double dsd;
int trackMultHL;
double probLL;
double effLL;
double ptoaLL;
double ptoaHL;

TObjArray draw_is(int idx_particle, TH1 *hist132, TH1 *hist108, double num_tracks_cut=0); 
TObjArray draw_ci(vector<TH1 *>hists, double ds, int isys);
TGraphErrors *draw_ratio(TGraphErrors *graph1, TGraphErrors *graph2, double max);

void draw_from_summary()
{
  ejungwoo::gcvspos(700);
  ejungwoo::gstat(0);
  ejungwoo::gsave(0);
  //ejungwoo::gdummytp();
  ejungwoo::gfast();
  ejungwoo::gsetvmark(0);

  bool withTestCuts = 0;
  int overlabTestCuts = 2;

  bool anaIS = 1;
  bool anaCI = 1;
  bool anaMain = 1;
  bool anaAngle = 1;
  bool anaSingle = 1;
  bool anaMult = 1;
  bool anaPid = 1;
  bool anaAll = 1;
  bool anaAlltta = false;

  //double num_tracks_per_event_cut = 0.02;
  double num_tracks_per_event_cut = 0.01;
  int nbinsIS = 20;
  double dKECI = 5.;

  auto setversion = [](int iversion) {
    if (iversion==-1) {
      spversion = "NewAna.2070.0231288"; 
      vshort = TString((ejungwoo::tok(spversion,".",0))[0])+ejungwoo::tok(spversion,".",1);
      vout = vshort+".raw"; 
      trackMultHL = 50;
      probLL = 0.1;
      effLL = 0.01;
      ptoaLL = 0;
      ptoaHL = 10000;
      phiLL = 160;
      phiHL = 200;
      ttaLL = 0;
      ttaHL = 180;
      fvar_ny_cm.binn = binning(100,-.5,.5);
      dsdinpi = 2 * (phiHL-phiLL)/360 * (new TF1("satta","sin(x)",0,TMath::Pi())) -> Integral(ttaLL*TMath::DegToRad(),ttaHL*TMath::DegToRad());
      dsd = dsdinpi * TMath::Pi();
      setpar("cut0", "$$(cut_pe)");
      setpar("cutn", "");
    }
    if (iversion==0) {
      spversion = "NewAna.2070.0231288"; 
      vshort = TString((ejungwoo::tok(spversion,".",0))[0])+ejungwoo::tok(spversion,".",1);
      vout = vshort+".tta711"; 
      trackMultHL = 50;
      probLL = 0.2;
      effLL = 0.08;
      ptoaLL = 0;
      ptoaHL = 10000;
      phiLL = 160;
      phiHL = 200;
      ttaLL = 70;
      ttaHL = 110;
      fvar_ny_cm.binn = binning(100,-.5,.5);
      dsdinpi = 2 * (phiHL-phiLL)/360 * (new TF1("satta","sin(x)",0,TMath::Pi())) -> Integral(ttaLL*TMath::DegToRad(),ttaHL*TMath::DegToRad());
      dsd = dsdinpi * TMath::Pi();
      setpar("cut0", "$$(cut_pe)");
      setpar("cutn", "");
    }
    else if (iversion==1) {
      spversion = "NewAna.2070.0231288"; 
      vshort = TString((ejungwoo::tok(spversion,".",0))[0])+ejungwoo::tok(spversion,".",1);
      vout = vshort+".nocut"; 
      trackMultHL = 50;
      probLL = 0.2;
      effLL = 0.08;
      ptoaLL = 0;
      ptoaHL = 10000;
      phiLL = 160;
      phiHL = 200;
      ttaLL = 0;
      ttaHL = 180;
      dsdinpi = 2 * (phiHL-phiLL)/360 * (new TF1("satta","sin(x)",0,TMath::Pi())) -> Integral(ttaLL*TMath::DegToRad(),ttaHL*TMath::DegToRad());
      dsd = dsdinpi * TMath::Pi();
      setpar("cut0", "$$(cut_pe)");
      setpar("cutn", "");
    }
    else if (iversion==2) {
      spversion = "NewAna.2070.0231288"; 
      vshort = TString((ejungwoo::tok(spversion,".",0))[0])+ejungwoo::tok(spversion,".",1);
      vout = vshort+".corrp_xe"; 
      trackMultHL = 50;
      probLL = 0.2;
      effLL = 0.05;
      ptoaLL = 0;
      ptoaHL = 10000;
      phiLL = 160;
      phiHL = 200;
      ttaLL = 0;
      ttaHL = 180;
      dsdinpi = 2 * (phiHL-phiLL)/360 * (new TF1("satta","sin(x)",0,TMath::Pi())) -> Integral(ttaLL*TMath::DegToRad(),ttaHL*TMath::DegToRad());
      dsd = dsdinpi * TMath::Pi();
      setpar("cut0", "$$(cut_p)");
      setpar("cutn", "");
    }

    vector<TString> entry_array;
    entry_array.push_back(Form("version=%s",spversion.Data()));
    entry_array.push_back(Form("N_{track}>%d",trackMultHL));
    entry_array.push_back(Form("probability>%.2f",probLL));
    entry_array.push_back(Form("efficiency>%.2f",effLL));
    if (ptoaLL<=0&&ptoaHL>=5000) entry_array.push_back("p_{T} cut: X");
    else                         entry_array.push_back(Form("%.2f<#p{T}/A_{CM}<%.2f",ptoaLL,ptoaHL));
    entry_array.push_back(Form("%.f<#phi_{CM}<%.f",phiLL,phiHL));
    entry_array.push_back(Form("%.f<#theta_{CM}<%.f",ttaLL,ttaHL));
    entry_array.push_back(Form("#Delta#Omega=%.3f#pi",dsdinpi));

    setpar("angle_cut",Form("($$(tta_cm)>%f&&$$(tta_cm)<%f)",ttaLL,ttaHL));
    setpar("poz_cut", "1");
    setpar("pes_cut",Form("$$(eff)>%f",effLL));

    ejungwoo::gversionin(spversion);
    ejungwoo::gversionout(vout);

    for (auto entry : entry_array) {
      cout << entry << endl;
    }
  };

  int systems[] = {132   , 108   , 112   , 124   };
  int targetA[] = {124   , 112   , 124   , 112   };
  double yAAs[] = {0.3822, 0.3647, 0.3538, 0.3902};
  double yNNs[] = {0.3696, 0.3697, 0.3705, 0.3706};
  Long64_t num_events[4] = {0};
  TString particleTreeNames[6] = {"p","d","t","he3","he4","he6"};
  TString particleNames[6] = {"p","d","t","he3","he4","he6"};
  double particleMass[6] = {938.272, 1871.06, 2809.41, 2809.41, 3728.4, 5606.55};
  double particlePozCut[6] = {100,300,300,100,200,200};
  int particlePDGs[6] = {2212, 1000010030, 1000010020, 1000020040, 1000020030,};
  int particleA[6] = {1,2,3,3,4,6};
  int particleZ[6] = {1,1,1,2,2,2};

  int numpids = 0;
  int numsyss = 0;
  for (auto i : ipids) numpids++;
  for (auto i : isyss) numsyss++;

  auto setpar_syspid = [
    &systems, &yAAs, &yNNs, &num_events, 
    &particleNames, &particleMass, &particlePozCut, &particlePDGs, 
    &particleA, &particleZ
  ] (int isys, int ipid=0)
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

  ejungwoo::cutg("/Users/ejungwoo/spirit/the3/mary/data__NewAna.2070.0231288/pefftest1.NewAna.2070.0231288.root","pefftest1","$$(tta_cm)","eff");
  ejungwoo::cutg("/Users/ejungwoo/spirit/the3/mary/data__NewAna.2070.0231288/pefftest2.NewAna.2070.0231288.root","pefftest2","$$(tta_cm)","eff");
  ejungwoo::cutg("/Users/ejungwoo/spirit/the3/mary/data__NewAna.2070.0231288/pefftest3.NewAna.2070.0231288.root","pefftest3","$$(tta_cm)","eff");
  ejungwoo::cutg("/Users/ejungwoo/spirit/the3/mary/data__NewAna.2070.0231288/pefftest4.NewAna.2070.0231288.root","pefftest4","$$(tta_cm)","eff");
  ejungwoo::cutg("/Users/ejungwoo/spirit/the3/mary/data__NewAna.2070.0231288/pefftest5.NewAna.2070.0231288.root","pefftest5","$$(tta_cm)","eff");
  //vector<TString> testCutArray = {"pefftest1", "pefftest2", "pefftest3", "pefftest4", "pefftest5",};

  ejungwoo::cutg("/Users/ejungwoo/spirit/the3/mary/data__NewAna.2070.0231288/pefftest_f1.NewAna.2070.0231288.root","pefftest_f1","$$(tta_cm)","eff");
  ejungwoo::cutg("/Users/ejungwoo/spirit/the3/mary/data__NewAna.2070.0231288/pefftest_f2.NewAna.2070.0231288.root","pefftest_f2","$$(tta_cm)","eff");
  //vector<TString> testCutArray = {"pefftest_f1", "pefftest_f2"};

  ejungwoo::cutg("/Users/ejungwoo/spirit/the3/mary/data__NewAna.2070.0231288/tpytest1.NewAna.2070.0231288.root","tpytest1","$$(ny_cm)","pt_cm/3");
  ejungwoo::cutg("/Users/ejungwoo/spirit/the3/mary/data__NewAna.2070.0231288/tpytest2.NewAna.2070.0231288.root","tpytest2","$$(ny_cm)","pt_cm/3");
  vector<TString> testCutArray = {"tpytest1", "tpytest2"};

  testCutArray.push_back("1");
 
  setpar("cut_n",                "1./$$(num_events)*($$(pes_cut))*($$(poz_cut))");
  setpar("cut_p",          "$$(prob)/$$(num_events)*($$(pes_cut))*($$(poz_cut))*($$(angle_cut))");
  setpar("cut_e",        "1./$$(eff)/$$(num_events)*($$(pes_cut))*($$(poz_cut))*($$(angle_cut))");
  setpar("cut_pe", "$$(prob)/$$(eff)/$$(num_events)*($$(pes_cut))*($$(poz_cut))*($$(angle_cut))");

  setpar("cut_prob6", "$$(prob)>.6");

  setpar("cut_ptest1", "(ptest1)*$$(cut_pe)");
  setpar("cut_ttest1", "(ttest1)*$$(cut_pe)");

  //setpar("cut0", "$$(cut_ptest1)");
  setpar("cut_compare", "$$(cut_p)");

  //setpar("cut0", "$$(cut_compare)");

  vector<variable> varListAll = {
    fvar_prob, fvar_eff, fvar_sd,
    fvar_ptoa_cm, fvar_pt_cm, fvar_fy_cm, fvar_p_cm, fvar_ke_cm,
    fvar_p_lab, fvar_dedx,
    fvar_dpoca, fvar_nr, fvar_nl,
    fvar_tta_cm, fvar_phi_cm,
    //fvar_tta_lab, fvar_phi_lab,
    //fvar_p_kab+fvar_dedx,
  };

  vector<variable> varListSingle = { fvar_keoa_cm+fvar_ny_cm };
  vector<variable> varListAngle  = { fvar_phi_lab+fvar_tta_lab, fvar_phi_cm+fvar_tta_cm, fvar_phi_lab+fvar_phi_cm, fvar_tta_lab+fvar_tta_cm, };
  vector<variable> varListMain   = { fvar_ny_cm+fvar_ptoa_cm, fvar_ny_cm+fvar_eff, fvar_tta_cm+fvar_eff, };
  vector<variable> varListIS     = { fvar_ny_cm };
  vector<variable> varListCI     = { fvar_keoa_cm };
  vector<variable> varListN1     = { fvar_np, fvar_nd, fvar_nt, fvar_nhe3, fvar_nhe4, };
  vector<variable> varListN2     = { fvar_ngp, fvar_ngd, fvar_ngt, fvar_nghe3, fvar_nghe4, };
  vector<variable> varListN      = { fvar_n+fvar_ng, fvar_np+fvar_ngp, fvar_nd+fvar_ngd, fvar_nt+fvar_ngt, fvar_nhe3+fvar_nghe3, fvar_nhe4+fvar_nghe4, };

  vector<variable> varListAlltta;
  for (auto var : varListAll)
    if (var.name==fvar_tta_cm.name)
      varListAlltta.push_back(fvar_tta_cm);
    else
      varListAlltta.push_back(fvar_tta_cm+var);

  auto fvar_pid = fvar_p_lab+fvar_dedx+"z";
  fvar_pid.cut = ejungwoo::getpar("cut_prob6");

  /******************************************************************************************/

  auto filem = new TFile("data/particle_cutg_s20_p2.pidcut.root");
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

  for (auto iversion : iversions)
  {
    setversion(iversion);
    setpar("vshort", vshort);

    TTree *tree_pid[4][5] = {0};
    TTree *tree_mult[4] = {0};
    for (auto isys : {0,1,2,3}) {
      auto sys = systems[isys];
      //TString fileName = Form("data__%s/sys%d.%s.ana.particle.root",ejungwoo::versioncc(),sys,ejungwoo::versioncc());
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
      if (anaMult) {
        setpar_syspid(isys);
        for (auto var : varListN)  { var.print(); var.drawaddnext(tree,ejungwoo::fHeader+"nn"+sys); }
        for (auto var : varListN1) { var.print(); var.drawadd(tree,TString("no")+sys,0,"hist"); }
        for (auto var : varListN2) { var.print(); var.drawadd(tree,TString("ng")+sys,0,"hist"); }
      }

      if (anaAll) {
        for (int ipid : ipids) {
          setpar_syspid(isys,ipid);
          auto tree = tree_pid[isys][ipid];
          for (auto var : varListAll) {
            if (withTestCuts) {
              auto cut0 = var.cut;
              TString cutName0 = var.cut.GetName();
              var.name = var.name + "_all";
              auto itest = 0;
              for (auto testCut : testCutArray) {
                var.cut = cut0*TCut(testCut.Data());
                var.cut.SetName(cutName0+"_"+testCut);
                var.drawaddsame(tree,TString("all__")+var.name+ipid,"hist",testCut);
                itest++;
              }
            }
            else{
              TString ename = TString("all__")+particleNames[ipid]+sys;
              var.cut = "$$(cut_pe)";
              auto hist1 = ejungwoo::norm_max(var.draw(tree));
              ejungwoo::addnext(ename, hist1, "hist addx", "pe");
              //var.name = var.name+2;
              //var.cut = "$$(cut_p)";
              //auto hist2 = ejungwoo::norm_max(var.draw(tree));
              //ejungwoo::addsame(ename, hist2, "hist", "p");
            }
          }

        }
      }

      if (anaAlltta) {
        for (int ipid : ipids) {
          setpar_syspid(isys,ipid);
          auto tree = tree_pid[isys][ipid];
          for (auto var : varListAlltta) {
            TString ename = TString("alltta__")+particleNames[ipid]+sys;
            ejungwoo::addnext(ename, var.draw(tree), "colz");
          }
        }
      }

      for (int ipid : ipids)
      {
        auto tree = tree_pid[isys][ipid];
        setpar_syspid(isys,ipid);

        if (anaAngle)
          for (auto var : varListAngle) {
            if (withTestCuts) {
              auto cut0 = var.cut;
              TString cutName0 = var.cut.GetName();
              TString varName0 = var.name;
              auto itest = 0;
              for (auto testCut : testCutArray) {
                var.cut = cut0*TCut(testCut.Data());
                var.cut.SetName(cutName0+"_"+testCut);
                var.name = TString("witCut_") + varName0;
                auto histCut = var.draw(tree);
                if (overlabTestCuts>0) {
                  var.cut = cut0;
                  var.cut.SetName(cutName0+"_for_"+testCut);
                  var.name = TString("witCut_") + varName0;
                  auto hist0 = var.draw(tree);
                  if (overlabTestCuts==1)
                    ejungwoo::set_vor0(hist0,histCut->GetMaximum()*0.001);
                  histCut -> Add(hist0);
                  //ejungwoo::add(TString("addtoCut_")+varName0+"_"+ipid,itest,hist0);
                }
                ejungwoo::add(varName0+"_"+ipid,itest,histCut);
                itest++;
              }
            }
            else
              var.drawadd(tree,var.name+sys,ipid);
          }
        if (anaMain)
          for (auto var : varListMain) {
            if (withTestCuts) {
              auto cut0 = var.cut;
              TString cutName0 = var.cut.GetName();
              TString varName0 = var.name;
              auto itest = 0;
              for (auto testCut : testCutArray) {
                var.cut = cut0*TCut(testCut.Data());
                var.cut.SetName(cutName0+"_"+testCut);
                var.name = TString("witCut_") + varName0;
                auto histCut = var.draw(tree);
                if (overlabTestCuts>0) {
                  var.cut = cut0;
                  var.cut.SetName(cutName0+"_for_"+testCut);
                  var.name = TString("witCut_") + varName0;
                  auto hist0 = var.draw(tree);
                  if (overlabTestCuts==1)
                    ejungwoo::set_vor0(hist0,histCut->GetMaximum()*0.001);
                  histCut -> Add(hist0);
                  //ejungwoo::add(TString("addtoCut_")+varName0+"_"+ipid,itest,hist0);
                }
                ejungwoo::add(varName0+"_"+ipid,itest,histCut);
                itest++;
              }
            }
            else
              var.drawadd(tree,var.name+sys,ipid);
          }

        if (anaSingle)
          for (auto var : varListSingle) {
            if (withTestCuts) {
              auto cut0 = var.cut;
              TString cutName0 = var.cut.GetName();
              //TString varName0 = var.name;
              var.name = var.name + "_single";
              auto itest = 0;
              for (auto testCut : testCutArray) {
                var.cut = cut0*TCut(testCut.Data());
                var.cut.SetName(cutName0+"_"+testCut);
                var.drawadd(tree,TString("single_")+var.name+ipid,itest);
                itest++;
              }
            }
            else {
              var.drawadd(tree,TString("single_")+var.name + ipid);
            }
          }

        if (anaPid) {
          auto var = fvar_pid;

          var.cut = ejungwoo::getpar("cut_prob6");
          auto hist_pid = var.draw(tree);
          ejungwoo::set_vor0(hist_pid,ipid+1);
          ejungwoo::addhist(TString("pid")+sys, hist_pid, "");
          ejungwoo::add(TString("pid")+sys, cutgs[ipid],"addx colorx samel");
          ejungwoo::add(TString("pid")+sys, means[ipid],"addx colorx samel");

          var.cut = ejungwoo::getpar("cut_pe");
          hist_pid = var.draw(tree);
          ejungwoo::addhist(TString("pid_pe_")+sys, hist_pid, "logz");
          ejungwoo::add(TString("pid_pe_")+sys, cutgs[ipid],"addx colorx samel");
          ejungwoo::add(TString("pid_pe_")+sys, means[ipid],"addx colorx samel");

          var.cut = ejungwoo::getpar("cut_compare");
          hist_pid = var.draw(tree);
          ejungwoo::addhist(TString("pid_compare_")+sys, hist_pid, "logz");
          ejungwoo::add(TString("pid_compare_")+sys, cutgs[ipid],"addx colorx samel");
          ejungwoo::add(TString("pid_compare_")+sys, means[ipid],"addx colorx samel");

          var.cut = ejungwoo::getpar("cut_p");
          hist_pid = var.draw(tree);
          ejungwoo::addhist(TString("pid_p_")+sys, hist_pid, "logz");
          ejungwoo::add(TString("pid_p_")+sys, cutgs[ipid],"addx colorx samel");
          ejungwoo::add(TString("pid_p_")+sys, means[ipid],"addx colorx samel");
        }
      }
    }

    if (anaIS)
    {
      for (auto var : varListIS) {
        setpar("pname","all");
        TString cname_nn = TString("NN_")+var.name+"__"+vshort;
        TString cname_is = TString("IS_")+var.name+"__"+vshort;

        var.binn.n=nbinsIS;
        var.title.main = "$$(pname) multiplicity 132 vs 108";
        ejungwoo::add(cname_is,0,new_h(var.histname,titles("isoscaling plot 132/108",var.title.x,"p(132+124)/p(108+112)"),var.binn,binning(100,0,2)),"addx");
        double ymax1 = 0, ymax2 = 0;
        for (int ipid : ipids) {
          setpar_syspid(0,ipid); auto hist1 = var.draw(tree_pid[0][ipid]);
          setpar_syspid(1,ipid); auto hist2 = var.draw(tree_pid[1][ipid]);
          ejungwoo::add(cname_nn,ipid,hist1,"hist","132");
          ejungwoo::add(cname_nn,ipid,hist2,"hist","108");
          auto graphs = draw_is(ipid, hist1, hist2, num_tracks_per_event_cut);
          ejungwoo::add(cname_is, 0, graphs.At(0), "pl", particleNames[ipid]);
          ejungwoo::add(cname_is, 1, graphs.At(1), "pl", particleNames[ipid]);
          ejungwoo::add(cname_is, 2, graphs.At(2), "pl", particleNames[ipid]);
          double ymax1_ = ejungwoo::y2_g((TGraph *) graphs.At(1)); if (ymax1 < ymax1_) ymax1 = ymax1_;
          double ymax2_ = ejungwoo::y2_g((TGraph *) graphs.At(2)); if (ymax2 < ymax2_) ymax2 = ymax2_;
        }

        auto line1 = new TLine(var.binn.min,num_tracks_per_event_cut,var.binn.max,num_tracks_per_event_cut); line1 -> SetLineStyle(2);
        auto line2 = new TLine(var.binn.min,num_tracks_per_event_cut,var.binn.max,num_tracks_per_event_cut); line2 -> SetLineStyle(2);
        ejungwoo::add(cname_is,1,line1);
        ejungwoo::add(cname_is,2,line1);

        titles ttlIS132("132+124",var.title.x,fttly_ntne);
        titles ttlIS108("108+112",var.title.x,fttly_ntne);

        if (0) {
          ejungwoo::add(cname_is,1,new_h(var.histname+"_frame132",ttlIS132,var.binn,binning(100,ymax1*0.0001,ymax1*5)),"addx rangex logy");
          ejungwoo::add(cname_is,2,new_h(var.histname+"_frame108",ttlIS108,var.binn,binning(100,ymax2*0.0001,ymax2*5)),"addx rangex logy");
        } else {
          ejungwoo::add(cname_is,1,new_h(var.histname+"frame132",ttlIS132,var.binn,binning(100,0,ymax1*1.05)),"addx rangex");
          ejungwoo::add(cname_is,2,new_h(var.histname+"frame108",ttlIS108,var.binn,binning(100,0,ymax2*1.05)),"addx rangex");
        }
      }
    }

    if (anaCI)
    {
      titles ttlCIM("","KE_{CM} (MeV)","#frac{dM_{n,CI}}{dKE_{CM} #Delta#Omega_{CM}}");
      titles ttlCIR("","KE_{CM} (MeV)","R = #frac{dM_{n,CI}}{dKE_{CM} #Delta#Omega_{CM}} / #frac{dM_{p,CI}}{dKE_{CM} #Delta#Omega_{CM}}");
      binning binnCIx(40,0,200);
      binning binnCIy(100,0.00005,1);
      binning binnCIy2(100,0.001,1);
      binning binnCIr(100,0,1.2);
      binning binnCIdr(100,0,2);

      for (auto var : varListCI) {
        var.binn.setw(dKECI);
        var.name = TString("ci_")+var.name; 
        var.title.main = "";
        TString cname_ci = var.name+"__"+vshort;
        ejungwoo::add(cname_ci, 0, new_h(var.histname+"_cin_allframe", ttlCIM, binnCIx, binnCIy2), "addx rangex logy ");
        ejungwoo::add(cname_ci, 1, new_h(var.histname+"_cip_allframe", ttlCIM, binnCIx, binnCIy2), "addx rangex logy ");
        ejungwoo::add(cname_ci, 2, new_h(var.histname+"_cir_allframe", ttlCIR, binnCIx, binnCIr), "addx rangex ");

        TGraphErrors *graphs_rnp[4];
        for (auto isys : {0,1,2,3}) {
          auto sys = systems[isys];
          auto tga = targetA[isys];
          var.title.main = "($$(sys)) particle multiplicity";
          titles ttlCIMN("($$(sys)) CI neutrons","KE_{CM} (MeV)",Form("#frac{dM(%d)_{n,CI}}{dKE_{CM} #Delta#Omega_{CM}}",sys));
          titles ttlCIMP("($$(sys)) CI protons","KE_{CM} (MeV)",Form("#frac{dM(%d)_{p,CI}}{dKE_{CM} #Delta#Omega_{CM}}",sys));
          titles ttlCIRs("($$(sys)) CI n/p","KE_{CM} (MeV)",Form("R(%d) = #frac{dM(%d)_{n,CI}}{dKE_{CM} #Delta#Omega_{CM}} / #frac{dM(%d)_{p,CI}}{dKE_{CM} #Delta#Omega_{CM}}",sys,sys,sys));
          vector<TH1*> hists;
          TString cname_ci_sys = var.name+"_"+sys+"__"+vshort;
          for (auto ipid : ipids) {
            setpar_syspid(isys,ipid);
            auto hist = var.draw(tree_pid[isys][ipid]);
            hists.push_back(hist);
            ejungwoo::add(cname_ci_sys, 1, hist, "histl logy grid", particleNames[ipid]);
          /****************************************************************************************************************/
          }
          setpar("pname","CI");
          auto graphs = draw_ci(hists, dsd, isys);
          /****************************************************************************************************************/
          ejungwoo::add(cname_ci_sys, 0, new_h(var.histname+"_cin_frame", ttlCIMN, binnCIx, binnCIy), "addx rangex logy gridy");
          ejungwoo::add(cname_ci_sys, 0, graphs.At(0), "colorx p", "CI neutrons");
          /****************************************************************************************************************/
          ejungwoo::add(cname_ci_sys, 2, new_h(var.histname+"_cip_frame", ttlCIMP, binnCIx, binnCIy), "addx rangex logy gridy");
          ejungwoo::add(cname_ci_sys, 2, graphs.At(1), "colorx p", "CI protons");
          /****************************************************************************************************************/
          ejungwoo::add(cname_ci_sys, 3, new_h(var.histname+"_cir_frame", ttlCIRs, binnCIx, binnCIr), "addx rangex gridy");
          ejungwoo::add(cname_ci_sys, 3, graphs.At(2), "colorx p", "CI n/p");
          /****************************************************************************************************************/

          ejungwoo::add(cname_ci, 0, graphs.At(0), "colorx l", Form("(%d+%d) CI neutrons",sys,tga));
          ejungwoo::add(cname_ci, 1, graphs.At(1), "colorx l", Form("(%d+%d) CI protons",sys,tga));
          ejungwoo::add(cname_ci, 2, graphs.At(2), "colorx l", Form("(%d+%d) CI n/p",sys,tga));

          graphs_rnp[isys] = (TGraphErrors *) graphs.At(2);
        }

        titles ttlCID("Double ratio","KE_{CM} (MeV)","DR(n/p)_{#frac{132}{108}} = #frac{R(132)}{R(108)}");
        ejungwoo::add(var.name+"_r132o108__"+vshort, new_h(var.histname+"_cidr_allframe", ttlCID, binnCIx, binnCIdr), "addx rangex gridx");
        ejungwoo::add(var.name+"_r132o108__"+vshort, draw_ratio(graphs_rnp[0], graphs_rnp[1], binnCIx.max), "p addx");
      }
    }
  }

  /******************************************************************************************/

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

TObjArray draw_ci(vector<TH1 *>hists, double ds, int isys)
{
  auto graph_np = new TGraphErrors();
  auto graph_nn = new TGraphErrors();
  auto graph_rr = new TGraphErrors();

  binning binn(hists[0]);
  binn.resetb();
  while (binn.nextb()) {
    binn.idx;
    double vpsum = 0;
    double vnsum = 0;
    for (auto ipid : ipids) {
      double ntracks = hists[ipid] -> GetBinContent(binn.idx);
      auto np = numProtons[ipid];
      auto nn = numNeutrons[ipid];
      vpsum += numProtons[ipid] * ntracks;
      vnsum += numNeutrons[ipid] * ntracks;
    }
    graph_np -> SetPoint(graph_np->GetN(), binn.value, vpsum/binn.w/ds);
    graph_nn -> SetPoint(graph_nn->GetN(), binn.value, vnsum/binn.w/ds);
    graph_rr -> SetPoint(graph_rr->GetN(), binn.value, vnsum/vpsum);
  }

  TObjArray array;
  int idx = 0;
  for (auto graph : {graph_nn, graph_np, graph_rr}) {
    graph -> SetMarkerStyle(ejungwoo::markeri(idx));
    graph -> SetMarkerColor(ejungwoo::colori(idx));
    graph -> SetLineColor(ejungwoo::colori(isys));
    graph -> SetLineWidth(2);
    array.Add(graph);
    idx++;
  }

  return array;
}

TGraphErrors *draw_ratio(TGraphErrors *graph1, TGraphErrors *graph2, double max)
{
  auto graph_dr = new TGraphErrors();

  double x1, y1, x2, y2;
  auto nn = graph1 -> GetN();
  for (auto i=0; i<nn; ++i) {
    graph1 -> GetPoint(i,x1,y1);
    graph2 -> GetPoint(i,x2,y2);
    if (x1 > max)
      break;
    if (y2>0) {
      graph_dr -> SetPoint(graph_dr->GetN(),x1,y1/y2);
    }
  }

  return graph_dr;
}
