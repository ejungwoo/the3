#include "/Users/ejungwoo/config/ejungwoo.h"
#include "init_variables.h"

int fParIndex1[2] = {0,2};
int fParIndex2[2] = {1,2};

struct GlobalChi2
{
  const ROOT::Math::IMultiGenFunction *fChi21;
  const ROOT::Math::IMultiGenFunction *fChi22;

  double operator() (const double *parIn) const {
    double par1[2] = {parIn[0], parIn[2]};
    double par2[2] = {parIn[1], parIn[2]};
    return (*fChi21)(par1) + (*fChi22)(par2);
  }

  GlobalChi2(ROOT::Math::IMultiGenFunction &chi21, ROOT::Math::IMultiGenFunction &chi22) : fChi21(&chi21), fChi22(&chi22) {}
};

TH1D *hist_dgraph(TString name, TTree *tree, TGraph *graph, double x1, double x2, double tta1, double tta2, int num_events);

double fPol1Function(double *x, double *par) {
  double value = TMath::Exp(par[0] + par[1]*x[0]);
  return value;
}

bool aaaaaaaa = 0;
bool docorr = 1;

void draw_all()
{
 ejungwoo::gfixcvsx(100);
 ejungwoo::gfixcvsy(100);

  //ejungwoo::gcvspos(1000);
  //ejungwoo::gcvspos(1300);
  ejungwoo::gstat(0);
  ejungwoo::gsave(0);
  ejungwoo::gsetvmark(0);
  //ejungwoo::gshortprint();
  //ejungwoo::gdummytp();

  //auto chooseVersion = "right_55";
  //auto chooseVersion = "right_50";
  //auto chooseVersion = "all_55";
  //auto chooseVersion = "all_45";

  //auto chooseVersion = "y21_left_45";
  //auto chooseVersion = "y21_left_55";
  //auto chooseVersion = "y21_right_45";
  //auto chooseVersion = "y21_right_55";

  //auto chooseVersion = "fix5_left_45";
  //auto chooseVersion = "fix5_left_55";
  //auto chooseVersion = "fix5_right_45";
  //auto chooseVersion = "fix5_right_55";

  auto chooseVersion = "fix6_left_45";
  //auto chooseVersion = "fix6_left_55";
  //auto chooseVersion = "fix6_right_45";
  //auto chooseVersion = "fix6_right_55";

  auto pidFileName = "PIDSigma_Sn132KanekoMult50.root";
  auto filePID = new TFile(pidFileName);

  int  pidmode = 3;
  bool addanddraw = 1;
  bool guideline = 1;

  bool draw_pydistpid = 0;

  int draw_plab = 0;
  int draw_pydist = 0;
  int draw_pndist = 0;
  int draw_pid = 1;
  int draw_like = 0;
  int draw_npratio = 1;
  int draw_dbratio = 0;
  int draw_r21 = 1;
  int draw_r21like = 0;
  int draw_r21nz = 0;
  int draw_srn0 = 0;
  int draw_psdon = 0;
  int draw_cici = 0;
  int draw_dist = 0;
  int draw_temp = 1;

  auto remove_draw = [&draw_pydist, &draw_pndist, &draw_pid, &draw_like, &draw_npratio, &draw_dbratio, &draw_r21, &draw_srn0, &draw_psdon, &draw_cici, &draw_dist]() {
    draw_pydist = 0; draw_pndist = 0; draw_pid = 0; draw_like = 0; draw_npratio = 0; draw_dbratio = 0; draw_r21 = 0; draw_srn0 = 0; draw_psdon = 0; draw_cici = 0; draw_dist = 0;
  };

  //remove_draw(); draw_temp = 1;

  double num_tracks_per_event_cut = 0.001;
  //double num_tracks_per_event_cut = 0.000000000001;

  init();

  auto filegraph1 = new TFile("data4/graph1.y21_right_45_yy_0_5.root");
  auto filegraph2 = new TFile("data4/graph2.y21_right_45_yy_0_5.root");
  auto filegraph3 = new TFile("data4/graph3.y21_right_45_yy_0_5.root");
  auto filegraph4 = new TFile("data4/graph4.y21_right_45_yy_0_5.root");
  auto filegraph5 = new TFile("data4/graph5.y21_right_45_yy_0_5.root");
  auto filegraph6 = new TFile("data4/graph6.y21_right_45_yy_0_5.root");

  auto graph1 = (TGraph *) filegraph1 -> Get("graph1");
  auto graph2 = (TGraph *) filegraph2 -> Get("graph2");
  auto graph3 = (TGraph *) filegraph3 -> Get("graph3");
  auto graph4 = (TGraph *) filegraph4 -> Get("graph4");
  auto graph5 = (TGraph *) filegraph5 -> Get("graph5");
  auto graph6 = (TGraph *) filegraph6 -> Get("graph6");

  //auto hist12 = new TH2D("hist","",100,0,2500,100,0,1500);
  auto graph_cut1 = ejungwoo::make_cutg("gcut0",graph1, graph2); graph_cut1 -> SetLineStyle(1);
  auto graph_cut2 = ejungwoo::make_cutg("gcut1",graph2, graph3); graph_cut2 -> SetLineStyle(1);
  auto graph_cut3 = ejungwoo::make_cutg("gcut2",graph3, graph4); graph_cut3 -> SetLineStyle(1);
  auto graph_cut4 = ejungwoo::make_cutg("gcut3",graph4, graph5); graph_cut4 -> SetLineStyle(1);
  auto graph_cut5 = ejungwoo::make_cutg("gcut4",graph5, graph6); graph_cut5 -> SetLineStyle(1);

  for (auto graph_cut : {graph_cut1,graph_cut2,graph_cut3,graph_cut4,graph_cut5}) {
    graph_cut -> SetVarX("p_lab");
    graph_cut -> SetVarY("dedx");
  }



    //if (!docorr) setpar("cut_pe", "($$(mult_cut))*($$(pes_cut))*($$(poz_cut))*($$(angle_cut))*($$(rap_cut))");

  //setpar("cut1","$$(cut0)");
  setpar("cut1","$$(cut_pek)");
  setpar("cut2","$$(cut_noy)");
  //setpar("cut2","$$(prob)/$$(eff)/$$(num_events)*($$(mult_cut))*($$(pes_cut))*($$(poz_cut))*($$(angle_cut))*($$(rap_cut))");
  //setpar("cut2","$$(prob)/$$(eff)/$$(num_events)*($$(mult_cut))*($$(pes_cut))*($$(poz_cut))*($$(angle_cut))*($$(rap_cut))");
  //setpar("cut_pek","$$(prob)/effk/$$(num_events)*($$(mult_cut))*($$(rap_cut))*($$(pes_cut))*($$(poz_cut))*($$(angle_cut))");
  setpar("cut_pek","$$(prob)/$$(eff)/$$(num_events)*($$(mult_cut))*($$(rap_cut))*($$(pes_cut))*($$(poz_cut))*($$(angle_cut))");
  //setpar("cut_pek","$$(gcutpid)*($$(tta_lab)>40&&$$(tta_lab)<60)");

  const char *particleNames[] = {"p","d","t","he3","he4","he6","1"};

  int ipid_nl[3][2] = {{4,3},{1,0},{2,1}}; int num_nl = 3;
  int ipid_pl[3][2] = {{0,6},{4,2},{3,1}}; int num_pl = 3;
  int ipid_dl[3][2] = {{1,6},{3,0},{4,1}}; int num_dl = 3;
  int ipid_tl[2][2] = {{2,6},{4,0}}; int num_tl = 2;

  //int isystemsAll[] = {0,1,2,3};
  int isystemsAll[] = {0,1};
  //int isystemsAna[] = {0,1,2,3};
  int isystemsAna[] = {0,1};
  auto findisys = [isystemsAna](int isys) {
    for (auto isys0 : isystemsAna)
      if (isys0==isys)
        return true;
     return false;
  };

  const int ipidsAna[] = {0,1,2,3,4};
  //const int ipidsAna[] = {3};
  //const int ipidsAna[] = {0,1,2};
  auto findipid = [ipidsAna](int ipid) {
    if (ipid==6)
      return true;
    for (auto ipid0 : ipidsAna)
      if (ipid0==ipid)
        return true;
     return false;
  };

  auto useSingleXBinning = false;
  int yy1, yy2, pt1, pt2;

  if (0) {
    yy1 = -10; yy2 = 20;
    //yy1 = 0; yy2 = 4;
    //yy1 = -4; yy2 = 0;
    //pt1 = 100; pt2 = 200;
    //pt1 = 200; pt2 = 400;
    //pt1 = 400; pt2 = 600;
    pt1 = 600; pt2 = 800;
    //remove_draw(); draw_pydist = 1; draw_r21nz = 2;
    useSingleXBinning = true;
  }
  else {
    //yy1 = 0; yy2 = 4;
    //yy1 = 6; yy2 = 10;
    //yy1 = -2; yy2 = 2;
    yy1 = 0; yy2 = 5;
    //yy1 = 5; yy2 = 10;
    pt1 = -1; pt2 = -1;
    useSingleXBinning = false;
  }

  if (pt1<0) {
    pt1 = 0;
    pt2 = 200;
  }

  //int pt1 = -1; int pt2 = -1;
  bool useRapidityBinning = 0;
  //const char *yBeam = "by_cm0";
  const char *yBeam = "by_cm/2";

  //remove_draw(); draw_pydist = 1; draw_pid = 1;

  for (auto aversion : {chooseVersion})
  {
    auto condition_array = setversion(aversion);

    auto vpozlab = variable("poz_lab", "p_lab", "$$(cut2)" ,titles(fttl_main, "p_{Lab}/Z (MeV/c^{2})", fttly_ntne), binning(400,0,2000));
    auto vpydist = variable("pydist", Form("pt_cm/$$(a):fy_cm/(%s)",yBeam), "$$(cut2)", titles(" ", "y_{0}", "p_{T}/A (MeV/c)"), binning(100,-1.,2.), binning(100,0,1000));
    //auto vpydist = variable("pydist", Form("pt_cm/$$(a):fy_cm/(%s)",yBeam), "$$(cut2)", titles(" ", "y_{0}", "p_{T}/A (MeV/c)"), binning(100,-4.,4), binning(100,0,1000));
    auto vpydist2 = variable("pydist2", Form("fy_cm/(%s)",yBeam), "$$(cut2)", titles(" ", "y_{0}", "#frac{d^{2}M}{#Delta#Omega_{CM}dy}"), binning(100,-1.,2.));
    //auto vpndist = variable("pndist", "$$(nr):pt_cm/$$(a)", "($$(nr)!=0)*($$(cut_x))", titles(" ", "p_{T}/A (MeV/c)", "nc"), binning(100,0,400), binning(100,0,100));
    //auto vpndist = variable("pndist", "$$(nl):pt_cm/$$(a)", "($$(nl)!=0)*($$(cut_x))", titles(" ", "p_{T}/A (MeV/c)", "nc"), binning(100,0,400), binning(100,0,100));
    auto vpndist = variable("pndist", "($$(nr)+$$(nl)):pt_cm/$$(a)", "$$(cut_x)", titles(" ", "p_{T}/A (MeV/c)", "nc"), binning(100,0,400), binning(100,0,100));
    //auto vpid2 = variable("pid", "dedx:p_lab", "$$(cut1)*(dedx<600)", titles(" ","p_{Lab}/Z (MeV/c^{2})","dE/dx"), binning(400,0,2500),binning(400,0,600)); 
    //auto vpid2 = variable("pid", "dedx:p_lab", "$$(cut2)", titles(" ","p_{Lab}/Z (MeV/c^{2})","dE/dx"), binning(400,0,2500),binning(400,0,600)); 
    auto vpid2 = variable("pid", "dedx:p_lab", "$$(cut0)", titles(" ","p_{Lab}/Z (MeV/c^{2})","dE/dx"), binning(200,0,3000),binning(200,0,1000)); 
    //auto vpid2 = variable("pid", "dedx:p_lab", "$$(cut_pe)", titles(" ","p_{Lab}/Z (MeV/c^{2})","dE/dx"), binning(200,0,3000),binning(200,0,1000)); 

    const char *selpid = "";
         if (pidmode==1) selpid = "$$(prob)/$$(eff)/$$(num_events)";
    else if (pidmode==2) selpid = "(prob>.5)";
    else if (pidmode==3)  {
      selpid = "$$(prob)/$$(eff)/$$(num_events)*($$(tta_lab)>$$(tta1)&&$$(tta_lab)<$$(tta2))";
      if (!docorr) selpid = "($$(tta_lab)>$$(tta1)&&$$(tta_lab)<$$(tta2))";
    }
    //else if (pidmode==4) selpid = "(prob>.5)*($$(tta_lab)>$$(tta1)&&$$(tta_lab)<$$(tta2))";
    else if (pidmode==4) selpid = "$$(gcutpid)";
      vpid2 = variable("pid", "dedx:p_lab", selpid, titles(" ","p_{Lab}/Z (MeV/c^{2})","dE/dx"), binning(200,0,3000),binning(200,0,1000)); 
    //auto vpid2 = variable("pid", "dedx:pt_cm/$$(a)", "$$(cut0)*(dedx<600)", titles(" ","p_{T}/A (MeV/c)","dE/dx"), binning(400,0,1000),binning(400,0,600)); 

    binning binn0;
    const char *xTitle;
    const char *xVar;
    variable vxvar;

    if (useRapidityBinning) {
      binn0 = binning(8,0.0,1.5);
      if (useSingleXBinning)
        binn0 = binning(1,yy1,yy2);
      xTitle = "y_{0} = 2y_{z}/y_{beam Lab.}";
      xVar = "y_{0}";
      vxvar = variable("y0", Form("fy_cm/(%s)",yBeam), "$$(cut1)", titles("",xTitle,""), binn0);
      if (useSingleXBinning) {
        vxvar = variable("y0", Form("fy_cm/(%s)",yBeam), "$$(cut1)", titles("",xTitle,""), binn0);
      }
      yy1 = -10;
      yy2 = 20;
    } else {
      //binn0 = binning(10,0,1000);
      //binn0 = binning(5,0,250);
      binn0 = binning(8,0,400);
      if (useSingleXBinning)
        binn0 = binning(1,pt1,pt2);
      xTitle = "p_{T}/A (MeV/c)";
      xVar = "p_{T}/A";
      vxvar = variable("ptoa_cm", "pt_cm/$$(a)", "$$(cut1)", titles("","p_{T}/A (MeV/c)",""), binn0);
    }

    auto ccc = 1.;
    setpar("rap_cut",Form("$$(foby_cm)>%.1f && $$(foby_cm)<%.1f",.1*yy1,.1*yy2));
    if (pt1!=0&&pt2!=400)
      setpar("rap_cut",Form("$$(foby_cm)>%.1f && $$(foby_cm)<%.1f && $$(ptoa_cm)>%d &&  $$(ptoa_cm)<%d",.1*yy1,.1*yy2,pt1,pt2));
    TString maintitle = Form("%s, Mult=%d~%d,  y_{0} = %.1f ~ %.1f",fHeadName.Data(), fMultLLCut[1],fMultHLCut[1], .1*yy1,.1*yy2);
    if (pt1!=0&&pt2!=400) maintitle = maintitle + Form(",  p_{T}/A = %d ~ %d",pt1, pt2);
    if (useSingleXBinning)
      fVOutShort = fVOutShort + "_single";
    TString vout2 = fVOutShort+"_yy_"+yy1+"_"+yy2;
    if (pt1!=0&&pt2!=400)
      vout2 = fVOutShort+Form("_pt%03d%03d_yy%02d%02d",pt1,pt2,yy1,yy2);

    cout << vout2 << endl;

    auto fnpdtlike = variable("fnpdtlike", "", "", titles(" ",xTitle,"Y(particle-0)/Y(particle-1)"), binn0, binning(100,0,5));
    auto fnpratio = variable("fnpratio", "", "", titles(" ",xTitle,"n-like/p-like"), binn0, binning(100,0,20));
    auto fdbratio = variable("fdbratio", "", "", titles(" ",xTitle,"DR(n-like/p-like)$$(systga21)"), binn0, binning(100,0,2.5));
    auto fdbratio2 = variable("fdbratio2", "", "", titles(" ",xTitle,"DR(#frac{CI-n}{CI-p})"), binn0, binning(100,0,2.0));
    //auto fptoar21 = variable("fptoar21", "", "", titles(" ",xTitle,"R_{21}"), binn0, binning(100,0.5,2.));
    auto fptoar21 = variable("fptoar21", "", "", titles(" ",xTitle,"R_{21}"), binn0, binning(100,0,2.));
    //auto fnr21 = variable("fnr21", "", "", titles(" ","n","R21($$(sys1)/$$(sys2))"), binning(4,-1,3), binning(100,0.6,4));
    //auto fzr21 = variable("fzr21", "", "", titles(" ","z","R21($$(sys1)/$$(sys2))"), binning(4,0,3), binning(100,0.6,4)); bool setfznzlog = true;
    auto fnr21 = variable("fnr21", "", "", titles(" ","N","R_{21}($$(sys1)+$$(tar1)/$$(sys2)+$$(tar2))"), binning(4,-1,3), binning(100,0.5,2));
    auto fzr21 = variable("fzr21", "", "", titles(" ","Z","R_{21}($$(sys1)+$$(tar1)/$$(sys2)+$$(tar2))"), binning(4,-1,3), binning(100,0.5,2)); bool setfznzlog = false;
    auto fptoasr = variable("fptoasr", "", "", titles(" ",xTitle,"SR(t/he3)"), binn0, binning(100,0,7));
    auto fptoan0 = variable("fptoan0", "", "", titles(" ",xTitle,"pseudo neutron  #frac{d^{2}M}{dyd(p_{T}/A)}"), binn0, binning(100,0,0.2));
    auto fcidist0 = variable("cidist", "", "", titles(" ",xTitle,"#frac{d^{2}M}{dyd(p_{T}/A)}"), binn0, binning(100,0,ccc*0.055));
    auto fcidist = variable("cidist", "", "", titles(" ",xTitle,"#frac{d^{2}M}{dyd(p_{T}/A)}"), binn0, binning(100,0,ccc*0.15));
    auto fciratio = variable("fciratio", "", "", titles(" ",xTitle,"#frac{CI-n}{CI-p}"), binn0, binning(100,0,2.0));
    auto ftemp = variable("temp", "", "", titles(" ",xTitle,"T = #frac{14.3}{log[1.59*R_{He-H}]}"), binn0, binning(100,0,25));
    auto yrheh = Form("R_{He-H} = #scale[0.8]{#frac{dM_{d}/d#Omegad(%s) #times dM_{he4}/d#Omegad(%s)}{dM_{t}/d#Omegad(%s) #times dM_{he3}/d#Omegad(%s)}}",xVar,xVar,xVar,xVar);
    auto frheh = variable("rheh", "", "", titles(" ",xTitle,yrheh), binn0, binning(100,0,5.));
    auto ftratio = variable("tratio", "", "", titles(" ",xTitle,"T_{2} / T_{1}"), binn0, binning(100,0,1.5));

    if (draw_plab) {
      for (auto isys : {0}) {
        TString ename = TString("plab_")+fSystems[isys];
        for (int ipid : ipidsAna)
        {
          setpar_syspid(isys,ipid);
          vpozlab.setmaint(Form("sys%d %s",fSystems[isys],fParticleNames[ipid].Data()));
          auto hist = vpozlab.draw(fTreePID[isys][ipid]);
          ejungwoo::addnext(ename,hist,"gridx gridy rangex hist");
        }
      }
    }
    if (draw_plab==2) break;

    if (draw_pydist) {
      cout << "draw_pydist" << endl;
      for (auto isys : isystemsAna)
      //for (auto isys : {0})
      {
        TString ename = TString("pydist_sysCM")+fSystems[isys];
        if (draw_pydistpid) 
          ename = TString("pydistpid_sysCM")+fSystems[isys];
        TString ename2 = TString("ydist_sysCM")+fSystems[isys];
        //for (int ipid : ipidsAna)
        for (int ipid : {0})
        {
          setpar_syspid(isys,ipid);
          vpydist.setmaint(Form("%s  (%d+%d)",fParticleNames2[ipid].Data(),fSystems[isys],fTargetA[isys]));
          if (draw_pydist==3) ejungwoo::gdummytp();
          auto hist = vpydist.draw(fTreePID[isys][ipid]);
          if (draw_pydist==3) ejungwoo::gdummytp(false);
          ejungwoo::addnext(ename,hist,"gridx gridy rangex");
          if (1) {
            auto graph = new TGraph();
            graph -> SetLineColor(kRed);
            graph -> SetLineWidth(2);
            graph -> SetPoint(graph->GetN(), .1*yy1, pt1);
            graph -> SetPoint(graph->GetN(), .1*yy1, pt2);
            graph -> SetPoint(graph->GetN(), .1*yy2, pt2);
            graph -> SetPoint(graph->GetN(), .1*yy2, pt1);
            graph -> SetPoint(graph->GetN(), .1*yy1, pt1);
            ejungwoo::addsame(ename,graph,"colorx samel addx");
          }

          vpydist2.setmaint(Form("sys%d %s",fSystems[isys],fParticleNames[ipid].Data()));
          //auto hist2 = vpydist2.draw(fTreePID[isys][ipid]);
          //ejungwoo::addnext(ename2,hist2,"gridx gridy rangex hist");
        }
      }
    }
    if (draw_pydist==2) break;

    if (draw_pndist) {
      cout << "draw_pndist" << endl;
      for (auto isys : isystemsAna)
      //for (auto isys : {0})
      {
        TString ename = TString("pn")+fSystems[isys];
        for (int ipid : ipidsAna)
        //for (int ipid : {0})
        {
          setpar_syspid(isys,ipid);
          vpndist.setmaint(Form("sys%d %s",fSystems[isys],fParticleNames[ipid].Data()));
          auto hist = vpndist.draw(fTreePID[isys][ipid]);
          ejungwoo::addnext(ename,hist,"gridx gridy");
        }
      }
    }


    //for (auto ipid : ipidsAna) fSDHLCut[ipid] = 2.0;
    //vout2 = vout2 + "_sd2";
    if (useRapidityBinning)
      vout2 = vout2 + "_by";
    ejungwoo::gversionout(vout2);
    if (draw_pndist==2) break;

    if (draw_pid) {
      cout << "draw_pid" << endl;

      //for (auto isys : isystemsAna)
      //for (auto isys : {0})
      for (auto isys : {0,1})
      {
        int itheta1 = 0;
        int itheta2 = 4;
        if (pidmode<3)
          itheta2 = 1;
        if (pidmode==4) {
          itheta1 = 2;
          itheta2 = 3;
        }

        for (auto itheta=itheta1; itheta<itheta2; ++itheta)
        {
          TString ename = "pid";
          if (draw_pydistpid) 
            ename = TString("pydistpid_sysCM")+fSystems[isys];
          TH2D *hist2 = nullptr;

          int tta1 = itheta*20;
          int tta2 = (itheta+1)*20;

          setpar("tta1",tta1);
          setpar("tta2",tta2);


          if (pidmode==1) {
            TString maintitle2 = Form("%s, Mult=%d~%d", fHeadName.Data(), fMultLLCut[1],fMultHLCut[1]);
            vpid2.setmaint(TString("<corr.> sys") + fSystems[isys] + " " + maintitle2);
            //ename = TString("pidcor1_sys") + fSystems[isys];
            ename = TString("pidcor1_sys");// + fSystems[isys];
          } else if (pidmode==2) {
            TString maintitle2 = Form("%s, Mult=%d~%d", fHeadName.Data(), fMultLLCut[1],fMultHLCut[1]);
            vpid2.setmaint(TString("<raw> sys") + fSystems[isys] + " " + maintitle2);
            ename = TString("pidraw2_sys") + fSystems[isys];
          } else if (pidmode==3) {
            TString maintitle2 = Form("%s, Mult=%d~%d, #theta_{Lab}=%d~%d", fHeadName.Data(), fMultLLCut[1],fMultHLCut[1], tta1, tta2);
            vpid2.setmaint(TString("<corr.> sys") + fSystems[isys] + " " + maintitle2);
            ename = TString("pidraw3_sys") + fSystems[isys] + "_tta" + itheta;
            //ename = TString("pidcor1_sys")+ "_tta" + itheta;
          } else if (pidmode==4) {
            TString maintitle2 = Form("%s, Mult=%d~%d, #theta_{Lab}=%d~%d", fHeadName.Data(), fMultLLCut[1],fMultHLCut[1], tta1, tta2);
            vpid2.setmaint(TString("<raw> sys") + fSystems[isys] + " " + maintitle2);
            ename = TString("pidraw4_sys") + fSystems[isys] + "_tta" + itheta;
          } else {
            vpid2.setmaint(TString("sys") + fSystems[isys] + " " + maintitle);
            ename = TString("pid_sys") + fSystems[isys];
          }

          for (int ipid : ipidsAna) {
            setpar_syspid(isys,ipid);
            TTree *tree_pid = fTreePID[isys][ipid];
            auto histpid = (TH2D *) vpid2.draw(tree_pid);
            if (addanddraw)  {
              //if (ipid==0) hist2 = histpid;
              if (hist2==nullptr) hist2 = histpid;
              //else hist2 -> Add(vpid2.draw(fTreePID[isys][ipid]));
              else hist2 -> Add(histpid);
            }
            else {
              if (pidmode>0) ejungwoo::addnext(TString("pidraw")+isys,histpid);
              else ejungwoo::addnext(TString("pid")+isys,histpid);
            }

            if (guideline) {
              auto bee = (TF1 *) filePID -> Get(Form("BEE%d",fParticlePDGs[ipid]));

              auto graph = new TGraph();
              double mom1 = 200;
              double mom2 = 2500;
              double x1, x2;
                   if (ipid==0) { mom1 = 110; mom2 = 1300; x1 =  300; x2 = 1000; }
              else if (ipid==1) { mom1 = 250; mom2 = 1800; x1 =  500; x2 = 1500; }
              else if (ipid==2) { mom1 = 365; mom2 = 2500; x1 =  700; x2 = 2000; }
              //else if (ipid==3) { mom1 = 450; mom2 = 1500; x1 =  700; x2 = 1200; }
              else if (ipid==3) { mom1 = 450; mom2 = 1500; x1 =  700; x2 = 1000; }
              else if (ipid==4) { mom1 = 500; mom2 = 2000; x1 =  700; x2 = 1500; }

              for (double mom = mom1; mom <= mom2; mom+=10) {
                auto dedx = bee -> Eval(mom);
                graph -> SetPoint(graph->GetN(), mom, dedx);
              }

              auto hist = hist_dgraph(Form("dist_%s_%s",histpid->GetName(),fParticleNames[ipid].Data()),tree_pid,graph,x1,x2,tta1,tta2,fNumEvents[isys]);
              auto fit1 = ejungwoo::fitg(hist,1);
              ejungwoo::add(ename,ipid+1,hist,"colorx addx hist");
              ejungwoo::add(ename,ipid+1,fit1,"colorx hist",Form("peak = %f",fit1->GetParameter(1)));

              //ejungwoo::addsame(ename,graph,"samel colorx addx");
              ejungwoo::add(ename,0,graph,"samel colorx addx");
            }
          }

          if (addanddraw) {
            if (pidmode==0||pidmode==1) {
              hist2 -> SetMaximum(0.02);
              hist2 -> SetMinimum(0.00005);
            }
            //ejungwoo::addnext(ename,hist2,"logz");
            //ejungwoo::addsame(ename,hist2,"logz");
            ejungwoo::add(ename,0,hist2,"logz");
            //ejungwoo::addsame(ename,hist2->ProjectionX(),"hist");
            if (guideline) {
              if (0)
              for (int ipid : ipidsAna)
              {
                auto bee = (TF1 *) filePID -> Get(Form("BEE%d",fParticlePDGs[ipid]));

                auto graph = new TGraph();
                double mom1 = 200;
                double mom2 = 2500;
                if (ipid==0) { mom1 = 110; mom2 = 1300; }
                else if (ipid==1) { mom1 = 250; mom2 = 1800; }
                else if (ipid==2) { mom1 = 365; mom2 = 2500; }
                else if (ipid==3) { mom1 = 450; mom2 = 1500; }
                else if (ipid==4) { mom1 = 500; mom2 = 2000; }

                for (double mom = mom1; mom <= mom2; mom+=10) {
                  auto dedx = bee -> Eval(mom);
                  graph -> SetPoint(graph->GetN(), mom, dedx);
                }

                ejungwoo::addsame(ename,graph,"samel colorx addx");
              }
              //ejungwoo::addsame(ename,graph_cut1,"samel colorx addx");
              //ejungwoo::addsame(ename,graph_cut2,"samel colorx addx");
              //ejungwoo::addsame(ename,graph_cut3,"samel colorx addx");
              //ejungwoo::addsame(ename,graph_cut4,"samel colorx addx");
              //ejungwoo::addsame(ename,graph_cut5,"samel colorx addx");
            }
          }
          if (aaaaaaaa) break;
        }
          if (aaaaaaaa) break;
      }
          if (aaaaaaaa) break;
    }
    if (draw_pid==2) break;


    TH1 *histYield[4][5];
    for (auto isys : isystemsAll) {
      for (auto ipid : ipidsAna) {
        setpar_syspid(isys,ipid);
        histYield[isys][ipid] = vxvar.draw(fTreePID[isys][ipid]);
      }
    }

    TH1 *histYieldLike[4][4][3]; // isys, inpdt, icomb
    TGraphErrors *graphdrmean[4] = {0};
    TH1D *histdrmean[4] = {0};
    TF1 *fitydrmean[4] = {0};
    while (draw_r21 || draw_r21like || draw_r21nz || draw_like || draw_npratio || draw_dbratio)
    {
      cout << "draw_like" << endl;
      for (auto isys : isystemsAna)
      {
        for (auto inpdt : {0,1,2,3})
        //for (auto inpdt : {3})
        {
          TString ename = Form("npdtlike_sys%d",fSystems[isys]);

          fnpdtlike.setmaint(maintitle);
          if (draw_like && isys==0) {
            auto ecvs = ejungwoo::add(ename,fnpdtlike.draw());
            //ecvs -> legendlt();
          }

          for (auto icomb=0; icomb<3; ++icomb)
          {
            if (inpdt==3 && icomb==2)
              continue;
               
            int ipid1, ipid2;
                 if (inpdt==0) { ipid1 = ipid_nl[icomb][0]; ipid2 = ipid_nl[icomb][1]; }
            else if (inpdt==1) { ipid1 = ipid_pl[icomb][0]; ipid2 = ipid_pl[icomb][1]; }
            else if (inpdt==2) { ipid1 = ipid_dl[icomb][0]; ipid2 = ipid_dl[icomb][1]; }
            else if (inpdt==3) { ipid1 = ipid_tl[icomb][0]; ipid2 = ipid_tl[icomb][1]; }
            if (!findipid(ipid1)||!findipid(ipid2))
              continue;

            auto histl = (TH1 *) histYield[isys][ipid1] -> Clone();

            TString titleTail = Form(" %s/%s",particleNames[ipid1],particleNames[ipid2]);
            int mstyle, mcolor;
                 if (inpdt==0) { mstyle = 29; histl -> SetTitle(TString("n-like ")+titleTail); mcolor = (icomb==0?kBlack:(icomb==1?kGray+2:kGray)); }
            else if (inpdt==1) { mstyle = 20; histl -> SetTitle(TString("p-like ")+titleTail); mcolor = (icomb==0?kRed+1:(icomb==1?kOrange-3:kPink+9)); }
            else if (inpdt==2) { mstyle = 21; histl -> SetTitle(TString("d-like ")+titleTail); mcolor = (icomb==0?kBlue+1:(icomb==1?kAzure-3:kAzure+8)); }
            else if (inpdt==3) { mstyle = 22; histl -> SetTitle(TString("t-like ")+titleTail); mcolor = (icomb==0?kGreen+3:(icomb==1?kSpring+5:kTeal+2)); }
            histl -> SetMarkerStyle(mstyle);
            histl -> SetMarkerColor(mcolor);
            histl -> SetLineColor(mcolor);
            if (mstyle==29) histl -> SetMarkerSize(1.2);

            if (ipid2!=6) histl -> Divide((TH1 *) histYield[isys][ipid2]);
            if (draw_like && isys==0) ejungwoo::addsame(ename,histl,"plhist colorx");

            histYieldLike[isys][inpdt][icomb] = histl;
          }
        }
      }
      if (draw_like==2) break;

      if (!(draw_like || draw_npratio || draw_dbratio)) break;


      cout << "draw_npratio" << endl;
      int countNP[4] = {0};
      TH1D *histnp[4][9];
      for (auto isys : isystemsAna)
      {
        TString ename = Form("sys%d_n-like/p-like",fSystems[isys]);
        fnpratio.setmaint(maintitle);

        if (draw_npratio) ejungwoo::add(ename,fnpratio.draw());
        for (auto inl=0; inl<3; ++inl)
          for (auto ipl=0; ipl<3; ++ipl) {
            int ipid1 = ipid_nl[inl][0];
            int ipid2 = ipid_nl[inl][1];
            int ipid3 = ipid_pl[ipl][0];
            int ipid4 = ipid_pl[ipl][1];
            if (!findipid(ipid1)||!findipid(ipid2)||!findipid(ipid3)||!findipid(ipid4))
              continue;

            TString npname = Form("sys%d_nplike%d%d",fSystems[isys],inl,ipl);
            histnp[isys][countNP[isys]] = (TH1D *) histYieldLike[isys][0][inl] -> Clone(npname);
            histnp[isys][countNP[isys]] -> SetTitle(Form("%d. (%s/%s) / (%s/%s)",countNP[isys]+1,
                  particleNames[ipid_nl[inl][0]],
                  particleNames[ipid_nl[inl][1]],
                  particleNames[ipid_pl[ipl][0]],
                  particleNames[ipid_pl[ipl][1]]));
            histnp[isys][countNP[isys]] -> Divide(histYieldLike[isys][1][ipl]);
                 if (inl==0) histnp[isys][countNP[isys]] -> SetMarkerStyle(5);
            else if (inl==1) histnp[isys][countNP[isys]] -> SetMarkerStyle(27);
            else if (inl==2) histnp[isys][countNP[isys]] -> SetMarkerStyle(42);
            if (draw_npratio) ejungwoo::addsame(ename,histnp[0][countNP[isys]],"hist");
            countNP[isys]++;
          }
      }
      //if (draw_npratio==2) break;


      cout << "draw_dbratio" << endl;
      for (auto isyscom : {0,1,2,3})
      {
        int isys1, isys2;
        TString drawoption = "";
             if (isyscom==0) { isys1=0; isys2=1; drawoption=""; }
        else if (isyscom==1) { isys1=3; isys2=2; }
        else if (isyscom==2) { isys1=0; isys2=2; }
        else if (isyscom==3) { isys1=3; isys2=1; }
        if (!findisys(isys1)||!findisys(isys2))
          continue;

        TString ename = TString("DR")+isyscom;
        fdbratio.setmaint(maintitle);
        ejungwoo::setpar("systga21",Form("_{#frac{%d+%d}{%d+%d}}", fSystems[isys1],fTargetA[isys1],fSystems[isys2],fTargetA[isys2]));
        if (draw_dbratio) {
          auto ecvs = fdbratio.drawadd(ename,"gridxgridy");
          ecvs -> legendxx();
        }
        vector<TH1 *> histdrarray;
        for (auto idx=0; idx<countNP[0]; ++idx) {
          TString drname = Form("DR%d",idx);
          auto histdr = (TH1D *) histnp[isys1][idx] -> Clone(drname);
          histdr -> Divide(histnp[isys2][idx]);
          if (draw_dbratio) ejungwoo::addsame(ename,histdr,drawoption+"plhist");
          if (draw_dbratio && idx==0) ejungwoo::addsame(ename,histdr,drawoption+"plhist");
          histdrarray.push_back(histdr);
        }
        //if (draw_dbratio==2) break;


        graphdrmean[isyscom] = ejungwoo::new_ge();
        auto binn = binning(histdrarray.at(0));
        binn.resetb();
        histdrmean[isyscom] = (TH1D *) ejungwoo::new_h(TString("histdrmean")+isyscom,"",binn);
        while (binn.nextb()) {
          double yvalue = 0;
          for (auto hist : histdrarray) {
            auto value = hist -> GetBinContent(binn.idx);
            yvalue += value;
          }
          yvalue = yvalue/9.;
          graphdrmean[isyscom] -> SetPoint(graphdrmean[isyscom]->GetN(),binn.value,yvalue);
          histdrmean[isyscom] -> Fill(binn.value,yvalue);
        }
        fitydrmean[isyscom] = new TF1(Form("fitydrmean%d",isyscom),"[0]",binn.min,binn.max);
        graphdrmean[isyscom] -> Fit(fitydrmean[isyscom],"QR0N");

        graphdrmean[isyscom] -> SetLineColor(ejungwoo::colori(isyscom+4));
        graphdrmean[isyscom] -> SetMarkerColor(kGray);
        graphdrmean[isyscom] -> SetMarkerStyle(21);
        graphdrmean[isyscom] -> SetMarkerSize(2);
        graphdrmean[isyscom] -> SetLineStyle(9);
        graphdrmean[isyscom] -> SetLineWidth(2);
        histdrmean[isyscom] -> SetLineColor(ejungwoo::colori(isyscom+4));
        histdrmean[isyscom] -> SetMarkerColor(kGray);
        histdrmean[isyscom] -> SetMarkerStyle(21);
        histdrmean[isyscom] -> SetMarkerSize(2);
        histdrmean[isyscom] -> SetLineStyle(9);
        histdrmean[isyscom] -> SetLineWidth(2);
        fitydrmean[isyscom] -> SetLineColor(ejungwoo::colori(isyscom+4));
        fitydrmean[isyscom] -> SetLineStyle(9);
        fitydrmean[isyscom] -> SetLineWidth(2);

        if (draw_dbratio) ejungwoo::addsame(ename,histdrmean[isyscom],drawoption+"phist colorx second","average");
      }
      break;
    }
    if (draw_dbratio==2) break;

    if (draw_r21 || draw_r21like || draw_r21nz)
    {
      cout << "draw_r21" << endl;
      //for (auto isyscom : {0,1,2,3})
      for (auto isyscom : {0})
      {
        int isys1, isys2;
        if (isyscom==0) { isys1=0; isys2=1; }
        else if (isyscom==1) { isys1=3; isys2=2; }
        else if (isyscom==2) { isys1=0; isys2=2; }
        else if (isyscom==3) { isys1=3; isys2=1; }
        if (!findisys(isys1)||!findisys(isys2))
          continue;

        TString systga1 = Form("%d+%d",fSystems[isys1],fTargetA[isys1]);
        TString systga2 = Form("%d+%d",fSystems[isys2],fTargetA[isys2]);

        TString ename = Form("r21_%d_%d",fSystems[isys1],fSystems[isys2]);
        fptoar21.setmaint(maintitle);
        fptoar21.setyt(Form("R_{21}( #frac{%s}{%s} )",systga1.Data(),systga2.Data()));
        if (draw_r21) fptoar21.drawadd(ename,"gridxgridy");
        TGraph *gfree[6] = {0};
        for (int ipid : ipidsAna) {
          auto graphs = draw_is(ipid, histYield[isys1][ipid], histYield[isys2][ipid], num_tracks_per_event_cut);
          gfree[ipid] = ((TGraph *) graphs.At(0));
          gfree[ipid] -> SetMarkerSize(1);
          if (draw_r21) ejungwoo::add(ename, 0, gfree[ipid], "pl", fParticleNames[ipid]);
        }

        if (draw_r21nz) {
          TString ename2 = Form("r21nz_%d_%d",fSystems[isys1],fSystems[isys2]);

          bool zzzz = 1;

          setpar("sys1",fSystems[isys1]);
          setpar("sys2",fSystems[isys2]);

          setpar("tar1",fTargetA[isys1]);
          setpar("tar2",fTargetA[isys2]);

          if (setfznzlog) {
            fnr21.setmaint(maintitle); ejungwoo::add(ename2,0,fnr21.draw(),"logy gridx gridy") -> legendlt();
            if (zzzz) fzr21.setmaint(maintitle); ejungwoo::add(ename2,1,fzr21.draw(),"logy gridx gridy");
          }
          else {
            fnr21.setmaint(maintitle); ejungwoo::add(ename2,0,fnr21.draw(),"gridx gridy") -> legendlt();
            if (zzzz) fzr21.setmaint(maintitle); ejungwoo::add(ename2,1,fzr21.draw(),"gridx gridy");
          }

          auto graph_hn = new TGraph();
          graph_hn -> SetMarkerSize(1.5);
          graph_hn -> SetMarkerStyle(20);
          graph_hn -> SetMarkerColor(kBlack);
          auto graph_hen = new TGraph();
          graph_hen -> SetMarkerSize(1.5);
          graph_hen -> SetMarkerStyle(25);
          graph_hen -> SetMarkerColor(kRed);

          for (int ipid : {0,1,2}) {
            Double_t xdummy, r21Value;
            gfree[ipid] -> GetPoint(0,xdummy,r21Value);
            auto znumber = fNumProtons[ipid];
            auto nnumber = fNumNeutrons[ipid];
            graph_hn -> SetPoint(graph_hn->GetN(),nnumber,r21Value);
          }

          for (int ipid : {3,4}) {
            Double_t xdummy, r21Value;
            gfree[ipid] -> GetPoint(0,xdummy,r21Value);
            auto znumber = fNumProtons[ipid];
            auto nnumber = fNumNeutrons[ipid];
            graph_hen -> SetPoint(graph_hen->GetN(),nnumber,r21Value);
          }

          ejungwoo::add(ename2, 0, graph_hn, "colorx p", "Z=1");
          ejungwoo::add(ename2, 0, graph_hen, "colorx p", "Z=2");

          auto graph_n0z = new TGraph();
          graph_n0z -> SetMarkerSize(1.5);
          graph_n0z -> SetMarkerStyle(20);
          graph_n0z -> SetMarkerColor(kBlack);
          auto graph_n1z = new TGraph();
          graph_n1z -> SetMarkerSize(1.5);
          graph_n1z -> SetMarkerStyle(25);
          graph_n1z -> SetMarkerColor(kRed);
          auto graph_n2z = new TGraph();
          graph_n2z -> SetMarkerSize(1.5);
          graph_n2z -> SetMarkerStyle(22);
          graph_n2z -> SetMarkerColor(kBlue);

          for (int ipid : {0}) {
            Double_t xdummy, r21Value;
            gfree[ipid] -> GetPoint(0,xdummy,r21Value);
            auto znumber = fNumProtons[ipid];
            auto nnumber = fNumNeutrons[ipid];
            graph_n0z -> SetPoint(graph_n0z->GetN(),znumber,r21Value);
          }

          for (int ipid : {1,3}) {
            Double_t xdummy, r21Value;
            gfree[ipid] -> GetPoint(0,xdummy,r21Value);
            auto znumber = fNumProtons[ipid];
            auto nnumber = fNumNeutrons[ipid];
            graph_n1z -> SetPoint(graph_n1z->GetN(),znumber,r21Value);
          }

          for (int ipid : {2,4}) {
            Double_t xdummy, r21Value;
            gfree[ipid] -> GetPoint(0,xdummy,r21Value);
            auto znumber = fNumProtons[ipid];
            auto nnumber = fNumNeutrons[ipid];
            graph_n2z -> SetPoint(graph_n2z->GetN(),znumber,r21Value);
          }

          ejungwoo::add(ename2, 1, graph_n0z, "colorx p", "N=0");
          ejungwoo::add(ename2, 1, graph_n1z, "colorx p", "N=1");
          ejungwoo::add(ename2, 1, graph_n2z, "colorx p", "N=2");

          if (1) {
            {
              auto graph1 = graph_hn;
              auto graph2 = graph_hen;
              //if (inz==1) { graph1 = graphZ1; graph2 = graphZ2; }
              TF1 *fit1 = new TF1("fit1",fPol1Function,0,100,2);
              TF1 *fit2 = new TF1("fit2",fPol1Function,0,100,2);
              double fitRange1 = -.5;
              double fitRange2 = 2.5;
              ROOT::Math::WrappedMultiTF1 wfit1(*fit1,1);
              ROOT::Math::WrappedMultiTF1 wfit2(*fit2,1);
              ROOT::Fit::DataOptions option;
              ROOT::Fit::DataRange range1(fitRange1,fitRange2);
              ROOT::Fit::DataRange range2(fitRange1,fitRange2);
              ROOT::Fit::BinData data1(option, range1);
              ROOT::Fit::BinData data2(option, range2);
              ROOT::Fit::FillData(data1, graph1);
              ROOT::Fit::FillData(data2, graph2);
              ROOT::Fit::Chi2Function chi21(data1, wfit1);
              ROOT::Fit::Chi2Function chi22(data2, wfit2);
              GlobalChi2 globalChi2(chi21, chi22);
              ROOT::Fit::Fitter fitter;
              vector<ROOT::Fit::ParameterSettings> parSetting = {
                ROOT::Fit::ParameterSettings("intercepti1", 0, 0.001, -100, 100),
                ROOT::Fit::ParameterSettings("intercepti2", 0, 0.001, -100, 100),
                ROOT::Fit::ParameterSettings("slope"      , 0, 0.001, -100, 100)};
              fitter.Config().SetParamsSettings(parSetting);
              fitter.Config().MinimizerOptions().SetPrintLevel(0);
              fitter.Config().SetMinimizer("Minuit","Minimize");
              fitter.FitFCN(3, globalChi2, 0, data1.Size()+data2.Size(), true);
              ROOT::Fit::FitResult result = fitter.Result();
              //result.Print(std::cout);
              fit1 -> SetFitResult(result, fParIndex1);
              fit1 -> SetRange(range1().first, range1().second);
              fit1 -> SetLineColor(kRed);
              graph1 -> GetListOfFunctions() -> Add(fit1);
              fit2 -> SetFitResult(result, fParIndex2);
              fit2 -> SetRange(range2().first, range2().second);
              fit2 -> SetLineColor(kRed);
              graph2 -> GetListOfFunctions() -> Add(fit2);
              //if (inz==0) { beta  = fit2 -> GetParameter(1); fitN1 = fit2; }
              //if (inz==1) { alpha = fit2 -> GetParameter(1); fitZ1 = fit2; }
              ejungwoo::add(ename2, 0, fit1, "colorx l addx");
              ejungwoo::add(ename2, 0, fit2, "colorx l", Form("#alpha=%.2f",fit2->GetParameter(1)));
              cout << "alpha = " << fit2 -> GetParameter(1) << endl;
            }
            {
              auto graph3 = graph_n1z;
              auto graph4 = graph_n2z;
              //if (inz==1) { graph3 = graphZ1; graph4 = graphZ2; }
              TF1 *fit3 = new TF1("fit3",fPol1Function,0,100,2);
              TF1 *fit4 = new TF1("fit4",fPol1Function,0,100,2);
              double fitRange1 = -.5;
              double fitRange2 = 2.5;
              ROOT::Math::WrappedMultiTF1 wfit1(*fit3,1);
              ROOT::Math::WrappedMultiTF1 wfit2(*fit4,1);
              ROOT::Fit::DataOptions option;
              ROOT::Fit::DataRange range1(fitRange1,fitRange2);
              ROOT::Fit::DataRange range2(fitRange1,fitRange2);
              ROOT::Fit::BinData data1(option, range1);
              ROOT::Fit::BinData data2(option, range2);
              ROOT::Fit::FillData(data1, graph3);
              ROOT::Fit::FillData(data2, graph4);
              ROOT::Fit::Chi2Function chi21(data1, wfit1);
              ROOT::Fit::Chi2Function chi22(data2, wfit2);
              GlobalChi2 globalChi2(chi21, chi22);
              ROOT::Fit::Fitter fitter;
              vector<ROOT::Fit::ParameterSettings> parSetting = {
                ROOT::Fit::ParameterSettings("intercepti1", 0, 0.001, -100, 100),
                ROOT::Fit::ParameterSettings("intercepti2", 0, 0.001, -100, 100),
                ROOT::Fit::ParameterSettings("slope"      , 0, 0.001, -100, 100)};
              fitter.Config().SetParamsSettings(parSetting);
              fitter.Config().MinimizerOptions().SetPrintLevel(0);
              fitter.Config().SetMinimizer("Minuit","Minimize");
              fitter.FitFCN(3, globalChi2, 0, data1.Size()+data2.Size(), true);
              ROOT::Fit::FitResult result = fitter.Result();
              //result.Print(std::cout);
              fit3 -> SetFitResult(result, fParIndex1);
              fit3 -> SetRange(range1().first, range1().second);
              fit3 -> SetLineColor(kRed);
              graph3 -> GetListOfFunctions() -> Add(fit3);
              fit4 -> SetFitResult(result, fParIndex2);
              fit4 -> SetRange(range2().first, range2().second);
              fit4 -> SetLineColor(kRed);
              graph4 -> GetListOfFunctions() -> Add(fit4);
              //if (inz==0) { beta  = fit4 -> GetParameter(1); fitN1 = fit4; }
              //if (inz==1) { alpha = fit4 -> GetParameter(1); fitZ1 = fit4; }
              ejungwoo::add(ename2, 1, fit3, "colorx l addx");
              ejungwoo::add(ename2, 1, fit4, "colorx l", Form("#beta=%.2f",fit4->GetParameter(1)));
              cout << "beta = " << fit4 -> GetParameter(1) << endl;
            }
          }

        }
        if (draw_r21nz==2) continue;

        if (draw_r21like)
        {
          ename = Form("r21_like_%d_%d",fSystems[isys1],fSystems[isys2]);
          fptoar21.setmaint(maintitle);
          fptoar21.setyt(Form("(npdt)-like R_{21}( #frac{%s}{%s} )",systga1.Data(),systga2.Data()));
          fptoar21.drawadd(ename,"gridxgridy");

          for (int inpdt : {0,1,2,3}) {
            //if (inpdt==1) { gfree[0]->SetLineColor(kGray); gfree[0]->SetMarkerColor(kRed); gfree[0]->SetMarkerStyle(24); gfree[0]->SetMarkerSize(2); ejungwoo::add(ename,0,gfree[0],"pl colorx", "p"); }
            //if (inpdt==2) { gfree[1]->SetLineColor(kGray); gfree[1]->SetMarkerColor(kRed); gfree[1]->SetMarkerStyle(25); gfree[1]->SetMarkerSize(2); ejungwoo::add(ename,0,gfree[1],"pl colorx", "d"); }
            //if (inpdt==3) { gfree[2]->SetLineColor(kGray); gfree[2]->SetMarkerColor(kRed); gfree[2]->SetMarkerStyle(26); gfree[2]->SetMarkerSize(2); ejungwoo::add(ename,0,gfree[2],"pl colorx", "t"); }
            for (auto icomb : {0,1,2}) {
              if (inpdt==3 && icomb==2)
                continue;

              int ipid1, ipid2;
              if (inpdt==0) { ipid1 = ipid_nl[icomb][0]; ipid2 = ipid_nl[icomb][1]; }
              else if (inpdt==1) { ipid1 = ipid_pl[icomb][0]; ipid2 = ipid_pl[icomb][1]; }
              else if (inpdt==2) { ipid1 = ipid_dl[icomb][0]; ipid2 = ipid_dl[icomb][1]; }
              else if (inpdt==3) { ipid1 = ipid_tl[icomb][0]; ipid2 = ipid_tl[icomb][1]; }
              if (!findipid(ipid1)||!findipid(ipid2))
                continue;

              //auto graphs = draw_is(0, histYieldLike[isys1][inpdt][icomb], histYieldLike[isys2][inpdt][icomb], num_tracks_per_event_cut);
              auto graphs = draw_is(0, histYieldLike[isys1][inpdt][icomb], histYieldLike[isys2][inpdt][icomb], num_tracks_per_event_cut);
              auto graph0 = (TGraphErrors *) graphs.At(0);

              TString graphTitle = Form("%s/%s",particleNames[ipid1],particleNames[ipid2]);
              int mstyle, mcolor;
              if (inpdt==0) { mstyle = 29; graphTitle = TString("n-like ")+graphTitle; mcolor = (icomb==0?kBlack:(icomb==1?kGray+2:kGray)); }
              else if (inpdt==1) { mstyle = 20; graphTitle = TString("p-like ")+graphTitle; mcolor = (icomb==0?kRed+1:(icomb==1?kOrange-3:kPink+9)); }
              else if (inpdt==2) { mstyle = 21; graphTitle = TString("d-like ")+graphTitle; mcolor = (icomb==0?kBlue+1:(icomb==1?kAzure-3:kAzure+8)); }
              else if (inpdt==3) { mstyle = 22; graphTitle = TString("t-like ")+graphTitle; mcolor = (icomb==0?kGreen+3:(icomb==1?kSpring+5:kTeal+2)); }
              if (graphTitle=="p-like p/1") graphTitle = "p";
              if (graphTitle=="d-like d/1") graphTitle = "d";
              if (graphTitle=="t-like t/1") graphTitle = "t";
              graph0 -> SetMarkerSize(1);
              if (mstyle==29) graph0 -> SetMarkerSize(1.2);

              if (inpdt==1&&icomb==0) { mstyle = 24; graph0 -> SetMarkerSize(1.5); }
              if (inpdt==2&&icomb==0) { mstyle = 25; graph0 -> SetMarkerSize(1.5); }
              if (inpdt==3&&icomb==0) { mstyle = 26; graph0 -> SetMarkerSize(1.5); }
              graph0 -> SetMarkerStyle(mstyle);
              graph0 -> SetMarkerColor(mcolor);
              graph0 -> SetLineColor(mcolor);

              auto ecvs = ejungwoo::add(ename, 0, graph0, "pl colorx", graphTitle);
              ecvs -> legendrr();
            }
          }
        }
      }
    }
    if (draw_r21==2) break;
    if (draw_r21nz==2) break;


    TH1 *hist_n0[4];
    if (draw_cici || draw_psdon)
    {
      if (!findipid(2)||!findipid(3)||!findipid(0))
        continue;

      cout << "draw_psdon" << endl;
      for (auto isys : isystemsAna) {
        TString systga = Form("%d+%d",fSystems[isys],fTargetA[isys]);

        setpar_syspid(isys,2);
        auto hist_sr = (TH1 *) histYield[isys][2] -> Clone();
        hist_sr -> SetMarkerStyle(20);
        hist_sr -> SetMarkerColor(ejungwoo::colori(isys));
        hist_sr -> SetMarkerSize(1);
        hist_sr -> SetLineColor(ejungwoo::colori(isys));
        hist_sr -> Divide(histYield[isys][3]);
        fptoasr.setmaint(maintitle);
        if (draw_srn0) fptoasr.drawadd("sr");
        if (draw_srn0) ejungwoo::addsame("sr",hist_sr,"colorx plhist",systga);

        setpar_syspid(isys,2);
        hist_n0[isys] = (TH1 *) histYield[isys][2] -> Clone();
        hist_n0[isys] -> SetMarkerStyle(20);
        hist_n0[isys] -> SetMarkerColor(ejungwoo::colori(isys));
        hist_n0[isys] -> SetLineColor(ejungwoo::colori(isys));
        hist_n0[isys] -> SetMarkerSize(1);
        hist_n0[isys] -> Divide(histYield[isys][3]);
        hist_n0[isys] -> Multiply(histYield[isys][0]);
        hist_n0[isys] = ejungwoo::dndx(hist_n0[isys]);
        fptoan0.setmaint(maintitle);
        if (draw_psdon) fptoan0.drawadd("n0");
        if (draw_psdon) ejungwoo::addsame("n0",hist_n0[isys],"colorx plhist",systga);
      }
    }
    if (draw_psdon==2) break;


    if (draw_cici)
    {
      cout << "draw_cici" << endl;
      //ejungwoo::add(ename, 0, new_h(vxvar.getHistName()+"_cin_allframe", ttlCIM, binnCIx, binnCIy2), "addx rangex logy ");
      //ejungwoo::add(ename, 1, new_h(vxvar.getHistName()+"_cip_allframe", ttlCIM, binnCIx, binnCIy2), "addx rangex logy ");
      //ejungwoo::add(ename, 2, new_h(vxvar.getHistName()+"_cir_allframe", ttlCIR, binnCIx, binnCIr), "addx rangex ");

      TString ename;

      TGraphErrors *graphs_rnp[4];
      for (auto isys : isystemsAna) {
        auto sys = fSystems[isys];
        auto tga = fTargetA[isys];
        vxvar.setmaint("($$(sys)) particle multiplicity");
        titles ttlCIMN("($$(sys)) CI neutrons","KE_{CM} (MeV)",Form("#frac{dM(%d)_{n,CI}}{dKE_{CM} #Delta#Omega_{CM}}",sys));
        titles ttlCIMP("($$(sys)) CI protons","KE_{CM} (MeV)",Form("#frac{dM(%d)_{p,CI}}{dKE_{CM} #Delta#Omega_{CM}}",sys));
        titles ttlCIRs("($$(sys)) CI n/p","KE_{CM} (MeV)",Form("R(%d) = #frac{dM(%d)_{n,CI}}{dKE_{CM} #Delta#Omega_{CM}} / #frac{dM(%d)_{p,CI}}{dKE_{CM} #Delta#Omega_{CM}}",sys,sys,sys));
        vector<TH1*> hists;
        ename = TString("ci_pdist_")+fSystems[isys];
        fcidist0.setmaint(maintitle);
        if (draw_dist && isys==0) fcidist0.drawadd(ename,"gridx");
        auto hist_n0c = (TH1 *) hist_n0[isys] -> Clone();
        hist_n0c -> SetMarkerStyle(29);
        hist_n0c -> SetMarkerSize(1.2);
        hist_n0c -> SetMarkerColor(ejungwoo::colori(0));
        hist_n0c -> SetLineColor(ejungwoo::colori(0));
        if (draw_dist && isys==0) ejungwoo::add(ename, 0, hist_n0c, "plhist gridx colorx", "pseudo-n");

        for (auto ipid : fIPIDAll) {
          setpar_syspid(isys,ipid);
          auto hist = histYield[isys][ipid];
          hist = ejungwoo::dndx(hist);
          if (ipid==0) hist -> SetMarkerStyle(20);
          if (ipid==1) hist -> SetMarkerStyle(21);
          if (ipid==2) hist -> SetMarkerStyle(22);
          if (ipid==3) hist -> SetMarkerStyle(33);
          if (ipid==4) hist -> SetMarkerStyle(34);
          hist -> SetMarkerColor(ejungwoo::colori(ipid+1));
          hist -> SetLineColor(ejungwoo::colori(ipid+1));
          hists.push_back(hist);
          if (draw_dist && isys==0) ejungwoo::add(ename, 0, hist, "plhist colorx", fParticleNames[ipid]);
        }

        auto graphs = draw_ci(hists, fSolidAngle, isys, hist_n0c);
        ((TGraph *) graphs.At(0)) -> SetMarkerStyle(20);
        ((TGraph *) graphs.At(0)) -> SetMarkerColor(ejungwoo::colori(isys));
        ((TGraph *) graphs.At(0)) -> SetLineColor(ejungwoo::colori(isys));
        ((TGraph *) graphs.At(1)) -> SetMarkerStyle(25);
        ((TGraph *) graphs.At(1)) -> SetMarkerColor(ejungwoo::colori(isys));
        ((TGraph *) graphs.At(1)) -> SetLineColor(ejungwoo::colori(isys));
        ((TGraph *) graphs.At(2)) -> SetMarkerStyle(22);
        ((TGraph *) graphs.At(2)) -> SetMarkerColor(ejungwoo::colori(isys));
        ((TGraph *) graphs.At(2)) -> SetLineColor(ejungwoo::colori(isys));

        ename = "CIN_CIP";
        fcidist.setmaint(maintitle);
        fcidist.drawadd(ename);
        ejungwoo::addsame(ename, graphs.At(0), "colorx pl", Form("CI-n (%d+%d)",sys,tga));
        ejungwoo::addsame(ename, graphs.At(1), "colorx pl", Form("CI-p (%d+%d)",sys,tga));

        ename = "CIN_ov_CIP";
        fciratio.setmaint(maintitle);
        if (ejungwoo::findc(ename)==nullptr) fciratio.drawadd(ename, "gridx");
        ((TGraph *) graphs.At(2)) -> SetMarkerStyle(20);
        ((TGraph *) graphs.At(2)) -> SetMarkerSize(1);
        ejungwoo::addsame(ename, graphs.At(2), "colorx pl", Form("CI-n/p (%d+%d)",sys,tga));

        graphs_rnp[isys] = (TGraphErrors *) graphs.At(2);
      }

      if (1) {
        ename = "DR_CI";
        fdbratio2.setmaint(maintitle);
        fdbratio2.drawadd(ename,"gridx");

        for (auto isyscom : {0,1,2,3}) {
          int isys1, isys2;
          Color_t color;
               if (isyscom==0) { isys1=0; isys2=1; color = kBlack; }
          else if (isyscom==1) { isys1=3; isys2=2; color = kRed; }
          else if (isyscom==2) { isys1=0; isys2=2; color = kBlue; }
          else if (isyscom==3) { isys1=3; isys2=1; color = kGreen+1; }
          if (!findisys(isys1)||!findisys(isys2))
            continue;

          auto graphr = draw_ratio(graphs_rnp[isys1], graphs_rnp[isys2],vxvar.getBinn().max);
          graphr -> SetMarkerColor(color);
          graphr -> SetLineColor(color);
          auto title1 = Form("DR(CI-n/CI-p)[(%d+%d)/(%d+%d)]", fSystems[isys1],fTargetA[isys1],fSystems[isys2],fTargetA[isys2]);
          auto ecvs = ejungwoo::add(ename, graphr, "pl colorx", title1);

          auto title2 = Form("DR(n-like/p-like)[(%d+%d)/(%d+%d)]", fSystems[isys1],fTargetA[isys1],fSystems[isys2],fTargetA[isys2]);
          if (fitydrmean[isyscom]!=nullptr) {
            fitydrmean[isyscom] -> SetLineColor(color);
            ejungwoo::addsame(ename, fitydrmean[isyscom], "colorx l addx", title2);
          }
          ecvs -> legendfb();
        }
      }
    }

    if (draw_temp)
    {
      cout << "draw_temp" << endl;

      TString ename1 = "temp";
      TString ename2 = "rheh";

      ftemp.setmaint(maintitle);
      frheh.setmaint(maintitle);

      frheh.drawadd(ename2) -> legendlt();
      ftemp.drawadd(ename1) -> legendlt();

      TH1D *histT[4] = {0};

      for (auto isys : isystemsAna) {
        auto histd = histYield[isys][1];
        auto histt = histYield[isys][2];
        auto hist3 = histYield[isys][3];
        auto hist4 = histYield[isys][4];

        auto histR = (TH1D *) histd -> Clone();
        histR -> Multiply(hist4);
        histR -> Divide(histt);
        histR -> Divide(hist3);

        auto binn = binning(histR);
        histT[isys] = (TH1D *) ejungwoo::new_h(TString("temperature")+isys,"",binn);
        histT[isys] -> SetMinimum(0);
        binn.resetb();
        while (binn.nextb()) {
          auto rvalue = histR -> GetBinContent(binn.idx);
          auto yvalue = 14.3 / (TMath::Log(1.59*rvalue));
          histT[isys] -> Fill(binn.value,yvalue);
        }

        histT[isys] -> SetLineColor(ejungwoo::colori(isys));
        histT[isys] -> SetMarkerColor(ejungwoo::colori(isys));
        histT[isys] -> SetMarkerStyle(ejungwoo::markeri(isys));

        histR -> SetLineColor(ejungwoo::colori(isys));
        histR -> SetMarkerColor(ejungwoo::colori(isys));
        histR -> SetMarkerStyle(ejungwoo::markeri(isys));

        ejungwoo::addsame(ename2,histR,"plhist gridx gridy",TString("R_{He-H}  sys")+fSystems[isys]);
        ejungwoo::addsame(ename1,histT[isys],"plhist gridx gridy",TString("temp.  sys")+fSystems[isys]);
      }

      TString ename3 = "tratio";
      ftratio.setmaint(maintitle);
      ftratio.drawadd(ename3);

      for (auto iii : {0,1})
      {
        int isys1 = 0; int isys2 = 1;
        if (iii==1) { isys1 = 3; isys2 = 2; }

        if (histT[isys1]!=nullptr||histT[isys2]!=nullptr) {
          auto histTRatio = (TH1D *) histT[isys1] -> Clone();
          histTRatio -> Divide(histT[isys2]);
          histTRatio -> SetLineColor(ejungwoo::colori(isys1));
          histTRatio -> SetMarkerColor(ejungwoo::colori(isys1));
          histTRatio -> SetMarkerStyle(ejungwoo::markeri(isys1));
          //ejungwoo::addsame(ename3,histTRatio,"plhist gridx gridy",Form("T_{%d}/T_{%d}",fSystems[isys1],fSystems[isys2]));
          ejungwoo::addsame(ename3,histTRatio,"plhist gridx gridy",Form("T%d / T%d",fSystems[isys1],fSystems[isys2]));
        }
      }
    }
  }

  //ejungwoo::drawsaveall("cvsl","pdf");
  ejungwoo::drawsaveall("cvsl","png");
}

TH1D *hist_dgraph(TString name, TTree *tree, TGraph *graph, double x1, double x2, double tta1, double tta2, int num_events)
{
  auto hist = new TH1D(name,name+";d(dedx)",200,-100,100);

  double dedx, poz, prob, eff, theta_lab;
  tree -> SetBranchAddress("dedx",&dedx);
  tree -> SetBranchAddress("p_lab",&poz);
  tree -> SetBranchAddress("prob",&prob);
  tree -> SetBranchAddress("eff",&eff);
  tree -> SetBranchAddress("theta_lab",&theta_lab);

  auto entries = tree -> GetEntries();
  for (auto i=0; i<entries; ++i)
  //for (auto i=0; i<100; ++i)
  {
    tree -> GetEntry(i);

    if (poz>x1 && poz<x2)
    {
      if (theta_lab*57.295780>tta1 && theta_lab*57.295780<tta2) {
        double scale = prob/eff/num_events;
        if (!docorr)
          scale = 1;
        auto ddedx = dedx - graph -> Eval(poz);
        hist -> Fill(ddedx,scale);
      }
    }
  }

  return hist;
}
