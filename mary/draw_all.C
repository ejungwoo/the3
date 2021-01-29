#include "/Users/ejungwoo/config/ejungwoo.h"
#include "init_variables.h"

void draw_all()
{
  ejungwoo::gcvspos(1000);
  //ejungwoo::gcvspos(1300);
  ejungwoo::gstat(0);
  ejungwoo::gsave(1);
  ejungwoo::gsetvmark(0);
  //ejungwoo::gshortprint();
  //ejungwoo::gdummytp();

  int draw_pydist = 1;
  int draw_pndist = 2;
  int draw_pid = 0;
  int draw_like = 0;
  int draw_npratio = 0;
  int draw_dbratio = 0;
  int draw_r21 = 0;
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

  double num_tracks_per_event_cut = 0.01;

  init();

  //setpar("cut1","$$(cut0)");
  setpar("cut1","$$(cut_pek)");
  setpar("cut2","$$(cut_noy)");

  //setpar("cut_pek","$$(prob)/effk/$$(num_events)*($$(mult_cut))*($$(rap_cut))*($$(pes_cut))*($$(poz_cut))*($$(angle_cut))");
  setpar("cut_pek","$$(prob)/$$(eff)/$$(num_events)*($$(mult_cut))*($$(rap_cut))*($$(pes_cut))*($$(poz_cut))*($$(angle_cut))");

  const char *particleNames[] = {"p","d","t","he3","he4","he6","1"};

  int ipid_nl[3][2] = {{4,3},{1,0},{2,1}}; int num_nl = 3;
  int ipid_pl[3][2] = {{0,6},{4,2},{3,1}}; int num_pl = 3;
  int ipid_dl[3][2] = {{1,6},{3,0},{4,1}}; int num_dl = 3;
  int ipid_tl[2][2] = {{2,6},{4,0}}; int num_tl = 2;

  int isystemsAll[] = {0,1,2,3};
  int isystemsAna[] = {0,1,2,3};
  auto findisys = [isystemsAna](int isys) {
    for (auto isys0 : isystemsAna)
      if (isys0==isys)
        return true;
     return false;
  };

  const int ipidsAna[] = {0,1,2,3,4};
  //const int ipidsAna[] = {0,1,2};
  auto findipid = [ipidsAna](int ipid) {
    if (ipid==6)
      return true;
    for (auto ipid0 : ipidsAna)
      if (ipid0==ipid)
        return true;
     return false;
  };

  int yy1 = 0; int yy2 = 4;
  //int yy1 = -4; int yy2 = 0;
  int pt1 = 100; int pt2 = 200;
  remove_draw(); draw_pydist = 1; draw_r21nz = 2;
  auto useSingleXBinning = true;

  //int yy1 = -2; int yy2 = 2;
  //int yy1 = 0; int yy2 = 5;
  //int yy1 = 5; int yy2 = 10;
  //int pt1 = -1; int pt2 = -1;
  //auto useSingleXBinning = false;

  if (pt1<0) {
    pt1 = 0;
    pt2 = 400;
  }

  //int pt1 = -1; int pt2 = -1;
  //int yy1 = -30; int yy2 = 30;
  //int yy1 = 0; int yy2 = 10;
  //int yy1 = 0; int yy2 = 5;
  //int yy1 = 5; int yy2 = 10;
  bool useRapidityBinning = 0;
  //const char *yBeam = "by_cm0";
  const char *yBeam = "by_cm/2";

  //remove_draw(); draw_pydist = 1; draw_pid = 1;

  //for (auto aversion : {"right_55"})
  for (auto aversion : {"right_50"})
  {
    auto condition_array = setversion(aversion);

    auto vpydist = variable("pydist", Form("pt_cm/$$(a):fy_cm/(%s)",yBeam), "$$(cut2)", titles(" ", "y_{0}", "p_{T}/A (MeV/c)"), binning(100,-1.,2.), binning(100,0,1000));
    //auto vpndist = variable("pndist", "$$(nr):pt_cm/$$(a)", "($$(nr)!=0)*($$(cut_x))", titles(" ", "p_{T}/A (MeV/c)", "nc"), binning(100,0,400), binning(100,0,100));
    //auto vpndist = variable("pndist", "$$(nl):pt_cm/$$(a)", "($$(nl)!=0)*($$(cut_x))", titles(" ", "p_{T}/A (MeV/c)", "nc"), binning(100,0,400), binning(100,0,100));
    auto vpndist = variable("pndist", "($$(nr)+$$(nl)):pt_cm/$$(a)", "$$(cut_x)", titles(" ", "p_{T}/A (MeV/c)", "nc"), binning(100,0,400), binning(100,0,100));
    auto vpid2 = variable("pid", "dedx:p_lab", "$$(cut1)*(dedx<600)", titles(" ","p_{Lab}/Z (MeV/c^{2})","dE/dx"), binning(400,0,2500),binning(400,0,600)); 
    //auto vpid2 = variable("pid", "dedx:pt_cm/$$(a)", "$$(cut0)*(dedx<600)", titles(" ","p_{T}/A (MeV/c)","dE/dx"), binning(400,0,1000),binning(400,0,600)); 

    binning binn0;
    const char *xTitle;
    const char *xVar;
    variable vxvar = variable("");

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
    if (pt1>=0)
      setpar("rap_cut",Form("$$(foby_cm)>%.1f && $$(foby_cm)<%.1f && $$(ptoa_cm)>%d &&  $$(ptoa_cm)<%d",.1*yy1,.1*yy2,pt1,pt2));
    TString maintitle = Form("Mult= %d,  y_{0} = %.1f ~ %.1f",fMultLLCut[1], .1*yy1,.1*yy2);
    if (pt1>=0)
      maintitle = Form("Mult= %d,  y_{0} = %.1f ~ %.1f, p_{T}/A = %d ~ %d",fMultLLCut[1], .1*yy1,.1*yy2, pt1, pt2);
    if (useSingleXBinning)
      fVOutShort = fVOutShort + "_single";
    TString vout2 = fVOutShort+"_yy_"+yy1+"_"+yy2;
    if (pt1>=0)
      vout2 = fVOutShort+Form("_pt%03d%03d_yy%02d%02d",pt1,pt2,yy1,yy2);

    cout << vout2 << endl;

    auto fnpdtlike = variable("fnpdtlike", "", "", titles(" ",xTitle,"Y(particle-0)/Y(particle-1)"), binn0, binning(100,0,5));
    auto fnpratio = variable("fnpratio", "", "", titles(" ",xTitle,"n-like/p-like"), binn0, binning(100,0,20));
    auto fdbratio = variable("fdbratio", "", "", titles(" ",xTitle,"DR(n-like/p-like)$$(systga21)"), binn0, binning(100,0,2.5));
    auto fdbratio2 = variable("fdbratio2", "", "", titles(" ",xTitle,"DR(#frac{CI-n}{CI-p})"), binn0, binning(100,0,2.0));
    auto fptoar21 = variable("fptoar21", "", "", titles(" ",xTitle,"R21"), binn0, binning(100,0.5,2.));
    //auto fnr21 = variable("fnr21", "", "", titles(" ","n","R21($$(sys1)/$$(sys2))"), binning(4,-1,3), binning(100,0.6,4));
    //auto fzr21 = variable("fzr21", "", "", titles(" ","z","R21($$(sys1)/$$(sys2))"), binning(4,0,3), binning(100,0.6,4)); bool setfznzlog = true;
    auto fnr21 = variable("fnr21", "", "", titles(" ","n","R21($$(sys1)/$$(sys2))"), binning(4,-1,3), binning(100,0.5,2));
    auto fzr21 = variable("fzr21", "", "", titles(" ","z","R21($$(sys1)/$$(sys2))"), binning(4,-1,3), binning(100,0.5,2)); bool setfznzlog = false;
    auto fptoasr = variable("fptoasr", "", "", titles(" ",xTitle,"SR(t/he3)"), binn0, binning(100,0,7));
    auto fptoan0 = variable("fptoan0", "", "", titles(" ",xTitle,"pseudo neutron  #frac{d^{2}M}{dyd(p_{T}/A)}"), binn0, binning(100,0,0.2));
    auto fcidist0 = variable("cidist", "", "", titles(" ",xTitle,"#frac{d^{2}M}{dyd(p_{T}/A)}"), binn0, binning(100,0,ccc*0.055));
    auto fcidist = variable("cidist", "", "", titles(" ",xTitle,"#frac{d^{2}M}{dyd(p_{T}/A)}"), binn0, binning(100,0,ccc*0.15));
    auto fciratio = variable("fciratio", "", "", titles(" ",xTitle,"#frac{CI-n}{CI-p}"), binn0, binning(100,0,2.0));
    auto ftemp = variable("temp", "", "", titles(" ",xTitle,"T = #frac{14.3}{log[1.59*R_{He-H}]}"), binn0, binning(100,0,60));
    auto yrheh = Form("R_{He-H} = #scale[0.8]{#frac{dM_{d}/d#Omegad(%s) #times dM_{he4}/d#Omegad(%s)}{dM_{t}/d#Omegad(%s) #times dM_{he3}/d#Omegad(%s)}}",xVar,xVar,xVar,xVar);
    auto frheh = variable("rheh", "", "", titles(" ",xTitle,yrheh), binn0, binning(100,0,5.));
    auto ftratio = variable("tratio", "", "", titles(" ",xTitle,"T_{2} / T_{1}"), binn0, binning(100,0,2.5));

    if (draw_pydist) {
      cout << "draw_pydist" << endl;
      //for (auto isys : isystemsAna)
      for (auto isys : {0})
      {
        TString ename = TString("sysCM")+fSystems[isys];
        for (int ipid : ipidsAna)
        //for (int ipid : {0})
        {
          setpar_syspid(isys,ipid);
          vpydist.setmaint(Form("sys%d %s",fSystems[isys],fParticleNames[ipid].Data()));
          auto hist = vpydist.draw(fTreePID[isys][ipid]);
          hist -> SetMinimum(0);
          ejungwoo::addnext(ename,hist,"gridx gridy");
          {
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
        }
      }
    }

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


    for (auto ipid : ipidsAna)
      fSDHLCut[ipid] = 2.0;
    vout2 = vout2 + "_sd2";
    if (useRapidityBinning)
      vout2 = vout2 + "_by";
    ejungwoo::gversionout(vout2);
    if (draw_pydist==2) break;
    if (draw_pndist==2) break;

    if (draw_pid) {
      cout << "draw_pid" << endl;
      bool addanddraw = 1;
      for (auto isys : isystemsAna) {
        TH2D *hist2 = nullptr;
        vpid2.setmaint(TString("sys") + fSystems[isys] + " " + maintitle);
        for (int ipid : ipidsAna) {
          setpar_syspid(isys,ipid);
          auto histpid = (TH2D *) vpid2.draw(fTreePID[isys][ipid]);
          if (addanddraw)  {
            if (ipid==0) hist2 = histpid;
            else hist2 -> Add(vpid2.draw(fTreePID[isys][ipid]));
          }
          else
            ejungwoo::addnext(TString("pid")+isys,histpid);
        }
        if (addanddraw) {
          hist2 -> SetMaximum(0.02);
          hist2 -> SetMinimum(0.00005);
          ejungwoo::addnext(TString("pid"),hist2,"logz");
        }
      }
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
    while (draw_r21 || draw_r21nz || draw_like || draw_npratio || draw_dbratio)
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

    if (draw_r21 || draw_r21nz)
    {
      cout << "draw_r21" << endl;
      for (auto isyscom : {0,1,2,3})
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

          setpar("sys1",fSystems[isys1]);
          setpar("sys2",fSystems[isys2]);

          if (setfznzlog) {
            fnr21.setmaint(maintitle); ejungwoo::add(ename2,0,fnr21.draw(),"gridx gridy logy") -> legendlt();
            fzr21.setmaint(maintitle); ejungwoo::add(ename2,1,fzr21.draw(),"gridx gridy logy");
          }
          else {
            fnr21.setmaint(maintitle); ejungwoo::add(ename2,0,fnr21.draw(),"gridx gridy") -> legendlt();
            fzr21.setmaint(maintitle); ejungwoo::add(ename2,1,fzr21.draw(),"gridx gridy");
          }

          auto graph_hz = new TGraph();
          auto graph_hn = new TGraph();
          graph_hz -> SetMarkerSize(1.5);
          graph_hn -> SetMarkerSize(1.5);
          graph_hz -> SetMarkerStyle(20);
          graph_hn -> SetMarkerStyle(20);
          graph_hz -> SetMarkerColor(kBlack);
          graph_hn -> SetMarkerColor(kBlack);
          for (int ipid : {0,1,2}) {
            Double_t xdummy, r21Value;
            gfree[ipid] -> GetPoint(0,xdummy,r21Value);
            auto znumber = fNumProtons[ipid];
            auto nnumber = fNumNeutrons[ipid];
            graph_hz -> SetPoint(graph_hz->GetN(),znumber,r21Value);
            graph_hn -> SetPoint(graph_hn->GetN(),nnumber,r21Value);
          }
          ejungwoo::add(ename2, 0, graph_hn, "colorx p", "H");
          ejungwoo::add(ename2, 1, graph_hz, "colorx p", "H");

          auto graph_hez = new TGraph();
          auto graph_hen = new TGraph();
          graph_hez -> SetMarkerSize(1.5);
          graph_hen -> SetMarkerSize(1.5);
          graph_hez -> SetMarkerStyle(25);
          graph_hen -> SetMarkerStyle(25);
          graph_hez -> SetMarkerColor(kRed);
          graph_hen -> SetMarkerColor(kRed);

          for (int ipid : {3,4}) {
            Double_t xdummy, r21Value;
            gfree[ipid] -> GetPoint(0,xdummy,r21Value);
            auto znumber = fNumProtons[ipid];
            auto nnumber = fNumNeutrons[ipid];
            graph_hez -> SetPoint(graph_hez->GetN(),znumber,r21Value);
            graph_hen -> SetPoint(graph_hen->GetN(),nnumber,r21Value);
          }
          ejungwoo::add(ename2, 1, graph_hez, "colorx p", "He");
          ejungwoo::add(ename2, 0, graph_hen, "colorx p", "He");
        }
        if (draw_r21nz==2) continue;

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
      ftemp.drawadd(ename1);

      TH1D *histT[4] = {0};

      for (auto isys : isystemsAna) {
        auto histd = histYield[isys][1];
        auto histt = histYield[isys][2];
        auto hist3 = histYield[isys][3];
        auto hist4 = histYield[isys][4];

        auto histR = (TH1D *) histd -> Clone();
        histR -> Multiply(hist4);
        histR -> Divide(histt);
        histR -> Divide(hist4);

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

  ejungwoo::drawsaveall("cvsl","pdf");
}
