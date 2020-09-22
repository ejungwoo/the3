#include "/Users/ejungwoo/config/ejungwoo.h"
#include "init_variables.h"
#include "functions.h"

const int fIPIDs[] = {0,1,2,3,4};

void draw_j()
{
  ejungwoo::gcvspos(1300);
  ejungwoo::gstat(0);
  ejungwoo::gsave(1);
  ejungwoo::gsetvmark(0);
  ejungwoo::gshortprint();

  int draw_pydist = 1;
  int draw_pid = 1;
  int draw_like = 1;
  int draw_npratio = 0;
  int draw_dbratio = 0;
  int draw_r21 = 1;
  int draw_srn0 = 0;
  int draw_psdon = 0;
  int draw_cici = 1;
  int draw_dist = 1;

  auto remove_draw = [&draw_pydist, &draw_pid, &draw_like, &draw_npratio, &draw_dbratio, &draw_r21, &draw_srn0, &draw_psdon, &draw_cici, &draw_dist]() {
    draw_pydist = 0; draw_pid = 0; draw_like = 0; draw_npratio = 0; draw_dbratio = 0; draw_r21 = 0; draw_srn0 = 0; draw_psdon = 0; draw_cici = 0; draw_dist = 0;
  };

  double num_tracks_per_event_cut = 0.01;

  init();

  const char *particleNames[] = {"p","d","t","he3","he4","he6","1"};

  int ipid_nl[3][2] = {{4,3},{1,0},{2,1}}; int num_nl = 3;
  int ipid_pl[3][2] = {{0,6},{4,2},{3,1}}; int num_pl = 3;
  int ipid_dl[3][2] = {{1,6},{3,0},{4,1}}; int num_dl = 3;
  int ipid_tl[2][2] = {{2,6},{4,0}}; int num_tl = 2;

  int isystems[] = {0,1};

  int yy1 = -0;
  int yy2 = 10;
  bool useRapidityBinning = 0;
  const char *yBeam = "by_cm0";
  //const char *yBeam = "by_cm/2";

  //remove_draw(); draw_pydist = 1; draw_pid = 1;

  //for (auto aversion : {"jright"})
  for (auto aversion : {"jrighti132"})
  //for (auto aversion : {"jkright"})
  {
    auto condition_array = setversion(aversion);

    auto vpydist = variable("pydist", Form("pt_cm/$$(a):fy_cm/(%s)",yBeam), "$$(cut_pe)", titles(" ", "y_{0}", "p_{T}/A (MeV/c)"), binning(100,-1.,2.), binning(100,0,1000));
    auto vpid2 = variable("pid", "dedx:p_lab", "$$(cut0)*(dedx<600)", titles(" ","p_{Lab}/Z (MeV/c^{2})","dE/dx"), binning(400,0,2500),binning(400,0,600)); 
    //auto vpid2 = variable("pid", "dedx:pt_cm/$$(a)", "$$(cut0)*(dedx<600)", titles(" ","p_{T}/A (MeV/c)","dE/dx"), binning(400,0,1000),binning(400,0,600)); 

    binning binn0;
    const char *xTitle;
    variable vxvar = variable("");

    if (useRapidityBinning) {
      binn0 = binning(8,0.5,1.5);
      xTitle = "y_{0} = 2y_{z}/y_{beam Lab.}";
      vxvar = variable("y0", Form("fy_cm/(%s)",yBeam), "$$(cut0)", titles("",xTitle,""), binn0);
      yy1 = -10;
      yy2 = 20;
    } else {
      //binn0 = binning(10,0,1000);
      binn0 = binning(8,0,400);
      //binn0 = binning(5,0,250);
      xTitle = "p_{T}/A (MeV/c)";
      vxvar = variable("ptoa_cm", "pt_cm/$$(a)", "$$(cut0)", titles("","p_{T}/A (MeV/c)",""), binn0);
    }

    auto ccc = 1.;
    setpar("rap_cut",Form("$$(foby_cm)>%.1f && $$(foby_cm)<%.1f",.1*yy1,.1*yy2));
    TString maintitle = Form("y_{0} = -%.1f ~ %.1f",.1*yy1,.1*yy2);
    TString vout2 = fVOutShort+"_yy_"+yy1+"_"+yy2;

    auto fnpdtlike = variable("fnpdtlike", "", "", titles(" ",xTitle,"Y(particle-0)/Y(particle-1)"), binn0, binning(100,0,5));
    auto fnpratio = variable("fnpratio", "", "", titles(" ",xTitle,"n-like/p-like"), binn0, binning(100,0,20));
    auto fdbratio = variable("fdbratio", "", "", titles(" ",xTitle,"DR(n-like/p-like)$$(systga21)"), binn0, binning(100,0,2.5));
    auto fdbratio2 = variable("fdbratio2", "", "", titles(" ",xTitle,"DR(#frac{CI-n}{CI-p})"), binn0, binning(100,0,2.0));
    auto fptoar21 = variable("fptoar21", "", "", titles(" ",xTitle,"R21"), binn0, binning(100,0.5,2.));
    auto fptoasr = variable("fptoasr", "", "", titles(" ",xTitle,"SR(t/he3)"), binn0, binning(100,0,7));
    auto fptoan0 = variable("fptoan0", "", "", titles(" ",xTitle,"pseudo neutron  #frac{d^{2}M}{dyd(p_{T}/A)}"), binn0, binning(100,0,0.2));
    auto fcidist0 = variable("cidist", "", "", titles(" ",xTitle,"#frac{d^{2}M}{dyd(p_{T}/A)}"), binn0, binning(100,0,ccc*0.06));
    auto fcidist = variable("cidist", "", "", titles(" ",xTitle,"#frac{d^{2}M}{dyd(p_{T}/A)}"), binn0, binning(100,0,ccc*0.15));
    auto fciratio = variable("fciratio", "", "", titles(" ",xTitle,"#frac{CI-n}{CI-p}"), binn0, binning(100,0,2.0));

    if (draw_pydist) {
      //for (auto isys : {0,1,2,3}) {
      for (auto isys : {0,1}) {
        TString ename = TString("sysCM")+fSystems[isys];
        for (int ipid : fIPIDs) {
          setpar_syspid(isys,ipid);
          vpydist.setmaint(Form("sys%d %s",fSystems[isys],fParticleNames[ipid].Data()));
          auto hist = vpydist.draw(fTreePID[isys][ipid]);
          hist -> SetMinimum(0);
          ejungwoo::addnext(ename,hist);
        }
      }
    }
    if (draw_pydist==2) break;


    for (auto ipid : fIPIDs)
      fSDHLCut[ipid] = 2.0;
    vout2 = vout2 + "_sd2";
    ejungwoo::gversionout(vout2);

    if (draw_pid) {
      bool addanddraw = 1;
      for (auto isys : {0,1,2,3}) {
        cout << isys << endl;
        TH2D *hist2 = nullptr;
        vpid2.setmaint(TString("sys") + fSystems[isys] + " " + maintitle);
        for (int ipid : fIPIDs) {
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
    for (auto isys : {0,1,2,3}) {
      for (auto ipid : fIPIDs) {
        setpar_syspid(isys,ipid);
        histYield[isys][ipid] = vxvar.draw(fTreePID[isys][ipid]);
      }
    }

    TH1 *histYieldLike[4][4][3]; // isys, inpdt, icomb
    TGraphErrors *graphdrmean[4] = {0};
    TH1D *histdrmean[4] = {0};
    TF1 *fitydrmean[4] = {0};
    while (draw_r21 || draw_like || draw_npratio || draw_dbratio)
    {
      for (auto isys : {0,1,2,3})
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

            if (ipid2!=6) histl -> Divide((TH1 *) histYield[isys][ipid2]);
            if (draw_like && isys==0) ejungwoo::addsame(ename,histl,"plhist colorx");

            histYieldLike[isys][inpdt][icomb] = histl;
          }
        }
      }
      if (draw_like==2) break;

      if (!(draw_like || draw_npratio || draw_dbratio)) break;


      TH1D *histnp[4][9];
      for (auto isys : {0,1,2,3})
      {
        TString ename = Form("sys%d_n-like/p-like",fSystems[isys]);
        fnpratio.setmaint(maintitle);

        if (draw_npratio) ejungwoo::add(ename,fnpratio.draw());
        int idx = 0;
        for (auto inl=0; inl<3; ++inl)
          for (auto ipl=0; ipl<3; ++ipl) {
            TString npname = Form("sys%d_nplike%d%d",fSystems[isys],inl,ipl);
            histnp[isys][idx] = (TH1D *) histYieldLike[isys][0][inl] -> Clone(npname);
            histnp[isys][idx] -> SetTitle(Form("%d. (%s/%s) / (%s/%s)",idx+1,particleNames[ipid_nl[inl][0]],particleNames[ipid_nl[inl][1]],particleNames[ipid_pl[ipl][0]],particleNames[ipid_pl[ipl][1]]));
            histnp[isys][idx] -> Divide(histYieldLike[isys][1][ipl]);
                 if (inl==0) histnp[isys][idx] -> SetMarkerStyle(5);
            else if (inl==1) histnp[isys][idx] -> SetMarkerStyle(27);
            else if (inl==2) histnp[isys][idx] -> SetMarkerStyle(42);
            if (draw_npratio) ejungwoo::addsame(ename,histnp[0][idx],"hist");
            idx++;
          }
      }
      if (draw_npratio==2) break;


      for (auto isyscom : {0,1,2,3})
      {
        int isys1, isys2;
        //TString drawoption = "addx";
        TString drawoption = "";
        //if (isyscom==0) { isys1=0; isys2=1; drawoption=""; }
        if (isyscom==0) { isys1=0; isys2=1; drawoption=""; }
        else if (isyscom==1) { isys1=3; isys2=2; }
        else if (isyscom==2) { isys1=0; isys2=2; }
        else if (isyscom==3) { isys1=3; isys2=1; }

        TString ename = TString("DR")+isyscom;
        fdbratio.setmaint(maintitle);
        //fdbratio.setmaint(Form("(%d+%d) / (%d+%d)", fSystems[isys1],fTargetA[isys1],fSystems[isys2],fTargetA[isys2]));
        ejungwoo::setpar("systga21",Form("_{#frac{%d+%d}{%d+%d}}", fSystems[isys1],fTargetA[isys1],fSystems[isys2],fTargetA[isys2]));
        if (draw_dbratio) {
          auto ecvs = fdbratio.drawadd(ename,"gridxgridy");
          ecvs -> legendxx();
        }
        vector<TH1 *> histdrarray;
        for (auto idx=0; idx<9; ++idx) {
          TString drname = Form("DR%d",idx);
          auto histdr = (TH1D *) histnp[isys1][idx] -> Clone(drname);
          histdr -> Divide(histnp[isys2][idx]);
          if (draw_dbratio) ejungwoo::addsame(ename,histdr,drawoption+"plhist");
          if (draw_dbratio && idx==0) ejungwoo::addsame(ename,histdr,drawoption+"plhist");
          histdrarray.push_back(histdr);
        }
        if (draw_dbratio==2) break;


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
          //graphdrmean[isyscom] -> SetPointError(graphdrmean[isyscom]->GetN()-1,0,yvalue*.1); //XXX
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

    if (draw_r21) {
      for (auto isyscom : {0,1,2,3})
      {
        int isys1, isys2;
        if (isyscom==0) { isys1=0; isys2=1; }
        else if (isyscom==1) { isys1=3; isys2=2; }
        else if (isyscom==2) { isys1=0; isys2=2; }
        else if (isyscom==3) { isys1=3; isys2=1; }
        //132, 108, 112, 124

        TString systga1 = Form("%d+%d",fSystems[isys1],fTargetA[isys1]);
        TString systga2 = Form("%d+%d",fSystems[isys2],fTargetA[isys2]);

        TString ename = Form("r21_%d_%d",fSystems[isys1],fSystems[isys2]);
        fptoar21.setmaint(maintitle);
        fptoar21.setyt(Form("R_{21}( #frac{%s}{%s} )",systga1.Data(),systga2.Data()));
        fptoar21.drawadd(ename,"gridxgridy");
        TGraph *gfree[6] = {0,};
        for (int ipid : fIPIDs) {
          auto graphs = draw_is(ipid, histYield[isys1][ipid], histYield[isys2][ipid], num_tracks_per_event_cut);
          gfree[ipid] = ((TGraph *) graphs.At(0));
          gfree[ipid] -> SetMarkerSize(1);
          ejungwoo::add(ename, 0, gfree[ipid], "pl", fParticleNames[ipid]);
        }

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

            auto graphs = draw_is(0, histYieldLike[isys1][inpdt][icomb], histYieldLike[isys2][inpdt][icomb], num_tracks_per_event_cut);
            auto graph0 = (TGraphErrors *) graphs.At(0);

            TString graphTitle = Form("%s/%s",particleNames[ipid1],particleNames[ipid2]);
            int mstyle, mcolor;
                 if (inpdt==0) { mstyle = 29; graphTitle = TString("n-like ")+graphTitle; mcolor = (icomb==0?kBlack:(icomb==1?kGray+2:kGray)); }
            else if (inpdt==1) { mstyle = 20; graphTitle = TString("p-like ")+graphTitle; mcolor = (icomb==0?kRed+1:(icomb==1?kOrange-3:kPink+9)); }
            else if (inpdt==2) { mstyle = 21; graphTitle = TString("d-like ")+graphTitle; mcolor = (icomb==0?kBlue+1:(icomb==1?kAzure-3:kAzure+8)); }
            else if (inpdt==3) { mstyle = 22; graphTitle = TString("t-like ")+graphTitle; mcolor = (icomb==0?kGreen+3:(icomb==1?kSpring+5:kTeal+2)); }
            graph0 -> SetMarkerSize(1);

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


    TH1 *hist_n0[4];
    if (draw_cici || draw_psdon)
    {
      for (auto isys : {0,1,2,3}) {
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
      //ejungwoo::add(ename, 0, new_h(vxvar.getHistName()+"_cin_allframe", ttlCIM, binnCIx, binnCIy2), "addx rangex logy ");
      //ejungwoo::add(ename, 1, new_h(vxvar.getHistName()+"_cip_allframe", ttlCIM, binnCIx, binnCIy2), "addx rangex logy ");
      //ejungwoo::add(ename, 2, new_h(vxvar.getHistName()+"_cir_allframe", ttlCIR, binnCIx, binnCIr), "addx rangex ");

      TString ename;

      TGraphErrors *graphs_rnp[4];
      for (auto isys : {0,1,2,3}) {
        auto sys = fSystems[isys];
        auto tga = fTargetA[isys];
        vxvar.setmaint("($$(sys)) particle multiplicity");
        titles ttlCIMN("($$(sys)) CI neutrons","KE_{CM} (MeV)",Form("#frac{dM(%d)_{n,CI}}{dKE_{CM} #Delta#Omega_{CM}}",sys));
        titles ttlCIMP("($$(sys)) CI protons","KE_{CM} (MeV)",Form("#frac{dM(%d)_{p,CI}}{dKE_{CM} #Delta#Omega_{CM}}",sys));
        titles ttlCIRs("($$(sys)) CI n/p","KE_{CM} (MeV)",Form("R(%d) = #frac{dM(%d)_{n,CI}}{dKE_{CM} #Delta#Omega_{CM}} / #frac{dM(%d)_{p,CI}}{dKE_{CM} #Delta#Omega_{CM}}",sys,sys,sys));
        vector<TH1*> hists;
        ename = TString("ci_pdist_")+fSystems[isys];
        fcidist0.setmaint(maintitle);
        //if (draw_dist && isys==0) fcidist0.drawadd(ename,"gridx");
        auto hist_n0c = (TH1 *) hist_n0[isys] -> Clone();
        hist_n0c -> SetMarkerStyle(29);
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
        //if (ejungwoo::findc(ename)==nullptr) {
        //  auto ecvs = fcidist.drawadd(ename,"gridx");
        //  ecvs -> legendlt();
        //}
        ejungwoo::addsame(ename, graphs.At(0), "colorx pl", Form("CI-n (%d+%d)",sys,tga));
        ejungwoo::addsame(ename, graphs.At(1), "colorx pl", Form("CI-p (%d+%d)",sys,tga));

        ename = "CIN_ov_CIP";
        //fciratio.setmaint("CI-n / CI-p");
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
  }

  ejungwoo::drawsaveall("cvsl","pdf");
}
