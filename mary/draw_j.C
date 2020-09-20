#include "/Users/ejungwoo/config/ejungwoo.h"
#include "init_variables.h"
#include "functions.h"

const int fIPIDs[] = {0,1,2,3,4};

void draw_j()
{
  ejungwoo::gcvspos(900);
  ejungwoo::gstat(0);
  ejungwoo::gsave(0);
  ejungwoo::gsetvmark(0);
  //ejungwoo::gshortprint();

  bool draw_pydist = 0;
  bool draw_like = 1;
  bool draw_npratio = 0;
  bool draw_dbratio = 1;
  bool draw_r21 = 1;
  bool draw_srn0 = 0;
  bool draw_psdon = 1;
  bool draw_cici = 1;
  bool draw_dist = 0;

  if (draw_cici) draw_psdon = 1;

  double num_tracks_per_event_cut = 0.01;

  init();

  const char *particleNames[] = {"p","d","t","he3","he4","he6","1"};

  //auto vptoa = variable("ptoa_cm", "pt_cm/$$(a)", "$$(cut0)", titles("","p_{T}/A (MeV/c)",""), binning(15,100,400));
  //auto vptoa = variable("ptoa_cm", "pt_cm/$$(a)", "$$(cut0)", titles("","p_{T}/A (MeV/c)",""), binning(10,100,400));
  //auto vptoa = variable("ptoa_cm", "pt_cm/$$(a)", "$$(cut0)", titles("","p_{T}/A (MeV/c)",""), binning(10,100,250));
  auto vptoa = variable("ptoa_cm", "pt_cm/$$(a)", "$$(cut0)", titles("","p_{T}/A (MeV/c)",""), binning(5,0,250));

  auto fnplike = variable("fnplike", "", "", titles(" ","p_{T}/A (MeV/c)","p0/p1"), 100,0,250,100,0,4);
  auto fnpratio = variable("fnpratio", "", "", titles(" ","p_{T}/A (MeV/c)","n-like/p-like"), 100,0,250,100,0,20);
  auto fdbratio = variable("fdbratio", "", "", titles(" ","p_{T}/A (MeV/c)","DR(n-like/p-like)$$(systga21)"), 100,0,250,100,0,2.5);
  auto fdbratio2 = variable("fdbratio2", "", "", titles(" ","p_{T}/A (MeV/c)","DR(CI-n/CI-p)"), 100,0,250,100,0,2.0);
  auto fptoar21 = variable("fptoar21", "", "", titles(" ","p_{T}/A (MeV/c)","R21"), 100,0,250,100,0.5,2);
  auto fptoasr = variable("fptoasr", "", "", titles(" ","p_{T}/A (MeV/c)","SR(t/he3)"), 100,0,250,100,0,7);
  //auto fptoan0 = variable("fptoan0", "", "", titles(" ","p_{T}/A (MeV/c)","#frac{d^{2}M}{dyd(p_{T}/A)}"), 100,0,250,100,0,0.06);
  auto fptoan0 = variable("fptoan0", "", "", titles(" ","p_{T}/A (MeV/c)","pseudo neutron  #frac{d^{2}M}{dyd(p_{T}/A)}"), 100,0,250,100,0,0.2);
  //auto fcidist = variable("cidist", "", "", titles(" ","p_{T}/A (MeV/c)","#frac{d^{2}M}{dyd(p_{T}/A)}"), 100,0,250,100,0.0003,0.0025);
  auto fcidist = variable("cidist", "", "", titles(" ","p_{T}/A (MeV/c)","#frac{d^{2}M}{dyd(p_{T}/A)}"), 100,0,250,100,0.0003,0.006);
  auto fciratio = variable("cidist", "", "", titles(" ","p_{T}/A (MeV/c)","CI-n/CI-p"), 100,0,250,100,0,2.0);

  int ipidnl[3][2] = {{4,3},{1,0},{2,1}};
  int ipidpl[3][2] = {{0,6},{4,2},{3,1}};
  //int ipiddl[3][2] = {{1,6},{3,0},{4,1}};
  //int ipidpl[3][2] = {{1,6},{3,0},{4,1}};

  int isystems[] = {0,1};

  for (auto aversion : {"jright"})
  {
    auto condition_array = setversion(aversion);
    settrees();

    //setpar("rap_cut","$$(foby_cm)>-.5&&$$(foby_cm)<.5"); TString maintitle = "y_{0} = -0.5 ~ 0.5"; TString vout2 = "jright_m50_y0";
    //setpar("rap_cut","$$(foby_cm)>-.5&&$$(foby_cm)<0."); TString maintitle = "y_{0} = -0.5 ~ 0."; TString vout2 = "jright_m50_y1";
    //setpar("rap_cut","$$(foby_cm)>0&&$$(foby_cm)<.5"); TString maintitle = "y_{0} = 0 ~ 0.5"; TString vout2 = "jright_m50_y02;
    setpar("rap_cut","$$(foby_cm)>0.5&&$$(foby_cm)<1."); TString maintitle = "y_{0} = 0.5 ~ 1."; TString vout2 = "jright_m50_y3";

    ejungwoo::gversionout(vout2);

    if (draw_pydist)
      for (auto isys : {0,1,2,3}) {
        //TString ename = TString("sysLab")+fSystems[isys];
        //for (int ipid : fIPIDs) {
        //  setpar_syspid(isys,ipid);
        //  fvar_ypt_lab.drawaddnext(fTreePID[isys][ipid],ename);
        //}
        TString ename = TString("sysCM")+fSystems[isys];
        for (int ipid : fIPIDs) {
          setpar_syspid(isys,ipid);
          fvar_ypt_cm.drawaddnext(fTreePID[isys][ipid],ename);
        }
      }

    TGraphErrors *graphdrmean[4];
    TH1D *histdrmean[4];
    TF1 *fitydrmean[4];
    if (draw_like || draw_npratio || draw_dbratio)
    {
      TH1 *histnlike[4][3]; // isys, inpl
      TH1 *histplike[4][3]; // isys, inpl
      TH1 *histdlike[4][3]; // isys, inpl
      TH1 *histndist[4][3][3]; // isys, inpl, idx
      TH1 *histpdist[4][3][3]; // isys, inpl, idx
      TH1 *histddist[4][3][3]; // isys, inpl, idx

      for (auto isys : {0,1,2,3})
      {
        cout << isys << endl;
        TString ename = Form("sys%d_n-like",fSystems[isys]);
        fnplike.setmaint(maintitle);
        if (draw_like) ejungwoo::add(ename,fnplike.draw());
        for (auto inl=0; inl<3; ++inl) {
          const char *pname[2] = {particleNames[ipidnl[inl][0]],particleNames[ipidnl[inl][1]]};
          TString ename3 = Form("sys%d_n-like_%s/%s",fSystems[isys],pname[0],pname[1]);
          for (auto idx : {0,1}) {
            auto ipid = ipidnl[inl][idx];
            setpar_syspid(isys,ipid);
            histndist[isys][inl][idx] = vptoa.draw(fTreePID[isys][ipid]);
            histndist[isys][inl][idx] -> SetTitle(pname[idx]);
            //ejungwoo::addsame(ename3,histndist[isys][inl][idx]);
            if (idx==0)
              histnlike[isys][inl] = vptoa.draw(fTreePID[isys][ipid]);
          }
          auto ecvs = ejungwoo::findc(ename3);
          if (ecvs!=nullptr)
            ecvs-> setTitle(Form("sys%d n-like %s/%s",fSystems[isys],pname[0],pname[1]));

          histnlike[isys][inl] -> Divide(histndist[isys][inl][1]);
          histnlike[isys][inl] -> SetTitle(Form("n-like_%s/%s",pname[0],pname[1]));
          if (draw_like) ejungwoo::addsame(ename,histnlike[isys][inl],"hist");
        }

        ename = Form("sys%d_p-like",fSystems[isys]);
        fnplike.setmaint(ename);
        if (draw_like) ejungwoo::add(ename,fnplike.draw());
        for (auto ipl=0; ipl<3; ++ipl) {
          const char *pname[2] = {particleNames[ipidpl[ipl][0]],particleNames[ipidpl[ipl][1]]};
          TString ename3 = Form("sys%d_p-like_%s/%s",fSystems[isys],pname[0],pname[1]);
          for (auto idx : {0,1}) {
            auto ipid = ipidpl[ipl][idx];
            if (ipid==6)
              continue;
            setpar_syspid(isys,ipid);
            histpdist[isys][ipl][idx] = vptoa.draw(fTreePID[isys][ipid]);
            histpdist[isys][ipl][idx] -> SetTitle(pname[idx]);
            //ejungwoo::addsame(ename3,histpdist[isys][ipl][idx]);
            if (idx==0)
              histplike[isys][ipl] = vptoa.draw(fTreePID[isys][ipid]);
          }
          auto ecvs = ejungwoo::findc(ename3);
          if (ecvs!=nullptr)
            ecvs-> setTitle(Form("sys%d n-like %s/%s",fSystems[isys],pname[0],pname[1]));

          if (ipidpl[ipl][1]!=6)
            histplike[isys][ipl] -> Divide(histpdist[isys][ipl][1]);

          histplike[isys][ipl] -> SetTitle(Form("p-like_%s/%s",pname[0],pname[1]));
          if (draw_like) ejungwoo::addsame(ename,histplike[isys][ipl],"hist");
        }
      }

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
            histnp[isys][idx] = (TH1D *) histnlike[isys][inl] -> Clone(npname);
            histnp[isys][idx] -> SetTitle(Form("%d. (%s/%s) / (%s/%s)",idx+1,particleNames[ipidnl[inl][0]],particleNames[ipidnl[inl][1]],particleNames[ipidpl[ipl][0]],particleNames[ipidpl[ipl][1]]));
            histnp[isys][idx] -> Divide(histplike[isys][ipl]);
                 if (inl==0) histnp[isys][idx] -> SetMarkerStyle(5);
            else if (inl==1) histnp[isys][idx] -> SetMarkerStyle(27);
            else if (inl==2) histnp[isys][idx] -> SetMarkerStyle(24);
            if (draw_npratio) ejungwoo::addsame(ename,histnp[0][idx],"hist");
            idx++;
          }
      }

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
          //auto ecvs = ejungwoo::add(ename,fdbratio.draw());
          auto ecvs = fdbratio.drawadd(ename,"gridxgridy");
          //fdbratio.drawaddsame(ename,"gridxgridy");
          //if (isyscom==0) ecvs -> legendlb();
          ecvs -> legendxx();
        }
        vector<TH1 *> histdrarray;
        for (auto idx=0; idx<9; ++idx) {
          TString drname = Form("DR%d",idx);
          auto histdr = (TH1D *) histnp[isys1][idx] -> Clone(drname);
          histdr -> Divide(histnp[isys2][idx]);
          if (draw_dbratio) ejungwoo::addsame(ename,histdr,drawoption+"plhist");
          if (idx==0) ejungwoo::addsame(ename,histdr,drawoption+"plhist");
          histdrarray.push_back(histdr);
        }


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
        //if (draw_dbratio) ejungwoo::addsame(ename,fitydrmean[isyscom],drawoption+"colorxl","average fit");
      }
    }

    if (draw_r21)
      for (auto isyscom : {0,1,2,3})
      {
        int isys1, isys2;
        if (isyscom==0) { isys1=0; isys2=1; }
        else if (isyscom==1) { isys1=3; isys2=2; }
        else if (isyscom==2) { isys1=0; isys2=2; }
        else if (isyscom==3) { isys1=3; isys2=1; }
        //132, 108, 112, 124

        TString systga1 = Form("(%d+%d)",fSystems[isys1],fTargetA[isys1]);
        TString systga2 = Form("(%d+%d)",fSystems[isys2],fTargetA[isys2]);
        TString ename = Form("r21_%d_%d",fSystems[isys1],fSystems[isys2]);
        fptoar21.setmaint(maintitle);
        //fptoar21.setmaint(systga1+" / "+systga2);
        fptoar21.setyt(Form("R_{21}(%s / %s)",systga1.Data(),systga2.Data()));
        fptoar21.drawadd(ename,"gridxgridy");
        double ymax1 = 0, ymax2 = 0;
        for (int ipid : fIPIDs) {
          setpar_syspid(isys1,ipid); auto hist1 = vptoa.draw(fTreePID[isys1][ipid]);
          setpar_syspid(isys2,ipid); auto hist2 = vptoa.draw(fTreePID[isys2][ipid]);
          auto graphs = draw_is(ipid, hist1, hist2, num_tracks_per_event_cut);
          ejungwoo::add(ename, 0, graphs.At(0), "pl", fParticleNames[ipid]);
        }
      }

    TH1 *hist_n0[4];
    if (draw_psdon)
    {
      for (auto isys : {0,1,2,3}) {
        TString systga = Form("%d+%d",fSystems[isys],fTargetA[isys]);

        setpar_syspid(isys,0);
        //auto hist_p = (TH1 *) vptoa.drawaddsame(fTreePID[isys][0],TString("yield"),"plhist",TString("p")  +fSystems[isys])->fLastObject;
        //hist_p -> SetMarkerStyle(20);
        auto hist_p = (TH1 *) vptoa.draw(fTreePID[isys][0]);

        setpar_syspid(isys,2);
        //auto hist_t = (TH1 *) vptoa.drawaddsame(fTreePID[isys][2],TString("yield"),"plhist",TString("t")  +fSystems[isys])->fLastObject;
        //hist_t -> SetMarkerStyle(21);
        auto hist_t = (TH1 *) vptoa.draw(fTreePID[isys][2]);

        setpar_syspid(isys,3);
        //auto hist_he3 = (TH1 *) vptoa.drawaddsame(fTreePID[isys][3],TString("yield"),"plhist",TString("he3")+fSystems[isys])->fLastObject;
        //auto ecvs = vptoa.drawaddsame(fTreePID[isys][3],TString("yield"),"plhist",TString("he3")+fSystems[isys]);
        //auto hist_he3 = (TH1 *) ecvs -> fLastObject;
        //hist_he3 -> SetMarkerStyle(22);
        //ecvs -> legendlt();
        auto hist_he3 = (TH1 *) vptoa.draw(fTreePID[isys][3]);

        
        setpar_syspid(isys,2);
        auto hist_sr = vptoa.draw(fTreePID[isys][2]);
        hist_sr -> SetMarkerStyle(20);
        hist_sr -> SetMarkerColor(ejungwoo::colori(isys));
        hist_sr -> SetLineColor(ejungwoo::colori(isys));
        hist_sr -> Divide(hist_he3);
        fptoasr.setmaint(maintitle);
        if (draw_srn0) fptoasr.drawadd("sr");
        if (draw_srn0) ejungwoo::addsame("sr",hist_sr,"colorx plhist",systga);
        //ejungwoo::add("sr",hist_sr,"",systga);

        setpar_syspid(isys,2);
        hist_n0[isys] = vptoa.draw(fTreePID[isys][2]);
        hist_n0[isys] -> SetMarkerStyle(20);
        hist_n0[isys] -> SetMarkerColor(ejungwoo::colori(isys));
        hist_n0[isys] -> SetLineColor(ejungwoo::colori(isys));
        hist_n0[isys] -> Divide(hist_he3);
        hist_n0[isys] -> Multiply(hist_p);
        hist_n0[isys] = ejungwoo::dndx(hist_n0[isys]);
        fptoan0.setmaint(maintitle);
        if (draw_srn0) fptoan0.drawadd("n0");
        if (draw_srn0) ejungwoo::addsame("n0",hist_n0[isys],"colorx plhist",systga);
      }
    }

    if (draw_cici)
    {
      //ejungwoo::add(ename, 0, new_h(vptoa.getHistName()+"_cin_allframe", ttlCIM, binnCIx, binnCIy2), "addx rangex logy ");
      //ejungwoo::add(ename, 1, new_h(vptoa.getHistName()+"_cip_allframe", ttlCIM, binnCIx, binnCIy2), "addx rangex logy ");
      //ejungwoo::add(ename, 2, new_h(vptoa.getHistName()+"_cir_allframe", ttlCIR, binnCIx, binnCIr), "addx rangex ");

      TString ename;

      TGraphErrors *graphs_rnp[4];
      for (auto isys : {0,1,2,3}) {
        auto sys = fSystems[isys];
        auto tga = fTargetA[isys];
        vptoa.setmaint("($$(sys)) particle multiplicity");
        titles ttlCIMN("($$(sys)) CI neutrons","KE_{CM} (MeV)",Form("#frac{dM(%d)_{n,CI}}{dKE_{CM} #Delta#Omega_{CM}}",sys));
        titles ttlCIMP("($$(sys)) CI protons","KE_{CM} (MeV)",Form("#frac{dM(%d)_{p,CI}}{dKE_{CM} #Delta#Omega_{CM}}",sys));
        titles ttlCIRs("($$(sys)) CI n/p","KE_{CM} (MeV)",Form("R(%d) = #frac{dM(%d)_{n,CI}}{dKE_{CM} #Delta#Omega_{CM}} / #frac{dM(%d)_{p,CI}}{dKE_{CM} #Delta#Omega_{CM}}",sys,sys,sys));
        vector<TH1*> hists;
        ename = TString("ci_pdist")+isys;
        //fcidist.setmaint("CI-n, CI-p");
        fcidist.setmaint(maintitle);
        for (auto ipid : fIPIDAll) {
          setpar_syspid(isys,ipid);
          auto hist = vptoa.draw(fTreePID[isys][ipid]);
          hist = ejungwoo::dndx(hist);
          hist -> SetMarkerStyle(20);
          hists.push_back(hist);
          if (draw_dist) ejungwoo::add(ename, 0, hist, "plhist grid", fParticleNames[ipid]);
        }
        auto hist_n0c = (TH1 *) hist_n0[isys] -> Clone();
        if (draw_dist) ejungwoo::add(ename, 0, hist_n0c, "plhist grid", "pseudo-n");

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
        fcidist.setmaint("CI-n,CI-p");
        if (ejungwoo::findc(ename)==nullptr) {
          auto ecvs = fcidist.drawaddnext(ename);
          ecvs -> legendlt();
        }
        ejungwoo::addsame(ename, graphs.At(0), "colorx pl", Form("CI-n (%d+%d)",sys,tga));
        ejungwoo::addsame(ename, graphs.At(1), "colorx pl", Form("CI-p (%d+%d)",sys,tga));

        ename = "CIN_ov_CIP";
        //fciratio.setmaint("CI-n / CI-p");
        fciratio.setmaint(maintitle);
        if (ejungwoo::findc(ename)==nullptr) fciratio.drawaddnext(ename);
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

          auto graphr = draw_ratio(graphs_rnp[isys1], graphs_rnp[isys2],vptoa.getBinn().max);
          graphr -> SetMarkerColor(color);
          graphr -> SetLineColor(color);
          auto title1 = Form("DR(CI-n/CI-p)[(%d+%d)/(%d+%d)]", fSystems[isys1],fTargetA[isys1],fSystems[isys2],fTargetA[isys2]);
          ejungwoo::add(ename, graphr, "pl colorx", title1);

          auto title2 = Form("DR(n-like/p-like)[(%d+%d)/(%d+%d)]", fSystems[isys1],fTargetA[isys1],fSystems[isys2],fTargetA[isys2]);
          fitydrmean[isyscom] -> SetLineColor(color);
          auto ecvs = ejungwoo::addsame(ename, fitydrmean[isyscom], "colorx l addx", title2);
          ecvs -> legendlb();
        }

        //graphr = draw_ratio(graphs_rnp[2], graphs_rnp[1],vptoa.getBinn().max);
        //graphr -> SetMarkerStyle(20);
        //ejungwoo::add(ename, graphr, "pl addx");
      }
    }

  }

  ejungwoo::drawsaveall("cvsl","pdf");
}
