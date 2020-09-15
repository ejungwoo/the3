#include "/Users/ejungwoo/config/ejungwoo.h"
#include "init_variables.h"
#include "functions.h"

const int fIPIDs[] = {0,1,2,3,4};

void draw_j()
{
  ejungwoo::gcvspos(1300);
  ejungwoo::gstat(0);
  ejungwoo::gsave(0);
  ejungwoo::gsetvmark(0);
  //ejungwoo::gshortprint();

  bool draw_pydist = 0;
  bool draw_nplike = 0;
  bool draw_npratio = 0;
  bool draw_dbratio = 1;
  bool draw_r21 = 1;
  bool draw_psdon = 1;
  bool draw_cici = 1;
  bool draw_dist = 0;

  if (draw_cici) draw_psdon = 1;

  double num_tracks_per_event_cut = 0.01;

  init();

  const char *particleNames[] = {"p","d","t","he3","he4","he6","1"};

  //auto vptoa = variable("ptoa_cm", "pt_cm/$$(a)", "$$(cut0)", titles("","p_{T}/A (MeV/c)",""), binning(15,100,400));
  auto vptoa = variable("ptoa_cm", "pt_cm/$$(a)", "$$(cut0)", titles("","p_{T}/A (MeV/c)",""), binning(10,100,400));

  auto fnplike = variable("fnplike", "", "", titles(" ","p_{T}/A (MeV/c)","p0/p1"), 100,100,400,100,0,4);
  auto fnpratio = variable("fnpratio", "", "", titles(" ","p_{T}/A (MeV/c)","n-like/p-like"), 100,100,400,100,0,20);
  auto fdbratio = variable("fdbratio", "", "", titles(" ","p_{T}/A (MeV/c)","DR(n-like/p-like)$$(systga21)"), 100,100,400,100,0,2.0);
  auto fptoar21 = variable("fptoar21", "", "", titles(" ","p_{T}/A (MeV/c)","R21"), 100,100,400,100,0.5,2);
  auto fptoasr = variable("fptoasr", "", "", titles(" ","p_{T}/A (MeV/c)","SR(t/he3)"), 100,100,400,100,0,7);
  auto fptoan0 = variable("fptoan0", "", "", titles(" ","p_{T}/A (MeV/c)","#frac{d^{2}M}{dyd(p_{T}/A)}"), 100,100,400,100,0,0.06);
  auto fcidist = variable("cidist", "", "", titles(" ","p_{T}/A (MeV/c)","#frac{d^{2}M}{dyd(p_{T}/A)}"), 100,100,400,100,0.0003,0.0025);
  auto fciratio = variable("cidist", "", "", titles(" ","p_{T}/A (MeV/c)","ratio"), 100,100,400,100,0,2.0);

  int ipidnl[3][2] = {{4,3},{1,0},{2,1}};
  int ipidpl[3][2] = {{0,6},{4,2},{3,1}};

  int isystems[] = {0,1};

  for (auto aversion : {"kright"})
  {
    auto condition_array = setversion(aversion);
    settrees();

    if (draw_pydist)
      for (auto isys : {0,1,2,3}) {
        TString ename = TString("sysCM")+fSystems[isys];
        ename = TString("sysLab")+fSystems[isys];
        for (int ipid : fIPIDs) {
          setpar_syspid(isys,ipid);
          fvar_ypt_lab.drawaddnext(fTreePID[isys][ipid],ename);
        }
      }

    setpar("rap_cut","$$(foby_lab)>.4&&$$(foby_lab)<.6");

    TGraphErrors *graphdr[4];
    TF1 *fitydr[4];
    if (draw_nplike || draw_npratio || draw_dbratio)
    {
      TH1 *histnlike[4][3]; // isys, inpl
      TH1 *histplike[4][3]; // isys, inpl
      TH1 *histndist[4][3][3]; // isys, inpl, idx
      TH1 *histpdist[4][3][3]; // isys, inpl, idx

      for (auto isys : {0,1,2,3})
      {
        cout << isys << endl;
        TString ename = Form("sys%d_n-like",fSystems[isys]);
        fnplike.setmaint(ename);
        if (draw_nplike) ejungwoo::add(ename,fnplike.draw());
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
          if (draw_nplike) ejungwoo::addsame(ename,histnlike[isys][inl],"hist");
        }

        ename = Form("sys%d_p-like",fSystems[isys]);
        fnplike.setmaint(ename);
        if (draw_nplike) ejungwoo::add(ename,fnplike.draw());
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
          if (draw_nplike) ejungwoo::addsame(ename,histplike[isys][ipl],"hist");
        }
      }

      TH1D *histnp[4][9];
      for (auto isys : {0,1,2,3})
      {
        TString ename = Form("sys%d_n-like/p-like",fSystems[isys]);
        fnpratio.setmaint(ename);

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
        TString drawoption = "addx";
        if (isyscom==0) { isys1=0; isys2=1; drawoption=""; }
        else if (isyscom==1) { isys1=3; isys2=2; }
        else if (isyscom==2) { isys1=0; isys2=2; }
        else if (isyscom==3) { isys1=3; isys2=1; }

        TString ename = TString("DR")+isyscom;
        fdbratio.setmaint(Form("(%d+%d) / (%d+%d)", fSystems[isys1],fTargetA[isys1],fSystems[isys2],fTargetA[isys2]));
        ejungwoo::setpar("systga21",Form("_{#frac{%d+%d}{%d+%d}}", fSystems[isys1],fTargetA[isys1],fSystems[isys2],fTargetA[isys2]));
        if (draw_dbratio) {
          auto ecvs = ejungwoo::add(ename,fdbratio.draw());
          if (isyscom==0) ecvs -> legendlb();
        }
        vector<TH1 *> histdrarray;
        for (auto idx=0; idx<9; ++idx) {
          TString drname = Form("DR%d",idx);
          auto histdr = (TH1D *) histnp[isys1][idx] -> Clone(drname);
          histdr -> Divide(histnp[isys2][idx]);
          if (draw_dbratio) ejungwoo::addsame(ename,histdr,drawoption+"plhist");
          histdrarray.push_back(histdr);
        }


        graphdr[isyscom] = ejungwoo::new_ge();
        auto binn = binning(histdrarray.at(0));
        binn.resetb();
        while (binn.nextb()) {
          double yvalue = 0;
          for (auto hist : histdrarray) {
            auto value = hist -> GetBinContent(binn.idx);
            yvalue += value;
          }
          yvalue = yvalue/9.;
          graphdr[isyscom] -> SetPoint(graphdr[isyscom]->GetN(),binn.value,yvalue);
        }
        fitydr[isyscom] = new TF1(Form("fitydr%d",isyscom),"[0]",binn.min,binn.max);
        graphdr[isyscom] -> Fit(fitydr[isyscom],"QR0N");

        graphdr[isyscom] -> SetLineColor(ejungwoo::colori(isyscom+4));
        graphdr[isyscom] -> SetMarkerStyle(25);
        graphdr[isyscom] -> SetMarkerSize(2);
        graphdr[isyscom] -> SetLineStyle(9);
        graphdr[isyscom] -> SetLineWidth(2);
        fitydr[isyscom] -> SetLineColor(ejungwoo::colori(isyscom+4));
        fitydr[isyscom] -> SetLineStyle(2);
        fitydr[isyscom] -> SetLineWidth(2);

        if (draw_dbratio) ejungwoo::addsame(ename,graphdr[isyscom],drawoption+"colorxp","average");
        //if (draw_dbratio) ejungwoo::addsame(ename,fitydr[isyscom],drawoption+"colorxl","average fit");
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
        fptoar21.setmaint(systga1+" / "+systga2);
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

        setpar_syspid(isys,0); auto hist_p   = (TH1 *) vptoa.drawaddsame(fTreePID[isys][0],TString("yield"),"plhist",TString("p")  +fSystems[isys])->fLastObject; hist_p   -> SetMarkerStyle(20);
        setpar_syspid(isys,2); auto hist_t   = (TH1 *) vptoa.drawaddsame(fTreePID[isys][2],TString("yield"),"plhist",TString("t")  +fSystems[isys])->fLastObject; hist_t   -> SetMarkerStyle(20);
        setpar_syspid(isys,3); auto hist_he3 = (TH1 *) vptoa.drawaddsame(fTreePID[isys][3],TString("yield"),"plhist",TString("he3")+fSystems[isys])->fLastObject; hist_he3 -> SetMarkerStyle(20);
        
        setpar_syspid(isys,2); auto hist_sr = vptoa.draw(fTreePID[isys][2]);
        hist_sr -> SetMarkerStyle(20);
        hist_sr -> Divide(hist_he3);
        fptoasr.drawadd("sr");
        ejungwoo::addsame("sr",hist_sr,"plhist",systga);
        //ejungwoo::add("sr",hist_sr,"",systga);

        setpar_syspid(isys,2);
        hist_n0[isys] = vptoa.draw(fTreePID[isys][2]);
        hist_n0[isys] -> SetMarkerStyle(20);
        hist_n0[isys] -> Divide(hist_he3);
        hist_n0[isys] -> Multiply(hist_p);
        hist_n0[isys] = ejungwoo::dndx(hist_n0[isys]);
        fptoan0.setmaint("pseudo neutron");
        fptoan0.drawadd("n0");
        ejungwoo::addsame("n0",hist_n0[isys],"plhist",systga);
      }
    }

    if (draw_cici)
    {
      //ejungwoo::add(ename, 0, new_h(vptoa.getHistName()+"_cin_allframe", ttlCIM, binnCIx, binnCIy2), "addx rangex logy ");
      //ejungwoo::add(ename, 1, new_h(vptoa.getHistName()+"_cip_allframe", ttlCIM, binnCIx, binnCIy2), "addx rangex logy ");
      //ejungwoo::add(ename, 2, new_h(vptoa.getHistName()+"_cir_allframe", ttlCIR, binnCIx, binnCIr), "addx rangex ");

      TGraphErrors *graphs_rnp[4];
      for (auto isys : {0,1,2,3}) {
        auto sys = fSystems[isys];
        auto tga = fTargetA[isys];
        vptoa.setmaint("($$(sys)) particle multiplicity");
        titles ttlCIMN("($$(sys)) CI neutrons","KE_{CM} (MeV)",Form("#frac{dM(%d)_{n,CI}}{dKE_{CM} #Delta#Omega_{CM}}",sys));
        titles ttlCIMP("($$(sys)) CI protons","KE_{CM} (MeV)",Form("#frac{dM(%d)_{p,CI}}{dKE_{CM} #Delta#Omega_{CM}}",sys));
        titles ttlCIRs("($$(sys)) CI n/p","KE_{CM} (MeV)",Form("R(%d) = #frac{dM(%d)_{n,CI}}{dKE_{CM} #Delta#Omega_{CM}} / #frac{dM(%d)_{p,CI}}{dKE_{CM} #Delta#Omega_{CM}}",sys,sys,sys));
        vector<TH1*> hists;
        TString ename = TString("ci_pdist")+isys;
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

        ename = "cinp";
        fcidist.setmaint("CI-neutrons,CI-protons");
        if (ejungwoo::findc(ename)==nullptr) fcidist.drawaddnext(ename);
        ejungwoo::addsame(ename, graphs.At(0), "colorx pl", Form("CI-neutrons (%d+%d)",sys,tga));
        ejungwoo::addsame(ename, graphs.At(1), "colorx pl", Form("CI-protons (%d+%d)",sys,tga));

        ename = "cinop";
        fciratio.setmaint("CI-neutrons / CI-protons");
        if (ejungwoo::findc(ename)==nullptr) fciratio.drawaddnext(ename);
        ejungwoo::addsame(ename, graphs.At(2), "colorx pl", Form("CI-n/p (%d+%d)",sys,tga));

        graphs_rnp[isys] = (TGraphErrors *) graphs.At(2);
      }

      for (auto isyscom : {0,1,2,3}) {
        int isys1, isys2;
             if (isyscom==0) { isys1=0; isys2=1; }
        else if (isyscom==1) { isys1=3; isys2=2; }
        else if (isyscom==2) { isys1=0; isys2=2; }
        else if (isyscom==3) { isys1=3; isys2=1; }

        auto title = Form("DR(n/p)[(%d+%d)/(%d+%d)]", fSystems[isys1],fTargetA[isys1],fSystems[isys2],fTargetA[isys2]);
        //auto ecvs = ejungwoo::addsame("cinop", graphdr[isyscom], "colorx l", title);
        auto ecvs = ejungwoo::addsame("cinop", fitydr[isyscom], "colorx l", title);
        ecvs -> legendlb();
      }


      if (0) {
        auto graphr = draw_ratio(graphs_rnp[0], graphs_rnp[1],vptoa.getBinn().max);
        graphr -> SetMarkerStyle(20);
        ejungwoo::add("ci_r132o108", graphr, "pl addx");
      }
    }

  }

  ejungwoo::drawsaveall("cvsl","pdf");
}
