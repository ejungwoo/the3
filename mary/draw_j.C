#include "/Users/ejungwoo/config/ejungwoo.h"
#include "init_variables.h"
#include "functions.h"

const int fIPIDs[] = {0,1,2,3,4};

void draw_j()
{
  ejungwoo::gcvspos(1000);
  ejungwoo::gstat(0);
  ejungwoo::gsave(0);
  ejungwoo::gsetvmark(0);
  //ejungwoo::gshortprint();

  bool draw_pydist = 0;
  bool draw_nplike = 0;
  bool draw_npratio = 0;
  bool draw_dbratio = 0;
  bool draw_r21 = 1;
  bool draw_psdon = 1;

  double num_tracks_per_event_cut = 0.01;

  init();

  const char *particleNames[] = {"p","d","t","he3","he4","he6","1"};

  //auto vptoa = variable("ptoa_cm", "pt_cm/$$(a)", "$$(cut0)", titles("","p_{T}/A (MeV/c)",""), binning(15,100,400));
  auto vptoa = variable("ptoa_cm", "pt_cm/$$(a)", "$$(cut0)", titles("","p_{T}/A (MeV/c)",""), binning(10,100,400));

  auto fnplike = variable("fnplike", "", "", titles(" ","p_{T}/A (MeV/c)","p0/p1"), 100,100,400,100,0,4);
  auto fnpratio = variable("fnpratio", "", "", titles(" ","p_{T}/A (MeV/c)","n-like/p-like"), 100,100,400,100,0,20);
  auto fdbratio = variable("fdbratio", "", "", titles(" ","p_{T}/A (MeV/c)","DR(n-like/p-like)"), 100,100,400,100,0,2.0);
  auto fptoar21 = variable("fptoar21", "", "", titles(" ","p_{T}/A (MeV/c)","R21"), 100,100,400,100,0.5,2);
  auto fptoasr = variable("fptoasr", "", "", titles(" ","p_{T}/A (MeV/c)","SR(t/he3)"), 100,100,400,100,0,7);
  auto fptoan0 = variable("fptoan0", "", "", titles(" ","p_{T}/A (MeV/c)","#frac{d^{2}M}{dyd(p_{t}/A)}"), 100,100,400,100,0,0.06);

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
            ecvs-> SetTitle(Form("sys%d n-like %s/%s",fSystems[isys],pname[0],pname[1]));

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
            ecvs-> SetTitle(Form("sys%d n-like %s/%s",fSystems[isys],pname[0],pname[1]));

          if (ipidpl[ipl][1]!=6)
            histplike[isys][ipl] -> Divide(histpdist[isys][ipl][1]);

          histplike[isys][ipl] -> SetTitle(Form("p-like_%s/%s",pname[0],pname[1]));
          if (draw_nplike) ejungwoo::addsame(ename,histplike[isys][ipl],"hist");
        }
      }

      TH1D *histnp[4][9];
      for (auto isys : {0,1,2,3})
      {
        cout << isys << endl;
        TString ename = Form("sys%d_n-like/p-like",fSystems[isys]);
        fnpratio.setmaint(ename);

        if (draw_npratio) ejungwoo::add(ename,fnpratio.draw());
        int idx = 0;
        for (auto inl=0; inl<3; ++inl) 
          for (auto ipl=0; ipl<3; ++ipl) {
            TString npname = Form("sys%d_nplike%d%d",fSystems[isys],inl,ipl);
            histnp[isys][idx] = (TH1D *) histnlike[isys][inl] -> Clone(npname);
            histnp[isys][idx] -> SetTitle(Form("%d. (%s/%s) / (%s/%s)",idx,particleNames[ipidnl[inl][0]],particleNames[ipidnl[inl][1]],particleNames[ipidpl[ipl][0]],particleNames[ipidpl[ipl][1]]));
            histnp[isys][idx] -> Divide(histplike[isys][ipl]);
            if (draw_npratio) ejungwoo::addsame(ename,histnp[0][idx],"hist");
            idx++;
          }
      }

      for (auto ii : {0,1,2,3})
      {
        cout << ii << endl;
        int isys1, isys2;
        if (ii==0) { isys1=0; isys2=1; }
        else if (ii==1) { isys1=3; isys2=2; }
        else if (ii==2) { isys1=0; isys2=2; }
        else if (ii==3) { isys1=3; isys2=1; }

        TString ename = TString("DR")+ii;
        fdbratio.setmaint(Form("(%d+%d) / (%d+%d)", fSystems[isys1],fTargetA[isys1],fSystems[isys2],fTargetA[isys2]));
        if (draw_dbratio) ejungwoo::add(ename,fdbratio.draw());
        for (auto idx=0; idx<9; ++idx) {
          TString drname = Form("DR%d",idx);
          auto histdr = (TH1D *) histnp[isys1][idx] -> Clone(drname);
          histdr -> Divide(histnp[isys2][idx]);
          histdr -> SetMarkerStyle(5);
          if (draw_dbratio) ejungwoo::addsame(ename,histdr,"plhist");
        }
      }
    }

    if (draw_r21)
      for (auto ii : {0,1,2,3})
      {
        int isys1, isys2;
        if (ii==0) { isys1=0; isys2=1; }
        else if (ii==1) { isys1=3; isys2=2; }
        else if (ii==2) { isys1=0; isys2=2; }
        else if (ii==3) { isys1=3; isys2=1; }
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

        setpar_syspid(isys,2); auto hist_n0 = vptoa.draw(fTreePID[isys][2]);
        hist_n0 -> SetMarkerStyle(20);
        hist_n0 -> Divide(hist_he3);
        hist_n0 -> Multiply(hist_p);
        hist_n0 = ejungwoo::dndx(hist_n0);
        fptoan0.drawadd("n0");
        ejungwoo::addsame("n0",hist_n0,"plhist",systga);
      }
    }
  }

  ejungwoo::drawsaveall("cvsl","pdf");
}
