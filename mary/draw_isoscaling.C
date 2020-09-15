#include "/Users/ejungwoo/config/ejungwoo.h"
#include "init_variables.h"
#include "functions.h"

const int fIPIDs[] = {0,1,2,3,4};
const int fISystems[] = {0,1};

void draw_isoscaling()
{
  ejungwoo::gstat(0);
  ejungwoo::gsave(0);
  ejungwoo::gsetvmark(0);
  ejungwoo::gshortprint();

  init();

  //for (auto aversion : {"kleft"})
  for (auto aversion : {"kright132"})
  {
    auto condition_array = setversion(aversion);
    settrees();

    for (auto ii : {0,1,2,3})
    {
      int isys1, isys2;
           if (ii==0) { isys1=0; isys2=1; }
      else if (ii==1) { isys1=3; isys2=2; }
      else if (ii==2) { isys1=0; isys2=2; }
      else if (ii==3) { isys1=3; isys2=1; }

      auto systgA1 = Form("%d+%d",fSystems[isys1],fTargetA[isys1]);
      auto systgA2 = Form("%d+%d",fSystems[isys2],fTargetA[isys2]);

      auto var = &fvar_foby_cm;
      setpar("pname","all");
      TString cname_is = Form("IS%d_",ii)+var -> getName()+"__"+aversion;

      auto name00 = cname_is + "frame";
      auto title00 = titles(Form("isoscaling plot %d/%d",fSystems[isys1],fSystems[isys2]),var->getTitle().x,Form("p(%s)/p(%s)",systgA1,systgA2));
      ejungwoo::add(cname_is,0,new_h(name00,title00,var->getBinn(),binning(100,0,2)),"addx");
      double ymax1 = 0, ymax2 = 0;
      for (int ipid : fIPIDs) {
        setpar_syspid(isys1,ipid); auto hist1 = var -> draw(fTreePID[isys1][ipid]); hist1 -> SetMarkerColor(ejungwoo::colori(ipid)); hist1 -> SetMarkerStyle(20);
        setpar_syspid(isys2,ipid); auto hist2 = var -> draw(fTreePID[isys2][ipid]); hist2 -> SetMarkerColor(ejungwoo::colori(ipid)); hist2 -> SetMarkerStyle(24); hist2 -> SetMarkerSize(1.5);
        auto graphs = draw_is(ipid, hist1, hist2, 0.01);
        ejungwoo::add(cname_is, 0, graphs.At(0), "pl", fParticleNames[ipid]);
        ejungwoo::add(cname_is, 1, hist1, "pl colorx", fParticleNames[ipid]+fSystems[isys1]);
        ejungwoo::add(cname_is, 1, hist2, "pl colorx", fParticleNames[ipid]+fSystems[isys2]);
      }

      auto var2 = &fvar_pid;
      TH2D *histpid1 = nullptr;
      TH2D *histpid2 = nullptr;
      var2 -> setmaint("$$(sys)");
      for (int ipid : fIPIDs) {
        setpar_syspid(isys1,ipid); auto hist1 = var2 -> draw(fTreePID[isys1][ipid]); if (histpid1==nullptr) histpid1 = (TH2D *) hist1; else histpid1 -> Add(hist1);
        setpar_syspid(isys2,ipid); auto hist2 = var2 -> draw(fTreePID[isys2][ipid]); if (histpid2==nullptr) histpid2 = (TH2D *) hist2; else histpid2 -> Add(hist2);
      }
      ejungwoo::add(cname_is, 2, histpid1, "gridx gridy logz", TString("pid_")+fSystems[isys1]);
      ejungwoo::add(cname_is, 3, histpid2, "gridx gridy logz", TString("pid_")+fSystems[isys2]);

      auto var3 = fvar_foby_cm.add(fvar_poz_lab);
      for (auto ipid : {3} ){
        //auto hist1 = var3 -> draw(fTreePID[isys1][ipid]); ejungwoo::add(cname_is, 4, hist1, "", fSystems[isys1]);
        //auto hist2 = var3 -> draw(fTreePID[isys2][ipid]); ejungwoo::add(cname_is, 5, hist1, "", fSystems[isys2]);
      }
    }
  }

  /******************************************************************************************/

  ejungwoo::drawsaveall("cvsl","png");
}
