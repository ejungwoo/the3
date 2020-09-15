#include "/Users/ejungwoo/config/ejungwoo.h"
#include "init_variables.h"
#include "functions.h"

const int fIPIDs[] = {0,1,2,3,4};
//const int fIPIDs[] = {0};
//const int fIPIDs[] = {0};

void draw_yp()
{
  ejungwoo::gstat(0);
  ejungwoo::gsave(0);
  ejungwoo::gsetvmark(0);

  init();

  auto varyp = ejungwoo::variable("yp;yp$$(sys)$$(pname)"
      ,"p_lab:fy_cm/(by_cm/2)"
      ,""
      ,titles(fttl_main, "y_{0} = y_{Frag.CM}/y_{NN}", "p/Z (MeV/c)")
      ,binning(200,-2,2)
      ,binning(200,0,3000.));

  for (auto aversion : {"kleft"})
  {
    auto condition_array = setversion(aversion);
    settrees();

    for (auto isys : {0})
    {
      TString ename = TString("yp")+isys;
      for (int ipid : fIPIDs) {
        setpar_syspid(isys,ipid);
        auto hist = (TH2D *) varyp.drawaddnext(fTreePID[isys][ipid],ename) -> fLastObject;
        auto hana2 = ejungwoo::histana2(hist,1,20,"m");
        auto graph = hana2.graphana();
        ejungwoo::addsame(ename,graph,"p colorx", "graph");
        auto fit = ejungwoo::fitgraph(graph,"pol2");
        ejungwoo::addsame(ename,fit, "l colorx", "fit");

        auto dummygraph = new TGraph();
        dummygraph -> SetMarkerStyle(25);
        dummygraph -> SetMarkerSize(kGreen);

        TString fileName = TString("data/")+fit->GetName()+".txt";
        cout << fileName << endl;
        ofstream file(fileName);
        for (auto mom : {200,500,1000,1500,2000,2500})
        {
          if (mom > fit->Eval(-2) && mom < fit->Eval(2)) {
            auto x = fit -> GetX(mom,-2,2);
            dummygraph -> SetPoint(dummygraph->GetN(),x,mom);
            file << mom << " " << x << " " << hist -> Integral() << endl;
          }
          else
            file << mom << " " << -999 << " " << hist -> Integral() << endl;
        }
        ejungwoo::addsame(ename,dummygraph,"colorx");
      }
    }
  }

  ejungwoo::drawsaveall("cvsl", "pdf");
}
