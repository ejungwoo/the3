#include "stp.h"

void check_prefit()
{
  ejungwoo::gcvspos(1300,10);
  ejungwoo::gstat("nm");
  ejungwoo::gzcolor(1);
  ejungwoo::gversion("check_prefit");
  ejungwoo::gstatleft(0);
  ejungwoo::gverbose(1);
  ejungwoo::gwrite(0);

  double p_hist_min = 0;
  double p_hist_max = 2500;

  bool draw_p5 = true;

  auto file = new TFile("/Users/ejungwoo/spirit/the3/mary/pid/data__sys108/hpid108_nn1000x1000_p2500_dedx1800.sys108.root");
  auto hist = (TH2D *) file -> Get("hpid108_nn1000x1000_p2500_dedx1800");
  hist -> SetName("hpid");
  hist -> SetTitle("pid;p;dedx");

  auto cvs = ejungwoo::cc4();
  cvs -> SetLogz();
  hist -> Draw("colz");

  //int idx_particles[] = {1,2,3,4,5,6,7,8,9,12};
  //int idx_particles[] = {12};
  int idx_particles[] = {stp::kp,stp::kd,stp::kt};
  for (auto idx : idx_particles)
  {
    TString name = stp::fNameShort[idx];

    if (0) {
      auto file = new TFile(TString("data__prefit/graph_pdedx_region_")+name+".prefit.root");
      auto graph = (TGraphErrors *) file -> Get(TString("graph_pdedx_region_")+name);
      graph -> Draw("samepz1");
    }

    if (draw_p5)
    {
      file = new TFile(TString("data__prefit/fit_dedx_par5_region_")+name+".prefit.root");
      auto f1_with_parameter = (TF1 *) file -> Get(TString("fit_dedx_par5_region_")+name);
      TF1 *fit_dedx = stp::f1_dedx_p5(idx,p_hist_min, p_hist_max);
      fit_dedx -> SetParameters(f1_with_parameter->GetParameter(0), f1_with_parameter->GetParameter(1), f1_with_parameter->GetParameter(2), f1_with_parameter->GetParameter(3), f1_with_parameter->GetParameter(4));
      fit_dedx -> SetLineColor(kBlack);
      fit_dedx -> SetNpx(1000);
      fit_dedx -> Draw("samel");
    }
  }
}
