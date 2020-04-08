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

  bool draw_bb = false;
  bool draw_p5 = true;

  //TString region_fit_names[] = { "pin", "pip", "p", "d", "t", "he3", "he4", "he6", "li6", "li7", "ele", "elp", };

  auto file = new TFile("/Users/ejungwoo/spirit/the3/mary/pid/data__sys108/hpid108_nn1000x1000_p2500_dedx1800.sys108.root");
  auto hist = (TH2D *) file -> Get("hpid108_nn1000x1000_p2500_dedx1800");
  hist -> SetName("hpid");
  hist -> SetTitle("pid;p;dedx");

  auto cvs = ejungwoo::cc4();
  cvs -> SetLogz();
  hist -> Draw("colz");

  //for (auto name : region_fit_names)
  int idx_particles[] = {1,2,3,4,5,6,7,8,9};
  for (auto idx : idx_particles)
  {
    TString name = stp::fNameShort[idx];

    if (0) {
      auto file = new TFile(TString("data__prefit/graph_pdedx_region_")+name+".prefit.root");
      auto graph = (TGraphErrors *) file -> Get(TString("graph_pdedx_region_")+name);
      graph -> Draw("samepz1");
    }

    if (draw_bb) {
      file = new TFile(TString("data__prefit/fit_dedx_region_")+name+".prefit.root");
      auto f1_with_parameter = (TF1 *) file -> Get(TString("fit_dedx_region_")+name);
      if (1) {
        f1_with_parameter -> SetNpx(1000);
        f1_with_parameter -> SetLineColor(kGray);
        f1_with_parameter -> Draw("samel");
      }
      else {
        TF1 *fit_dedx;
        /**/ if (idx== 0) { fit_dedx = new TF1(TString("fit_dedx_")+name, stp::bb_dedx_pin, p_hist_min, p_hist_max, 2); }
        else if (idx== 1) { fit_dedx = new TF1(TString("fit_dedx_")+name, stp::bb_dedx_pip, p_hist_min, p_hist_max, 2); fit_dedx -> SetParameters(0.815692,-48.1502);}
        else if (idx== 2) { fit_dedx = new TF1(TString("fit_dedx_")+name, stp::bb_dedx_p  , p_hist_min, p_hist_max, 2); fit_dedx -> SetParameters(0.815692,-48.1502);}
        else if (idx== 3) { fit_dedx = new TF1(TString("fit_dedx_")+name, stp::bb_dedx_d  , p_hist_min, p_hist_max, 2); fit_dedx -> SetParameters(0.815692,-48.1502);}
        else if (idx== 4) { fit_dedx = new TF1(TString("fit_dedx_")+name, stp::bb_dedx_t  , p_hist_min, p_hist_max, 2); fit_dedx -> SetParameters(0.815692,-48.1502);}
        else if (idx== 5) { fit_dedx = new TF1(TString("fit_dedx_")+name, stp::bb_dedx_he3, p_hist_min, p_hist_max, 2); fit_dedx -> SetParameters(0.0350209,-1397.76);}
        else if (idx== 6) { fit_dedx = new TF1(TString("fit_dedx_")+name, stp::bb_dedx_he4, p_hist_min, p_hist_max, 2); fit_dedx -> SetParameters(0.0350209,-1397.76);}
        else if (idx== 7) { fit_dedx = new TF1(TString("fit_dedx_")+name, stp::bb_dedx_he6, 1000      , p_hist_max, 2); fit_dedx -> SetParameters(0.0350209,-1397.76);}
        else if (idx== 8) { fit_dedx = new TF1(TString("fit_dedx_")+name, stp::bb_dedx_li6, p_hist_min, p_hist_max, 2); fit_dedx -> SetParameters(0.00112263,-45209.6);}
        else if (idx== 9) { fit_dedx = new TF1(TString("fit_dedx_")+name, stp::bb_dedx_li7, p_hist_min, p_hist_max, 2); fit_dedx -> SetParameters(0.00112263,-45209.6);}
        else if (idx==10) { fit_dedx = new TF1(TString("fit_dedx_")+name, stp::bb_dedx_ele, p_hist_min, p_hist_max, 2); }
        else if (idx==11) { fit_dedx = new TF1(TString("fit_dedx_")+name, stp::bb_dedx_pos, p_hist_min, p_hist_max, 2); }
        cout << name << " " << f1_with_parameter->GetParameter(0) << " " << f1_with_parameter->GetParameter(1) << endl;
        fit_dedx -> SetParameters(f1_with_parameter->GetParameter(0), f1_with_parameter->GetParameter(1));
        fit_dedx -> SetNpx(1000);
        fit_dedx -> Draw("samel");
      }
    }

    if (draw_p5)
    {
      file = new TFile(TString("data__prefit/fit_dedx_par5_region_")+name+".prefit.root");
      auto f1_with_parameter = (TF1 *) file -> Get(TString("fit_dedx_par5_region_")+name);
      if (0) {
        f1_with_parameter -> SetNpx(1000);
        f1_with_parameter -> Draw("samel");
      }
      else {
        TF1 *fit_dedx;
        /**/ if (idx== 0) { fit_dedx = new TF1(TString("fit_dedx_par5_")+name, stp::p5_dedx_pin, p_hist_min, p_hist_max, 5); }
        else if (idx== 1) { fit_dedx = new TF1(TString("fit_dedx_par5_")+name, stp::p5_dedx_pip, p_hist_min, p_hist_max, 5); }
        else if (idx== 2) { fit_dedx = new TF1(TString("fit_dedx_par5_")+name, stp::p5_dedx_p  , p_hist_min, p_hist_max, 5); }
        else if (idx== 3) { fit_dedx = new TF1(TString("fit_dedx_par5_")+name, stp::p5_dedx_d  , p_hist_min, p_hist_max, 5); }
        else if (idx== 4) { fit_dedx = new TF1(TString("fit_dedx_par5_")+name, stp::p5_dedx_t  , p_hist_min, p_hist_max, 5); }
        else if (idx== 5) { fit_dedx = new TF1(TString("fit_dedx_par5_")+name, stp::p5_dedx_he3, p_hist_min, p_hist_max, 5); }
        else if (idx== 6) { fit_dedx = new TF1(TString("fit_dedx_par5_")+name, stp::p5_dedx_he4, p_hist_min, p_hist_max, 5); }
        else if (idx== 7) { fit_dedx = new TF1(TString("fit_dedx_par5_")+name, stp::p5_dedx_he6, 1000      , p_hist_max, 5); }
        else if (idx== 8) { fit_dedx = new TF1(TString("fit_dedx_par5_")+name, stp::p5_dedx_li6, p_hist_min, p_hist_max, 5); }
        else if (idx== 9) { fit_dedx = new TF1(TString("fit_dedx_par5_")+name, stp::p5_dedx_li7, p_hist_min, p_hist_max, 5); }
        else if (idx==10) { fit_dedx = new TF1(TString("fit_dedx_par5_")+name, stp::p5_dedx_ele, p_hist_min, p_hist_max, 5); }
        else if (idx==11) { fit_dedx = new TF1(TString("fit_dedx_par5_")+name, stp::p5_dedx_pos, p_hist_min, p_hist_max, 5); }
        fit_dedx -> SetParameters(f1_with_parameter->GetParameter(0), f1_with_parameter->GetParameter(1), f1_with_parameter->GetParameter(2), f1_with_parameter->GetParameter(3), f1_with_parameter->GetParameter(4));
        fit_dedx -> SetNpx(1000);
        fit_dedx -> Draw("samel");
      }
    }
  }
}
