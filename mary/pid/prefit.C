#include "stp.h"

void prefit()
{
  ejungwoo::gcvspos(1300,10);
  ejungwoo::gstat("nm");
  ejungwoo::gversion("prefit");
  ejungwoo::gstatleft(0);

  ejungwoo::gverbose(0);
  ejungwoo::gwrite(1);

  //ejungwoo::gfixcvsxy(1);

  int idx_particles[] = {12};
  //int idx_particles[] = {1,2,3,4,5,6,7,8,9,12};
  //int idx_particles[] = {5,6,7,8,9};
  //int idx_particles[] = {1};
  //int idx_particles[] = {8,9};

  bool fit_projection = true;
  bool draw_projection_fits = false;
  bool draw_pid = true;

  double p_hist_min = 0;
  double p_hist_max = 2500;

  double s_hist_min = 0;
  double s_hist_max = 2500;

  auto sum_bin_s = 20;
  auto sum_bin_p = 20;

  Double_t limits[13][2] = {
    /* pin" */ {100,100},
    /* pip" */ {50,200},
    /* p"   */ {500,500},
    /* d"   */ {800,600},
    /* t"   */ {400,1000},
    /* he3" */ {800,800},
    /* he4" */ {800,800},
    /* he6" */ {800,800},
    /* li6" */ {800,1500},
    /* li7" */ {800,1500},
    /* ele" */ {500,500},
    /* pos" */ {500,500},
    /* li " */ {800,1500},
  };

  auto file = new TFile("/Users/ejungwoo/spirit/the3/mary/pid/data__sys108/hpid108_nn1000x1000_p2500_dedx1800.sys108.root");
  auto hist = (TH2D *) file -> Get("hpid108_nn1000x1000_p2500_dedx1800");
  hist -> SetName("hpid");
  hist -> SetTitle("pid;p;dedx");

  //for (auto name : cut_names)
  //for (auto idx : {0,1})
  for (auto idx : idx_particles)
  {
    auto name_cut = TString("region_") + stp::fNameShort[idx];
    //cout << " ======= " << name_cut << endl;
    TH2D *hist_cut;
    if (name_cut.IsNull()) {}
    else {
      auto cut = ejungwoo::cutg(TString("data__region/")+name_cut+".region.root",name_cut,"p","dedx");
      hist_cut = (TH2D *) ejungwoo::cutg(hist,cut);
    }

    auto graph_fitproj = ejungwoo::new_ge(TString("graph_pdedx_")+name_cut);

    TF1 *fit_dedx = nullptr;
    if (0) {
           if (idx== 0) { fit_dedx = new TF1(TString("fit_dedx_")+name_cut, stp::bb_dedx_pin, p_hist_min, p_hist_max, 2); fit_dedx -> SetParameters(0.730379, -51.818967);}
      else if (idx== 1) { fit_dedx = new TF1(TString("fit_dedx_")+name_cut, stp::bb_dedx_pip, p_hist_min, p_hist_max, 2); fit_dedx -> SetParameters(0.945242, -42.6304);  }
      else if (idx== 2) { fit_dedx = new TF1(TString("fit_dedx_")+name_cut, stp::bb_dedx_p  , p_hist_min, p_hist_max, 2); fit_dedx -> SetParameters(0.945242, -42.6304);  }
      else if (idx== 3) { fit_dedx = new TF1(TString("fit_dedx_")+name_cut, stp::bb_dedx_d  , p_hist_min, p_hist_max, 2); fit_dedx -> SetParameters(0.758891, -54.0868);  }
      else if (idx== 4) { fit_dedx = new TF1(TString("fit_dedx_")+name_cut, stp::bb_dedx_t  , p_hist_min, p_hist_max, 2); fit_dedx -> SetParameters(1.37889, -29.8536);   }
      else if (idx== 5) { fit_dedx = new TF1(TString("fit_dedx_")+name_cut, stp::bb_dedx_he3, p_hist_min, p_hist_max, 2); fit_dedx -> SetParameters(0.867347, -54.3779);  }
      else if (idx== 6) { fit_dedx = new TF1(TString("fit_dedx_")+name_cut, stp::bb_dedx_he4, p_hist_min, p_hist_max, 2); fit_dedx -> SetParameters(0.19517, -247.207);   }
      else if (idx== 7) { fit_dedx = new TF1(TString("fit_dedx_")+name_cut, stp::bb_dedx_he6, p_hist_min, p_hist_max, 2); fit_dedx -> SetParameters(9.46557, -5.98054);   }
      else if (idx== 8) { fit_dedx = new TF1(TString("fit_dedx_")+name_cut, stp::bb_dedx_li6, p_hist_min, p_hist_max, 2); fit_dedx -> SetParameters(0.00328109, -15279.4); }
      else if (idx== 9) { fit_dedx = new TF1(TString("fit_dedx_")+name_cut, stp::bb_dedx_li7, p_hist_min, p_hist_max, 2); fit_dedx -> SetParameters(0.00471953, -11279.9); }
      else if (idx==10) { fit_dedx = new TF1(TString("fit_dedx_")+name_cut, stp::bb_dedx_ele, p_hist_min, p_hist_max, 2); fit_dedx -> SetParameters(0.945242, -42.6304);  }
      else if (idx==11) { fit_dedx = new TF1(TString("fit_dedx_")+name_cut, stp::bb_dedx_pos, p_hist_min, p_hist_max, 2); fit_dedx -> SetParameters(0.945242, -42.6304);  }
    } else {
           if (idx== 0) { fit_dedx = new TF1(TString("fit_dedx_")+name_cut, stp::expression_bb_dedx_pin(), p_hist_min, p_hist_max); /*cout << stp::expression_bb_dedx_pin() << endl;*/ fit_dedx -> SetParameters(0.730379, -51.818967);}
      else if (idx== 1) { fit_dedx = new TF1(TString("fit_dedx_")+name_cut, stp::expression_bb_dedx_pip(), p_hist_min, p_hist_max); /*cout << stp::expression_bb_dedx_pip() << endl;*/ fit_dedx -> SetParameters(0.945242, -42.6304);  }
      else if (idx== 2) { fit_dedx = new TF1(TString("fit_dedx_")+name_cut, stp::expression_bb_dedx_p()  , p_hist_min, p_hist_max); /*cout << stp::expression_bb_dedx_p()   << endl;*/ fit_dedx -> SetParameters(0.945242, -42.6304);  }
      else if (idx== 3) { fit_dedx = new TF1(TString("fit_dedx_")+name_cut, stp::expression_bb_dedx_d()  , p_hist_min, p_hist_max); /*cout << stp::expression_bb_dedx_d()   << endl;*/ fit_dedx -> SetParameters(0.758891, -54.0868);  }
      else if (idx== 4) { fit_dedx = new TF1(TString("fit_dedx_")+name_cut, stp::expression_bb_dedx_t()  , p_hist_min, p_hist_max); /*cout << stp::expression_bb_dedx_t()   << endl;*/ fit_dedx -> SetParameters(1.37889, -29.8536);   }
      else if (idx== 5) { fit_dedx = new TF1(TString("fit_dedx_")+name_cut, stp::expression_bb_dedx_he3(), p_hist_min, p_hist_max); /*cout << stp::expression_bb_dedx_he3() << endl;*/ fit_dedx -> SetParameters(0.867347, -54.3779);  }
      else if (idx== 6) { fit_dedx = new TF1(TString("fit_dedx_")+name_cut, stp::expression_bb_dedx_he4(), p_hist_min, p_hist_max); /*cout << stp::expression_bb_dedx_he4() << endl;*/ fit_dedx -> SetParameters(0.19517, -247.207);   }
      else if (idx== 7) { fit_dedx = new TF1(TString("fit_dedx_")+name_cut, stp::expression_bb_dedx_he6(), p_hist_min, p_hist_max); /*cout << stp::expression_bb_dedx_he6() << endl;*/ fit_dedx -> SetParameters(9.46557, -5.98054);   }
      else if (idx== 8) { fit_dedx = new TF1(TString("fit_dedx_")+name_cut, stp::expression_bb_dedx_li6(), p_hist_min, p_hist_max); /*cout << stp::expression_bb_dedx_li6() << endl;*/ fit_dedx -> SetParameters(0.00328109, -15279.4); }
      else if (idx== 9) { fit_dedx = new TF1(TString("fit_dedx_")+name_cut, stp::expression_bb_dedx_li7(), p_hist_min, p_hist_max); /*cout << stp::expression_bb_dedx_li7() << endl;*/ fit_dedx -> SetParameters(0.00471953, -11279.9); }
      else if (idx==10) { fit_dedx = new TF1(TString("fit_dedx_")+name_cut, stp::expression_bb_dedx_ele(), p_hist_min, p_hist_max); /*cout << stp::expression_bb_dedx_ele() << endl;*/ fit_dedx -> SetParameters(0.945242, -42.6304);  }
      else if (idx==11) { fit_dedx = new TF1(TString("fit_dedx_")+name_cut, stp::expression_bb_dedx_pos(), p_hist_min, p_hist_max); /*cout << stp::expression_bb_dedx_pos() << endl;*/ fit_dedx -> SetParameters(0.945242, -42.6304);  }
      else if (idx==12) { fit_dedx = new TF1(TString("fit_dedx_")+name_cut, stp::expression_bb_dedx_li (), p_hist_min, p_hist_max); /*cout << stp::expression_bb_dedx_li () << endl;*/ fit_dedx -> SetParameters(0.945242, -42.6304);  }
    }

    TF1 *fit_dedx_par5 = nullptr;
    {
           if (idx== 0) { fit_dedx_par5 = new TF1(TString("fit_dedx_par5_")+name_cut, stp::p5_dedx_pin, p_hist_min, p_hist_max, 5); }
      else if (idx== 1) { fit_dedx_par5 = new TF1(TString("fit_dedx_par5_")+name_cut, stp::p5_dedx_pip, p_hist_min, p_hist_max, 5); fit_dedx_par5 -> SetParameters(-0.838172, -17.197, 2.74259, 1.1186, -2.03426);   }
      else if (idx== 2) { fit_dedx_par5 = new TF1(TString("fit_dedx_par5_")+name_cut, stp::p5_dedx_p  , p_hist_min, p_hist_max, 5); fit_dedx_par5 -> SetParameters(-5.82913, 3.4971, 488.69, 0.972273, 1.48717);     }
      else if (idx== 3) { fit_dedx_par5 = new TF1(TString("fit_dedx_par5_")+name_cut, stp::p5_dedx_d  , p_hist_min, p_hist_max, 5); fit_dedx_par5 -> SetParameters(-5.81648, 8.19503, 65740.6, 0.975845, 2.58973);   }
      else if (idx== 4) { fit_dedx_par5 = new TF1(TString("fit_dedx_par5_")+name_cut, stp::p5_dedx_t  , p_hist_min, p_hist_max, 5); fit_dedx_par5 -> SetParameters(-0.561179, -24.1703, 22146.2, 0.963454, 2.53092); }
      else if (idx== 5) { fit_dedx_par5 = new TF1(TString("fit_dedx_par5_")+name_cut, stp::p5_dedx_he3, p_hist_min, p_hist_max, 5); fit_dedx_par5 -> SetParameters(-20.7459, 4.11136, 3845.17, 0.902391, 3.60306);   }
      else if (idx== 6) { fit_dedx_par5 = new TF1(TString("fit_dedx_par5_")+name_cut, stp::p5_dedx_he4, p_hist_min, p_hist_max, 5); fit_dedx_par5 -> SetParameters(-15.2432, 4.09085, 26656.4, 0.863061, 4.75901);   }
      else if (idx== 7) { fit_dedx_par5 = new TF1(TString("fit_dedx_par5_")+name_cut, stp::p5_dedx_he6, p_hist_min, p_hist_max, 5); fit_dedx_par5 -> SetParameters(-1.53396, -57.5896, 20306.2, 0.816746, -17.5053); }
      else if (idx== 8) { fit_dedx_par5 = new TF1(TString("fit_dedx_par5_")+name_cut, stp::p5_dedx_li6, p_hist_min, p_hist_max, 5); fit_dedx_par5 -> SetParameters(-1.53396, -57.5896, 20306.2, 0.816746, -17.5053); }
      else if (idx== 9) { fit_dedx_par5 = new TF1(TString("fit_dedx_par5_")+name_cut, stp::p5_dedx_li7, p_hist_min, p_hist_max, 5); fit_dedx_par5 -> SetParameters(-1.53396, -57.5896, 20306.2, 0.816746, -17.5053); }
      else if (idx==10) { fit_dedx_par5 = new TF1(TString("fit_dedx_par5_")+name_cut, stp::p5_dedx_ele, p_hist_min, p_hist_max, 5); }
      else if (idx==11) { fit_dedx_par5 = new TF1(TString("fit_dedx_par5_")+name_cut, stp::p5_dedx_pos, p_hist_min, p_hist_max, 5); }
      else if (idx==12) { fit_dedx_par5 = new TF1(TString("fit_dedx_par5_")+name_cut, stp::p5_dedx_li , p_hist_min, p_hist_max, 5); fit_dedx_par5 -> SetParameters(-1.53396, -57.5896, 20306.2, 0.816746, -17.5053); }

      fit_dedx_par5 -> SetLineColor(kRed);
      //fit_dedx_par5 -> SetParameters(1,1,1,1,1);
    }


    if (fit_projection)
    {
      double s_max = 2000.;
      double s_min = limits[idx][0];
      double p_max = 2000;
      double p_min = limits[idx][1];

      auto axis_p = hist_cut -> GetXaxis();
      auto axis_s = hist_cut -> GetYaxis();

      auto bin_s_min = axis_s -> FindBin(s_min);
      auto bin_s_max = axis_s -> FindBin(s_max);

      if (s_min < s_hist_max)
      for (auto bin_s = bin_s_max; bin_s >= bin_s_min; bin_s=bin_s-sum_bin_s)
      {
        auto bin_s_low = bin_s - sum_bin_s + 1;
        auto bin_s_up = bin_s;

        auto s_bin_up = axis_s -> GetBinLowEdge(bin_s_up);
        auto s_bin_low = axis_s -> GetBinUpEdge(bin_s_low);
        auto dbin_s = s_bin_up - s_bin_low;
        auto s_center = (s_bin_up + s_bin_low)/2.;

        auto hist_proj = (TH1D *) hist_cut -> ProjectionX(Form("bin_s%d_x%.2f",bin_s,s_center),bin_s_low,bin_s,"o");
        if (hist_proj -> Integral() < 1000)
          continue;

        hist_proj -> SetLineColor(kBlack);
        auto fit = ejungwoo::fitg(hist_proj,1);

        auto mean = fit -> GetParameter(1);
        auto sigma = fit -> GetParameter(2);
        auto hist_range1 = mean-5*sigma; if (hist_range1 < 0) hist_range1 = 0;
        auto hist_range2 = mean+5*sigma; if (hist_range2 > p_max) hist_range2 = p_max;
        hist_proj -> GetXaxis() -> SetRangeUser(hist_range1, hist_range2);

        auto idxPoint = graph_fitproj -> GetN();
        graph_fitproj -> SetPoint(idxPoint, mean, s_center);
        graph_fitproj -> SetPointError(idxPoint, sigma, dbin_s/2.);

        if (draw_projection_fits) {
          ejungwoo::cv();
          hist_proj -> Draw("hist");
          fit -> Draw("samel");
        }
      }

      auto bin_p_min = axis_p -> FindBin(p_min);
      auto bin_p_max = axis_p -> FindBin(p_max);
      for (auto bin_p = bin_p_min; bin_p <= bin_p_max; bin_p=bin_p+sum_bin_p) {
        auto bin_p_up = bin_p + sum_bin_p - 1;
        auto bin_p_low = bin_p;

        auto p_bin_low = axis_p -> GetBinLowEdge(bin_p_low);
        auto p_bin_up = axis_p -> GetBinUpEdge(bin_p_up);
        auto dbin_p = p_bin_up - p_bin_low;
        auto p_center = (p_bin_low + p_bin_up)/2.;

        auto hist_proj = (TH1D *) hist_cut -> ProjectionY(Form("bin_p%d_x%.2f",bin_p,p_center),bin_p_low,bin_p_up,"o");
        if (hist_proj -> Integral() < 1000)
          continue;

        hist_proj -> SetLineColor(kBlack);
        auto fit = ejungwoo::fitg(hist_proj);

        auto mean = fit -> GetParameter(1);
        auto sigma = fit -> GetParameter(2);
        auto hist_range1 = mean-5*sigma; if (hist_range1 < 0) hist_range1 = 0;
        auto hist_range2 = mean+5*sigma; if (hist_range2 > p_max) hist_range2 = p_max;
        hist_proj -> GetXaxis() -> SetRangeUser(hist_range1, hist_range2);

        auto idxPoint = graph_fitproj -> GetN();
        graph_fitproj -> SetPoint(graph_fitproj->GetN(), p_center, mean);
        graph_fitproj -> SetPointError(graph_fitproj->GetN()-1, dbin_p/2., sigma);
      }

      graph_fitproj -> Fit(fit_dedx,"Q0");
      ejungwoo::write(graph_fitproj);
      ejungwoo::write(fit_dedx);


      graph_fitproj -> Fit(fit_dedx_par5,"Q0");
      ejungwoo::write(fit_dedx_par5);

      cout << name_cut << 
      " (" << fit_dedx_par5 -> GetParameter(0) <<
      ", " << fit_dedx_par5 -> GetParameter(1) <<
      ", " << fit_dedx_par5 -> GetParameter(2) <<
      ", " << fit_dedx_par5 -> GetParameter(3) <<
      ", " << fit_dedx_par5 -> GetParameter(4) <<
      ");" << endl;
    }

    TCanvas *cvs;
    if (draw_pid) {
      cvs = ejungwoo::cc3();
      cvs -> SetLogz();
      hist_cut -> Draw("colz");
      graph_fitproj -> Draw("1samepz");
      fit_dedx -> SetLineColor(kGray+1);
      fit_dedx -> SetNpx(2000);
      fit_dedx -> Draw("samel");
      fit_dedx_par5 -> SetLineColor(kRed);
      fit_dedx_par5 -> SetNpx(2000);
      fit_dedx_par5 -> Draw("samel");
      ejungwoo::make_c(cvs);
    }
  }
}
