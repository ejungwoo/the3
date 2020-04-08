#include "stp.h"

double c_limit_p = .4;
double c_limit_oli = .15;
double c_limit_all = .10;

void get_par_limits(double p_current, int idx_particle, TH1D * hist, double &a, double &m, double &s, double &a1, double &a2, double &m1, double &m2, double &s1, double &s2);
void get_par_limits2(double s_current, int idx_particle, TH1D * hist, double &a, double &m, double &s, double &a1, double &a2, double &m1, double &m2, double &s1, double &s2);

void fit_sn108(bool verbose = 1, int pit = -1, int yaw = -1)
{
  int system = 108;
  ejungwoo::gzcolor(0);
  ejungwoo::gstat(0);
  ejungwoo::gsave(0);
  ejungwoo::gcvspos(1300,0);

  bool eeeeeeeee = true;

  bool draw_proj = false;
  bool draw_fit_status = false;

  double p_start = 100.;
  double p_end = 2000.;
  int bin_size_p = 5;

  double s_start = 100.;
  double s_end = 1300.;
  //double s_end = 988.;
  int bin_size_s = 10;

  int idx_particles[] = {
    stp::kpip,
    stp::kp,
    stp::kd,
    stp::kt,
    stp::khe3,
    stp::khe4,
    stp::khe6,
    stp::kli,
    };

  double p_hist_min = 0;
  double p_hist_max = 2500;
  int nbins_p = 500;
  int nbins_s = 500;
  double p_min = 0;
  double p_max = 2500;
  double s_min = 0;
  double s_max = 1800;
  double dbin_p = (p_max-p_min)/nbins_p;
  double dbin_s = (s_max-s_min)/nbins_s;

  double plimit1[] = {
  /*stp::kpin*/ p_hist_min, 
  /*stp::kpip*/ 100, 
  /*stp::kp  */ 350, 
  /*stp::kd  */ 400, 
  /*stp::kt  */ 600       , 
  /*stp::khe3*/ 600       , 
  /*stp::khe4*/ 700       , 
  /*stp::khe6*/ 1700      , 
  /*stp::kli6*/ 1000      , 
  /*stp::kli7*/ 1500      , 
  /*stp::kele*/ p_hist_min, 
  /*stp::kpos*/ p_hist_min,
  /*stp::kli */ 1000    , 
  };

  double plimit2[] = {
  /*stp::kpin*/ p_hist_max, 
  /*stp::kpip*/ 500       , 
  /*stp::kp  */ 1700   , 
  /*stp::kd  */ p_hist_max, 
  /*stp::kt  */ p_hist_max, 
  /*stp::khe3*/ 1400 , 
  /*stp::khe4*/ p_hist_max, 
  /*stp::khe6*/ p_hist_max, 
  /*stp::kli6*/ 2000 , 
  /*stp::kli7*/ 2300      , 
  /*stp::kele*/ p_hist_max, 
  /*stp::kpos*/ p_hist_max,
  /*stp::kli */ 2000 , 
  };

  double slimit1[] = {
  /*stp::kpin*/ 0,
  /*stp::kpip*/ 100,
  /*stp::kp  */ 1000, 
  /*stp::kd  */ 1000, 
  /*stp::kt  */ 1000, 
  /*stp::khe3*/ 1000, 
  /*stp::khe4*/ 1000, 
  /*stp::khe6*/ 0, 
  /*stp::kli6*/ 0, 
  /*stp::kli7*/ 0, 
  /*stp::kele*/ 0, 
  /*stp::kpos*/ 0,
  /*stp::kli */ 0, 
  };

  double slimit2[] = {
  /*stp::kpin*/ 0, 
  /*stp::kpip*/ 200, 
  /*stp::kp  */ 1500, 
  /*stp::kd  */ 1500, 
  /*stp::kt  */ 1500, 
  /*stp::khe3*/ 1500, 
  /*stp::khe4*/ 1500, 
  /*stp::khe6*/ 0, 
  /*stp::kli6*/ 0, 
  /*stp::kli7*/ 0, 
  /*stp::kele*/ 0, 
  /*stp::kpos*/ 0,
  /*stp::kli */ 0, 
  };

  TString nameVersion = Form("sys108_p%dy%d",pit,yaw);
  if (pit<0&&yaw<0) nameVersion = "sys108";
  ejungwoo::gversion(nameVersion);
  TString name_hist = TString("hist_") + ejungwoo::version();
  TString name_cvs = TString("cvs_") + ejungwoo::version();
  TString name_pid_file = "";

  TH2D *hist_pid = nullptr;
  TString name_ana = Form("hpid%d_nn%dx%d_p%.0f_dedx%.0f",system,nbins_p,nbins_s,p_max,s_max);
  if (!name_pid_file.IsNull()) {
    auto file_pid = new TFile(name_pid_file);
    hist_pid = (TH2D *) file_pid -> Get(name_hist);
    name_ana = name_hist;
  }
  else {
    name_pid_file = Form("data__%s/%s.%s.root",ejungwoo::version().Data(),name_ana.Data(),ejungwoo::version().Data());
    cout << name_pid_file << endl;
    auto file_pid = new TFile(name_pid_file);
    if (file_pid -> IsOpen()) {
      hist_pid = (TH2D *) file_pid -> Get(name_ana);
    } else {
      auto fileTree = new TFile(Form("data__in/pdedx%d.root",system));
      auto tree = (TTree *) fileTree -> Get("data");
      hist_pid = (TH2D *) ejungwoo::tp(Form("%s",name_ana.Data()),tree,"dedx:p","qrun","",nbins_p,p_min,p_max,nbins_s,s_min,s_max);
      ejungwoo::write(hist_pid);
    }
  }

  TCanvas *cvs_pid;
  if (1) {
    cvs_pid = ejungwoo::cc4();
    cvs_pid -> SetLogz();
    hist_pid -> Draw("colz");
  }

  double parameters[50] = {0};
  double parlimits1[50] = {0};
  double parlimits2[50] = {0};
  TGraphErrors *graph_dedx_onpid[3][20];
  TGraphErrors *graph_particle_ams[2][3][20];

  TF1 *fit_dedx_array[20] = {0};
  for (auto idx_particle : idx_particles)
  {
    TString name_particle = stp::fNameShort[idx_particle];
    graph_dedx_onpid[0][idx_particle] = ejungwoo::new_ge(TString("dedx_") + name_particle);
    graph_dedx_onpid[1][idx_particle] = ejungwoo::new_ge(TString("dedx_") + name_particle + "_1");
    graph_dedx_onpid[2][idx_particle] = ejungwoo::new_ge(TString("dedx_") + name_particle + "_2");

    graph_particle_ams[0][0][idx_particle] = ejungwoo::new_ge("amp_"  +name_particle+"_s");
    graph_particle_ams[0][1][idx_particle] = ejungwoo::new_ge("mean_" +name_particle+"_s");
    graph_particle_ams[0][2][idx_particle] = ejungwoo::new_ge("sigma_"+name_particle+"_s");
    graph_particle_ams[1][0][idx_particle] = ejungwoo::new_ge("amp_"  +name_particle);
    graph_particle_ams[1][1][idx_particle] = ejungwoo::new_ge("mean_" +name_particle);
    graph_particle_ams[1][2][idx_particle] = ejungwoo::new_ge("sigma_"+name_particle);

    auto file_fit = new TFile(TString("data__prefit/fit_dedx_par5_region_")+name_particle+".prefit.root");
    auto fit_par5 = (TF1 *) file_fit -> Get(TString("fit_dedx_par5_region_")+name_particle);
    TF1 *fit_dedx;
    auto name_par5 = TString("fpre_dedx_")+name_particle;
    /**/ if (idx_particle == stp::kpin) { fit_dedx = new TF1(name_par5, stp::p5_dedx_pin, p_hist_min, p_hist_max, 5); }
    else if (idx_particle == stp::kpip) { fit_dedx = new TF1(name_par5, stp::p5_dedx_pip, p_hist_min, p_hist_max, 5); }
    else if (idx_particle == stp::kp  ) { fit_dedx = new TF1(name_par5, stp::p5_dedx_p  , p_hist_min, p_hist_max, 5); }
    else if (idx_particle == stp::kd  ) { fit_dedx = new TF1(name_par5, stp::p5_dedx_d  , p_hist_min, p_hist_max, 5); }
    else if (idx_particle == stp::kt  ) { fit_dedx = new TF1(name_par5, stp::p5_dedx_t  , p_hist_min, p_hist_max, 5); }
    else if (idx_particle == stp::khe3) { fit_dedx = new TF1(name_par5, stp::p5_dedx_he3, p_hist_min, p_hist_max, 5); }
    else if (idx_particle == stp::khe4) { fit_dedx = new TF1(name_par5, stp::p5_dedx_he4, p_hist_min, p_hist_max, 5); }
    else if (idx_particle == stp::khe6) { fit_dedx = new TF1(name_par5, stp::p5_dedx_he6, p_hist_min, p_hist_max, 5); }
    else if (idx_particle == stp::kli6) { fit_dedx = new TF1(name_par5, stp::p5_dedx_li6, p_hist_min, p_hist_max, 5); }
    else if (idx_particle == stp::kli7) { fit_dedx = new TF1(name_par5, stp::p5_dedx_li7, p_hist_min, p_hist_max, 5); }
    else if (idx_particle == stp::kele) { fit_dedx = new TF1(name_par5, stp::p5_dedx_ele, p_hist_min, p_hist_max, 5); }
    else if (idx_particle == stp::kpos) { fit_dedx = new TF1(name_par5, stp::p5_dedx_pos, p_hist_min, p_hist_max, 5); }
    else if (idx_particle == stp::kli ) { fit_dedx = new TF1(name_par5, stp::p5_dedx_li , p_hist_min, p_hist_max, 5); }

    fit_dedx -> SetParameters(
        fit_par5->GetParameter(0),
        fit_par5->GetParameter(1),
        fit_par5->GetParameter(2),
        fit_par5->GetParameter(3),
        fit_par5->GetParameter(4));

    fit_dedx -> SetNpx(2000);
    fit_dedx -> SetLineColor(kGray);
    //fit_dedx -> Draw("samel");
    fit_dedx_array[idx_particle] = fit_dedx;
  }

  int bin_s_start = hist_pid -> GetYaxis() -> FindBin(s_start);
  int bin_s_end = hist_pid -> GetYaxis() -> FindBin(s_end);

  for (auto bin1 = bin_s_end; bin1 > bin_s_start; bin1-=bin_size_s)
  {
    bool idx_is_good[20] = {0};

    auto bin2 = bin1-bin_size_p+1;
    auto s_bin_up = hist_pid -> GetXaxis() -> GetBinLowEdge(bin1);
    auto s_bin_low = hist_pid -> GetXaxis() -> GetBinUpEdge(bin2);
    auto ds_bin = s_bin_up - s_bin_low;

    auto s_current1 = hist_pid -> GetYaxis() -> GetBinCenter(bin1);
    auto s_current2 = hist_pid -> GetYaxis() -> GetBinCenter(bin2);
    auto s_current = (s_current1 + s_current2)/2.;
    //TString name_s = Form("proj_bin_s%.2f_b%d",s_current,bin1);
    TString name_s = Form("proj_bin_s%.2f_b%d_%d",s_current,bin1,bin2);
    auto hist_proj = hist_pid -> ProjectionX(name_s,bin2,bin1);
    hist_proj -> SetTitle(name_s);

    TCanvas *cvs_proj = nullptr;

    double p_amp_array[20] = {0};
    double p_mean_array[20] = {0};
    double p_sigma_array[20] = {0};

    auto graph_init = ejungwoo::new_ge(name_s+"_inits");
    graph_init -> SetMarkerColor(kRed);
    graph_init -> SetLineColor(kRed);

    bool aaaaaaaaaaaaa = false;

    for (auto idx_particle : idx_particles)
    {
      auto fit = fit_dedx_array[idx_particle];
      double srange1 = slimit1[idx_particle];
      double srange2 = slimit2[idx_particle];
      if (s_current > srange1 && s_current < srange2)
      {
        idx_is_good[idx_particle] = true;
        aaaaaaaaaaaaa = true;
        auto p_mean_init = fit -> GetX(s_current);

        p_amp_array[idx_particle] = hist_proj -> GetBinContent(hist_proj->GetXaxis()->FindBin(p_mean_init));
        p_mean_array[idx_particle] = p_mean_init;
      }
    }

    if (!aaaaaaaaaaaaa)
      continue;

    if (draw_proj) {
      cvs_proj = ejungwoo::cv(name_s);
      cvs_proj -> SetLogy();
      hist_proj -> Draw();
    }

    TString expression = "0";
    int count_s = 0;
    for (auto idx_particle : idx_particles)
    {
      TString name_particle = stp::fNameShort[idx_particle];

      if (!idx_is_good[idx_particle])
        continue;

      double p_mean_init = p_mean_array[idx_particle];
      if (p_mean_init <= 0)
        continue;

      double p_amp_init = hist_proj -> GetBinContent(hist_proj->GetXaxis()->FindBin(p_mean_init));

      if (p_amp_init < 2.) {
        idx_is_good[idx_particle] = false;
        continue;
      }

      expression = expression + Form("+gaus(%d)",3*count_s);

      double p_sigma_init = (idx_particle+1)*2;

      double a1, a2, m1, m2, s1, s2;

      get_par_limits2(s_current, idx_particle, hist_proj, p_amp_init, p_mean_init, p_sigma_init, a1, a2, m1, m2, s1, s2);

      parameters[3*idx_particle+0] = p_amp_init;
      parameters[3*idx_particle+1] = p_mean_init;
      parameters[3*idx_particle+2] = p_sigma_init;
      parlimits1[3*idx_particle+0] = a1;
      parlimits2[3*idx_particle+0] = a2;
      parlimits1[3*idx_particle+1] = m1;
      parlimits2[3*idx_particle+1] = m2;
      parlimits1[3*idx_particle+2] = s1;
      parlimits2[3*idx_particle+2] = s2;

      if (draw_proj && draw_fit_status) {
        cvs_proj -> cd();
        graph_init -> SetPoint(idx_particle, p_mean_init,100);
        graph_init -> SetPointError(idx_particle, (m2-m1)/2., 0);
      }

      count_s++;
    }

    if (draw_proj && draw_fit_status) {
      cvs_proj -> cd();
      graph_init -> Sort();
      graph_init -> Draw("samepz1");
    }

    int count_particle = 0;
    TF1 *fit_p_total = new TF1("fit_p_total",expression,0,1000);
    for (auto idx_particle : idx_particles)
    {
      if (!idx_is_good[idx_particle])
        continue;

      TString name_particle = stp::fNameShort[idx_particle];

      fit_p_total -> SetParameter(3*count_particle+0,parameters[3*idx_particle+0]);
      fit_p_total -> SetParameter(3*count_particle+1,parameters[3*idx_particle+1]);
      fit_p_total -> SetParameter(3*count_particle+2,parameters[3*idx_particle+2]);

      fit_p_total -> SetParLimits(3*count_particle+0,parlimits1[3*idx_particle+0],parlimits2[3*idx_particle+0]);
      fit_p_total -> SetParLimits(3*count_particle+1,parlimits1[3*idx_particle+1],parlimits2[3*idx_particle+1]);
      fit_p_total -> SetParLimits(3*count_particle+2,parlimits1[3*idx_particle+2],parlimits2[3*idx_particle+2]);

      count_particle++;
    }

    hist_proj -> Fit(fit_p_total,"Q0");

    count_particle = 0;
    auto legend = new TLegend();
    for (auto idx_particle : idx_particles)
    {
      if (!idx_is_good[idx_particle])
        continue;

      TString name_particle = stp::fNameShort[idx_particle];

      auto p_amp = fit_p_total -> GetParameter(3*count_particle+0);
      auto p_mean = fit_p_total -> GetParameter(3*count_particle+1);
      auto p_sigma = fit_p_total -> GetParameter(3*count_particle+2);

      auto marker = new TMarker(p_mean,p_amp,30);
      count_particle++;
      marker -> SetMarkerSize(1.5);
      marker -> SetMarkerColor(kBlack);
      if (draw_proj && draw_fit_status) {
        cvs_proj -> cd();
        marker -> Draw("samep");
        auto fit_p_single = new TF1(Form("fit_p_single_%s",name_particle.Data()),"gaus(0)",0,1000);
        fit_p_single -> SetParameters(p_amp,p_mean,p_sigma);
        fit_p_single -> SetLineColor(kGray);
        fit_p_single -> Draw("samel");

        legend -> AddEntry(fit_p_single,Form("%s: %.2e, %.2e, %.2e",name_particle.Data(), p_amp,p_mean,p_sigma),"l");
      }

      auto idxn = graph_dedx_onpid[0][idx_particle] -> GetN();
      graph_dedx_onpid[0][idx_particle] -> SetPoint(idxn, p_mean, s_current);
      graph_dedx_onpid[0][idx_particle] -> SetPointError(idxn, p_sigma, ds_bin/2.);
      graph_dedx_onpid[1][idx_particle] -> SetPoint(idxn, p_mean-p_sigma, s_current);
      graph_dedx_onpid[2][idx_particle] -> SetPoint(idxn, p_mean+p_sigma, s_current);

      graph_particle_ams[0][0][idx_particle] -> SetPoint(graph_particle_ams[0][0][idx_particle] -> GetN(), s_current, p_amp);
      graph_particle_ams[0][1][idx_particle] -> SetPoint(graph_particle_ams[0][1][idx_particle] -> GetN(), s_current, p_mean);
      graph_particle_ams[0][2][idx_particle] -> SetPoint(graph_particle_ams[0][2][idx_particle] -> GetN(), s_current, p_sigma);
    }

    if (draw_proj && draw_fit_status) {
      cvs_proj -> cd();
      fit_p_total -> SetNpx(2000);
      fit_p_total -> SetLineColor(kGray+2);
      fit_p_total -> Draw("samel");
      ejungwoo::make_l(cvs_proj,legend,0,0,0.6,0.4) -> Draw("same");

      auto npar = fit_p_total -> GetNpar();
      auto p_mean_last = fit_p_total -> GetParameter(npar-2);
      auto p_sigma_last = fit_p_total -> GetParameter(npar-1);

      hist_proj -> GetXaxis() -> SetRangeUser(0,p_mean_last+2*p_sigma_last);

      //auto hist_sub = ejungwoo::subtract(hist_proj, fit_p_total);
      //hist_sub -> SetLineColor(kGray+2);
      //hist_sub -> Draw("same hist");
    }
  }

  int bin_p_start = hist_pid -> GetXaxis() -> FindBin(p_start);
  int bin_p_end = hist_pid -> GetXaxis() -> FindBin(p_end);

  if (eeeeeeeee)
  for (auto bin1 = bin_p_start; bin1 < bin_p_end; bin1+=bin_size_p)
  {
    bool idx_is_good[20] = {0};

    auto bin2 = bin1+bin_size_p-1;
    auto p_bin_low = hist_pid -> GetXaxis() -> GetBinLowEdge(bin1);
    auto p_bin_up = hist_pid -> GetXaxis() -> GetBinUpEdge(bin2);
    auto dp_bin = p_bin_up - p_bin_low;
    auto p_current1 = hist_pid -> GetXaxis() -> GetBinCenter(bin1);
    auto p_current2 = hist_pid -> GetXaxis() -> GetBinCenter(bin2);
    auto p_current = (p_current1 + p_current2)/2.;
    TString name_p = Form("proj_bin_p%.2f_b%d_%d",p_current,bin1,bin2);
    auto hist_proj = hist_pid -> ProjectionY(name_p,bin1,bin2);
    hist_proj -> SetTitle(name_p);

    TCanvas *cvs_proj = nullptr;
    if (draw_proj) {
      cvs_proj = ejungwoo::cv(name_p);
      cvs_proj -> SetLogy();
      hist_proj -> Draw();
    }

    double s_amp_array[20] = {0};
    double s_mean_array[20] = {0};
    double s_sigma_array[20] = {0};

    auto graph_init = ejungwoo::new_ge(name_p+"_inits");
    graph_init -> SetMarkerColor(kRed);
    graph_init -> SetLineColor(kRed);

    for (auto idx_particle : idx_particles)
    {
      auto fit = fit_dedx_array[idx_particle];
      double prange1 = plimit1[idx_particle];
      double prange2 = plimit2[idx_particle];
      if (p_current > prange1 && p_current < prange2)
      {
        idx_is_good[idx_particle] = true;
        auto s_mean_init = fit -> Eval(p_current);

        s_amp_array[idx_particle] = hist_proj -> GetBinContent(hist_proj->GetXaxis()->FindBin(s_mean_init));
        s_mean_array[idx_particle] = s_mean_init;
      }
    }

    TString expression = "0";
    int count_s = 0;
    for (auto idx_particle : idx_particles)
    {
      TString name_particle = stp::fNameShort[idx_particle];

      if (!idx_is_good[idx_particle])
        continue;

      double s_mean_init = s_mean_array[idx_particle];
      if (s_mean_init <= 0)
        continue;

      double s_amp_init = hist_proj -> GetBinContent(hist_proj->GetXaxis()->FindBin(s_mean_init));

      if (s_amp_init < 2.) {
        idx_is_good[idx_particle] = false;
        continue;
      }

      expression = expression + Form("+gaus(%d)",3*count_s);

      double s_sigma_init = (idx_particle+1)*2;

      double a1, a2, m1, m2, s1, s2;

      get_par_limits(p_current, idx_particle, hist_proj, s_amp_init, s_mean_init, s_sigma_init, a1, a2, m1, m2, s1, s2);

      parameters[3*idx_particle+0] = s_amp_init;
      parameters[3*idx_particle+1] = s_mean_init;
      parameters[3*idx_particle+2] = s_sigma_init;
      parlimits1[3*idx_particle+0] = a1;
      parlimits2[3*idx_particle+0] = a2;
      parlimits1[3*idx_particle+1] = m1;
      parlimits2[3*idx_particle+1] = m2;
      parlimits1[3*idx_particle+2] = s1;
      parlimits2[3*idx_particle+2] = s2;

      if (draw_proj && draw_fit_status) {
        cvs_proj -> cd();
        graph_init -> SetPoint(idx_particle, s_mean_init,100);
        graph_init -> SetPointError(idx_particle, (m2-m1)/2., 0);
      }

      count_s++;
    }

    if (draw_proj && draw_fit_status) {
      cvs_proj -> cd();
      graph_init -> Sort();
      graph_init -> Draw("samepz1");
    }

    int count_particle = 0;
    TF1 *fit_s_total = new TF1("fit_s_total",expression,0,1000);
    for (auto idx_particle : idx_particles)
    {
      if (!idx_is_good[idx_particle])
        continue;

      TString name_particle = stp::fNameShort[idx_particle];

      fit_s_total -> SetParameter(3*count_particle+0,parameters[3*idx_particle+0]);
      fit_s_total -> SetParameter(3*count_particle+1,parameters[3*idx_particle+1]);
      fit_s_total -> SetParameter(3*count_particle+2,parameters[3*idx_particle+2]);

      fit_s_total -> SetParLimits(3*count_particle+0,parlimits1[3*idx_particle+0],parlimits2[3*idx_particle+0]);
      fit_s_total -> SetParLimits(3*count_particle+1,parlimits1[3*idx_particle+1],parlimits2[3*idx_particle+1]);
      fit_s_total -> SetParLimits(3*count_particle+2,parlimits1[3*idx_particle+2],parlimits2[3*idx_particle+2]);

      count_particle++;
    }

    hist_proj -> Fit(fit_s_total,"Q0");

    count_particle = 0;
    auto legend = new TLegend();
    for (auto idx_particle : idx_particles)
    {
      if (!idx_is_good[idx_particle])
        continue;

      TString name_particle = stp::fNameShort[idx_particle];

      auto s_amp = fit_s_total -> GetParameter(3*count_particle+0);
      auto s_mean = fit_s_total -> GetParameter(3*count_particle+1);
      auto s_sigma = fit_s_total -> GetParameter(3*count_particle+2);

      auto marker = new TMarker(s_mean,s_amp,30);
      count_particle++;
      marker -> SetMarkerSize(1.5);
      marker -> SetMarkerColor(kBlack);
      if (draw_proj && draw_fit_status) {
        cvs_proj -> cd();
        marker -> Draw("samep");
        auto fit_s_single = new TF1(Form("fit_s_single_%s",name_particle.Data()),"gaus(0)",0,1000);
        fit_s_single -> SetParameters(s_amp,s_mean,s_sigma);
        fit_s_single -> SetLineColor(kGray);
        fit_s_single -> Draw("samel");

        legend -> AddEntry(fit_s_single,Form("%s: %.2e, %.2e, %.2e",name_particle.Data(), s_amp,s_mean,s_sigma),"l");
      }

      auto idxn = graph_dedx_onpid[0][idx_particle] -> GetN();
      graph_dedx_onpid[0][idx_particle] -> SetPoint(idxn, p_current, s_mean);
      graph_dedx_onpid[0][idx_particle] -> SetPointError(idxn, dp_bin/2., s_sigma);
      graph_dedx_onpid[1][idx_particle] -> SetPoint(idxn, p_current, s_mean-s_sigma);
      graph_dedx_onpid[2][idx_particle] -> SetPoint(idxn, p_current, s_mean+s_sigma);

      graph_particle_ams[1][0][idx_particle] -> SetPoint(graph_particle_ams[1][0][idx_particle] -> GetN(), p_current, s_amp);
      graph_particle_ams[1][1][idx_particle] -> SetPoint(graph_particle_ams[1][1][idx_particle] -> GetN(), p_current, s_mean);
      graph_particle_ams[1][2][idx_particle] -> SetPoint(graph_particle_ams[1][2][idx_particle] -> GetN(), p_current, s_sigma);
    }

    if (draw_proj && draw_fit_status) {
      cvs_proj -> cd();
      fit_s_total -> SetNpx(2000);
      fit_s_total -> SetLineColor(kGray+2);
      fit_s_total -> Draw("samel");
      ejungwoo::make_l(cvs_proj,legend,0,0,0.6,0.4) -> Draw("same");

      auto npar = fit_s_total -> GetNpar();
      auto s_mean_last = fit_s_total -> GetParameter(npar-2);
      auto s_sigma_last = fit_s_total -> GetParameter(npar-1);

      hist_proj -> GetXaxis() -> SetRangeUser(0,s_mean_last+2*s_sigma_last);

      //auto hist_sub = ejungwoo::subtract(hist_proj, fit_s_total);
      //hist_sub -> SetLineColor(kGray+2);
      //hist_sub -> Draw("same hist");
    }
  }


  auto name_summary = ejungwoo::name_full("summary");
  auto file_summary = new TFile(name_summary,"recreate");
  vector<TGraphErrors *> gglist;
  vector<TF1 *> fflist;

  for (auto idx_particle : idx_particles)
  {
    TString name_particle = stp::fNameShort[idx_particle];
    cout << name_particle << endl;

    TF1 *fit_dedx_par5[3];
    auto name_par5 = TString("fit_dedx_")+name_particle+"_";
         if (idx_particle == stp::kpin) { for (auto ifit : {0,1,2}) fit_dedx_par5[ifit] = new TF1(name_par5+ifit, stp::p5_dedx_pin, p_hist_min, p_hist_max, 5); }
    else if (idx_particle == stp::kpip) { for (auto ifit : {0,1,2}) fit_dedx_par5[ifit] = new TF1(name_par5+ifit, stp::p5_dedx_pip, p_hist_min, p_hist_max, 5); }
    else if (idx_particle == stp::kp  ) { for (auto ifit : {0,1,2}) fit_dedx_par5[ifit] = new TF1(name_par5+ifit, stp::p5_dedx_p  , p_hist_min, p_hist_max, 5); }
    else if (idx_particle == stp::kd  ) { for (auto ifit : {0,1,2}) fit_dedx_par5[ifit] = new TF1(name_par5+ifit, stp::p5_dedx_d  , p_hist_min, p_hist_max, 5); }
    else if (idx_particle == stp::kt  ) { for (auto ifit : {0,1,2}) fit_dedx_par5[ifit] = new TF1(name_par5+ifit, stp::p5_dedx_t  , p_hist_min, p_hist_max, 5); }
    else if (idx_particle == stp::khe3) { for (auto ifit : {0,1,2}) fit_dedx_par5[ifit] = new TF1(name_par5+ifit, stp::p5_dedx_he3, p_hist_min, p_hist_max, 5); }
    else if (idx_particle == stp::khe4) { for (auto ifit : {0,1,2}) fit_dedx_par5[ifit] = new TF1(name_par5+ifit, stp::p5_dedx_he4, p_hist_min, p_hist_max, 5); }
    else if (idx_particle == stp::khe6) { for (auto ifit : {0,1,2}) fit_dedx_par5[ifit] = new TF1(name_par5+ifit, stp::p5_dedx_he6, p_hist_min, p_hist_max, 5); }
    else if (idx_particle == stp::kli6) { for (auto ifit : {0,1,2}) fit_dedx_par5[ifit] = new TF1(name_par5+ifit, stp::p5_dedx_li6, p_hist_min, p_hist_max, 5); }
    else if (idx_particle == stp::kli7) { for (auto ifit : {0,1,2}) fit_dedx_par5[ifit] = new TF1(name_par5+ifit, stp::p5_dedx_li7, p_hist_min, p_hist_max, 5); }
    else if (idx_particle == stp::kele) { for (auto ifit : {0,1,2}) fit_dedx_par5[ifit] = new TF1(name_par5+ifit, stp::p5_dedx_pos, p_hist_min, p_hist_max, 5); }
    else if (idx_particle == stp::kpos) { for (auto ifit : {0,1,2}) fit_dedx_par5[ifit] = new TF1(name_par5+ifit, stp::p5_dedx_li , p_hist_min, p_hist_max, 5); }
    else if (idx_particle == stp::kli ) { for (auto ifit : {0,1,2}) fit_dedx_par5[ifit] = new TF1(name_par5+ifit, stp::p5_dedx_ele, p_hist_min, p_hist_max, 5); }

    for (auto ifit : {0,1,2})
    {
      auto gg = graph_dedx_onpid[ifit][idx_particle];
      if (gg -> GetN()<2)
        continue;

      if (ifit!=0) {
        gg -> SetMarkerStyle(25);
        gg -> SetMarkerSize(0.2);
      }

      cvs_pid -> cd();
      //gg -> SetLineColor(kOrange-2);
      //gg -> Draw("samepz1");

      auto ff = fit_dedx_par5[ifit];
      if (ifit==0) {
        auto fx = fit_dedx_array[idx_particle];
        ff -> SetParameters(fx->GetParameter(0), fx->GetParameter(1), fx->GetParameter(2), fx->GetParameter(3), fx->GetParameter(4));
        ff -> SetLineColor(kBlack);
      } else {
        auto f0 = fit_dedx_par5[0];
        ff -> SetParameters(f0->GetParameter(0), f0->GetParameter(1), f0->GetParameter(2), f0->GetParameter(3), f0->GetParameter(4));
        ff -> SetLineColor(kOrange+8);
      }

      gg -> Fit(ff,"Q0");

      if (idx_particle == stp::kp || idx_particle == stp::kd || idx_particle == stp::kt || idx_particle == stp::khe3 || idx_particle == stp::khe4) {
        auto x1 = ff -> GetX(1500.);
        ff -> SetRange(x1, ejungwoo::x2_g(gg)+500);
      } else
        ff -> SetRange(ejungwoo::x1_g(gg)-10, ejungwoo::x2_g(gg)+10);

      ff -> SetNpx(2000);
      if (ifit !=0) ff -> Draw("samel");

      gglist.push_back(gg);
      fflist.push_back(ff);
    }
  }

  for (auto idx_particle : idx_particles) {
    ejungwoo::addto("amp0"  ,graph_particle_ams[0][0][idx_particle]);
    ejungwoo::addto("mean0" ,graph_particle_ams[0][1][idx_particle]);
    ejungwoo::addto("sigma0",graph_particle_ams[0][2][idx_particle]);
    ejungwoo::addto("amp1"  ,graph_particle_ams[1][0][idx_particle]);
    ejungwoo::addto("mean1" ,graph_particle_ams[1][1][idx_particle]);
    ejungwoo::addto("sigma1",graph_particle_ams[1][2][idx_particle]);
  }

  ejungwoo::drawc("amp0");
  ejungwoo::drawc("mean0");
  ejungwoo::drawc("sigma0");
  ejungwoo::drawc("mean1");
  ejungwoo::drawc("sigma1");

  auto cvs_amp = ejungwoo::drawc("amp1");
  for (auto idx_particle : idx_particles)
  {
    cvs_amp -> cd();
    auto gg = graph_particle_ams[1][0][idx_particle];
    auto hh = ejungwoo::tohist(gg,0,2500);
    auto ff = ejungwoo::fitgg(hh,5);
    ff -> SetName(TString("fit_")+gg->GetName());
    gg -> Fit(ff,"Q0");
    ff -> SetLineColor(gg->GetLineColor());
    ff -> SetNpx(1000);
    ff -> Draw("samel");

    gglist.push_back(gg);
    fflist.push_back(ff);
  }

  file_summary -> cd();
  for (auto gg : gglist) gg -> Write();
  for (auto ff : fflist) ff -> Write();
  cout << name_summary << endl;
}

void get_par_limits(double p_current, int idx_particle, TH1D * hist, double &a, double &m, double &s, double &a1, double &a2, double &m1, double &m2, double &s1, double &s2)
{
       if (idx_particle == stp::kpin) { }
  else if (idx_particle == stp::kpip) { }
  else if (idx_particle == stp::kp  ) { auto s_sigma0 = 4.;  auto s_mean0 = 40.; s = s_sigma0 / s_mean0 * m; }
  else if (idx_particle == stp::kd  ) { auto s_sigma0 = 10.; auto s_mean0 = 80.; s = s_sigma0 / s_mean0 * m; }
  else if (idx_particle == stp::kt  ) { auto s_sigma0 = 25.; auto s_mean0 = 180.; s = s_sigma0 / s_mean0 * m; }
  else if (idx_particle == stp::khe3) { auto s_sigma0 = 22.; auto s_mean0 = 250.; s = s_sigma0 / s_mean0 * m; }
  else if (idx_particle == stp::khe4) { auto s_sigma0 = 35.; auto s_mean0 = 350.; s = s_sigma0 / s_mean0 * m; }
  else if (idx_particle == stp::khe6) { auto s_sigma0 = 21.; auto s_mean0 = 240.; s = s_sigma0 / s_mean0 * m; }
  else if (idx_particle == stp::kli6) { }
  else if (idx_particle == stp::kli7) { }
  else if (idx_particle == stp::kele) { }
  else if (idx_particle == stp::kpos) { }
  else if (idx_particle == stp::kli ) { auto s_sigma0 = 60.; auto s_mean0 = 400.; s = s_sigma0 / s_mean0 * m; }

  auto dm = c_limit_all * m;

  a1 = 10;
  a2 = hist -> GetMaximum();
  m1 = m - dm;
  m2 = m + dm;
  s1 = 1;
  s2 = 100;

  auto ccc = c_limit_all;
       if (idx_particle==stp::kp)   ccc = c_limit_p;
  else if (idx_particle>=stp::kli6) ccc = c_limit_oli;

  /**/ if (idx_particle == stp::kpin)  { }
  else if (idx_particle == stp::kpip)  { }
  else if (idx_particle == stp::kp  )  { s1 = (1.-ccc)*s; s2 = (1.+ccc)*s; }
  else if (idx_particle == stp::kd  )  { s1 = (1.-ccc)*s; s2 = (1.+ccc)*s; }
  else if (idx_particle == stp::kt  )  { s1 = (1.-ccc)*s; s2 = (1.+ccc)*s; }
  else if (idx_particle == stp::khe3)  { s1 = (1.-ccc)*s; s2 = (1.+ccc)*s; }
  else if (idx_particle == stp::khe4)  { s1 = (1.-ccc)*s; s2 = (1.+ccc)*s; }
  else if (idx_particle == stp::khe6)  { s1 = (1.-ccc)*s; s2 = (1.+ccc)*s; }
  else if (idx_particle == stp::kli6)  { s1 = (1.-ccc)*s; s2 = (1.+ccc)*s; }
  else if (idx_particle == stp::kli7)  { }
  else if (idx_particle == stp::kele)  { }
  else if (idx_particle == stp::kpos)  { }
  else if (idx_particle == stp::kli )  { s1 = (1.-ccc)*s; s2 = (1.+ccc)*s; }

  if (a < a1 || a > a2) a = (a1 + a2)/2.;
  if (m < m1 || m > m2) m = (m1 + m2)/2.;
  if (s < s1 || s > s2) s = (s1 + s2)/2.;
}

void get_par_limits2(double s_current, int idx_particle, TH1D * hist, double &a, double &m, double &s, double &a1, double &a2, double &m1, double &m2, double &s1, double &s2)
{
  auto dm = c_limit_all * m;

  a1 = 10;
  a2 = hist -> GetMaximum();
  m1 = m - dm;
  m2 = m + dm;
  s1 = 10;
  s2 = 60;

  if (a < a1 || a > a2) a = (a1 + a2)/2.;
  if (m < m1 || m > m2) m = (m1 + m2)/2.;
  if (s < s1 || s > s2) s = (s1 + s2)/2.;
}
