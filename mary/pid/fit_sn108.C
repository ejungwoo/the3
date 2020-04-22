#include "stp.h"

double c_limit_p = .4;
double c_limit_oli = .20;
double c_limit_all = .20;

void get_par_limits(double p_current, int idx_particle, TH1D * hist, double &a, double &m, double &s, double &a1, double &a2, double &m1, double &m2, double &s1, double &s2);
void get_par_limits2(double s_current, int idx_particle, TH1D * hist, double &a, double &m, double &s, double &a1, double &a2, double &m1, double &m2, double &s1, double &s2);

void fit_sn108(bool verbose = 1, int index_py = -1)
{
  int system = 108;
  ejungwoo::gzcolor(0);
  ejungwoo::gstat("ne");
  ejungwoo::gsave(0);
  ejungwoo::gcvspos(1300,0);
  ejungwoo::gsetjumpcvs();

  bool write_summary = true;
  bool draw_proj1 = true;
  bool draw_proj2 = false;
  bool draw_fit_status = true;
  bool draw_refit = false;
  bool draw_pid = true;

  double p_start = 100.;
  double p_end = 2500.;
  int bin_size_p = 8;
  double s_start = 100.;
  double s_end = 1300.;
  int bin_size_s = 20;

  int idx_particles[] = { stp::kpip, stp::kp, stp::kd, stp::kt, stp::khe3, stp::khe4, stp::khe6, stp::kli, };

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

  double p_limit1[] = {
  /*stp::kpin*/ p_hist_min, 
  /*stp::kpip*/ 100,
  /*stp::kp  */ 350,
  /*stp::kd  */ 400,
  /*stp::kt  */ 500,
  /*stp::khe3*/ 500,
  /*stp::khe4*/ 500,
  /*stp::khe6*/ 1700,
  /*stp::kli6*/ 1000,
  /*stp::kli7*/ 1500,
  /*stp::kele*/ p_hist_min, 
  /*stp::kpos*/ p_hist_min,
  /*stp::kli */ 1000,
  };

  double p_limit2[] = {
  /*stp::kpin*/ p_hist_max, 
  /*stp::kpip*/ 500,
  /*stp::kp  */ 1700,
  /*stp::kd  */ 2200,
  /*stp::kt  */ p_hist_max, 
  ///*stp::khe3*/ 1500,
  /*stp::khe3*/ 2000,
  /*stp::khe4*/ p_hist_max, 
  /*stp::khe6*/ p_hist_max, 
  /*stp::kli6*/ 2000,
  /*stp::kli7*/ 2300,
  /*stp::kele*/ p_hist_max, 
  /*stp::kpos*/ p_hist_max,
  /*stp::kli */ 2000,
  };

  double s_limit1[] = {
  /*stp::kpin*/ 0,
  /*stp::kpip*/ 100,
  /*stp::kp  */ 500,
  /*stp::kd  */ 700,
  /*stp::kt  */ 700,
  /*stp::khe3*/ 700,
  /*stp::khe4*/ 700,
  /*stp::khe6*/ 0,
  /*stp::kli6*/ 0,
  /*stp::kli7*/ 0,
  /*stp::kele*/ 0,
  /*stp::kpos*/ 0,
  /*stp::kli */ 0,
  };

  double s_limit2[] = {
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

  //TString nameVersion = Form("sys108_p%dy%d",pit,yaw);
  TString nameVersion = Form("sys108_ipy%d",index_py);
  if (index_py<0) nameVersion = "sys108";
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
      auto fileTree = new TFile(Form("data__in/pdedx%d2.root",system));
      auto tree = (TTree *) fileTree -> Get("data");
      auto cut = TCut();
      //if (pit>0&&yaw>0)
      if (index_py==0)
        cut = cut + TCut("ipy==1515||ipy==1516||ipy==1615||ipy==1616");
        //cut = cut + TCut("ipy==1515");
      cout << cut << endl;
      hist_pid = (TH2D *) ejungwoo::tp(Form("%s",name_ana.Data()),tree,"dedx:p",cut,"",nbins_p,p_min,p_max,nbins_s,s_min,s_max);
      ejungwoo::write(hist_pid);
    }
  }

  TCanvas *cvs_pid;
  if (draw_pid) {
    cvs_pid = ejungwoo::cc4();
    cvs_pid -> SetLogz();
    hist_pid -> Draw("colz");
  }

  double parameters[50] = {0};
  double parlimits1[50] = {0};
  double parlimits2[50] = {0};
  TGraphErrors *graph_ams[2][4][stp::fNumParticles];

  TF1 *prefit_dedx[stp::fNumParticles] = {0};
  for (auto idx_particle : idx_particles)
  {
    TString name_particle = stp::fNameShort[idx_particle];

    graph_ams[0][0][idx_particle] = ejungwoo::new_ge("ampl_"+name_particle+"_s");
    graph_ams[0][1][idx_particle] = ejungwoo::new_ge("mean_"+name_particle+"_s");
    graph_ams[0][2][idx_particle] = ejungwoo::new_ge("sigm_"+name_particle+"_s");
    graph_ams[0][3][idx_particle] = ejungwoo::new_ge("rslt_"+name_particle+"_s");
    graph_ams[1][0][idx_particle] = ejungwoo::new_ge("ampl_"+name_particle);
    graph_ams[1][1][idx_particle] = ejungwoo::new_ge("mean_"+name_particle);
    graph_ams[1][2][idx_particle] = ejungwoo::new_ge("sigm_"+name_particle);
    graph_ams[1][3][idx_particle] = ejungwoo::new_ge("rslt_"+name_particle);

    auto file_fit = new TFile(TString("data__prefit/fit_dedx_par5_region_")+name_particle+".prefit.root");
    prefit_dedx[idx_particle] = (TF1 *) file_fit -> Get(TString("fit_dedx_par5_region_")+name_particle);
    prefit_dedx[idx_particle] -> Draw("samel");
  }

  TGraphErrors *graph_ampl[stp::fNumParticles];
  TGraphErrors *graph_mean[stp::fNumParticles];
  TGraphErrors *graph_sigm[stp::fNumParticles];
  TGraphErrors *graph_rslt[stp::fNumParticles];

  for (auto idx_particle : idx_particles) {
    TString name_particle = stp::fNameShort[idx_particle];
    graph_ampl[idx_particle] = ejungwoo::new_ge(name_particle+"_refit_ampl");
    graph_mean[idx_particle] = ejungwoo::new_ge(name_particle+"_refit_mean");
    graph_sigm[idx_particle] = ejungwoo::new_ge(name_particle+"_refit_sigm");
    graph_rslt[idx_particle] = ejungwoo::new_ge(name_particle+"_refit_rslt");

    graph_ampl[idx_particle] -> SetMarkerColor(stp::fColor[idx_particle]);
    graph_mean[idx_particle] -> SetMarkerColor(stp::fColor[idx_particle]);
    graph_sigm[idx_particle] -> SetMarkerColor(stp::fColor[idx_particle]);
    graph_rslt[idx_particle] -> SetMarkerColor(stp::fColor[idx_particle]);
  }


  int bin_s_start = hist_pid -> GetYaxis() -> FindBin(s_start);
  int bin_s_end = hist_pid -> GetYaxis() -> FindBin(s_end);

  for (auto bin1 = bin_s_end; bin1 > bin_s_start; bin1-=bin_size_s)
  {
    bool idx_is_good[stp::fNumParticles] = {0};

    auto bin2 = bin1-bin_size_p+1;
    auto s_bin_up = hist_pid -> GetXaxis() -> GetBinLowEdge(bin1);
    auto s_bin_low = hist_pid -> GetXaxis() -> GetBinUpEdge(bin2);
    auto ds_bin = s_bin_up - s_bin_low;

    auto s_current1 = hist_pid -> GetYaxis() -> GetBinCenter(bin1);
    auto s_current2 = hist_pid -> GetYaxis() -> GetBinCenter(bin2);
    auto s_current = (s_current1 + s_current2)/2.;
    TString name_s = Form("proj_bin_s%.2f_b%d_%d",s_current,bin1,bin2);
    auto hist_proj = hist_pid -> ProjectionX(name_s,bin2,bin1);
    if (hist_proj -> GetEntries() < 200)
      continue;
    hist_proj -> SetTitle(name_s);

    TCanvas *cvs_proj = nullptr;
    TCanvas *cvs_proj1 = nullptr;
    TCanvas *cvs_proj2 = nullptr;

    double p_ampl_array[stp::fNumParticles] = {0};
    double p_mean_array[stp::fNumParticles] = {0};
    //double p_sigm_array[stp::fNumParticles] = {0};

    auto graph_init = ejungwoo::new_ge(name_s+"_inits");
    graph_init -> SetMarkerColor(kRed);
    graph_init -> SetLineColor(kRed);

    bool particle_is_in_range = false;

    for (auto idx_particle : idx_particles)
    {
      auto fit = prefit_dedx[idx_particle];
      double srange1 = s_limit1[idx_particle];
      double srange2 = s_limit2[idx_particle];
      if (s_current > srange1 && s_current < srange2)
      {
        idx_is_good[idx_particle] = true;
        particle_is_in_range = true;
        auto p_mean_init = fit -> GetX(s_current);

        p_ampl_array[idx_particle] = hist_proj -> GetBinContent(hist_proj->GetXaxis()->FindBin(p_mean_init));
        p_mean_array[idx_particle] = p_mean_init;
      }
    }

    if (!particle_is_in_range)
      continue;

    if (draw_proj1) {
      cvs_proj = ejungwoo::cv(name_s, 500,800);
      ejungwoo::div(cvs_proj,1,2);
      cvs_proj1 = (TCanvas *) cvs_proj -> cd(1);
      cvs_proj2 = (TCanvas *) cvs_proj -> cd(2);

      cvs_proj1 -> cd();
      cvs_proj1 -> SetLogy();
      hist_proj -> Draw();

      cvs_proj2 -> cd();
      hist_proj -> Draw();
    }

    TString expression_sel = "0";
    int count_s = 0;
    for (auto idx_particle : idx_particles)
    {
      TString name_particle = stp::fNameShort[idx_particle];

      if (!idx_is_good[idx_particle])
        continue;

      double p_mean_init = p_mean_array[idx_particle];
      if (p_mean_init <= 0)
        continue;

      double p_ampl_init = hist_proj -> GetBinContent(hist_proj->GetXaxis()->FindBin(p_mean_init));

      if (p_ampl_init < 2.) {
        idx_is_good[idx_particle] = false;
        continue;
      }

      expression_sel = expression_sel + Form("+gaus(%d)",3*count_s);

      double p_sigm_init = (idx_particle+1)*2;

      double a1, a2, m1, m2, s1, s2;

      get_par_limits2(s_current, idx_particle, hist_proj, p_ampl_init, p_mean_init, p_sigm_init, a1, a2, m1, m2, s1, s2);

      parameters[3*idx_particle+0] = p_ampl_init;
      parameters[3*idx_particle+1] = p_mean_init;
      parameters[3*idx_particle+2] = p_sigm_init;
      parlimits1[3*idx_particle+0] = a1;
      parlimits2[3*idx_particle+0] = a2;
      parlimits1[3*idx_particle+1] = m1;
      parlimits2[3*idx_particle+1] = m2;
      parlimits1[3*idx_particle+2] = s1;
      parlimits2[3*idx_particle+2] = s2;

      if (draw_proj1 && draw_fit_status) {
        cvs_proj1 -> cd();
        graph_init -> SetPoint(idx_particle, p_mean_init,100);
        graph_init -> SetPointError(idx_particle, (m2-m1)/2., 0);
      }

      count_s++;
    }

    if (draw_proj1 && draw_fit_status) {
      cvs_proj1 -> cd();
      graph_init -> Sort();
      graph_init -> Draw("samepz1");
    }

    int count_particle = 0;
    TF1 *fit_p_all = new TF1("fit_p_all",expression_sel,0,1000);
    TF1 *fit_p_sel = new TF1("fit_p_sel",expression_sel,0,1000);
    for (auto idx_particle : idx_particles)
    {
      if (!idx_is_good[idx_particle])
        continue;

      TString name_particle = stp::fNameShort[idx_particle];

      fit_p_sel -> SetParameter(3*count_particle+0,parameters[3*idx_particle+0]);
      fit_p_sel -> SetParameter(3*count_particle+1,parameters[3*idx_particle+1]);
      fit_p_sel -> SetParameter(3*count_particle+2,parameters[3*idx_particle+2]);

      fit_p_sel -> SetParLimits(3*count_particle+0,parlimits1[3*idx_particle+0],parlimits2[3*idx_particle+0]);
      fit_p_sel -> SetParLimits(3*count_particle+1,parlimits1[3*idx_particle+1],parlimits2[3*idx_particle+1]);
      fit_p_sel -> SetParLimits(3*count_particle+2,parlimits1[3*idx_particle+2],parlimits2[3*idx_particle+2]);

      count_particle++;
    }

    hist_proj -> Fit(fit_p_sel,"Q0");

    count_particle = 0;
    auto legend = new TLegend();
    for (auto idx_particle : idx_particles)
    {
      if (!idx_is_good[idx_particle])
        continue;

      TString name_particle = stp::fNameShort[idx_particle];

      auto p_ampl = fit_p_sel -> GetParameter(3*count_particle+0);
      auto p_mean = fit_p_sel -> GetParameter(3*count_particle+1);
      auto p_sigm = fit_p_sel -> GetParameter(3*count_particle+2);

      auto marker = new TMarker(p_mean,p_ampl,30);
      count_particle++;
      marker -> SetMarkerSize(1.5);
      marker -> SetMarkerColor(kBlack);
      if (draw_proj1 && draw_fit_status) {
        cvs_proj1 -> cd();
        marker -> Draw("samep");
        auto fit_p_single = new TF1(Form("fit_p_single_%s",name_particle.Data()),"gaus(0)",0,1000);
        fit_p_single -> SetParameters(p_ampl,p_mean,p_sigm);
        fit_p_single -> SetLineColor(stp::fColor[idx_particle]);
        fit_p_single -> Draw("samel");

        legend -> AddEntry(fit_p_single,Form("%s: %.2e, %.2e, %.2e",name_particle.Data(), p_ampl,p_mean,p_sigm),"l");

        cvs_proj2 -> cd();
        marker -> Draw("samep");
        fit_p_single -> Draw("samel");

        //legend -> AddEntry(fit_p_single,Form("%s: %.2e, %.2e, %.2e",name_particle.Data(), p_ampl,p_mean,p_sigm),"l");
      }

      auto idxn = graph_ams[0][0][idx_particle] -> GetN();
      graph_ams[0][0][idx_particle] -> SetPoint(idxn, s_current, p_ampl);
      graph_ams[0][1][idx_particle] -> SetPoint(idxn, s_current, p_mean);
      graph_ams[0][1][idx_particle] -> SetPointError(idxn, ds_bin/2., p_sigm);
      graph_ams[0][2][idx_particle] -> SetPoint(idxn, s_current, p_sigm);
      graph_ams[0][3][idx_particle] -> SetPoint(idxn, s_current, p_mean/p_sigm);

      graph_ams[1][1][idx_particle] -> SetPoint(idxn, p_mean, s_current);
      graph_ams[1][1][idx_particle] -> SetPointError(idxn, p_sigm, ds_bin/2.);

      idxn = graph_mean[idx_particle] -> GetN();
      graph_mean[idx_particle] -> SetPoint(idxn, p_mean, s_current);
      graph_mean[idx_particle] -> SetPointError(idxn, p_sigm, ds_bin/2.);
    }

    if (draw_proj1 && draw_fit_status) {
      cvs_proj1 -> cd();
      fit_p_sel -> SetNpx(2000);
      fit_p_sel -> SetLineColor(kGray+2);
      fit_p_sel -> Draw("samel");
      ejungwoo::make_l(cvs_proj1,legend,0,0,0.6,0.4) -> Draw("same");

      cvs_proj2 -> cd();
      //fit_p_sel -> SetNpx(2000);
      //fit_p_sel -> SetLineColor(kGray+2);
      fit_p_sel -> Draw("samel");
      //ejungwoo::make_l(cvs_proj2,legend,0,0,0.6,0.4) -> Draw("same");

      auto npar = fit_p_sel -> GetNpar();
      auto p_mean_last = fit_p_sel -> GetParameter(npar-2);
      auto p_sigm_last = fit_p_sel -> GetParameter(npar-1);

      hist_proj -> GetXaxis() -> SetRangeUser(0,p_mean_last+2*p_sigm_last);

      ejungwoo::make_c(cvs_proj1);
      ejungwoo::make_c(cvs_proj2);
    }
  }

  int bin_p_start = hist_pid -> GetXaxis() -> FindBin(p_start);
  int bin_p_end = hist_pid -> GetXaxis() -> FindBin(p_end);

  ////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////
  // XXX
  vector<double> list_mom;
  vector<double> list_dp;
  vector<TH1D *> list_hist;
  vector<TF1 *> list_fsel;
  ////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////

  for (auto bin1 = bin_p_start; bin1 < bin_p_end; bin1+=bin_size_p)
  {
    bool idx_is_good[stp::fNumParticles] = {0};

    auto bin2 = bin1+bin_size_p-1;
    auto p_bin_low = hist_pid -> GetXaxis() -> GetBinLowEdge(bin1);
    auto p_bin_up = hist_pid -> GetXaxis() -> GetBinUpEdge(bin2);
    auto dp_bin = p_bin_up - p_bin_low;
    auto p_current1 = hist_pid -> GetXaxis() -> GetBinCenter(bin1);
    auto p_current2 = hist_pid -> GetXaxis() -> GetBinCenter(bin2);
    auto p_current = (p_current1 + p_current2)/2.;
    TString name_p = Form("proj_bin_p%.2f_b%d_%d",p_current,bin1,bin2);
    auto hist_proj = hist_pid -> ProjectionY(name_p,bin1,bin2);
    list_hist.push_back(hist_proj);
    list_mom.push_back(p_current);
    list_dp.push_back(dp_bin);
    if (hist_proj -> GetEntries() < 200)
      continue;
    hist_proj -> SetTitle(name_p);

    TCanvas *cvs_proj = nullptr;
    TCanvas *cvs_proj1 = nullptr;
    TCanvas *cvs_proj2 = nullptr;

    if (draw_proj2) {
      cvs_proj = ejungwoo::cv(name_p, 500,800);
      ejungwoo::div(cvs_proj,1,2);
      cvs_proj1 = (TCanvas *) cvs_proj -> cd(1);
      cvs_proj2 = (TCanvas *) cvs_proj -> cd(2);

      cvs_proj1 -> cd();
      hist_proj -> Draw();

      cvs_proj2 -> cd();
      cvs_proj2 -> SetLogy();
      hist_proj -> Draw();
    }

    double s_ampl_array[stp::fNumParticles] = {0};
    double s_mean_array[stp::fNumParticles] = {0};

    auto graph_init = ejungwoo::new_ge(name_p+"_inits");
    graph_init -> SetMarkerColor(kRed);
    graph_init -> SetLineColor(kRed);

    for (auto idx_particle : idx_particles)
    {
      auto fit = prefit_dedx[idx_particle];
      double prange1 = p_limit1[idx_particle];
      double prange2 = p_limit2[idx_particle];
      if (p_current > prange1 && p_current < prange2)
      {
        idx_is_good[idx_particle] = true;
        auto s_mean_init = fit -> Eval(p_current);

        s_ampl_array[idx_particle] = hist_proj -> GetBinContent(hist_proj->GetXaxis()->FindBin(s_mean_init));
        s_mean_array[idx_particle] = s_mean_init;
      }
    }

    TString expression_sel = "0";
    int count_s = 0;
    for (auto idx_particle : idx_particles)
    {
      TString name_particle = stp::fNameShort[idx_particle];

      if (!idx_is_good[idx_particle])
        continue;

      double s_mean_init = s_mean_array[idx_particle];
      if (s_mean_init <= 0)
        continue;

      double s_ampl_init = hist_proj -> GetBinContent(hist_proj->GetXaxis()->FindBin(s_mean_init));

      if (s_ampl_init < 2.) {
        idx_is_good[idx_particle] = false;
        continue;
      }

      expression_sel = expression_sel + Form("+gaus(%d)",3*count_s);

      double s_sigm_init = (idx_particle+1)*2;

      double a1, a2, m1, m2, s1, s2;

      get_par_limits(p_current, idx_particle, hist_proj, s_ampl_init, s_mean_init, s_sigm_init, a1, a2, m1, m2, s1, s2);

      parameters[3*idx_particle+0] = s_ampl_init;
      parameters[3*idx_particle+1] = s_mean_init;
      parameters[3*idx_particle+2] = s_sigm_init;
      parlimits1[3*idx_particle+0] = a1;
      parlimits2[3*idx_particle+0] = a2;
      parlimits1[3*idx_particle+1] = m1;
      parlimits2[3*idx_particle+1] = m2;
      parlimits1[3*idx_particle+2] = s1;
      parlimits2[3*idx_particle+2] = s2;

      if (draw_proj2 && draw_fit_status) {
        cvs_proj1 -> cd();
        graph_init -> SetPoint(idx_particle, s_mean_init,100);
        graph_init -> SetPointError(idx_particle, (m2-m1)/2., 0);
      }

      count_s++;
    }

    if (draw_proj2 && draw_fit_status) {
      cvs_proj1 -> cd();
      graph_init -> Sort();
      graph_init -> Draw("samepz1");
    }

    int count_particle = 0;
    TF1 *fit_s_sel = new TF1("fit_s_sel",expression_sel,0,1000);
    list_fsel.push_back(fit_s_sel);
    for (auto idx_particle : idx_particles)
    {
      if (!idx_is_good[idx_particle])
        continue;

      TString name_particle = stp::fNameShort[idx_particle];

      fit_s_sel -> SetParameter(3*count_particle+0,parameters[3*idx_particle+0]);
      fit_s_sel -> SetParameter(3*count_particle+1,parameters[3*idx_particle+1]);
      fit_s_sel -> SetParameter(3*count_particle+2,parameters[3*idx_particle+2]);

      fit_s_sel -> SetParLimits(3*count_particle+0,parlimits1[3*idx_particle+0],parlimits2[3*idx_particle+0]);
      fit_s_sel -> SetParLimits(3*count_particle+1,parlimits1[3*idx_particle+1],parlimits2[3*idx_particle+1]);
      fit_s_sel -> SetParLimits(3*count_particle+2,parlimits1[3*idx_particle+2],parlimits2[3*idx_particle+2]);

      count_particle++;
    }

    hist_proj -> Fit(fit_s_sel,"Q0");

    count_particle = 0;
    auto legend = new TLegend();
    for (auto idx_particle : idx_particles)
    {
      if (!idx_is_good[idx_particle]) {
        continue;
      }

      TString name_particle = stp::fNameShort[idx_particle];

      auto s_ampl = fit_s_sel -> GetParameter(3*count_particle+0);
      auto s_mean = fit_s_sel -> GetParameter(3*count_particle+1);
      auto s_sigm = fit_s_sel -> GetParameter(3*count_particle+2);

      auto marker = new TMarker(s_mean,s_ampl,30);
      count_particle++;
      marker -> SetMarkerSize(1.5);
      marker -> SetMarkerColor(kBlack);
      if (draw_proj2 && draw_fit_status) {
        cvs_proj1 -> cd();
        marker -> Draw("samep");
        auto fit_s_single = new TF1(Form("fit_s_single_%s",name_particle.Data()),"gaus(0)",0,1000);
        fit_s_single -> SetParameters(s_ampl,s_mean,s_sigm);
        fit_s_single -> SetLineColor(stp::fColor[idx_particle]);
        fit_s_single -> Draw("samel");

        legend -> AddEntry(fit_s_single,Form("%s: %.2e, %.2e, %.2e",name_particle.Data(), s_ampl,s_mean,s_sigm),"l");

        cvs_proj2 -> cd();
        marker -> Draw("samep");
        fit_s_single -> Draw("samel");

        legend -> AddEntry(fit_s_single,Form("%s: %.2e, %.2e, %.2e",name_particle.Data(), s_ampl,s_mean,s_sigm),"l");
      }

      auto idxn = graph_ams[1][1][idx_particle] -> GetN();
      graph_ams[1][1][idx_particle] -> SetPoint(idxn, p_current, s_mean);
      graph_ams[1][1][idx_particle] -> SetPointError(idxn, dp_bin/2., s_sigm);

      idxn = graph_ams[1][0][idx_particle] -> GetN();
      graph_ams[1][0][idx_particle] -> SetPoint(idxn, p_current, s_ampl);
      graph_ams[1][2][idx_particle] -> SetPoint(idxn, p_current, s_sigm);
      graph_ams[1][3][idx_particle] -> SetPoint(idxn, p_current, s_mean/s_sigm);
    }

    if (draw_proj2 && draw_fit_status) {
      cvs_proj1 -> cd();
      fit_s_sel -> SetNpx(2000);
      fit_s_sel -> SetLineColor(kGray+2);
      fit_s_sel -> Draw("samel");
      ejungwoo::make_l(cvs_proj1,legend,0,0,0.6,0.4) -> Draw("same");

      cvs_proj2 -> cd();
      fit_s_sel -> Draw("samel");

      auto npar = fit_s_sel -> GetNpar();
      auto s_mean_last = fit_s_sel -> GetParameter(npar-2);
      auto s_sigm_last = fit_s_sel -> GetParameter(npar-1);

      hist_proj -> GetXaxis() -> SetRangeUser(0,s_mean_last+2*s_sigm_last);

      ejungwoo::make_c(cvs_proj1);
      ejungwoo::make_c(cvs_proj2);
    }
  }

  ////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////

  auto name_summary = ejungwoo::name_full("summary");
  auto file_summary = new TFile(name_summary,"recreate");

  TF1 *fit_ampl[stp::fNumParticles] = {0};
  TF1 *fit_mean[stp::fNumParticles] = {0};
  TF1 *fit_sigm[stp::fNumParticles] = {0};

  for (auto idx_particle : idx_particles)
  {
    TString name_particle = stp::fNameShort[idx_particle];
    cout << name_particle << endl;

    TString name_amp = Form("fit_ampl_%s",name_particle.Data());
    auto gg = graph_ams[1][0][idx_particle];
    gg -> SetMarkerColor(stp::fColor[idx_particle]);
    auto hh = ejungwoo::tohist(gg,0,2500);
    fit_ampl[idx_particle] = ejungwoo::fitgg(hh,5);
    fit_ampl[idx_particle] -> SetName(name_amp);
    gg -> Fit(fit_ampl[idx_particle],"Q0");
    fit_ampl[idx_particle] -> SetRange(ejungwoo::x1_g(gg)-100,ejungwoo::x2_g(gg)+100);
    fit_ampl[idx_particle] -> SetLineColor(stp::fColor[idx_particle]);
    fit_ampl[idx_particle] -> SetNpx(2000);
    ejungwoo::addto("ampl",gg,"samep colorx");
    ejungwoo::addto("ampl",fit_ampl[idx_particle],"samepl addx colorx");



    TString name_mean = Form("fit_mean_%s",name_particle.Data());
    fit_mean[idx_particle] = stp::f1_dedx_p5(idx_particle,p_hist_min,p_hist_max);
    fit_mean[idx_particle] -> SetParameters(
        prefit_dedx[idx_particle] -> GetParameter(0),
        prefit_dedx[idx_particle] -> GetParameter(1),
        prefit_dedx[idx_particle] -> GetParameter(2),
        prefit_dedx[idx_particle] -> GetParameter(3),
        prefit_dedx[idx_particle] -> GetParameter(4));
    fit_mean[idx_particle] -> SetName(name_mean);
    gg = graph_ams[1][1][idx_particle];
    gg -> SetMarkerColor(stp::fColor[idx_particle]);
    gg -> Fit(fit_mean[idx_particle],"Q0");
    if (idx_particle == stp::kp || idx_particle == stp::kd || idx_particle == stp::kt || idx_particle == stp::khe3 || idx_particle == stp::khe4)
      fit_mean[idx_particle] -> SetRange(fit_mean[idx_particle] -> GetX(1500.), ejungwoo::x2_g(gg)+500);
    else
      fit_mean[idx_particle] -> SetRange(ejungwoo::x1_g(gg)-10, ejungwoo::x2_g(gg)+10);
    fit_mean[idx_particle] -> SetLineColor(stp::fColor[idx_particle]);
    fit_mean[idx_particle] -> SetNpx(2000);
    ejungwoo::addto("mean",gg,"samepz colorx");
    ejungwoo::addto("mean",fit_mean[idx_particle],"samel addx colorx");



    TString name_sigm = Form("fit_sigm_%s",name_particle.Data());
    fit_sigm[idx_particle] = stp::f1_dedx_p5(idx_particle,p_hist_min,p_hist_max);
    fit_sigm[idx_particle] -> SetParameters(
        prefit_dedx[idx_particle] -> GetParameter(0)*0.12,
        prefit_dedx[idx_particle] -> GetParameter(1),
        prefit_dedx[idx_particle] -> GetParameter(2),
        prefit_dedx[idx_particle] -> GetParameter(3),
        prefit_dedx[idx_particle] -> GetParameter(4));
    fit_sigm[idx_particle] -> SetName(name_sigm);
    gg = graph_ams[1][2][idx_particle];
    gg -> SetMarkerColor(stp::fColor[idx_particle]);
    gg -> Fit(fit_sigm[idx_particle],"Q0");
    if (idx_particle == stp::kp || idx_particle == stp::kd || idx_particle == stp::kt || idx_particle == stp::khe3 || idx_particle == stp::khe4)
      fit_sigm[idx_particle] -> SetRange(fit_sigm[idx_particle] -> GetX(200.), ejungwoo::x2_g(gg)+500);
    else
      fit_sigm[idx_particle] -> SetRange(ejungwoo::x1_g(gg)-10, ejungwoo::x2_g(gg)+10);
    fit_sigm[idx_particle] -> SetLineColor(stp::fColor[idx_particle]);
    fit_sigm[idx_particle] -> SetNpx(2000);
    ejungwoo::addto("sigm",gg,"samepz colorx");
    ejungwoo::addto("sigm",fit_sigm[idx_particle],"samel addx colorx");

    gg = graph_ams[1][3][idx_particle];
    gg -> SetMarkerColor(stp::fColor[idx_particle]);
    ejungwoo::addto("rslt",gg,"samepz colorx");
  }

  auto cvs_ampl = ejungwoo::drawc("ampl");
  auto cvs_mean = ejungwoo::drawc("mean");
  auto cvs_sigm = ejungwoo::drawc("sigm");
  auto cvs_rslt = ejungwoo::drawc("rslt");

  //return;

  ////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////

  double pr_limit1[] = {
  /*stp::kpin*/ p_hist_min, 
  /*stp::kpip*/ 100,
  /*stp::kp  */ 350,
  /*stp::kd  */ 400,
  /*stp::kt  */ 700,
  /*stp::khe3*/ 700,
  /*stp::khe4*/ 700,
  /*stp::khe6*/ 700,
  /*stp::kli6*/ 700,
  /*stp::kli7*/ 700,
  /*stp::kele*/ p_hist_min, 
  /*stp::kpos*/ p_hist_min,
  /*stp::kli */ 1000,
  };

  double pr_limit2[] = {
  /*stp::kpin*/ 2500,
  /*stp::kpip*/ 1000,
  /*stp::kp  */ 2500,
  /*stp::kd  */ 2500,
  /*stp::kt  */ 2500,
  /*stp::khe3*/ 2500,
  /*stp::khe4*/ 2500,
  /*stp::khe6*/ 2500,
  /*stp::kli6*/ 2500,
  /*stp::kli7*/ 2500,
  /*stp::kele*/ 2500,
  /*stp::kpos*/ 2500,
  /*stp::kli */ 2500,
  };

  auto numHists = list_hist.size();
  for (auto idxp=0; idxp<numHists; ++idxp)
  {
    double p_center = list_mom.at(idxp);
    double dp = list_dp.at(idxp);
    if (p_center < 350)
      continue;

    auto hist = list_hist.at(idxp);
    auto fit_sel = list_fsel.at(idxp);

    TString expression_all = "0";
    int count = 0;
    for (auto idx_particle : idx_particles)
      expression_all = expression_all + Form("+gaus(%d)",3*count++);

    auto fit_refit = new TF1(Form("fit_refit_p%.2f",p_center),expression_all,0,1000);

    count = 0;
    for (auto idx_particle : idx_particles)
    {
      double prange1 = pr_limit1[idx_particle];
      double prange2 = pr_limit2[idx_particle];

      auto a0 = fit_ampl[idx_particle] -> Eval(p_center);
      auto m0 = fit_mean[idx_particle] -> Eval(p_center);
      auto s0 = fit_sigm[idx_particle] -> Eval(p_center);

      if (!(p_center > prange1 && p_center < prange2)) {
        fit_refit -> FixParameter(3*count+0,0);
        fit_refit -> FixParameter(3*count+1,m0);
        fit_refit -> FixParameter(3*count+2,s0);
        count++;
        continue;
      }

      fit_refit -> SetParameter(3*count+0,a0);
      fit_refit -> SetParLimits(3*count+0,a0-0.5*a0,a0+0.2*a0);
      fit_refit -> SetParameter(3*count+1,m0);
      fit_refit -> SetParLimits(3*count+1,m0-0.1*m0,m0+0.1*m0);
      fit_refit -> SetParameter(3*count+2,s0);
      fit_refit -> SetParLimits(3*count+2,s0-0.1*s0,s0+0.1*s0);

      if (idx_particle==stp::khe3) {
        fit_refit -> SetParLimits(3*count+0,a0-0.2*a0,a0+3*a0);
      }
      //if (p_center <= 1500 && idx_particle==stp::khe6) {
      if (p_center <= 1400 && idx_particle==stp::khe6) {
        fit_refit -> FixParameter(3*count+0,0);
        fit_refit -> FixParameter(3*count+1,m0);
        fit_refit -> FixParameter(3*count+2,s0);
      }
      count++;
    }
    hist -> Fit(fit_refit,"Q0");
    fit_refit -> SetNpx(1000);
    fit_refit -> SetLineColor(kBlack);

    count = 0;
    for (auto idx_particle : idx_particles)
    {
      auto a1 = fit_refit -> GetParameter(3*count+0);
      auto m1 = fit_refit -> GetParameter(3*count+1);
      auto s1 = fit_refit -> GetParameter(3*count+2);

      if (a1!=0) {
        graph_ampl[idx_particle] -> SetPoint(graph_ampl[idx_particle]->GetN(), p_center, a1);
        graph_mean[idx_particle] -> SetPoint(graph_mean[idx_particle]->GetN(), p_center, m1);
        graph_mean[idx_particle] -> SetPointError(graph_mean[idx_particle]->GetN()-1, dp/2., s1);
        graph_sigm[idx_particle] -> SetPoint(graph_sigm[idx_particle]->GetN(), p_center, s1);
        graph_rslt[idx_particle] -> SetPoint(graph_rslt[idx_particle]->GetN(), p_center, m1/s1);
      }
      count++;
    }


    auto npar = fit_sel -> GetNpar();
    auto s_mean_last = fit_sel -> GetParameter(npar-2);
    auto s_sigm_last = fit_sel -> GetParameter(npar-1);

    auto tttttttttt = s_mean_last+2*s_sigm_last;
    if (tttttttttt < 800)
      tttttttttt = 800;
    hist -> GetXaxis() -> SetRangeUser(0,tttttttttt);

    if (draw_refit)
    {
      auto cvsall = ejungwoo::div(ejungwoo::cv(TString("call_")+hist -> GetName(), 500,800),1,2);
      {
        cvsall -> cd(1);
        hist -> Draw();
        int ig = 0;
        for (auto idx_particle : idx_particles) {
          auto g1 = ejungwoo::gg(fit_refit,ig++);
          g1 -> SetLineColor(stp::fColor[idx_particle]);
          g1 -> Draw("samel");
        }
        fit_refit -> Draw("samel");
      }

      {
        cvsall -> cd(2) -> SetLogy();
        hist -> Draw();
        int ig = 0;
        for (auto idx_particle : idx_particles) {
          auto g1 = ejungwoo::gg(fit_refit,ig++);
          g1 -> SetLineColor(stp::fColor[idx_particle]);
          g1 -> Draw("samel");
        }
        fit_refit -> Draw("samel");
      }
    }
  }

  TF1 *refit_ampl[stp::fNumParticles] = {0};
  TF1 *refit_mean[stp::fNumParticles] = {0};
  TF1 *refit_sigm[stp::fNumParticles] = {0};

  vector<TGraphErrors *> gglist;
  vector<TF1 *> fflist;
  vector<TF1 *> fflist_mean;

  for (auto idx_particle : idx_particles)
  {
    //if (idx_particle==stp::kpip) continue;

    //TString name_particle = TString(stp::fNameShort[idx_particle]) + "_refit";
    TString name_particle = TString(stp::fNameShort[idx_particle]) + "_refit";

    TString name_amp = Form("fit_%s_ampl",name_particle.Data());
    auto gg = graph_ampl[idx_particle];
    gg -> SetMarkerColor(stp::fColor[idx_particle]);
    auto hh = ejungwoo::tohist(gg,0,2500);
    refit_ampl[idx_particle] = ejungwoo::fitgg(hh,5);
    auto ff = refit_ampl[idx_particle];
    ff -> SetName(name_amp);
    gg -> Fit(ff,"Q0");
    ff -> SetRange(ejungwoo::x1_g(gg)-100,ejungwoo::x2_g(gg)+100);
    ff -> SetLineColor(stp::fColor[idx_particle]);
    ff -> SetNpx(2000);
    ejungwoo::addto("refit_ampl",gg,"samep colorx");
    ejungwoo::addto("refit_ampl",ff,"samepl addx colorx");
    gglist.push_back(gg);
    fflist.push_back(ff);



    TString name_mean = Form("fit_%s_mean",name_particle.Data());
    refit_mean[idx_particle] = stp::f1_dedx_p5(idx_particle,p_hist_min,p_hist_max);
    ff = refit_mean[idx_particle];
    ff -> SetParameters(
        fit_mean[idx_particle] -> GetParameter(0),
        fit_mean[idx_particle] -> GetParameter(1),
        fit_mean[idx_particle] -> GetParameter(2),
        fit_mean[idx_particle] -> GetParameter(3),
        fit_mean[idx_particle] -> GetParameter(4));
    ff -> SetName(name_mean);
    gg = graph_mean[idx_particle];
    gg -> SetMarkerColor(stp::fColor[idx_particle]);
    gg -> Fit(ff,"Q0");
    if (idx_particle == stp::kp || idx_particle == stp::kd || idx_particle == stp::kt || idx_particle == stp::khe3 || idx_particle == stp::khe4)
      ff -> SetRange(ff -> GetX(1500.), ejungwoo::x2_g(gg)+500);
    else
      ff -> SetRange(ejungwoo::x1_g(gg)-10, ejungwoo::x2_g(gg)+10);
    ff -> SetLineColor(stp::fColor[idx_particle]);
    ff -> SetNpx(2000);
    ejungwoo::addto("refit_mean",gg,"samepz colorx");
    ejungwoo::addto("refit_mean",ff,"samel addx colorx");
    gglist.push_back(gg);
    fflist.push_back(ff);
    fflist_mean.push_back(ff);



    TString name_sigm = Form("fit_%s_sigm",name_particle.Data());
    refit_sigm[idx_particle] = stp::f1_dedx_p5(idx_particle,p_hist_min,p_hist_max);
    ff = refit_sigm[idx_particle];
    ff -> SetParameters(
        fit_sigm[idx_particle] -> GetParameter(0),
        fit_sigm[idx_particle] -> GetParameter(1),
        fit_sigm[idx_particle] -> GetParameter(2),
        fit_sigm[idx_particle] -> GetParameter(3),
        fit_sigm[idx_particle] -> GetParameter(4));
    ff -> SetName(name_sigm);
    gg = graph_sigm[idx_particle];
    gg -> SetMarkerColor(stp::fColor[idx_particle]);
    gg -> Fit(ff,"Q0");
    if (idx_particle == stp::kp || idx_particle == stp::kd || idx_particle == stp::kt || idx_particle == stp::khe3 || idx_particle == stp::khe4)
      ff -> SetRange(ff -> GetX(200.), ejungwoo::x2_g(gg)+500);
    else
      ff -> SetRange(ejungwoo::x1_g(gg)-10, ejungwoo::x2_g(gg)+10);
    ff -> SetLineColor(stp::fColor[idx_particle]);
    ff -> SetNpx(2000);
    ejungwoo::addto("refit_sigm",gg,"samepz colorx");
    ejungwoo::addto("refit_sigm",ff,"samel addx colorx");
    gglist.push_back(gg);
    fflist.push_back(ff);



    gg = graph_rslt[idx_particle];
    gg -> SetMarkerColor(stp::fColor[idx_particle]);
    ejungwoo::addto("refit_rslt",gg,"samepz colorx");

    TString name_rslt = Form("fit_%s_rslt",name_particle.Data());
    TString expression_rslt = Form("(%s)/(%s)",ejungwoo::expf1(refit_mean[idx_particle]).Data(),ejungwoo::expf1(refit_sigm[idx_particle]).Data());
    auto eval_rslt = new TF1(name_rslt ,expression_rslt,p_hist_min,p_hist_max);
    eval_rslt -> SetLineColor(stp::fColor[idx_particle]);
    ejungwoo::addto("refit_rslt",eval_rslt,"samel colorx addx");
  }

  ejungwoo::drawc("refit_ampl");
  ejungwoo::drawc("refit_mean");
  ejungwoo::drawc("refit_sigm");
  ejungwoo::drawc("refit_rslt");

  if (draw_pid) {
    cvs_pid -> cd();
    for (auto idx_particle : idx_particles) {
      graph_mean[idx_particle] -> Draw("samepz");
    }

    for(auto ff : fflist_mean) {
      ff -> SetLineColor(kBlack);
      ff -> Draw("samel");
    }
  }

  if (write_summary) {
    file_summary -> cd();
    hist_pid -> Write("hist_pid");
    for (auto gg : gglist) gg -> Write();
    for (auto ff : fflist) ff -> Write();
    cout << name_summary << endl;
  }
}


void get_par_limits(double p_current, int idx_particle, TH1D * hist, double &a, double &m, double &s, double &a1, double &a2, double &m1, double &m2, double &s1, double &s2)
{
       if (idx_particle == stp::kpin) { }
  else if (idx_particle == stp::kpip) { }
  else if (idx_particle == stp::kp  ) { auto s_sigm0 = 4.;  auto s_mean0 = 40.; s = s_sigm0 / s_mean0 * m; }
  else if (idx_particle == stp::kd  ) { auto s_sigm0 = 10.; auto s_mean0 = 80.; s = s_sigm0 / s_mean0 * m; }
  else if (idx_particle == stp::kt  ) { auto s_sigm0 = 25.; auto s_mean0 = 180.; s = s_sigm0 / s_mean0 * m; }
  else if (idx_particle == stp::khe3) { auto s_sigm0 = 22.; auto s_mean0 = 250.; s = s_sigm0 / s_mean0 * m; }
  else if (idx_particle == stp::khe4) { auto s_sigm0 = 35.; auto s_mean0 = 350.; s = s_sigm0 / s_mean0 * m; }
  else if (idx_particle == stp::khe6) { auto s_sigm0 = 21.; auto s_mean0 = 240.; s = s_sigm0 / s_mean0 * m; }
  else if (idx_particle == stp::kli6) { }
  else if (idx_particle == stp::kli7) { }
  else if (idx_particle == stp::kele) { }
  else if (idx_particle == stp::kpos) { }
  else if (idx_particle == stp::kli ) { auto s_sigm0 = 60.; auto s_mean0 = 400.; s = s_sigm0 / s_mean0 * m; }

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
