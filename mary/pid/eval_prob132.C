#include "stp.h"

void eval_prob132()
{
  ejungwoo::gstat(0);
  ejungwoo::gcvspos(1300,0);

  ejungwoo::gversion("eval_prob132");

  int idx_particles[] = { stp::kpip, stp::kp, stp::kd, stp::kt, stp::khe3, stp::khe4, stp::khe6, stp::kli, };
  //int idx_particles[] = { stp::kp, stp::kd, stp::kt, stp::khe3, stp::khe4, stp::khe6, stp::kli, };

  auto file = new TFile("/Users/ejungwoo/spirit/the3/mary/pid/data__sys132/summary.sys132.root");
  file -> ls();

  auto hist = (TH2D *) file -> Get("hist_pid");
  auto cvs = ejungwoo::cc4("pid");
  cvs -> SetLogz();
  hist -> Draw("colz");

  TF1 *fit_dedx[3][20];
  for (auto idx_particle : idx_particles)
  {
    TString name_particle = stp::fNameShort[idx_particle];
    cout << TString("fit_")+name_particle+"_refit_ampl" << endl;
    cout << TString("fit_")+name_particle+"_refit_mean" << endl;
    cout << TString("fit_")+name_particle+"_refit_sigm" << endl;
    fit_dedx[0][idx_particle] = (TF1 *) file -> Get(TString("fit_")+name_particle+"_refit_ampl");
    fit_dedx[1][idx_particle] = (TF1 *) file -> Get(TString("fit_")+name_particle+"_refit_mean");
    fit_dedx[2][idx_particle] = (TF1 *) file -> Get(TString("fit_")+name_particle+"_refit_sigm");

    //fit_dedx[0][idx_particle] -> SetRange(10,3000);
    //fit_dedx[1][idx_particle] -> SetRange(10,3000);
    //fit_dedx[2][idx_particle] -> SetRange(10,3000);

    ejungwoo::addto("ampl",fit_dedx[0][idx_particle],"l",name_particle);
    ejungwoo::addto("mean",fit_dedx[1][idx_particle],"l",name_particle);
    ejungwoo::addto("sigm",fit_dedx[2][idx_particle],"l",name_particle);

    fit_dedx[1][idx_particle] -> Draw("samel");
  }

  ejungwoo::save(cvs,"png");

  ejungwoo::save(ejungwoo::drawc("ampl"));
  //ejungwoo::save(ejungwoo::drawc("mean"));
  //ejungwoo::save(ejungwoo::drawc("sigm"));

  // dummy run
  for (auto idx_particle : idx_particles) {
    auto amp = fit_dedx[0][idx_particle] -> Eval(50);
    auto mean = fit_dedx[1][idx_particle] -> Eval(50);
    auto sigma = fit_dedx[2][idx_particle] -> Eval(50);
  }

  auto s0 = 400.;
  for (auto p0 = 10.; p0<4000; p0+=200.)
  {
    int idxMax = 0;
    double valMax = 0.;
    double valSum = 0.;

    for (auto idx_particle : idx_particles)
    {
      auto amp = fit_dedx[0][idx_particle] -> Eval(p0);
      auto mean = fit_dedx[1][idx_particle] -> Eval(p0);
      auto sigma = fit_dedx[2][idx_particle] -> Eval(p0);

      //cout << stp::fNameShort[idx_particle] << " " << amp << " " << mean << " " << sigma << endl;

      auto val = amp * TMath::Gaus(s0,mean,sigma);

      valSum += val;
      if (val > valMax) {
        valMax = val;
        idxMax = idx_particle;
      }
    }

    cout << p0 << " " << stp::fNameShort[idxMax] << " " << 100*valMax/valSum << endl;
  }
}
