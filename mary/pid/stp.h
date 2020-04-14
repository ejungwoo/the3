#ifndef SPIRIT_PARTICLE_CONFIGURATION
#define SPIRIT_PARTICLE_CONFIGURATION

namespace stp/*article_configuration*/
{
  const int fNumParticles = 13;
  const int fNumZs = 3;

  const int kpin = 0;
  const int kpip = 1;
  const int kp   = 2;
  const int kd   = 3;
  const int kt   = 4;
  const int khe3 = 5;
  const int khe4 = 6;
  const int khe6 = 7;
  const int kli6 = 8;
  const int kli7 = 9;
  const int kele = 10;
  const int kpos = 11;
  const int kli  = 12;
  const int kall = 99;

  // MeV/c^2
  const double kmpi  = 139.57018; //
  const double kmpip = 139.57018;
  const double kmpin = 139.57018;
  const double kmp   = 938.2720813;
  const double kmd   = 1875.612762;
  const double kmt   = 2808.921112;
  const double kmhe3 = 2808.39132;
  const double kmhe4 = 3727.379378;
  const double kmhe6 = 5606.556709674;
  const double kmli6 = 5603.051494945;
  const double kmli7 = 6535.365824124;
  const double kmele = 0.5109989461;
  const double kmpos = 0.5109989461;
  const double kmli  = kmli6;

  const char *fNameShort[fNumParticles] = {"pin","pip","p","d","t","he3","he4","he6","li6","li7","ele","pos","li"};
  const char *fNameLong[fNumParticles] = {"pin","pip","proton","deuteron","triton","he3","he4","he6","li6","li7","electron","positron","li"};
  const int fZParticle[fNumParticles] = {-1,1,1,1,1,2,2,2,3,3,1,-1,3};
  const double fMass[fNumParticles] = { kmpip,kmpin,kmp,kmd,kmt,kmhe3,kmhe4,kmhe6,kmli6,kmli7,kmele,kmpos,kmli};
  int fColor[fNumParticles] = {kBlack,kOrange+8,kOrange,kSpring+7,kSpring-6,kAzure+7,kBlue,kViolet-5,kPink+7,kGray+1,kGray+1,kGray+1,kPink+7};
  int fMarkerStyle[fNumParticles] = {42,42,24,25,32,28,30,5,42,44,46,46,42};

  const double klogI_p10 = -8.66662;
  const double kzova_p10 = 0.457446809; // = 17.2 / 37.6 (Z/A medium)

  double bethebloch_dedx(double z_reco, double m_particle, double *x, double *par);
  double parameter5_dedx(double z_reco, double m_particle, double *x, double *par);
  double mostprob_dedx(double z_reco, double m_particle, double *x, double *par);
  double bethebloch_dedx_diff(double *x, double *par);

  TString expression_bb(double z_reco, double m_particle);
  TString expression_p5(double z_reco, double m_particle);

  TF1 *f1_dedx_bb(int idx_particle, double p1=0, double p2=3000);
  TF1 *f1_dedx_p5(int idx_particle, double p1=0, double p2=3000);

  double p5_dedx_pin (double *x, double *par) { return parameter5_dedx(-1, kmpin, x, par); }
  double p5_dedx_pip (double *x, double *par) { return parameter5_dedx( 1, kmpip, x, par); }
  double p5_dedx_p   (double *x, double *par) { return parameter5_dedx( 1, kmp,   x, par); }
  double p5_dedx_d   (double *x, double *par) { return parameter5_dedx( 1, kmd,   x, par); }
  double p5_dedx_t   (double *x, double *par) { return parameter5_dedx( 1, kmt,   x, par); }
  double p5_dedx_he3 (double *x, double *par) { return parameter5_dedx( 2, kmhe3, x, par); }
  double p5_dedx_he4 (double *x, double *par) { return parameter5_dedx( 2, kmhe4, x, par); }
  double p5_dedx_he6 (double *x, double *par) { return parameter5_dedx( 2, kmhe6, x, par); }
  double p5_dedx_li6 (double *x, double *par) { return parameter5_dedx( 3, kmli6, x, par); }
  double p5_dedx_li7 (double *x, double *par) { return parameter5_dedx( 3, kmli7, x, par); }
  double p5_dedx_ele (double *x, double *par) { return parameter5_dedx( 1, kmele, x, par); }
  double p5_dedx_pos (double *x, double *par) { return parameter5_dedx(-1, kmpos, x, par); }
  double p5_dedx_li  (double *x, double *par) { return parameter5_dedx( 3, kmli , x, par); }
};

/////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////
TString stp::expression_bb(double z_reco, double m_particle) {
  auto val1 = z_reco*z_reco*kzova_p10;
  auto val2 = m_particle/z_reco;
  auto val3 = 4*kmele*m_particle*m_particle/(z_reco*z_reco);
  auto val4 = klogI_p10;
  auto expression_dedx = Form("[0]*%.4e*(1.+((%.4e)/x)**2)*(.5*log(%.4e/(x**2))-%.4e-1./(1.+((%.4e)/x)**2)-[1])", val1, val2, val3, val4, val2);

  return expression_dedx;
}

/////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////
TString stp::expression_p5(double z_reco, double m_particle) {
  TString beta_square = Form("1./(1+(%f/(x*%f))**2)",m_particle,z_reco);
  TString beta_gamma = Form("%f/(x*%f)",m_particle,z_reco);
  TString expression_beta4 = Form("TMath::Power(%s,[3])",beta_square.Data());
  TString expression_beta5 = Form("TMath::Power(%s,2*[4])",beta_gamma.Data());
  TString expression_dedx = Form("[0] / %s * ([1] - %s - log([2] + %s))",expression_beta4.Data(),expression_beta4.Data(),expression_beta5.Data());
  return expression_dedx;
}

/////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////
TF1 *stp::f1_dedx_bb(int idx_particle, double p1, double p2)
{
  if (idx_particle >=0 && idx_particle < fNumParticles)
    return new TF1(Form("f1_dedx_bb_%s",stp::fNameShort[idx_particle]),expression_bb(fZParticle[idx_particle], fMass[idx_particle]), p1, p2);
  return (TF1 *) nullptr;
}

/////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////
TF1 *stp::f1_dedx_p5(int idx_particle, double p1, double p2)
{
  if (idx_particle >=0 && idx_particle < fNumParticles)
    return new TF1(Form("f1_dedx_p5_%s",stp::fNameShort[idx_particle]),expression_p5(fZParticle[idx_particle], fMass[idx_particle]), p1, p2);
  return (TF1 *) nullptr;
}

/////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////
double stp::bethebloch_dedx(double z_reco, double m_particle, double *x, double *par) {
  auto p_reco = x[0]; // momentum reco
  auto fit_par1 = par[0]; // fit parameter 1
  auto fit_par2 = par[1]; // fit parameter 2
  auto val1 = z_reco*z_reco*kzova_p10;
  auto val2 = m_particle/z_reco;
  auto val3 = 4*kmele*m_particle*m_particle/(z_reco*z_reco);
  auto val4 = klogI_p10;
  auto dedx = fit_par1*val1 * (1.+(val2/p_reco)*(val2/p_reco)) * (0.5*log(val3/(p_reco*p_reco)) - klogI_p10 - 1./(1.+(val2/p_reco)*(val2/p_reco)) - fit_par2);
  return dedx;
}

/////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////
double stp::parameter5_dedx(double z_reco, double m_particle, double *x, double *par) {
  auto pz_reco = x[0]*z_reco;
  auto beta_square = pz_reco*pz_reco/(pz_reco*pz_reco+m_particle*m_particle);
  auto gamma_square = 1./(1.-beta_square);
  auto beta4 = TMath::Power(beta_square,par[3]);
  auto beta5 = TMath::Power(1/(beta_square*gamma_square),par[4]);

  auto dedx = (par[0]) / beta4 * (par[1] - beta4 - log(par[2] + beta5));
  return dedx;
}

/////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////
double stp::mostprob_dedx(double z_reco, double m_particle, double *x, double *par) {
  auto p_reco = x[0]; // momentum reco
  auto fit_par1 = par[0]; // fit parameter 1
  auto fit_par2 = par[1]; // fit parameter 2

  auto pz_reco = p_reco*z_reco;
  auto beta_square = pz_reco*pz_reco/(pz_reco*pz_reco+m_particle*m_particle);
  auto gamma_square = 1./(1.-beta_square);

  auto density = 14.422205; // g/cm^-2
  auto xi = (fit_par1/2)*kzova_p10*z_reco*z_reco*(density/beta_square);
  auto dedx = -xi * (TMath::Log(2*m_particle*beta_square*gamma_square) - klogI_p10 + TMath::Log(xi) - klogI_p10 + 0.02 - beta_square + fit_par2);
  return dedx;
}

/////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////
double stp::bethebloch_dedx_diff(double *x, double *par) {
  double m_reco = x[2];
  double mom[] = {par[2]};
  auto z_reco = par[3]; // z reco
  auto dedx_reco = par[4];

  auto dedx_bb = bethebloch_dedx(z_reco, m_reco, mom, par);
  auto ddedx = dedx_reco - dedx_bb;
  return ddedx;
}

#endif
