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
  const double kmpi  = 139.57018;

  const int kzpin = -1;
  const int kzpip =  1;
  const int kzp   =  1;
  const int kzd   =  1;
  const int kzt   =  1;
  const int kzhe3 =  2;
  const int kzhe4 =  2;
  const int kzhe6 =  2;
  const int kzli6 =  3;
  const int kzli7 =  3;
  const int kzli  =  3;
  const int kzele =  1;
  const int kzpos = -1;

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
  const double kmli  = kmli6;
  const double kmele = 0.5109989461;
  const double kmpos = 0.5109989461;

  const double klogI_p10 = -8.66662;
  const double kzova_p10 = 0.457446809; // = 17.2 / 37.6 (Z/A medium)

  TString expression_bethebloch_dedx(double z_reco, double m_particle) {
    auto val1 = z_reco*z_reco*kzova_p10;
    auto val2 = m_particle/z_reco;
    auto val3 = 4*kmele*m_particle*m_particle/(z_reco*z_reco);
    auto val4 = klogI_p10;
    auto expression_dedx = Form("[0]*%.4e*(1.+((%.4e)/x)**2)*(.5*log(%.4e/(x**2))-%.4e-1./(1.+((%.4e)/x)**2)-[1])", val1, val2, val3, val4, val2);

    return expression_dedx;
  }

  double bethebloch_dedx(double z_reco, double m_particle, double *x, double *par) {
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

  double parameter5_dedx(double z_reco, double m_particle, double *x, double *par) {

    //auto bbbb = (m_particle/z_reco/x[0]);
    //auto beta = sqrt(1./(1.-(bbbb*bbbb)));
    //auto gamma = sqrt(1./(1.-beta*beta));
    //auto beta4 = TMath::Power(beta,par[3]);
    //auto beta5 = TMath::Power(1/(beta*gamma),par[4]);

    auto pz_reco = x[0]*z_reco;
    auto beta_square = pz_reco*pz_reco/(pz_reco*pz_reco+m_particle*m_particle);
    auto gamma_square = 1./(1.-beta_square);
    auto beta4 = TMath::Power(beta_square,par[3]);
    auto beta5 = TMath::Power(1/(beta_square*gamma_square),par[4]);

    auto dedx = (par[0]) / beta4 * (par[1] - beta4 - log(par[2] + beta5));
    return dedx;
  }

  /*
  double bethebloch_dedx(double z_reco, double m_particle, double *x, double *par) {
    auto p_reco = x[0]; // momentum reco
    auto fit_par1 = par[0]; // fit parameter 1
    auto fit_par2 = par[1]; // fit parameter 2

    //auto beta_square = 1./(1.-((m_particle/z_reco/p_reco)**2);
    auto pz_reco = p_reco*z_reco;
    auto beta_square = pz_reco*pz_reco/(pz_reco*pz_reco+m_particle*m_particle);
    auto gamma_square = 1./(1.-beta_square);

    auto beta_gamma_square = (m_particle/p_reco/z_reco);

    auto w_max1 = 2*kmele*((p_reco*z_reco/m_particle)*(p_reco*z_reco/m_particle));
    //auto w_max2 = (1 + 2*kmele*sqrt(gamma_square)/m_particle + (kmele/m_particle)*(kmele/m_particle));
    auto w_max2 = 1;
    auto w_max = w_max1 / w_max2;

    //auto dedx = fit_par1*z_reco*z_reco*kzova_p10/beta_square*(0.5*log(2*kmele*beta_square*gamma_square*w_max)-klogI_p10-beta_square-fit_par2);
    auto dedx = fit_par1*z_reco*z_reco*kzova_p10/beta_square*(0.5*log(2*kmele*beta_gamma_square*w_max)-klogI_p10-beta_square-fit_par2);
    return dedx;
  }
  */

  double mostprob_dedx(double z_reco, double m_particle, double *x, double *par) {
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

  TString expression_bb_dedx_pin() { return expression_bethebloch_dedx(-1, kmpin); }
  TString expression_bb_dedx_pip() { return expression_bethebloch_dedx( 1, kmpip); }
  TString expression_bb_dedx_p()   { return expression_bethebloch_dedx( 1, kmp  ); }
  TString expression_bb_dedx_d()   { return expression_bethebloch_dedx( 1, kmd  ); }
  TString expression_bb_dedx_t()   { return expression_bethebloch_dedx( 1, kmt  ); }
  TString expression_bb_dedx_he3() { return expression_bethebloch_dedx( 2, kmhe3); }
  TString expression_bb_dedx_he4() { return expression_bethebloch_dedx( 2, kmhe4); }
  TString expression_bb_dedx_he6() { return expression_bethebloch_dedx( 2, kmhe6); }
  TString expression_bb_dedx_li6() { return expression_bethebloch_dedx( 3, kmli6); }
  TString expression_bb_dedx_li7() { return expression_bethebloch_dedx( 3, kmli7); }
  TString expression_bb_dedx_li () { return expression_bethebloch_dedx( 3, kmli ); }
  TString expression_bb_dedx_ele() { return expression_bethebloch_dedx( 1, kmele); }
  TString expression_bb_dedx_pos() { return expression_bethebloch_dedx(-1, kmpos); }

  double bb_dedx_pin (double *x, double *par) { return bethebloch_dedx(-1, kmpin, x, par); }
  double bb_dedx_pip (double *x, double *par) { return bethebloch_dedx( 1, kmpip, x, par); }
  double bb_dedx_p   (double *x, double *par) { return bethebloch_dedx( 1, kmp,   x, par); }
  double bb_dedx_d   (double *x, double *par) { return bethebloch_dedx( 1, kmd,   x, par); }
  double bb_dedx_t   (double *x, double *par) { return bethebloch_dedx( 1, kmt,   x, par); }
  double bb_dedx_he3 (double *x, double *par) { return bethebloch_dedx( 2, kmhe3, x, par); }
  double bb_dedx_he4 (double *x, double *par) { return bethebloch_dedx( 2, kmhe4, x, par); }
  double bb_dedx_he6 (double *x, double *par) { return bethebloch_dedx( 2, kmhe6, x, par); }
  double bb_dedx_li6 (double *x, double *par) { return bethebloch_dedx( 3, kmli6, x, par); }
  double bb_dedx_li7 (double *x, double *par) { return bethebloch_dedx( 3, kmli7, x, par); }
  double bb_dedx_li  (double *x, double *par) { return bethebloch_dedx( 3, kmli , x, par); }
  double bb_dedx_ele (double *x, double *par) { return bethebloch_dedx( 1, kmele, x, par); }
  double bb_dedx_pos (double *x, double *par) { return bethebloch_dedx(-1, kmpos, x, par); }

  double mp_dedx_pin (double *x, double *par) { return mostprob_dedx(-1, kmpin, x, par); }
  double mp_dedx_pip (double *x, double *par) { return mostprob_dedx( 1, kmpip, x, par); }
  double mp_dedx_p   (double *x, double *par) { return mostprob_dedx( 1, kmp,   x, par); }
  double mp_dedx_d   (double *x, double *par) { return mostprob_dedx( 1, kmd,   x, par); }
  double mp_dedx_t   (double *x, double *par) { return mostprob_dedx( 1, kmt,   x, par); }
  double mp_dedx_he3 (double *x, double *par) { return mostprob_dedx( 2, kmhe3, x, par); }
  double mp_dedx_he4 (double *x, double *par) { return mostprob_dedx( 2, kmhe4, x, par); }
  double mp_dedx_he6 (double *x, double *par) { return mostprob_dedx( 2, kmhe6, x, par); }
  double mp_dedx_li6 (double *x, double *par) { return mostprob_dedx( 3, kmli6, x, par); }
  double mp_dedx_li7 (double *x, double *par) { return mostprob_dedx( 3, kmli7, x, par); }
  double mp_dedx_li  (double *x, double *par) { return mostprob_dedx( 3, kmli , x, par); }
  double mp_dedx_ele (double *x, double *par) { return mostprob_dedx( 1, kmele, x, par); }
  double mp_dedx_pos (double *x, double *par) { return mostprob_dedx(-1, kmpos, x, par); }

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
  double p5_dedx_li  (double *x, double *par) { return parameter5_dedx( 3, kmli , x, par); }
  double p5_dedx_ele (double *x, double *par) { return parameter5_dedx( 1, kmele, x, par); }
  double p5_dedx_pos (double *x, double *par) { return parameter5_dedx(-1, kmpos, x, par); }

  double bethebloch_dedx_diff(double *x, double *par) {
    double m_reco = x[2];
    double mom[] = {par[2]};
    auto z_reco = par[3]; // z reco
    auto dedx_reco = par[4];

    auto dedx_bb = bethebloch_dedx(z_reco, m_reco, mom, par);
    auto ddedx = dedx_reco - dedx_bb;
    return ddedx;
  }

  const double fMass[fNumParticles] = {
    kmpip,
    kmpin,
    kmp  ,
    kmd  ,
    kmt  ,
    kmhe3,
    kmhe4,
    kmhe6,
    kmli6,
    kmli7,
    kmele,
    kmpos};

  const int fZParticle[fNumParticles] = {
    kzpin,
    kzpip,
    kzp  ,
    kzd  ,
    kzt  ,
    kzhe3,
    kzhe4,
    kzhe6,
    kzli6,
    kzli7,
    kzele,
    kzpos,};

  const char *fNameLong[fNumParticles] = {
    /*pi*/ "pin", "pip",
    /*Z1*/ "proton", "deuteron", "triton",
    /*Z2*/ "he3", "he4", "he6",
    /*Z3*/ "li6", "li7",
    /*el*/ "electron", "positron", "li"};

  const char *fNameShort[fNumParticles] = {
    /*pi*/ "pin", "pip",
    /*Z1*/ "p"  , "d"  , "t", 
    /*Z2*/ "he3", "he4", "he6",
    /*Z3*/ "li6", "li7",
    /*el*/ "ele", "pos", "li"};

  int fMarkerStyle[fNumParticles] = {
    /*pi*/ 42,42,
    /*Z1*/ 24,25,32,
    /*Z2*/ 28,30,5,
    /*Z3*/ 42,44,
    /*el*/ 46,46};

  int fColor[fNumParticles] = {
    /*pi*/ kBlack, kBlack,
    /*Z1*/ kRed+1, kOrange+7, kOrange-2,
    /*Z2*/ kBlue+1, kAzure-4, kCyan+1,
    /*Z3*/ kGreen+2, kGreen,
    /*el*/ kGray+1, kGray+1};

  Double_t fSigmaInitParticles900[] = {
    /*pi*/ 2,2,
    /*Z1*/ 3,5,5,
    /*Z2*/ 12,12,30,
    /*Z3*/ 60,60,
    /*el*/ 40,40};

  Double_t fSigmaInitParticles2000[] = {
    /*pi*/ 2,2,
    /*Z1*/ 3,5,5,
    /*Z2*/ 10,10,30,
    /*Z3*/ 60,60,
    /*el*/ 40,40};

  int fZIdx[] = {
    /*pi*/ 0,0,
    /*Z1*/ 0,0,0,
    /*Z2*/ 1,1,
    /*he6*/2,
    /*Z3*/ 2,2,
    /*el*/ 0,0}; 

  double fPLowCut[] = {
    /*pi*/ 50,50,
    /*Z1*/ 200,400,
           600,
    /*Z2*/ 600,600,
    /*he6*/1000,
    /*Z3*/ 1000,1000,
    /*el*/ 0,0
  };

  const int kSpl3TestG = 0;
  const int kTrueGuide = 1;
  const int kTestGuide = 2;
  const int kHandGuide = 3;
  const int kSpl3TrueG = 4;

  const char *fNameMethods[5] = {"Spl3-Guideline", "True-Guideline", "Test-Guideline", "Hand-Guideline", "Spl3-Guideline"};

  const int fNumFuncPars = 4;
  const char *fNamePar[fNumFuncPars] = {"Amplitude","Mean","Sigma","SLength"};

  const double fParMin[fNumZs][fNumFuncPars]  = {{        5,  40,  2,100}, {        5,  40,  2,100},{        5,  40,  2,100}};
  const double fParMax[fNumZs][fNumFuncPars]  = {{100000000,1500,250,100}, {100000000,1500, 80,100},{100000000,1500,280,100}};
  const double fDParMin[fNumZs][fNumFuncPars] = {{     1000,  10,  5,  1}, {     1000,  10,  5,  1},{      50,   10,  5,  1}};
  const double fDParMax[fNumZs][fNumFuncPars] = {{  1000000, 100, 30,  1}, {  1000000, 100, 30,  1},{  1000000, 100, 30,  1}};

  double fMaxD[fNumZs][3] = {{0,0,0},{60,60,60},{20,20,20}};

  double fMeanBR = 1;
  double fMeanRR = 3.;
  double fParC = 1.;

  double fEntryLow = 40;

  double fCutSignalToNoise = 1;
};

#endif
