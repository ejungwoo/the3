const int kall = -1;

const int      kp = 0, kd = 1, kt = 2, khe3 = 3, khe4 = 4;
const int      fNumParticles = 5;
const int      fParticleIdx      [fNumParticles] = {kp, kd, kt, khe3, khe4};
const int      fParticleNumPs    [fNumParticles] = {1, 1, 1, 2, 2};
const int      fParticleNumNs    [fNumParticles] = {0, 1, 2, 1, 2};
const int      fParticleA        [fNumParticles] = {1, 2, 3, 3, 4};
const int      fParticleZ        [fNumParticles] = {1, 1, 1, 2, 2};
const TString  fParticleNames    [fNumParticles] = {"p", "d", "t", "he3", "he4"};
const TString  fParticleTitles   [fNumParticles] = {"p", "d", "t", "^{3}He", "^{4}He"};
const double   fParticlePozLLCut [fNumParticles] = {100, 100, 800, 400, 400};
const double   fParticleMass     [fNumParticles] = {938.272, 1871.06, 2809.41, 2809.41, 3728.4};
const int      fParticlePDGs     [fNumParticles] = {2212, 1000010020, 1000010030, 1000020030, 1000020040};
      double   fParticleSDHL     [fNumParticles] = {2.2, 2.0, 1.8, 1.8, 1.8};
      double   fParticleSDLL     [fNumParticles] = {2.2, 2.0, 1.8, 1.8, 1.8};

const int      k132 = 0, k108 = 1, k112 = 2, k124 = 3;
const int      fNumSystems = 4;
const int      fSystemIdx         [fNumSystems] = {k132, k108, k112, k124};
const int      fSystems           [fNumSystems] = {132, 108, 112, 124};
const int      fSystemTargets     [fNumSystems] = {124, 112, 124, 112};
const double   fSystemYAAs        [fNumSystems] = {0.3822, 0.3647, 0.3538, 0.3902};
const double   fSystemYNNs        [fNumSystems] = {0.3696, 0.3697, 0.3705, 0.3706};
const TString  fSystemNames       [fNumSystems] = {"sys132", "sys108", "sys124", "sys112"};
const TString  fSystemTitles      [fNumSystems] = {"system(132+124)", "system(108+112)", "system(124+112)", "system(112+124)"};
const TString  fSysEnergyLossFile [fNumSystems] = {"data/Sn132Sn124.txt","data/Sn108Sn112.txt","data/Sn112Sn124.txt","data/Sn112Sn124.txt"};
      Long64_t fSystemNumEvents   [fNumSystems] = {0};

const int      kf6 = 0, kx6 = 1;
const int      fNumAna = 2;
const int      fAnaIdx     [fNumAna] = {kf6, kx6};
const TString  fAnaNames   [fNumAna] = {"fix6","x0_f6set"};
const TString  fAnaONames  [fNumAna] = {"after","before"};
const TString  fAnaShort   [fNumAna] = {"f6","x6"};
const TString  fAnaTitles  [fNumAna] = {"After Fix","Before Fix"};
const TString  fAnaVersion [fNumAna] = {"NewAna.2107.4fd2bca","NewAna.2107.4fd2bca"};

const int      klr = 0, kleft = 1, kright = 2;
const int      fNumLRs = 3;
const int      fLRIdx    [fNumLRs] = {0, 1, 2};
const TString  fLRNames  [fNumLRs] = {"all","left","right"};
const TString  fLRFNames [fNumLRs] = {"*","left","right"};
const TString  fLRTitles [fNumLRs] = {"TPC-Left&Right","TPC-Left","TPC-Right"};

const int      knAll = 0, kn45 = 1, kn55 = 2;
const int      fNumMultOption = 3;
const int      fMultIdx     [fNumMultOption] = {knAll, kn45, kn55};
const TString  fMultFNames  [fNumMultOption] = {"*","45_54", "55_100"};
const TString  fMultNames   [fNumMultOption] = {"45_100","45_54", "55_100"};
const TString  fMultShort   [fNumMultOption] = {"ml", "mh", "mlh"};
const TString  fMultNames2  [fNumMultOption] = {"m45to100", "m45to54", "m55to100"};
const TString  fMultTitles  [fNumMultOption] = {"mult=45~100", "mult=45~54", "mult=55~100"};
const int      fMultLL      [fNumMultOption] = {45,  45,  55};
const int      fMultHL      [fNumMultOption] = {100, 54, 100};

const int      fNumSysComb = 4;
const int      fSysCombIdxSameTarget [] = {1,2};
const int      fSysCombIndx          [fNumSysComb] = {0,1,2,3};
const int      fSysCombIdx           [fNumSysComb][2] = {{0,1},{0,3},{2,1},{2,3}};
const TString  fSysCombNames         [fNumSysComb] = {"comb_132_108", "comb_132_112", "comb_124_108", "comb_124_112"};
const TString  fSysCombNames2        [fNumSysComb] = {"132 / 108", "132 / 112", "124 / 108", "124 / 112"};
const TString  fSysCombTitles        [fNumSysComb] = {"(132 / 108)", "(132 / 112)", "(124 / 108)", "(124 / 112)"};

const int      kttaAll = 0, ktta0 = 1, ktta20 = 2, ktta40 = 3, ktta60 = 4, kttaST60 = 5, kttaLT60 = 6;
const int      fNumCutTTAs = 7;
const int      fCutTTAIdx    [fNumCutTTAs] = {0,1,2,3,4,5,6};
const TString  fCutTTANames  [fNumCutTTAs] = {"ttaAll","ttaRG0_20","ttaRG20_40","ttaRG40_60","ttaRG60_80","ttaST60","ttaLT60"};
const TString  fCutTTATitles [fNumCutTTAs] = {"All #theta_{lab}","0<=#theta_{lab}<20","20<=#theta_{lab}<40","40<=#theta_{lab}<60","60<=#theta_{lab}<80","#theta_{lab}<60","#theta_{lab}>=60"};
const TCut     fCutTTAValues [fNumCutTTAs] = {
  "",
  "theta_lab>= 0*TMath::DegToRad() && theta_lab<20*TMath::DegToRad()",
  "theta_lab>=20*TMath::DegToRad() && theta_lab<40*TMath::DegToRad()",
  "theta_lab>=40*TMath::DegToRad() && theta_lab<60*TMath::DegToRad()",
  "theta_lab>=60*TMath::DegToRad() && theta_lab<80*TMath::DegToRad()",
  "theta_lab<60*TMath::DegToRad()",
  "theta_lab>=60*TMath::DegToRad()",
};
