const int      kall = -1;

const int      kp = 0, kd = 1, kt = 2, khe3 = 3, khe4 = 4;
const int      fNumParticles = 5;
const int      fParticleIdx      [fNumParticles] = {kp, kd, kt, khe3, khe4};
const int      fParticleNumPs    [fNumParticles] = {1, 1, 1, 2, 2};
const int      fParticleNumNs    [fNumParticles] = {0, 1, 2, 1, 2};
const int      fParticleA        [fNumParticles] = {1, 2, 3, 3, 4};
const int      fParticleZ        [fNumParticles] = {1, 1, 1, 2, 2};
const char*    fParticleNames    []              = {"p", "d", "t", "he3", "he4", ""};
const char*    fParticleTitles   []              = {"p", "d", "t", "^{3}He", "^{4}He", ""};
const double   fParticlePozLLCut [fNumParticles] = {100, 100, 800, 400, 400};
const double   fParticleMass     [fNumParticles] = {938.272, 1871.06, 2809.41, 2809.41, 3728.4};
const int      fParticlePDGs     [fNumParticles] = {2212, 1000010020, 1000010030, 1000020030, 1000020040};
const double   fParticleSDHL     [fNumParticles] = {2.2, 2.0, 1.8, 1.8, 1.8};
const double   fParticleSDLL     [fNumParticles] = {2.2, 2.0, 1.8, 1.8, 1.8};

const int      k132 = 0, k108 = 1, k112 = 2, k124 = 3;
const int      fNumSystems = 4;
const int      fSystemIdx         [fNumSystems] = {k132, k108, k112, k124};
const int      fSystems           [fNumSystems] = {132, 108, 112, 124};
const int      fSystemTargets     [fNumSystems] = {124, 112, 124, 112};
const double   fSystemYAAs        [fNumSystems] = {0.3822, 0.3647, 0.3538, 0.3902};
const double   fSystemYNNs        [fNumSystems] = {0.3696, 0.3697, 0.3705, 0.3706};
const char*    fSystemNames       [fNumSystems] = {"132", "108", "124", "112"};
const char*    fSystemTitles      [fNumSystems] = {"system(132+124)", "system(108+112)", "system(124+112)", "system(112+124)"};
const char*    fSysEnergyLossFile [fNumSystems] = {"data/Sn132Sn124.txt","data/Sn108Sn112.txt","data/Sn112Sn124.txt","data/Sn112Sn124.txt"};
int            fSystemNumEvents   [fNumSystems] = {0};

const int      fNumSysComb = 4;
const int      fSysCombIdxTest       [] = {0};
const int      fSysCombIdxSameTarget [] = {1,2};
const int      fSysCombIndx          [fNumSysComb] = {0,1,2,3};
const int      fSysCombIdx           [fNumSysComb][2] = {{0,1},{0,3},{2,1},{2,3}};
const char*    fSysCombNames         [fNumSysComb] = {"comb_132_108", "comb_132_112", "comb_124_108", "comb_124_112"};
const char*    fSysCombNames2        [fNumSysComb] = {"132 / 108", "132 / 112", "124 / 108", "124 / 112"};
const char*    fSysCombTitles        [fNumSysComb] = {"(132 / 108)", "(132 / 112)", "(124 / 108)", "(124 / 112)"};

const int      kf7 = 0, kf6 = 1, kx6 = 2;
const int      fNumAna = 3;
const int      fAnaIdx     [fNumAna] = {kf7, kf6, kx6};
const char*    fAnaFNames  [fNumAna] = {"f7", "fix6","x0_f6set"};
const char*    fAnaNames   [fNumAna] = {"f7", "f","x"};
const char*    fAnaONames  [fNumAna] = {"f7", "after","before"};
const char*    fAnaShort   [fNumAna] = {"f7", "f6","x6"};
const char*    fAnaTitles  [fNumAna] = {"After Fix", "After Fix","Before Fix"};
const char*    fAnaVersion [fNumAna] = {"NewAna.2107.4fd2bca","NewAna.2107.4fd2bca","NewAna.2107.4fd2bca"};

const int      klr = 0, kleft = 1, kright = 2;
const int      fNumLRs = 3;
const int      fLRIdx    [fNumLRs] = {0, 1, 2};
const char*    fLRNames  [fNumLRs] = {"all","left","right"};
const char*    fLRFNames [fNumLRs] = {"*","left","right"};
const char*    fLRTitles [fNumLRs] = {"TPC-All","TPC-Left","TPC-Right"};

const int      knAll = 0, kn45 = 1, kn55 = 2;
const int      fNumMultOption = 3;
//const int      fMultIdx     [fNumMultOption] = {knAll, kn45, kn55};
const int      fMultIdx     [fNumMultOption] = {kn45, kn55};
const char*    fMultFNames  [fNumMultOption] = {"*","45_54", "55_100"};
const char*    fMultNames   [fNumMultOption] = {"45100","4554", "55100"};
const char*    fMultShort   [fNumMultOption] = {"ml", "mh", "mlh"};
const char*    fMultNames2  [fNumMultOption] = {"m45to100", "m45to54", "m55to100"};
const char*    fMultTitles  [fNumMultOption] = {"mult=45~100", "mult=45~54", "mult=55~100"};
const int      fMultLL      [fNumMultOption] = {45,  45,  55};
const int      fMultHL      [fNumMultOption] = {100, 54, 100};

const int      kttaAll = 0, ktta0 = 1, ktta20 = 2, ktta40 = 3, ktta60 = 4, kttaLT60 = 5, kttaGT60 = 6;
const int      fNumCutTTAs = 7;
const int      fCutTTAIdx    [fNumCutTTAs] = {0,1,2,3,4,5,6};
const char*    fCutTTANames  [fNumCutTTAs] = {"ttaAll","ttaRG0_20","ttaRG20_40","ttaRG40_60","ttaRG60_80","ttaLT60","ttaGT60"};
const char*    fCutTTATitles [fNumCutTTAs] = {"All-#theta_{lab}","0<=#theta_{lab}<20","20<=#theta_{lab}<40","40<=#theta_{lab}<60","60<=#theta_{lab}<80","#theta_{lab}<60","#theta_{lab}>=60"};
const TCut     fCutTTAValues [fNumCutTTAs] = {
  "",
  "theta_lab>= 0*TMath::DegToRad()&&theta_lab<20*TMath::DegToRad()",
  "theta_lab>=20*TMath::DegToRad()&&theta_lab<40*TMath::DegToRad()",
  "theta_lab>=40*TMath::DegToRad()&&theta_lab<60*TMath::DegToRad()",
  "theta_lab>=60*TMath::DegToRad()&&theta_lab<80*TMath::DegToRad()",
  "theta_lab <60*TMath::DegToRad()",
  "theta_lab>=60*TMath::DegToRad()",
};

const int      fCutTTA0Idx    [] = {0,1,2,3,4,5,6,7,8};
const TCut     fCutTTA0Values [] = {
  TCut("theta_lab>=(0)*TMath::DegToRad()&&theta_lab<(1)*TMath::DegToRad()"),
  TCut("theta_lab>=(10-1)*TMath::DegToRad()&&theta_lab<(10+1)*TMath::DegToRad()"),
  TCut("theta_lab>=(20-1)*TMath::DegToRad()&&theta_lab<(20+1)*TMath::DegToRad()"),
  TCut("theta_lab>=(30-1)*TMath::DegToRad()&&theta_lab<(30+1)*TMath::DegToRad()"),
  TCut("theta_lab>=(40-1)*TMath::DegToRad()&&theta_lab<(40+1)*TMath::DegToRad()"),
  TCut("theta_lab>=(50-1)*TMath::DegToRad()&&theta_lab<(50+1)*TMath::DegToRad()"),
  TCut("theta_lab>=(60-1)*TMath::DegToRad()&&theta_lab<(60+1)*TMath::DegToRad()"),
  TCut("theta_lab>=(70-1)*TMath::DegToRad()&&theta_lab<(70+1)*TMath::DegToRad()"),
  TCut("theta_lab>=(80-1)*TMath::DegToRad()&&theta_lab<(80+1)*TMath::DegToRad()"),
};

const int      kya = 0, ky02 = 1, ky0 = 2, ky04 = 3, ky0610 = 4;
const int      fNumCutY0s = 5;
const int      fCutY0Idx[] = {kya, ky02, ky0, ky04, ky0610};
const char*    fCutY0Names[] = {"yAll","yGT02","yGT0","y0004","y0610"};
const char*    fCutY0Titles[] = {"", "y_{0}>0.2", "y_{0}>0", "0<y_{0}<.4", ".6<y_{0}<1."};
const TCut     fCutY0Values[] = {
  "",
  "fy_cm/(by_cm/2)>0.2",
  "fy_cm/(by_cm/2)>0",
  "fy_cm/(by_cm/2)>0&&fy_cm/(by_cm/2)<0.4",
  "fy_cm/(by_cm/2)>0.6&&fy_cm/(by_cm/2)<1"
};
