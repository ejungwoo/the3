const int      kall = -1;

const int      kp = 0, kd = 1, kt = 2, khe3 = 3, khe4 = 4;
const int      fNumParticles = 5;
const int      fParticleIdx      [fNumParticles] = {kp, kd, kt, khe3, khe4};
const int      fParticleNumPs    [fNumParticles] = {1, 1, 1, 2, 2};
const int      fParticleN        [fNumParticles] = {0, 1, 2, 1, 2};
const int      fParticleA        [fNumParticles] = {1, 2, 3, 3, 4};
const int      fParticleZ        [fNumParticles] = {1, 1, 1, 2, 2};
const char*    fParticleNames0   []              = {"p", "d", "t", "he3", "he4"};
const char*    fParticleNames    []              = {"p", "d", "t", "he3", "he4", ""};
const char*    fParticleTitles   []              = {"p", "d", "t", "^{3}He", "^{4}He", ""};
const double   fParticlePozLLCut [fNumParticles] = {100, 100, 800, 400, 400};
const double   fParticleMass     [fNumParticles] = {938.272, 1871.06, 2809.41, 2809.41, 3728.4};
const int      fParticlePDGs     [fNumParticles] = {2212, 1000010020, 1000010030, 1000020030, 1000020040};
const double   fParticleSDHL     [fNumParticles] = {2.2, 2.0, 1.8, 1.8, 1.8};
const double   fParticleSDLL     [fNumParticles] = {2.2, 2.0, 1.8, 1.8, 1.8};
const TCut     fParticlePozCut   [fNumParticles] = {
  "p_lab>100",
  "p_lab>200",
  "p_lab>400",
  "p_lab>400",
  "p_lab>400",
};

const int      k132 = 0, k108 = 1, k112 = 2, k124 = 3;
const int      fNumSystems = 4;
const int      fSystemIdx         [fNumSystems] = {k132, k108, k112, k124};
const int      fSystems           [fNumSystems] = {132, 108, 112, 124};
const int      fSystemTargets     [fNumSystems] = {124, 112, 124, 112};
const double   fSystemYAAs        [fNumSystems] = {0.3822, 0.3647, 0.3538, 0.3902};
const double   fSystemYNNs        [fNumSystems] = {0.3696, 0.3697, 0.3705, 0.3706};
const char*    fSystemNames       [fNumSystems] = {"132", "108", "112", "124"};
const char*    fSystemTitles      [fNumSystems] = {"system(132+124)", "system(108+112)", "system(124+112)", "system(112+124)"};
const char*    fSysEnergyLossFile [fNumSystems] = {"data/Sn132Sn124.txt","data/Sn108Sn112.txt","data/Sn112Sn124.txt","data/Sn112Sn124.txt"};
int            fSystemNumEvents   [fNumSystems] = {0};

const int      fNumSysComb = 4;
const int      fSysCombIdxTest       [] = {0};
const int      fSysCombIdxSameTarget [] = {1,2};
const int      fSysCombIndx          [fNumSysComb] = {0,1,2,3};
const int      fSysCombIdx           [fNumSysComb][2] = {{0,1},{0,2},{3,1},{2,3}};
const char*    fSysCombNames         [fNumSysComb] = {"comb_132_108", "comb_132_112", "comb_124_108",   "comb_124_112"};
const char*    fSysCombNames2        [fNumSysComb] = {"132 / 108",     "132 / 112",     "124 / 108",     "124 / 112"};
const char*    fSysCombTitles        [fNumSysComb] = {"(132 / 108)",  "(132 / 112)",   "(124 / 108)",   "(124 / 112)"};

const int      kf7 = 0, kx0 = 1, kf6 = 2, kx6 = 3;
const int      fNumAna = 4;
const int      fAnaIdx     [fNumAna] = {kf7, kx0, kf6, kx6};
const char*    fAnaFNames  [fNumAna] = {"f7", "x0", "fix6","x0_f6set"};
const char*    fAnaNames   [fNumAna] = {"f7", "x0", "f","x"};
const char*    fAnaONames  [fNumAna] = {"f7", "x0", "after","before"};
const char*    fAnaShort   [fNumAna] = {"f7", "x0", "f6","x6"};
const char*    fAnaTitles  [fNumAna] = {"After Fix", "Before Fix", "After Fix", "Before Fix"};
const char*    fAnaVersion [fNumAna] = {"NewAna.2107.4fd2bca", "NewAna.2107.4fd2bca","NewAna.2107.4fd2bca","NewAna.2107.4fd2bca"};

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
const int      kptoa0 = 5, kptoa50 = 6, kptoa100 = 7, kptoa150 = 8, kptoa200 = 9, kptoa250 = 10, kptoa300 = 11, kptoa350 = 12;
const int      kyF=13, kyH = 14;
const int      fNumCutY0s = 7;
const int      fNumCutPtoas = 8;
const int      fNumCutYPs = fNumCutY0s + fNumCutPtoas;
const int      fCutYPIdx[] = {kya, ky02, ky0, ky04, ky0610,kyF, kyH};
const int      fCutPtoaIdx[] = {kptoa0, kptoa50, kptoa100, kptoa150, kptoa200, kptoa250, kptoa300, kptoa350};
const int      fCutY0Idx[] = {kya, ky02, ky0, ky04, ky0610, kptoa0, kptoa50, kptoa100, kptoa150, kptoa200, kptoa250, kptoa300, kptoa350, kyF, kyH};
const char*    fCutYPNames[] = {"yAll","yGT02","yGT0","y0004","y0610", "ptoa0", "ptoa50", "ptoa100", "ptoa150", "ptoa200", "ptoa250", "ptoa300", "ptoa350", "yF", "yH"};
const char*    fCutYPTitles[] = {"", "y_{0}>0.2", "y_{0}>0", "0<y_{0}<.4", ".6<y_{0}<1.",
  "p_{T}/A=0~50", "p_{T}/A=50~100", "p_{T}/A=100~150", "p_{T}/A=150~200", "p_{T}/A=200~250", "p_{T}/A=250~300", "p_{T}/A=300~350", "p_{T}/A=350~400",
  //"-.025<y_{0}<1",};
  "-.025<y_{0}<0.5",
  ".7<y_{0}<1.5"};
const TCut     fCutYPValues[] = {
  "",
  "fy_cm/(by_cm/2)>0.2",
  "fy_cm/(by_cm/2)>0",
  "fy_cm/(by_cm/2)>0&&fy_cm/(by_cm/2)<0.4",
  "fy_cm/(by_cm/2)>0.6&&fy_cm/(by_cm/2)<1",

  "pt_cm/PARTICLEA>=0&&pt_cm/PARTICLEA<50",
  "pt_cm/PARTICLEA>=50&&pt_cm/PARTICLEA<100",
  "pt_cm/PARTICLEA>=100&&pt_cm/PARTICLEA<150",
  "pt_cm/PARTICLEA>=150&&pt_cm/PARTICLEA<200",
  "pt_cm/PARTICLEA>=200&&pt_cm/PARTICLEA<250",
  "pt_cm/PARTICLEA>=250&&pt_cm/PARTICLEA<300",
  "pt_cm/PARTICLEA>=300&&pt_cm/PARTICLEA<350",
  "pt_cm/PARTICLEA>=350&&pt_cm/PARTICLEA<400",

  //"fy_cm/(by_cm/2)>-.25&&fy_cm/(by_cm/2)<1",
  "fy_cm/(by_cm/2)>-.25&&fy_cm/(by_cm/2)<0.5",
  "fy_cm/(by_cm/2)>.7&&fy_cm/(by_cm/2)<1.5",
};
