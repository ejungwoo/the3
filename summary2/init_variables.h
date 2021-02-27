const bool     kDraw = true, kHold = false, kSet = true, kUnSet = false;
const int      kAll = -1;

const int      kP = 0, kD = 1, kT = 2, kHe3 = 3, kHe4 = 4;
const int      fNumParticles = 5;
const int      fParticleIdx      [fNumParticles] = {kP, kD, kT, kHe3, kHe4};
const int      fParticleNumPs    [fNumParticles] = {1, 1, 1, 2, 2};
const int      fParticleN        [fNumParticles] = {0, 1, 2, 1, 2};
const int      fParticleA        [fNumParticles] = {1, 2, 3, 3, 4};
const int      fParticleZ        [fNumParticles] = {1, 1, 1, 2, 2};
const char*    fParticleNames    []              = {"p", "d", "t", "he3", "he4", ""};
const char*    fParticleTitles   []              = {"p", "d", "t", "^{3}He", "^{4}He", ""};
const double   fParticlePozLLCut [fNumParticles] = {100, 100, 800, 400, 400};
const double   fParticleMass     [fNumParticles] = {938.272, 1871.06, 2809.41, 2809.41, 3728.4};
const int      fParticlePDGs     [fNumParticles] = {2212, 1000010020, 1000010030, 1000020030, 1000020040};
const double   fParticleSDHL     [fNumParticles] = {2.2, 2.0, 1.8, 1.8, 1.8};
const double   fParticleSDLL     [fNumParticles] = {2.2, 2.0, 1.8, 1.8, 1.8};

const TCut     fParticlePozCut   [fNumParticles] = {"p_lab>100", "p_lab>200", "p_lab>200", "p_lab>400", "p_lab>400"};

const int      kSDAll = 0, kSD_0_x = 1, kSD_x_0 = 2, kSD_xx_l3 = 3, kSD_xx = 4;
const char*    fSDNames[]  =  {"sdall", "sd_0_SDVALUE ", "sd_SDVALUE_0", "sd_SDVALUE_1_l3", "asdSDVALUE"};
const char*    fSDTitles[] =  {"sdall", "sd<SDVALUE,sd>0 ", "sd<0,sd>SDVALUE", "|sd|<SDVALUE_L3", "|sd|<SDVALUE"};
const TCut     fSDParticleCut[][fNumParticles] = {
  {"","","","",""},
  {"sd<SDVALUE&&sd>0", "sd<SDVALUE&&sd>0", "sd<SDVALUE&&sd>0", "sd<SDVALUE&&sd>0", "sd<SDVALUE&&sd>0"},
  {"sd<0&&sd>-SDVALUE", "sd<0&&sd>-SDVALUE", "sd<0&&sd>-SDVALUE", "sd<0&&sd>-SDVALUE", "sd<0&&sd>-SDVALUE"},
  {"sd<SDVALUE&&sd>-3","abs(sd)<SDVALUE",  "abs(sd)<SDVALUE",  "abs(sd)<SDVALUE",  "sd<3&&sd>-SDVALUE"},
  {"abs(sd)<SDVALUE",  "abs(sd)<SDVALUE",  "abs(sd)<SDVALUE",  "abs(sd)<SDVALUE",  "abs(sd)<SDVALUE"},
};

const int      kDOP = 0, kTOP = 1;
const int      fNumPR = 2;
const int      fPRIdx[fNumPR]            = {kDOP, kTOP};
const int      fPRParticleIdx[fNumPR][2] = {{1,0},{2,0}};
const char*    fPRNames[fNumPR]          = {"dop","top"};
const char*    fPRTitles[fNumPR]         = {"d/p","t/p"};

const int      k132 = 0, k108 = 1, k112 = 2, k124 = 3, kSysAll=4;
const int      fNumSyss = 4;
const int      fSysIdx            [fNumSyss] = {k132, k108, k112, k124};
const int      fSysBeams          [fNumSyss] = {132, 108, 112, 124};
const int      fSysTargets        [fNumSyss] = {124, 112, 124, 112};
const double   fSysYAAs           [fNumSyss] = {0.3822, 0.3647, 0.3538, 0.3902};
const double   fSysYNNs           [fNumSyss] = {0.3696, 0.3697, 0.3705, 0.3706};
const char*    fSysNames          []            = {"132", "108", "112", "124", "sysAll"};
const char*    fSysTitles         []            = {"system(132+124)", "system(108+112)", "system(124+112)", "system(112+124)", "all-systems"};
const char*    fSysEnergyLossFile [fNumSyss] = {"data/Sn132Sn124.txt","data/Sn108Sn112.txt","data/Sn112Sn124.txt","data/Sn112Sn124.txt"};
int            fSysNumEvents      [fNumSyss] = {0};
const int      fSysColor          [fNumSyss] = {kRed, kBlue, kSpring-6, kOrange-3};
const int      fSysMStyle         [fNumSyss] = {25,26,27,28};
const double   fSysMSize          [fNumSyss] = {1.3,1.3,1.3,1.3};

const int      fNumSysComb = 4;
const int      fSysCombIdxTest       [] = {0};
const int      fSysCombIdxSameTarget [] = {1,2};
const int      fSysCombIndx          [fNumSysComb]    = {0,1,2,3};
const int      fSysCombIdx           [fNumSysComb][2] = {{0,1},{0,2},{3,1},{2,3}};
const char*    fSysCombNames         [fNumSysComb]    = {"comb_132_108", "comb_132_112", "comb_124_108",   "comb_124_112"};
const char*    fSysCombNames2        [fNumSysComb]    = {"132 / 108",     "132 / 112",     "124 / 108",     "124 / 112"};
const char*    fSysCombTitles        [fNumSysComb]    = {"(132 / 108)",  "(132 / 112)",   "(124 / 108)",   "(124 / 112)"};

const int      kf7 = 0, kx0 = 1, kfp=2;
const int      fNumAna = 3;
const int      fAnaIdx     [fNumAna] = {kf7, kx0, kfp};
const char*    fAnaFNames  [fNumAna] = {"f7", "x0", "fixPID"};
const char*    fAnaNames   [fNumAna] = {"f7", "x0", "fp"};
const char*    fAnaONames  [fNumAna] = {"f7", "x0", "fp"};
const char*    fAnaShort   [fNumAna] = {"f7", "x0", "fp"};
const char*    fAnaTitles  [fNumAna] = {"After Fix", "Before Fix", "After Fix+1"};
const char*    fAnaVersion [fNumAna] = {"NewAna.2107.4fd2bca", "NewAna.2107.4fd2bca", "NewAna.2107.4fd2bca"};

const int      kLR = 0, kLeft = 1, kRight = 2;
const int      fNumLRs = 3;
const int      fLRIdx    [fNumLRs] = {0, 1, 2};
const char*    fLRNames  [fNumLRs] = {"all","left","right"};
const char*    fLRFNames [fNumLRs] = {"*","left","right"};
const char*    fLRTitles [fNumLRs] = {"TPC-All","TPC-Left","TPC-Right"};

const int      kMultAll = 0, kMult45 = 1, kMult55 = 2;
const int      fNumMultOption = 3;
const int      fMultIdx     []               = {kMultAll, kMult45, kMult55};
const char*    fMultFNames  [fNumMultOption] = {"*","45_54", "55_100"};
const char*    fMultNames   [fNumMultOption] = {"45100","4554", "55100"};
const char*    fMultShort   [fNumMultOption] = {"ml", "mh", "mlh"};
const char*    fMultNames2  [fNumMultOption] = {"m45to100", "m45to54", "m55to100"};
const char*    fMultTitles  [fNumMultOption] = {"mult=45~100", "mult=45~54", "mult=55~100"};
const int      fMultLL      [fNumMultOption] = {45,  45,  55};
const int      fMultHL      [fNumMultOption] = {100, 54, 100};

const int      kThetaAll = 0, kTheta0 = 1, kTheta20 = 2, kTheta40 = 3, kTheta60 = 4, kThetaLT60 = 5, kThetaGT60 = 6;
const int      fNumCutThetas = 7;
const int      fCutThetaIdx    [fNumCutThetas] = {kThetaAll, kTheta0, kTheta20, kTheta40, kTheta60, kThetaLT60, kThetaGT60};
const char*    fCutThetaNames  [fNumCutThetas] = {"ttaAll","ttaRG0_20","ttaRG20_40","ttaRG40_60","ttaRG60_80","ttaLT60","ttaGT60"};
const char*    fCutThetaTitles [fNumCutThetas] = {"All-#theta_{lab}","0<=#theta_{lab}<20","20<=#theta_{lab}<40","40<=#theta_{lab}<60","60<=#theta_{lab}<80","#theta_{lab}<60","#theta_{lab}>=60"};
const double   fCutThetaRanges [fNumCutThetas][2] = {{0,90},{0,20},{20,40},{40,60},{60,80},{0,60},{60,90}};
const TCut     fCutThetaValues [fNumCutThetas] = {
  "",
  "theta_lab>= 0*TMath::DegToRad()&&theta_lab<20*TMath::DegToRad()",
  "theta_lab>=20*TMath::DegToRad()&&theta_lab<40*TMath::DegToRad()",
  "theta_lab>=40*TMath::DegToRad()&&theta_lab<60*TMath::DegToRad()",
  "theta_lab>=60*TMath::DegToRad()&&theta_lab<80*TMath::DegToRad()",
  "theta_lab <60*TMath::DegToRad()",
  "theta_lab>=60*TMath::DegToRad()",
};

const int      kYPAll = 0, kYGT02 = 1, kYGT0 = 2, kYGT04 = 3, kY0610 = 4;
const int      kPtoa0 = 5, kPtoa50 = 6, kPtoa100 = 7, kPtoa150 = 8, kPtoa200 = 9, kPtoa250 = 10, kPtoa300 = 11, kPtoa350 = 12;
const int      kYF = 13, kPF = 14, kYPF = 15;
const int      fNumCutYPs = 15;
const int      fCutY0Idx   [] = {kYPAll, kYGT02, kYGT0, kYGT04, kY0610, kYF, kPF};
const int      fCutPtoaIdx [] = {kPtoa0, kPtoa50, kPtoa100, kPtoa150, kPtoa200, kPtoa250, kPtoa300, kPtoa350};
const int      fCutYPIdx   [] = {kYPAll, kYGT02, kYGT0, kYGT04, kY0610, kPtoa0, kPtoa50, kPtoa100, kPtoa150, kPtoa200, kPtoa250, kPtoa300, kPtoa350, kYF, kPF, kYPF};
const char*    fCutYPNames [] = {"yAll","yGT02","yGT0","y0004","y0610", "ptoa0", "ptoa50", "ptoa100", "ptoa150", "ptoa200", "ptoa250", "ptoa300", "ptoa350", "yF", "pF", "ypF"};
const char*    fCutYPTitles[] = {"ypAll", "y_{0}>0.2", "y_{0}>0", "0<y_{0}<.4", ".6<y_{0}<1.",
                                 "p_{T}/A=0~50", "p_{T}/A=50~100", "p_{T}/A=100~150", "p_{T}/A=150~200", "p_{T}/A=200~250", "p_{T}/A=250~300", "p_{T}/A=300~350", "p_{T}/A=350~400", "-.025<y_{0}<1", "0<p_{T}/A<300", "-0.25<y_{0}<1,0<p_{T}/A<300"};
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

  "fy_cm/(by_cm/2)>-.25&&fy_cm/(by_cm/2)<1",
  "pt_cm/PARTICLEA>=0&&pt_cm/PARTICLEA<400",
  "fy_cm/(by_cm/2)>-.25&&fy_cm/(by_cm/2)<1&&pt_cm/PARTICLEA>=0&&pt_cm/PARTICLEA<400",
};

const int      fDrawColor    [] = {kBlack, kRed, kBlue, kSpring-6, kOrange-3,kViolet-5,kAzure-1,kPink+7};
const int      fDrawMStyle   [] = {24,25,26,27,28,30,42,46};
const double   fDrawMSize    [] = {1.3,1.3,1.3,1.3,1.3,1.3,1.5,1.3};
