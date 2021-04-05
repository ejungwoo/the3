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
const double   fParticleMass     [fNumParticles] = {938.272, 1871.06, 2809.41, 2809.41, 3728.4};
const int      fParticlePDGs     [fNumParticles] = {2212, 1000010020, 1000010030, 1000020030, 1000020040};
const double   fParticleSDHL     [fNumParticles] = {2.2, 2.0, 1.8, 1.8, 1.8};
const double   fParticleSDLL     [fNumParticles] = {2.2, 2.0, 1.8, 1.8, 1.8};

//const TCut     fParticlePozCut   [fNumParticles] = {"p_lab>100", "p_lab>200", "p_lab>200", "p_lab>400", "p_lab>400"};
const TCut     fParticlePozCut   [fNumParticles] = {"p_lab>100", "p_lab>200", "p_lab>200", "p_lab>350", "p_lab>400"};

const int      kSDAll = 0, kSD_0_x = 1, kSD_x_0 = 2, kSD_xx_l3 = 3, kSD_xx = 4;
const char*    fSDNames[]  =  {"sdall", "sd_0_SDVALUE ", "sd_SDVALUE_0", "sd_SDVALUE_1_l3", "asdSDVALUE"};
const char*    fSDTitles[] =  {"sd_{ANATTL2}All", "0<sd_{ANATTL2}<SDVALUE ", "SDVALUE<sd_{ANATTL2}<0", "|sd_{ANATTL2}|<SDVALUE_L3", "|sd_{ANATTL2}|<SDVALUE"};
const TCut     fSDParticleCut[][fNumParticles] = {
  {"","","","",""},
  {"sd<SDVALUE&&sd>0", "sd<SDVALUE&&sd>0", "sd<SDVALUE&&sd>0", "sd<SDVALUE&&sd>0", "sd<SDVALUE&&sd>0"},
  {"sd<0&&sd>-SDVALUE", "sd<0&&sd>-SDVALUE", "sd<0&&sd>-SDVALUE", "sd<0&&sd>-SDVALUE", "sd<0&&sd>-SDVALUE"},
  {"sd<SDVALUE&&sd>-1.5","abs(sd)<SDVALUE",  "abs(sd)<SDVALUE",  "abs(sd)<SDVALUE",  "sd<1.5&&sd>-SDVALUE"},
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
const char*    fSysTitles         []            = {"s(132+124)", "s(108+112)", "s(112+124)", "s(124+112)", "all-systems"};
const char*    fSysEnergyLossFile [fNumSyss] = {"data/Sn132Sn124.txt","data/Sn108Sn112.txt","data/Sn112Sn124.txt","data/Sn112Sn124.txt"};
int            fSysNumEvents      [fNumSyss] = {0};
const int      fSysColor          [fNumSyss] = {kRed, kBlue, kSpring-6, kOrange-3};
const int      fSysMStyle         [fNumSyss] = {25,26,27,28};
const double   fSysMSize          [fNumSyss] = {1.3,1.3,1.3,1.3};

const int      fSysNEventsAll     []         = {1849598, 1212050, 716286, 266066};
const int      fSysNEvents55      []         = {355036, 189764, 100080, 56008};

//const TCut     fSysCut            [fNumSyss] = {"n/z>0.81978&&n/z<0.860337", "n/z>0.711004&&n/z<0.749269", "n/z>0.751898&&n/z<0.790893", "n/z>0.76741&&n/z<0.806511"};
//const TCut     fSysCut            [fNumSyss] = {"n/z>0.799501&&n/z<0.880615", "n/z>0.691872&&n/z<0.768402", "n/z>0.732401&&n/z<0.81039", "n/z>0.74786&&n/z<0.826061"};
//const TCut     fSysCut            [fNumSyss] = { "n/z>0.76175&&n/z<0.909377", "n/z>0.657071&&n/z<0.801513", "n/z>0.697387&&n/z<0.845191", "n/z>0.710847&&n/z<0.858042"};
const TCut     fSysCut            [fNumSyss] = { "","","","" };

const int      fNumSysComb = 4;
const int      fSysCombIdxTest       [] = {0};
const int      fSysCombIdxSameTarget [] = {1,2};
const int      fSysCombIndx          [fNumSysComb]    = {0,1,2,3};
const int      fSysCombIdx           [fNumSysComb][2] = {{0,1},{0,2},{3,1},{3,2}};
const int      fSysCombBeam          [fNumSysComb][2] = {{132,108},{132,112},{124,108},{124,112}};
const char*    fSysCombNames         [fNumSysComb]    = {"comb_132_108", "comb_132_112", "comb_124_108",   "comb_124_112"};
const char*    fSysCombNames2        [fNumSysComb]    = {"132 / 108",     "132 / 112",     "124 / 108",     "124 / 112"};
const char*    fSysCombTitles        [fNumSysComb]    = {"(132 / 108)",  "(132 / 112)",   "(124 / 108)",   "(124 / 112)"};

const int      kF132 = 0, kB132 = 1, kFSys=2, kFNN=3, kFNN2=4, kFNN50=5;
const int      fNumAna = 6;
const int      fAnaIdx     [fNumAna] = {kF132, kB132, kFSys, kFNN, kFNN2, kFNN50};
const char*    fAnaFNames  [fNumAna] = {"f7", "x0", "fixPID", "nn", "nn2", "nn50"};
const char*    fAnaNames   [fNumAna] = {"f132", "x0", "fsys", "nn", "nn2", "nn50"};
const char*    fAnaONames  [fNumAna] = {"f132", "x0", "fsys", "nn", "nn2", "nn50"};
const char*    fAnaShort   [fNumAna] = {"f132", "x0", "fsys", "nn", "nn2", "nn50"};
const char*    fAnaTitles2 [fNumAna] = {"132", "b132", "sys", "nn", "nn2", "nn50"};
const char*    fAnaTitles  [fNumAna] = {"After_Fix(132)", "Before_Fix(132)", "After_Fix(Sys)", "NN C.M.", "NN C.M. (vadpoca cut)", "central cut 50"};
const char*    fAnaVersion [fNumAna] = {"NewAna.2107.4fd2bca", "NewAna.2107.4fd2bca", "NewAna.2107.4fd2bca","NewAna.2107.4fd2bca","NewAna.2107.4fd2bca","NewAna.2107.4fd2bca"};

const int      kLR = 0, kLeft = 1, kRight = 2;
const int      fNumLRs = 3;
const int      fLRIdx    [fNumLRs] = {0, 1, 2};
const char*    fLRNames  [fNumLRs] = {"all","left","right"};
const char*    fLRFNames [fNumLRs] = {"*","left","right"};
const char*    fLRTitles [fNumLRs] = {"TPC-All","TPC-Left","TPC-Right"};

const int      kMultAll = 0, kMult45 = 1, kMult55 = 2, kMult50=3;
const int      fNumMultOption = 4;
const int      fMultIdx     []               = {kMultAll, kMult45, kMult55, kMult50};
const char*    fMultFNames  [fNumMultOption] = {"*","45_54", "55_100", "50_100"};
const char*    fMultNames   [fNumMultOption] = {"45100","4554", "55100", "50100"};
const char*    fMultShort   [fNumMultOption] = {"ml", "mh", "mlh", "m50"};
const char*    fMultNames2  [fNumMultOption] = {"m45to100", "m45to54", "m55to100", "m50to100"};
const char*    fMultTitles  [fNumMultOption] = {"mult=45~100", "mult=45~54", "mult=55~100", "mult=50~100"};
const int      fMultLL      [fNumMultOption] = {45,  45,  55, 50};
const int      fMultHL      [fNumMultOption] = {100, 54, 100, 55};

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

const int      kYPAll = 0, kY1 = 1, kY2 = 2, kY3 = 3, kY4 = 4;
const int      kPtoa0 = 5, kPtoa50 = 6, kPtoa100 = 7, kPtoa150 = 8, kPtoa200 = 9, kPtoa250 = 10, kPtoa300 = 11, kPtoa350 = 12;
const int      kYF = 13, kPF3 = 14, kYPFI = 16, kPF4 = 15, kYPUD = 17;
const int      kP0 = 18, kP100 = 19, kP200 = 20, kP300 = 21;
const int      fNumCutYPs = 22;
const int      fCutY0Idx   [] = {kYPAll, kY1, kY2, kY3, kY4, kYF};
const int      fCutPtoaIdx [] = {kPtoa0, kPtoa50, kPtoa100, kPtoa150, kPtoa200, kPtoa250, kPtoa300, kPtoa350};
const int      fCutYPIdx   [] = {kYPAll, kY1, kY2, kY3, kY4, kPtoa0, kPtoa50, kPtoa100, kPtoa150, kPtoa200, kPtoa250, kPtoa300, kPtoa350, kYF, kPF3, kPF4, kYPFI, kYPUD, kP0, kP100, kP200, kP300};
const int      fCutPtoaIdx2[] = {kP0, kP100, kP200, kP300};
const char*    fCutYPNames [] = {"yAll","y1", "y2", "y3", "y4",
                                 "ptoa0", "ptoa50", "ptoa100", "ptoa150", "ptoa200", "ptoa250", "ptoa300", "ptoa350",
                                 "yF3", "pF3", "pF4", "ypFI", "ypUD",
                                 "p0", "p100", "p200", "p300"};
const char*    fCutYPTitles[] = {"ypAll", "y1", "y_{0}=-.25~.25","y_{0}=.25~.75","y_{0}=.75~1.25",
                                 "p_{T}/A=0~50", "p_{T}/A=50~100", "p_{T}/A=100~150", "p_{T}/A=150~200", "p_{T}/A=200~250", "p_{T}/A=250~300", "p_{T}/A=300~350", "p_{T}/A=350~400",
                                 "y_{0}=-0.25~1", "p_{T}/A=0~300", "p_{T}/A=0~400", "yp=UserDefinedFull", "yp=UserDefined",
                                 "p_{T}/A=0~100", "p_{T}/A=100~200", "p_{T}/A=200~300", "p_{T}/A=300~400"};
const TCut     fCutYPValues[] = {
  "",
  "fy_cm/(by_cm/2)>0",
  "fy_cm/(by_cm/2)>=-.25&&fy_cm/(by_cm/2)<0.25",
  "fy_cm/(by_cm/2)>=0.25&&fy_cm/(by_cm/2)<0.75",
  "fy_cm/(by_cm/2)>=0.75&&fy_cm/(by_cm/2)<1.25",

  "pt_cm/PARTICLEA>=0&&pt_cm/PARTICLEA<50",
  "pt_cm/PARTICLEA>=50&&pt_cm/PARTICLEA<100",
  "pt_cm/PARTICLEA>=100&&pt_cm/PARTICLEA<150",
  "pt_cm/PARTICLEA>=150&&pt_cm/PARTICLEA<200",
  "pt_cm/PARTICLEA>=200&&pt_cm/PARTICLEA<250",
  "pt_cm/PARTICLEA>=250&&pt_cm/PARTICLEA<300",
  "pt_cm/PARTICLEA>=300&&pt_cm/PARTICLEA<350",
  "pt_cm/PARTICLEA>=350&&pt_cm/PARTICLEA<400",

  "fy_cm/(by_cm/2)>-.25&&fy_cm/(by_cm/2)<1",

  "pt_cm/PARTICLEA>=0&&pt_cm/PARTICLEA<300",
  "pt_cm/PARTICLEA>=0&&pt_cm/PARTICLEA<400",

  // ypfi
  "(sqrt((p_cm)*(p_cm)+PARTICLEM*PARTICLEM)-PARTICLEM)/PARTICLEA>30&&(sqrt((p_cm)*(p_cm)+PARTICLEM*PARTICLEM)-PARTICLEM)/PARTICLEA<60",
  //"fy_cm/(by_cm/2)>-.15&&fy_cm/(by_cm/2)<.15&&pt_cm/PARTICLEA>=0&&pt_cm/PARTICLEA<100",
  //"fy_cm/(by_cm/2)>-0.25&&fy_cm/(by_cm/2)<0.25&&pt_cm/PARTICLEA>=0&&pt_cm/PARTICLEA<100",

  "",

  "pt_cm/PARTICLEA>=0&&pt_cm/PARTICLEA<100",
  "pt_cm/PARTICLEA>=100&&pt_cm/PARTICLEA<200",
  "pt_cm/PARTICLEA>=200&&pt_cm/PARTICLEA<300",
  "pt_cm/PARTICLEA>=300&&pt_cm/PARTICLEA<400",
};

const int      fDrawMStyle   [] = {24,25,26,27,28,30,42,46,3};
const double   fDrawMSize    [] = {1.7,1.7,1.7,2.0,2.0, 1.3,1.3,1.3,1.3};
const int      fDrawColor    [] = {kBlack, kRed, kBlue, kSpring-6, kOrange-3,kViolet-5, kGray+2, kPink+6, kAzure-1};

const int      fDrawMStyle2  [] = {24,25,26,27,28,30, 26,25,24};
const double   fDrawMSize2   [] = {1.5,1.5,1.5,2.0,1.8,1.9,1.3,1.3,1.3};
const int      fDrawColor2   [] = {kGray+1, kRed-8, kBlue-8, kGreen-8, kOrange-8, kViolet-8, kGray+1, kPink+2};

const int      fDrawMStyle3  [] = {20,21,22,33,3,28,26,25,24};
const double   fDrawMSize3   [] = {1.7,1.7,1.7,2.0,2.0, 1.3,1.3,1.3,1.3};
const int      fDrawColor3   [] = {kGray, kRed-10, kBlue-10, kGreen-10, kOrange-9, kViolet-9, kGray, kPink+1};
