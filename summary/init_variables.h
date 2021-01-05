const int      fNumParticles = 5;
const int      fParticleNumPs    [fNumParticles] = {1, 1, 1, 2, 2};
const int      fParticleNumNs    [fNumParticles] = {0, 1, 2, 1, 2};
const int      fParticleIdx      [fNumParticles] = {0, 1, 2, 3, 4};
const int      fParticleA        [fNumParticles] = {1, 2, 3, 3, 4}; const int      fParticleZ        [fNumParticles] = {1, 1, 1, 2, 2};
const TString  fParticleNames    [fNumParticles] = {"p", "d", "t", "he3", "he4"};
const TString  fParticleNames2   [fNumParticles] = {"p", "d", "t", "^{3}He", "^{4}He"};
const double   fParticlePozLLCut [fNumParticles] = {100, 100, 800, 400, 400};
const double   fParticleMass     [fNumParticles] = {938.272, 1871.06, 2809.41, 2809.41, 3728.4};
const int      fParticlePDGs     [fNumParticles] = {2212, 1000010030, 1000010020, 1000020040};
      double   fParticleSDHL     [fNumParticles] = {2.2, 2.0, 1.8, 1.8, 1.8};
      double   fParticleSDLL     [fNumParticles] = {2.2, 2.0, 1.8, 1.8, 1.8};

const int      fNumSystems = 4;
const int      fSystems           [fNumSystems] = {132, 108, 112, 124};
const int      fSystemTargets     [fNumSystems] = {124, 112, 124, 112};
const double   fSystemYAAs        [fNumSystems] = {0.3822, 0.3647, 0.3538, 0.3902};
const double   fSystemYNNs        [fNumSystems] = {0.3696, 0.3697, 0.3705, 0.3706};
const TString  fSystemNames       [fNumSystems] = {"sys132", "sys108", "sys124", "sys112"};
const TString  fSystemNames2      [fNumSystems] = {"132+124", "108+112", "124+112", "112+124"};
const TString  fSysEnergyLossFile [fNumSystems] = {"data/Sn132Sn124.txt","data/Sn108Sn112.txt","data/Sn112Sn124.txt","data/Sn112Sn124.txt"};
      Long64_t fSystemNumEvents   [fNumSystems] = {0};

const int      fNumMultOption = 2;
const int      fMultOptions [fNumMultOption] = {0, 1};
const int      fMultLL      [fNumMultOption] = {45, 55};
const int      fMultHL      [fNumMultOption] = {54, 100};

const int      fNumSysComb = 4;
const int      fSysCombIdxSameTarget [] = {1,2};
const int      fSysCombIndx          [fNumSysComb] = {0,1,2,3};
const int      fSysCombIdx           [fNumSysComb][2] = {{0,1},{0,3},{2,1},{2,3}};
const TString  fSysCombNames         [fNumSysComb] = {"comb_132_108", "comb_132_112", "comb_124_108", "comb_124_112"};
const TString  fSysCombNames2        [fNumSysComb] = {"132 / 108", "132 / 112", "124 / 108", "124 / 112"};
const TString  fSysCombNames3        [fNumSysComb] = {"(132 / 108)", "(132 / 112)", "(124 / 108)", "(124 / 112)"};

auto hist_prob      = new ehist("prob", ";PID probability",             "prob", "",  200, 0, 1, "logy=1");
auto hist_eff       = new ehist("eff",  ";Track efficiency",            "eff",  "",  200, 0, 1, "logy=1");
auto hist_sd        = new ehist("sd",   ";PID sigma-distance (#sigma)", "sd",   "",  200, -5, 5);

auto hist_y0        = new ehist("y0",   ";Rapidity y_{0}",   "fy_cm/(by_cm/2)", "", 200, -1, 2);
auto hist_ptoa      = new ehist("ptoa", ";p_{T}/A (MeV/c)",  "pt_cm/[0]",       "", 200, 0, 1000);

auto hist_ke_cm     = new ehist("ke_cm",   ";KE_{CM} (MeV)",   "ke_cm",  "", 200, 0, 500);
auto hist_p_lab     = new ehist("p_lab",   ";p_{Lab} (MeV/c)", "p_lab",  "", 200, 0, 2500);
auto hist_dedx      = new ehist("dedx",    ";dE/dx",           "dedx",   "", 200, 0, 500);

auto hist_theta_cm  = new ehist("theta_cm",  ";#theta_{CM}",  "TMath::RadToDeg()*theta_cm",  "", 200, 0, 180);
auto hist_phi_cm    = new ehist("phi_cm",    ";#phi_{CM}",    "TMath::RadToDeg()*phi_cm",    "", 200, -180, 180);
auto hist_theta_lab = new ehist("theta_lab", ";#theta_{Lab}", "TMath::RadToDeg()*theta_lab", "", 200, 0, 180);
auto hist_phi_lab   = new ehist("phi_lab",   ";#phi_{Lab}",   "TMath::RadToDeg()*phi_lab",   "", 200, -180, 180);

auto hist_pid = new ehist(hist_p_lab, hist_dedx);

vector<ehist *> fHists1;
vector<ehist *> fHists2;

void Init()
{
  fHists1.push_back(hist_prob     );
  fHists1.push_back(hist_eff      );
  fHists1.push_back(hist_sd       );
  fHists1.push_back(hist_y0       );
  fHists1.push_back(hist_ptoa     );
  fHists1.push_back(hist_ke_cm    );
  fHists1.push_back(hist_p_lab    );
  fHists1.push_back(hist_dedx     );
  fHists1.push_back(hist_theta_cm );
  fHists1.push_back(hist_phi_cm   );
  fHists1.push_back(hist_theta_lab);
  fHists1.push_back(hist_phi_lab  );

  fHists2.push_back(hist_pid);
}
