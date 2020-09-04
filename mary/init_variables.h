#ifndef INIT_VARIABLES_H
#define INIT_VARIABLES_H

using ejungwoo::variable;
using ejungwoo::binning;
using ejungwoo::titles;
using ejungwoo::setpar;
using ejungwoo::dbstr;

const int fIPIDAll[] = {0,1,2,3,4};
const int fNumProtons[6]  = {1,1,1,2,2,2};
const int fNumNeutrons[6] = {0,1,2,1,2,4};
const int fSystems[] = {132   , 108   , 112   , 124   };
const int fTargetA[] = {124   , 112   , 124   , 112   };
const double fYAAs[] = {0.3822, 0.3647, 0.3538, 0.3902};
const double fYNNs[] = {0.3696, 0.3697, 0.3705, 0.3706};
int fMultLLCut[6] = {56,55,50,50};
//double fSDHLCut[6] = {5,5,5,5,5,5,};
double fSDHLCut[6] = {2.5,2,2,2,2,2,};
Long64_t fNumEvents[4] = {0};
TString fParticleTreeNames[6] = {"p","d","t","he3","he4","he6"};
TString fParticleNames[6] = {"p","d","t","he3","he4","he6"};
//const double fParticlePozCut[6] = {100,300,300,100,200,200};
const double fParticlePozCut[6] = {100,100,100,100,100,100};
double fParticleMass[6] = {938.272, 1871.06, 2809.41, 2809.41, 3728.4, 5606.55};
int fParticlePDGs[6] = {2212, 1000010030, 1000010020, 1000020040, 1000020030,};
int fParticleA[6] = {1,2,3,3,4,6};
int fParticleZ[6] = {1,1,1,2,2,2};

TString fttl_main = "$$(sys)$$(pname)";
TString fttly_ntne = "N_{tracks} / N_{events}";
TString fttly_n1 = "N_{prob > 0.5}";
TString fttly_n2 = "N_{prob > 0.2, eff > 0.05, abs(sd) < 5}";

auto fvar_prob     = variable("prob" , "prob" , "$$(cut0)", titles(fttl_main, "probability"                 , fttly_ntne), binning(200,0,1.01));
auto fvar_eff      = variable("eff"  , "eff"  , "$$(cut0)", titles(fttl_main, "efficiency"                  , fttly_ntne), binning(200,0,1.01));
auto fvar_sd       = variable("sd"   , "sd"   , "$$(cut0)", titles(fttl_main, "sd = distance from pid_mean (#sigma)", fttly_ntne), binning(200,-5,5));
auto fvar_pt_cm    = variable("pt_cm", "pt_cm", "$$(cut0)", titles(fttl_main, "p_{T,CM} (MeV/c^{2})"        , fttly_ntne), binning(400,0,2000));
auto fvar_by_cm    = variable("by_cm", "by_cm", "$$(cut0)", titles(fttl_main, "y_{NN}"                      , fttly_ntne), binning(40,-2.,2.));
auto fvar_fy_cm    = variable("fy_cm", "fy_cm", "$$(cut0)", titles(fttl_main, "y_{Frag,CM}"                 , fttly_ntne), binning(100,-2.,2.));
auto fvar_p_cm     = variable("p_cm" , "p_cm" , "$$(cut0)", titles(fttl_main, "p_{CM} (MeV/c^{2})"          , fttly_ntne), binning(400,0,3000));
auto fvar_ke_cm    = variable("ke_cm", "ke_cm", "$$(cut0)", titles(fttl_main, "K.E._{CM} (MeV)"             , fttly_ntne), binning(40,0,400));
//auto fvar_p_lab    = variable("p_lab", "p_lab", "$$(cut0)", titles(fttl_main, "p_{Lab} (MeV/c^{2})"         , fttly_ntne), binning(400,0,2000));
auto fvar_dedx     = variable("dedx" , "dedx" , "$$(cut0)", titles(fttl_main, "dE/dx"                       , fttly_ntne), binning(400,0,1500));
auto fvar_dpoca    = variable("dpoca", "dpoca", "$$(cut0)", titles(fttl_main, "distance_{POCA-PV} (mm)"     , fttly_ntne), binning(200,0,10));
auto fvar_nr       = variable("nr"   , "nc"   , "$$(cut0)", titles(fttl_main, "N_{row-cluster}"             , fttly_ntne), binning(150,1,150));
auto fvar_nl       = variable("nl"   , "nl"   , "$$(cut0)", titles(fttl_main, "N_{layer-cluster}"           , fttly_ntne), binning(150,1,150));
auto fvar_ntk      = variable("ntk"  , "nt"   , ""        , titles(fttl_main, "N_{tracks}"                  , fttly_ntne), binning(150,40,80));
auto fvar_tta_cm   = variable("tta_cm" ,  "theta_cm*57.295780", "$$(cut0)", titles(fttl_main, "#theta_{CM}" , fttly_ntne), binning(200,   0,180));
auto fvar_phi_cm   = variable("phi_cm" ,    "phi_cm*57.295780", "$$(cut0)", titles(fttl_main, "#phi_{CM}"   , fttly_ntne), binning(200,-180,180));
auto fvar_tta_lab  = variable("tta_lab", "theta_lab*57.295780", "$$(cut0)", titles(fttl_main, "#theta_{lab}", fttly_ntne), binning(200,   0,180));
auto fvar_phi_lab  = variable("phi_lab",   "phi_lab*57.295780", "$$(cut0)", titles(fttl_main, "#phi_{lab}"  , fttly_ntne), binning(200,-180,180));

auto fvar_foby_cm  = variable("foby_cm", "y0_$$(sys)$$(pname)", "fy_cm/(by_cm/2)", "$$(cut0)", titles(fttl_main, "y_{0} = y_{Frag.CM}/y_{NN}", fttly_ntne), binning(40,-2.,2.));
auto fvar_asd      = variable("asd"    , "abs(sd)"      , "$$(cut0)", titles(fttl_main, "|distance_{PID mean}} (#sigma)", fttly_ntne), binning(200,0,5));
auto fvar_poz_lab  = variable("poz_lab", "p_lab"        , "$$(cut0)" ,titles(fttl_main, "p_{Lab}/Z (MeV/c^{2})", fttly_ntne), binning(400,0,2000));
auto fvar_keoa_cm  = variable("keoa_cm", "ke_cm/$$(a)"  , "$$(cut0)", titles(fttl_main, "KE_{CM}/A (MeV)"      , fttly_ntne), binning(40,0,400));
auto fvar_ptoa_cm  = variable("ptoa_cm", "pt_cm/$$(a)"  , "$$(cut0)", titles(fttl_main, "p_{T}/A (MeV/c)"      , fttly_ntne), binning(100,0,1000));

auto fvar_n      = variable("nall"     , "n"         , "$$(cutn)", titles("", "all_$$(sys)", fttly_n1), binning(100,0,100));
auto fvar_np     = variable("np"       , "np"        , "$$(cutn)", titles("",   "p_$$(sys)", fttly_n1), binning(100,0,100));
auto fvar_nd     = variable("nd"       , "nd"        , "$$(cutn)", titles("",   "d_$$(sys)", fttly_n1), binning(100,0,100));
auto fvar_nt     = variable("nt"       , "nt"        , "$$(cutn)", titles("",   "t_$$(sys)", fttly_n1), binning(100,0,100));
auto fvar_nhe3   = variable("nhe3"     , "nhe3"      , "$$(cutn)", titles("", "he3_$$(sys)", fttly_n1), binning(100,0,100));
auto fvar_nhe4   = variable("nhe4"     , "nhe4"      , "$$(cutn)", titles("", "he4_$$(sys)", fttly_n1), binning(100,0,100));
auto fvar_ng     = variable("n_good"   , "n_good"    , "$$(cutn)", titles("", "all_$$(sys)", fttly_n2), binning(100,0,100));
auto fvar_ngp    = variable("np_good"  , "np_good"   , "$$(cutn)", titles("",   "p_$$(sys)", fttly_n2), binning(100,0,100));
auto fvar_ngd    = variable("nd_good"  , "nd_good"   , "$$(cutn)", titles("",   "d_$$(sys)", fttly_n2), binning(100,0,100));
auto fvar_ngt    = variable("nt_good"  , "nt_good"   , "$$(cutn)", titles("",   "t_$$(sys)", fttly_n2), binning(100,0,100));
auto fvar_nghe3  = variable("nhe3_good", "nhe3_good" , "$$(cutn)", titles("", "he3_$$(sys)", fttly_n2), binning(100,0,100));
auto fvar_nghe4  = variable("nhe4_good", "nhe4_good" , "$$(cutn)", titles("", "he4_$$(sys)", fttly_n2), binning(100,0,100));

TString fSTVersion;
TString fVOutShort;
TString fTrackMultHL;
TString fCut0String;
TString fSDHLString;
double fPhiLL;
double fPhiHL;
double fPhiLL2 = 0;
double fPhiHL2 = 0;
double fTtaLL;
double fTtaHL;
double fProbLL;
double fEffLL;
double fPtoaLL;
double fPtoaHL;
double fNyLL;
//
double fSolidAngleInPi;
double fSolidAngle;


vector<TString> setversion(TString version)
{
  TString confName = Form("data/conf_%s.C",version.Data());
  cout << confName << endl;
  gROOT -> Macro(confName);

  fSolidAngleInPi = 2 * ((fPhiHL-fPhiLL)+(fPhiHL2-fPhiLL2))/360 * (new TF1("satta","sin(x)",0,TMath::Pi())) -> Integral(fTtaLL*TMath::DegToRad(),fTtaHL*TMath::DegToRad());
  fSolidAngle = fSolidAngleInPi * TMath::Pi();

  vector<TString> entry_array;
  entry_array.push_back(Form("N_{track}>=%s",fTrackMultHL.Data()));
  if (fPhiLL2==0&&fPhiHL2==0) entry_array.push_back(Form("#phi_{cm}=%s~%s",dbstr(fPhiLL).Data(),dbstr(fPhiHL).Data()));
  else entry_array.push_back(Form("#phi_{cm}=%s~%s,%s~%s",dbstr(fPhiLL).Data(),dbstr(fPhiHL).Data(),dbstr(fPhiLL2).Data(),dbstr(fPhiHL2).Data()));
  entry_array.push_back(Form("#theta_{cm}=%s~%s",dbstr(fTtaLL).Data(),dbstr(fTtaHL).Data()));
  entry_array.push_back(Form("prob>%s",dbstr(fProbLL).Data()));
  entry_array.push_back(Form("eff>%s",dbstr(fEffLL).Data()));
  entry_array.push_back(Form("|sd|<%s",fSDHLString.Data()));
  cout << "===================================== conf" << endl;
  for (auto entry : entry_array) cout << entry << endl;
  cout << "==========================================" << endl;

  ejungwoo::gversionin(fSTVersion);
  ejungwoo::gversionout(fVOutShort);
  setpar("vshort", fVOutShort);
  setpar("mult_cut",Form("$$(ntk)>=%s",fTrackMultHL.Data()));
  setpar("mult_cut2",Form("$$(nall)>%s",fTrackMultHL.Data()));
  setpar("angle_cut",Form("$$(tta_cm)>%s&&$$(tta_cm)<%s",dbstr(fTtaLL).Data(),dbstr(fTtaHL).Data()));
  setpar("poz_cut", "$$(poz_lab)>$$(pozll)");
  setpar("pes_cut",Form("$$(prob)>%s&&$$(eff)>%s&&$$(asd)<%s",dbstr(fProbLL).Data(),dbstr(fEffLL).Data(),fSDHLString.Data()));
  setpar("rap_cut",Form("$$(foby_cm)>%s",dbstr(fNyLL).Data()));
  setpar("cut0",fCut0String);

  return entry_array;
}

void setpar_syspid(int isys, int ipid=0)
{
  setpar("sys",fSystems[isys]);
  setpar("yaa",fYAAs[isys]);
  setpar("ynn",fYNNs[isys]);
  setpar("num_events",fNumEvents[isys]);
  setpar("multll",fMultLLCut[isys]);
  setpar("pname",fParticleNames[ipid]);
  setpar("pozll",fParticlePozCut[ipid]);
  setpar("pdg",fParticlePDGs[ipid]);
  setpar("m",fParticleMass[ipid]);
  setpar("a",fParticleA[ipid]);
  setpar("z",fParticleZ[ipid]);
  setpar("nump",fNumProtons[ipid]);
  setpar("numn",fNumNeutrons[ipid]);
  setpar("asdcut",fSDHLCut[ipid]);
};

#endif
