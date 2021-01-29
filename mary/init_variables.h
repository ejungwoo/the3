#ifndef INIT_VARIABLES_H
#define INIT_VARIABLES_H

using ejungwoo::variable;
using ejungwoo::binning;
using ejungwoo::titles;
using ejungwoo::setpar;
using ejungwoo::dbstr;

void init();
void settrees();
vector<TString> setversion(TString version);
void setpar_syspid(int isys, int ipid=0);

const int fIPIDAll[] = {0,1,2,3,4};
const int fNumProtons[6]  = {1,1,1,2,2,0};
const int fNumNeutrons[6] = {0,1,2,1,2,0};
const int fSystems[] = {132   , 108   , 112   , 124   };
const int fTargetA[] = {124   , 112   , 124   , 112   };
const double fYAAs[] = {0.3822, 0.3647, 0.3538, 0.3902};
const double fYNNs[] = {0.3696, 0.3697, 0.3705, 0.3706};
int fMultLLCut[6] = {56,55,55,55};
int fMultHLCut[6] = {100,100,100,100};
//double fSDHLCut[6] = {5,5,5,5,5,5,};
double fSDHLCut[6] = {2.2,2.0,1.8,1.8,1.8,1.8,};
//double fSDHLCut[6] = {2.0,2.0,2.0,2.0,2.0,2.0,};
//double fSDHLCut[6] = {3.0,3.0,3.0,3.0,3.0,3.0,};
//double fSDHLCut[6] = {5.0,5.0,5.0,5.0,5.0,5.0,};
Long64_t fNumEvents[4] = {0};
TString fParticleTreeNames[6] = {"p","d","t","he3","he4","he6"};
TString fParticleNames[6] = {"p","d","t","he3","he4","he6"};
TString fParticleNames2[6] = {"p","d","t","^{3}He","^{4}He","^{6}He"};
//const double fParticlePozCut[6] = {100,300,300,100,200,200};
//const double fParticlePozCut[6] = {100,500,500,500,500,500};
const double fParticlePozCut[6] = {100,100,800,400,400,100};
double fParticleMass[6] = {938.272, 1871.06, 2809.41, 2809.41, 3728.4, 5606.55};
int fParticlePDGs[6] = {2212, 1000010020, 1000010030, 1000020030, 1000020040,};
int fParticleA[6] = {1,2,3,3,4,6};
int fParticleZ[6] = {1,1,1,2,2,2};

TChain *fTreePID[4][5] = {0};
TChain *fTreeMult[4] = {0};
bool fAddEffk = false;

TString fttl_main = "$$(sys)$$(pname)";
TString fttly_ntne = "N_{tracks} / N_{events}";
TString fttly_n1 = "N_{prob > 0.5}";
TString fttly_n2 = "N_{prob > 0.2, eff > 0.05, abs(sd) < 5}";

auto fvar_prob     = variable("prob" , "prob" , "$$(cut0)", titles(fttl_main, "probability"                 , fttly_ntne), binning(200,0,1.01));
auto fvar_eff      = variable("eff"  , "eff"  , "$$(cut0)", titles(fttl_main, "efficiency"                  , fttly_ntne), binning(200,0,1.01));
auto fvar_sd       = variable("sd"   , "sd"   , "$$(cut0)", titles(fttl_main, "sd = distance from pid_mean (#sigma)", fttly_ntne), binning(200,-5,5));
auto fvar_pt_cm    = variable("pt_cm", "pt_cm", "$$(cut0)", titles(fttl_main, "p_{T,CM} (MeV/c^{2})"        , fttly_ntne), binning(400,0,2000));
auto fvar_by_cm    = variable("by_cm", "by_cm", "$$(cut0)", titles(fttl_main, "2*y_{Beam,NN}"               , fttly_ntne));
auto fvar_by_cm0   = variable("by_cm0", "by_cm0", "$$(cut0)", titles(fttl_main, "y_{Beam,AA?}"              , fttly_ntne));
auto fvar_fy_cm    = variable("fy_cm", "fy_cm", "$$(cut0)", titles(fttl_main, "y_{Frag,CM}"                 , fttly_ntne), binning(100,-2.,2.));
auto fvar_fy_lab    = variable("fy_lab", "fy_lab", "$$(cut0)", titles(fttl_main, "y_{Frag,Lab}"                 , fttly_ntne), binning(100,-2.,2.));
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

auto fvar_foby_cm  = variable("foby_cm", "y0_$$(sys)$$(pname)", "fy_cm/(by_cm/2)", "$$(cut0)", titles(fttl_main, "y_{0} = y_{Frag.CM}/y_{NN}", fttly_ntne), binning(100,-2.,2.));
auto fvar_foby_lab = variable("foby_lab", "y0L_$$(sys)$$(pname)", "fy_lab/(by_cm)", "$$(cut0)", titles(fttl_main, "y_{0,Lab} = y_{Frag,Lab}/y_{beam}", fttly_ntne), binning(100,0,2.));
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

auto fvar_pid = fvar_poz_lab+fvar_dedx+"z";
auto fvar_ypt_cm = fvar_foby_cm+fvar_ptoa_cm;
auto fvar_ypt_lab = fvar_foby_lab+fvar_ptoa_cm;

TString fSTVersion;
TString fVOutShort;
TString fTrackMultHL;
TString fCut0String;
TString fSDHLString;
TString fHeadName;
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

void init()
{
  setpar("cut_x",                               "1.*($$(mult_cut))*($$(pes_cut))*($$(poz_cut))*($$(angle_cut))*($$(rap_cut))");
  setpar("cut_p",          "$$(prob)/$$(num_events)*($$(mult_cut))*($$(pes_cut))*($$(poz_cut))*($$(angle_cut))*($$(rap_cut))");
  setpar("cut_e",        "1./$$(eff)/$$(num_events)*($$(mult_cut))*($$(pes_cut))*($$(poz_cut))*($$(angle_cut))*($$(rap_cut))");
  setpar("cut_pe", "$$(prob)/$$(eff)/$$(num_events)*($$(mult_cut))*($$(pes_cut))*($$(poz_cut))*($$(angle_cut))*($$(rap_cut))");
  setpar("cut_noy","$$(prob)/$$(eff)/$$(num_events)*($$(mult_cut))*($$(pes_cut))*($$(poz_cut))*($$(angle_cut))");
}

TString fFileVersion = "";
TString fFileVersion2 = "";
TString fFileVersion3 = "";

bool tree_is_set = false;

void settrees()
{
  if (tree_is_set) return;
  tree_is_set = true;

  for (auto isys : {0,1,2,3}) {
    setpar_syspid(isys);
    auto sys = fSystems[isys];

    if (fTreeMult[isys]==nullptr) fTreeMult[isys] = new TChain("mult");
    auto fileName = Form("data/sys%d_%s.ana.particle.root",sys,fFileVersion.Data());
    cout << " ++ " << fileName << endl;
    fTreeMult[isys] -> Add(fileName);
    if (!fFileVersion2.IsNull()) {
      auto fileName2 = Form("data/sys%d_%s.ana.particle.root",sys,fFileVersion2.Data());
      cout << " ++2 " << fileName2 << endl;
      fTreeMult[isys] -> Add(fileName2);
    }
    if (!fFileVersion3.IsNull()) {
      auto fileName2 = Form("data/sys%d_%s.ana.particle.root",sys,fFileVersion3.Data());
      cout << " ++3 " << fileName2 << endl;
      fTreeMult[isys] -> Add(fileName2);
    }
    fNumEvents[isys] = fTreeMult[isys] -> GetEntries(ejungwoo::getpar("mult_cut2"));
    cout << isys << " " << fNumEvents[isys] << endl;

    if (0) {
      for (int ipid : fIPIDAll)
      {
        if (fTreePID[isys][ipid]==nullptr)
          fTreePID[isys][ipid] = new TChain("data");

        auto fileName = Form("data/sys%d_prob05.root",fSystems[isys]);
        if (ipid==0) cout << " ++ " << fileName << endl;
        fTreePID[isys][ipid] -> Add(fileName);
      }
    }
    else {
      for (int ipid : fIPIDAll) {
        bool addeffk = ((fAddEffk&&(sys==132||sys==108)&&(ipid==0||ipid==1||ipid==2))?true:false);
        if (fTreePID[isys][ipid]==nullptr) fTreePID[isys][ipid] = new TChain(fParticleTreeNames[ipid]);
        auto fileName = Form("data/sys%d_%s.ana.particle.root",sys,fFileVersion.Data());
        if (ipid==0) cout << " ## " << fileName << endl;
        fTreePID[isys][ipid] -> Add(fileName);
        if (addeffk) {
          auto fileName1 = Form("data/sys%d_%s.ana.particle.effk_%s.root",sys,fFileVersion.Data(),fParticleNames[ipid].Data());
          cout << fileName1 << endl;
          auto file1 = new TFile(fileName1,"read");
          auto tree1 = (TTree *) file1 -> Get(fParticleTreeNames[ipid]);
          fTreePID[isys][ipid] -> AddFriend(tree1);
        }
        if (!fFileVersion2.IsNull()) {
          auto fileName2 = Form("data/sys%d_%s.ana.particle.root",sys,fFileVersion2.Data());
          if (ipid==0) cout << " ## " << fileName2 << endl;
          fTreePID[isys][ipid] -> Add(fileName2);
          //if (addeffk) fTreePID[isys][ipid] -> AddFriend(Form("data/sys%d.%s.ana.particle.effk_%s.root",sys,fFileVersion2.Data(),fParticleNames[ipid].Data()));
        }
        if (!fFileVersion3.IsNull()) {
          auto fileName2 = Form("data/sys%d_%s.ana.particle.root",sys,fFileVersion3.Data());
          if (ipid==0) cout << " ## " << fileName2 << endl;
          fTreePID[isys][ipid] -> Add(fileName2);
          //if (addeffk) fTreePID[isys][ipid] -> AddFriend(Form("data/sys%d.%s.ana.particle.effk_%s.root",sys,fFileVersion3.Data(),fParticleNames[ipid].Data()));
        }
      }
    }
  }
}


vector<TString> setversion(TString version)
{
  TString confName = Form("data/conf_%s.C",version.Data());
  cout << confName << endl;
  gROOT -> Macro(confName);

  if (fFileVersion.IsNull()) fFileVersion = fSTVersion;

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
  //setpar("mult_cut",Form("$$(ntk)>=%s",fTrackMultHL.Data()));
  setpar("mult_cut","1");
  //setpar("mult_cut2",Form("$$(nall)>%s",fTrackMultHL.Data()));
  setpar("mult_cut2","1");
  //setpar("angle_cut",Form("$$(tta_cm)>%s&&$$(tta_cm)<%s",dbstr(fTtaLL).Data(),dbstr(fTtaHL).Data()));
  setpar("angle_cut","1");
  setpar("poz_cut", "$$(poz_lab)>$$(pozll)");
  setpar("pes_cut",Form("$$(prob)>%s&&$$(eff)>%s&&$$(asd)<%s",dbstr(fProbLL).Data(),dbstr(fEffLL).Data(),fSDHLString.Data()));
  setpar("rap_cut",Form("$$(foby_cm)>%s",dbstr(fNyLL).Data()));
  setpar("cut0",fCut0String);

  settrees();

  return entry_array;
}

void setpar_syspid(int isys, int ipid=0)
{
  setpar("sys",fSystems[isys]);
  setpar("tar",fTargetA[isys]);
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
  setpar("gcutpid",Form("(prob>.5)*gcut%d",ipid));
};

TObjArray draw_is(int idx_particle, TH1 *hist132, TH1 *hist108, double num_tracks_cut)
{
  auto graph_is = new TGraphErrors();
  auto graph_v132 = new TGraphErrors();
  auto graph_v108 = new TGraphErrors();
  binning binn(hist132);

  for (auto bin=1; bin<=binn.n; ++bin)
  {
    auto binc = binn.getc(bin);
    auto v132 = hist132 -> GetBinContent(bin);
    auto v108 = hist108 -> GetBinContent(bin);
    auto idxvv = graph_v132 -> GetN();
    graph_v132 -> SetPoint(idxvv, binc, v132);
    graph_v108 -> SetPoint(idxvv, binc, v108);
    if (v132 < num_tracks_cut || v108 < num_tracks_cut) {
      cout << "X  " << idxvv << " " << binc << " 132:" << v132 << " / 108:" << v108 << "  =  " << v132/v108 << endl;
      continue;
    }
    cout << "O  " << idxvv << " " << binc << " 132:" << v132 << " / 108:" << v108 << "  =  " << v132/v108 << endl;
    auto idxis = graph_is -> GetN();
    graph_is -> SetPoint(idxis, binc, v132/v108);
  }

  TObjArray array;
  for (auto graph : {graph_is, graph_v132, graph_v108}) {
    graph -> SetMarkerStyle(ejungwoo::markeri(idx_particle));
    graph -> SetLineColor(ejungwoo::colori(idx_particle));
    //graph -> SetLineWidth(2);
    array.Add(graph);
  }

  return array;
}

TObjArray draw_ci(vector<TH1 *>hists, double ds, int isys, TH1 *histn0=nullptr)
{
  auto graph_np = new TGraphErrors();
  auto graph_nn = new TGraphErrors();
  auto graph_rr = new TGraphErrors();

  binning binn(hists[0]);
  binn.resetb();
  while (binn.nextb()) {
    binn.idx;
    double vpsum = 0;
    double vnsum = 0;
    for (auto ipid : fIPIDAll) {
      double ntracks = hists[ipid] -> GetBinContent(binn.idx);
      vpsum += fNumProtons[ipid] * ntracks;
      vnsum += fNumNeutrons[ipid] * ntracks;
    }
    if (histn0!=nullptr) {
      double ntracks = histn0 -> GetBinContent(binn.idx);
      vnsum += ntracks;
    }

    //graph_np -> SetPoint(graph_np->GetN(), binn.value, vpsum/binn.w/ds);
    //graph_nn -> SetPoint(graph_nn->GetN(), binn.value, vnsum/binn.w/ds);
    graph_np -> SetPoint(graph_np->GetN(), binn.value, vpsum);
    graph_nn -> SetPoint(graph_nn->GetN(), binn.value, vnsum);
    graph_rr -> SetPoint(graph_rr->GetN(), binn.value, vnsum/vpsum);
  }

  TObjArray array;
  int idx = 0;
  for (auto graph : {graph_nn, graph_np, graph_rr}) {
    graph -> SetMarkerStyle(ejungwoo::markeri(idx));
    graph -> SetMarkerColor(ejungwoo::colori(idx));
    //graph -> SetLineColor(ejungwoo::colori(isys));
    graph -> SetLineWidth(2);
    array.Add(graph);
    idx++;
  }

  return array;
}

TGraphErrors *draw_ratio(TGraphErrors *graph1, TGraphErrors *graph2, double max)
{
  auto graph_dr = new TGraphErrors();

  double x1, y1, x2, y2;
  auto nn = graph1 -> GetN();
  for (auto i=0; i<nn; ++i) {
    graph1 -> GetPoint(i,x1,y1);
    graph2 -> GetPoint(i,x2,y2);
    if (x1 > max)
      break;
    if (y2>0) {
      graph_dr -> SetPoint(graph_dr->GetN(),x1,y1/y2);
    }
  }

  return graph_dr;
}

#endif
