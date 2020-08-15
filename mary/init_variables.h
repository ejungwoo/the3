#ifndef INIT_VARIABLES_H
#define INIT_VARIABLES_H

using ejungwoo::variable;
using ejungwoo::binning;
using ejungwoo::titles;
using ejungwoo::setpar;

const int fIPIDAll[] = {0,1,2,3,4};
const int fNumProtons[6]  = {1,1,1,2,2,2};
const int fNumNeutrons[6] = {0,1,2,1,2,4};

TString fttl_main = "$$(pname)_$$(sys)_$$(vshort)";
TString fttly_ntne = "N_{tracks} / N_{events}";
TString fttly_n1 = "N_{prob > 0.5}";
TString fttly_n2 = "N_{prob > 0.2, eff > 0.05, abs(sd) < 5}";

auto fvar_prob     = variable("prob" , "prob" , "$$(cut0)", titles(fttl_main, "probability"                 , fttly_ntne), binning(200,0,1.01));
auto fvar_eff      = variable("eff"  , "eff"  , "$$(cut0)", titles(fttl_main, "efficiency"                  , fttly_ntne), binning(200,0,1.01));
auto fvar_sd       = variable("sd"   , "sd"   , "$$(cut0)", titles(fttl_main, "distance_{PID mean} (#sigma)", fttly_ntne), binning(200,-5,5));
auto fvar_pt_cm    = variable("pt_cm", "pt_cm", "$$(cut0)", titles(fttl_main, "p_{T,CM} (MeV/c^{2})"        , fttly_ntne), binning(400,0,2000));
auto fvar_fy_cm    = variable("fy_cm", "fy_cm", "$$(cut0)", titles(fttl_main, "y_{Frag,CM}"                 , fttly_ntne), binning(100,-2.,2.));
auto fvar_p_cm     = variable("p_cm" , "p_cm" , "$$(cut0)", titles(fttl_main, "p_{CM} (MeV/c^{2})"          , fttly_ntne), binning(400,0,3000));
auto fvar_ke_cm    = variable("ke_cm", "ke_cm", "$$(cut0)", titles(fttl_main, "K.E._{CM} (MeV)"             , fttly_ntne), binning(40,0,400));
auto fvar_p_lab    = variable("p_lab", "p_lab", "$$(cut0)", titles(fttl_main, "p_{Lab} (MeV/c^{2})"         , fttly_ntne), binning(400,0,2000));
auto fvar_dedx     = variable("dedx" , "dedx" , "$$(cut0)", titles(fttl_main, "dE/dx"                       , fttly_ntne), binning(400,0,1500));
auto fvar_dpoca    = variable("dpoca", "dpoca", "$$(cut0)", titles(fttl_main, "distance_{POCA-PV} (mm)"     , fttly_ntne), binning(200,0,10));
auto fvar_nr       = variable("nr"   , "nc"   , "$$(cut0)", titles(fttl_main, "N_{row-cluster}"             , fttly_ntne), binning(150,1,150));
auto fvar_nl       = variable("nl"   , "nl"   , "$$(cut0)", titles(fttl_main, "N_{layer-cluster}"           , fttly_ntne), binning(150,1,150));
auto fvar_tta_cm   = variable("tta_cm" ,  "theta_cm*57.295780", "$$(cut0)", titles(fttl_main, "#theta_{CM}" , fttly_ntne), binning(200,   0,180));
auto fvar_phi_cm   = variable("phi_cm" ,    "phi_cm*57.295780", "$$(cut0)", titles(fttl_main, "#phi_{CM}"   , fttly_ntne), binning(200,-180,180));
auto fvar_tta_lab  = variable("tta_lab", "theta_lab*57.295780", "$$(cut0)", titles(fttl_main, "#theta_{lab}", fttly_ntne), binning(200,   0,180));
auto fvar_phi_lab  = variable("phi_lab",   "phi_lab*57.295780", "$$(cut0)", titles(fttl_main, "#phi_{lab}"  , fttly_ntne), binning(200,-180,180));

auto fvar_asd      = variable("asd"    , "abs(sd)"      , "$$(cut0)", titles(fttl_main, "|distance_{PID mean}} (#sigma)", fttly_ntne), binning(200,0,5));
auto fvar_poz_lab  = variable("poz_lab", "p_lab/$$(z)"  , "$$(cut0)" ,titles(fttl_main, "p_{Lab}/Z (MeV/c^{2})", fttly_ntne), binning(400,0,2000));
auto fvar_keoa_cm  = variable("keoa_cm", "ke_cm/$$(a)"  , "$$(cut0)", titles(fttl_main, "KE_{CM}/A (MeV)"      , fttly_ntne), binning(40,0,400));
auto fvar_ptoa_cm  = variable("ptoa_cm", "pt_cm/$$(a)"  , "$$(cut0)", titles(fttl_main, "p_{T}/A (MeV/c)"      , fttly_ntne), binning(100,0,1000));
auto fvar_ny_cm    = variable("ny_cm"  , "fy_cm/$$(ynn)", "$$(cut0)", titles(fttl_main, "y_{CM}/y_{NN}"        , fttly_ntne), binning(100,-2.,2.));

auto fvar_n      = variable("n"        , "np+nd+nt+nhe3+nhe4", "$$(cutn)", titles("", "all_$$(sys)", fttly_n1), binning(100,0,100));
auto fvar_np     = variable("np"       , "np"                , "$$(cutn)", titles("",   "p_$$(sys)", fttly_n1), binning(100,0,100));
auto fvar_nd     = variable("nd"       , "nd"                , "$$(cutn)", titles("",   "d_$$(sys)", fttly_n1), binning(100,0,100));
auto fvar_nt     = variable("nt"       , "nt"                , "$$(cutn)", titles("",   "t_$$(sys)", fttly_n1), binning(100,0,100));
auto fvar_nhe3   = variable("nhe3"     , "nhe3"              , "$$(cutn)", titles("", "he3_$$(sys)", fttly_n1), binning(100,0,100));
auto fvar_nhe4   = variable("nhe4"     , "nhe4"              , "$$(cutn)", titles("", "he4_$$(sys)", fttly_n1), binning(100,0,100));
auto fvar_ng     = variable("n_good"   , "n_good"            , "$$(cutn)", titles("", "all_$$(sys)", fttly_n2), binning(100,0,100));
auto fvar_ngp    = variable("np_good"  , "np_good"           , "$$(cutn)", titles("",   "p_$$(sys)", fttly_n2), binning(100,0,100));
auto fvar_ngd    = variable("nd_good"  , "nd_good"           , "$$(cutn)", titles("",   "d_$$(sys)", fttly_n2), binning(100,0,100));
auto fvar_ngt    = variable("nt_good"  , "nt_good"           , "$$(cutn)", titles("",   "t_$$(sys)", fttly_n2), binning(100,0,100));
auto fvar_nghe3  = variable("nhe3_good", "nhe3_good"         , "$$(cutn)", titles("", "he3_$$(sys)", fttly_n2), binning(100,0,100));
auto fvar_nghe4  = variable("nhe4_good", "nhe4_good"         , "$$(cutn)", titles("", "he4_$$(sys)", fttly_n2), binning(100,0,100));

TString fSTVersion;
TString fVShort;
TString fVOutShort;
TString fVOut;
double fPhiLL;
double fPhiHL;
double fTtaLL;
double fTtaHL;
double fSolidAngleInPi;
double fSolidAngle;
int fTrackMultHL;
double fProbLL;
double fEffLL;
double fPtoaLL;
double fPtoaHL;

void setversion(int iversion)
{
  fSTVersion = "NewAna.2070.0231288";
  fVShort = TString((ejungwoo::tok(fSTVersion,".",0))[0])+ejungwoo::tok(fSTVersion,".",1);
  fTrackMultHL = 50;
  fProbLL = 0.2;
  fEffLL = 0.08;
  fPtoaLL = 0;
  fPtoaHL = 10000;
  fPhiLL = 160;
  fPhiHL = 200;
  fTtaLL = 0;
  fTtaHL = 180;
  (&fvar_eff)->binn.setmax(0.3);
  setpar("cut0", "$$(cut_pe)");
  setpar("cutn", "");

  if (iversion==-1) {
    fVOutShort = "rawp";
    fProbLL = 0.1;
    fEffLL = 0.01;
    setpar("cut0", "$$(cut_p)");
  }
  else if (iversion==0) {
    fVShort = TString((ejungwoo::tok(fSTVersion,".",0))[0])+ejungwoo::tok(fSTVersion,".",1);
    fVOutShort = "tta711";
    fTtaLL = 70;
    fTtaHL = 110;
    (&fvar_eff)->binn.setmax(0.15);
    (&fvar_fy_cm)->binn.setminmax(-.5,.5);
    (&fvar_ny_cm)->binn.setminmax(-.5,.5);
  }
  else if (iversion==1) { fVOutShort = "nocut"; }
  else if (iversion==2) { fVOutShort = "corrp_xe"; setpar("cut0", "$$(cut_p)"); }

  fVOut = fVShort + "." + fVOutShort;
  fSolidAngleInPi = 2 * (fPhiHL-fPhiLL)/360 * (new TF1("satta","sin(x)",0,TMath::Pi())) -> Integral(fTtaLL*TMath::DegToRad(),fTtaHL*TMath::DegToRad());
  fSolidAngle = fSolidAngleInPi * TMath::Pi();

  vector<TString> entry_array;
  entry_array.push_back(Form("version=%s",fSTVersion.Data()));
  entry_array.push_back(Form("N_{track}>%d",fTrackMultHL));
  entry_array.push_back(Form("probability>%.2f",fProbLL));
  entry_array.push_back(Form("efficiency>%.2f",fEffLL));
  if (fPtoaLL<=0&&fPtoaHL>=5000) entry_array.push_back("p_{T} cut: X");
  else                         entry_array.push_back(Form("%.2f<#p{T}/A_{CM}<%.2f",fPtoaLL,fPtoaHL));
  entry_array.push_back(Form("%.f<#phi_{CM}<%.f",fPhiLL,fPhiHL));
  entry_array.push_back(Form("%.f<#theta_{CM}<%.f",fTtaLL,fTtaHL));
  entry_array.push_back(Form("#Delta#Omega=%.3f#pi",fSolidAngleInPi));
  setpar("angle_cut",Form("$$(tta_cm)>%f&&$$(tta_cm)<%f",fTtaLL,fTtaHL));
  setpar("poz_cut", "$$(poz_lab)>$$(pozll)");
  setpar("vshort", fVOut);
  setpar("pes_cut",Form("$$(eff)>%f",fEffLL));
  ejungwoo::gversionin(fSTVersion);
  ejungwoo::gversionout(fVOut);
  for (auto entry : entry_array) cout << entry << endl;
}

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
    if (v132 < num_tracks_cut || v108 < num_tracks_cut)
      continue;
    auto idxis = graph_is -> GetN();
    graph_is -> SetPoint(idxis, binc, v132/v108);
  }

  TObjArray array;
  for (auto graph : {graph_is, graph_v132, graph_v108}) {
    graph -> SetMarkerStyle(ejungwoo::markeri(idx_particle));
    graph -> SetLineColor(ejungwoo::colori(idx_particle));
    graph -> SetLineWidth(2);
    array.Add(graph);
  }

  return array;
}

TObjArray draw_ci(vector<TH1 *>hists, double ds, int isys)
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
      auto np = fNumProtons[ipid];
      auto nn = fNumNeutrons[ipid];
      vpsum += fNumProtons[ipid] * ntracks;
      vnsum += fNumNeutrons[ipid] * ntracks;
    }
    graph_np -> SetPoint(graph_np->GetN(), binn.value, vpsum/binn.w/ds);
    graph_nn -> SetPoint(graph_nn->GetN(), binn.value, vnsum/binn.w/ds);
    graph_rr -> SetPoint(graph_rr->GetN(), binn.value, vnsum/vpsum);
  }

  TObjArray array;
  int idx = 0;
  for (auto graph : {graph_nn, graph_np, graph_rr}) {
    graph -> SetMarkerStyle(ejungwoo::markeri(idx));
    graph -> SetMarkerColor(ejungwoo::colori(idx));
    graph -> SetLineColor(ejungwoo::colori(isys));
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
