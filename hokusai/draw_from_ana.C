#include "/home/ejungwoo/config/ejungwoo.h"

void draw_from_ana(
  int system=132,
  int iSplit=1,
  int version=1
  )
{
  ejungwoo::gstats(0);

  TString file_name_in  = TString("/home/ejungwoo/data/pid4/Sn")+system+"_"+iSplit+"_ana_v"+version+".NewAna.2034.45b9400.root";
  TString file_name_out = TString("/home/ejungwoo/the3/hokusai/data_xml/summary_hist_")+system+"_s"+iSplit+"_v"+version+".root";

  auto yBeamCM = 0.36;
  auto probCut = 0.2;
  auto effCut = 0.1;

  TString particle_name[] = {"p","d","t","he3","he4"};
  double particle_mass[] = {938.272 ,1871.06 ,2809.41 ,2809.41 ,3728.4};
  double particle_a[] = {1,2,3,3,4,};

  ejungwoo::binning binning_pt{200,0,1500};
  ejungwoo::binning binning_yyCM{200,-1.5,3.0};
  ejungwoo::binning binning_keCM{200,0,500};
  ejungwoo::binning binning_phiCM{200,-180,180};
  ejungwoo::binning binning_ttaCM{200,0,180};
  ejungwoo::binning binning_pLab{400,0,3000};
  ejungwoo::binning binning_dedx{400,0,1500};

  ejungwoo::titles ttl_yypt{"","y_{CM}/y_{beam,CM}","p_{T} (MeV/c)",""};
  ejungwoo::titles ttl_keCM{"","KE_{CM}/A (MeV)","",""};
  ejungwoo::titles ttl_phiCM{"","#phi (Deg.)","",""};
  ejungwoo::titles ttl_ttaCM{"","#theta (Deg.)","",""};
  ejungwoo::titles ttl_pLab{"","p/A (MeV/c)","",""};
  ejungwoo::titles ttl_dedx{"","dE/dx","",""};

  TString ttl_x = " ";
  TString ttl_p = " (prob. cor.) ";
  TString ttl_a = " (prob., eff. cor.) ";
  TString ttl_r = " (raw) ";
  TString ttl_t = " (prob., eff. cut) ";

  auto www_p = [](int idx) { return TCut(Form("Prob[%d].fElements",idx)); };
  auto www_a = [](int idx) { return TCut(Form("Prob[%d].fElements/Eff[%d].fElements",idx,idx)); };
  auto cut_a = [probCut,effCut](int idx) { return TCut(Form("(Prob[%d].fElements > %f && Eff[%d].fElements > %f)",idx,probCut,idx,effCut)); };

  auto tree = new TChain("cbmsim");
  tree -> Add(file_name_in);
  TFile *file_hist = new TFile(file_name_out,"recreate");

  for (auto idx=0; idx<5; ++idx)
  {
    TString pname = particle_name[idx];
    TString mttl = pname + " sn" + system;
    auto mass = particle_mass[idx];

    ejungwoo::gheader(particle_name[idx]+"_");

    TString val_pt    = Form("%f*CMVector[%d].fElements.Perp()",1./particle_a[idx],idx);
    TString val_yyCM  = Form("FragRapidity[%d].fElements/%f",idx,yBeamCM);
    TString val_keCM  = Form("(sqrt(CMVector[%d].fElements.Mag2()+%f)-%f)/%f",idx,mass*mass,mass,particle_a[idx]);
    TString val_phiCM = Form("CMVector[%d].fElements.Phi()*TMath::RadToDeg()",idx);
    TString val_ttaCM = Form("CMVector[%d].fElements.Theta()*TMath::RadToDeg()",idx);
    TString val_pLab  = Form("LabVector[%d].fElements.Mag()",idx);
    TString val_dedx  = "dEdx[0].fElements";

    auto hist1 = ejungwoo::tp("yypt"       , tree, val_pt+":"+val_yyCM   ,            cut_a(idx), mttl+ttl_x + ttl_yypt.xyz()              , binning_yyCM, binning_pt);
    auto hist2 = ejungwoo::tp("yypt_cP"    , tree, val_pt+":"+val_yyCM   , www_p(idx)*cut_a(idx), mttl+ttl_p + ttl_yypt.xyz()              , binning_yyCM, binning_pt);
    auto hist3 = ejungwoo::tp("yypt_cPE"   , tree, val_pt+":"+val_yyCM   , www_a(idx)*cut_a(idx), mttl+ttl_a + ttl_yypt.xyz()              , binning_yyCM, binning_pt);
    auto hist6 = ejungwoo::tp("keCM"       , tree, val_keCM              ,            cut_a(idx), mttl+ttl_x + ttl_keCM.xyz()              , binning_keCM);
    auto hist7 = ejungwoo::tp("keCM_cP"    , tree, val_keCM              , www_p(idx)*cut_a(idx), mttl+ttl_p + ttl_keCM.xyz()              , binning_keCM);
    auto hist8 = ejungwoo::tp("keCM_cPE"   , tree, val_keCM              , www_a(idx)*cut_a(idx), mttl+ttl_a + ttl_keCM.xyz()              , binning_keCM);
    auto histb = ejungwoo::tp("phikeCM"    , tree, val_keCM+":"+val_phiCM,            cut_a(idx), mttl+ttl_x + ttl_phiCM.tx()+ttl_keCM.tx(), binning_phiCM, binning_keCM);
    auto histc = ejungwoo::tp("phikeCM_cPE", tree, val_keCM+":"+val_phiCM, www_a(idx)*cut_a(idx), mttl+ttl_a + ttl_phiCM.tx()+ttl_keCM.tx(), binning_phiCM, binning_keCM);
    auto histd = ejungwoo::tp("ttakeCM"    , tree, val_keCM+":"+val_ttaCM,            cut_a(idx), mttl+ttl_x + ttl_ttaCM.tx()+ttl_keCM.tx(), binning_ttaCM, binning_keCM);
    auto histe = ejungwoo::tp("ttakeCM_cPE", tree, val_keCM+":"+val_ttaCM, www_a(idx)*cut_a(idx), mttl+ttl_a + ttl_ttaCM.tx()+ttl_keCM.tx(), binning_ttaCM, binning_keCM);
    auto histf = ejungwoo::tp("pdedx_raw"  , tree, val_dedx+":"+val_pLab ,                TCut(), mttl+ttl_r + ttl_pLab.tx() +ttl_dedx.tx(), binning_pLab, binning_dedx);
    auto histg = ejungwoo::tp("pdedx_cut"  , tree, val_dedx+":"+val_pLab ,            cut_a(idx), mttl+ttl_t + ttl_pLab.tx() +ttl_dedx.tx(), binning_pLab, binning_dedx);

    auto hist4 = (TH2D *) hist2 -> Clone(ejungwoo::fHeader+"yypt_effi"+ejungwoo::fFooter); hist4 -> Divide(hist3); hist4 -> SetTitle(mttl+" efficiency " +ttl_yypt.xyz());
    auto hist5 = (TH2D *) hist2 -> Clone(ejungwoo::fHeader+"yypt_prob"+ejungwoo::fFooter); hist5 -> Divide(hist1); hist5 -> SetTitle(mttl+" probability "+ttl_yypt.xyz());
    auto hist9 = (TH2D *) hist7 -> Clone(ejungwoo::fHeader+"keCM_effi"+ejungwoo::fFooter); hist9 -> Divide(hist8); hist9 -> SetTitle(mttl+" efficiency " +ttl_keCM.xyz());
    auto hista = (TH2D *) hist7 -> Clone(ejungwoo::fHeader+"keCM_prob"+ejungwoo::fFooter); hista -> Divide(hist6); hista -> SetTitle(mttl+" probability "+ttl_keCM.xyz());

    file_hist -> cd();
    for (auto hist : {hist1, hist2, hist3, hist4, hist5, hist6, hist7, hist8, hist9, hista, histb, histc, histd, histe, histf, histg})
      hist -> Write();
  }

  file_hist -> Print();
}
