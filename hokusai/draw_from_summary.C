#include "/home/ejungwoo/config/ejungwoo.h"

void draw_from_summary(
  int system = 108,
  int version = 1
  )
{
  ejungwoo::gsave(0);
  ejungwoo::gversion("vvv");
  
  bool do_draw_angle = true;
  bool do_draw_yypt = true;
  bool do_draw_keCM = true;

  const char *pnames[] = {"p","d","t","he3","he4"};

  auto n1 = 20;
  auto n2 = 20;
  auto n3 = 20;


  TFile *files[20];
  for (auto iSplit=0; iSplit<n1; ++iSplit) {
    //files[iSplit] = new TFile(Form("data_from_hokusai_xml2/summary2_hist_%d_s%d.root",system,iSplit));
    //files[iSplit] = new TFile(Form("data_xml/summary3_hist_%d_s%d.root",system,iSplit));
    //files[iSplit] -> ls();
    files[iSplit] = new TFile(Form("data_xml/summary_hist_%d_s%d_v%d.root",system,iSplit,version));
  }

  TString file_name = Form("data_xml/summary_hist_%d_all_v%d.root",system,version);
  auto file_out = new TFile(file_name,"recreate");

  if (do_draw_yypt) {
    TH2D *hsum1[5] = {nullptr,nullptr,nullptr,nullptr,nullptr};
    TH2D *hsum2[5] = {nullptr,nullptr,nullptr,nullptr,nullptr};
    TH2D *hsum3[5] = {nullptr,nullptr,nullptr,nullptr,nullptr};
    TH2D *hsum4[5] = {nullptr,nullptr,nullptr,nullptr,nullptr};
    TH2D *hsum5[5] = {nullptr,nullptr,nullptr,nullptr,nullptr};

    for (auto iSplit=0; iSplit<n1; ++iSplit) {
      auto file = files[iSplit];
      for (auto idx : {0,1,2,3,4})
      {
        auto hist1 = (TH2D *) file -> Get(Form("yypt_cPE_%s" ,pnames[idx])); if (hsum1[idx]==nullptr) hsum1[idx] = hist1; else hsum1[idx] -> Add(hist1);
        auto hist2 = (TH2D *) file -> Get(Form("yypt_cP_%s"  ,pnames[idx])); if (hsum2[idx]==nullptr) hsum2[idx] = hist2; else hsum2[idx] -> Add(hist2);
        auto hist3 = (TH2D *) file -> Get(Form("yypt_%s"     ,pnames[idx])); if (hsum3[idx]==nullptr) hsum3[idx] = hist3; else hsum3[idx] -> Add(hist3);
        auto hist4 = (TH2D *) file -> Get(Form("yypt_eff_%s" ,pnames[idx])); if (hsum4[idx]==nullptr) hsum4[idx] = hist4; //else hsum4[idx] -> Add(hist4);
        auto hist5 = (TH2D *) file -> Get(Form("yypt_prob_%s",pnames[idx])); if (hsum5[idx]==nullptr) hsum5[idx] = hist5; //else hsum5[idx] -> Add(hist5);
      }
    }

    for (auto idx : {0,1,2,3,4})
    {
      auto cvs1 = ejungwoo::cc(TString("canvas1_")+hsum1[idx]->GetName()); hsum1[idx] -> Draw("colz"); ejungwoo::make_c(cvs1); ejungwoo::save(cvs1,"png"); file_out -> cd(); hsum1[idx] -> Write();
      auto cvs2 = ejungwoo::cc(TString("canvas1_")+hsum2[idx]->GetName()); hsum2[idx] -> Draw("colz"); ejungwoo::make_c(cvs2); ejungwoo::save(cvs2,"png"); file_out -> cd(); hsum2[idx] -> Write();
      auto cvs3 = ejungwoo::cc(TString("canvas1_")+hsum3[idx]->GetName()); hsum3[idx] -> Draw("colz"); ejungwoo::make_c(cvs3); ejungwoo::save(cvs3,"png"); file_out -> cd(); hsum3[idx] -> Write();
      auto cvs4 = ejungwoo::cc(TString("canvas1_")+hsum4[idx]->GetName()); hsum4[idx] -> Draw("colz"); ejungwoo::make_c(cvs4); ejungwoo::save(cvs4,"png"); file_out -> cd(); hsum4[idx] -> Write();
      auto cvs5 = ejungwoo::cc(TString("canvas1_")+hsum5[idx]->GetName()); hsum5[idx] -> Draw("colz"); ejungwoo::make_c(cvs5); ejungwoo::save(cvs5,"png"); file_out -> cd(); hsum5[idx] -> Write();
    }
  }

  if (do_draw_keCM) {
    TH1D *hsum1[5] = {nullptr,nullptr,nullptr,nullptr,nullptr};
    TH1D *hsum2[5] = {nullptr,nullptr,nullptr,nullptr,nullptr};
    TH1D *hsum3[5] = {nullptr,nullptr,nullptr,nullptr,nullptr};
    TH1D *hsum4[5] = {nullptr,nullptr,nullptr,nullptr,nullptr};
    TH1D *hsum5[5] = {nullptr,nullptr,nullptr,nullptr,nullptr};

    for (auto iSplit=0; iSplit<n2; ++iSplit) {
      auto file = files[iSplit];
      for (auto idx : {0,1,2,3,4})
      {
        auto hist1 = (TH1D *) file -> Get(Form("keCM_cPE_%s" ,pnames[idx])); if (hsum1[idx]==nullptr) hsum1[idx] = hist1; else hsum1[idx] -> Add(hist1);
        auto hist2 = (TH1D *) file -> Get(Form("keCM_cP_%s"  ,pnames[idx])); if (hsum2[idx]==nullptr) hsum2[idx] = hist2; else hsum2[idx] -> Add(hist2);
        auto hist3 = (TH1D *) file -> Get(Form("keCM_%s"     ,pnames[idx])); if (hsum3[idx]==nullptr) hsum3[idx] = hist3; else hsum3[idx] -> Add(hist3);
        auto hist4 = (TH1D *) file -> Get(Form("keCM_eff_%s" ,pnames[idx])); if (hsum4[idx]==nullptr) hsum4[idx] = hist4; //else hsum4[idx] -> Add(hist4);
        auto hist5 = (TH1D *) file -> Get(Form("keCM_prob_%s",pnames[idx])); if (hsum5[idx]==nullptr) hsum5[idx] = hist5; //else hsum5[idx] -> Add(hist5);
      }
    }

    for (auto idx : {0,1,2,3,4})
    {
      auto cvs1 = ejungwoo::cv(TString("canvas2_")+hsum1[idx]->GetName()); hsum1[idx] -> Draw(); ejungwoo::make_c(cvs1); ejungwoo::save(cvs1,"png"); file_out -> cd(); hsum1[idx] -> Write();
      auto cvs2 = ejungwoo::cv(TString("canvas2_")+hsum2[idx]->GetName()); hsum2[idx] -> Draw(); ejungwoo::make_c(cvs2); ejungwoo::save(cvs2,"png"); file_out -> cd(); hsum2[idx] -> Write();
      auto cvs3 = ejungwoo::cv(TString("canvas2_")+hsum3[idx]->GetName()); hsum3[idx] -> Draw(); ejungwoo::make_c(cvs3); ejungwoo::save(cvs3,"png"); file_out -> cd(); hsum3[idx] -> Write();
      auto cvs4 = ejungwoo::cv(TString("canvas2_")+hsum4[idx]->GetName()); hsum4[idx] -> Draw(); ejungwoo::make_c(cvs4); ejungwoo::save(cvs4,"png"); file_out -> cd(); hsum4[idx] -> Write();
      auto cvs5 = ejungwoo::cv(TString("canvas2_")+hsum5[idx]->GetName()); hsum5[idx] -> Draw(); ejungwoo::make_c(cvs5); ejungwoo::save(cvs5,"png"); file_out -> cd(); hsum5[idx] -> Write();
    }
  }

  if (do_draw_angle) {
    TH1D *hsum1[5] = {nullptr,nullptr,nullptr,nullptr,nullptr};
    TH1D *hsum2[5] = {nullptr,nullptr,nullptr,nullptr,nullptr};
    TH1D *hsum3[5] = {nullptr,nullptr,nullptr,nullptr,nullptr};
    TH1D *hsum4[5] = {nullptr,nullptr,nullptr,nullptr,nullptr};

    for (auto iSplit=0; iSplit<n3; ++iSplit) {
      auto file = files[iSplit];
      for (auto idx : {0,1,2,3,4})
      {
        auto hist1 = (TH1D *) file -> Get(Form("phi_cPE_%s"  ,pnames[idx])); if (hsum1[idx]==nullptr) hsum1[idx] = hist1; else hsum1[idx] -> Add(hist1);
        auto hist2 = (TH1D *) file -> Get(Form("phi_%s"      ,pnames[idx])); if (hsum2[idx]==nullptr) hsum2[idx] = hist2; else hsum2[idx] -> Add(hist2);
        auto hist3 = (TH1D *) file -> Get(Form("theta_cPE_%s",pnames[idx])); if (hsum3[idx]==nullptr) hsum3[idx] = hist3; else hsum3[idx] -> Add(hist3);
        auto hist4 = (TH1D *) file -> Get(Form("theta_%s"    ,pnames[idx])); if (hsum4[idx]==nullptr) hsum4[idx] = hist4; else hsum4[idx] -> Add(hist4);
      }
    }

    for (auto idx : {0,1,2,3,4})
    {
      auto cvs1 = ejungwoo::cv(TString("canvas3_")+hsum1[idx]->GetName()); hsum1[idx] -> Draw("hist"); ejungwoo::make_c(cvs1); ejungwoo::save(cvs1,"png"); file_out -> cd(); hsum1[idx] -> Write();
      auto cvs2 = ejungwoo::cv(TString("canvas3_")+hsum2[idx]->GetName()); hsum2[idx] -> Draw("hist"); ejungwoo::make_c(cvs2); ejungwoo::save(cvs2,"png"); file_out -> cd(); hsum2[idx] -> Write();
      auto cvs3 = ejungwoo::cv(TString("canvas3_")+hsum3[idx]->GetName()); hsum3[idx] -> Draw("hist"); ejungwoo::make_c(cvs3); ejungwoo::save(cvs3,"png"); file_out -> cd(); hsum3[idx] -> Write();
      auto cvs4 = ejungwoo::cv(TString("canvas3_")+hsum4[idx]->GetName()); hsum4[idx] -> Draw("hist"); ejungwoo::make_c(cvs4); ejungwoo::save(cvs4,"png"); file_out -> cd(); hsum4[idx] -> Write();
    }
  }
  
  cout << file_name  << endl;
}
