#include "init_variables.h"

void draw_kk()
{
  auto file = new TFile("data2/nn50/compact_sys132_nn50.root");
  auto cvs_kk = new TCanvas("kt","",2000,800);
  cvs_kk -> Divide(5,2);

  TCut cut("keoac>10");

  for (auto iParticle : fParticleIdx)
  {
    auto name_particle = fParticleNames[iParticle];
    auto tree = (TTree *) file -> Get(name_particle);

    cvs_kk -> cd(iParticle+1);
    auto name_kt = Form("kt_%s",name_particle);
    auto kt = new TH2D(name_kt,"",100,0,200,100,0,180);
    tree -> Draw(Form("ttc:keoac>>%s",name_kt), cut, "colz");
    //auto kt = new TH1D(name_kt,"",100,0,100);
    //tree -> Draw(Form("keoac>>%s",name_kt), "", "colz");

    cvs_kk -> cd(5+iParticle+1);
    auto name_yp = Form("yp_%s",name_particle);
    auto yp = new TH2D(name_yp,"",100,-1,2,100,0,1000);
    tree -> Draw(Form("ptoac:y0>>%s",name_yp), cut, "colz");
    //auto yp = new TH1D(name_yp,"",100,0,500);
    //tree -> Draw(Form("ptoac>>%s",name_yp), "", "colz");
  }
}
