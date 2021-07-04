#include "init_variables.h"
using namespace ejungwoo;

void draw_quick_summary()
{
  gStyle -> SetOptStat(0);

  //const char *version = "vb3";
  const char *version = "rb3";

  vector<binning> binnings;
  //binning bn_prob(100,0,1.01, "probability", "prob"); binnings.push_back(bn_prob);
  //binning bn_eff(100,0,1.01, "efficienty", "eff"); binnings.push_back(bn_eff);
  binning bn_sd(100,-3,3, "PID residual (#sigma);count/N_{event}", "sd"); binnings.push_back(bn_sd);
  binning bn_dpoca(100,0,10, "distance to vertex (mm);count/N_{event}", "dpoca"); binnings.push_back(bn_dpoca);
  binning bn_nr(150,0,150, "number of row-clusters;count/N_{event}", "nr", "nr!=0"); binnings.push_back(bn_nr);
  binning bn_nl(150,0,150, "number of layer-clusters;count/N_{event}", "nl", "nl!=0"); binnings.push_back(bn_nl);
  binning bn_nc(150,0,150, "number of clusters;count/N_{event}", "nr+nl"); binnings.push_back(bn_nc);
  binning bn_keoa_cm(100,0,400, "KE^{c.m.} (MeV);count/N_{event}", "ke_cm/[A]"); binnings.push_back(bn_keoa_cm);
  binning bn_ptoa_cm(100,0,1000, "p_{T}^{c.m.}/A (MeV/c);count/N_{event}","pt_cm/[A]"); binnings.push_back(bn_ptoa_cm);
  binning bn_y0(100,-1,2, "normalized rapidity y_{0};count/N_{event}", "fy_cm/(by_cm/2)"); binnings.push_back(bn_y0);
  binning bn_phi_cm(100,-180,180, "#phi^{c.m.} (deg);count/N_{event}", "phi_cm*TMath::RadToDeg()"); binnings.push_back(bn_phi_cm);
  binning bn_theta_cm(100,0,180, "#theta^{c.m.} (deg);count/N_{event}", "theta_cm*TMath::RadToDeg()"); binnings.push_back(bn_theta_cm);

  vector<binning2> binnings2;
  //auto bn_phitta = bn_theta_cm*bn_phi_cm; binnings2.push_back(bn_phitta);
  auto bn_y0ptoa = bn_y0*bn_ptoa_cm; binnings2.push_back(bn_y0ptoa);
  auto bn_keoatheta = bn_keoa_cm*bn_theta_cm; binnings2.push_back(bn_keoatheta);

  for (auto iSys : fSysIdx)
  {
    auto sys = fSysBeams[iSys];
    auto name_file_mult = Form("data2/%s/sys%d_%s_*.mult.root",version,sys,version);
    auto chain_mult = new TChain("mult");
    chain_mult -> Add(name_file_mult);
    auto numEvents = chain_mult -> GetEntries();
    for (auto iParticle : fParticleIdx)
    {
      const char *name_particle = fParticleNames[iParticle];
      cout << sys << " " << name_particle << endl;

      auto name_file_in = Form("data2/%s/sys%d_%s_*.%s.root",version,sys,version,name_particle);
      auto chain = new TChain(name_particle);
      chain -> Add(name_file_in);

      auto cvs = canvas(Form("qs_%s_%d_%s",version,sys,name_particle),3,4,"snn");

      int particle_a = fParticleA[iParticle];

      int idx = 0;
      for (auto bnx : binnings) {
        idx++;
        const char *name_hist = Form("hist_%d_%s_%d",sys,name_particle,idx);
        TH1D *hist = bnx.newHist(name_hist);
        TString expression = bnx.getExpression();
        expression.ReplaceAll("[A]",Form("%d",particle_a));
        TCut selection = bnx.getSelection();
        selection = selection * TCut(Form("1/%lld",numEvents));
        chain -> Project(name_hist,expression,selection);
        draw(hist,cvs->cd(idx),"hist");
        if (expression=="prob"||expression=="eff")
          cvs -> cd(idx) -> SetLogy();
      }

      for (auto bnxy : binnings2) {
        idx++;
        const char *name_hist = Form("hist_%d_%s_%d",sys,name_particle,idx);
        TH2D *hist = bnxy.newHist(name_hist);
        TString expression = bnxy.getExpression();
        expression.ReplaceAll("[A]",Form("%d",particle_a));
        TCut selection = bnxy.getSelection();
        selection = selection * TCut(Form("1/%lld",numEvents));
        chain -> Project(name_hist,expression,selection);
        draw(hist,cvs->cd(idx),"col");
        if (expression=="prob"||expression=="eff")
          cvs -> cd(idx) -> SetLogy();
      }
    }
  }

  saveAll(Form("qs_%s",version));
}
