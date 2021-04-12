#include "init_variables.h"

void make_compact_root_files()
{
  const char *ana_version = "nn50";
  const char *path_to_data = Form("data2/%s",ana_version);

  for (auto iSystem : fSysIdx)
  {
    auto name_file_out = Form("%s/compact_sys%d_%s.root",path_to_data,fSysBeams[iSystem],ana_version);
    cout << name_file_out << endl;
    auto file_out = new TFile(name_file_out,"recreate");

    auto name_file_mult = Form("%s/sys%d_%s_*_50_100.NewAna.2107.4fd2bca.ana.mult.root",path_to_data,fSysBeams[iSystem],ana_version);
    auto tree_mult = new TChain("mult");
    tree_mult -> Add(name_file_mult);
    auto n_events = tree_mult -> GetEntries("1");

    for (auto iParticle : fParticleIdx)
    {
      auto name_particle = fParticleNames[iParticle];
      cout << name_particle << endl;
      auto name_file_in = Form("%s/sys%d_%s_*_50_100.NewAna.2107.4fd2bca.ana.%s.root",path_to_data,fSysBeams[iSystem],ana_version,name_particle);
      auto chain = new TChain(name_particle);
      chain -> Add(name_file_in);

      int particle_a = fParticleA[iParticle];

      short br_nt;
      int br_nr, br_nl, br_nt2;
      double br_prob, br_eff, br_sd, br_ke_cm, br_pt_cm, br_by_cm0, br_by_cm, br_fy_cm, br_fy_lab, br_p_lab, br_dedx;
      double br_p_cm, br_theta_cm, br_phi_cm, br_theta_lab, br_phi_lab, br_dpoca, br_dpoca_r;

      chain -> SetBranchAddress("prob"     , &br_prob     );
      chain -> SetBranchAddress("eff"      , &br_eff      );
      chain -> SetBranchAddress("sd"       , &br_sd       );
      chain -> SetBranchAddress("ke_cm"    , &br_ke_cm    );
      chain -> SetBranchAddress("pt_cm"    , &br_pt_cm    );
      chain -> SetBranchAddress("by_cm0"   , &br_by_cm0   );
      chain -> SetBranchAddress("by_cm"    , &br_by_cm    );
      chain -> SetBranchAddress("fy_cm"    , &br_fy_cm    );
      chain -> SetBranchAddress("fy_lab"   , &br_fy_lab   );
      chain -> SetBranchAddress("p_lab"    , &br_p_lab    );
      chain -> SetBranchAddress("dedx"     , &br_dedx     );
      chain -> SetBranchAddress("p_cm"     , &br_p_cm     );
      chain -> SetBranchAddress("theta_cm" , &br_theta_cm );
      chain -> SetBranchAddress("phi_cm"   , &br_phi_cm   );
      chain -> SetBranchAddress("theta_lab", &br_theta_lab);
      chain -> SetBranchAddress("phi_lab"  , &br_phi_lab  );
      chain -> SetBranchAddress("dpoca"    , &br_dpoca    );
      chain -> SetBranchAddress("dpoca_r"  , &br_dpoca_r  );
      chain -> SetBranchAddress("nr"       , &br_nr       );
      chain -> SetBranchAddress("nl"       , &br_nl       );
      chain -> SetBranchAddress("nt"       , &br_nt       );
      chain -> SetBranchAddress("nt2"      , &br_nt2      );

      file_out -> cd();
      double br_keoac, br_ptoac, br_y0, br_ttc, br_corr;
      auto tree_out = new TTree(name_particle,"");
      tree_out -> Branch("prob"     , &br_prob  );
      tree_out -> Branch("corr"     , &br_corr  );
      tree_out -> Branch("keoac"    , &br_keoac );
      tree_out -> Branch("ptoac"    , &br_ptoac );
      tree_out -> Branch("y0"       , &br_y0    );
      tree_out -> Branch("ttc"      , &br_ttc   );
      tree_out -> Branch("pozl"     , &br_p_lab );
      tree_out -> Branch("dedx"     , &br_dedx  );

      auto entries = chain -> GetEntries();
      auto i0 = entries; i0=0;

      for (auto i=i0; i<entries; ++i)
      {
        chain -> GetEntry(i);

        if (br_prob < 0.5) continue;
        if (br_eff < 0.05) continue;
        if (abs(br_sd) > 3) continue;
        if (br_dpoca > 15) continue;
        if (br_nr+br_nl < 15) continue;

             if (iParticle==0) { if (br_p_lab<100) continue; }
        else if (iParticle==1) { if (br_p_lab<200) continue; }
        else if (iParticle==2) { if (br_p_lab<200) continue; }
        else if (iParticle==3) { if (br_p_lab<350) continue; }
        else if (iParticle==4) { if (br_p_lab<400) continue; }

        br_keoac = br_ke_cm / particle_a;
        br_ptoac = br_pt_cm / particle_a;
        br_y0 = br_fy_cm / (br_by_cm/2);
        br_ttc = br_theta_cm * TMath::RadToDeg();
        br_corr = br_prob/br_eff/n_events;

        tree_out -> Fill();
      }
      tree_out -> Write();
    }
  }
}
