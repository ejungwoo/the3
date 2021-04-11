#include "init_variables.h"

void draw_quick_summary()
{
  const char *name_version = "nn50";
  const char *path_to_data = Form("data2/%s",name_version);

  for (auto iSystem : fSysIdx)
  {
    for (auto iParticle : fParticleIdx)
    {
      auto name_particle = fParticleNames[iParticle];
      auto name_file_in = Form("%s/sys%d_%s_*_50_100.NewAna.2107.4fd2bca.ana.%s.root",path_to_data,fSysBeams[iSystem],name_version,name_particle);
      auto chain = new TChain(name_particle);
      chain -> Add(name_file_in);

      auto cvs = new TCanvas(name_particle,"",2000,800);
      cvs -> Divide(3,2);

      auto idx = 1;
      cvs -> cd(idx); chain -> Draw(Form("prob>>%s_%d",name_particle,idx) ); idx++;
      cvs -> cd(idx); chain -> Draw(Form("eff>>%s_%d",name_particle,idx)  ); idx++;
      cvs -> cd(idx); chain -> Draw(Form("sd>>%s_%d",name_particle,idx)   ); idx++;
      cvs -> cd(idx); chain -> Draw(Form("dpoca>>%s_%d",name_particle,idx)); idx++;
      cvs -> cd(idx); chain -> Draw(Form("nr+nl>>%s_%d",name_particle,idx)); idx++;
    }
    return;
  }
}
