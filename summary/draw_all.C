#include "init_variables.h"

void draw_all()
{
  Init();

  for (auto multOption : fMultOptions) {
    for (auto system : fSystems) {
      for (auto iParticle : fParticleIdx)
      {
        TString particleName = fParticleNames[iParticle];

        auto fileName = Form("data/sys%d.all_%d_%d.NewAna.2100.de29163.ana.%s.root", system, fMultLL[multOption], fMultHL[multOption], particleName.Data());
        auto file = new TFile(fileName,"read");
        auto tree = (TTree *) file -> Get("data");

        ecanvas *cvs = new ecanvas();

        for (auto hist : fHists1)
        {
          hist -> addOption("setstat=0,addtol=0");
          hist -> setMainTitle(particleName + " (" + system +")");
          hist -> setPar(0, fParticleA[iParticle]);

          hist -> make(tree);
          cvs -> addnext(hist);
        }

        for (auto hist : fHists2)
        {
          hist -> addOption("setstat=0,addtol=0");
          hist -> setMainTitle(particleName + " (" + system +")");
          hist -> setPar(0, fParticleA[iParticle]);

          hist -> make(tree);
          cvs -> addnext(hist);
        }

        cvs -> draw();
      }
      return;
    }
    return;
  }
}
