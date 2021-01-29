#include "init_variables.h"

void create_trees()
{
  if (0)
  for (auto fileName : {
      //"data/sys108.all_45_49.NewAna.2100.de29163.ana.particle.root",
      //"data/sys108.all_50_54.NewAna.2100.de29163.ana.particle.root",
      //"data/sys108.all_55_100.NewAna.2100.de29163.ana.particle.root",
      //"data/sys112.all_45_49.NewAna.2100.de29163.ana.particle.root",
      //"data/sys112.all_50_54.NewAna.2100.de29163.ana.particle.root",
      //"data/sys112.all_55_100.NewAna.2100.de29163.ana.particle.root",

      "data/sys124.all_45_49.NewAna.2100.de29163.ana.particle.root",
      "data/sys124.all_50_54.NewAna.2100.de29163.ana.particle.root",
      "data/sys124.all_55_100.NewAna.2100.de29163.ana.particle.root",

      //"data/sys132.all_45_49.NewAna.2100.de29163.ana.particle.root",
      //"data/sys132.all_50_54.NewAna.2100.de29163.ana.particle.root",
      //"data/sys132.all_55_100.NewAna.2100.de29163.ana.particle.root"
      })
  {
    auto file = new TFile(fileName);

    for (auto treeName : {"p","d","t","he3","he4"})
    {
      auto tree = (TTree *) file -> Get(treeName);
      tree -> SetBranchStatus("nc",0);
      tree -> SetBranchStatus("nl",0);
      tree -> SetBranchStatus("nt",0);

      TString fileName0 = fileName;
      fileName0.ReplaceAll("particle",treeName);
      tree -> SetName("data");

      auto file0 = new TFile(fileName0,"recreate");
      auto tree0 = tree -> CopyTree("");
      tree0 -> Write("data",TObject::kOverwrite);

      cout << fileName0 << endl;
      file0 -> ls();
      file0 -> Close();
    }
  }


  else
  for (auto fileName : {"data/sysSYSTEM.all_MULT.NewAna.2100.de29163.ana.PARTICLE.root"})
  {
    //for (auto system : fSystems) {
    for (auto system : {124}) {
      for (auto particleName : fParticleNames)
      {
        TString fileName0 = fileName;
        TString fileName1 = fileName;
        TString fileName2 = fileName;

        fileName0.ReplaceAll("SYSTEM",Form("%d",system));
        fileName0.ReplaceAll("MULT",Form("%d_%d",45,54));
        fileName0.ReplaceAll("PARTICLE",particleName.Data());

        fileName1.ReplaceAll("SYSTEM",Form("%d",system));
        fileName1.ReplaceAll("MULT",Form("%d_%d",45,49));
        fileName1.ReplaceAll("PARTICLE",particleName.Data());

        fileName2.ReplaceAll("SYSTEM",Form("%d",system));
        fileName2.ReplaceAll("MULT",Form("%d_%d",50,54));
        fileName2.ReplaceAll("PARTICLE",particleName.Data());

        auto chain = new TChain("data");
        chain -> Add(fileName1);
        chain -> Add(fileName2);

        chain -> Merge(fileName0);
        cout << fileName0 << endl;
      }
    }
  }
}
