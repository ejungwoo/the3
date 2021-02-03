void create_trees()
{
  TString treeCut = "prob>.5";
  TString cutTag = "p05_";

  bool separate_particles = 1;
  bool merge_45_54 = 1;
  bool merge_all_particles = 1;

  if (separate_particles)
    for (TString fileNameIn : {
        "data/sys108.TAG_45_49.NewAna.2100.de29163.ana.particle.root",
        "data/sys108.TAG_50_54.NewAna.2100.de29163.ana.particle.root",
        "data/sys108.TAG_55_100.NewAna.2100.de29163.ana.particle.root",

        "data/sys112.TAG_45_49.NewAna.2100.de29163.ana.particle.root",
        "data/sys112.TAG_50_54.NewAna.2100.de29163.ana.particle.root",
        "data/sys112.TAG_55_100.NewAna.2100.de29163.ana.particle.root",

        "data/sys124.TAG_45_49.NewAna.2100.de29163.ana.particle.root",
        "data/sys124.TAG_50_54.NewAna.2100.de29163.ana.particle.root",
        "data/sys124.TAG_55_100.NewAna.2100.de29163.ana.particle.root",

        "data/sys132.TAG_45_49.NewAna.2100.de29163.ana.particle.root",
        "data/sys132.TAG_50_54.NewAna.2100.de29163.ana.particle.root",
        "data/sys132.TAG_55_100.NewAna.2100.de29163.ana.particle.root"
        })
  {
    for (auto tag : {"y21_left","y21_right"})
    {
      TString fileName = fileNameIn;
      fileName.ReplaceAll("TAG",tag);
      auto file = new TFile(fileName);

      for (auto treeName : {"p","d","t","he3","he4"})
      {
        auto tree = (TTree *) file -> Get(treeName);
        tree -> SetBranchStatus("nc",0);
        tree -> SetBranchStatus("nl",0);
        tree -> SetBranchStatus("nt",0);

        TString fileName0 = fileName;
        fileName0.ReplaceAll("particle",cutTag+treeName);
        tree -> SetName("data");

        auto file0 = new TFile(fileName0,"recreate");
        auto tree0 = tree -> CopyTree(treeCut);
        tree0 -> Write("data",TObject::kOverwrite);

        cout << fileName0 << endl;
        file0 -> Close();
      }

      for (auto treeName : {"mult"})
      {
        auto tree = (TTree *) file -> Get(treeName);

        TString fileName0 = fileName;
        fileName0.ReplaceAll("particle",treeName);
        tree -> SetName("data");

        auto file0 = new TFile(fileName0,"recreate");
        auto tree0 = tree -> CopyTree("");
        tree0 -> Write("data",TObject::kOverwrite);

        cout << fileName0 << endl;
        file0 -> Close();
      }
    }
  }

  if (merge_45_54)
    for (auto fileNameIn : {TString("data/sysSYSTEM.TAG_MULT.NewAna.2100.de29163.ana.")+cutTag+"PARTICLE.root"})
    {
      for (auto tag : {"y21_left","y21_right"})
      {
        TString fileName = fileNameIn;
        fileName.ReplaceAll("TAG",tag);

        for (auto system : {132,108,124,112})
        {
          for (auto particleName : {"p","d","t","he3","he4","mult"})
          {
            TString fileName0 = fileName;
            TString fileName1 = fileName;
            TString fileName2 = fileName;

            fileName0.ReplaceAll("SYSTEM",Form("%d",system));
            fileName0.ReplaceAll("MULT",Form("%d_%d",45,54));
            fileName0.ReplaceAll("PARTICLE",particleName);

            fileName1.ReplaceAll("SYSTEM",Form("%d",system));
            fileName1.ReplaceAll("MULT",Form("%d_%d",45,49));
            fileName1.ReplaceAll("PARTICLE",particleName);

            fileName2.ReplaceAll("SYSTEM",Form("%d",system));
            fileName2.ReplaceAll("MULT",Form("%d_%d",50,54));
            fileName2.ReplaceAll("PARTICLE",particleName);

            auto chain = new TChain("data");
            chain -> Add(fileName1);
            chain -> Add(fileName2);

            chain -> Merge(fileName0);
            cout << fileName0 << endl;
          }
        }
      }
    }

  if (merge_all_particles)
    for (auto fileNameIn : {TString("data/sysSYSTEM.TAG_MULT.NewAna.2100.de29163.ana.")+cutTag+"PARTICLE.root"})
    {
      for (auto tag : {"y21_left","y21_right"})
      {
        TString fileName = fileNameIn;
        fileName.ReplaceAll("TAG",tag);

        for (auto multString : {"45_54", "55_100"})
        {
          for (auto system : {132,108,124,112})
          {
            TString fileName0 = fileName;
            fileName0.ReplaceAll("SYSTEM",Form("%d",system));
            fileName0.ReplaceAll("MULT",multString);
            fileName0.ReplaceAll("PARTICLE","all");

            auto chain = new TChain("data");

            for (auto particleName : {"p","d","t","he3","he4"})
            {
              TString fileName1 = fileName;

              fileName1.ReplaceAll("SYSTEM",Form("%d",system));
              fileName1.ReplaceAll("MULT",multString);
              fileName1.ReplaceAll("PARTICLE",particleName);

              cout << " ++ " << fileName1 << endl;
              chain -> Add(fileName1);
            }

            chain -> Merge(fileName0);
            cout << fileName0 << endl;

          }
        }
      }
  }
}
