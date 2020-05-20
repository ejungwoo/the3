void read_ana(
    int idx=0,
    bool is132or108=1,
    bool debug=false
)
{
  const char *fileNames[] = { // 23
    "/home/ejungwoo/data/pid2/beam132_run2841_2846_ana.root",
    "/home/ejungwoo/data/pid2/beam132_run2848_2852_ana.root",
    "/home/ejungwoo/data/pid2/beam132_run2855_2859_ana.root",
    "/home/ejungwoo/data/pid2/beam132_run2860_2878_ana.root",
    "/home/ejungwoo/data/pid2/beam132_run2879_2883_ana.root",
    "/home/ejungwoo/data/pid2/beam132_run2884_2890_ana.root",
    "/home/ejungwoo/data/pid2/beam132_run2891_2896_ana.root",
    "/home/ejungwoo/data/pid2/beam132_run2898_2902_ana.root",
    "/home/ejungwoo/data/pid2/beam132_run2903_2914_ana.root",
    "/home/ejungwoo/data/pid2/beam132_run2916_2921_ana.root", //
    "/home/ejungwoo/data/pid2/beam132_run2922_2927_ana.root",
    "/home/ejungwoo/data/pid2/beam132_run2929_2933_ana.root",
    "/home/ejungwoo/data/pid2/beam132_run2934_2940_ana.root",
    "/home/ejungwoo/data/pid2/beam132_run2941_2945_ana.root",
    "/home/ejungwoo/data/pid2/beam132_run2946_2958_ana.root",
    "/home/ejungwoo/data/pid2/beam132_run2959_2964_ana.root",
    "/home/ejungwoo/data/pid2/beam132_run2965_2970_ana.root",
    "/home/ejungwoo/data/pid2/beam132_run2971_2976_ana.root",
    "/home/ejungwoo/data/pid2/beam132_run2977_2981_ana.root",
    "/home/ejungwoo/data/pid2/beam132_run2982_2986_ana.root",
    "/home/ejungwoo/data/pid2/beam132_run2988_2992_ana.root",
    "/home/ejungwoo/data/pid2/beam132_run2993_3002_ana.root",
    "/home/ejungwoo/data/pid2/beam132_run3003_3039_ana.root",
  };

  const char *fileNames2[] = { // 18
    "/home/ejungwoo/data/pid2/beam108_run2272_2276_ana.root",
    "/home/ejungwoo/data/pid2/beam108_run2283_2288_ana.root",
    "/home/ejungwoo/data/pid2/beam108_run2289_2314_ana.root",
    "/home/ejungwoo/data/pid2/beam108_run2315_2324_ana.root",
    "/home/ejungwoo/data/pid2/beam108_run2325_2334_ana.root",
    "/home/ejungwoo/data/pid2/beam108_run2335_2341_ana.root",
    "/home/ejungwoo/data/pid2/beam108_run2362_2370_ana.root",
    "/home/ejungwoo/data/pid2/beam108_run2371_2375_ana.root",
    "/home/ejungwoo/data/pid2/beam108_run2378_2382_ana.root",
    "/home/ejungwoo/data/pid2/beam108_run2383_2387_ana.root",
    "/home/ejungwoo/data/pid2/beam108_run2388_2393_ana.root",
    "/home/ejungwoo/data/pid2/beam108_run2394_2398_ana.root",
    "/home/ejungwoo/data/pid2/beam108_run2399_2429_ana.root",
    "/home/ejungwoo/data/pid2/beam108_run2432_2438_ana.root",
    "/home/ejungwoo/data/pid2/beam108_run2439_2461_ana.root",
    "/home/ejungwoo/data/pid2/beam108_run2462_2503_ana.root",
    "/home/ejungwoo/data/pid2/beam108_run2505_2509_ana.root",
    "/home/ejungwoo/data/pid2/beam108_run2505_2509_ana.root",
  };




  auto numPDGs = 6;
  const std::vector<int> listPDGs{2212, 1000010020, 1000010030, 1000020030, 1000020040, 1000020060};

  TFile *file = nullptr;

  if (is132or108)
    file = new TFile(fileNames[idx],"read");
  else
    file = new TFile(fileNames2[idx],"read");

  auto tree = (TTree *) file -> Get("cbmsim");

  TClonesArray *dataArray = nullptr;
  TClonesArray *probArray = nullptr;
  TClonesArray *effiArray = nullptr;
  TClonesArray *momcArray = nullptr;
  TClonesArray *frapArray = nullptr;
  TClonesArray *cmkeArray = nullptr;
  STVectorF *brapArray = nullptr;
  tree -> SetBranchAddress("STData",&dataArray);
  tree -> SetBranchAddress("Prob",&probArray);
  tree -> SetBranchAddress("Eff",&effiArray);
  tree -> SetBranchAddress("CMVector",&momcArray);
  tree -> SetBranchAddress("FragRapidity",&frapArray);
  tree -> SetBranchAddress("BeamRapidity",&brapArray);
  tree -> SetBranchAddress("CMKE",&cmkeArray);

  std::map<int, STVectorF*> probMap;
  std::map<int, STVectorF*> effiMap;
  std::map<int, STVectorVec3*> momcMap;
  std::map<int, STVectorF*> frapMap;
  std::map<int, STVectorVec3*> cmkeMap;

  int pdgMP = 0;
  bool isGood = false;
  double probMP = 0.;
  double effiMP = 0.;
  double frapMP = 0.;
  double ptMP = 0.;
  double keMP = 0.;

  double brap = 0.;

  TFile *fileOut = nullptr;
  if (is132or108)
    fileOut = new TFile(Form("data/summary_132_%d.root",idx),"recreate");
  else
    fileOut = new TFile(Form("data/summary_108_%d.root",idx),"recreate");

  auto treeOut = new TTree("data","");
  treeOut -> Branch("good",&isGood);
  treeOut -> Branch("pdg",&pdgMP);
  treeOut -> Branch("prob",&probMP);
  treeOut -> Branch("eff",&effiMP);
  treeOut -> Branch("y",&frapMP);
  treeOut -> Branch("yBeam",&brap);
  treeOut -> Branch("pt",&ptMP);
  treeOut -> Branch("ke",&keMP);

  auto numEvents = tree -> GetEntries();
  cout << "Number of events: " << numEvents << endl;
  for (auto event=0; event<numEvents; ++event)
  {
    tree -> GetEntry(event);

    brap = brapArray -> fElements[0];
    //brap = brapArray -> fElements[1];

    if (event%1000==0)
      cout << "event:" << event << "/" << numEvents << ",  filled:" << treeOut -> GetEntries() << endl;

    for (auto iPDG=0; iPDG<numPDGs; ++iPDG) {
      auto pdg = listPDGs[iPDG];
      probMap[pdg] = static_cast<STVectorF*>(probArray->At(iPDG));
      effiMap[pdg] = static_cast<STVectorF*>(effiArray->At(iPDG));
      momcMap[pdg] = static_cast<STVectorVec3*>(momcArray->At(iPDG));
      frapMap[pdg] = static_cast<STVectorF*>(frapArray->At(iPDG));
      cmkeMap[pdg] = static_cast<STVectorVec3*>(cmkeArray->At(iPDG));
    }

    auto data = (STData *) dataArray -> At(0);
    auto numTracks = data -> multiplicity;
    //cout << "numTracks in event-" << event << ": " << numTracks << endl;

    for (auto iTrack=0; iTrack<numTracks; ++iTrack)
    {
      pdgMP = 0;
      probMP = 0.;
      effiMP = 0.;
      frapMP = 0.;
      ptMP = 0.;
      keMP = 0.;

      for (auto pdg : listPDGs)
      {
        auto prob = probMap[pdg] -> fElements.at(iTrack);
        auto effi = effiMap[pdg] -> fElements.at(iTrack);

        auto frap = frapMap[pdg] -> fElements[iTrack];
        auto momc = momcMap[pdg] -> fElements[iTrack];
        auto cmke = cmkeMap[pdg] -> fElements[iTrack].z();
        double pt = momc.Perp();
        //double z = momc.z();

        if (probMP < prob) {
          pdgMP = pdg;
          probMP = prob;
          effiMP = effi;
          frapMP = frap;
          ptMP = pt;
          keMP = cmke;
        }
      }

      if (probMP > 0.8 && effiMP != 0 )
        isGood = true;
      else
        isGood = false;

      if (debug) {
        if (isGood) cout << right << "O " << setw(5) << iTrack << " " << setw(12) << pdgMP << left << " p:" << setw(12) << probMP << " e:" << setw(12) << effiMP << " y:" << setw(12) << frapMP << " beam_y:" << setw(12) << brap << endl;
        else        cout << right << "X " << setw(5) << iTrack << " " << setw(12) << pdgMP << left << " p:" << setw(12) << probMP << " e:" << setw(12) << effiMP << " y:" << setw(12) << frapMP << " beam_y:" << setw(12) << brap << endl;
      }

      if (probMP!=0)
        treeOut -> Fill();
    }

    if (debug) if (event>5) return;
  }

  fileOut -> cd();
  cout << fileOut -> GetName() << " " << treeOut -> GetEntries() << endl;
  treeOut -> Write();
}
