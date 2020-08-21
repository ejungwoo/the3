#include "/home/ejungwoo/config/ejungwoo.h"
#include "STConcReader.C"

void make_submit_embed()
{
  int multiplicityCut = 50;
  int numEventsEmbed = 200;
  int numEventsForTest = 10;

  TString workDir = gSystem -> Getenv("VMCWORKDIR");
  TString parDir = workDir + "/parameters/";
  TString pwdDir = TString(gSystem -> pwd()) + "/";
  TString confDir = pwdDir + "analysisInputFiles/embedConfig/";
  TString mcDir = "/home/ejungwoo/data/mc/";
  TString digiDir = "/home/ejungwoo/data/digi/";
  TString embedDir = "/home/ejungwoo/data/embed/";

  auto tree = new TChain("spirit");
  ifstream listFile("analysisInputFiles/listFiles/list_reco_20200529_132.txt");
  TString concName;
  vector<int> concRunArray;
  for (auto i : {0,1,2}) {
    listFile >> concName;
    cout << concName << endl;
    tree -> Add(concName);
    TString justRun = ejungwoo::tok(ejungwoo::justname(concName),"_",0);
    justRun.ReplaceAll("run","");
    if (concRunArray.size()==0) concRunArray.push_back(justRun.Atoi());
    //else if (concRunArray.back()!=justRun.Atoi()) concRunArray.push_back(justRun.Atoi());
    else if (concRunArray.back()!=justRun.Atoi()) break;
  }

  TString concFilesName;
       if (concRunArray.size()==0) concFilesName = "x";
  else if (concRunArray.size()==1) concFilesName = Form("run%d",concRunArray[0]);
  else if (concRunArray.size() >1) concFilesName = Form("run%dTo%d",concRunArray[0],concRunArray.back());

  TString listEVName = TString("eventVertexList__") + numEventsEmbed + "__" + concFilesName + ".txt";
  TString listEVFull = Form("%s/%s",parDir.Data(),listEVName.Data());
  ofstream eventVertexListFile(listEVFull);
  eventVertexListFile << "#RunNum EventNum x(mm) y(mm) z(mm), " << ejungwoo::lastname(concName) << " with multiplicity cut:" << multiplicityCut << endl;

  TString listSEFull = confDir + TString("selectedEvents__") + numEventsEmbed + "__" + concFilesName + ".txt"; 
  ofstream selectedEventListFile(listSEFull);
  vector<int> selectedEventPlus1Array;

  auto reader = new STConcReader(tree);
  auto numEvents = tree -> GetEntries();
  int countGoodEvents = 0;
  for (auto iEvent=0; iEvent<numEvents; ++iEvent) {
    reader -> GetEntry(iEvent);
    auto pos = reader -> tpcVertex;
    if (reader -> vaMultiplicity >= multiplicityCut) {
      eventVertexListFile << concRunArray[0] << " " << iEvent << " " << .1*pos.X() << " " << .1*pos.Y() << " " << .1*pos.Z() << endl;
      selectedEventListFile << iEvent << endl;
      selectedEventPlus1Array.push_back(iEvent+1);
      countGoodEvents++;
    }
    if (countGoodEvents >= numEventsEmbed)
      break;
  }

  TString eventArrayString = "\\{";
  for (auto event : selectedEventPlus1Array)
    eventArrayString = eventArrayString + event + ",";
  eventArrayString = eventArrayString + "\\}";
  cout << eventArrayString << endl;

  cout << listEVFull << " " << countGoodEvents << endl;

  int numMomVaules = 8.;
  vector<vector<int>> momRangeArray = {
    vector<int>{ 100,1500},
    vector<int>{ 200,2000},
    vector<int>{ 400,3000},
    vector<int>{ 400,1500},
    vector<int>{ 500,2500},
    vector<int>{1500,2500},
  };

  int dtta = 15.;
  vector<int> ttaRange = {0,90};

  int dphi = 10.;
  vector<vector<int>> phiRangeArray = {
    vector<int>{160,180},
    vector<int>{0,20}
  }; 

  bool makeTest = true;

  const int numParticles = 5;
  const int particlePDGs[] = {2212, 1000010030, 1000010020, 1000020040, 1000020030,};
  TString particleNames[] = {"p","d","t","he3","he4","he6"};

  TString runName = "submit_embed";
  TString runTag = "se";

  TString subDir = pwdDir + runName + "/";
  TString outDir = subDir + "job_out/";
  TString logDir = subDir + "job_log/";
  TString allFull = pwdDir + runName+".sh";
  TString endFull = pwdDir + runName+".end";

  std::ofstream submit_all(allFull);
  submit_all << "set +x" << endl;

  std::ofstream submit_end(endFull);

  int countConfFiles = 0;
  for (auto ipdg=0; ipdg<numParticles; ipdg++) {
    auto pdg = particlePDGs[ipdg];
    const char *pname = particleNames[ipdg];
    auto momRange =  momRangeArray[ipdg];
    auto mom0 = momRange[0];
    auto mom1 = momRange[1];
    auto dmom = double(momRange[1]-momRange[0])/numMomVaules;
    for (auto imom=0; imom<numMomVaules; imom++) {
      auto mom = mom0 + .5*dmom + imom*dmom;
      for (auto tta=ttaRange[0]+.5*dtta; tta<=ttaRange[1]; tta+=dtta) {
        for (auto phiRange : phiRangeArray) {
          for (auto phi=phiRange[0]+.5*dphi; phi<=phiRange[1]; phi+=dphi)
          {
            TString confName = Form("embed_%s_mom%d_phi%d_tta%d",pname,imom,int(phi-.5*dphi),int(tta-.5*dtta));
            TString embedMCDigiFull = digiDir + confName + ".digi.root";
            TString confFull = Form("%s/%s.conf.par",confDir.Data(),confName.Data());
            ofstream confFile(confFull);
            confFile << "VertexFile " << listEVName << endl;
            confFile << "Particle " << pdg << endl;
            confFile << "Momentum " << mom*0.001 << endl;
            confFile << "Theta " << tta << endl;
            confFile << "Phi " << phi << endl;

            TString subName = Form("%s%d",runTag.Data(),countConfFiles);
            TString macFull = subDir + subName + ".sh";
            TString outFull = outDir + subName + ".out";
            TString lmcFull = logDir + subName + ".mc.log";
            TString tailmc9 = TString("$(tail -9 ") + lmcFull + " | head -1)";
            TString ldgFull = logDir + subName + ".digi.log";
            TString taildg4 = TString("$(tail -4 ") + ldgFull + " | head -1)";
            TString lrcFull = logDir + subName + ".reco.log";
            TString tailrc1 = TString("$(tail -1 ") + lrcFull + " | head -1)";

            submit_end << subName << " " << lmcFull << endl;
            //cout << macFull << " " << lmcFull << endl;

            std::ofstream submit_macro(macFull);
            submit_all << "pjsub --bulk --sparam 0-0 " << macFull << endl;

            submit_macro << "#!/bin/bash" << endl;
            submit_macro << "#------ pjsub option -------- #" << endl;
            submit_macro << "#PJM -L rscunit=bwmpc" << endl;
            submit_macro << "#PJM -L rscgrp=batch" << endl;
            submit_macro << "#PJM -L vnode=1" << endl;
            submit_macro << "#PJM -L vnode-core=10" << endl;
            submit_macro << "#PJM -L vnode-mem=10Gi" << endl;
            submit_macro << "#PJM -g Q20393" << endl;
            submit_macro << "#PJM -j" << endl;
            submit_macro << "#PJM -o " << outFull << endl;
            submit_macro << "#------- Program execution ------- #" << endl;
            submit_macro << "export OMP_NUM_THREADS=1" << endl;
            submit_macro << "source /home/ejungwoo/environment.spiritroot.bwmpc.sh" << endl;
            submit_macro << "cd " << subDir << endl;

            submit_macro << "root -q -b -l run_mc.C\\("
                         << "\\\"" << mcDir << "\\\","
                         << "\\\"" << confName << "\\\","
                         << "-1,"
                         << "\\\"" << confFull << "\\\"" 
                         << "\\) > " << lmcFull << " 2>&1 " << endl;
            submit_macro << "echo '"<< confName << " ' " << tailmc9 << " >> " << endFull << endl;

            submit_macro << "root -q -b -l run_digi.C\\("
                         << "\\\"" << mcDir << "\\\","
                         << "\\\"" << digiDir << "\\\","
                         << "\\\"" << confName << "\\\""
                         << "\\) > " << ldgFull << " 2>&1 " << endl;
            submit_macro << "echo '"<< confName << " ' " << taildg4 << " >> " << endFull << endl;

            submit_macro << "root -q -b -l run_reco_embed.C\\("
                         << "\\\"" << confName << "\\\","
                         //<< concRunArray[0] << ",0," << selectedEventPlus1Array.back() << ","
                         << concRunArray[0] << ",0," << numEventsEmbed << ","
                         << eventArrayString << ","
                         << "\\\"" << embedMCDigiFull << "\\\","
                         << "\\\"" << embedDir << "\\\""
                         << "\\) > " << lrcFull << " 2>&1 " << endl;
            submit_macro << "echo '"<< confName << " ' " << tailrc1 << " >> " << endFull << endl;

            if (makeTest) {
              TString testFull = pwdDir + runName+".test.sh";
              std::ofstream submit_test(testFull);

              submit_test << "root -q -b -l run_mc.C\\("
                          << "\\\"" << mcDir << "\\\","
                          << "\\\"" << confName << "\\\","
                          << numEventsForTest << ","
                          << "\\\"" << confFull << "\\\"" 
                          << "\\);" << endl;

              submit_test << "root -q -b -l run_digi.C\\("
                          << "\\\"" << mcDir << "\\\","
                          << "\\\"" << digiDir << "\\\","
                          << "\\\"" << confName << "\\\""
                          << "\\);" << endl;

              submit_test << "root -q -b -l run_reco_embed.C\\("
                          << "\\\"" << confName << "\\\","
                          //<< concRunArray[0] << ",0," << selectedEventPlus1Array[numEventsForTest-1] << ","
                          << concRunArray[0] << ",0," << numEventsForTest << ","
                          << eventArrayString << ","
                          << "\\\"" << embedMCDigiFull << "\\\","
                          << "\\\"" << embedDir << "\\\""
                          << "\\);" << endl;

              makeTest = false;
            }

            countConfFiles++;
          }
        }
      }
    }
  }

  cout << endl;
  cout << allFull << endl;
  cout << endFull << endl;
}
