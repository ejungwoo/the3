void make_submit_reco()
{
  //TString anaName = "fix4"; TString readme = "with spiritroot install Tommy's hokusai instruction of Tommy";
  //TString anaName = "fix5"; TString readme = "after checking it works fine";
  TString anaName = "fix6"; TString readme = "vertex fix parameter was not set in fix5";

  bool useDB = 1;
  int numEventsMax = 800000;
  int numEventsInSplit = 5000;
  //auto fSelectRun = 2962; // set 0 to run all
  auto fSelectRun = 0; // set 0 to run all
  auto formSelectRunEvent = "Pick_PiEvt/Sn%d_KanekoEvt/list_all";
  auto nameRunDB = "/home/ejungwoo/spiritroot/parameters/runDB.csv";

  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////

  //int iSystemsForAna[] = {0,1,2,3};
  int iSystemsForAna[] = {0,3};
  int systems[] = {108,112,124,132};

  TString runName = "submit_reco";
  TString runTag = "sr";

  TString pwdDir = TString(gSystem -> pwd()) + "/";
  TString subDir = pwdDir + runName + "/";
  TString outDir = subDir + "job_out/";
  TString logDir = subDir + "job_log/";
  TString allFull = pwdDir + runName + "_" + anaName + ".sh";
  TString nowFull = pwdDir + runName + "_" + anaName + ".now.sh";
  TString endFull = pwdDir + runName + "_" + anaName + ".end";

  TString recoDir = Form("/home/ejungwoo/data/reco/%s/",anaName.Data());
  gSystem -> mkdir(recoDir);
  std::ofstream fileReadme(recoDir+"Readme");
  fileReadme << readme << endl;

  std::ofstream submit_now(nowFull);
  std::ofstream submit_all(allFull);
  std::ofstream submit_end(endFull);
  submit_all << "set +x" << endl;

  // 108 112 124 132
  vector<int> runArray[4];
  vector<int> numEventsArray[4];
  int numEventsTotal[4] = {0};

  if (useDB)
  {
    ifstream fileRunDB(nameRunDB);

    std::string line0;
    std::getline(fileRunDB, line0);
    while (std::getline(fileRunDB, line0))
    {
      int run, sys, systar, numEvents, iSystem;
      stringstream((TString(line0).ReplaceAll(","," ")).Data()) >> run >> systar >> numEvents;
      //cout << run << " " << int(systar/1000) << " " << numEvents << endl;

      sys = int(systar/1000);

      if (sys==108) iSystem = 0;
      else if (sys==112) iSystem = 1;
      else if (sys==124) iSystem = 2;
      else if (sys==132) iSystem = 3;

      if (fSelectRun>0) {
        if (fSelectRun==run) {
          runArray[iSystem].push_back(run);
          numEventsArray[iSystem].push_back(numEvents);
          numEventsTotal[iSystem] += numEventsArray[iSystem].back();
        }
      }
      else {
        if (numEventsTotal[iSystem] < numEventsMax) {
          runArray[iSystem].push_back(run);
          numEventsArray[iSystem].push_back(numEvents);
          numEventsTotal[iSystem] += numEventsArray[iSystem].back();
        }
      }
    }
  }
  else
  {
    for (auto iSystem : iSystemsForAna)
    {
      auto sys = systems[iSystem];
      auto nameSelectRunEvent = Form(formSelectRunEvent ,sys);
      std::ifstream fileListSelRE(nameSelectRunEvent);

      TString fileNameSelRE;
      while (fileListSelRE >> fileNameSelRE) {
        std::ifstream fileSelRE(fileNameSelRE);

        int numEvents = 0;
        int run, event;
        while (fileSelRE >> run >> event)
          numEvents++;

        runArray[iSystem].push_back(run);
        numEventsArray[iSystem].push_back(numEvents);
      }
    }
  }

  for (auto iSystem : iSystemsForAna) {
    submit_now << "# " << systems[iSystem] << ": Number of Runs = " << runArray[iSystem].size() << ", Number of Events = " << numEventsTotal[iSystem] << endl;
    submit_all << "# " << systems[iSystem] << ": Number of Runs = " << runArray[iSystem].size() << ", Number of Events = " << numEventsTotal[iSystem] << endl;
    submit_end << "# " << systems[iSystem] << ": Number of Runs = " << runArray[iSystem].size() << ", Number of Events = " << numEventsTotal[iSystem] << endl;
  }

  int countSubmit = 0;

  for (auto iSystem : iSystemsForAna)
  {
    auto sys = systems[iSystem];
    TString recoSysDir = recoDir + Form("Sn%d/",sys);
    gSystem -> mkdir(recoSysDir);

    for (auto iRun=0; iRun<int(runArray[iSystem].size()); ++iRun) {
      auto run = runArray[iSystem].at(iRun);
      auto numEvents = numEventsArray[iSystem].at(iRun);

      auto numSplits = int(numEvents/numEventsInSplit)+1;

      for (auto split=0; split<numSplits; ++split)
      {
        TString subName = Form("%s%s_%d",runTag.Data(),anaName.Data(),countSubmit);
        //TString subName = Form("%s%d%s_%d_%d",runTag.Data(),sys,anaName.Data(),run,split);
        TString macFull = subDir + subName + ".sh";
        TString outFull = outDir + subName + ".out";
        TString logFull = logDir + subName + "_" + sys + "_" + run + "_" + split + ".log";
        //TString logFull = logDir + subName + ".log";
        TString tailLog1 = TString("$(tail -1 ") + logFull + " | head -1)";
        TString tailLog2 = TString("$(tail -2 ") + logFull + " | head -1)";
        TString tailLog3 = TString("$(tail -3 ") + logFull + " | head -1)";
        TString tailLog4 = TString("$(tail -4 ") + logFull + " | head -1)";

        submit_end << subName << " " << logFull << endl;
        cout << macFull << " " << logFull << endl;

        std::ofstream submit_macro(macFull);
        submit_all << "pjsub --bulk --sparam 0-0 " << macFull << " # " << logFull << endl;

        submit_macro << "#!/bin/bash" << endl;
        submit_macro << "#------ pjsub option -------- #" << endl;
        submit_macro << "#PJM -L rscunit=bwmpc" << endl;
        submit_macro << "#PJM -L rscgrp=batch" << endl;
        submit_macro << "#PJM -L vnode=1" << endl;
        submit_macro << "#PJM -L vnode-core=6" << endl;
        submit_macro << "#PJM -L vnode-mem=6Gi" << endl;
        submit_macro << "#PJM -g Q20393" << endl;
        submit_macro << "#PJM -j" << endl;
        submit_macro << "#PJM -o " << outFull << endl;
        submit_macro << "#------- Program execution ------- #" << endl;
        submit_macro << "export OMP_NUM_THREADS=1" << endl;
        submit_macro << "source /home/ejungwoo/environment.spiritroot.bwmpc.sh" << endl;
        submit_macro << "cd " << subDir << endl;
        TString commandRun = Form("root -q -b -l run_reco_experiment_auto.C\\(%d,%d,%d,%d,\\\"%s\\\"\\)", run, split, numEventsInSplit, sys, recoSysDir.Data());
        TString commandLog = Form(" > %s 2>&1", logFull.Data());
        submit_macro << commandRun << commandLog << endl;
        submit_now   << commandRun << endl;
        submit_macro << "echo '"<< sys << " ' " << tailLog4 << " >> " << endFull << endl;

        countSubmit++;
      }
    }
  }

  submit_end << endl;
  submit_end << "########################################" << endl;
  submit_end << endl;

  cout << endl;
  cout << allFull << endl;
  cout << nowFull << endl;
  cout << endFull << endl;
  cout << endl;

  gApplication -> Terminate();
}
