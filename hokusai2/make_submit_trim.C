void make_submit_trim()
{
  const char *anaName = "fix6";
  int systems[] = {108,132};
  TString runName = "submit_trim";
  TString runTag = "st";
  int numGroup = 10;

  gSystem -> Exec(Form("mkdir -p %s",runName.Data()));

  TString pwdDir = TString(gSystem -> pwd()) + "/";
  TString subDir = pwdDir + runName + "/";
  TString outDir = subDir + "job_out/";
  TString logDir = subDir + "job_log/";
  TString nowFull = pwdDir + runName + "_"+anaName + ".now.sh";
  TString allFull = pwdDir + runName + "_"+anaName + ".sh";
  TString endFull = pwdDir + runName + "_"+anaName + ".end";
  TString nameRecoEnd = Form("submit_reco_%s.end",anaName);

  std::ofstream submit_now(nowFull);
  std::ofstream submit_all(allFull);
  std::ofstream submit_end(endFull);
  submit_all << "set +x" << endl;

  int countSubmit = 0;
  for (auto system : systems)
  {
    TString nameListConc = Form("analysisInputFiles/listFiles/list_reco_%s_%d.txt",anaName,system);
    TString nameListRS = nameListConc;
    nameListRS.ReplaceAll("list_reco","list_run_split");

    ifstream fileRecoEnd(nameRecoEnd);
    ofstream fileListConc(nameListConc);
    int countFiles = 0;
    std::string line0;
    while (std::getline(fileRecoEnd, line0))
    {
      TString line(line0);
      if (line.Index(".conc.root")>0) {
        auto index0 = line.Index("/home/");
        if (line.Index((TString("Sn")+system).Data())>=0) {
          cout << line0 << endl;
          fileListConc << line(index0,line.Sizeof()) << endl;
          countFiles++;
        }
      }
    }
    fileListConc.close();

    ifstream fileListConc2(nameListConc);
    ofstream fileListRS(nameListRS);
    TObjArray array;
    TString recoName, pathName;
    while(fileListConc2 >> recoName) {
      auto indexAtRun = recoName.Index("/run");
      if (pathName.IsNull()) {
        pathName = recoName(0,indexAtRun+1);
      }
      TString runNo = recoName(indexAtRun+4,4);
      TString sptName = recoName(indexAtRun+10,3);
      if (!sptName.IsOct()) sptName = recoName(indexAtRun+10,2);
      if (!sptName.IsOct()) sptName = recoName(indexAtRun+10,1);
      auto named = (TNamed *) array.FindObject(runNo);
      if (named!=nullptr) {
        TString title = named -> GetTitle();
        named -> SetTitle(title+" "+sptName);
      }
      else
        array.Add(new TNamed(runNo,sptName));
    }
    fileListRS << pathName << endl;
    for (auto i=0; i<array.GetEntries(); ++i) {
      auto named = (TNamed *) array.At(i);
      TString title = named -> GetTitle();
      auto title2 = title;
      title2.ReplaceAll(" ","  ");
      auto num = title2.Sizeof() - title.Sizeof() + 1;
      fileListRS << named -> GetName() << " " << Form("-%d", num) << " " << title <<endl;
    }
    fileListConc2.close();
    fileListRS.close();

    cout << nameRecoEnd << " " << nameListConc << " " << nameListRS << " #(" << countFiles << ")" << endl;
    cout << endl;

    ifstream fileListRS2(nameListRS);
    int run, numSplits;
    TString pathToReco;
    fileListRS2 >> pathToReco;

    TString pathToTrim = pathToReco;
    pathToTrim.ReplaceAll("reco","trim");
    gSystem -> Exec(TString("mkdir -p ")+pathToTrim);

    while (fileListRS2 >> run >> numSplits)
    {
      TString subName = Form("%s%s_%d",runTag.Data(),anaName,countSubmit);
      TString macFull = subDir + subName + ".sh";
      TString outFull = outDir + subName + ".out";
      TString logFull = logDir + subName + "_" + system + "_" + run + "_$PJM_BULKNUM.log";
      TString tailLog1 = TString("$(tail -1 ") + logFull + " | head -1)";
      TString tailLog2 = TString("$(tail -2 ") + logFull + " | head -1)";
      TString tailLog3 = TString("$(tail -3 ") + logFull + " | head -1)";
      TString tailLog4 = TString("$(tail -4 ") + logFull + " | head -1)";

      int numBulk;

      if (numSplits<0)
      {
        numGroup = 1;
        numSplits = -numSplits;
        int splitNo;
        numBulk = 0;
        for (auto iSplit=0; iSplit<numSplits; ++iSplit)
        {
          fileListRS2 >> splitNo;

          subName = Form("%s%s_%d",runTag.Data(),anaName,countSubmit);
          macFull = subDir + subName + ".sh";
          outFull = outDir + subName + ".out";
          logFull = logDir + subName + "_" + system + "_" + run + ".log";
          tailLog1 = TString("$(tail -1 ") + logFull + " | head -1)";
          tailLog2 = TString("$(tail -2 ") + logFull + " | head -1)";
          tailLog3 = TString("$(tail -3 ") + logFull + " | head -1)";
          tailLog4 = TString("$(tail -4 ") + logFull + " | head -1)";

          submit_end << subName << " " << logFull << endl;
          cout << macFull << " " << logFull << endl;

          std::ofstream submit_macro(macFull);
          submit_all << "pjsub --bulk --sparam 0-" << numBulk << " " << macFull << " # " << logFull << endl;

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
          TString commandRun = Form("root -q -b -l run_trim_data.C\\(%d,%d,%d,\\\"%s\\\",\\\"%s\\\",\\\"%s\\\"\\)",run,10000+splitNo,numGroup,anaName,pathToReco.Data(),pathToTrim.Data());
          TString commandLog = Form(" > %s 2>&1", logFull.Data());
          submit_macro << commandRun << commandLog << endl;
          submit_macro << "echo '"<< run << " ' " << tailLog3 << " >> " << endFull << endl;

          submit_now << commandRun << endl;

          ++countSubmit;
        }

        continue; // continue run scope
      }

      ++countSubmit;

      submit_end << subName << " " << logFull << endl;
      cout << macFull << " " << logFull << endl;

      if (numSplits<=numGroup)
        numBulk = 0;
      else {
        numBulk = std::ceil(double(numSplits)/numGroup);
      }

      std::ofstream submit_macro(macFull);
      submit_all << "pjsub --bulk --sparam 0-" << numBulk << " " << macFull << endl;

      submit_macro << "#!/bin/bash" << endl;
      submit_macro << "#------ pjsub option -------- #" << endl;
      submit_macro << "#PJM -L rscunit=bwmpc" << endl;
      submit_macro << "#PJM -L rscgrp=batch" << endl;
      submit_macro << "#PJM -L vnode=1" << endl;
      submit_macro << "#PJM -L vnode-core=10" << endl;
      submit_macro << "#PJM -L vnode-mem=16Gi" << endl;
      submit_macro << "#PJM -g Q20393" << endl;
      submit_macro << "#PJM -j" << endl;
      submit_macro << "#PJM -o " << outFull << endl;
      submit_macro << "#------- Program execution ------- #" << endl;
      submit_macro << "export OMP_NUM_THREADS=16" << endl;
      submit_macro << "source /home/ejungwoo/environment.spiritroot.bwmpc.sh" << endl;
      submit_macro << "cd " << subDir << endl;
      TString commandRun = Form("root -q -b -l run_trim_data.C\\(%d,$PJM_BULKNUM,%d,\\\"%s\\\",\\\"%s\\\",\\\"%s\\\"\\)",run,numGroup,anaName,pathToReco.Data(),pathToTrim.Data());
      TString commandLog = Form(" > %s 2>&1", logFull.Data());
      submit_macro << commandRun << commandLog << endl;
      submit_macro << "echo '"<< run << " ' " << tailLog3 << " >> " << endFull << endl;

      submit_now << commandRun << endl;
    }
  }

  submit_end << endl;
  submit_end << "########################################" << endl;
  submit_end << endl;

  cout << endl;
  cout << nowFull << endl;
  cout << allFull << endl;
  cout << endFull << endl;
  cout << endl;

  gApplication -> Terminate();
}
