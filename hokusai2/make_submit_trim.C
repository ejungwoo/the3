void make_submit_trim()
{
  int recoDate = 20191214;

  int numGroup = 10;

  TString runName = "submit_trim";
  TString subTag = "st";

  TString pwdDir = TString(gSystem -> pwd()) + "/";
  TString subDir = pwdDir + runName + "/";
  TString outDir = subDir + "job_out/";
  TString logDir = subDir + "job_log/";
  TString allFull = pwdDir + runName+".sh";
  TString endFull = pwdDir + runName+".end";

  gSystem -> Exec(Form("mkdir -p %s",runName.Data()));

  std::ofstream submit_all(allFull);
  submit_all << "set +x" << endl;
  cout << allFull << endl;
  cout << endl;

  std::ofstream submit_end(endFull);

  for (auto sys : {108,132,112,124})
  {
    ifstream list_run_split(Form("analysisInputFiles/listFiles/list_run_split_%d_%d.txt",recoDate,sys));
    int run, numSplits;
    TString pathToData;
    list_run_split >> pathToData;

    while (list_run_split >> run >> numSplits)
    {
      TString subName = Form("%s%d%d",subTag.Data(),run,sys);
      TString macFull = subDir + subName + ".sh";
      TString outFull = outDir + subName + ".out";
      TString logFull = logDir + subName + "_$PJM_BULKNUM.log";
      TString tailLog1 = TString("$(tail -1 ") + logFull + " | head -1)";
      TString tailLog2 = TString("$(tail -2 ") + logFull + " | head -1)";
      TString tailLog3 = TString("$(tail -3 ") + logFull + " | head -1)";
      TString tailLog4 = TString("$(tail -4 ") + logFull + " | head -1)";

      submit_end << subName << " " << logFull << endl;
      cout << macFull << " " << logFull << endl;

      int numBulk;
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
      submit_macro << Form("root -q -b -l run_trim_data.C\\(%d,$PJM_BULKNUM,%d,%d,\\\"%s\\\"\\) > %s 2>&1 ",run,numGroup,recoDate,pathToData.Data(),logFull.Data()) << endl;
      submit_macro << "echo '"<< run << " ' " << tailLog3 << " >> " << endFull << endl;
    }
  }

  submit_end << endl;
  submit_end << "====" << endl;
  submit_end << endl;

  cout << endl;
  cout << endFull << endl;

  gApplication -> Terminate();
}
