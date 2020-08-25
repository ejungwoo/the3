void make_submit_summary()
{
  //TString anaName = "";
  TString anaName = ".Tommy";

  TString runName = "submit_summary";
  TString runTag = "ss";

  TString pwdDir = TString(gSystem -> pwd()) + "/";
  TString subDir = pwdDir + runName + "/";
  TString outDir = subDir + "job_out/";
  TString logDir = subDir + "job_log/";
  TString allFull = pwdDir + runName+".sh";
  TString endFull = pwdDir + runName+".end";

  std::ofstream submit_all(allFull);
  submit_all << "set +x" << endl;
  cout << allFull << endl;
  cout << endl;

  std::ofstream submit_end(endFull);

  for (auto sys : {108,112,124,132})
  //for (auto sys : {108,132})
  {
    TString subName = Form("%s%d",runTag.Data(),sys);
    TString macFull = subDir + subName + ".sh";
    TString outFull = outDir + subName + ".out";
    TString logFull = logDir + subName + ".log";
    TString tailLog1 = TString("$(tail -1 ") + logFull + " | head -1)";
    TString tailLog2 = TString("$(tail -2 ") + logFull + " | head -1)";
    TString tailLog3 = TString("$(tail -3 ") + logFull + " | head -1)";
    TString tailLog4 = TString("$(tail -4 ") + logFull + " | head -1)";

    submit_end << subName << " " << logFull << endl;
    cout << macFull << " " << logFull << endl;

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
    submit_macro << "root -q -b -l run_analysis_xml.C\\(" << sys << ",\\\"" << anaName << "\\\"" << "\\)" << " > " << logFull.Data() << " 2>&1 " << endl;
    submit_macro << "echo '"<< sys << " ' " << tailLog4 << " >> " << endFull << endl;
    submit_macro << "echo '"<< sys << " ' " << tailLog3 << " >> " << endFull << endl;
    submit_macro << "echo '"<< sys << " ' " << tailLog2 << " >> " << endFull << endl;
    submit_macro << "echo '"<< sys << " ' " << tailLog1 << " >> " << endFull << endl;
  }

  submit_end << endl;
  submit_end << "====" << endl;
  submit_end << endl;

  cout << endl;
  cout << endFull << endl;

  gApplication -> Terminate();
}
