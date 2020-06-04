void make_exe()
{
  auto numSplits = 20;

  ofstream file132("exe_132.sh"); for (auto iSplit=0; iSplit<numSplits; ++iSplit) file132 << Form("root -q -b -l run_analysis_xml.C\\(132,%d,%d\\) > log/log_132_%d 2>&1",numSplits,iSplit,iSplit) << endl;
  ofstream file108("exe_108.sh"); for (auto iSplit=0; iSplit<numSplits; ++iSplit) file108 << Form("root -q -b -l run_analysis_xml.C\\(108,%d,%d\\) > log/log_108_%d 2>&1 ",numSplits,iSplit,iSplit) << endl;

  for (auto version : {0,1}) {
    std::ofstream submit_all(Form("submit_all_v%d.sh",version));
    submit_all << "set +x" << endl;

    //for (auto sys : {108,132})
    for (auto sys : {132})
      for (auto idx=0; idx<numSplits; ++idx)
      {
        TString macName = Form("mac_%d_%d_v%d.sh",sys,idx,version);
        std::ofstream submit_macro(macName);
        submit_all << "pjsub --bulk --sparam 0-0 " << macName << endl;

        submit_macro << "#!/bin/bash" << endl;
        submit_macro << "#------ pjsub option -------- #" << endl;
        submit_macro << "#PJM -L rscunit=bwmpc" << endl;
        submit_macro << "#PJM -L rscgrp=batch" << endl;
        submit_macro << "#PJM -L vnode=1" << endl;
        submit_macro << "#PJM -L vnode-core=1" << endl;
        submit_macro << "#PJM -L vnode-mem=1Gi" << endl;
        submit_macro << "#PJM -g Q20393" << endl;
        submit_macro << "#PJM -j" << endl;
        submit_macro << "#------- Program execution ------- #" << endl;
        submit_macro << "export OMP_NUM_THREADS=1" << endl;
        submit_macro << "source /home/ejungwoo/environment.spiritroot.bwmpc.sh" << endl;
        submit_macro << "cd /home/ejungwoo/the3/hokusai/submit_xml/" << endl;
        submit_macro << Form("root -q -b -l run_analysis_xml.C\\(%d,%d,%d,%d\\) > log/log_%s 2>&1 ",sys,numSplits,idx,version,macName.Data()) << endl;
      }
  }
}
