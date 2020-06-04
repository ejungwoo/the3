void make_submit_macros()
{
  int version = 1;

  std::ofstream submit_all(Form("submit_all_v%d.sh",version));
  submit_all << "set +x" << endl;

  for (auto sys : {108,132})
    for (auto idx=0; idx<20; ++idx)
    {
      TString macName = Form("mac_%d_%d_%d.sh",sys,idx,version);
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
      submit_macro << "cd /home/ejungwoo/the3/hokusai/submit_ana/" << endl;
      submit_macro << "root draw_from_ana.C\\(" << sys << "," << idx << "," << version << "\\) -b -q -l > log/log_" << macName << " 2>&1" << endl;
    }
}
