void make_submit_macros()
{
  std::ofstream submit_all("submit_all.sh");
  submit_all << "set +x" << endl;

  for (auto prob_cut : {0.1,0.8,0.9})
    for (auto version : {0,1})
      for (auto sys : {108,132})
        for (auto idx=0; idx<20; ++idx)
        {
          TString macName = Form("m_%d_%d_%d_%d.sh",sys,idx,version,int(10*prob_cut));
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
          submit_macro << "root draw_from_ana.C\\(" << sys << "," << idx << "," << version << "," << prob_cut << "\\) -b -q -l > log/log_" << macName << " 2>&1" << endl;
        }
}
