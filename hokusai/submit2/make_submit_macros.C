void make_submit_macros()
{
  int runs108[] = {2272, 2273, 2274, 2275, 2276, 2283, 2284, 2285, 2286, 2288, 2289, 2291, 2310, 2311, 2314, 2315, 2320, 2322, 2323, 2324, 2325, 2331, 2332, 2333, 2334, 2335, 2336, 2337, 2340, 2341, 2362, 2363, 2368, 2369, 2370, 2371, 2372, 2373, 2374, 2375, 2378, 2379, 2380, 2381, 2382, 2383, 2384, 2385, 2386, 2387, 2388, 2389, 2391, 2392, 2393, 2394, 2395, 2396, 2397, 2398, 2399, 2400, 2401, 2402, 2429, 2432, 2433, 2434, 2437, 2438, 2439, 2440, 2442, 2453, 2461, 2462, 2463, 2501, 2502, 2503, 2505, 2506, 2507, 2508, 2509};

  int runs132[] = {2841, 2843, 2844, 2845, 2846, 2848, 2849, 2850, 2851, 2852, 2855, 2856, 2857, 2858, 2859, 2860, 2861, 2875, 2877, 2878, 2879, 2880, 2881, 2882, 2883, 2884, 2887, 2888, 2889, 2890, 2891, 2892, 2893, 2894, 2896, 2898, 2899, 2900, 2901, 2902, 2903, 2904, 2905, 2907, 2914, 2916, 2917, 2919, 2920, 2921, 2922, 2924, 2925, 2926, 2927, 2929, 2930, 2931, 2932, 2933, 2934, 2935, 2936, 2939, 2940, 2941, 2942, 2943, 2944, 2945, 2946, 2948, 2955, 2956, 2958, 2959, 2960, 2961, 2962, 2964, 2965, 2966, 2968, 2969, 2970, 2971, 2972, 2973, 2975, 2976, 2977, 2978, 2979, 2980, 2981, 2982, 2983, 2984, 2985, 2986, 2988, 2989, 2990, 2991, 2992, 2993, 2997, 2999, 3000, 3002, 3003, 3007, 3039};

  {
    int sys1 = 108;
    int sys2 = 112;
    TString strSys1 = Form("%d",sys1);
    TString strSys2 = Form("%d",sys2);

    std::ofstream submit_all_108("submit_all_108.sh");
    submit_all_108 << "set +x" << endl;

    for (auto run : runs108)
    {
      int run1 = run;
      int run2 = run;
      TString strRun1 = Form("%d",run1);
      TString strRun2 = Form("%d",run2);

      TString macName = Form("mac_%d.sh",run);
      std::ofstream submit_macro(macName);
      submit_all_108 << "pjsub --bulk --sparam 0-0 " << macName << endl;

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
      submit_macro << "cd /home/ejungwoo/reco/submit_pid/" << endl;
      submit_macro << "root run_analysis.C\\("+strRun1+","+strRun2+","+strSys1+","+strSys2+"\\) -b -q -l > log/log_run"+strRun1+" 2>&1" << endl;
    }

    TString exeName = "exe_108.sh";
    std::ofstream execute_macro(exeName);
    auto numRuns = 0;
    auto run1 = 0;
    auto run2 = 0;

    for (auto run : runs108) {
      numRuns++;
      if (numRuns == 1)
        run1 = run;
      else if (numRuns == 5) {
        run2 = run;

        TString strRun1 = Form("%d",run1);
        TString strRun2 = Form("%d",run2);

        numRuns = 0;
        execute_macro << "root run_analysis.C\\("+strRun1+","+strRun2+","+strSys1+","+strSys2+"\\) -b -q -l > log/log_run"+strRun1+" 2>&1" << endl;
      }

      run2 = run;
    }

    TString strRun1 = Form("%d",run1);
    TString strRun2 = Form("%d",run2);

    execute_macro << "root run_analysis.C\\("+strRun1+","+strRun2+","+strSys1+","+strSys2+"\\) -b -q -l > log/log_run"+strRun1+" 2>&1" << endl;
  }

  {
    int sys1 = 132;
    int sys2 = 124;
    TString strSys1 = Form("%d",sys1);
    TString strSys2 = Form("%d",sys2);

    std::ofstream submit_all_132("submit_all_132.sh");
    submit_all_132 << "set +x" << endl;

    for (auto run : runs132)
    {
      int run1 = run;
      int run2 = run;
      TString strRun1 = Form("%d",run1);
      TString strRun2 = Form("%d",run2);

      TString macName = Form("mac_%d.sh",run);
      std::ofstream submit_macro(macName);
      submit_all_132 << "pjsub --bulk --sparam 0-0 " << macName << endl;

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
      submit_macro << "cd /home/ejungwoo/reco/submit_pid/" << endl;
      submit_macro << "root run_analysis.C\\("+strRun1+","+strRun2+","+strSys1+","+strSys2+"\\) -b -q -l > log/log_run"+strRun1+" 2>&1" << endl;
    }

    TString exeName = "exe_132.sh";
    std::ofstream execute_macro(exeName);
    auto numRuns = 0;
    auto run1 = 0;
    auto run2 = 0;

    for (auto run : runs132) {
      numRuns++;
      if (numRuns == 1)
        run1 = run;
      else if (numRuns == 5) {
        run2 = run;

        TString strRun1 = Form("%d",run1);
        TString strRun2 = Form("%d",run2);

        numRuns = 0;
        execute_macro << "root run_analysis.C\\("+strRun1+","+strRun2+","+strSys1+","+strSys2+"\\) -b -q -l > log/log_run"+strRun1+" 2>&1" << endl;
      }

      run2 = run;
    }

    TString strRun1 = Form("%d",run1);
    TString strRun2 = Form("%d",run2);

    execute_macro << "root run_analysis.C\\("+strRun1+","+strRun2+","+strSys1+","+strSys2+"\\) -b -q -l > log/log_run"+strRun1+" 2>&1" << endl;

  }
}
