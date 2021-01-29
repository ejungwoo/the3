void make_submit_summary()
{
  //TString anaName = "fix6";
  TString anaName = "fix5";
  const char *trimPath = "/home/ejungwoo/data/trim";
  int systemSummary[] = {108,132};

  TString runName = "submit_summary";
  TString runTag = "ss";
  TString pwdDir = TString(gSystem -> pwd()) + "/";
  TString subDir = pwdDir + runName + "/";
  TString outDir = subDir + "job_out/";
  TString logDir = subDir + "job_log/";
  TString allFull = pwdDir + runName + "_" + anaName + ".sh";
  TString endFull = pwdDir + runName + "_" + anaName + ".end";
  TString nowFull = pwdDir + runName + "_" + anaName + ".now.sh";


  for (int idxRightLeft : {0,1}) {
    for (int iMultRange : {0,1})
    {
      TString confName, phiRange;
           if (idxRightLeft==0) { confName = anaName + "_right";  phiRange = "160-220"; }
      else if (idxRightLeft==1) { confName = anaName + "_left";   phiRange = "0-20,320-360"; }

      int min = 45;
      int max = 54;
           if (iMultRange==0) { min = 45; max =  54; }
      else if (iMultRange==1) { min = 55; max = 100; }
      confName = confName + "_" + min + "_" + max;

      const int numParameters = 6;
      TString parameters[numParameters][5] = { // 108, 132, 112, 124
        {"DataDir"          , Form("%s/%s/Sn108/",trimPath,anaName.Data())
                            , Form("%s/%s/Sn132/",trimPath,anaName.Data())
                            , Form("%s/%s/Sn112/",trimPath,anaName.Data())
                            , Form("%s/%s/Sn124/",trimPath,anaName.Data())},
        {"MultiplicityMin"  , TString()+(min-1), "","",""},
        {"MultiplicityMax"  , TString()+(max+1), "","",""},
        {"MultiplicityDPOCA", "20",  "","",""},
        {"NClus"            , "15",  "","",""},
        {"Phi"              , phiRange ,"","",""},
      };

      int systemsAll[] = {108, 132, 112, 124};

      ofstream newConfFile[4];
      for (auto isys : {0,1,2,3}) {
        TString newConfName = Form("analysisInputFiles/analysisConfig/analysis.%s.Sn%dCM.xml",confName.Data(),systemsAll[isys]);
        TString dummyConfName = Form("analysisInputFiles/analysisConfig/analysis.dummy.Sn%dCM.xml",systemsAll[isys]);
        cout << dummyConfName << " " << newConfName << endl;
        newConfFile[isys].open(newConfName);
      }

      for (auto isys : {0,1,2,3})
      {
        std::string ssline;
        TString dummyConfName = Form("analysisInputFiles/analysisConfig/analysis.dummy.Sn%dCM.xml",systemsAll[isys]);
        ifstream dummy(dummyConfName);

        while (std::getline(dummy, ssline)) {
          auto line = TString(ssline);

          for (auto ipar=0; ipar<numParameters; ++ipar) {
            TString parName = parameters[ipar][0];
            TString parValue = parameters[ipar][isys+1];
            if (parValue.IsNull())
              parValue = parameters[ipar][1];

            if (line.Index(parName.Data())>=0) {
              auto nnn = line.First('<');
              TString spacings; for (auto i=0; i<nnn; ++i) spacings += " ";
              line = Form("%s<%s>%s</%s>",spacings.Data(),parName.Data(),parValue.Data(),parName.Data());
            }
          }

          newConfFile[isys] << line << endl;
        }
      }
    }
  }

  std::ofstream submit_all(allFull);
  std::ofstream submit_end(endFull);
  std::ofstream submit_now(nowFull);

  submit_all << "set +x" << endl;

  int countSubmit = 0;
  for (TString anaName0 : { "right_55_100", "right_45_54", "left_55_100", "left_45_54" })
  {
    TString anaName1 = Form("%s_%s",anaName.Data(),anaName0.Data());

    for (auto system : systemSummary)
    {
      //TString subName = Form("%s%d%s",runTag.Data(),system,anaName1.Data());
      TString subName = Form("%s%s_%d",runTag.Data(),anaName.Data(),countSubmit);
      TString macFull = subDir + subName + ".sh";
      TString outFull = outDir + subName + ".out";
      //TString logFull = logDir + subName + ".log";
      TString logFull = logDir + subName + "_" + system + "_" + anaName0 + ".log";
      TString tailLog1 = TString("$(tail -1 ") + logFull + " | head -1)";
      TString tailLog2 = TString("$(tail -2 ") + logFull + " | head -1)";
      TString tailLog3 = TString("$(tail -3 ") + logFull + " | head -1)";
      TString tailLog4 = TString("$(tail -4 ") + logFull + " | head -1)";

      submit_end << subName << " " << logFull << endl;
      cout << macFull << " " << logFull << endl;

      std::ofstream submit_macro(macFull);
      submit_all << "pjsub --bulk --sparam 0-0 " << macFull << "  # " << logFull << endl;

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
      TString commandRun = Form("root -q -b -l run_analysis_xml.C\\(%d,\\\"%s\\\",\\\"%s\\\"\\)",system,anaName.Data(),anaName1.Data());
      TString commandLog = Form(" > %s 2>&1",logFull.Data());
      submit_macro << commandRun << commandLog << endl;
      submit_macro << "echo '"<< system << " ' " << tailLog4 << " >> " << endFull << endl;
      submit_macro << "echo '"<< system << " ' " << tailLog3 << " >> " << endFull << endl;
      submit_macro << "echo '"<< system << " ' " << tailLog2 << " >> " << endFull << endl;
      submit_macro << "echo '"<< system << " ' " << tailLog1 << " >> " << endFull << endl;

      submit_now << commandRun << endl;

      ++countSubmit;
    }

    submit_end << endl;
    submit_end << "====" << endl;
    submit_end << endl;
  }

  cout << endl;
  cout << allFull << endl;
  cout << nowFull << endl;
  cout << endFull << endl;
  cout << endl;

  gApplication -> Terminate();
}
