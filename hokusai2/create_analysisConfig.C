void create_analysisConfig()
{
  /*
  TString anaName = ".Kaneko";
  const int numParameters = 6;
  TString parameters[numParameters][5] = { // 108, 132, 112, 124
    {"DataDir"          , "/home/ejungwoo/data/trim/NewAna.2070.0231288/Sn108/"
                        , "/home/ejungwoo/data/trim/NewAna.2070.0231288/Sn132/"
                        , "/home/ejungwoo/data/trim/NewAna.2070.0231288/Sn112/"
                        , "/home/ejungwoo/data/trim/NewAna.2070.0231288/Sn124/"},
    {"MultiplicityMin"  , "55" ,"56","100","100"},
    {"MultiplicityMax"  , "100"          ,"","",""},
    {"MultiplicityDPOCA", "20"           ,"","",""},
    {"NClus"            , "15"           ,"","",""},
    {"Phi"              , "0-20,330-360" ,"","",""},
  };
  */

  TString anaName = ".Tommy";
  const int numParameters = 6;
  TString parameters[numParameters][5] = { // 108, 132, 112, 124
    {"DataDir"          , "/home/ejungwoo/data/trim/NewAna.2070.0231288/Sn108/"
                        , "/home/ejungwoo/data/trim/NewAna.2070.0231288/Sn132/"
                        , "/home/ejungwoo/data/trim/NewAna.2070.0231288/Sn112/"
                        , "/home/ejungwoo/data/trim/NewAna.2070.0231288/Sn124/"},
    {"MultiplicityMin"  , "50"   ,"","",""},
    {"MultiplicityMax"  , "100"          ,"","",""},
    {"MultiplicityDPOCA", "20"           ,"","",""},
    {"NClus"            , "15"           ,"","",""},
    {"Phi"              , "160-220" ,"","",""},
  };



  /////////////////////
  /////////////////////
  /////////////////////


  int systems[] = {108, 132, 112, 124};
  ofstream newConfFile[4];
  for (auto isys : {0,1,2,3}) {
    TString newConfName = Form("analysisInputFiles/analysisConfig/analysisSn%dCM%s.xml",systems[isys],anaName.Data());
    cout << newConfName << endl;
    newConfFile[isys].open(newConfName);
  }


  for (auto isys : {0,1,2,3})
  {
    std::string ssline;
    ifstream dummy(Form("analysisInputFiles/analysisConfig/analysisSn%dCM.dummy.xml",systems[isys]));
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
