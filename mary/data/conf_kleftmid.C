void conf_kleftmid()
{
  fSTVersion     = "kleft.NewAna.2100.de29163";
  fVOutShort     = "kleftmid";
  fTrackMultHL   = "$$(multll)";
  fCut0String    = "$$(cut_pe)";
  fPhiLL         = 0;
  fPhiHL         = 20;
  fPhiLL2        = 330;
  fPhiHL2        = 360;
  fTtaLL         = 70;
  fTtaHL         = 110;
  fProbLL        = 0.2;
  fEffLL         = 0.05;
  fSDHLString    = "$$(asdcut)";
  fPtoaLL        = 0;
  fPtoaHL        = 10000;
  fNyLL          = -1;

  fSDHLCut[0] = 2.2;
  fSDHLCut[1] = 2.0;
  fSDHLCut[2] = 1.8;
  fSDHLCut[3] = 2.0;
  fSDHLCut[4] = 2.0;
  fSDHLCut[5] = 2.0;

//  fSDHLCut[0] = 3.0;
//  fSDHLCut[1] = 2.5;
//  fSDHLCut[2] = 2.0;
//  fSDHLCut[3] = 2.0;
//  fSDHLCut[4] = 2.0;
//  fSDHLCut[5] = 2.0;

  fvar_foby_cm.setnmm(20,-0.5,0.5);
}
