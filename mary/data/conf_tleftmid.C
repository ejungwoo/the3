void conf_tleftmid()
{
  fSTVersion     = "tleft.NewAna.2100.de29163";
  fVOutShort     = "tleftmid";
  fTrackMultHL   = "50";
  fCut0String    = "$$(cut_pe)";
  fPhiLL         = 0;
  fPhiHL         = 20;
  fPhiLL2        = 330;
  fPhiHL2        = 360;
  fTtaLL         = 70;
  fTtaHL         = 110;
  fProbLL        = 0.2;
  fEffLL         = 0.05;
  fSDHLString    = "3";
  fPtoaLL        = 0;
  fPtoaHL        = 10000;
  fNyLL          = -1;

  fvar_foby_cm.setnmm(20,-0.5,0.5);
}
