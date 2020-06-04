void print_mass()
{
  const Double_t kAu2Gev = 0.9314943228;
  const Double_t khSlash = 1.0545726663e-27;
  const Double_t kErg2Gev = 1/1.6021773349e-3;
  const Double_t khShGev = khSlash*kErg2Gev;
  const Double_t kYear2Sec = 3600*24*365.25;

  TDatabasePDG *db = TDatabasePDG::Instance();
  db -> AddParticle("Deuteron","Deuteron",2*kAu2Gev+8.071e-3,kTRUE, 0,3,"Ion",1000010020);
  db -> AddParticle("Triton","Triton",3*kAu2Gev+14.931e-3,kFALSE, khShGev/(12.33*kYear2Sec),3,"Ion",1000010030);
  db -> AddParticle("Alpha","Alpha",4*kAu2Gev+2.424e-3,kTRUE, khShGev/(12.33*kYear2Sec),6,"Ion",1000020040);
  db -> AddParticle("HE3","HE3",3*kAu2Gev+14.931e-3,kFALSE, 0,6,"Ion",1000020030);

  for (auto pdg : {2212, 1000010020, 1000010030, 1000020030, 1000020040})
  {
    auto particle = TDatabasePDG::Instance() -> GetParticle(pdg);
    auto particleMass = particle -> Mass()*1000; // GeV to MeV
    cout << pdg << " " << particleMass << endl;
  }
}
