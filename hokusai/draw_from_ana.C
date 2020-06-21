#include "/home/ejungwoo/config/ejungwoo.h"

void draw_from_ana(
  int system=132,
  int iSplit=1,
  int version=1,
  double probCut = 0.9,
  double effCut = 0.1
  )
{
  bool testmode = true;

  ejungwoo::gstat("ien");
  ejungwoo::gfill(0);

  const char *cut_string = Form("e%d_p%d",int(10*effCut),int(10*probCut));

  TString file_name_in = Form("/home/ejungwoo/data/ana/NewAna.2034.45b9400/Sn%d/Sn%d_%d_ana_v%d.NewAna.2034.45b9400.root",system,system,iSplit,version);
  TString file_name_out = Form("/home/ejungwoo/the3/hokusai/submit_ana/data_with_trim/summary_hist_%d_s%d_v%d_%s.root",system,iSplit,version,cut_string);
  if (testmode)
    file_name_out = Form("/home/ejungwoo/the3/hokusai/submit_ana/data_with_trim/test_summary_hist_%d_s%d_v%d_%s.root",system,iSplit,version,cut_string);

  auto yBeamCM = 0.36;
  ejungwoo::setpar("ybeam",yBeamCM);
  ejungwoo::setpar("rd","TMath::RadToDeg()");

  TString particle_name[] = {"p","d","t","he3","he4"};
  double particle_mass[] = {938.272 ,1871.06 ,2809.41 ,2809.41 ,3728.4};
  double particle_a[] = {1,2,3,3,4,};
  double particle_z[] = {1,1,1,2,2,};

  ejungwoo::variable var_keCM   ("keCM"   ,"(sqrt(CMVector[$$(i)].fElements.Mag2()+$$(m2))-$$(m))"        ,"$$(WCUT)" ,"KE_{CM} (MeV)  "        ,200, 0, 500);
  ejungwoo::variable var_keoaCM ("keoaCM" ,"(sqrt(CMVector[$$(i)].fElements.Mag2()+$$(m2))-$$(m))/$$(a)"  ,"$$(WCUT)" ,"KE_{CM}/A (MeV)"        ,200, 0, 500);
  ejungwoo::variable var_keozCM ("keozCM" ,"(sqrt(CMVector[$$(i)].fElements.Mag2()+$$(m2))-$$(m))/$$(z)"  ,"$$(WCUT)" ,"KE_{CM}/Z (MeV)"        ,200, 0, 500);
  ejungwoo::variable var_ptoaCM ("ptoaCM" ,"CMVector[$$(i)].fElements.Perp()/$$(a)"                       ,"$$(WCUT)" ,"p_{T}/A (MeV/c)"        ,200, 0, 1500);
  ejungwoo::variable var_yyCM   ("yyCM"   ,"FragRapidity[$$(i)].fElements/$$(ybeam)"                      ,"$$(WCUT)" ,"y_{CM}/y_{beam,CM}"     ,200, -1.5, 3.0);
  ejungwoo::variable var_ttaCM  ("ttaCM"  ,"CMVector[$$(i)].fElements.Theta()*$$(rd)"                     ,"$$(WCUT)" ,"#theta_{CM} (Deg.)"     ,200, 0, 180);
  ejungwoo::variable var_phiCM  ("phiCM"  ,"CMVector[$$(i)].fElements.Phi()*$$(rd)"                       ,"$$(WCUT)" ,"#phi_{CM} (Deg.)"       ,200, -180, 180);
  ejungwoo::variable var_ttaLab ("ttaLab" ,"LabVector[$$(i)].fElements.Theta()*$$(rd)"                    ,"$$(WCUT)" ,"#theta_{Lab} (Deg.)"    ,200, 0, 180);
  ejungwoo::variable var_phiLab ("phiLab" ,"LabVector[$$(i)].fElements.Phi()*$$(rd)"                      ,"$$(WCUT)" ,"#phi_{Lab} (Deg.)"      ,200, -180, 180);
  ejungwoo::variable var_eff    ("eff"    ,"Eff[$$(i)].fElements"                                         ,"$$(WCUT)" ,"efficiency"             ,200, 0, 1);
  ejungwoo::variable var_prob   ("prob"   ,"Prob[$$(i)].fElements"                                        ,"$$(WCUT)" ,"probability"            ,200, 0, 1);
  ejungwoo::variable var_pLab   ("pLab"   ,"LabVector[$$(i)].fElements.Mag()"                             ,"$$(WCUT)" ,"p_{Lab} (MeV/c^{2})"    ,400, 0, 3000);
  ejungwoo::variable var_pozLab ("pozLab" ,"LabVector[$$(i)].fElements.Mag()/$$(z)"                       ,"$$(WCUT)" ,"p_{Lab}/Z (MeV/c^{2})"  ,400, 0, 3000);
  ejungwoo::variable var_dedx   ("dedx"   ,"dEdx[0].fElements"                                            ,"$$(WCUT)" ,"dE/dx"                  ,400, 0, 2000);

  ejungwoo::variable variables1[] = {
    var_keoaCM,
    var_ptoaCM,
    var_yyCM,
    var_yyCM + var_ptoaCM,
    var_phiCM + var_keoaCM,
    var_ttaCM + var_keoaCM,
    var_phiCM + var_ttaCM,
    var_phiLab + var_ttaLab,
    var_keoaCM + var_pozLab,
    var_eff + "y",
    var_prob + "y",
    var_pozLab + var_dedx + "z",
    var_keoaCM + var_dedx + "z",

    //var_keCM,
    //var_keozCM,
    //var_keCM + var_pozLab,
    //var_pLab + var_dedx + "z",
    //var_keCM + var_dedx + "z",
    //var_keozCM + var_dedx + "z",
  };

  auto tree = new TChain("cbmsim");
  tree -> Add(file_name_in);
  auto num_events = tree -> GetEntries();
  cout << "number of events: " << num_events << endl;
  ejungwoo::setpar("n",num_events);

  TFile *file_hist = new TFile(file_name_out,"recreate");
  (new TParameter<int>("num_events",num_events)) -> Write();

  auto www_x = TCut("x1",      "1./$$(n)");
  auto www_p = TCut("xP",      "$$(prob)/$$(n)");
  auto www_a = TCut("xPoE",    "$$(prob)/$$(eff)/$$(n)");
  auto cut_a = TCut("cPE",     Form("$$(prob)>%f&&$$(eff)>%f",probCut,effCut));
  auto to_70 = TCut("to70",    "$$(ttaCM)>70");
  auto fm_70 = TCut("70to110", "$$(ttaCM)>70&&$$(ttaCM)<110");

  for (auto idx : {0,1,2,3,4})
  {
    ejungwoo::gheader(particle_name[idx]+"_");
    ejungwoo::setpar("i",idx);
    ejungwoo::setpar("m",particle_mass[idx]);
    ejungwoo::setpar("m2",particle_mass[idx]*particle_mass[idx]);
    ejungwoo::setpar("a",particle_a[idx]);
    ejungwoo::setpar("z",particle_z[idx]);

    for (auto weight : {www_x,www_p,www_a}) {
      for (auto cut_ttl : {to_70,fm_70}) {
        TCut wcut = weight * (cut_ttl+cut_a);
        wcut.SetName(Form("%s_%s",weight.GetName(),cut_ttl.GetName()));
        ejungwoo::gwcut(wcut);
        int count_testmode = 0;
        for (auto var : variables1) {
          auto chist = ejungwoo::drawv(var,tree);
          file_hist -> cd();
          chist.hist -> Write();
          if (testmode && count_testmode++>5)
            break;
        }
        if (testmode) break;
      }
      if (testmode) break;
    }
    if (testmode) break;
  }

  file_hist -> ls();

  cout << file_name_out << endl;
}
