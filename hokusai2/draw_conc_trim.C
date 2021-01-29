void draw_conc_trim()
{
  //const char *anaName = "fix5";
  const char *anaName = "fix6";

  TCut cutNClDPoca = "vaNRowClusters + vaNLayerClusters > 15 && recodpoca.Mag() < 15";
  TCut cutVertexZ = "fabs(tpcVertex.z() + 15) < 5";
  TCut cutTheta2 = "vaMom.Theta()*TMath::RadToDeg() > 60";
  TCut cutTheta1 = "vaMom.Theta()*TMath::RadToDeg() <= 60";

  for (int sys : {108,132})
  //for (int sys : {108})
  {
    TString name = Form("pid_%s_%d",anaName,sys);
    auto tree = new TChain("cbmsim");
    tree -> Add(Form("/home/ejungwoo/data/trim/%s/Sn%d/run*_s*.reco.*.conc.trimmed.root",anaName,sys));

    int idx = 0;
    for (auto cut :  {cutNClDPoca && cutTheta2, cutNClDPoca && cutTheta1}) {
      TString name1 = name + "_" + idx;
      cout << name1 << endl;
      auto cvs = new TCanvas(name1);
      tree -> Draw(Form("vadedx:vaMom.Mag()>>%s_%d(200,0,3000,200,0,1000)",name1.Data(),idx), cut, "colz");
      cvs -> SetLogz();
      cvs -> SaveAs(Form("figures/%s.png",name1.Data()));
      idx++;
    }

  }
} 

