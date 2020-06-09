void check_efficiency_data(
  //const char *anaFile = "~/data/pid2/beam108_run2272_2276_ana.root"
  //const char *anaFile = "/home/ejungwoo/data/pid3/beam132_run2841_2846_ana2.root"
  //const char *anaFile = "/home/ejungwoo/data/pid4/Sn132_0_ana.root"
  const char *anaFile = "/home/ejungwoo/data/pid4/Sn108_0_ana.NewAna.2034.45b9400.root"
  )
{
  gStyle -> SetOptStat(0);

  int pdgs[] = {2212 ,1000010020 ,1000010030 ,1000020030 ,1000020040};

  TCanvas *cvs = new TCanvas("cvs","",0,0,2400,1200);
  cvs -> SetLeftMargin(0.01);
  cvs -> Divide(5,2,0.01,0.001);

  Int_t idx = 0;
  auto file = new TFile(anaFile);
  for (auto pdg : pdgs)
  {
    auto teff = (TEfficiency *) file -> Get(Form("Efficiency%d",pdg));

    idx++;

    auto hist1 = (TH3D *) teff -> GetTotalHistogram();
    hist1 -> SetTitle(Form("%s %d;pz;pt",hist1->GetTitle(),pdg));
    cvs -> cd(idx);
    hist1 -> Draw("colz");

    auto hist2 = (TH3D *) teff -> GetPassedHistogram();
    hist2 -> SetTitle(Form("%s %d;pz;pt",hist2->GetTitle(),pdg));
    cvs -> cd(5+idx);
    hist2 -> Draw("colz");
  }

  //cvs -> SaveAs("figures_eff/eff132.pdf");
  cvs -> SaveAs("figures_eff/eff108.png");
}
