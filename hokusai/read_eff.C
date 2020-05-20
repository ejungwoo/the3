void read_eff()
{
  gStyle -> SetOptStat(0);

  int pdgs[] = {2212 ,1000010020 ,1000010030 ,1000020030 ,1000020040 ,1000020060};

  TCanvas *cvs = new TCanvas("cvs","",0,0,1000,700);
  cvs -> SetLeftMargin(0.01);
  cvs -> Divide(6,2,0.01,0.001);

  Int_t idx = 0;
  auto file = new TFile("~/data/pid2/beam108_run2272_2276_ana.root");
  for (auto pdg : pdgs)
  {
    auto teff = (TEfficiency *) file -> Get(Form("Efficiency%d",pdg));

    idx++;

    auto hist1 = (TH3D *) teff -> GetTotalHistogram();
    hist1 -> SetTitle(Form("%s %d;z;pt",hist1->GetTitle(),pdg));
    cvs -> cd(idx);
    hist1 -> Draw("colz");

    auto hist2 = (TH3D *) teff -> GetPassedHistogram();
    hist2 -> SetTitle(Form("%s %d;z;pt",hist2->GetTitle(),pdg));
    cvs -> cd(6+idx);
    hist2 -> Draw("colz");
  }
}
