void draw_rapidity_pt()
{
  bool test = false;
  int system = 108;

  ejungwoo::gzcolor(1);
  ejungwoo::gstat(0);
  ejungwoo::gcvspos(1300,0);
  ejungwoo::gversion("ty");
  ejungwoo::gversion(Form("sys%d",system));

  ejungwoo::binning binningNy{200, -1.5, 3.0};
  ejungwoo::binning binningPt{200, 0, 1000};
  ejungwoo::binning binningKe{200, 0, 200};

  auto tree = new TChain("data");
  bool printTree = true;
  //for (auto iSummary=0; iSummary<15; ++iSummary) {
  //for (auto iSummary : {0}) {

  if (system==132)
  for (auto iSummary=0; iSummary<23; ++iSummary) {
    tree -> AddFile(Form("data_from_hokusai/summary_132_%d.root",iSummary));
    if (printTree) {
      tree -> Print("top");
      printTree = false;
    }
    if (test)
      break;
  }
  else if (system==108)
  for (auto iSummary=0; iSummary<10; ++iSummary) {
    tree -> AddFile(Form("data_from_hokusai/summary_108_%d.root",iSummary));
    if (printTree) {
      tree -> Print("top");
      printTree = false;
    }
    if (test)
      break;
  }

  bool isGood;
  int pdg;
  double prob, eff, y, ybeam;

  tree -> SetBranchAddress("good",&isGood);
  tree -> SetBranchAddress("pdg",&pdg);
  tree -> SetBranchAddress("prob",&prob);
  tree -> SetBranchAddress("eff",&eff);
  tree -> SetBranchAddress("y",&y);
  tree -> SetBranchAddress("yBeam",&ybeam);

  auto line0 = new TLine(0,0,0,binningPt.max);
  line0 -> SetLineColor(kBlack);
  line0 -> SetLineStyle(9);

  auto numPDGs = 6;
  const std::vector<int> listPDGs{2212, 1000010020, 1000010030, 1000020030, 1000020040, 1000020060};
  const char *listNames[] = {"p", "d", "t", "he3", "he4", "he6"};
  const int listA[] = {1,2,3,3,4,6};

  //for (auto ipdg=0; ipdg<numPDGs;
  for (auto ipdg=0; ipdg<5; ++ipdg)
  //for (auto ipdg : {0})
  {
    auto pdg = listPDGs[ipdg];
    TString pdgname = TString(listNames[ipdg]) + Form("_%d",system);
    int particleA = listA[ipdg];

    TString name = pdgname;
    //auto cvs = ejungwoo::div(ejungwoo::cc(name,600,1200),1,2);
    auto cvs = ejungwoo::div(ejungwoo::cc(name,1200,1200),2,2);

    cvs -> cd(1);
    ejungwoo::gfooter("_nypt");
    ejungwoo::titles titles1{pdgname,"y_{CM}/y_{beam,CM}","p_{T}/A"};
    //TCut cut1(Form("(eff>0.01&&pdg==%d&&prob>0.95)",pdg));
    TCut cut1(Form("(pdg==%d&&prob>0.95)",pdg));
    auto hist1 = ejungwoo::tp(name,tree,Form("pt/%d:y/yBeam",particleA),cut1,titles1.data(),binningNy,binningPt);
    hist1 -> Draw("colz");
    line0 -> Draw("samel");

    ejungwoo::gfooter("_nypt1");
    TCut cutx(Form("eff>0.01&&pdg==%d&&prob>0.95",pdg));
    auto histx = ejungwoo::tp(name,tree,Form("pt/%d:y/yBeam",particleA),cutx,titles1.data(),binningNy,binningPt);

    cvs -> cd(2);// -> SetLogz();
    ejungwoo::gfooter("_nypt_Ceff");
    ejungwoo::titles titles2{pdgname+" efficiency-corrected","y_{CM}/y_{beam,CM}","p_{T}/A"};
    TCut cut2 = TCut("1./eff")*cutx;
    auto hist2 = ejungwoo::tp(name,tree,Form("pt/%d:y/yBeam",particleA),cut2,titles2.data(),binningNy,binningPt);
    hist2 -> Draw("colz"); 
    line0 -> Draw("samel");
    ejungwoo::gfooter("");

    cvs -> cd(3);
    ejungwoo::gfooter("_eff");
    TString name3 = TString(histx->GetName());
    name3.ReplaceAll("nypt1","eff");
    ejungwoo::titles titles3{pdgname+" efficiency plot","y_{CM}/y_{beam,CM}","p_{T}/A"};
    auto hist3 = (TH2D *) histx -> Clone(name3);
    hist3 -> Divide(hist2);
    hist3 -> SetTitle(titles3.data());
    hist3 -> Draw("colz");

    cvs -> cd(4);
    ejungwoo::titles titles4{pdgname,"KE_{CM}/A",""};
    TCut cut5 = cutx;
    TCut cut6 = TCut("1./(eff*3)")*cutx;

    auto legend = new TLegend();

    ejungwoo::gfooter("_ke");
    auto hist4 = ejungwoo::tp(name,tree,Form("ke/%d",particleA),cut1,titles4.data(),binningKe);
    legend -> AddEntry(hist4,"prob>0.95","l");

    ejungwoo::gfooter("_ke");
    auto hist5 = ejungwoo::tp(name,tree,Form("ke/%d",particleA),cut5,titles4.data(),binningKe);
    legend -> AddEntry(hist5,"prob>0.95 eff>0.01","l");

    ejungwoo::gfooter("_ke_Ceff");
    auto hist6 = ejungwoo::tp(name,tree,Form("ke/%d",particleA),cut6,titles4.data(),binningKe);
    legend -> AddEntry(hist6,"prob>0.95 eff>0.01, eff.#times3 cor.","l");

    hist4 -> SetLineColor(kGreen);
    hist5 -> SetLineColor(kBlue);
    hist6 -> SetLineColor(kRed);

    hist4 -> Draw("hist"); 
    hist6 -> Draw("same hist"); 
    hist5 -> Draw("same hist"); 

    ejungwoo::gfooter("");
    ejungwoo::make_c(cvs);
    ejungwoo::make_l(legend) -> Draw();

    ejungwoo::gfooter("");
  }

  if (!test)
    ejungwoo::savecvsall("png");
}
