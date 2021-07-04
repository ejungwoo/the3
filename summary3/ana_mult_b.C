#include "init_variables.h"
//using namespace ejungwoo;

TH1D *make_sumy(TH1D *hist, const char *name, double scale=1);
TGraph *copy_marker(TH1D *hist, double msize);

void ana_mult_b()
{
  gStyle -> SetOptStat(0);

  const char *version = "f7MA";
  bool reco_or_va = true;

  //const char *version2 = Form("%s_%s",version,(reco_or_va?"reco":"va"));
  TString version2 = Form("%s_%s",version,(reco_or_va?"reco":"va"));

  double bmax[] = {7.52,7.13,7.33,7.31}; // 132, 108, 112, 124
  double bmax_error[] = {0.0254, 0.0223, 0.0195, 0.0283};

  const char *attName = "sysmult";

  binning bnx(80,0,80);
  //binning bnx(100,0,100);

  double mult_max = 0.04;
  double redb_max = 1;
  double realb_max = 8;

  auto cvs_mult = canvas("cvs_mult","multb");
  cvs_mult -> SetGrid();
  auto frame_mult = bnx.newHist("mult_frame",";Charged particle multiplicity #it{N};#it{P}(#it{N})");
  frame_mult -> SetMaximum(0.04);
  //frame_mult -> SetMaximum(0.0372);
  draw(frame_mult,cvs_mult);
  auto legend_mult = new TLegend();

  //reduced impact parameter
  auto cvs_redb = canvas("cvs_redb","multb");
  cvs_redb -> SetGrid();
  //auto frame_redb = bnx.newHist("redb_frame",";Multiplicity cut #it{N}_{C};#hat{#it{b}} = #Sigma_{#it{N} > #it{N}_{C}} #it{P}(#it{N})");
  //auto frame_redb = bnx.newHist("redb_frame",";Multiplicity cut #it{N}_{C};#hat{#it{b}} = #Sigma_{#it{N} #geq #it{N}_{C}} #it{P}(#it{N})");
  auto frame_redb = bnx.newHist("redb_frame",";Multiplicity cut #it{N}_{C};#hat{#it{b}} = #Sigma_{#it{N} = #it{N}_{C}}^{#infty} #it{P}(#it{N})");
  frame_redb -> SetMaximum(redb_max);
  draw(frame_redb,cvs_redb);
  auto legend_redb = new TLegend();

  //real b
  auto cvs_realb = canvas("cvs_realb","multb");
  cvs_realb -> SetGrid();
  auto frame_realb = bnx.newHist("realb_frame",";Multiplicity cut #it{N}_{C};#it{b} = #hat{#it{b}} #times #it{b}_{max} (fm)");
  frame_realb -> SetMaximum(realb_max);
  //frame_realb -> GetXaxis() -> SetRangeUser(10,100);
  draw(frame_realb,cvs_realb);
  auto legend_realb = new TLegend();
  TH1D *hist_realb[4];

  auto cvs_all = canvas("cvs_all",3,4,"snn2");

  for (auto iSys : {0,1,2,3}) {
    auto sys = fSysBeams[iSys];
    auto target = fSysTargets[iSys];
    auto tree = new TChain("mult");
    tree -> Add(Form("data2/%s/sys%d_%s_all_*_100.NewAna.2107.4fd2bca.ana.mult.root",version,sys,version));
    auto num_events = tree -> GetEntries();

    auto name_mult = Form("hist_mult_%d",sys);
    auto hist_mult = bnx.newHist(name_mult,";Charged particle multiplicity #it{N};#it{P}(#it{N})");
    if (reco_or_va) tree -> Project(name_mult,"nt_mult", Form("1./%lld",num_events));
    else            tree -> Project(name_mult,"nt_va"  , Form("1./%lld",num_events));
    hist_mult -> SetMaximum(0.04);
    att(hist_mult,iSys,attName);
    draw(hist_mult,cvs_mult,"histplsame");
    legend_mult -> AddEntry(copy_marker(hist_mult,1.2),Form("^{%d}Sn+^{%d}Sn",sys,target),"pl");

    auto name_redb = Form("hist_redb_%d",sys);
    auto hist_redb = make_sumy(hist_mult,name_redb);
    hist_redb -> SetTitle(";Multiplicity cut #it{N}_{C};#hat{#it{b}} = #Sigma_{#it{N} = #it{N}_{C}}^{#infty} #it{P}(#it{N})");

    hist_redb -> SetMaximum(1);
    att(hist_redb,iSys,attName);
    draw(hist_redb,cvs_redb,"plsame");
    legend_redb -> AddEntry(copy_marker(hist_redb,1.2),Form("^{%d}Sn+^{%d}Sn",sys,target),"pl");

    auto name_realb = Form("hist_realb_%d",sys);
    hist_realb[iSys] = make_sumy(hist_mult,name_realb,bmax[iSys]);
    hist_realb[iSys] -> SetTitle(";Multiplicity cut #it{N}_{C};#it{b} = #hat{#it{b}} #times #it{b}_{max} (fm)");
    hist_realb[iSys] -> SetMaximum(8);
    att(hist_realb[iSys],iSys,attName);
    draw(hist_realb[iSys],cvs_realb,"plsame");
    legend_realb -> AddEntry(copy_marker(hist_realb[iSys],1.2),Form("^{%d}Sn+^{%d}Sn",sys,target),"pl");

    {
      TVirtualPad *cvs_cur;

      cvs_cur = cvs_all -> cd(iSys*3+1);
      if (iSys==0)  {
        draw(hist_mult,cvs_all->cd(1*3+1),"histl");
        draw(hist_mult,cvs_all->cd(2*3+1),"histl");
        draw(hist_mult,cvs_all->cd(3*3+1),"histl");
        draw(hist_mult,cvs_all->cd(0*3+1),"histpl");
      }
      else
        draw(hist_mult,cvs_cur,"histplsame");

      cvs_cur -> SetGrid();
      TString sysFull = Form("^{%d}Sn+^{%d}Sn",sys,target);
      auto tt0 = (new TLatex(bnx.xByRatio(0.55),binning(2,0,mult_max).xByRatio(0.85),sysFull));
      tt0 -> SetTextFont(133);
      tt0 -> SetTextSize(24);
      tt0 -> Draw();

      cvs_cur = cvs_all -> cd(iSys*3+2);
      draw(hist_redb,cvs_cur,"pl");
      //if (iSys==0) draw(hist_redb,cvs_all,iSys*3+2,"histpl");
      //else draw(hist_redb,cvs_cur,0,"histplsame");
      cvs_cur -> SetGrid();
      auto tt1 = (new TLatex(bnx.xByRatio(0.55),binning(2,0,redb_max).xByRatio(0.85),sysFull));
      tt1 -> SetTextFont(133);
      tt1 -> SetTextSize(24);
      tt1 -> Draw();

      cvs_cur = cvs_all -> cd(iSys*3+3);
      draw(hist_realb[iSys],cvs_cur,"pl");
      //if (iSys==0) draw(hist_realb[iSys],cvs_all,iSys*3+2,"histpl");
      //else draw(hist_realb[iSys],cvs_cur,0,"histplsame");
      cvs_cur -> SetGrid();
      auto tt2 = (new TLatex(bnx.xByRatio(0.55),binning(2,0,realb_max).xByRatio(0.85),sysFull));
      tt2 -> SetTextFont(133);
      tt2 -> SetTextSize(24);
      tt2 -> Draw();
    }
  }

  cout << version2 << endl;
  const char *name_mult_b = Form("data_mult_b/%s_mult_b.csv",version2.Data());
  cout << name_mult_b << endl;
  ofstream file_mult_b(name_mult_b);
  file_mult_b << version2 << " 132 108 112 124" << endl;
  bnx.reset();
  while (bnx.next()) {
    file_mult_b << bnx.low() << " "
      << hist_realb[0] -> GetBinContent(bnx.bi()) << " "
      << hist_realb[1] -> GetBinContent(bnx.bi()) << " "
      << hist_realb[2] -> GetBinContent(bnx.bi()) << " "
      << hist_realb[3] -> GetBinContent(bnx.bi()) << endl;

    continue;

    cout << bnx.low() << " "
      << hist_realb[0] -> GetBinContent(bnx.bi()) << " "
      << hist_realb[1] -> GetBinContent(bnx.bi()) << " "
      << hist_realb[2] -> GetBinContent(bnx.bi()) << " "
      << hist_realb[3] -> GetBinContent(bnx.bi()) << endl;
  }

  draw(legend_mult,cvs_mult,0.05,-1, .35,.4,0.3);
  draw(legend_redb,cvs_redb,-1,-1,   .35,.4,0.3);
  draw(legend_realb,cvs_realb,-1,-1, .35,.4,0.3);

  //saveAll(version2);
}

TH1D *make_sumy(TH1D *hist, const char *name, double scale)
{
  auto bnx = binning(hist);
  auto hist_sumy = bnx.newHist(name);
  bnx.reset();
  while(bnx.next()) {
    double sumy = 0;
    for (auto bin=bnx.bi(); bin<=bnx.fN; ++bin)
    //for (auto bin=bnx.bi()+1; bin<=bnx.fN; ++bin)
      sumy += hist -> GetBinContent(bin);
    hist_sumy -> SetBinContent(bnx.bi(),sqrt(sumy)*scale);
  }
  return hist_sumy;
}


TGraph *copy_marker(TH1D *hist, double msize) {
  auto graph = new TGraph();
  graph -> SetMarkerStyle(hist->GetMarkerStyle());
  graph -> SetMarkerColor(hist->GetMarkerColor());
  graph -> SetLineColor(hist->GetLineColor());
  graph -> SetMarkerSize(msize);
  return graph;
}
