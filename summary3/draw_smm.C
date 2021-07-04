#include "init_variables.h"
#include "draw_r21.C"

void draw_smm()
{
  gStyle -> SetOptStat(0);

  double dataValues[5][2] =
  /*
  { {0.84625,0.725867},
    {1.13812,1.02502},
    {1.54882,1.38822},
    {0.93565,0.655999},
    {1.21135,0.920377} }; // pt/a < 250
  */
  { {0.835549, 0.707272},
    {1.10853,  1.01954},
    {1.53436,  1.41199},
    {0.923894, 0.646844},
    {1.19937,  0.925081}, }; // pt/a < 200

  double refValues[5];
  TString smmTitle0;

  int selSMM = 7;
  //auto fNumSMM = 11;
  //for (auto idxSMM=3; idxSMM<fNumSMM; ++idxSMM)

  if (0) {
    binning bnx(100,0,4);
    binning bny(100,0.5,2.0,"R_{21}");
    auto bnbn = bnx*bny;

    auto cvs = canvas("smm_a1", "smm_compare2");
    auto hist = bnbn.newHist("hist"); 
    draw(hist,cvs);
    TGraphErrors *graphc[5];
    graphc[0] = new TGraphErrors();
    graphc[1] = new TGraphErrors();
    graphc[2] = new TGraphErrors();
    graphc[3] = new TGraphErrors();
    graphc[4] = new TGraphErrors();

    int idx = 0;
    for (auto idxSMM : {3,4,5,6})
    {
      TString fileName;
      TString smmTitle;
      TString smmChange;
      if (idxSMM==3)  { fileName = "data2/smm/t7_nz157_vv3.dat";  smmTitle = "T=7.0,  #frac{V_{bk}}{V_{0}} = 3,  #frac{N_{2}}{Z_{2}} = 1.57"; smmChange = "Decrease temperature"; }
      else if (idxSMM==4)  { fileName = "data2/smm/t7_nz157_vv9.dat";  smmTitle = "T=7.0,  #frac{V_{bk}}{V_{0}} = 9,  #frac{N_{2}}{Z_{2}} = 1.57"; smmChange = ""; }
      else if (idxSMM==5)  { fileName = "data2/smm/t7_nz138_vv3.dat";  smmTitle = "T=7.0,  #frac{V_{bk}}{V_{0}} = 3,  #frac{N_{2}}{Z_{2}} = 1.38"; smmChange = ""; }
      else if (idxSMM==6)  { fileName = "data2/smm/t7_nz138_vv9.dat";  smmTitle = "T=7.0,  #frac{V_{bk}}{V_{0}} = 9,  #frac{N_{2}}{Z_{2}} = 1.38"; smmChange = ""; }

      ifstream smmfile(fileName);
      double x0, r21Value;
      for (auto iParticle : {0,1,2,3,4}) {
        smmfile >> x0 >> r21Value;
        if (idxSMM==3) {
          auto tt = new TLatex(0.1,r21Value,fParticleTitles[iParticle]);
          tt -> SetTextAlign(22);
          tt -> SetTextFont(133);
          tt -> SetTextSize(15);
          tt -> SetTextColor(fDrawColor[iParticle]);
          tt -> Draw("samel");
        }
        graphc[iParticle] -> SetPoint(2*idx,  idx+0.2, r21Value);
        graphc[iParticle] -> SetPoint(2*idx+1,idx+0.8, r21Value);
        attGraph(graphc[iParticle],iParticle);
        graphc[iParticle] -> SetLineWidth(2);
      }
      idx++;

      auto paveText = newpt(smmTitle, cvs, idx+0.5, 0.75, .5, .25);
      paveText -> Draw("same");
    }

    graphc[0] -> Draw("samel");
    graphc[1] -> Draw("samel");
    graphc[2] -> Draw("samel");
    graphc[3] -> Draw("samel");
    graphc[4] -> Draw("samel");
  return;
  }

  binning bnx(100,0,1);
  binning bny(100,0.2,2.0,"R_{21}");
  auto bnbn = bnx*bny;

  for (auto idxSMM : {7,3,8,9})
  {
    TString fileName;
    TString smmTitle;
    TString smmChange;
         if (idxSMM==0)  { fileName = "data2/smm/srs_1.dat";         smmTitle = "T=8.7,  participant"; }
    else if (idxSMM==1)  { fileName = "data2/smm/srs_2.dat";         smmTitle = "T=8.7,  0.2#timespart + proj"; }
    else if (idxSMM==2)  { fileName = "data2/smm/srs_3.dat";         smmTitle = "T=8.7,  0.1#timespart + proj"; }
    else if (idxSMM==3)  { fileName = "data2/smm/t7_nz157_vv3.dat";  smmTitle = "T=7.0,  #frac{V_{bk}}{V_{0}} = 3,  #frac{N_{2}}{Z_{2}} = 1.57"; smmChange = "Decrease temperature"; }
    else if (idxSMM==4)  { fileName = "data2/smm/t7_nz157_vv9.dat";  smmTitle = "T=7.0,  #frac{V_{bk}}{V_{0}} = 9,  #frac{N_{2}}{Z_{2}} = 1.57"; smmChange = ""; }
    else if (idxSMM==5)  { fileName = "data2/smm/t7_nz138_vv3.dat";  smmTitle = "T=7.0,  #frac{V_{bk}}{V_{0}} = 3,  #frac{N_{2}}{Z_{2}} = 1.38"; smmChange = ""; }
    else if (idxSMM==6)  { fileName = "data2/smm/t7_nz138_vv9.dat";  smmTitle = "T=7.0,  #frac{V_{bk}}{V_{0}} = 9,  #frac{N_{2}}{Z_{2}} = 1.38"; smmChange = ""; }
    else if (idxSMM==7)  { fileName = "data2/smm/t8_nz157_vv3.dat";  smmTitle = "T=8.0,  #frac{V_{bk}}{V_{0}} = 3,  #frac{N_{2}}{Z_{2}} = 1.57"; smmChange = ""; }
    else if (idxSMM==8)  { fileName = "data2/smm/t8_nz157_vv9.dat";  smmTitle = "T=8.0,  #frac{V_{bk}}{V_{0}} = 9,  #frac{N_{2}}{Z_{2}} = 1.57"; smmChange = "Increase volume"; }
    else if (idxSMM==9)  { fileName = "data2/smm/t8_nz138_vv3.dat";  smmTitle = "T=8.0,  #frac{V_{bk}}{V_{0}} = 3,  #frac{N_{2}}{Z_{2}} = 1.38"; smmChange = "Decrease N/Z"; }
    else if (idxSMM==10) { fileName = "data2/smm/t8_nz138_vv9.dat";  smmTitle = "T=8.0,  #frac{V_{bk}}{V_{0}} = 9,  #frac{N_{2}}{Z_{2}} = 1.38"; smmChange = "Change volume and N/Z"; }
    if (idxSMM==selSMM) smmTitle0 = smmTitle;

    auto cvs = canvas(Form("smm_v%d",idxSMM), "smm_compare");
    auto hist0 = bnbn.newHist(Form("smm_1_v%d",idxSMM));
    hist0 -> SetTitle(smmChange);
    draw(hist0,cvs);
    auto line1 = new TLine(0,1,1,1);
    auto line2 = new TLine(0,1.5,1,1.5);
    auto line3 = new TLine(.5,.2,.5,2.0);
    for (auto line : {line1, line2, line3}) {
      line -> SetLineColor(kGray+1);
      line -> SetLineStyle(2);
      //line -> Draw("same");
    }

    ifstream smmfile(fileName);
    double x0, r21Value;
    for (auto iParticle : {0,1,2,3,4}) {
      smmfile >> x0 >> r21Value;
      if (idxSMM==selSMM) refValues[iParticle] = r21Value;

      auto graphc = new TGraphErrors();
      auto tt = new TLatex(0.04,refValues[iParticle],fParticleTitles[iParticle]);
      tt -> SetTextAlign(22);
      tt -> SetTextFont(133);
      tt -> SetTextSize(15);
      tt -> SetTextColor(fDrawColor[iParticle]);
      tt -> Draw("samel");
      graphc -> SetPoint(0,.08,refValues[iParticle]);
      graphc -> SetPoint(1,.48,refValues[iParticle]);
      graphc -> SetPoint(2,.52,r21Value);
      graphc -> SetPoint(3,0.92,r21Value);
      attGraph(graphc,iParticle);
      graphc -> SetLineWidth(2);
      graphc -> Draw("samel");

      double data0 = dataValues[iParticle][0];
      double data1 = dataValues[iParticle][1];
      auto zValue = fParticleZ[iParticle];
      zValue = 0;
      TGraph *graphyy1 = new TGraph();
      graphyy1 -> SetPoint(graphyy1->GetN(),0.15-0.02*zValue,data1);
      graphyy1 -> SetPoint(graphyy1->GetN(),0.35-0.02*zValue,data0);
      graphyy1 -> SetMarkerStyle(fDrawMStyle[iParticle]);
      graphyy1 -> SetMarkerColor(fDrawColor[iParticle]);
      graphyy1 -> SetLineColor(fDrawColor[iParticle]);
      graphyy1 -> SetMarkerSize(1.8);
      graphyy1 -> SetLineStyle(3);
      graphyy1 -> Draw("samepl");

      TGraph *graphyy2 = new TGraph();
      graphyy2 -> SetPoint(graphyy2->GetN(),0.65-0.02*zValue,data1);
      graphyy2 -> SetPoint(graphyy2->GetN(),0.85-0.02*zValue,data0);
      graphyy2 -> SetMarkerStyle(fDrawMStyle[iParticle]);
      graphyy2 -> SetMarkerColor(fDrawColor[iParticle]);
      graphyy2 -> SetLineColor(fDrawColor[iParticle]);
      graphyy2 -> SetMarkerSize(1.8);
      graphyy2 -> SetLineStyle(3);
      graphyy2 -> Draw("samepl");
      /*
      TMarker *marker1 = new TMarker(0.15-0.02*zValue,data1,fDrawMStyle[iParticle]);
      TMarker *marker0 = new TMarker(0.35-0.02*zValue,data0,fDrawMStyle[iParticle]);
      TMarker *marker3 = new TMarker(0.65-0.02*zValue,data1,fDrawMStyle[iParticle]);
      TMarker *marker2 = new TMarker(0.85-0.02*zValue,data0,fDrawMStyle[iParticle]);
      for (auto marker : {marker0, marker1, marker2, marker3}) {
        marker -> SetMarkerColor(fDrawColor[iParticle]);
        marker -> SetMarkerSize(1.8);
        marker -> Draw("samep");
      }
      */
    }

    auto paveText1 = newpt(smmTitle0, cvs, 0, 0.75, .5, .25);
    auto paveText2 = newpt(smmTitle, cvs, 0.5, 0.75, .5, .25);
    paveText1 -> SetTextAlign(22);
    paveText2 -> SetTextAlign(22);
    auto paveText_y1 = newpt("y_{0}~1", cvs, 0.15-0.1, 0.0, .2, .15); paveText_y1 -> SetTextAlign(22);
    auto paveText_y2 = newpt("y_{0}~0", cvs, 0.35-0.1, 0.0, .2, .15); paveText_y2 -> SetTextAlign(22);
    auto paveText_y3 = newpt("y_{0}~1", cvs, 0.65-0.1, 0.0, .2, .15); paveText_y3 -> SetTextAlign(22);
    auto paveText_y4 = newpt("y_{0}~0", cvs, 0.85-0.1, 0.0, .2, .15); paveText_y4 -> SetTextAlign(22);
    auto paveText_y5 = newpt("y_{0}~1  #rightarrow  y_{0}~0", cvs, 0.25-0.2, 0.0, .4, .15); paveText_y5 -> SetTextAlign(22);
    auto paveText_y6 = newpt("y_{0}~1  #rightarrow  y_{0}~0", cvs, 0.75-0.2, 0.0, .4, .15); paveText_y6 -> SetTextAlign(22);

    //for (auto ptt : {paveText1,paveText2,paveText_y1,paveText_y2,paveText_y3,paveText_y4}) {
    for (auto ptt : {paveText1,paveText2,paveText_y5,paveText_y6}) {
      ptt -> SetTextSize(20);
      ptt -> Draw("same");
    }
  }
}
