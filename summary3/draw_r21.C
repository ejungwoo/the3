#include "KBGlobal.hh"
#include "init_variables.h"
#include "binning.h"

TString fNameV;

int fiSys, fiSys2, fiSys1, fiComb = -1, fiParticle;
int fCountCvs = 0;
int fXOffCvs = 1200;
int fdxCvsAB = 150;
int fdyCvsAB = 200;
vector<TCanvas *> fCvsArray;

TH2D *fHYValue[4][5][3] = {0};
TH2D *fHYLoose[4][5][3] = {0};

double fLMargin = 0.12, fRMargin = 0.055, fBMargin = 0.10, fTMargin = 0.11;

const char *make_name(const char *name0, int idx=0, const char *tag="");
const char *make_title(binning bn_range=binning(), int idx=-1);
void project(TTree *tree, TH1 *hist, const char *expression, const char *selection="");
TF2 *fit_r21(TGraph2DErrors *graph_r21);
TCanvas *make_cvs(const char *name, int opt=0, int w=630, int h=550, int nx=0, int ny=0);
void save_cvs();
TGraphErrors *attGraph(TGraphErrors *graph, int idx);
TMarker *attMarker(TMarker *marker, int idx);
TH2D *make_r21_frame(TH2D *hist);
TF1 *get_fit1_r21(int nValue, int zValue, double alpha, double beta, double cnorm, double x1, double x2);
void draw_ab_fit(binning bnx, binning bny, double r21Values[8][8][5][2], double abcValues[8][8][3][2], double tptValues[8][8][2][2]);
void draw_r21_x(binning bnx, binning bny, double r21Values[8][8][5][2], double abcValues[8][8][3][2], bool do_for_y=true);
void draw_apmb(binning bnx, binning bny, double r21Values[8][8][5][2], double abcValues[8][8][3][2], double tptValues[8][8][2][2], bool multt=0);
void draw_avb(binning bnx, binning bny, double r21Values[8][8][5][2], double abcValues[8][8][3][2]);
TLegend *make_legend(TVirtualPad *cvs, TLegend *legend, TString opt = "", double x_offset=0, double y_offset=0, double width_fixed=0, double height_fixed=0);

void draw_r21(int fnn=8, int recreate_hist=0, int hold_all=0, int add_all=0)
{
  bool set_recreate_hist = recreate_hist;

  bool set_draw1 = kSet;
  bool set_draw2 = kUnSet;

  bool set_draw_yield = kHold;
  bool set_draw_r21_x = kHold;
  bool set_draw_apmb = kHold;
  bool set_draw_apmbt = kHold;
  bool set_draw_ab_fit = kDraw;
  bool set_draw_avb = kHold;

  if (hold_all) {
    set_draw1 = kUnSet;
    set_draw2 = kUnSet;
    set_draw_yield = kHold;
    set_draw_r21_x = kHold;
    set_draw_apmb = kHold;
    set_draw_apmbt = kHold;
    set_draw_ab_fit = kHold;
    set_draw_avb = kHold;
  }

  const char *anaV = "nn50";

  binning bn_pozl(200,0,2000,"p_{Lab.}/Z","pozl");
  binning bn_dedx(200,0,1500,"dE/dx","dedx");

  binning bn_y0(fnn,0,1,"y_{0}","y0");
  binning bn_pt(fnn,0,400,"p_{T}/A (MeV/c)","ptoac");
  binning bn_ke(fnn,0,120,"KE_{CM}/A (MeV)","keoac");
  binning bn_tt(fnn,0,90,"#theta_{CM} (deg)","ttc");

  vector<int> selSysIdx = {0,1,2,3};
  //vector<int> selSysCombIdx = {3};
  vector<int> selSysCombIdx = {0,1,2,3,4,5};

  //vector<binning2> var2Array = {(bn_ke*bn_tt)};
  //vector<binning2> var2Array = {(bn_pt*bn_y0)};
  vector<binning2> var2Array = {(bn_ke*bn_tt), (bn_pt*bn_y0)};

  if (add_all) {
    var2Array.clear();
    var2Array.push_back(bn_ke*bn_tt);
    var2Array.push_back(bn_pt*bn_y0);
    selSysIdx.clear(); for (auto i : {0,1,2,3}) selSysIdx.push_back(i);
    selSysCombIdx.clear(); for (auto i : {0,1,2,3,4,5}) selSysCombIdx.push_back(i);
  }

  vector<binning2> var2ArrayPP;
  for (auto var2 : var2Array) var2ArrayPP.push_back(var2);
  var2ArrayPP.push_back((bn_pozl*bn_dedx));

  auto findv = [&var2ArrayPP](binning2 var2) {
    for (auto idx=0; idx<var2ArrayPP.size(); ++idx) {
      if (strcmp(var2ArrayPP[idx].name(),var2.name())==0)
        return idx;
    }
    return -1;
  };

  fNameV = Form("%s_%d_%d_%d_%d",anaV,bn_y0.fN,bn_pt.fN,bn_ke.fN,bn_tt.fN);
  gSystem -> mkdir(fNameV);
  gSystem -> mkdir(Form("%s/figures",fNameV.Data()));

  gStyle -> SetOptStat(0);

  for (auto iSys : selSysIdx)
  {
    fiSys = iSys;

    auto nameV = Form("%s.hist_%s_%s_%s_%s",anaV,bn_y0.nnmm(),bn_pt.nnmm(),bn_ke.nnmm(),bn_tt.nnmm());
    auto name_file_hist = Form("data2/%s/compact_sys%d_%s.root",anaV,fSysBeams[fiSys],nameV);
    TFile *file_hist = new TFile(name_file_hist,"read");
    cout_info << "hist file: " << name_file_hist << endl;

    if (set_recreate_hist || file_hist->IsZombie())
    {
      auto name_file_in = Form("data2/%s/compact_sys%d_%s.root",anaV,fSysBeams[fiSys],anaV);
      auto file_tree = new TFile(name_file_in);
      cout_info << "compact file: " << name_file_in  << endl;
      file_hist = new TFile(name_file_hist,"recreate");

      for (auto iParticle : fParticleIdx)
      {
        fiParticle = iParticle;

        auto tree = (TTree *) file_tree -> Get(fParticleNames[fiParticle]);

        for (auto var2 : var2ArrayPP)
        {
          auto hist_vv0 = var2.newHist(make_name(var2.name(),0)); project(tree,hist_vv0,var2.ee(),"corr*(prob>.7)");
          auto hist_vv1 = var2.newHist(make_name(var2.name(),1)); project(tree,hist_vv1,var2.ee(),"corr");

          fHYValue[fiSys][fiParticle][findv(var2)] = hist_vv0;
          fHYLoose[fiSys][fiParticle][findv(var2)] = hist_vv1;

          cout_info << fiSys2 << " " <<  fiParticle << " " << findv(var2) << " " << hist_vv0 -> GetName() << endl;
          cout_info << fiSys2 << " " <<  fiParticle << " " << findv(var2) << " " << hist_vv1 -> GetName() << endl;

          file_hist -> cd();
          hist_vv0 -> Write();
          hist_vv1 -> Write();
        }
      }
      file_tree -> Close();
    }
    else
    {
      for (auto iParticle : fParticleIdx)
      {
        fiParticle = iParticle;

        for (auto var2 : var2ArrayPP)
        {
          auto nameValue = make_name(var2.name(),0);
          auto nameLoose = make_name(var2.name(),1);
          fHYValue[fiSys][fiParticle][findv(var2)] = (TH2D *) file_hist -> Get(nameValue);
          fHYLoose[fiSys][fiParticle][findv(var2)] = (TH2D *) file_hist -> Get(nameLoose);
          auto histValue = fHYValue[fiSys][fiParticle][findv(var2)];
          auto histLoose = fHYLoose[fiSys][fiParticle][findv(var2)]; 
          if (histValue!=nullptr) cout_info << fiSys << " " << fiParticle << " " << findv(var2) << " " << histValue->GetName() << endl; else cout_error << "nullptr " << nameValue << endl;
          if (histLoose!=nullptr) cout_info << fiSys << " " << fiParticle << " " << findv(var2) << " " << histLoose->GetName() << endl; else cout_error << "nullptr " << nameLoose << endl;
        }
      }
    }
  }

  if (set_draw_yield)
  {
    for (auto iSys : selSysIdx)
    {
      for (auto var2 : var2ArrayPP)
      {
        fiSys = iSys;

        auto addAll = false;
        if (strcmp(var2.name(),"dedxpozl")==0)
          addAll = true;

        TH2D *histAll1 = nullptr;
        TH2D *histAll2 = nullptr;

        for (auto iParticle : fParticleIdx)
        {
          fiParticle = iParticle;
          auto hist1 = fHYValue[fiSys][fiParticle][findv(var2)];
          if (addAll) {
            if (histAll1 == nullptr) histAll1 = hist1;
            else histAll1 -> Add(hist1);
          }
          else {
            make_cvs(hist1->GetName(),2);
            hist1 -> Draw("colz");
          }

          auto hist2 = fHYLoose[fiSys][fiParticle][findv(var2)];
          if (addAll) {
            if (histAll2 == nullptr) histAll2 = hist2;
            else histAll2 -> Add(hist2);
          }
          else {
            make_cvs(hist2->GetName(),2);
            hist2 -> Draw("colz");
          }
        }
        if (addAll) {
          make_cvs(histAll1->GetName(),2) -> SetLogz();
          histAll1 -> Draw("colz");

          make_cvs(histAll2->GetName(),2) -> SetLogz();
          histAll2 -> Draw("colz");
        }
      }
    }
  }

  /////////////////////////////////////////////////////
  /////////////////////////////////////////////////////
  /////////////////////////////////////////////////////
  /////////////////////////////////////////////////////

  if (set_draw1)
  {
    for (int ivar : {0,1})
    {
      for (auto var2 : var2Array)
      {
        auto bnx = var2.bx();
        auto bny = var2.by();

        for (auto iComb : selSysCombIdx)
        {
          fiComb = iComb;
          fiSys2 = fSysCombIdx[fiComb][0];
          fiSys1 = fSysCombIdx[fiComb][1];

          auto bnv = bnx;
          auto bns = bny;
          if (ivar==1) {
            bnv = bny;
            bns = bnx;
          }

          double tptValues[8][8][2][2] = {{0}};
          double r21Values[8][8][5][2] = {{0}};
          double abcValues[8][8][3][2] = {{0}};

          bnv.reset();
          while (bnv.next())
          {
            auto graph_r21 = new TGraph2DErrors();

            double yieldValues[5][2] = {0};

            for (auto iParticle : fParticleIdx)
            {
              fiParticle = iParticle;

              double value2 = 0;
              double value1 = 0;
              double loose2 = 0;
              double loose1 = 0;

              if (ivar==0) {
                bns.reset(); while (bns.next()) { value2 += fHYValue[fiSys2][fiParticle][findv(var2)] -> GetBinContent(bnv.bi(),bns.bi()); }
                bns.reset(); while (bns.next()) { value1 += fHYValue[fiSys1][fiParticle][findv(var2)] -> GetBinContent(bnv.bi(),bns.bi()); }
                bns.reset(); while (bns.next()) { loose2 += fHYLoose[fiSys2][fiParticle][findv(var2)] -> GetBinContent(bnv.bi(),bns.bi()); }
                bns.reset(); while (bns.next()) { loose1 += fHYLoose[fiSys1][fiParticle][findv(var2)] -> GetBinContent(bnv.bi(),bns.bi()); }
              }
              else {
                bns.reset(); while (bns.next()) { value2 += fHYValue[fiSys2][fiParticle][findv(var2)] -> GetBinContent(bns.bi(),bnv.bi()); }
                bns.reset(); while (bns.next()) { value1 += fHYValue[fiSys1][fiParticle][findv(var2)] -> GetBinContent(bns.bi(),bnv.bi()); }
                bns.reset(); while (bns.next()) { loose2 += fHYLoose[fiSys2][fiParticle][findv(var2)] -> GetBinContent(bns.bi(),bnv.bi()); }
                bns.reset(); while (bns.next()) { loose1 += fHYLoose[fiSys1][fiParticle][findv(var2)] -> GetBinContent(bns.bi(),bnv.bi()); }
              }

              yieldValues[iParticle][0] = value2;
              yieldValues[iParticle][1] = value1;

              auto r21Value = value2/value1;
              auto r21Loose = loose2/loose1;
              auto r21Error = abs(r21Value - r21Loose);

              r21Values[bnv.ii()][0][fiParticle][kVal] = r21Value;
              r21Values[bnv.ii()][0][fiParticle][kErr] = r21Error;

              auto zValue = fParticleZ[fiParticle];
              auto nValue = fParticleN[fiParticle];

              auto idxR21 = graph_r21 -> GetN();
              graph_r21 -> SetPoint(idxR21, nValue, zValue, r21Value);
              graph_r21 -> SetPointError(idxR21, 0, 0, r21Value);
            }

            for (auto iss : {0,1}) {
              auto yieldd = yieldValues[kD][iss];
              auto yield4 = yieldValues[kHe4][iss];
              auto yieldt = yieldValues[kT][iss];
              auto yield3 = yieldValues[kHe3][iss];
              tptValues[bnv.ii()][0][iss][kVal] = 14.3/TMath::Log(1.59*(yieldd * yield4 / yieldt / yield3));
            }

            auto fitIsoscaling = fit_r21(graph_r21);
            double alpha = fitIsoscaling -> GetParameter(0);
            double beta = fitIsoscaling -> GetParameter(1);
            double cnorm = fitIsoscaling -> GetParameter(2);
            abcValues[bnv.ii()][0][0][kVal] = alpha;
            abcValues[bnv.ii()][0][1][kVal] = beta;
            abcValues[bnv.ii()][0][2][kVal] = cnorm;

            double alphaError = 0.;
            double betaError = 0.;

            for (auto iParticle : fParticleIdx)
            {
              fiParticle = iParticle;

              auto nValue = fParticleN[iParticle];
              auto zValue = fParticleZ[iParticle];

              auto r21Value = r21Values[bnv.ii()][0][fiParticle][kVal];
              auto r21Eval = fitIsoscaling -> Eval(nValue, zValue);
              auto alphaError1 = (r21Value - r21Eval) / r21Value / nValue;
              auto betaError1 = (r21Value - r21Eval) / r21Value / zValue;

              if (nValue!=0) alphaError += alphaError1*alphaError1;
              if (zValue!=0) betaError += betaError1*betaError1;
            }
            alphaError = sqrt(alphaError / 4);
            betaError = sqrt(betaError / 5);
            abcValues[bnv.ii()][0][0][kErr] = alphaError;
            abcValues[bnv.ii()][0][1][kErr] = betaError;
          }

          auto bn2 = bns;
          bn2.setN(1);
          bn2.setExpression("all");

          if (set_draw_ab_fit) draw_ab_fit(bnv, bn2, r21Values, abcValues, tptValues);
          if (set_draw_r21_x) draw_r21_x(bnv, bn2, r21Values, abcValues, 0);
          if (set_draw_apmb) draw_apmb(bnv, bn2, r21Values, abcValues, tptValues);
          if (set_draw_apmbt) draw_apmb(bnv, bn2, r21Values, abcValues, tptValues, 1);
        }
      }
    }
  }

  if (set_draw2)
  {
    for (auto var2 : var2Array)
    {
      auto bnx = var2.bx();
      auto bny = var2.by();

      for (auto iComb : selSysCombIdx)
      {
        fiComb = iComb;
        fiSys2 = fSysCombIdx[fiComb][0];
        fiSys1 = fSysCombIdx[fiComb][1];

        auto graph_ab = new TGraphErrors();

        double tptValues[8][8][2][2] = {{0}};
        double r21Values[8][8][5][2] = {{0}};
        double abcValues[8][8][3][2] = {{0}};

        bnx.reset();
        while (bnx.next())
        {
          bny.reset();
          while (bny.next())
          {
            auto graph_r21 = new TGraph2DErrors();

            double yieldValues[5][2] = {0};

            for (auto iParticle : fParticleIdx)
            {
              fiParticle = iParticle;

              auto value2 = fHYValue[fiSys2][fiParticle][findv(var2)] -> GetBinContent(bnx.bi(),bny.bi());
              auto value1 = fHYValue[fiSys1][fiParticle][findv(var2)] -> GetBinContent(bnx.bi(),bny.bi());
              auto loose2 = fHYLoose[fiSys2][fiParticle][findv(var2)] -> GetBinContent(bnx.bi(),bny.bi());
              auto loose1 = fHYLoose[fiSys1][fiParticle][findv(var2)] -> GetBinContent(bnx.bi(),bny.bi());

              yieldValues[iParticle][0] = value2;
              yieldValues[iParticle][1] = value1;

              auto r21Value = value2/value1;
              auto r21Loose = loose2/loose1;
              auto r21Error = abs(r21Value - r21Loose);

              r21Values[bnx.ii()][bny.ii()][fiParticle][kVal] = r21Value;
              r21Values[bnx.ii()][bny.ii()][fiParticle][kErr] = r21Error;

              auto zValue = fParticleZ[fiParticle];
              auto nValue = fParticleN[fiParticle];

              auto idxR21 = graph_r21 -> GetN();
              graph_r21 -> SetPoint(idxR21, nValue, zValue, r21Value);
              graph_r21 -> SetPointError(idxR21, 0, 0, r21Value);
            }

            for (auto iss : {0,1}) {
              auto yieldd = yieldValues[kD][iss];
              auto yield4 = yieldValues[kHe4][iss];
              auto yieldt = yieldValues[kT][iss];
              auto yield3 = yieldValues[kHe3][iss];
              tptValues[bnx.ii()][bny.ii()][iss][kVal] = 14.3/TMath::Log(1.59*(yieldd * yield4 / yieldt / yield3));
            }

            auto fitIsoscaling = fit_r21(graph_r21);
            double alpha = fitIsoscaling -> GetParameter(0);
            double beta = fitIsoscaling -> GetParameter(1);
            double cnorm = fitIsoscaling -> GetParameter(2);


            abcValues[bnx.ii()][bny.ii()][0][kVal] = alpha;
            abcValues[bnx.ii()][bny.ii()][1][kVal] = beta;
            abcValues[bnx.ii()][bny.ii()][2][kVal] = cnorm;

            double alphaError = 0.;
            double betaError = 0.;

            for (auto iParticle : fParticleIdx)
            {
              fiParticle = iParticle;

              auto nValue = fParticleN[fiParticle];
              auto zValue = fParticleZ[fiParticle];

              auto r21Value = r21Values[bnx.ii()][bny.ii()][fiParticle][kVal];
              auto r21Eval = fitIsoscaling -> Eval(nValue, zValue);
              auto alphaError1 = (r21Value - r21Eval) / r21Value / nValue;
              auto betaError1 = (r21Value - r21Eval) / r21Value / zValue;

              if (nValue!=0) alphaError += alphaError1*alphaError1;
              if (zValue!=0) betaError += betaError1*betaError1;

            }
            alphaError = sqrt(alphaError / 4);
            betaError = sqrt(betaError / 5);
            abcValues[bnx.ii()][bny.ii()][0][kErr] = alphaError;
            abcValues[bnx.ii()][bny.ii()][1][kErr] = betaError;
          }
        }

        if (set_draw_ab_fit) draw_ab_fit(bnx, bny, r21Values, abcValues, tptValues);
        if (set_draw_r21_x) draw_r21_x(bnx, bny, r21Values, abcValues);
        if (set_draw_avb) draw_avb(bnx, bny, r21Values, abcValues);
      }
    }
  }

}


TH2D *make_r21_frame(TH2D *hist)
{
  hist -> GetYaxis() -> SetTitleOffset(0.5);

  hist -> GetXaxis() -> SetTitleSize(0.10);
  hist -> GetYaxis() -> SetTitleSize(0.10);
  hist -> GetXaxis() -> SetLabelSize(0.10);
  hist -> GetYaxis() -> SetLabelSize(0.10);
  hist -> GetXaxis() -> SetNdivisions(506);
  hist -> GetYaxis() -> SetNdivisions(506);
  hist -> GetXaxis() -> CenterTitle();
  //hist -> GetYaxis() -> CenterTitle();
  hist -> GetXaxis() -> SetNdivisions(4,0,0,true);

  return hist;
}

TF1 *get_fit1_r21(int nValue, int zValue, double alpha, double beta, double cnorm, double x1, double x2)
{
  TF1 *fit = nullptr;
  if (nValue>=0) {
    fit = new TF1("fit","[2]*exp([0]+[1]*x)",x1,x2);
    fit -> SetParameters(alpha*nValue,beta,cnorm);
    fit -> SetLineColor(kGray+1);
  }
  else if (zValue>=0) {
    fit = new TF1("fit","[2]*exp([0]*x+[1])",x1,x2);
    fit -> SetParameters(alpha,beta*zValue,cnorm);
    fit -> SetLineColor(kGray+1);
  }
  return fit;
}

void draw_ab_fit(binning bnx, binning bny, double r21Values[8][8][5][2], double abcValues[8][8][3][2], double tptValues[8][8][2][2])
{
  binning bn_n(100,-.8,2.8,"N");
  binning bn_z(100,0.2,2.8,"Z");
  //binning bn_r21(100,0.5,1.8,"R_{21}");
  binning bn_r21(100,0.3,2.0,"R_{21}");

  auto var2 = (bnx*bny);

  auto cvs_alpha = make_cvs(make_name("alpha",0,var2.namen()), 0, 100+fdxCvsAB*bnx.fN, fdyCvsAB*bny.fN, bnx.fN, bny.fN);
  auto cvs_beta  = make_cvs(make_name("beta",0,var2.namen()),  0, 100+fdxCvsAB*bnx.fN, fdyCvsAB*bny.fN, bnx.fN, bny.fN);

  bnx.reset();
  while (bnx.next())
  {
    bny.reset();
    while (bny.next())
    {
      TObjArray draw_alpha_fit(10);
      TObjArray draw_beta_fit(10);

      auto color_r21_error = kBlack;

      auto graph_alpha_error = new TGraphErrors();
      graph_alpha_error -> SetMarkerStyle(10);
      graph_alpha_error -> SetLineColor(color_r21_error);

      auto graph_beta_error  = new TGraphErrors();
      graph_beta_error -> SetMarkerStyle(10);
      graph_beta_error -> SetLineColor(color_r21_error);

      for (auto iParticle : fParticleIdx)
      {
        fiParticle = iParticle;

        auto nValue = fParticleN[iParticle];
        auto zValue = fParticleZ[iParticle];
        auto r21Value = r21Values[bnx.ii()][bny.ii()][fiParticle][0];
        auto r21Error = r21Values[bnx.ii()][bny.ii()][fiParticle][1];

        graph_alpha_error -> SetPoint(graph_alpha_error->GetN(),nValue,r21Value);
        graph_alpha_error -> SetPointError(graph_alpha_error->GetN()-1,0,r21Error);

        graph_beta_error -> SetPoint(graph_beta_error->GetN(),zValue,r21Value);
        graph_beta_error -> SetPointError(graph_beta_error->GetN()-1,0,r21Error);

        draw_alpha_fit.Add(attMarker(new TMarker(nValue,r21Value,20),fiParticle));
        draw_beta_fit.Add(attMarker(new TMarker(zValue,r21Value,20),fiParticle));
      }

      auto alpha = abcValues[bnx.ii()][bny.ii()][0][kVal];
      auto beta = abcValues[bnx.ii()][bny.ii()][1][kVal];
      auto cnorm = abcValues[bnx.ii()][bny.ii()][2][kVal];
      auto iii = var2.icvs(bnx.ii(),bny.ii());
      auto tpt2 = tptValues[bnx.ii()][bny.ii()][0][kVal];
      auto tpt1 = tptValues[bnx.ii()][bny.ii()][1][kVal];
      int dtpt = floor(0.5+100*abs(tpt2 - tpt1 / tpt1));

      auto cvs_a = cvs_alpha -> cd(iii);
      cvs_a -> SetLogy();
      auto frame_alpha = (bn_n*bn_r21).newHist(make_name("nr21",iii,var2.namen()));
      make_r21_frame(frame_alpha);
      frame_alpha -> Draw();
      graph_alpha_error -> Draw("same e");
      draw_alpha_fit.Draw("same");
      get_fit1_r21(-1,1,alpha,beta,cnorm,-0.2,2.2) -> Draw("samel");
      get_fit1_r21(-1,2,alpha,beta,cnorm, 0.8,2.2) -> Draw("samel");

      {
        auto mma = new TMarker(-.4,bn_r21.fMax-0.20*bn_r21.getFullWidth(),fDrawMStyle2[bnx.ii()]);
        mma -> SetMarkerColor(fDrawColor[bny.ii()]);
        mma -> SetMarkerSize(fDrawMSize2[bnx.ii()]);
        mma -> Draw();

        auto aa1 = new TLatex(-.1,bn_r21.fMax-0.20*bn_r21.getFullWidth(),Form("#alpha = %.2f",alpha));
        aa1 -> SetTextAlign(12);
        aa1 -> SetTextSize(0.1);
        aa1 -> SetTextFont(42);
        aa1 -> Draw();

        auto tt1 = new TLatex(-.6,bn_r21.fMin+0.1*bn_r21.getFullWidth(),Form("T_{%d} = %.1f",fSysCombBeam[fiComb][0],tpt2));
        tt1 -> SetTextAlign(11);
        tt1 -> SetTextSize(0.1);
        tt1 -> SetTextFont(42);
        tt1 -> Draw();

        auto tt2 = new TLatex(-.6,bn_r21.fMin+0.03*bn_r21.getFullWidth(),Form("T_{%d} = %.1f (%d%s)",fSysCombBeam[fiComb][1],tpt1,dtpt,"%"));
        tt2 -> SetTextAlign(11);
        tt2 -> SetTextSize(0.1);
        tt2 -> SetTextFont(42);
        tt2 -> Draw();

        auto line_r21_1 = new TLine(bn_n.fMin,1,bn_n.fMax,1);
        line_r21_1 -> SetLineColor(kGray);
        line_r21_1 -> SetLineStyle(2);
        line_r21_1 -> Draw("samel");
      }

      auto cvs_b = cvs_beta -> cd(iii);
      cvs_b -> SetLogy();
      make_r21_frame((bn_z*bn_r21).newHist(make_name("zr21",iii,var2.namen()))) -> Draw();
      graph_beta_error -> Draw("samee");
      draw_beta_fit.Draw("same");
      get_fit1_r21(0,-1,alpha,beta,cnorm,0.8,1.2) -> Draw("samel");
      get_fit1_r21(1,-1,alpha,beta,cnorm,0.8,2.2) -> Draw("samel");
      get_fit1_r21(2,-1,alpha,beta,cnorm,0.8,2.2) -> Draw("samel");

      {
        double offb = .8;

        auto mmb = new TMarker(offb+.5,bn_r21.fMax-0.20*bn_r21.getFullWidth(),fDrawMStyle2[bnx.ii()]);
        mmb -> SetMarkerColor(fDrawColor[bny.ii()]);
        mmb -> SetMarkerSize(fDrawMSize2[bnx.ii()]);
        mmb -> Draw();

        auto bb1 = new TLatex(offb+.7,bn_r21.fMax-0.20*bn_r21.getFullWidth(),Form("#beta = %.2f",beta));
        bb1 -> SetTextAlign(12);
        bb1 -> SetTextSize(0.1);
        bb1 -> SetTextFont(42);
        bb1 -> Draw();

        auto tt1 = new TLatex(0.4,bn_r21.fMin+0.1*bn_r21.getFullWidth(),Form("T_{%d} = %.1f",fSysCombBeam[fiComb][0],tpt2));
        tt1 -> SetTextAlign(11);
        tt1 -> SetTextSize(0.1);
        tt1 -> SetTextFont(42);
        tt1 -> Draw();

        auto tt2 = new TLatex(0.4,bn_r21.fMin+0.03*bn_r21.getFullWidth(),Form("T_{%d} = %.1f (%d%s)",fSysCombBeam[fiComb][1],tpt1,dtpt,"%"));
        tt2 -> SetTextAlign(11);
        tt2 -> SetTextSize(0.1);
        tt2 -> SetTextFont(42);
        tt2 -> Draw();

        auto line_r21_1 = new TLine(bn_z.fMin,1,bn_z.fMax,1);
        line_r21_1 -> SetLineColor(kGray);
        line_r21_1 -> SetLineStyle(2);
        line_r21_1 -> Draw("samel");
      }
    }
  }
}

void draw_r21_x(binning bnx, binning bny, double r21Values[8][8][5][2], double abcValues[8][8][3][2], bool do_for_y)
{
  binning bn_r21(100,0.5,1.8,"R_{21}");

  bny.reset();
  while (bny.next())
  {
    auto name = make_name("r21_x",bny.ii(),(bnx*bny).namen());
    auto title = make_title(bny,bny.bi());
    auto cvs = make_cvs(name);
    (bnx*bn_r21).newHist(name,title) -> Draw();

    auto line = new TLine(bnx.fMin,1,bnx.fMax,1);
    line -> SetLineColor(kGray+1);
    line -> SetLineStyle(2);
    line -> Draw("samle");

    auto legend = new TLegend();
    for (auto iParticle : fParticleIdx)
    {
      fiParticle = iParticle;

      auto graph_r21_x = new TGraphErrors();
      attGraph(graph_r21_x, fiParticle);
      legend -> AddEntry(graph_r21_x,fParticleNames[fiParticle],"pl");

      bnx.reset();
      while (bnx.next())
      {
        auto r21Value = r21Values[bnx.ii()][bny.ii()][fiParticle][0];
        graph_r21_x -> SetPoint(bnx.ii(),bnx.getValue(),r21Value);
      }

      graph_r21_x -> Draw("pl");
    }
    make_legend(cvs,legend) -> Draw();
  }

  if (!do_for_y) return;

  bnx.reset();
  while (bnx.next())
  {
    auto name = make_name("r21_y",bnx.ii(),(bny*bnx).namen());
    auto title = make_title(bnx,bnx.bi());
    auto cvs = make_cvs(name);
    (bny*bn_r21).newHist(name,title) -> Draw();

    auto line = new TLine(bny.fMin,1,bny.fMax,1);
    line -> SetLineColor(kGray+1);
    line -> SetLineStyle(2);
    line -> Draw("samle");

    auto legend = new TLegend();
    for (auto iParticle : fParticleIdx)
    {
      fiParticle = iParticle;

      auto graph_r21_x = new TGraphErrors();
      attGraph(graph_r21_x, fiParticle);
      legend -> AddEntry(graph_r21_x,fParticleNames[fiParticle],"pl");

      bny.reset();
      while (bny.next())
      {
        auto r21Value = r21Values[bnx.ii()][bny.ii()][fiParticle][0];
        graph_r21_x -> SetPoint(bny.ii(),bny.getValue(),r21Value);
      }

      graph_r21_x -> Draw("pl");
    }
    make_legend(cvs,legend) -> Draw();
  }
}

void draw_avb(binning bnx, binning bny, double r21Values[8][8][5][2], double abcValues[8][8][3][2])
{
  auto var2 = (bnx*bny);

  auto name = make_name("avb",0,var2.namen());
  auto cvs = make_cvs(name);
  binning bn_alpha(100,-.1,.5,"#alpha");
  binning bn_beta(100,-.5,.1,"#beta");
  (bn_alpha*bn_beta).newHist(name) -> Draw();

  vector<TMarker *> markers;
  vector<TGraphErrors *> graph_errorls;

  bnx.reset();
  while (bnx.next()) {
    TGraphErrors *graph_errorl = new TGraphErrors(); 
    graph_errorl -> SetLineColor(fDrawColor2[bnx.ii()]);
    graph_errorl -> SetMarkerColor(fDrawColor2[bnx.ii()]);
    graph_errorl -> SetMarkerSize(fDrawMSize2[bnx.ii()]);
    bny.reset();
    while (bny.next())
    {
      auto alpha = abcValues[bnx.ii()][bny.ii()][0][kVal];
      auto beta = abcValues[bnx.ii()][bny.ii()][1][kVal];
      auto alphaError = abcValues[bnx.ii()][bny.ii()][0][kErr];
      auto betaError = abcValues[bnx.ii()][bny.ii()][1][kErr];
      if (!bn_alpha.isInside(alpha)||!bn_beta.isInside(beta)) continue;

      graph_errorl -> SetPoint(graph_errorl->GetN(),alpha,beta);
      graph_errorl -> SetPointError(graph_errorl->GetN()-1,alphaError,betaError);

      auto marker = new TMarker(alpha,beta,20);
      marker -> SetMarkerColor(fDrawColor[bnx.ii()]);
      marker -> SetMarkerStyle(fDrawMStyle2[bny.ii()]);
      marker -> SetMarkerSize(fDrawMSize2[bny.ii()]);
      markers.push_back(marker);
    }

    graph_errorls.push_back(graph_errorl);
  }

  auto legend = new TLegend();
  bny.reset();
  while (bny.next()) {
    auto marker = new TMarker(0,0,20);
    marker -> SetMarkerStyle(fDrawMStyle2[bny.ii()]);
    marker -> SetMarkerSize(fDrawMSize2[bny.ii()]);
    //legend -> AddEntry(marker,Form("%s = %.0f ~ %.0f",bny.getTitle(),bny.cLow(),bny.cHigh()),"p");
    legend -> AddEntry(marker,Form("%s = %.1f ~ %.1f",bny.getTitle(),bny.cLow(),bny.cHigh()),"p");
  }
  bnx.reset();
  while (bnx.next()) {
    auto graph = new TGraph();
    graph -> SetLineColor(fDrawColor[bnx.ii()]);
    legend -> AddEntry(graph,Form("%s = %.0f ~ %.0f",bnx.getTitle(),bnx.cLow(),bnx.cHigh()),"l");
  }
  make_legend(cvs,legend,"",0,0,.5,.4) -> Draw();

  auto refxy0 = new TGraph;
  refxy0 -> SetPoint(0,bn_alpha.fMin,bn_beta.fMax);
  refxy0 -> SetPoint(1,bn_alpha.fMax,bn_beta.fMin);
  refxy0 -> SetLineColor(kGray);
  refxy0 -> SetLineStyle(2);
  refxy0 -> Draw("samel");

  auto refx0 = new TGraph();
  refx0 -> SetPoint(0,0,bn_beta.fMin);
  refx0 -> SetPoint(1,0,bn_beta.fMax);
  refx0 -> SetLineColor(kGray);
  refx0 -> SetLineStyle(2);
  refx0 -> Draw("samel");

  auto refy0 = new TGraph();
  refy0 -> SetPoint(0,bn_alpha.fMin,0);
  refy0 -> SetPoint(1,bn_alpha.fMax,0);
  refy0 -> SetLineColor(kGray);
  refy0 -> SetLineStyle(2);
  refy0 -> Draw("samel");

  for (auto graph : graph_errorls) graph -> Draw("sameezl");
  for (auto marker : markers) marker -> Draw("samep");
}

void draw_apmb(binning bnx, binning bny, double r21Values[8][8][5][2], double abcValues[8][8][3][2], double tptValues[8][8][2][2], bool multt)
{
  auto legend = new TLegend();

  auto bnab = binning(100,-.5,1);

  auto graph_alpha   = new TGraphErrors();
  graph_alpha -> SetMarkerStyle(24);
  graph_alpha -> SetMarkerSize(1.5);
  auto graph_beta    = new TGraphErrors();
  graph_beta -> SetMarkerStyle(25);
  graph_beta -> SetMarkerSize(1.5);
  auto graph_aplusb  = new TGraphErrors();
  graph_aplusb -> SetLineColor(kBlue);
  auto graph_aminusb = new TGraphErrors();
  graph_aminusb -> SetLineColor(kRed);

  auto line0 = new TLine(bnx.fMin,0,bnx.fMax,0);
  line0 -> SetLineStyle(2);
  line0 -> SetLineColor(kGray+1);
  auto line5 = new TLine(bnx.fMin,.5,bnx.fMax,.5);
  line5 -> SetLineStyle(2);
  line5 -> SetLineColor(kGray+1);

  if (multt) {
    bnab.set(100,-5,10);
    legend -> AddEntry(graph_alpha  ,"#mu_{n}",         "p");
    legend -> AddEntry(graph_beta   ,"#mu_{p}",         "p");
    legend -> AddEntry(graph_aplusb ,"#mu_{n}+#mu_{p}", "l");
    legend -> AddEntry(graph_aminusb,"#mu_{n}-#mu_{p}", "l");
    line5 -> SetY1(5);
    line5 -> SetY2(5);
  }
  else {
    legend -> AddEntry(graph_alpha  ,"#alpha",       "p");
    legend -> AddEntry(graph_beta   ,"#beta",        "p");
    legend -> AddEntry(graph_aplusb ,"#alpha+#beta", "l");
    legend -> AddEntry(graph_aminusb,"#alpha-#beta", "l");
  }

  bnx.reset();
  while (bnx.next())
  {
    auto alpha = abcValues[bnx.ii()][0][0][kVal];
    auto beta =  abcValues[bnx.ii()][0][1][kVal];

    if (multt) {
      auto temperature = (tptValues[bnx.ii()][0][0][kVal] + tptValues[bnx.ii()][0][1][kVal])/2.;
      alpha = alpha*temperature;
      beta = beta*temperature;
    }

    graph_alpha   -> SetPoint(bnx.ii(),bnx.getValue(),alpha);
    graph_beta    -> SetPoint(bnx.ii(),bnx.getValue(),beta);
    graph_aplusb  -> SetPoint(bnx.ii(),bnx.getValue(),alpha+beta);
    graph_aminusb -> SetPoint(bnx.ii(),bnx.getValue(),alpha-beta);
  }

  auto name_apmp = make_name(Form("apmp_%s",bnx.getE()));
  if (multt) name_apmp = make_name(Form("apmpt_%s",bnx.getE()));
  auto cvs = make_cvs(name_apmp);
  auto hist = (bnx*bnab).newHist(name_apmp);
  auto iSys2 = fSysCombIdx[fiComb][0];
  auto iSys1 = fSysCombIdx[fiComb][1];
  if (TString(bny.getTitle()).Index("y_{0}")>=0)
    hist -> SetTitle(Form("%s/%s : %s=%.2f~%.2f",fSysTitles[iSys2],fSysTitles[iSys1],bny.getTitle(),bny.fMin,bny.fMax));
  else
    hist -> SetTitle(Form("%s/%s : %s=%.0f~%.0f",fSysTitles[iSys2],fSysTitles[iSys1],bny.getTitle(),bny.fMin,bny.fMax));
  hist -> SetYTitle("#alpha, #beta");
  if (multt) hist -> SetYTitle("#mu");

  hist -> Draw();
  line0 -> Draw("samel");
  line5 -> Draw("samel");
  graph_alpha   -> Draw("samep");
  graph_beta    -> Draw("samep");
  graph_aplusb  -> Draw("samel");
  graph_aminusb -> Draw("samel");
  make_legend(cvs, legend) -> Draw();
}

TLegend *make_legend(TVirtualPad *cvs, TLegend *legend, TString fLegendDrawStyle="rt", double x_offset, double y_offset, double width_fixed, double height_fixed)
{
  double fWidthPerLengthLegend = 0.012;
  double fHeightPerEntryLegend = 0.05;
  double fWidthDefaultLegend = 0.15;
  int fFillStyleLegend = 0;
  int fBorderSizeLegend = 0;

  auto lmargin_cvs = fLMargin;
  auto rmargin_cvs = fRMargin;
  auto bmargin_cvs = fBMargin;
  auto tmargin_cvs = fTMargin;
  if (cvs != nullptr) {
    lmargin_cvs = cvs -> GetLeftMargin();
    rmargin_cvs = cvs -> GetRightMargin();
    bmargin_cvs = cvs -> GetBottomMargin();
    tmargin_cvs = cvs -> GetTopMargin();
  }
  else if (gPad != nullptr) {
    lmargin_cvs = ((TCanvas *) gPad) -> GetLeftMargin();
    rmargin_cvs = ((TCanvas *) gPad) -> GetRightMargin();
    bmargin_cvs = ((TCanvas *) gPad) -> GetBottomMargin();
    tmargin_cvs = ((TCanvas *) gPad) -> GetTopMargin();
  }

  auto x1_box =     lmargin_cvs;
  auto x2_box = 1.- rmargin_cvs;
  auto y1_box =     bmargin_cvs;
  auto y2_box = 1.- tmargin_cvs;
  auto unit_x = x2_box - x1_box;
  auto unit_y = y2_box - y1_box;

  auto y1_stat = y2_box;
  auto y2_stat = y2_box;
  auto x1_stat = x1_box;
  auto x2_stat = x2_box;
  auto dx_stat = x2_stat - x1_stat;
  auto dy_stat = y2_stat - y1_stat;

  if (fLegendDrawStyle=="xx1")
  {
    auto legend0 = legend;
    legend = new TLegend();
    auto count = 0;
    TIter next_entry(legend0->GetListOfPrimitives());
    while (auto label=(TLegendEntry*)next_entry()) {
      auto name = label -> GetLabel();
      auto obj = label -> GetObject();
      auto opt = label -> GetOption();
      if (count<5)
        legend -> AddEntry(obj,name,opt);
      count++;
    }
  }
  else if (fLegendDrawStyle=="xx2")
  {
    auto legend0 = legend;
    legend = new TLegend();
    auto count = 0;
    TIter next_entry(legend0->GetListOfPrimitives());
    while (auto label=(TLegendEntry*)next_entry()) {
      auto name = label -> GetLabel();
      auto obj = label -> GetObject();
      auto opt = label -> GetOption();
      if (count>=5)
        legend -> AddEntry(obj,name,opt);
      count++;
    }
  }

  auto length_text_max = 0;
  TIter next_entry(legend->GetListOfPrimitives());
  while (auto label=(TLegendEntry*)next_entry()) {
    auto lenth_text = TString(label->GetLabel()).Length();
    if (length_text_max<lenth_text)
      length_text_max = lenth_text;
  }

  auto y2_legend = y1_stat;
  auto y1_legend = y2_legend - legend->GetNRows() * fHeightPerEntryLegend;
  if (height_fixed>0) y1_legend = y2_legend - height_fixed*unit_y;

  auto x2_legend = x2_box;
  auto x1_legend = x2_legend - (length_text_max * fWidthPerLengthLegend) - fWidthDefaultLegend*unit_x;
  if (width_fixed>0) x1_legend = x2_legend - width_fixed*unit_x;

  if (fLegendDrawStyle=="rr")
  {
    x1_legend = x2_box;
    auto www = (length_text_max * fWidthPerLengthLegend) + fWidthDefaultLegend*unit_x;
    x2_legend = x1_legend + www;
    //cvs -> SetRightMargin(www);
  }
  if (fLegendDrawStyle=="lt")
  {
    x1_legend = x1_box;
    x2_legend = x1_legend + (length_text_max * fWidthPerLengthLegend) + fWidthDefaultLegend*unit_x;
    if (width_fixed>0) x2_legend = x1_legend + width_fixed*unit_x;
    if (x_offset==0) x_offset = 0.02;
  }
  else if (fLegendDrawStyle=="lb")
  {
    x1_legend = x1_box;
    x2_legend = x1_legend + (length_text_max * fWidthPerLengthLegend) + fWidthDefaultLegend*unit_x;
    if (width_fixed>0) x2_legend = x1_legend + width_fixed*unit_x;

    y1_legend = y1_box;
    y2_legend = y1_legend + legend->GetNRows() * fHeightPerEntryLegend;
    if (height_fixed>0) y1_legend = y2_legend - height_fixed*unit_y;
    if (x_offset==0) x_offset = 0.02;
    if (y_offset==0) y_offset = 0.04;
  }
  else if (fLegendDrawStyle=="fb")
  {
    x1_legend = x1_box;
    x2_legend = x2_box;

    y1_legend = y1_box;
    y2_legend = y1_legend + legend->GetNRows() * fHeightPerEntryLegend*1.2;
    if (height_fixed>0) y1_legend = y2_legend - height_fixed*unit_y;
    if (x_offset==0) x_offset = 0.02;
    if (y_offset==0) y_offset = 0.04;
  }

  else if (fLegendDrawStyle=="rb")
  {
    x2_legend = x2_box;
    x1_legend = x2_legend - (length_text_max * fWidthPerLengthLegend) - fWidthDefaultLegend*unit_x;
    if (width_fixed>0) x1_legend = x2_legend - width_fixed*unit_x;

    y1_legend = y1_box;
    y2_legend = y1_legend + legend->GetNRows() * fHeightPerEntryLegend;
    if (height_fixed>0) y1_legend = y2_legend - height_fixed*unit_y;
    if (y_offset==0) y_offset = 0.04;
  }
  else if (fLegendDrawStyle=="xx1")
  {
    x1_legend = x1_box;
    x2_legend = x1_legend + (length_text_max * fWidthPerLengthLegend) + fWidthDefaultLegend*unit_x;
    if (width_fixed>0) x2_legend = x1_legend + width_fixed*unit_x;

    y1_legend = y1_box;
    y2_legend = y1_legend + legend->GetNRows() * fHeightPerEntryLegend;
    if (height_fixed>0) y1_legend = y2_legend - height_fixed*unit_y;

    x_offset = 0.02;
    y_offset = 0.04;
  }
  else if (fLegendDrawStyle=="xx2")
  {
    x2_legend = x2_box;
    x1_legend = x2_legend - (length_text_max * fWidthPerLengthLegend) - fWidthDefaultLegend*unit_x;
    if (width_fixed>0) x1_legend = x2_legend - width_fixed*unit_x;

    y1_legend = y1_box;
    y2_legend = y1_legend + legend->GetNRows() * fHeightPerEntryLegend;
    if (height_fixed>0) y1_legend = y2_legend - height_fixed*unit_y;

    x_offset = 0;
    y_offset = 0.04;
  }

  x1_legend += x_offset*unit_x;
  x2_legend += x_offset*unit_x;
  y1_legend += y_offset*unit_y;
  y2_legend += y_offset*unit_y;

  auto dx_legend = x2_legend - x1_legend;

  if (dx_stat > 0) {
    if (abs(dx_legend-dx_stat)/dx_stat < 0.4)
      dx_legend = dx_stat;

    x1_legend = x2_legend - dx_legend;
  }

  legend -> SetX1(x1_legend);
  legend -> SetX2(x2_legend);
  legend -> SetY1(y1_legend);
  legend -> SetY2(y2_legend);

  legend -> SetFillStyle(fFillStyleLegend);
  legend -> SetBorderSize(fBorderSizeLegend);

  return legend;
}


const char *make_name(const char *name0, int idx, const char *tag) {
  const char *nameSC = ( fiComb>=0 ? Form("c%d_%s",fSysCombIdx2[fiComb],fSysCombNames[fiComb]) : Form("%d",fSysBeams[fiSys]) );
  const char *name = Form("%s_%s_%s_%s_%d",name0,tag,nameSC,fParticleNames[fiParticle],idx);
  return name;
}

const char *make_title(binning bn_range, int idx)
{
  const char *tt_range = "";
  if (!bn_range.isNull())
  {
    fiSys2 = fSysCombIdx[fiComb][0];
    fiSys1 = fSysCombIdx[fiComb][1];
    auto x1 = bn_range.fMin;
    auto x2 = bn_range.fMax;
    if (idx>=0) {
      x1 = bn_range.lowEdge(idx);
      x2 = bn_range.highEdge(idx);
    }
    tt_range = Form(" : %s=%.2f~%.2f",bn_range.getTitle(),x1,x2);
  }
  auto tt_sys = Form("%s/%s",fSysTitles[fiSys2],fSysTitles[fiSys1]);
  auto title = Form("%s%s",tt_sys,tt_range);
  return title;
}

void project(TTree *tree, TH1 *hist, const char *expression, const char *selection) {
  tree -> Project(hist->GetName(),expression,selection);
}

TF2 *fit_r21(TGraph2DErrors *graph_r21)
{
  TF2 *fitIsoscaling = new TF2("fitIsoscaling","[2]*exp([0]*x+[1]*y)",1,2,0,2);
  fitIsoscaling -> SetParameters(1,.3,-.3);
  graph_r21 -> Fit(fitIsoscaling,"RQ0");
  return fitIsoscaling;
}

TCanvas *make_cvs(const char *name, int opt, int w, int h, int nx, int ny)
{
  auto cvs = new TCanvas(name,name, fXOffCvs+20*(fCountCvs+1), 20*(fCountCvs+1), w, h);

  if (nx>0||ny>0) {
    if (opt==0)
    {
      auto sMargin = 0;
      double lMargin = 0.15;
      double rMargin = 0.055;
      double bMargin = 0.20;
      double tMargin = 0.12;
      cvs -> Divide(nx,ny,0,0);
      for (auto ix=1; ix<=nx; ++ix) {
        for (auto iy=1; iy<=ny; ++iy)
        {
          auto i = ix + nx*(iy-1);
          if (iy==1&&iy==ny) {
            if (ix==1) cvs->cd(i)->SetMargin(lMargin,0,bMargin,0);
            else if (ix==nx) cvs->cd(i)->SetMargin(0,rMargin,bMargin,0);
            else cvs->cd(i)->SetMargin(0,0,bMargin,0);
          } else if (iy==1) {
            if (ix==1&&ix==nx) cvs->cd(i)->SetMargin(lMargin,rMargin,sMargin,tMargin);
            else if (ix==1) cvs->cd(i)->SetMargin(lMargin,0,sMargin,tMargin);
            else if (ix==nx) cvs->cd(i)->SetMargin(0,rMargin,sMargin,tMargin);
            else cvs->cd(i)->SetMargin(0,0,sMargin,tMargin);
          } else if (iy==ny) {
            if (ix==1&&ix==nx) cvs->cd(i)->SetMargin(lMargin,rMargin,bMargin,sMargin);
            else if (ix==1) cvs->cd(i)->SetMargin(lMargin,0,bMargin,sMargin);
            else if (ix==nx) cvs->cd(i)->SetMargin(0,rMargin,bMargin,sMargin);
            else cvs->cd(i)->SetMargin(0,0,bMargin,sMargin);
          } else {
            if (ix==1&&ix==nx) cvs->cd(i)->SetMargin(lMargin,rMargin,sMargin,sMargin);
            else if (ix==1) cvs->cd(i)->SetMargin(lMargin,0,sMargin,sMargin);
            else if (ix==nx) cvs->cd(i)->SetMargin(0,rMargin,sMargin,sMargin);
            else cvs->cd(i)->SetMargin(0,0,sMargin,sMargin);
          }
        }
      }
    }
    else {
      cvs -> Divide(nx,ny);
    }
  }
  else if (opt==2) {
    cvs -> SetMargin(fLMargin,0.15,fBMargin,fTMargin);
  }
  else {
    cvs -> SetMargin(fLMargin,fRMargin,fBMargin,fTMargin);
  }

  fCountCvs++;

  fCvsArray.push_back(cvs);

  return cvs;
}

void save_cvs() {
  for (auto cvs : fCvsArray) {
    cout_info << fNameV << " " << cvs -> GetName() << endl;
    cvs -> cd();
    cvs -> SaveAs(Form("%s/figures/%s.png",fNameV.Data(),cvs->GetName()));
  }
}

TGraphErrors *attGraph(TGraphErrors *graph, int idx)
{
  graph -> SetMarkerStyle(fDrawMStyle[idx]);
  graph -> SetMarkerColor(fDrawColor[idx]);
  graph -> SetMarkerSize(fDrawMSize[idx]);
  graph -> SetLineColor(fDrawColor[idx]);
  return graph;
}

TMarker *attMarker(TMarker *marker, int idx)
{
  marker -> SetMarkerStyle(fDrawMStyle[idx]);
  marker -> SetMarkerColor(fDrawColor[idx]);
  marker -> SetMarkerSize(fDrawMSize[idx]);
  return marker;
}
