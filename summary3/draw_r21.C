#include "init_variables.h"
#include "binning.h"

int fiSys;
int fiSys2;
int fiSys1;
int fiComb = -1;
int fiParticle;
int fCountCvs = 0;
int fdxCvsAB = 150;
int fdyCvsAB = 200;

const char *make_name(const char *name0, int idx=0, const char *tag="");
void project(TTree *tree, TH1 *hist, const char *expression, const char *selection="");
TF2 *fit_r21(TGraph2DErrors *graph_r21);
TCanvas *make_cvs(const char *name, int opt=0, int w=630, int h=550, int nx=0, int ny=0);
TMarker *attMarker(TMarker *marker, int idx);
TH2D *make_r21_frame(TH2D *hist);
TF1 *get_fit1_r21(int nValue, int zValue, double alpha, double beta, double cnorm, double x1, double x2);
void draw_ab_fit(binning bnx, binning bny, double r21ValueArray[10][10][5], double r21ErrorArray[10][10][5], double abcValueArray[10][10][3]);
void draw_apmb(binning bnx, binning bny, double r21ValueArray[10][10][5], double r21ErrorArray[10][10][5], double abcValueArray[10][10][3]);

void draw_r21()
{
  bool set_recreate_hist = false;
  bool set_draw1 = true;
  bool set_draw2 = false;
  bool set_draw_ab_fit = false;
  bool set_draw_apmb = true;

  const char *anaV = "nn50";
  binning bn_y0(4,-.15,1.05,"y_{0}","y0");
  binning bn_pt(4,0,400,"p_{T}/A (MeV/c)","ptoac");
  binning bn_ke(4,0,120,"KE_{CM}/A (MeV)","keoac");
  binning bn_tt(4,0,100,"#theta_{CM} (deg)","ttc");

  vector<binning2> var2Array = {(bn_ke*bn_tt), (bn_pt*bn_y0)};

  auto findv = [&var2Array](binning2 var2) {
    for (auto idx=0; idx<var2Array.size(); ++idx) {
      if (strcmp(var2Array[idx].name(),var2.name())==0)
        return idx;
    }
    return -1;
  };


  binning bn_alpha(100,-.1,.5,"#alpha");
  binning bn_beta(100,-.5,.1,"#beta");


  TH2D *hyValue[4][5][2] = {{0}};
  TH2D *hyLoose[4][5][2] = {{0}};

  gStyle -> SetOptStat(0);

  for (auto iSys : fSysIdx)
  {
    fiSys = iSys;

    auto name_file_hist = Form("data2/%s/compact_sys%d_%s.hist_%d_%d_%d_%d.root",anaV,fSysBeams[fiSys],anaV,bn_y0.fN,bn_pt.fN,bn_ke.fN,bn_tt.fN);
    TFile *file_hist = new TFile(name_file_hist,"read");
    cout << name_file_hist << endl;

    if (set_recreate_hist || file_hist->IsZombie())
    {
      auto name_file_in = Form("data2/%s/compact_sys%d_%s.root",anaV,fSysBeams[fiSys],anaV);
      auto file_tree = new TFile(name_file_in);
      cout << name_file_in  << endl;
      file_hist = new TFile(name_file_hist,"recreate");

      for (auto iParticle : fParticleIdx)
      {
        fiParticle = iParticle;

        auto tree = (TTree *) file_tree -> Get(fParticleNames[fiParticle]);

        for (auto var2 : var2Array)
        {
          auto hist_vv0 = var2.newHist(make_name(var2.name(),0)); project(tree,hist_vv0,var2.ee(),"prob/eff*(prob>.7)");
          auto hist_vv1 = var2.newHist(make_name(var2.name(),1)); project(tree,hist_vv1,var2.ee(),"prob/eff*(prob>.5)");

          hyValue[fiSys][fiParticle][findv(var2)] = hist_vv0;
          hyLoose[fiSys][fiParticle][findv(var2)] = hist_vv1;

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

        for (auto var2 : var2Array)
        {
          hyValue[fiSys][fiParticle][findv(var2)] = (TH2D *) file_hist -> Get(make_name(var2.name(),0));
          hyLoose[fiSys][fiParticle][findv(var2)] = (TH2D *) file_hist -> Get(make_name(var2.name(),1));
        }
      }
    }
  }

  /////////////////////////////////////////////////////
  /////////////////////////////////////////////////////
  /////////////////////////////////////////////////////
  /////////////////////////////////////////////////////

  for (auto var2 : var2Array)
  {
    auto bnx = var2.bx();
    auto bny = var2.by();

    for (auto iComb : fSysCombIndx)
    {
      fiComb = iComb;
      fiSys2 = fSysCombIdx[fiComb][0];
      fiSys1 = fSysCombIdx[fiComb][1];

      if (set_draw1)
      {
        for (int ivar : {0,1})
        {
          auto bn_var = bnx;
          auto bn_sum = bny;
          if (ivar==1) {
            bn_var = bny;
            bn_sum = bnx;
          }

          double r21ValueArray[10][10][5] = {{0}};
          double r21ErrorArray[10][10][5] = {{0}};
          double abcValueArray[10][10][3] = {{0}};
          double abcErrorArray[10][10][3] = {{0}};

          bn_var.reset();
          while (bn_var.next())
          {
            auto graph_r21 = new TGraph2DErrors();

            for (auto iParticle : fParticleIdx)
            {
              fiParticle = iParticle;

              double value2 = 0;
              double value1 = 0;
              double loose2 = 0;
              double loose1 = 0;

              if (ivar==0) {
                bn_sum.reset(); while (bn_sum.next()) { value2 += hyValue[fiSys2][fiParticle][findv(var2)] -> GetBinContent(bn_var.bi(),bn_sum.bi()); }
                bn_sum.reset(); while (bn_sum.next()) { value1 += hyValue[fiSys1][fiParticle][findv(var2)] -> GetBinContent(bn_var.bi(),bn_sum.bi()); }
                bn_sum.reset(); while (bn_sum.next()) { loose2 += hyLoose[fiSys2][fiParticle][findv(var2)] -> GetBinContent(bn_var.bi(),bn_sum.bi()); }
                bn_sum.reset(); while (bn_sum.next()) { loose1 += hyLoose[fiSys1][fiParticle][findv(var2)] -> GetBinContent(bn_var.bi(),bn_sum.bi()); }
              }
              else {
                bn_sum.reset(); while (bn_sum.next()) { value2 += hyValue[fiSys2][fiParticle][findv(var2)] -> GetBinContent(bn_sum.bi(),bn_var.bi()); }
                bn_sum.reset(); while (bn_sum.next()) { value1 += hyValue[fiSys1][fiParticle][findv(var2)] -> GetBinContent(bn_sum.bi(),bn_var.bi()); }
                bn_sum.reset(); while (bn_sum.next()) { loose2 += hyLoose[fiSys2][fiParticle][findv(var2)] -> GetBinContent(bn_sum.bi(),bn_var.bi()); }
                bn_sum.reset(); while (bn_sum.next()) { loose1 += hyLoose[fiSys1][fiParticle][findv(var2)] -> GetBinContent(bn_sum.bi(),bn_var.bi()); }
              }

              auto r21Value = value2/value1;
              auto r21Loose = loose2/loose1;
              auto r21Error = abs(r21Value - r21Loose);

              r21ValueArray[bn_var.ii()][0][fiParticle] = r21Value;
              r21ErrorArray[bn_var.ii()][0][fiParticle] = r21Error;

              auto zValue = fParticleZ[fiParticle];
              auto nValue = fParticleN[fiParticle];

              auto idxR21 = graph_r21 -> GetN();
              graph_r21 -> SetPoint(idxR21, nValue, zValue, r21Value);
              graph_r21 -> SetPointError(idxR21, 0, 0, r21Value);
            }

            auto fitIsoscaling = fit_r21(graph_r21);
            double alpha = fitIsoscaling -> GetParameter(0);
            double beta = fitIsoscaling -> GetParameter(1);
            double cnorm = fitIsoscaling -> GetParameter(2);
            abcValueArray[bn_var.ii()][0][0] = alpha;
            abcValueArray[bn_var.ii()][0][1] = beta;
            abcValueArray[bn_var.ii()][0][2] = cnorm;

            double alphaError = 0.;
            double betaError = 0.;

            for (auto iParticle : fParticleIdx)
            {
              fiParticle = iParticle;

              auto nValue = fParticleN[iParticle];
              auto zValue = fParticleZ[iParticle];

              auto r21Value = r21ValueArray[bn_var.ii()][0][fiParticle];
              auto r21Eval = fitIsoscaling -> Eval(nValue, zValue);
              auto alphaError1 = (r21Value - r21Eval) / r21Value / nValue;
              auto betaError1 = (r21Value - r21Eval) / r21Value / zValue;

              if (nValue!=0) alphaError += alphaError1*alphaError1;
              if (zValue!=0) betaError += betaError1*betaError1;
            }
            alphaError = sqrt(alphaError / 4);
            betaError = sqrt(betaError / 5);
            abcErrorArray[bn_var.ii()][0][0] = alphaError;
            abcErrorArray[bn_var.ii()][0][1] = betaError;
          }

          auto bn2 = bn_sum;
          bn2.setN(1);
          bn2.setExpression("all");

          if (set_draw_ab_fit) draw_ab_fit(bn_var, bn2, r21ValueArray, r21ErrorArray, abcValueArray);
          if (set_draw_apmb) draw_apmb(bn_var, bn2, r21ValueArray, r21ErrorArray, abcValueArray);

        }
      }

      if (set_draw2)
      {
        auto graph_ab = new TGraphErrors();

        double r21ValueArray[10][10][5] = {{0}};
        double r21ErrorArray[10][10][5] = {{0}};
        double abcValueArray[10][10][3] = {{0}};
        double abcErrorArray[10][10][3] = {{0}};

        bnx.reset();
        while (bnx.next())
        {
          bny.reset();
          while (bny.next())
          {
            auto graph_r21 = new TGraph2DErrors();

            for (auto iParticle : fParticleIdx)
            {
              fiParticle = iParticle;

              auto value2 = hyValue[fiSys2][fiParticle][findv(var2)] -> GetBinContent(bnx.bi(),bny.bi());
              auto value1 = hyValue[fiSys1][fiParticle][findv(var2)] -> GetBinContent(bnx.bi(),bny.bi());
              auto loose2 = hyLoose[fiSys2][fiParticle][findv(var2)] -> GetBinContent(bnx.bi(),bny.bi());
              auto loose1 = hyLoose[fiSys1][fiParticle][findv(var2)] -> GetBinContent(bnx.bi(),bny.bi());

              auto r21Value = value2/value1;
              auto r21Loose = loose2/loose1;
              auto r21Error = abs(r21Value - r21Loose);

              r21ValueArray[bnx.ii()][bny.ii()][fiParticle] = r21Value;
              r21ErrorArray[bnx.ii()][bny.ii()][fiParticle] = r21Error;

              auto zValue = fParticleZ[fiParticle];
              auto nValue = fParticleN[fiParticle];

              auto idxR21 = graph_r21 -> GetN();
              graph_r21 -> SetPoint(idxR21, nValue, zValue, r21Value);
              graph_r21 -> SetPointError(idxR21, 0, 0, r21Value);
            }

            auto fitIsoscaling = fit_r21(graph_r21);
            double alpha = fitIsoscaling -> GetParameter(0);
            double beta = fitIsoscaling -> GetParameter(1);
            double cnorm = fitIsoscaling -> GetParameter(2);
            abcValueArray[bnx.ii()][bny.ii()][0] = alpha;
            abcValueArray[bnx.ii()][bny.ii()][1] = beta;
            abcValueArray[bnx.ii()][bny.ii()][2] = cnorm;

            double alphaError = 0.;
            double betaError = 0.;

            for (auto iParticle : fParticleIdx)
            {
              fiParticle = iParticle;

              auto nValue = fParticleN[iParticle];
              auto zValue = fParticleZ[iParticle];

              auto r21Value = r21ValueArray[bnx.ii()][bny.ii()][fiParticle];
              auto r21Eval = fitIsoscaling -> Eval(nValue, zValue);
              auto alphaError1 = (r21Value - r21Eval) / r21Value / nValue;
              auto betaError1 = (r21Value - r21Eval) / r21Value / zValue;

              if (nValue!=0) alphaError += alphaError1*alphaError1;
              if (zValue!=0) betaError += betaError1*betaError1;
            }
            alphaError = sqrt(alphaError / 4);
            betaError = sqrt(betaError / 5);
            abcErrorArray[bnx.ii()][bny.ii()][0] = alphaError;
            abcErrorArray[bnx.ii()][bny.ii()][1] = betaError;
          }
        }

        if (set_draw_ab_fit)
          draw_ab_fit(bnx, bny, r21ValueArray, r21ErrorArray, abcValueArray);

        bool draw_avb = true;
        if (draw_avb)
        {
          make_cvs(make_name("avb"));
          (bn_alpha*bn_beta).newHist(make_name("avb")) -> Draw();

          bnx.reset();
          while (bnx.next()) {
            TGraphErrors *graph_avb = new TGraphErrors(); 
            bny.reset();
            while (bny.next())
            {
              auto alpha = abcValueArray[bnx.ii()][bny.ii()][0];
              auto beta = abcValueArray[bnx.ii()][bny.ii()][1];
              auto cnorm = abcValueArray[bnx.ii()][bny.ii()][2];

              graph_avb -> SetPoint(graph_avb->GetN(),alpha,beta);
            }

            graph_avb -> Draw("same*l");
          }
        }

        return;
      }
    }
  }

}

const char *make_name(const char *name0, int idx, const char *tag) {
  const char *nameSC = ( fiComb>=0 ? fSysCombNames[fiComb] : Form("%d",fSysBeams[fiSys]) );
  const char *name = Form("%s_%s_%s_%s_%d",name0,tag,nameSC,fParticleNames[fiParticle],idx);
  return name;
}

void project(TTree *tree, TH1 *hist, const char *expression, const char *selection) {
  tree -> Project(hist->GetName(),expression,selection);
}

TF2 *fit_r21(TGraph2DErrors *graph_r21)
{
  TF2 *fitIsoscaling = new TF2("fitIsoscaling","[2]*exp([0]*x+[1]*y)",1,2,0,2);
  fitIsoscaling -> SetParameters(1,.3,-.3);
  graph_r21 -> Fit(fitIsoscaling,"","RQ0");
  return fitIsoscaling;
}

TCanvas *make_cvs(const char *name, int opt, int w, int h, int nx, int ny)
{
  auto cvs = new TCanvas(name,name, 20+20*(fCountCvs+1), 20*(fCountCvs+1), w, h);
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
  fCountCvs++;
  return cvs;
}

TMarker *attMarker(TMarker *marker, int idx)
{
  marker -> SetMarkerStyle(fDrawMStyle[idx]);
  marker -> SetMarkerColor(fDrawColor[idx]);
  marker -> SetMarkerSize(fDrawMSize[idx]);
  return marker;
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

void draw_ab_fit(binning bnx, binning bny, double r21ValueArray[10][10][5], double r21ErrorArray[10][10][5], double abcValueArray[10][10][3])
{
  binning bn_n(100,-.8,2.8,"N");
  binning bn_z(100,0.2,2.8,"Z");
  binning bn_r21(100,0.5,2.8,"R_{21}");

  auto var2 = (bnx*bny);

  auto cvs_alpha = make_cvs(make_name("alpha",0,var2.name()), 0, 100+fdxCvsAB*bnx.fN, fdyCvsAB*bny.fN, bnx.fN, bny.fN);
  auto cvs_beta  = make_cvs(make_name("beta",0,var2.name()),  0, 100+fdxCvsAB*bnx.fN, fdyCvsAB*bny.fN, bnx.fN, bny.fN);

  bnx.reset();
  while (bnx.next()) {
    bny.reset();
    while (bny.next())
    {
      TObjArray draw_alpha_fit(10);
      TObjArray draw_beta_fit(10);

      auto color_r21_error = kGray+1;

      auto graph_alpha_fit = new TGraphErrors();
      graph_alpha_fit -> SetMarkerStyle(10);
      graph_alpha_fit -> SetLineColor(color_r21_error);

      auto graph_beta_fit  = new TGraphErrors();
      graph_beta_fit -> SetMarkerStyle(10);
      graph_alpha_fit -> SetLineColor(color_r21_error);

      for (auto iParticle : fParticleIdx)
      {
        fiParticle = iParticle;

        auto nValue = fParticleN[iParticle];
        auto zValue = fParticleZ[iParticle];
        auto r21Value = r21ValueArray[bnx.ii()][bny.ii()][fiParticle];
        auto r21Error = r21ErrorArray[bnx.ii()][bny.ii()][fiParticle];

        graph_alpha_fit -> SetPoint(graph_alpha_fit->GetN(),nValue,r21Value);
        graph_alpha_fit -> SetPointError(graph_alpha_fit->GetN()-1,0,r21Error);

        graph_beta_fit -> SetPoint(graph_beta_fit->GetN(),zValue,r21Value);
        graph_beta_fit -> SetPointError(graph_beta_fit->GetN()-1,0,r21Error);

        draw_alpha_fit.Add(attMarker(new TMarker(nValue,r21Value,20),fiParticle));
        draw_beta_fit.Add(attMarker(new TMarker(zValue,r21Value,20),fiParticle));
      }

      auto alpha = abcValueArray[bnx.ii()][bny.ii()][0];
      auto beta = abcValueArray[bnx.ii()][bny.ii()][1];
      auto cnorm = abcValueArray[bnx.ii()][bny.ii()][2];
      auto iii = var2.icvs(bnx.ii(),bny.ii());

      auto cvs_a = cvs_alpha -> cd(iii);
      cvs_a -> SetLogy();
      auto frame_alpha = (bn_n*bn_r21).newHist(make_name("nr21",iii,var2.name()));
      make_r21_frame(frame_alpha);
      frame_alpha -> Draw();
      graph_alpha_fit -> Draw("samee");
      draw_alpha_fit.Draw("same");
      get_fit1_r21(-1,1,alpha,beta,cnorm,-0.2,2.2) -> Draw("samel");
      get_fit1_r21(-1,2,alpha,beta,cnorm,-0.2,2.2) -> Draw("samel");

      auto cvs_b = cvs_beta -> cd(iii);
      cvs_b -> SetLogy();
      make_r21_frame((bn_z*bn_r21).newHist(make_name("zr21",iii,var2.name()))) -> Draw();
      graph_beta_fit -> Draw("samee");
      draw_beta_fit.Draw("same");
      get_fit1_r21(0,-1,alpha,beta,cnorm,0.8,2.2) -> Draw("samel");
      get_fit1_r21(1,-1,alpha,beta,cnorm,0.8,2.2) -> Draw("samel");
      get_fit1_r21(2,-1,alpha,beta,cnorm,0.8,2.2) -> Draw("samel");
    }
  }
}

void draw_apmb(binning bnx, binning bny, double r21ValueArray[10][10][5], double r21ErrorArray[10][10][5], double abcValueArray[10][10][3])
{
  if (1) {
    auto graph_alpha   = new TGraph();
    graph_alpha -> SetMarkerStyle(24);
    graph_alpha -> SetMarkerSize(1.5);
    graph_alpha -> SetMarkerColor(kRed);
    auto graph_beta    = new TGraph();
    graph_beta -> SetMarkerStyle(25);
    graph_beta -> SetMarkerSize(1.5);
    graph_beta -> SetMarkerColor(kRed);
    auto graph_aplusb  = new TGraph();
    auto graph_aminusb = new TGraph();
    graph_aminusb -> SetLineStyle(9);

    bnx.reset();
    while (bnx.next())
    {
      auto alpha = abcValueArray[bnx.ii()][0][0];
      auto beta =  abcValueArray[bnx.ii()][0][1];
      auto cnorm = abcValueArray[bnx.ii()][0][2];

      graph_alpha   -> SetPoint(bnx.ii(),bnx.getValue(),alpha);
      graph_beta    -> SetPoint(bnx.ii(),bnx.getValue(),beta);
      graph_aplusb  -> SetPoint(bnx.ii(),bnx.getValue(),alpha+beta);
      graph_aminusb -> SetPoint(bnx.ii(),bnx.getValue(),alpha-beta);
    }

    make_cvs(make_name(Form("apmp_%s",bnx.getE())));
    (bnx*binning(100,-.5,1)).newHist(make_name(Form("apmp_%s",bnx.getE()))) -> Draw();
    graph_alpha   -> Draw("samep");
    graph_beta    -> Draw("samep");
    graph_aplusb  -> Draw("samel");
    graph_aminusb -> Draw("samel");
  }
}
