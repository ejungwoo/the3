
#include "init_variables.h"

using namespace ejungwoo;

TString fNameV;
int fiSys, fiSys2, fiSys1, fiComb = -1, fiParticle;
vector<TCanvas *> fCvsArray;
TH2D *fHYValue[4][5][5] = {0};
TH2D *fHYLoose[4][5][5] = {0};

TString fTag = "";
const int fTemperatureCut = 10;
binning fbn_r21(100,0.5,1.8,"R_{21}");
binning fbn_dndy0(100,0,20,"dN/d(y_{0})");
//binning fbn_dndy0(100,0,8,"dN/d(y_{0})");
binning fbn_dndpt(100,0,.1,"dN/d(p_{T}/A)");
binning fbn_dndke(100,0,.1,"dN/d(KE_{CM}/A)");
binning fbn_dndtt(100,15,1,"dN/d(#theta_{CM})");

// SMM
bool fSetSMM = 0;
int fIdxSMM = 3;
int fLStyleSMM = 1;
auto fNumSMM = 11;
double fSMMValues[20][8] = {};
TString fSMMTitles[20];

// AMD
bool fSetAMD = 0;
double fAMDY0[20] = {-1.9, -1.7, -1.5, -1.3, -1.1, -0.9, -0.7, -0.5, -0.3, -0.1, 0.1, 0.3, 0.5, 0.7, 0.9, 1.1, 1.3, 1.5, 1.7, 1.9,};
double fAMDValues[20][2][3] = {
  {{0.097,0.001,0.000,  }, {0.102,0.000,0.000,}},
  {{0.244,0.012,0.000,  }, {0.327,0.017,0.000,}},
  {{0.762,0.128,0.010,  }, {0.976,0.129,0.008,}},
  {{1.985,0.713,0.157,  }, {2.795,0.749,0.148,}},
  {{4.769,2.670,0.967,  }, {6.653,2.944,0.871,}},
  {{8.612,6.081,2.933,  }, {11.476,6.555,2.611,}},
  {{11.997,9.261,4.706, }, {14.960,9.215,4.028,}},
  {{13.796,10.790,5.519,}, {16.778,10.242,4.405,}},
  {{14.788,11.311,5.636,}, {17.345,10.274,4.206,}},
  {{15.029,11.432,5.671,}, {17.118,9.948,4.056,}},
  {{15.001,11.499,5.617,}, {17.127,9.933,3.888,}},
  {{15.190,11.772,5.824,}, {17.026,10.171,4.012,}},
  {{14.305,11.630,6.023,}, {16.429,9.862,4.234,}},
  {{12.368,9.853,5.485, }, {14.785,8.858,3.697,}},
  {{9.045,6.959,3.521,  }, {11.076,5.996,2.318,}},
  {{5.200,3.111,1.209,  }, {6.347,2.748,0.778,}},
  {{2.312,0.846,0.190,  }, {2.724,0.676,0.119,}},
  {{0.855,0.135,0.010,  }, {1.005,0.106,0.006,}},
  {{0.271,0.018,0.001,  }, {0.327,0.013,0.000,}},
  {{0.099,0.001,0.000,  }, {0.110,0.003,0.000,}},
};

const char *make_name(const char *name0, int idx=0, const char *tag="");
const char *make_title(binning *bn_range, int idx=-1);
const char *make_title(binning bn_range=binning(), int idx=-1);
const char *make_title1(binning *bn_range, int idx, int isys);
void project(TTree *tree, TH1 *hist, const char *expression, const char *selection="");
TF2 *fit_r21(TGraph2DErrors *graph_r21);
void save_cvs();
TGraphErrors *attGraph(TGraphErrors *graph, int idx);
TH2D *make_r21_frame(TH2D *hist);

TF1 *get_fit1_r21(int nValue, int zValue, double alpha, double beta, double cnorm, double x1, double x2);
void draw_ab_fit(binning bnx, binning bny, double r21Values[20][20][5][4], double abcValues[20][20][3][2], double tptValues[20][20][2][2], bool useTitle = 0, TString titles[20][20] = {{0}});
void draw_r21_x(binning bnx, binning bny, double r21Values[20][20][5][4], double abcValues[20][20][3][2], bool draw_yield, bool do_for_y=true, int is_x_smm_amd=0);
void draw_apmb(binning bnx, binning bny, double r21Values[20][20][5][4], double abcValues[20][20][3][2], double tptValues[20][20][2][2], bool multt=0);
void draw_avb(binning bnx, binning bny, double r21Values[20][20][5][4], double abcValues[20][20][3][2], bool drawLegend = true, int is_x_smm_amd = 0);
TLegend *make_legend(TVirtualPad *cvs, TLegend *legend, TString opt = "", double x_offset=0, double y_offset=0, double width_fixed=0, double height_fixed=0);

void draw_r21(int fnn=100, int recreate_hist=0, int hold_all=0, int add_all=0, int iversion=3/*13*/, bool paper_version=0)
{
  bool set_recreate_hist = recreate_hist;

  bool set_draw1 = kSet;
  bool set_draw2 = kUnSet;

  bool set_draw_range_box = kHold;
  bool set_draw_yield = kHold;
  bool set_draw_pncount = kHold;
  bool set_draw_r21_x = kDraw;
  bool set_draw_yield_x = kDraw;
  bool set_draw_apmb = kHold;
  bool set_draw_apmbt = kHold;
  bool set_draw_ab_fit = kHold;
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
  if (iversion==1) anaV = "rb3";
  if (iversion==2) anaV = "vb3";
  if (iversion==3) anaV = "rb3_tight";
  if (iversion==4) anaV = "rb3_mid";
  if (iversion==5) anaV = "rb3_loose";
  if (iversion==6) anaV = "af7";
  if (iversion==7) anaV = "f7";
  const char *compact_name = "compact";

  fTag = Form("nn%d",fnn);

  binning bn_pozl(200,0,2000,"p_{Lab.}/Z","pozl");
  binning bn_dedx(200,0,1500,"dE/dx","dedx");

  binning bn_y0(fnn,0,1,"y_{0}","y0");
  binning bn_pt(fnn,0,250,"p_{T}/A (MeV/c)","ptoac");
  binning bn_ke(fnn,0,120,"KE_{CM}/A (MeV)","keoac");
  binning bn_tt(fnn,15,95,"#theta_{CM} (deg)","ttc");

  if (fnn==4) { bn_tt.setMin(20); bn_tt.setMax(100); }
  if (fnn==0) {
    fTag = "zoom";
    bn_y0.set(8,0,1,"y_{0}","y0");
    bn_pt.set(5,0,200,"p_{T}/A (MeV/c)","ptoac");
    bn_ke.set(8,0,100,"KE_{CM}/A (MeV)","keoac");
    bn_tt.set(8,0,90,"#theta_{CM} (deg)","ttc");
  }
  if (fnn==100) {
    fTag = "full";
    bn_y0.set(24,-1,2,"y_{0}","y0");
    bn_pt.set(24,0,600,"p_{T}/A (MeV/c)","ptoac");
    bn_ke.set(24,0,200,"KE_{CM}/A (MeV)","keoac");
    bn_tt.set(24,0,180,"#theta_{CM} (deg)","ttc");
  }

  binning bn_y0_full(100,-1,2,"y_{0}","y0");
  binning bn_pt_full(100,0,600,"p_{T}/A (MeV/c)","ptoac");
  binning bn_ke_full(100,0,200,"KE_{CM}/A (MeV)","keoac");
  binning bn_tt_full(100,0,180,"#theta_{CM} (deg)","ttc");

  if (fSetSMM)
    fbn_r21.setMax(2.0);

  //vector<int> selSysIdx = {0,1,2,3};
  vector<int> selSysIdx = {0,1};

  //vector<int> selSysCombIdx = {0,1,2};
  vector<int> selSysCombIdx = {0};
  //vector<int> selSysCombIdx = {0,3};
  //vector<int> selSysCombIdx = {0,1,2,3,4,5};
  //vector<int> selSysCombIdx = {0,1,4,3,2,5};

  auto bnbnl_kt = bn_ke * bn_tt;
  auto bnbnl_py = bn_y0 * bn_pt;
  auto bnbnf_pe = bn_pozl * bn_dedx;
  auto bnbnf_kt = bn_ke_full * bn_tt_full;
  auto bnbnf_py = bn_y0_full * bn_pt_full;

  //vector<binning2> var2Array = {bnbnl_kt};
  vector<binning2> var2Array = {bnbnl_py};
  //vector<binning2> var2Array = {bnbnl_kt, bnbnl_py};

  if (add_all) {
    var2Array.clear();
    var2Array.push_back(bnbnl_kt);
    var2Array.push_back(bnbnl_py);
    selSysIdx.clear(); for (auto i : {0,1,2,3}) selSysIdx.push_back(i);
    selSysCombIdx.clear(); for (auto i : {0,1,2,3,4,5}) selSysCombIdx.push_back(i);
  }

  vector<binning2> var2ArrayPP;
  //for (auto var2 : var2Array) var2ArrayPP.push_back(var2);
  var2ArrayPP.push_back(bnbnl_kt);
  var2ArrayPP.push_back(bnbnl_py);
  var2ArrayPP.push_back(bnbnf_pe);
  var2ArrayPP.push_back(bnbnf_kt);
  var2ArrayPP.push_back(bnbnf_py);

  auto findv = [&var2ArrayPP](binning2 var2) {
    for (auto idx=0; idx<var2ArrayPP.size(); ++idx) {
      //cout << idx << " " << var2ArrayPP[idx].namexynn() << " " << var2.namexynn() << endl;
      if (strcmp(var2ArrayPP[idx].namexynn(),var2.namexynn())==0)
        return idx;
    }
    return -1;
  };

  fNameV = Form("%s_%s_%d_%d_%d_%d",compact_name,anaV,bn_y0.fN,bn_pt.fN,bn_ke.fN,bn_tt.fN);
  if (paper_version)
    fNameV = Form("%s_%s_%d_%d_%d_%d","paper_version",anaV,bn_y0.fN,bn_pt.fN,bn_ke.fN,bn_tt.fN);

  gStyle -> SetOptStat(0);

  if (fSetSMM)
  {
    double r21Values[20][20][5][4] = {{0}};
    double abcValues[20][20][3][2] = {{0}};
    double tptValues[20][20][2][2] = {{0}};
    TString smmTitles[20][20];

    for (auto idxSMM=0; idxSMM<fNumSMM; ++idxSMM)
    {
      TString fileName;
      TString smmTitle;
      int ix, iy;
           if (idxSMM==0)  { iy = 0; ix = 0; fileName = "data2/smm/srs_1.dat";         smmTitle = "[SMM]  T=8.7,  participant"; }
      else if (idxSMM==1)  { iy = 0; ix = 1; fileName = "data2/smm/srs_2.dat";         smmTitle = "[SMM]  T=8.7,  0.2#timespart + proj"; }
      else if (idxSMM==2)  { iy = 0; ix = 2; fileName = "data2/smm/srs_3.dat";         smmTitle = "[SMM]  T=8.7,  0.1#timespart + proj"; }
      else if (idxSMM==3)  { iy = 1; ix = 0; fileName = "data2/smm/t7_nz157_vv3.dat";  smmTitle = "[SMM]  T=7.0,  #frac{V_{bk}}{V_{0}} = 3,  #frac{N_{2}}{Z_{2}} = 1.57"; }
      else if (idxSMM==4)  { iy = 1; ix = 1; fileName = "data2/smm/t7_nz157_vv9.dat";  smmTitle = "[SMM]  T=7.0,  #frac{V_{bk}}{V_{0}} = 9,  #frac{N_{2}}{Z_{2}} = 1.57"; }
      else if (idxSMM==5)  { iy = 1; ix = 2; fileName = "data2/smm/t7_nz138_vv3.dat";  smmTitle = "[SMM]  T=7.0,  #frac{V_{bk}}{V_{0}} = 3,  #frac{N_{2}}{Z_{2}} = 1.38"; }
      else if (idxSMM==6)  { iy = 1; ix = 3; fileName = "data2/smm/t7_nz138_vv9.dat";  smmTitle = "[SMM]  T=7.0,  #frac{V_{bk}}{V_{0}} = 9,  #frac{N_{2}}{Z_{2}} = 1.38"; }
      else if (idxSMM==7)  { iy = 2; ix = 0; fileName = "data2/smm/t8_nz157_vv3.dat";  smmTitle = "[SMM]  T=8.0,  #frac{V_{bk}}{V_{0}} = 3,  #frac{N_{2}}{Z_{2}} = 1.57"; }
      else if (idxSMM==8)  { iy = 2; ix = 1; fileName = "data2/smm/t8_nz157_vv9.dat";  smmTitle = "[SMM]  T=8.0,  #frac{V_{bk}}{V_{0}} = 9,  #frac{N_{2}}{Z_{2}} = 1.57"; }
      else if (idxSMM==9)  { iy = 2; ix = 2; fileName = "data2/smm/t8_nz138_vv3.dat";  smmTitle = "[SMM]  T=8.0,  #frac{V_{bk}}{V_{0}} = 3,  #frac{N_{2}}{Z_{2}} = 1.38"; }
      else if (idxSMM==10) { iy = 2; ix = 3; fileName = "data2/smm/t8_nz138_vv9.dat";  smmTitle = "[SMM]  T=8.0,  #frac{V_{bk}}{V_{0}} = 9,  #frac{N_{2}}{Z_{2}} = 1.38"; }
      smmTitles[ix][iy] = smmTitle;
      fSMMTitles[idxSMM] = smmTitle;

      ifstream smmfile(fileName);
      double x0, r21Value;
      auto graph_r21 = new TGraph2DErrors();
      for (auto iParticle : fParticleIdx)
      {
        smmfile >> x0 >> r21Value;
        fSMMValues[idxSMM][iParticle] = r21Value;
        r21Values[ix][iy][iParticle][kVal] = r21Value;
        auto zValue = fParticleZ[iParticle];
        auto nValue = fParticleN[iParticle];
        auto idxR21 = graph_r21 -> GetN();
        graph_r21 -> SetPoint(idxR21, nValue, zValue, r21Value);
        graph_r21 -> SetPointError(idxR21, 0, 0, r21Value);
      }

      auto fitIsoscaling = fit_r21(graph_r21);
      double alpha = fitIsoscaling -> GetParameter(0);
      double beta = fitIsoscaling -> GetParameter(1);
      double cnorm = fitIsoscaling -> GetParameter(2);
      abcValues[ix][iy][0][kVal] = alpha;
      abcValues[ix][iy][1][kVal] = beta;
      abcValues[ix][iy][2][kVal] = cnorm;
    }

    auto bn1 = binning(4,0,0,"","smm");
    auto bn2 = binning(3,0,0,"","par");
    if (set_draw_ab_fit) draw_ab_fit(bn1, bn2, r21Values, abcValues, tptValues, 1, smmTitles);
    if (set_draw_avb) draw_avb(binning(), binning(), r21Values, abcValues, 0, 1);
  }

  for (auto iSys : selSysIdx)
  {
    fiSys = iSys;

    auto nameV = Form("%s.hist_%s_%s_%s_%s",anaV,bn_y0.nnmm(),bn_pt.nnmm(),bn_ke.nnmm(),bn_tt.nnmm());
    auto name_file_hist = Form("data2/%s/%s_sys%d_%s.root",anaV,compact_name,fSysBeams[fiSys],nameV);
    TFile *file_hist = nullptr;
    if (!set_recreate_hist) {
      file_hist = new TFile(name_file_hist,"read");
      cout_info << "hist file: " << name_file_hist << endl;
    }

    if (set_recreate_hist || file_hist->IsZombie())
    {
      auto name_file_in = Form("data2/%s/%s_sys%d_%s.root",anaV,compact_name,fSysBeams[fiSys],anaV);
      auto file_tree = new TFile(name_file_in);
      cout_info << "compact file: " << name_file_in  << endl;
      file_hist = new TFile(name_file_hist,"recreate");

      for (auto iParticle : fParticleIdx)
      {
        fiParticle = iParticle;

        auto tree = (TTree *) file_tree -> Get(fParticleNames[fiParticle]);

        for (auto var2 : var2ArrayPP)
        {
          auto hist_vv0 = var2.newHist(make_name(var2.namexynn(),0)); project(tree,hist_vv0,var2.expr(),"corr*(prob>.7)");
          auto hist_vv1 = var2.newHist(make_name(var2.namexynn(),1)); project(tree,hist_vv1,var2.expr(),"corr");

          fHYValue[fiSys][fiParticle][findv(var2)] = hist_vv0;
          fHYLoose[fiSys][fiParticle][findv(var2)] = hist_vv1;

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
          auto nameValue = make_name(var2.namexynn(),0);
          auto nameLoose = make_name(var2.namexynn(),1);
          fHYValue[fiSys][fiParticle][findv(var2)] = (TH2D *) file_hist -> Get(nameValue);
          fHYLoose[fiSys][fiParticle][findv(var2)] = (TH2D *) file_hist -> Get(nameLoose);
          auto histValue = fHYValue[fiSys][fiParticle][findv(var2)];
          auto histLoose = fHYLoose[fiSys][fiParticle][findv(var2)]; 
          //if (histValue!=nullptr) cout_info << fiSys << " " << fiParticle << " " << findv(var2) << " " << histValue->GetName() << endl;
          //else cout_error << "nullptr " << nameValue << endl;
          //if (histLoose!=nullptr) cout_info << fiSys << " " << fiParticle << " " << findv(var2) << " " << histLoose->GetName() << endl;
          //else cout_error << "nullptr " << nameLoose << endl;
        }
      }
    }
  }

  if (set_draw_range_box)
  {
    for (auto iSys : selSysIdx)
    {
      fiSys = iSys;
      for (auto ibnbn : {0,1})
      {
        auto bnbnf = bnbnf_kt;
        auto bnbnl = bnbnl_kt;
        if (ibnbn==1) {
          bnbnf = bnbnf_py;
          bnbnl = bnbnl_py;
        }
        auto bnx = bnbnl.bx();
        auto bny = bnbnl.by();

        double tptValues[20][20] = {{0}};
        auto histH2  = fHYValue[fiSys][kD][findv(bnbnl)];
        auto histH3  = fHYValue[fiSys][kT][findv(bnbnl)];
        auto histHe3 = fHYValue[fiSys][kHe3][findv(bnbnl)];
        auto histHe4 = fHYValue[fiSys][kHe4][findv(bnbnl)];

        bnx.reset();
        while (bnx.next()) {
          bny.reset();
          while (bny.next()) {
            auto yieldd = histH2 -> GetBinContent(bnx.bi(), bny.bi());
            auto yieldt = histH3 -> GetBinContent(bnx.bi(), bny.bi());
            auto yield3 = histHe3 -> GetBinContent(bnx.bi(), bny.bi());
            auto yield4 = histHe4 -> GetBinContent(bnx.bi(), bny.bi());
            tptValues[bnx.ii()][bny.ii()] = 14.3/TMath::Log(1.59*(yieldd * yield4 / yieldt / yield3));
          }
        }

        TString confName = "nn_z";
        int nx=3, ny=2;
        if (paper_version) {
          confName = "paper_nn";
          nx = 2;
          ny = 3;
        }
        auto cvs = canvas(make_name(bnbnf.namexy(),0),nx,ny,confName);
        for (auto iParticle : fParticleIdx)
        {
          fiParticle = iParticle;
          auto histFull = fHYValue[fiSys][fiParticle][findv(bnbnf)];
          //histFull -> SetTitle(Form("%s,  %s,  temperature (MeV)",fSysTitles2[fiSys],fParticleTitles[iParticle]));
          if (paper_version) {
            histFull -> SetTitle("");
            draw(histFull,cvs->cd(iParticle+1),"col");
          }
          else {
            histFull -> SetTitle(Form("%s,  %s",fSysTitles2[fiSys],fParticleTitles[iParticle]));
            draw(histFull,cvs->cd(iParticle+1),"colz");
          }

          bnbnl.rangeBoxGrid() -> Draw("samel");

          bnx.reset();
          while (bnx.next()) {
            bny.reset();
            while (bny.next()) {
              auto temperature = tptValues[bnx.ii()][bny.ii()];
              auto text = new TText(bnx.center(), bny.center(), Form("%.1f",temperature));
              text -> SetTextFont(133);
              text -> SetTextSize(15);
              text -> SetTextAlign(22);
              text -> Draw("same");

              if (temperature>fTemperatureCut) {
                auto line = bnbnl.rangeBox(bnx.bi(),bny.bi());
                line -> SetFillColor(kBlack);
                line -> SetFillStyle(3005);
                line -> Draw("samef");
              }
            }
          }

        }
      }
    }
  }

  if (set_draw_yield)
  {
    for (auto iSys : selSysIdx)
    {
      fiSys = iSys;

      auto numVar2 = var2ArrayPP.size();
      auto cvs_kk = canvas(make_name("kk"),5,numVar2,"single");

      for (auto var2 : var2ArrayPP)
      {
        auto addAll = false;
        //if (strcmp(var2.namexynn(),"dedxpozl")==0) addAll = true;

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
            cvs_kk -> cd((iParticle+1)+5*findv(var2));
            hist1 -> Draw("colz");
          }

          auto hist2 = fHYLoose[fiSys][fiParticle][findv(var2)];
          if (addAll) {
            if (histAll2 == nullptr) histAll2 = hist2;
            else histAll2 -> Add(hist2);
          }
          else {
            cvs_kk -> cd((iParticle+1)+5*findv(var2));
            hist2 -> Draw("colz");
          }
        }
        if (addAll) {
          canvas(histAll1->GetName()) -> SetLogz();
          histAll1 -> Draw("colz");

          canvas(histAll2->GetName()) -> SetLogz();
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

          double tptValues[20][20][2][2] = {{0}};
          double r21Values[20][20][5][4] = {{0}};
          double abcValues[20][20][3][2] = {{0}};

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
              r21Values[bnv.ii()][0][fiParticle][kSys2] = value2;
              r21Values[bnv.ii()][0][fiParticle][kSys1] = value1;

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
          if (set_draw_r21_x) draw_r21_x(bnv, bn2, r21Values, abcValues, set_draw_yield_x, 0);
          if (set_draw_apmb) draw_apmb(bnv, bn2, r21Values, abcValues, tptValues);
          if (set_draw_apmbt) draw_apmb(bnv, bn2, r21Values, abcValues, tptValues, 1);
          if (set_draw_avb) draw_avb(bnv, bn2, r21Values, abcValues);
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

        double tptValues[20][20][2][2] = {{0}};
        double r21Values[20][20][5][4] = {{0}};
        double abcValues[20][20][3][2] = {{0}};

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
              r21Values[bnx.ii()][bny.ii()][fiParticle][kSys2] = value2;
              r21Values[bnx.ii()][bny.ii()][fiParticle][kSys1] = value1;

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
        if (set_draw_r21_x) draw_r21_x(bnx, bny, r21Values, abcValues, set_draw_yield_x);
        if (set_draw_avb) draw_avb(bnx, bny, r21Values, abcValues);
      }
    }
  }

  //save_cvs();
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
    //fit -> SetLineColor(kGray+1);
    att(fit,nValue,"r21fit");
  }
  else if (zValue>=0) {
    fit = new TF1("fit","[2]*exp([0]*x+[1])",x1,x2);
    fit -> SetParameters(alpha,beta*zValue,cnorm);
    //fit -> SetLineColor(kGray+1);
    att(fit,zValue,"r21fit");
  }
  return fit;
}

void draw_ab_fit(binning bnx, binning bny, double r21Values[20][20][5][4], double abcValues[20][20][3][2], double tptValues[20][20][2][2], bool useTitle, TString titles[20][20])
{
  binning bn_n(100,-.8,2.8,"N");
  binning bn_z(100,0.2,2.8,"Z");

  bool draw_avb_marker = false;
  bool draw_two_ttt = false;
  bool legend_particle = true;

  auto var2 = (bnx*bny);

  auto cvs_alpha = canvas(make_name("alpha",0,var2.namexynn()), bnx.fN, bny.fN, "nn_s");
  auto cvs_beta  = canvas(make_name("beta",0,var2.namexynn()),  bnx.fN, bny.fN, "nn_s");

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

      double r21ValueMax = 0;
      double r21ValueMin = 10;

      if (bnx.ii()==0&&bny.ii()==0)
        legend_particle = true;
      else
        legend_particle = false;

      for (auto iParticle : fParticleIdx)
      {
        fiParticle = iParticle;

        auto nValue = fParticleN[iParticle];
        auto zValue = fParticleZ[iParticle];
        auto r21Value = r21Values[bnx.ii()][bny.ii()][fiParticle][0];
        auto r21Error = r21Values[bnx.ii()][bny.ii()][fiParticle][1];

        if (r21ValueMax<r21Value) r21ValueMax = r21Value;
        if (r21ValueMin>r21Value) r21ValueMin = r21Value;

        graph_alpha_error -> SetPoint(graph_alpha_error->GetN(),nValue,r21Value);
        graph_alpha_error -> SetPointError(graph_alpha_error->GetN()-1,0,r21Error);

        graph_beta_error -> SetPoint(graph_beta_error->GetN(),zValue,r21Value);
        graph_beta_error -> SetPointError(graph_beta_error->GetN()-1,0,r21Error);

        draw_alpha_fit.Add(att(new TMarker(nValue,r21Value,20),fiParticle));
        draw_beta_fit.Add(att(new TMarker(zValue,r21Value,20),fiParticle));
      }

      auto alpha = abcValues[bnx.ii()][bny.ii()][0][kVal];
      auto beta = abcValues[bnx.ii()][bny.ii()][1][kVal];
      auto cnorm = abcValues[bnx.ii()][bny.ii()][2][kVal];
      auto tpt2 = tptValues[bnx.ii()][bny.ii()][0][kVal];
      auto tpt1 = tptValues[bnx.ii()][bny.ii()][1][kVal];
      int dtpt = floor(0.5+100*abs(tpt2 - tpt1) / tpt1);

      auto iii = var2.icvs(bnx.ii(),bny.ii());
      auto cvs_a = cvs_alpha -> cd(iii);
      cvs_a -> SetLogy();
      auto frame_alpha = (bn_n*fbn_r21).newHist(make_name("nr21",iii,var2.namexynn()));
      if (bnx.getMin()!=bnx.getMax()) frame_alpha -> SetTitle(Form("#downarrow  %s=%s~%s",bnx.getTitle(), toString(bnx.low()).Data(), toString(bnx.high()).Data()));
      draw(frame_alpha,cvs_a,"z");
      if (bny.getMin()!=bny.getMax()) drawz(frame_alpha,cvs_a,Form("#downarrow  %s=%s~%s",bny.getTitle(), toString(bny.low()).Data(), toString(bny.high()).Data()));
      graph_alpha_error -> Draw("same e");
      draw_alpha_fit.Draw("same");
      get_fit1_r21(-1,1,alpha,beta,cnorm,-0.2,2.2) -> Draw("samel");
      get_fit1_r21(-1,2,alpha,beta,cnorm, 0.8,2.2) -> Draw("samel");

      {
        double x1tt1 = bn_n.xByRatio(.08);
        double y1tt2 = fbn_r21.xByRatio(.1*1./3);
        double y1tt1 = fbn_r21.xByRatio(.1);

        double x1mma = bn_n.xByRatio(.10);
        double y1mma = fbn_r21.xByRatio(.87);
        double x1aa1 = bn_n.xByRatio(.17);
        double y1aa1 = y1mma;
        if (!draw_avb_marker) x1aa1 = x1mma;

        if (draw_avb_marker) {
          auto mma = new TMarker(x1mma,y1mma,fDrawMStyle2[bnx.ii()]);
          mma -> SetMarkerColor(fDrawColor[bny.ii()]);
          mma -> SetMarkerSize(fDrawMSize2[bnx.ii()]);
          mma -> Draw();
        }

        auto aa1 = new TLatex(x1aa1,y1aa1,Form("#alpha = %.2f",alpha));
        aa1 -> SetTextAlign(12);
        aa1 -> SetTextSize(23);
        aa1 -> SetTextFont(133);
        aa1 -> Draw();

        if (draw_two_ttt) {
          auto tt1 = new TLatex(x1tt1,y1tt1,Form("T_{%d} = %.1f",fSysCombBeam[fiComb][0],tpt2));
          tt1 -> SetTextAlign(11);
          tt1 -> SetTextSize(23);
          tt1 -> SetTextFont(133);
          tt1 -> Draw();

          auto tt2 = new TLatex(x1tt1,y1tt2,Form("T_{%d} = %.1f (%d%s)",fSysCombBeam[fiComb][1],tpt1,dtpt,"%"));
          tt2 -> SetTextAlign(11);
          tt2 -> SetTextSize(23);
          tt2 -> SetTextFont(133);
          tt2 -> Draw();
        }
        else {
          auto tt1 = new TLatex(x1tt1,y1tt2,Form("T_{%d,%d} = %.1f, %.1f",fSysCombBeam[fiComb][0],fSysCombBeam[fiComb][1],tpt2,tpt1));
          tt1 -> SetTextAlign(11);
          tt1 -> SetTextSize(23);
          tt1 -> SetTextFont(133);
          if (useTitle) {
            TString inTitle = titles[bnx.ii()][bny.ii()];
            inTitle.ReplaceAll("  "," ");
            inTitle.ReplaceAll("[SMM]","");
            tt1 -> SetTitle(inTitle);
            tt1 -> SetTextSize(18);
          }
          if (tpt1>fTemperatureCut||tpt2>fTemperatureCut)
            tt1 -> SetTextColor(kRed);
          tt1 -> Draw();
        }

        auto line_r21_1 = new TLine(bn_n.fMin,1,bn_n.fMax,1);
        line_r21_1 -> SetLineColor(kGray);
        line_r21_1 -> SetLineStyle(2);
        line_r21_1 -> Draw("samel");
      }

      auto cvs_b = cvs_beta -> cd(iii);
      cvs_b -> SetLogy();
      auto frame_beta = (bn_z*fbn_r21).newHist(make_name("zr21",iii,var2.namexynn()));
      if (bny.getMin()!=bny.getMax()) frame_beta -> SetTitle(Form("#downarrow  %s=%s~%s",bnx.getTitle(), toString(bnx.low()).Data(), toString(bnx.high()).Data()));
      draw(frame_beta,cvs_b);
      if (bny.getMin()!=bny.getMax()) drawz(frame_beta,cvs_b,Form("#downarrow  %s=%s~%s",bny.getTitle(), toString(bny.low()).Data(), toString(bny.high()).Data()));
      graph_beta_error -> Draw("samee");
      draw_beta_fit.Draw("same");
      get_fit1_r21(0,-1,alpha,beta,cnorm,0.8,1.2) -> Draw("samel");
      get_fit1_r21(1,-1,alpha,beta,cnorm,0.8,2.2) -> Draw("samel");
      get_fit1_r21(2,-1,alpha,beta,cnorm,0.8,2.2) -> Draw("samel");

      if (legend_particle)
        for (auto iParticle : fParticleIdx)
        {
          fiParticle = iParticle;

          auto zValue = fParticleZ[iParticle];
          auto r21Value = r21Values[bnx.ii()][bny.ii()][fiParticle][0];

          double xOffset = 0.5;
          if (iParticle<3) xOffset = -0.45;
          auto ttp = new TLatex(zValue+xOffset,r21Value,Form("%s",fParticleTitles[iParticle]));
          ttp -> SetTextAlign(22);
          ttp -> SetTextSize(19);
          ttp -> SetTextFont(133);
          ttp -> SetTextColor(kBlack);
          ttp -> Draw();
        }

      {
        double x1tt1 = bn_z.xByRatio(.08);
        double y1tt2 = fbn_r21.xByRatio(.1*1./3);
        double y1tt1 = fbn_r21.xByRatio(.1);

        double x1mmb = bn_z.xByRatio(.45); //offb+.5,
        double y1mmb = fbn_r21.xByRatio(.87);//fbn_r21.fMax-0.13*fbn_r21.getFullWidth();
        double x1bb1 = bn_z.xByRatio(.52);
        double y1bb1 = y1mmb;
        //if (!draw_avb_marker) x1bb1 = x1mmb;

        if (draw_avb_marker) {
          auto mmb = new TMarker(x1mmb,y1mmb,fDrawMStyle2[bnx.ii()]);
          mmb -> SetMarkerColor(fDrawColor[bny.ii()]);
          mmb -> SetMarkerSize(fDrawMSize2[bnx.ii()]);
          mmb -> Draw();
        }

        auto bb1 = new TLatex(x1bb1,y1bb1,Form("#beta = %.2f",beta));
        bb1 -> SetTextAlign(12);
        bb1 -> SetTextSize(23);
        bb1 -> SetTextFont(133);
        bb1 -> Draw();

        if (draw_two_ttt) {
          auto tt1 = new TLatex(x1tt1,y1tt1,Form("T_{%d} = %.1f",fSysCombBeam[fiComb][0],tpt2));
          tt1 -> SetTextAlign(11);
          tt1 -> SetTextSize(23);
          tt1 -> SetTextFont(133);
          tt1 -> Draw();

          auto tt2 = new TLatex(x1tt1,y1tt2,Form("T_{%d} = %.1f (%d%s)",fSysCombBeam[fiComb][1],tpt1,dtpt,"%"));
          tt2 -> SetTextAlign(11);
          tt2 -> SetTextSize(23);
          tt2 -> SetTextFont(133);
          tt2 -> Draw();
        }
        else {
          auto tt1 = new TLatex(x1tt1,y1tt2,Form("T_{%d,%d} = %.1f, %.1f",fSysCombBeam[fiComb][0],fSysCombBeam[fiComb][1],tpt2,tpt1));
          tt1 -> SetTextAlign(11);
          tt1 -> SetTextSize(23);
          tt1 -> SetTextFont(133);
          if (useTitle) {
            TString inTitle = titles[bnx.ii()][bny.ii()];
            inTitle.ReplaceAll("  "," ");
            inTitle.ReplaceAll("[SMM]","");
            tt1 -> SetTitle(inTitle);
            tt1 -> SetTextSize(18);
          }
          if (tpt1>fTemperatureCut||tpt2>fTemperatureCut)
            tt1 -> SetTextColor(kRed);
          tt1 -> Draw();
        }

        auto line_r21_1 = new TLine(bn_z.fMin,1,bn_z.fMax,1);
        line_r21_1 -> SetLineColor(kGray);
        line_r21_1 -> SetLineStyle(2);
        line_r21_1 -> Draw("samel");
      }
    }
  }
}

void draw_r21_x(binning bnx, binning bny, double r21Values[20][20][5][4], double abcValues[20][20][3][2], bool draw_yield, bool do_for_y, int is_x_smm_amd)
{
  bool printYR21 = true;

  double spaceForLegend = 0.1;
  double lowBoundLegend = fbn_r21.fMax - spaceForLegend;
  double marginMarker = 0.5;
  TString headerName = "r21";
  TString headerName2 = "yield_s2";
  TString headerName1 = "yield_s1";
  if (fSetSMM)
    headerName = Form("r21_smm_%d",fIdxSMM);

  bool setAMD = false;
  if (fSetAMD && TString(bnx.getExpression())=="y0")
    setAMD = true;

  for (auto ibnx : {0,1})
  {
    binning *bnvar = &bnx;
    binning *bnfix = &bny;
    if (ibnx==0) {}
    else if (ibnx==1) {
      if (!do_for_y)
        return;
      bnvar = &bny;
      bnfix = &bnx;
    } else
      return;

    binning bn_dndx;
    if (TString(bnvar->getExpression())=="y0")    bn_dndx = fbn_dndy0;
    if (TString(bnvar->getExpression())=="ptoac") bn_dndx = fbn_dndpt;
    if (TString(bnvar->getExpression())=="keoac") bn_dndx = fbn_dndke;
    if (TString(bnvar->getExpression())=="ttc")   bn_dndx = fbn_dndtt;
    
    int cnx = 1, cny = 1;
    TString cconf = "";
    if (bnfix->fN==1)  { cnx=1; cny=1; cconf="single"; }
    else if (bnfix->fN==2)  { cnx=2; cny=1; cconf="nn"; }
    else if (bnfix->fN==4)  { cnx=2; cny=2; cconf="nn"; }
    else if (bnfix->fN==8)  { cnx=4; cny=2; cconf="nn"; }
    else if (bnfix->fN==5)  { cnx=3; cny=2; cconf="nn"; }
    else if (bnfix->fN==10) { cnx=5; cny=2; cconf="nn"; }

    auto name = make_name(Form("%s_all",headerName.Data()),0,(bnfix->mult(bnvar)).namexynn());
    TCanvas *cvsAll = canvas(name,cnx,cny,cconf);

    TCanvas *cvsAll2 = nullptr;
    TCanvas *cvsAll1 = nullptr;
    if (draw_yield)  {
      auto names2 = make_name(Form("%s_ysys2",headerName.Data()),0,(bnfix->mult(bnvar)).namexynn());
      auto names1 = make_name(Form("%s_ysys1",headerName.Data()),0,(bnfix->mult(bnvar)).namexynn());
      cvsAll2 = canvas(names2,cnx,cny,cconf);
      cvsAll1 = canvas(names1,cnx,cny,cconf);
    }

    bnfix->reset();
    while (bnfix->next())
    {
      if (printYR21) {
        auto x1 = bnfix->fMin;
        auto x2 = bnfix->fMax;
        if (bnfix->bi()>=0) {
          x1 = bnfix->lowEdge(bnfix->bi());
          x2 = bnfix->highEdge(bnfix->bi());
        }
        const char *x1String = toString(x1);
        const char *x2String = toString(x2);
        cout << endl;
        cout << "========================" << endl;
        auto tt2 = Form("(%s)/(%s) : %s=%s~%s",fSysTitles[fiSys2],fSysTitles[fiSys1], bnfix->getTitle(),x1String,x2String);
        cout << tt2 << endl;
      }

      TPad *cvs_r21 = (TPad *) cvsAll;
      if (bnfix->fN>1)
        cvs_r21 = (TPad *) cvsAll -> cd(bnfix->bi());
      auto name = make_name(Form("%s",headerName.Data()),bnfix->ii(),(bnfix->mult(bnvar)).namexynn());
      auto title = make_title(bnfix,bnfix->bi());
      auto hist = (bnvar->mult(&fbn_r21)).newHist(name,title);
      draw(hist,cvs_r21);

      TPad *cvs2 = nullptr;
      TPad *cvs1 = nullptr;
      if (draw_yield)
      {
        cvs2 = (TPad *) cvsAll2;
        if (bnfix->fN>1)
          cvs2 = (TPad *) cvsAll2 -> cd(bnfix->bi());
        auto name2 = make_name(Form("%s",headerName2.Data()),bnfix->ii(),(bnfix->mult(bnvar)).namexynn());
        auto title2 = make_title1(bnfix,bnfix->bi(),fiSys2);
        auto hist2 = (bnvar->mult(&bn_dndx)).newHist(name2,title2);
        draw(hist2,cvs2);

        cvs1 = (TPad *) cvsAll1;
        if (bnfix->fN>1)
          cvs1 = (TPad *) cvsAll1 -> cd(bnfix->bi());
        auto name1 = make_name(Form("%s",headerName1.Data()),bnfix->ii(),(bnfix->mult(bnvar)).namexynn());
        auto title1 = make_title1(bnfix,bnfix->bi(),fiSys1);
        auto hist1 = (bnvar->mult(&bn_dndx)).newHist(name1,title1);
        draw(hist1,cvs1);
      }

      if (!fSetSMM && !setAMD) {
        auto line1 = new TLine(bnvar->fMin,1,bnvar->fMax,1);
        line1 -> SetLineColor(kGray+1);
        line1 -> SetLineStyle(2);
        cvs_r21 -> cd();
        line1 -> Draw("samel");
      }

      auto lineTop = new TLine(bnvar->fMin,lowBoundLegend,bnvar->fMax,lowBoundLegend);
      lineTop -> SetLineColor(kBlack);
      cvs_r21 -> cd(); lineTop -> Draw("samel");

      if (draw_yield) {
        lowBoundLegend = bn_dndx.xByRatio(0.92);
        lineTop = new TLine(bnvar->fMin,lowBoundLegend,bnvar->fMax,lowBoundLegend);
        cvs2 -> cd();
        lineTop -> Draw("samel");

        lineTop = new TLine(bnvar->fMin,lowBoundLegend,bnvar->fMax,lowBoundLegend);
        cvs1 -> cd();
        lineTop -> Draw("samel");

      }

      if (fSetSMM) {
        for (auto iParticle : fParticleIdx) {
          double calValue = fSMMValues[fIdxSMM][iParticle];
          auto lineSMM = new TGraphErrors();
          lineSMM -> SetPoint(0,bnvar->fMin,calValue);
          lineSMM -> SetPointError(0,0,.01);
          lineSMM -> SetPoint(1,bnvar->fMax,calValue);
          lineSMM -> SetPointError(1,0,.01);
          attGraph(lineSMM,iParticle);
          lineSMM -> SetLineStyle(fLStyleSMM);
          cvs_r21 -> cd();
          lineSMM -> Draw("same l ");
        }
      }

      bnvar->reset(); cout << bnvar->getExpression() << ": "; while (bnvar->next()) { cout << bnvar->val() << ", ";  } cout << endl;

      std::vector<TGraphErrors*> graphsr;
      std::vector<TGraphErrors*> graphs2;
      std::vector<TGraphErrors*> graphs1;

      for (auto iParticle : fParticleIdx)
      {
        fiParticle = iParticle;

        auto graph_r21_x = new TGraphErrors();
        auto graph_s2 = new TGraphErrors(); graph_s2 -> SetName(Form("dndx_%s",fParticleNames[iParticle]));
        auto graph_s1 = new TGraphErrors(); graph_s1 -> SetName(Form("dndx_%s",fParticleNames[iParticle]));
        graphsr.push_back(graph_r21_x);
        graphs2.push_back(graph_s2);
        graphs1.push_back(graph_s1);

        for (auto graph : {graph_r21_x, graph_s2, graph_s1}) {
          attGraph(graph, fiParticle);
          graph -> SetLineStyle(1);
        }

        if (printYR21) {
          cout << "dndx_" << fSysBeams[fiSys2] << "_" << fParticleNames[iParticle] << ": "; bnvar->reset(); while (bnvar->next()) { cout << r21Values[bnx.ii()][bny.ii()][fiParticle][kSys2] << ", "; } cout << endl;
          cout << "dndx_" << fSysBeams[fiSys1] << "_" << fParticleNames[iParticle] << ": "; bnvar->reset(); while (bnvar->next()) { cout << r21Values[bnx.ii()][bny.ii()][fiParticle][kSys1] << ", "; } cout << endl;
          cout << "yield_" << fSysBeams[fiSys2] << "_" << fParticleNames[iParticle] << ": "; bnvar->reset(); while (bnvar->next()) { cout << r21Values[bnx.ii()][bny.ii()][fiParticle][kSys2]/bnvar->getW() << ", "; } cout << endl;
          cout << "yield_" << fSysBeams[fiSys1] << "_" << fParticleNames[iParticle] << ": "; bnvar->reset(); while (bnvar->next()) { cout << r21Values[bnx.ii()][bny.ii()][fiParticle][kSys1]/bnvar->getW() << ", "; } cout << endl;
          cout << "r21_" << fParticleNames[iParticle] << ": "; bnvar->reset(); while (bnvar->next()) { cout << r21Values[bnx.ii()][bny.ii()][fiParticle][kVal] << ", "; } cout << endl;
        }

        bnvar->reset();
        while (bnvar->next())
        {
          auto r21Value = r21Values[bnx.ii()][bny.ii()][fiParticle][kVal];
          auto r21Error = r21Values[bnx.ii()][bny.ii()][fiParticle][kErr];
          if (fbn_r21.isInside(r21Value)) {
            graph_r21_x -> SetPoint(graph_r21_x->GetN(),bnvar->val(),r21Value);
            graph_r21_x -> SetPointError(graph_r21_x->GetN()-1,0,r21Error);
          }

          graph_s2 -> SetPoint(bnvar->ii(),bnvar->val(), r21Values[bnx.ii()][bny.ii()][fiParticle][kSys2]/bnvar->getW());
          graph_s1 -> SetPoint(bnvar->ii(),bnvar->val(), r21Values[bnx.ii()][bny.ii()][fiParticle][kSys1]/bnvar->getW());
        }

        cvs_r21 -> cd();
        if (fSetSMM || setAMD)
          graph_r21_x -> Draw("p");
        else
          graph_r21_x -> Draw("pl");

        if (draw_yield)
        {
          cvs2 -> cd();
          graph_s2 -> Draw("pl");
          cvs1 -> cd();
          graph_s1 -> Draw("pl");
        }

        if (fSetSMM) {
          auto lineSMM = new TGraphErrors();
          lineSMM -> SetPoint(0,1,1);
          lineSMM -> SetPointError(0,0,.01);
          lineSMM -> SetLineColor(fDrawColor[iParticle]);
          lineSMM -> SetLineStyle(fLStyleSMM);
          auto legendSMM = new TLegend();
          legendSMM -> AddEntry(lineSMM,"","l");
          double sfl = spaceForLegend/fbn_r21.getFullWidth();
          draw(legendSMM,cvs_r21,iParticle*(1./fNumParticles),1.-sfl,1./fNumParticles,sfl,marginMarker);
          if (draw_yield) {
            draw(legendSMM,cvs2,iParticle*(1./fNumParticles),1.-sfl,1./fNumParticles,sfl,marginMarker);
            draw(legendSMM,cvs1,iParticle*(1./fNumParticles),1.-sfl,1./fNumParticles,sfl,marginMarker);
          }
        }
        auto legendSingle = new TLegend();
        if (fSetSMM || setAMD) legendSingle -> AddEntry(graph_r21_x,fParticleTitles[fiParticle],"p");
        else legendSingle -> AddEntry(graph_r21_x,fParticleTitles[fiParticle],"pl");
        draw(legendSingle,cvs_r21,iParticle*(1./fNumParticles),-1.,1./fNumParticles,spaceForLegend/fbn_r21.getFullWidth(),marginMarker);
        if (draw_yield) {
          draw(legendSingle,cvs2,iParticle*(1./fNumParticles),-1.,1./fNumParticles,spaceForLegend/fbn_r21.getFullWidth(),marginMarker);
          draw(legendSingle,cvs1,iParticle*(1./fNumParticles),-1.,1./fNumParticles,spaceForLegend/fbn_r21.getFullWidth(),marginMarker);
        }
      }

      if (fSetSMM) {
        auto paveText = newpt(("lines : %s",fSMMTitles[fIdxSMM]), cvs_r21,
            0, 1.-3.0*spaceForLegend/fbn_r21.getFullWidth(),
            1., spaceForLegend/fbn_r21.getFullWidth());
        cvs_r21 -> cd();
        paveText -> Draw();
      }

      if (setAMD && TString(bnvar->getExpression())=="y0")
      {
        for (auto iParticle : {0,1,2})
        {
          fiParticle = iParticle;
          auto graph_amd_r21 = new TGraphErrors();
          attGraph(graph_amd_r21, fiParticle);
          graph_amd_r21 -> SetLineStyle(2);

          for (auto ix=0; ix<20; ++ix) {
            double xValue = fAMDY0[ix];
            double yValue132 = fAMDValues[ix][0][iParticle];
            double yValue108 = fAMDValues[ix][1][iParticle];
            double r21Value = yValue132 / yValue108;
            if (bnvar->isInside(xValue)) graph_amd_r21 -> SetPoint(graph_amd_r21 -> GetN(), xValue, r21Value);
          }
          cvs_r21 -> cd();
          graph_amd_r21 -> Draw("samel");
        }

        auto paveText = newpt("dotted lines : AMD", cvs_r21,
            0, 1.-2*spaceForLegend/fbn_r21.getFullWidth(),
            1., spaceForLegend/fbn_r21.getFullWidth());
        cvs_r21 -> cd();
        paveText -> Draw();
      }

      if (0)
      {
        if (TString(bnx.getExpression())=="y0"&&draw_yield)
        {
          auto filer = new TFile(Form("canvas_r21_%s_%d_%d_%s.root",bnx.getExpression(),fSysBeams[fiSys2],fSysBeams[fiSys1],fTag.Data()),"recreate");
          cvs_r21 -> Write("c1");
          for (auto graph : graphsr)
            graph -> Write();

          auto file2 = new TFile(Form("canvas_dnd%s_%d_%s.root",bnx.getExpression(),fSysBeams[fiSys2],fTag.Data()),"recreate");
          cvs2 -> Write("c1");
          for (auto graph : graphs2)
            graph -> Write();

          auto file1 = new TFile(Form("canvas_dnd%s_%d_%s.root",bnx.getExpression(),fSysBeams[fiSys1],fTag.Data()),"recreate");
          cvs1 -> Write("c1");
          for (auto graph : graphs1)
            graph -> Write();
        }
      }

    }
  }

}

void draw_avb(binning bnx, binning bny, double r21Values[20][20][5][4], double abcValues[20][20][3][2], bool drawLegend, int is_x_smm_amd)
{
  auto var2 = (bnx*bny);

  int isN1XY = 0; // false:0, x:1, y:2
       if (bnx.getN()==1) isN1XY = 1;
  else if (bny.getN()==1) isN1XY = 2;

  auto name = make_name("avb",0,var2.namexynn());
  auto cvs = canvas(name,"single_sq");
  binning bn_alpha(100,-.1,.5,"#alpha");
  binning bn_beta(100,-.5,.1,"#beta");
  auto hist = (bn_alpha*bn_beta).newHist(name);
  TString addTitle;
  if (is_x_smm_amd==1) addTitle = "Calculation from S.R.Souza";
  else if (isN1XY==1) addTitle = Form("%s=%s~%s",bnx.getTitle(),toString(bnx.fMin).Data(),toString(bnx.fMax).Data());
  else if (isN1XY==2) addTitle = Form("%s=%s~%s",bny.getTitle(),toString(bny.fMin).Data(),toString(bny.fMax).Data());
  auto tt_sys = Form("#frac{%s}{%s}  %s",fSysTitles[fiSys2],fSysTitles[fiSys1],addTitle.Data());

  hist -> SetTitle(tt_sys);
  draw(hist,cvs);

  vector<TMarker *> markers;
  vector<TGraphErrors *> graph_errorls;

  bnx.reset();
  if (is_x_smm_amd==0)
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
      if (is_x_smm_amd==1) marker -> SetMarkerStyle(20);
      markers.push_back(marker);
    }

    graph_errorls.push_back(graph_errorl);
  }

  if (fSetSMM) {
    auto marker = new TMarker(fSMMValues[fIdxSMM][5],fSMMValues[fIdxSMM][6],20);
    marker -> SetMarkerSize(1.2);
    markers.push_back(marker);
  }

  if (drawLegend) {
    auto legend = new TLegend();
    bny.reset();
    if (isN1XY==0) {
      while (bny.next()) {
        auto marker = new TMarker(0,0,20);
        marker -> SetMarkerStyle(fDrawMStyle2[bny.ii()]);
        marker -> SetMarkerSize(fDrawMSize2[bny.ii()]);
        //legend -> AddEntry(marker,Form("%s = %.0f ~ %.0f",bny.getTitle(),bny.low(),bny.high()),"p");
        legend -> AddEntry(marker,Form("%s = %.1f ~ %.1f",bny.getTitle(),bny.low(),bny.high()),"p");
      }
      bnx.reset();
      while (bnx.next()) {
        auto graph = new TGraph();
        graph -> SetLineColor(fDrawColor[bnx.ii()]);
        legend -> AddEntry(graph,Form("%s = %.1f ~ %.1f",bnx.getTitle(),bnx.low(),bnx.high()),"l");
      }
    }
    if (isN1XY==1) {
      bny.reset();
      while (bny.next()) {
        auto marker = new TMarker(0,0,20);
        marker -> SetMarkerStyle(fDrawMStyle2[0]);
        marker -> SetMarkerSize(fDrawMSize2[0]);
        marker -> SetMarkerColor(fDrawColor[bny.ii()]);
        legend -> AddEntry(marker,Form("%s = %.1f ~ %.1f",bny.getTitle(),bny.low(),bny.high()),"p");
      }
    }
    if (isN1XY==2) {
      bnx.reset();
      while (bnx.next()) {
        auto marker = new TMarker(0,0,20);
        marker -> SetMarkerStyle(fDrawMStyle2[0]);
        marker -> SetMarkerSize(fDrawMSize2[0]);
        marker -> SetMarkerColor(fDrawColor[bnx.ii()]);
        legend -> AddEntry(marker,Form("%s = %.1f ~ %.1f",bnx.getTitle(),bnx.low(),bnx.high()),"p");
      }
    }
    make_legend(cvs,legend,"",0,0,.6,.4) -> Draw();
  }

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

void draw_apmb(binning bnx, binning bny, double r21Values[20][20][5][4], double abcValues[20][20][3][2], double tptValues[20][20][2][2], bool multt)
{
  auto legend1 = new TLegend();
  auto legend2 = new TLegend();

  auto bnab = binning(100,-.5,1);

  auto graph_alpha   = new TGraphErrors();
  graph_alpha -> SetMarkerStyle(24);
  graph_alpha -> SetMarkerSize(1.5);
  auto graph_beta    = new TGraphErrors();
  graph_beta -> SetMarkerStyle(25);
  graph_beta -> SetMarkerSize(1.5);
  auto graph_aplusb  = new TGraphErrors();
  graph_aplusb -> SetLineColor(fDrawColor[2]);
  graph_aplusb -> SetFillColor(fDrawColor3[2]);
  auto graph_aminusb = new TGraphErrors();
  graph_aminusb -> SetLineColor(fDrawColor[1]);
  graph_aminusb -> SetFillColor(fDrawColor3[1]);

  auto line0 = new TLine(bnx.fMin,0,bnx.fMax,0);
  line0 -> SetLineStyle(2);
  line0 -> SetLineColor(kGray+1);
  auto line5 = new TLine(bnx.fMin,.5,bnx.fMax,.5);
  line5 -> SetLineStyle(2);
  line5 -> SetLineColor(kGray+1);

  if (multt) {
    bnab.set(100,-5,10);
    line5 -> SetY1(5);
    line5 -> SetY2(5);

    legend1 -> AddEntry(graph_alpha  ,"#mu_{n}",         "p");
    legend1 -> AddEntry(graph_beta   ,"#mu_{p}",         "p");
    legend1 -> AddEntry(graph_aminusb,"#mu_{n}-#mu_{p}", "l");

    legend2 -> AddEntry(graph_alpha  ,"#mu_{n}",         "p");
    legend2 -> AddEntry(graph_beta   ,"#mu_{p}",         "p");
    legend2 -> AddEntry(graph_aplusb ,"#mu_{n}+#mu_{p}", "l");
  }
  else {
    legend1 -> AddEntry(graph_alpha  ,"#alpha",       "p");
    legend1 -> AddEntry(graph_beta   ,"#beta",        "p");
    legend1 -> AddEntry(graph_aminusb,"#alpha-#beta", "l");

    legend2 -> AddEntry(graph_alpha  ,"#alpha",       "p");
    legend2 -> AddEntry(graph_beta   ,"#beta",        "p");
    legend2 -> AddEntry(graph_aplusb ,"#alpha+#beta", "l");
  }

  bnx.reset();
  while (bnx.next())
  {
    auto alpha = abcValues[bnx.ii()][0][0][kVal];
    auto beta =  abcValues[bnx.ii()][0][1][kVal];
    auto alphaError = abcValues[bnx.ii()][0][0][kErr];
    auto betaError =  abcValues[bnx.ii()][0][1][kErr];
    auto abError = sqrt(alphaError*alphaError + betaError*betaError);

    if (multt) {
      auto temperature = (tptValues[bnx.ii()][0][0][kVal] + tptValues[bnx.ii()][0][1][kVal])/2.;
      alpha = alpha*temperature;
      beta = beta*temperature;
    }

    graph_alpha   -> SetPoint(bnx.ii(),bnx.val(),alpha);
    graph_beta    -> SetPoint(bnx.ii(),bnx.val(),beta);
    graph_alpha   -> SetPointError(bnx.ii(),0,alphaError);
    graph_beta    -> SetPointError(bnx.ii(),0,betaError);

    graph_aplusb  -> SetPoint(bnx.ii(),bnx.val(),alpha+beta);
    graph_aminusb -> SetPoint(bnx.ii(),bnx.val(),alpha-beta);
    graph_aplusb  -> SetPointError(bnx.ii(),0,abError);
    graph_aminusb -> SetPointError(bnx.ii(),0,abError);
  }

  auto name_apmb1 = make_name(Form("aminusb_%s",bnx.getExpression()));
  if (multt) name_apmb1 = make_name(Form("aminusbt_%s",bnx.getExpression()));
  auto hist1 = (bnx*bnab).newHist(name_apmb1);

  auto name_apmb2 = make_name(Form("aplusb_%s",bnx.getExpression()));
  if (multt) name_apmb2 = make_name(Form("aplusbt_%s",bnx.getExpression()));
  auto hist2 = (bnx*bnab).newHist(name_apmb2);

  auto iSys2 = fSysCombIdx[fiComb][0];
  auto iSys1 = fSysCombIdx[fiComb][1];
  if (TString(bny.getTitle()).Index("y_{0}")>=0) {
    hist1 -> SetTitle(Form("#frac{%s}{%s} : %s=%s~%s",fSysTitles[iSys2],fSysTitles[iSys1],bny.getTitle(),toString(bny.fMin).Data(),toString(bny.fMax).Data()));
    hist2 -> SetTitle(Form("#frac{%s}{%s} : %s=%s~%s",fSysTitles[iSys2],fSysTitles[iSys1],bny.getTitle(),toString(bny.fMin).Data(),toString(bny.fMax).Data()));
  }
  else {
    hist1 -> SetTitle(Form("#frac{%s}{%s} : %s=%s~%s",fSysTitles[iSys2],fSysTitles[iSys1],bny.getTitle(),toString(bny.fMin).Data(),toString(bny.fMax).Data()));
    hist2 -> SetTitle(Form("#frac{%s}{%s} : %s=%s~%s",fSysTitles[iSys2],fSysTitles[iSys1],bny.getTitle(),toString(bny.fMin).Data(),toString(bny.fMax).Data()));
  }
  hist1 -> SetYTitle("#alpha - #beta");
  hist2 -> SetYTitle("#alpha + #beta");

  if (multt) hist1 -> SetYTitle("#mu");
  if (multt) hist2 -> SetYTitle("#mu");

  auto cvs1 = canvas(name_apmb1,1,1,"single");
  draw(hist1,cvs1);
  line0 -> Draw("samel");
  line5 -> Draw("samel");
  graph_aminusb -> Draw("same3");
  graph_aminusb -> Draw("samelx");
  graph_alpha   -> Draw("samep");
  graph_beta    -> Draw("samep");
  make(legend1, cvs1, -1, -1, 0.25, 1./3, 0.3) -> Draw();

  auto cvs2 = canvas(name_apmb2,1,1,"single");
  draw(hist2,cvs2);
  line0 -> Draw("samel");
  line5 -> Draw("samel");
  graph_aplusb  -> Draw("same3");
  graph_aplusb  -> Draw("samelx");
  graph_alpha   -> Draw("samep");
  graph_beta    -> Draw("samep");
  make(legend2, cvs2, -1, -1, 0.25, 1./3, 0.3) -> Draw();
}

TLegend *make_legend(TVirtualPad *cvs, TLegend *legend, TString fLegendDrawStyle="rt", double x_offset, double y_offset, double width_fixed, double height_fixed)
{
  return make(legend,(TCanvas *)cvs,-1,-1,width_fixed,height_fixed);
}


const char *make_name(const char *name0, int idx, const char *tag) {
  const char *nameSC = ( fiComb>=0 ? Form("c%d_%s",fSysCombIdx2[fiComb],fSysCombNames[fiComb]) : Form("%d",fSysBeams[fiSys]) );
  const char *name = Form("%s_%s_%s_%s_%d",name0,tag,nameSC,fParticleNames[fiParticle],idx);
  return name;
}

const char *make_title(binning bn_range, int idx) { return make_title(&bn_range, idx); }

const char *make_title(binning *bn_range, int idx)
{
  const char *tt_range;
  if (!bn_range->isNull())
  {
    fiSys2 = fSysCombIdx[fiComb][0];
    fiSys1 = fSysCombIdx[fiComb][1];
    auto x1 = bn_range->fMin;
    auto x2 = bn_range->fMax;
    if (idx>=0) {
      x1 = bn_range->lowEdge(idx);
      x2 = bn_range->highEdge(idx);
    }
    const char *x1String = toString(x1);
    const char *x2String = toString(x2);
    tt_range = Form(" : %s=%s~%s",bn_range->getTitle(),x1String,x2String);
  }
  auto tt_sys = Form("#frac{%s}{%s}",fSysTitles[fiSys2],fSysTitles[fiSys1]);
  auto title = Form("%s%s",tt_sys,tt_range);
  return title;
}

const char *make_title1(binning *bn_range, int idx, int isys)
{
  const char *tt_range;
  if (!bn_range->isNull())
  {
    //isys = fSysCombIdx[fiComb][0];
    auto x1 = bn_range->fMin;
    auto x2 = bn_range->fMax;
    if (idx>=0) {
      x1 = bn_range->lowEdge(idx);
      x2 = bn_range->highEdge(idx);
    }
    const char *x1String = toString(x1);
    const char *x2String = toString(x2);
    tt_range = Form(" : %s=%s~%s",bn_range->getTitle(),x1String,x2String);
  }
  auto tt_sys = Form("(%s)",fSysTitles[isys]);
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

void save_cvs() {
  savePDF(fNameV.Data());
  savePNG(fNameV.Data());
}

TGraphErrors *attGraph(TGraphErrors *graph, int idx)
{
  att(graph,idx,"asame");
  //graph -> SetMarkerStyle(fDrawMStyle[idx]);
  //graph -> SetMarkerColor(fDrawColor[idx]);
  //graph -> SetMarkerSize(fDrawMSize[idx]);
  //graph -> SetLineColor(fDrawColor[idx]);
  //graph -> SetFillColor(fDrawColor3[idx]);
  return graph;
}
