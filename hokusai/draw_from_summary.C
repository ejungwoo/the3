#include "/home/ejungwoo/config/ejungwoo.h"

void draw_from_summary(
  int system = 108,
  int version = 1
  )
{
  ejungwoo::gsave(0);
  ejungwoo::gversion("sum");
  
  const char *particle_names[] = {"p","d","t","he3","he4"};
  const int numSplits = 20;

  TFile *files[numSplits];
  for (auto iSplit=0; iSplit<numSplits; ++iSplit)
    files[iSplit] = new TFile(Form("data_with_trim/summary_hist_%d_s%d_v%d.root",system,iSplit,version));

  TString file_name_out = Form("data_with_trim/summary_hist_%d_all_v%d.root",system,version);
  auto file_out = new TFile(file_name_out,"recreate");

  TString values[] = {
    "yypt",
    "yypt_cP",
    "yypt_cPE",
    "yypt_effi",
    "yypt_prob",
    "keCM",
    "keCM_cP",
    "keCM_cPE",
    "keCM_effi",
    "keCM_prob",
    "phikeCM",
    "phikeCM_cPE",
    "ttakeCM",
    "ttakeCM_cPE",
    "phittaCM",
    "phittaCM_cP",
    "phittaCM_cPE",
    "phittaLab",
    "phittaLab_cP",
    "phittaLab_cPE",
    "pdedx_raw",
    "pdedx_cut"
  };

  vector<TH1 *> hsums;

  for (auto pidx : {0,1,2,3,4})
  {
    TString pname = particle_names[pidx];

    hsums.clear();
    for (auto value : values)
      hsums.push_back((TH1 *) nullptr);

    auto ivalue = 0;
    for (auto value : values) {
      auto name = pname + "_" + value;
      auto hsum = hsums.at(ivalue++);

      for (auto iSplit=0; iSplit<numSplits; ++iSplit) {
        auto file = files[iSplit];
        auto hist = (TH1 *) file -> Get(name);
        if (hsum == nullptr) hsum = hist;
        else                 hsum -> Add(hist);
      }

      if (value.Index("effi")>=0 || value.Index("prob")>=0) {
        if (hsum -> InheritsFrom("TH2")) {
          auto hdiv = (TH2 *) hsum -> Clone(Form("div_%s",hsum->GetName()));
          auto nx = hdiv -> GetNbinsX();
          auto ny = hdiv -> GetNbinsY();
          for (auto ix=1; ix<=nx; ++ix)
            for (auto iy=1; iy<=ny; ++iy)
              hdiv -> SetBinContent(ix,iy,numSplits);
          hsum -> Divide(hdiv);
        }
        else {
          auto hdiv = (TH1 *) hsum -> Clone(Form("div_%s",hsum->GetName()));
          auto nx = hdiv -> GetNbinsX();
          for (auto ix=1; ix<=nx; ++ix)
            hdiv -> SetBinContent(ix,numSplits);
          hsum -> Divide(hdiv);
        }
      }

      file_out -> cd();
      hsum -> Write();
      cout << name << endl;
    }
  }

  cout << file_name_out << endl;
}
