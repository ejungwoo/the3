#include "/home/ejungwoo/config/ejungwoo.h"

void draw_from_summary(
  int system = 132,
  int version = 1,
  int probv = 8
  )
{
  ejungwoo::gsave(0);
  
  const char *particle_names[] = {"p","d","t","he3","he4"};
  const int numSplits = 20;

  ///home/ejungwoo/the3/hokusai/submit_ana/data_with_trim/summary_hist_108_s0_v0_e5_p8.root

  TFile *files[numSplits];
  int num_events[numSplits];
  int num_total_events = 0;
  for (auto iSplit=0; iSplit<numSplits; ++iSplit) {
    TString file_name = Form("data_with_trim/summary_hist_%d_s%d_v%d_e1_p%d.root",system,iSplit,version,probv);
    files[iSplit] = new TFile(file_name);
    num_events[iSplit] = ((TParameter<int>*)files[iSplit]->Get("num_events"))->GetVal();
    cout << file_name << " " << iSplit << " " << num_events[iSplit] << endl;
    num_total_events += num_events[iSplit];
  }

  TString file_name_out = Form("data_with_trim/summary_hist_%d_all_v%d_e1_p%d.root",system,version,probv);
  auto file_out = new TFile(file_name_out,"recreate");
  (new TParameter<int>("num_events",num_total_events)) -> Write();

  vector<TString> values;
  auto keys = files[0] -> GetListOfKeys();
  auto next = TIter(keys);
  while (TObject *key = next()) {
    TString name = key -> GetName();
    if (name.Index("p_")==0) {
      name.ReplaceAll("p_","");
      values.push_back(name);
    }
  }

  vector<TH1 *> hsums;

  for (auto pidx : {0,1,2,3,4})
  {
    TString pname = particle_names[pidx];
    cout << endl;
    cout << pname << endl;

    hsums.clear();
    for (auto value : values)
      hsums.push_back((TH1 *) nullptr);

    auto ivalue = 0;
    for (auto value : values)
    {
      auto name = pname + "_" + value;
      auto hsum = hsums.at(ivalue++);

      for (auto iSplit=0; iSplit<numSplits; ++iSplit) {
        auto file = files[iSplit];
        auto hist = (TH1 *) file -> Get(name);
        if (hsum == nullptr) {
          hsum = ejungwoo::new_h(name+"_total",hist);
          hsum -> SetTitle(hist->GetTitle());
          hsum -> GetXaxis() -> SetTitle(hist -> GetXaxis() -> GetTitle());
          hsum -> GetYaxis() -> SetTitle(hist -> GetYaxis() -> GetTitle());
          hsum -> GetZaxis() -> SetTitle(hist -> GetZaxis() -> GetTitle());
          ejungwoo::make_h(hsum);
        }

        if (value.Index("eff")>=0 || value.Index("prob")>=0) {
          hsum -> Add(hist,1./numSplits);
        }
        else {
          hsum -> Add(hist,Double_t(num_events[iSplit])/num_total_events);
        }
      }

      cout << setw(30) << hsum -> GetName() 
        << setw(15) << hsum -> GetEntries() 
        << setw(15) << hsum -> Integral()
        << endl;

      file_out -> cd();
      hsum -> Write();
    }
  }

  cout << file_name_out << endl;
}

  /*
  TString values[] = {
    "eff_x1",
    "prob_x1",
    "pLab_dedx_x1",
    "pozLab_dedx_x1",
    "keCM_dedx_x1",
    "keoaCM_dedx_x1",
    "keozCM_dedx_x1",
    "keCM_x1",
    "keoaCM_x1",
    "keozCM_x1",
    "ptoaCM_x1",
    "yyCM_x1",
    "yyCM_ptoaCM_x1",
    "phiCM_keoaCM_x1",
    "ttaCM_keoaCM_x1",
    "phiCM_ttaCM_x1",
    "phiLab_ttaLab_x1",
    "keCM_pozLab_x1",
    "keoaCM_pozLab_x1",
    "keCM_xP",
    "keoaCM_xP",
    "keozCM_xP",
    "ptoaCM_xP",
    "yyCM_xP",
    "yyCM_ptoaCM_xP",
    "phiCM_keoaCM_xP",
    "ttaCM_keoaCM_xP",
    "phiCM_ttaCM_xP",
    "phiLab_ttaLab_xP",
    "keCM_pozLab_xP",
    "keoaCM_pozLab_xP",
    "keCM_xPoE",
    "keoaCM_xPoE",
    "keozCM_xPoE",
    "ptoaCM_xPoE",
    "yyCM_xPoE",
    "yyCM_ptoaCM_xPoE",
    "phiCM_keoaCM_xPoE",
    "ttaCM_keoaCM_xPoE",
    "phiCM_ttaCM_xPoE",
    "phiLab_ttaLab_xPoE",
    "keCM_pozLab_xPoE",
    "keoaCM_pozLab_xPoE",
  };
  */
