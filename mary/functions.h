#ifndef FUNCTION_H
#define FUNCTION_H

#include "init_variables.h"

TObjArray draw_is(int idx_particle, TH1 *hist132, TH1 *hist108, double num_tracks_cut)
{
  auto graph_is = new TGraphErrors();
  auto graph_v132 = new TGraphErrors();
  auto graph_v108 = new TGraphErrors();
  binning binn(hist132);

  for (auto bin=1; bin<=binn.n; ++bin)
  {
    auto binc = binn.getc(bin);
    auto v132 = hist132 -> GetBinContent(bin);
    auto v108 = hist108 -> GetBinContent(bin);
    auto idxvv = graph_v132 -> GetN();
    graph_v132 -> SetPoint(idxvv, binc, v132);
    graph_v108 -> SetPoint(idxvv, binc, v108);
    if (v132 < num_tracks_cut || v108 < num_tracks_cut)
      continue;
    auto idxis = graph_is -> GetN();
    graph_is -> SetPoint(idxis, binc, v132/v108);
  }

  TObjArray array;
  for (auto graph : {graph_is, graph_v132, graph_v108}) {
    graph -> SetMarkerStyle(ejungwoo::markeri(idx_particle));
    graph -> SetLineColor(ejungwoo::colori(idx_particle));
    graph -> SetLineWidth(2);
    array.Add(graph);
  }

  return array;
}

TObjArray draw_ci(vector<TH1 *>hists, double ds, int isys)
{
  auto graph_np = new TGraphErrors();
  auto graph_nn = new TGraphErrors();
  auto graph_rr = new TGraphErrors();

  binning binn(hists[0]);
  binn.resetb();
  while (binn.nextb()) {
    binn.idx;
    double vpsum = 0;
    double vnsum = 0;
    for (auto ipid : fIPIDAll) {
      double ntracks = hists[ipid] -> GetBinContent(binn.idx);
      auto np = fNumProtons[ipid];
      auto nn = fNumNeutrons[ipid];
      vpsum += fNumProtons[ipid] * ntracks;
      vnsum += fNumNeutrons[ipid] * ntracks;
    }
    graph_np -> SetPoint(graph_np->GetN(), binn.value, vpsum/binn.w/ds);
    graph_nn -> SetPoint(graph_nn->GetN(), binn.value, vnsum/binn.w/ds);
    graph_rr -> SetPoint(graph_rr->GetN(), binn.value, vnsum/vpsum);
  }

  TObjArray array;
  int idx = 0;
  for (auto graph : {graph_nn, graph_np, graph_rr}) {
    graph -> SetMarkerStyle(ejungwoo::markeri(idx));
    graph -> SetMarkerColor(ejungwoo::colori(idx));
    graph -> SetLineColor(ejungwoo::colori(isys));
    graph -> SetLineWidth(2);
    array.Add(graph);
    idx++;
  }

  return array;
}

TGraphErrors *draw_ratio(TGraphErrors *graph1, TGraphErrors *graph2, double max)
{
  auto graph_dr = new TGraphErrors();

  double x1, y1, x2, y2;
  auto nn = graph1 -> GetN();
  for (auto i=0; i<nn; ++i) {
    graph1 -> GetPoint(i,x1,y1);
    graph2 -> GetPoint(i,x2,y2);
    if (x1 > max)
      break;
    if (y2>0) {
      graph_dr -> SetPoint(graph_dr->GetN(),x1,y1/y2);
    }
  }

  return graph_dr;
}

#endif
