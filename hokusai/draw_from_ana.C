void draw_from_ana(
  int system=132,
  int iSplit=1
  //int numSplits=1
  )
{
  gStyle -> SetOptStat(0);

  bool write_to_file = true;

  auto tree = new TChain("cbmsim");
  TString path1 = "/home/ejungwoo/data/pid4/Sn";
  //for (auto iSplit=0; iSplit<numSplits; ++iSplit) {
    TString fileName = path1+system+"_"+iSplit+"_ana.root";
    tree -> Add(fileName);
  //}

  TFile *file_hist = new TFile(Form("data_xml/summary_hist_yypt_%d_s%d.root",system,iSplit),"recreate");

  auto yBeamCM = 0.36;
  auto probCut = 0.2;
  auto effCut = 0.1;

  double particle_a[] = {1,2,3,3,4,};
  TString particle_name[] = {"p","d","t","he3","he4"};

  for (auto idx : {0,1,2,3,4})
  //for (auto idx : {0})
  //for (auto idx : {0,1,2})
  {
    TH2D *hist_yypt_cPE;
    {
      TString name = Form("yypt_cPE_%s",particle_name[idx].Data());
      TString title = particle_name[idx]+Form(" Sn%d (prob., eff. corrected);y_{CM}/y_{beam,CM};p_{T}/A (MeV/c)",system);
      hist_yypt_cPE = new TH2D(name,title,200,-1.5,3.0,200,0,1500);

      TString ey = Form("%f*CMVector[%d].fElements.Perp()",1./particle_a[idx],idx);
      TString ex = Form("FragRapidity[%d].fElements/%f",idx,yBeamCM);
      TString xy = ey + ":" + ex + ">>" + hist_yypt_cPE -> GetName();
      TCut ww = Form("Prob[%d].fElements/Eff[%d].fElements",idx,idx);
      TCut cut = Form("(Prob[%d].fElements > %f && Eff[%d].fElements > %f)",idx,probCut,idx,effCut);
      TCut wcut = ww * cut;

      auto cvs = new TCanvas(name+"_canvas","",600,500);
      //cvs -> SetMargin(0.19,0.155,0.16,0.12);
      tree -> Draw(xy, wcut,"colz");

      if (write_to_file) {
        file_hist -> cd();
        hist_yypt_cPE -> Write();
      }
      cout << "end " << name << endl;
    }

    TH2D *hist_yypt_cP;
    {
      TString name = Form("yypt_cP_%s",particle_name[idx].Data());
      TString title = particle_name[idx]+Form(" Sn%d (prob. corrected);y_{CM}/y_{beam,CM};p_{T}/A (MeV/c)",system);
      hist_yypt_cP = new TH2D(name,title,200,-1.5,3.0,200,0,1500);

      TString ey = Form("%f*CMVector[%d].fElements.Perp()",1./particle_a[idx],idx);
      TString ex = Form("FragRapidity[%d].fElements/%f",idx,yBeamCM);
      TString xy = ey + ":" + ex + ">>" + hist_yypt_cP -> GetName();
      TCut ww = Form("Eff[%d].fElements",idx);
      TCut cut = Form("(Prob[%d].fElements > %f && Eff[%d].fElements > %f)",idx,probCut,idx,effCut);
      TCut wcut = ww * cut;

      auto cvs = new TCanvas(name+"_canvas","",600,500);
      //cvs -> SetMargin(0.19,0.155,0.16,0.12);
      tree -> Draw(xy, wcut,"colz");

      if (write_to_file) {
        file_hist -> cd();
        hist_yypt_cP -> Write();
      }
      cout << "end " << name << endl;
    }

    TH2D *hist_yypt;
    {
      TString name = Form("yypt_%s",particle_name[idx].Data());
      TString title = particle_name[idx]+Form(" Sn%d;y_{CM}/y_{beam,CM};p_{T}/A (MeV/c)",system);
      hist_yypt = new TH2D(name,title,200,-1.5,3.0,200,0,1500);

      TString ey = Form("%f*CMVector[%d].fElements.Perp()",1./particle_a[idx],idx);
      TString ex = Form("FragRapidity[%d].fElements/%f",idx,yBeamCM);
      TString xy = ey + ":" + ex + ">>" + hist_yypt -> GetName();
      TCut ww = "1";
      TCut cut = Form("(Prob[%d].fElements > %f && Eff[%d].fElements > %f)",idx,probCut,idx,effCut);
      TCut wcut = ww * cut;

      auto cvs = new TCanvas(name+"_canvas","",600,500);
      //cvs -> SetMargin(0.19,0.155,0.16,0.12);
      tree -> Draw(xy, wcut,"colz");

      if (write_to_file) {
        file_hist -> cd();
        hist_yypt -> Write();
      }
      cout << "end " << name << endl;
    }

    TH2D *hist_eff;
    {
      TString name = Form("eff_%s",particle_name[idx].Data());
      TString title = particle_name[idx]+Form(" Sn%d efficiency;y_{CM}/y_{beam,CM};p_{T}/A (MeV/c)",system);
      hist_eff = (TH2D *) hist_yypt_cP -> Clone(name);
      hist_eff -> Divide(hist_yypt_cPE);
      hist_eff -> SetTitle(title);

      auto cvs = new TCanvas(name+"_canvas","",600,500);
      //cvs -> SetMargin(0.19,0.155,0.16,0.12);
      hist_eff -> Draw("colz");

      if (write_to_file) {
        file_hist -> cd();
        hist_eff -> Write();
      }
      cout << "end " << name << endl;
    }

    TH2D *hist_prob;
    {
      TString name = Form("prob_%s",particle_name[idx].Data());
      TString title = particle_name[idx]+Form(" Sn%d probability;y_{CM}/y_{beam,CM};p_{T}/A (MeV/c)",system);
      hist_prob = (TH2D *) hist_yypt_cP -> Clone(name);
      hist_prob -> Divide(hist_yypt);
      hist_prob -> SetTitle(title);

      auto cvs = new TCanvas(name+"_canvas","",600,500);
      //cvs -> SetMargin(0.19,0.155,0.16,0.12);
      hist_prob -> Draw("colz");

      if (write_to_file) {
        file_hist -> cd();
        hist_prob -> Write();
      }
      cout << "end " << name << endl;
    }
  }
}
