void draw_ke_from_ana(
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

  TFile *file_hist = new TFile(Form("data_xml/summary_hist_keCM_%d_s%d.root",system,iSplit),"recreate");

  auto yBeamCM = 0.36;
  auto probCut = 0.2;
  auto effCut = 0.1;

  double particle_a[] = {1,2,3,3,4,};
  TString particle_name[] = {"p","d","t","he3","he4"};
  double particle_mass[] = {938.272 ,1871.06 ,2809.41 ,2809.41 ,3728.4};

  for (auto idx : {0,1,2,3,4})
  //for (auto idx : {0})
  //for (auto idx : {0,1,2})
  {
    auto mass = particle_mass[idx];
    TString ex = Form("sqrt(CMVector[%d].fElements.Mag2()+%f)-%f",idx,mass*mass, mass);

    TH1D *hist_keCM_cPE;
    {
      TString name = Form("keCM_cPE_%s",particle_name[idx].Data());
      TString title = particle_name[idx]+Form(" Sn%d (prob., eff. corrected);y_{CM}/y_{beam,CM};p_{T}/A (MeV/c)",system);
      hist_keCM_cPE = new TH1D(name,title,200,0,500);

      TString xy = ex + ">>" + hist_keCM_cPE -> GetName();
      TCut ww = Form("Prob[%d].fElements/Eff[%d].fElements",idx,idx);
      TCut cut = Form("(Prob[%d].fElements > %f && Eff[%d].fElements > %f)",idx,probCut,idx,effCut);
      TCut wcut = ww * cut;

      auto cvs = new TCanvas(name+"_canvas","",600,500);
      //cvs -> SetMargin(0.19,0.155,0.16,0.12);
      tree -> Draw(xy, wcut,"colz");

      if (write_to_file) {
        file_hist -> cd();
        hist_keCM_cPE -> Write();
      }
      cout << "end " << name << endl;
    }

    TH1D *hist_keCM_cP;
    {
      TString name = Form("keCM_cP_%s",particle_name[idx].Data());
      TString title = particle_name[idx]+Form(" Sn%d (prob. corrected);y_{CM}/y_{beam,CM};p_{T}/A (MeV/c)",system);
      hist_keCM_cP = new TH1D(name,title,200,0,500);

      TString xy = ex + ">>" + hist_keCM_cP -> GetName();
      TCut ww = Form("Eff[%d].fElements",idx);
      TCut cut = Form("(Prob[%d].fElements > %f && Eff[%d].fElements > %f)",idx,probCut,idx,effCut);
      TCut wcut = ww * cut;

      auto cvs = new TCanvas(name+"_canvas","",600,500);
      //cvs -> SetMargin(0.19,0.155,0.16,0.12);
      tree -> Draw(xy, wcut,"colz");

      if (write_to_file) {
        file_hist -> cd();
        hist_keCM_cP -> Write();
      }
      cout << "end " << name << endl;
    }

    TH1D *hist_keCM;
    {
      TString name = Form("keCM_%s",particle_name[idx].Data());
      TString title = particle_name[idx]+Form(" Sn%d;y_{CM}/y_{beam,CM};p_{T}/A (MeV/c)",system);
      hist_keCM = new TH1D(name,title,200,0,500);

      TString xy = ex + ">>" + hist_keCM -> GetName();
      TCut ww = "1";
      TCut cut = Form("(Prob[%d].fElements > %f && Eff[%d].fElements > %f)",idx,probCut,idx,effCut);
      TCut wcut = ww * cut;

      auto cvs = new TCanvas(name+"_canvas","",600,500);
      //cvs -> SetMargin(0.19,0.155,0.16,0.12);
      tree -> Draw(xy, wcut,"colz");

      if (write_to_file) {
        file_hist -> cd();
        hist_keCM -> Write();
      }
      cout << "end " << name << endl;
    }

    TH1D *hist_eff;
    {
      TString name = Form("eff_%s",particle_name[idx].Data());
      TString title = particle_name[idx]+Form(" Sn%d efficiency;y_{CM}/y_{beam,CM};p_{T}/A (MeV/c)",system);
      hist_eff = (TH1D *) hist_keCM_cP -> Clone(name);
      hist_eff -> Divide(hist_keCM_cPE);
      hist_eff -> SetMinimum(0);
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

    TH1D *hist_prob;
    {
      TString name = Form("prob_%s",particle_name[idx].Data());
      TString title = particle_name[idx]+Form(" Sn%d probability;y_{CM}/y_{beam,CM};p_{T}/A (MeV/c)",system);
      hist_prob = (TH1D *) hist_keCM_cP -> Clone(name);
      hist_prob -> Divide(hist_keCM);
      hist_prob -> SetTitle(title);
      hist_prob -> SetMinimum(0);

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
