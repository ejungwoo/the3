//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Mon Aug 17 20:18:05 2020 by ROOT version 6.13/02
// from TTree spirit/
// found on file: /data/Q20393/production/20200529/SpiRITROOT/macros/data/Sn132/run2841_s0.reco.develop.1988.bf2b00e.conc.root
//////////////////////////////////////////////////////////

#ifndef STConcReader_h
#define STConcReader_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include "TObject.h"

class STConcReader {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   Int_t           runID;
 //STData          *EvtData;
   UInt_t          fUniqueID;
   UInt_t          fBits;
   Double_t        aoq;
   Double_t        z;
   Double_t        a;
   Double_t        b;
   Double_t        proja;
   Double_t        projb;
   Double_t        projx;
   Double_t        projy;
   Double_t        beamEnergy;
   Double_t        beta;
   TVector3        tpcVertex;
   TVector3        bdcVertex;
   Int_t           multiplicity;
   vector<TVector3> recoMom;
   vector<TVector3> recoPosPOCA;
   vector<TVector3> recoPosTargetPlane;
   vector<TVector3> recodpoca;
   vector<int>     recoNRowClusters;
   vector<int>     recoNLayerClusters;
   vector<int>     recoCharge;
   vector<bool>    recoEmbedTag;
   vector<double>  recodedx;
   Int_t           vaMultiplicity;
   vector<TVector3> vaMom;
   vector<TVector3> vaPosPOCA;
   vector<TVector3> vaPosTargetPlane;
   vector<TVector3> vadpoca;
   vector<int>     vaNRowClusters;
   vector<int>     vaNLayerClusters;
   vector<int>     vaCharge;
   vector<bool>    vaEmbedTag;
   vector<double>  vadedx;
   TVector3        embedMom;
   Double_t        beamEnergyTargetPlane;
   Double_t        betaTargetPlane;
   Int_t           eventID;
   Int_t           eventType;

   // List of branches
   TBranch        *b_runID;   //!
   TBranch        *b_EvtData_fUniqueID;   //!
   TBranch        *b_EvtData_fBits;   //!
   TBranch        *b_EvtData_aoq;   //!
   TBranch        *b_EvtData_z;   //!
   TBranch        *b_EvtData_a;   //!
   TBranch        *b_EvtData_b;   //!
   TBranch        *b_EvtData_proja;   //!
   TBranch        *b_EvtData_projb;   //!
   TBranch        *b_EvtData_projx;   //!
   TBranch        *b_EvtData_projy;   //!
   TBranch        *b_EvtData_beamEnergy;   //!
   TBranch        *b_EvtData_beta;   //!
   TBranch        *b_EvtData_tpcVertex;   //!
   TBranch        *b_EvtData_bdcVertex;   //!
   TBranch        *b_EvtData_multiplicity;   //!
   TBranch        *b_EvtData_recoMom;   //!
   TBranch        *b_EvtData_recoPosPOCA;   //!
   TBranch        *b_EvtData_recoPosTargetPlane;   //!
   TBranch        *b_EvtData_recodpoca;   //!
   TBranch        *b_EvtData_recoNRowClusters;   //!
   TBranch        *b_EvtData_recoNLayerClusters;   //!
   TBranch        *b_EvtData_recoCharge;   //!
   TBranch        *b_EvtData_recoEmbedTag;   //!
   TBranch        *b_EvtData_recodedx;   //!
   TBranch        *b_EvtData_vaMultiplicity;   //!
   TBranch        *b_EvtData_vaMom;   //!
   TBranch        *b_EvtData_vaPosPOCA;   //!
   TBranch        *b_EvtData_vaPosTargetPlane;   //!
   TBranch        *b_EvtData_vadpoca;   //!
   TBranch        *b_EvtData_vaNRowClusters;   //!
   TBranch        *b_EvtData_vaNLayerClusters;   //!
   TBranch        *b_EvtData_vaCharge;   //!
   TBranch        *b_EvtData_vaEmbedTag;   //!
   TBranch        *b_EvtData_vadedx;   //!
   TBranch        *b_EvtData_embedMom;   //!
   TBranch        *b_EvtData_beamEnergyTargetPlane;   //!
   TBranch        *b_EvtData_betaTargetPlane;   //!
   TBranch        *b_eventID;   //!
   TBranch        *b_eventType;   //!

   STConcReader(TTree *tree=0);
   virtual ~STConcReader();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef STConcReader_cxx
STConcReader::STConcReader(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("/data/Q20393/production/20200529/SpiRITROOT/macros/data/Sn132/run2841_s0.reco.develop.1988.bf2b00e.conc.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("/data/Q20393/production/20200529/SpiRITROOT/macros/data/Sn132/run2841_s0.reco.develop.1988.bf2b00e.conc.root");
      }
      f->GetObject("spirit",tree);

   }
   Init(tree);
}

STConcReader::~STConcReader()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t STConcReader::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t STConcReader::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void STConcReader::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("runID", &runID, &b_runID);
   fChain->SetBranchAddress("fUniqueID", &fUniqueID, &b_EvtData_fUniqueID);
   fChain->SetBranchAddress("fBits", &fBits, &b_EvtData_fBits);
   fChain->SetBranchAddress("aoq", &aoq, &b_EvtData_aoq);
   fChain->SetBranchAddress("z", &z, &b_EvtData_z);
   fChain->SetBranchAddress("a", &a, &b_EvtData_a);
   fChain->SetBranchAddress("b", &b, &b_EvtData_b);
   fChain->SetBranchAddress("proja", &proja, &b_EvtData_proja);
   fChain->SetBranchAddress("projb", &projb, &b_EvtData_projb);
   fChain->SetBranchAddress("projx", &projx, &b_EvtData_projx);
   fChain->SetBranchAddress("projy", &projy, &b_EvtData_projy);
   fChain->SetBranchAddress("beamEnergy", &beamEnergy, &b_EvtData_beamEnergy);
   fChain->SetBranchAddress("beta", &beta, &b_EvtData_beta);
   fChain->SetBranchAddress("tpcVertex", &tpcVertex, &b_EvtData_tpcVertex);
   fChain->SetBranchAddress("bdcVertex", &bdcVertex, &b_EvtData_bdcVertex);
   fChain->SetBranchAddress("multiplicity", &multiplicity, &b_EvtData_multiplicity);
   fChain->SetBranchAddress("recoMom", &recoMom, &b_EvtData_recoMom);
   fChain->SetBranchAddress("recoPosPOCA", &recoPosPOCA, &b_EvtData_recoPosPOCA);
   fChain->SetBranchAddress("recoPosTargetPlane", &recoPosTargetPlane, &b_EvtData_recoPosTargetPlane);
   fChain->SetBranchAddress("recodpoca", &recodpoca, &b_EvtData_recodpoca);
   fChain->SetBranchAddress("recoNRowClusters", &recoNRowClusters, &b_EvtData_recoNRowClusters);
   fChain->SetBranchAddress("recoNLayerClusters", &recoNLayerClusters, &b_EvtData_recoNLayerClusters);
   fChain->SetBranchAddress("recoCharge", &recoCharge, &b_EvtData_recoCharge);
   fChain->SetBranchAddress("recoEmbedTag", &recoEmbedTag, &b_EvtData_recoEmbedTag);
   fChain->SetBranchAddress("recodedx", &recodedx, &b_EvtData_recodedx);
   fChain->SetBranchAddress("vaMultiplicity", &vaMultiplicity, &b_EvtData_vaMultiplicity);
   fChain->SetBranchAddress("vaMom", &vaMom, &b_EvtData_vaMom);
   fChain->SetBranchAddress("vaPosPOCA", &vaPosPOCA, &b_EvtData_vaPosPOCA);
   fChain->SetBranchAddress("vaPosTargetPlane", &vaPosTargetPlane, &b_EvtData_vaPosTargetPlane);
   fChain->SetBranchAddress("vadpoca", &vadpoca, &b_EvtData_vadpoca);
   fChain->SetBranchAddress("vaNRowClusters", &vaNRowClusters, &b_EvtData_vaNRowClusters);
   fChain->SetBranchAddress("vaNLayerClusters", &vaNLayerClusters, &b_EvtData_vaNLayerClusters);
   fChain->SetBranchAddress("vaCharge", &vaCharge, &b_EvtData_vaCharge);
   fChain->SetBranchAddress("vaEmbedTag", &vaEmbedTag, &b_EvtData_vaEmbedTag);
   fChain->SetBranchAddress("vadedx", &vadedx, &b_EvtData_vadedx);
   fChain->SetBranchAddress("embedMom", &embedMom, &b_EvtData_embedMom);
   fChain->SetBranchAddress("beamEnergyTargetPlane", &beamEnergyTargetPlane, &b_EvtData_beamEnergyTargetPlane);
   fChain->SetBranchAddress("betaTargetPlane", &betaTargetPlane, &b_EvtData_betaTargetPlane);
   fChain->SetBranchAddress("eventID", &eventID, &b_eventID);
   fChain->SetBranchAddress("eventType", &eventType, &b_eventType);
   Notify();
}

Bool_t STConcReader::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void STConcReader::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t STConcReader::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef STConcReader_cxx
