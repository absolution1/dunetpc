//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Fri May  2 15:13:16 2014 by ROOT version 5.34/03
// from TTree RawData/Raw Data Display
// found on file: duneEVDraw35t_tree_Cosmic_1000evt.root
//////////////////////////////////////////////////////////

#ifndef MakeHistograms_h
#define MakeHistograms_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TH2.h>

// Header file for the classes stored in the TTree if any.
#include <vector>

// Fixed size dimensions of array or collections stored in the TTree if any.

class MakeHistograms {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain
   UInt_t          fEvent;
	
   // Declaration of leaf types
   UInt_t          Run;
   UInt_t          SubRun;
   UInt_t          Event;
   vector<vector<int> > *ADC;
   vector<vector<int> > *TDC;
   UInt_t          Nticks;
   UInt_t          NofUChan;
   UInt_t          NofVChan;
   UInt_t          NofZ0Chan;
   UInt_t          NofZ1Chan;
   vector<unsigned int> *Chan;
   vector<unsigned int> *APA;
   vector<unsigned int> *Plane;

   // List of branches
   TBranch        *b_Run;   //!
   TBranch        *b_SubRun;   //!
   TBranch        *b_Event;   //!
   TBranch        *b_ADC;   //!
   TBranch        *b_TDC;   //!
   TBranch        *b_Nticks;   //!
   TBranch        *b_NofUChan;   //!
   TBranch        *b_NofVChan;   //!
   TBranch        *b_NofZ0Chan;   //!
   TBranch        *b_NofZ1Chan;   //!
   TBranch        *b_Chan;   //!
   TBranch        *b_APA;   //!
   TBranch        *b_Plane;   //!

//   MakeHistograms(TTree *tree=0);
   MakeHistograms(std::string dfname="duneEVDraw35t_tree_Cosmic_1000evt.root");
   virtual ~MakeHistograms();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     PrintPlots(UInt_t eventnb=1);
   virtual void     HistDef();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
	virtual void 	  Pal1();
	virtual void	  Pal2();
	virtual void	  PrintHistos();
	virtual void	  PrintGraphs(std::vector<TH2I*> h);
	virtual void 	  DrawUplane();
	virtual void 	  DrawVplane();
	virtual void 	  DrawZplane();

private :
   std::vector<TH2I*> fTimeChanU;
   std::vector<TH2I*> fTimeChanV;
   std::vector<TH2I*> fTimeChanZ0;
   std::vector<TH2I*> fTimeChanZ1;
   TH2I* fChargeSumU;
   TH2I* fChargeSumV;
   TH2I* fChargeSumZ;
	TCanvas* fCanvas;
};

#endif

#ifdef MakeHistograms_cxx
//MakeHistograms::MakeHistograms(TTree *tree) : fChain(0) 
MakeHistograms::MakeHistograms(std::string dfname) : fChain(0)
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
	TTree *tree=0; 
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject(dfname.c_str());
      if (!f || !f->IsOpen()) {
         f = new TFile(dfname.c_str());
      }
		dfname=dfname+":/rawdraw";
      TDirectory * dir = (TDirectory*)f->Get(dfname.c_str());
      dir->GetObject("RawData",tree);

   }
//	fCanvas = new TCanvas("Ecanvas","DUNE 35t detector event display",800,600);
   Init(tree);
	HistDef();
}

MakeHistograms::~MakeHistograms()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t MakeHistograms::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t MakeHistograms::LoadTree(Long64_t entry)
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

void MakeHistograms::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   ADC = 0;
   TDC = 0;
   Chan = 0;
   APA = 0;
   Plane = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("Run", &Run, &b_Run);
   fChain->SetBranchAddress("SubRun", &SubRun, &b_SubRun);
   fChain->SetBranchAddress("Event", &Event, &b_Event);
   fChain->SetBranchAddress("ADC", &ADC, &b_ADC);
   fChain->SetBranchAddress("TDC", &TDC, &b_TDC);
   fChain->SetBranchAddress("Nticks", &Nticks, &b_Nticks);
   fChain->SetBranchAddress("NofUChan", &NofUChan, &b_NofUChan);
   fChain->SetBranchAddress("NofVChan", &NofVChan, &b_NofVChan);
   fChain->SetBranchAddress("NofZ0Chan", &NofZ0Chan, &b_NofZ0Chan);
   fChain->SetBranchAddress("NofZ1Chan", &NofZ1Chan, &b_NofZ1Chan);
   fChain->SetBranchAddress("Chan", &Chan, &b_Chan);
   fChain->SetBranchAddress("APA", &APA, &b_APA);
   fChain->SetBranchAddress("Plane", &Plane, &b_Plane);
   Notify();
}

Bool_t MakeHistograms::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void MakeHistograms::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t MakeHistograms::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef MakeHistograms_cxx
