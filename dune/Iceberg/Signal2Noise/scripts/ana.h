//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Tue Sep  1 23:54:29 2020 by ROOT version 6.18/04
// from TTree Event/Event
// found on file: signal2noise.root
//////////////////////////////////////////////////////////

#ifndef ana_h
#define ana_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.

class ana {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   Int_t           event;
   Int_t           run;
   Int_t           subrun;
   Double_t        evttime;
   Int_t           year_month_date;
   Int_t           hour_min_sec;
   Int_t           ntrks;
   Int_t           trkid[336];   //[ntrks]
   Float_t         trkstart[336][3];   //[ntrks]
   Float_t         trkend[336][3];   //[ntrks]
   Float_t         trklen[336];   //[ntrks]
   Float_t         trkthetaxz[336];   //[ntrks]
   Float_t         trkthetayz[336];   //[ntrks]
   Float_t         trkstartcosxyz[336][3];   //[ntrks]
   Float_t         trkendcosxyz[336][3];   //[ntrks]
   Int_t           ntrkhits[336][3];   //[ntrks]
   Float_t         trkdqdx[336][3][1000];   //[ntrks]
   Float_t         trkx[336][3][1000];   //[ntrks]
   Float_t         trkt[336][3][1000];   //[ntrks]
   Double_t        trkhitx[336][3][1000];   //[ntrks]
   Double_t        trkhity[336][3][1000];   //[ntrks]
   Double_t        trkhitz[336][3][1000];   //[ntrks]
   Int_t           wireid[336][1000];   //[ntrks]
   Int_t           chid[336][1000];   //[ntrks]
   Int_t           tpcid[336][1000];   //[ntrks]
   Float_t         hit_plane[336][1000];   //[ntrks]
   Float_t         ped[336][1000];   //[ntrks]
   Float_t         amp[336][1000];   //[ntrks]
   Int_t           tamp[336][1000];   //[ntrks]
   Double_t        cosgma[336][1000];   //[ntrks]
   Float_t         noiserms[336][1000];   //[ntrks]
   Float_t         noisermsfit[336][1000];   //[ntrks]

   // List of branches
   TBranch        *b_event;   //!
   TBranch        *b_run;   //!
   TBranch        *b_subrun;   //!
   TBranch        *b_evttime;   //!
   TBranch        *b_year_month_date;   //!
   TBranch        *b_hour_min_sec;   //!
   TBranch        *b_ntrks;   //!
   TBranch        *b_trkid;   //!
   TBranch        *b_trkstart;   //!
   TBranch        *b_trkend;   //!
   TBranch        *b_trklen;   //!
   TBranch        *b_trkthetaxz;   //!
   TBranch        *b_trkthetayz;   //!
   TBranch        *b_trkstartcosxyz;   //!
   TBranch        *b_trkendcosxyz;   //!
   TBranch        *b_ntrkhits;   //!
   TBranch        *b_trkdqdx;   //!
   TBranch        *b_trkx;   //!
   TBranch        *b_trkt;   //!
   TBranch        *b_trkhitx;   //!
   TBranch        *b_trkhity;   //!
   TBranch        *b_trkhitz;   //!
   TBranch        *b_wireid;   //!
   TBranch        *b_chid;   //!
   TBranch        *b_tpcid;   //!
   TBranch        *b_hit_plane;   //!
   TBranch        *b_ped;   //!
   TBranch        *b_amp;   //!
   TBranch        *b_tamp;   //!
   TBranch        *b_cosgma;   //!
   TBranch        *b_noiserms;   //!
   TBranch        *b_noisermsfit;   //!

   ana(TTree *tree=0);
   virtual ~ana();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef ana_cxx
ana::ana(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("signal2noise.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("signal2noise.root");
      }
      TDirectory * dir = (TDirectory*)f->Get("signal2noise.root:/signal2noise");
      dir->GetObject("Event",tree);

   }
   Init(tree);
}

ana::~ana()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t ana::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t ana::LoadTree(Long64_t entry)
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

void ana::Init(TTree *tree)
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

   fChain->SetBranchAddress("event", &event, &b_event);
   fChain->SetBranchAddress("run", &run, &b_run);
   fChain->SetBranchAddress("subrun", &subrun, &b_subrun);
   fChain->SetBranchAddress("evttime", &evttime, &b_evttime);
   fChain->SetBranchAddress("year_month_date", &year_month_date, &b_year_month_date);
   fChain->SetBranchAddress("hour_min_sec", &hour_min_sec, &b_hour_min_sec);
   fChain->SetBranchAddress("ntrks", &ntrks, &b_ntrks);
   fChain->SetBranchAddress("trkid", trkid, &b_trkid);
   fChain->SetBranchAddress("trkstart", trkstart, &b_trkstart);
   fChain->SetBranchAddress("trkend", trkend, &b_trkend);
   fChain->SetBranchAddress("trklen", trklen, &b_trklen);
   fChain->SetBranchAddress("trkthetaxz", trkthetaxz, &b_trkthetaxz);
   fChain->SetBranchAddress("trkthetayz", trkthetayz, &b_trkthetayz);
   fChain->SetBranchAddress("trkstartcosxyz", trkstartcosxyz, &b_trkstartcosxyz);
   fChain->SetBranchAddress("trkendcosxyz", trkendcosxyz, &b_trkendcosxyz);
   fChain->SetBranchAddress("ntrkhits", ntrkhits, &b_ntrkhits);
   fChain->SetBranchAddress("trkdqdx", trkdqdx, &b_trkdqdx);
   fChain->SetBranchAddress("trkx", trkx, &b_trkx);
   fChain->SetBranchAddress("trkt", trkt, &b_trkt);
   fChain->SetBranchAddress("trkhitx", trkhitx, &b_trkhitx);
   fChain->SetBranchAddress("trkhity", trkhity, &b_trkhity);
   fChain->SetBranchAddress("trkhitz", trkhitz, &b_trkhitz);
   fChain->SetBranchAddress("wireid", wireid, &b_wireid);
   fChain->SetBranchAddress("chid", chid, &b_chid);
   fChain->SetBranchAddress("tpcid", tpcid, &b_tpcid);
   fChain->SetBranchAddress("hit_plane", hit_plane, &b_hit_plane);
   fChain->SetBranchAddress("ped", ped, &b_ped);
   fChain->SetBranchAddress("amp", amp, &b_amp);
   fChain->SetBranchAddress("tamp", tamp, &b_tamp);
   fChain->SetBranchAddress("cosgma", cosgma, &b_cosgma);
   fChain->SetBranchAddress("noiserms", noiserms, &b_noiserms);
   fChain->SetBranchAddress("noisermsfit", noisermsfit, &b_noisermsfit);
   Notify();
}

Bool_t ana::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void ana::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t ana::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef ana_cxx
