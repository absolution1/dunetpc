////////////////////////////////////////////////////////////////////////////////////////////////
// Class:       RegCVNAna
// Module Type: analyzer
// File:        RegCVNAna_module.cc
// Author:      
////////////////////////////////////////////////////////////////////////////////////////////////

#include <string> 
#include <vector>
#include "TTree.h"
#include "TBranch.h"

// Framework includes: 
#include "art/Framework/Core/ModuleMacros.h"  
#include "art/Framework/Principal/Event.h" 
#include "art/Framework/Principal/Handle.h"
#include "fhiclcpp/ParameterSet.h" 
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Optional/TFileDirectory.h"
#include "art/Framework/Core/EDAnalyzer.h"

#include "larsim/MCCheater/ParticleInventoryService.h"
#include "lardataobj/RecoBase/Hit.h"
#include "dune/FDSensOpt/FDSensOptData/EnergyRecoOutput.h"
#include "dune/RegCVN/func/RegCVNResult.h"

const int kMax = 1000;
namespace myana {
  class RegCVNAna;
}

class myana::RegCVNAna : public art::EDAnalyzer
{
  public:
    RegCVNAna(fhicl::ParameterSet const& pset);
    void analyze(art::Event const& evt);
    void reconfigure(fhicl::ParameterSet const& p);
    void reset();
  
  private:
    TTree* fTree; 
    int ievt;
    double trueEnergy;
    float  regcvn_energy;
    
    int    inDet;
    int    nhits;
    int    nu_truth_N;
    int    nupdg_truth[kMax];
    int    numode_truth[kMax];    
    int    nuccnc_truth[kMax];     
    double nueng_truth[kMax];

    double ErecoNue;
    double RecoLepEnNue; 
    double RecoHadEnNue; 
    int    RecoMethodNue; 

    std::string fMCGenModuleLabel;   
    std::string fRegCVNModuleLabel;   
    std::string fHitsModuleLabel;
    std::string fEnergyRecoNueLabel;

    art::ServiceHandle<art::TFileService> tfs;

};

myana::RegCVNAna::RegCVNAna(fhicl::ParameterSet const& pset) : EDAnalyzer(pset)
{
  this->reconfigure(pset);
  fTree = tfs->make<TTree>("anatree", "anatree");
  fTree->Branch("ievent",            &ievt,          "ievent/I");
  fTree->Branch("InDet",             &inDet,         "InDet/I");
  fTree->Branch("TrueEnergy",        &trueEnergy,    "TrueEnergy/D");
  fTree->Branch("CVNEnergy",         &regcvn_energy, "CVNEnergy/F");
  fTree->Branch("NuTruthN",          &nu_truth_N,    "NuTruthN/I");  
  fTree->Branch("NuEngTruth",        nueng_truth,    "NuEngTruth[NuTruthN]/D");   
  fTree->Branch("NuPDGTruth",        nupdg_truth,    "NuPDGTruth[NuTruthN]/I");   
  fTree->Branch("NuModeTruth",       numode_truth,   "NuModeTruth[NuTruthN]/I"); 
  fTree->Branch("NuCCNCTruth",       nuccnc_truth,   "NuCCNCTruth[NuTruthN]/I");   

  fTree->Branch("NHits",             &nhits,         "NHits/I");

  fTree->Branch("ErecoNue",          &ErecoNue,      "ErecoNue/D");
  fTree->Branch("RecoLepEnNue",      &RecoLepEnNue,  "RecoLepEnNue/D");
  fTree->Branch("RecoHadEnNue",      &RecoHadEnNue,  "RecoHadEnNue/D");
  fTree->Branch("RecoMethodNue",     &RecoMethodNue, "RecoMethodNue/I");


}

void myana::RegCVNAna::reconfigure(fhicl::ParameterSet const& pset)
{
  fMCGenModuleLabel  = pset.get<std::string>("MCGenModuleLabel");
  fRegCVNModuleLabel = pset.get<std::string>("RegCVNModuleLabel");
  fHitsModuleLabel   = pset.get<std::string>("HitsModuleLabel");
  fEnergyRecoNueLabel = pset.get<std::string>("EnergyRecoNueLabel");
}

void myana::RegCVNAna::analyze(art::Event const& evt)
{
  this->reset();
  ievt = evt.id().event();
  bool isMC = !evt.isRealData(); 

  // * MC truth information 
  art::Handle< std::vector<simb::MCTruth> > mctruthListHandle;             
  std::vector<art::Ptr<simb::MCTruth> > mclist;      
  if (isMC){    
    if (evt.getByLabel(fMCGenModuleLabel,mctruthListHandle))         
      art::fill_ptr_vector(mclist, mctruthListHandle);  
  }  

  // Get the hits out of the event record
  art::Handle<std::vector<recob::Hit> > hitHandle;
  std::vector<art::Ptr<recob::Hit> > hits;
  if (evt.getByLabel(fHitsModuleLabel,hitHandle))
    art::fill_ptr_vector(hits, hitHandle);

  // Get DUNE energy Reco
  art::Handle<dune::EnergyRecoOutput> engrecoHandle;
  evt.getByLabel(fEnergyRecoNueLabel,engrecoHandle);


  // Get RegCVN Results
  art::Handle<std::vector<cvn::RegCVNResult>> cvnresultListHandle;
  evt.getByLabel(fRegCVNModuleLabel, "regcvnresult", cvnresultListHandle);
  //std::vector<art::Ptr<cvn::Result> > cvnlist;
  //if (evt.getByLabel(fRegCVNModuleLabel, cvnresultListHandle))
  //  art::fill_ptr_vector(cvnlist, cvnresultListHandle);


  // Get Truth information
  if (mclist.size()>0)
  {
        int neutrino_i = 0;   
        for(size_t iList = 0; (iList < mclist.size()) && (neutrino_i < kMax) ; ++iList)
        {
          if (mclist[iList]->NeutrinoSet())
          {
            nueng_truth[neutrino_i]  = mclist[iList]->GetNeutrino().Nu().Momentum().E();
            nupdg_truth[neutrino_i]  = mclist[iList]->GetNeutrino().Nu().PdgCode(); 
            nuccnc_truth[neutrino_i] = mclist[iList]->GetNeutrino().CCNC(); 
            numode_truth[neutrino_i] = mclist[iList]->GetNeutrino().Mode();
            neutrino_i++; 
          }        
        }   
        nu_truth_N = neutrino_i;  
  }
  // Get Hit information
  nhits = hits.size();
  inDet = 0;
  if (hits.size()>0)
  {
      bool fid_flag = false;
      for(size_t iHit = 0; iHit < hits.size(); ++iHit)
      {
        art::Ptr<recob::Hit> hit = hits.at(iHit);
    	float peakT = hit->PeakTime();
    	unsigned int channel = hit->Channel();
	unsigned int plane = hit->WireID().Plane;
    	if (peakT > 4482.) fid_flag = true; 
        if (plane == 2 && channel < 1605) fid_flag = true;
	if (fid_flag) break; 
      }
      if (!fid_flag) inDet = 1; 
  }

  // Get RecoE from DUNE
  if (!engrecoHandle.failedToGet())
  {
      ErecoNue          = engrecoHandle->fNuLorentzVector.E();
      RecoLepEnNue      = engrecoHandle->fLepLorentzVector.E();
      RecoHadEnNue      = engrecoHandle->fHadLorentzVector.E();
      RecoMethodNue     = engrecoHandle->recoMethodUsed;
      std::cout<< ErecoNue << std::endl;
  }

  // Get RegCVN Results
  if (!cvnresultListHandle.failedToGet())
  {
	if (!cvnresultListHandle->empty())
	{
  	  const std::vector<float>& v = (*cvnresultListHandle)[0].fOutput;
	  regcvn_energy = v[0];
          //std::cout << v[0] << std::endl;
        }
  }
  // fill entry
  fTree->Fill();
}

void myana::RegCVNAna::reset()
{
    ievt = -9999;
    trueEnergy = -99999;
    regcvn_energy = -99999;
    ErecoNue = -99999;
    RecoLepEnNue = -99999;
    RecoHadEnNue = -99999;
    RecoMethodNue = -99999;
    
    nu_truth_N = 0;
    for (int ii = 0; ii < kMax; ++ii)
    {
        nupdg_truth[ii] = -99999;
        numode_truth[ii] = -99999;
        nuccnc_truth[ii] = -99999;
        nueng_truth[ii] = -99999;
    }

}

DEFINE_ART_MODULE(myana::RegCVNAna)
