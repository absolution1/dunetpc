////////////////////////////////////////////////////////////////////////
// Class:       NeutronDecayN2Ana
// Module Type: analyzer
// File:        NeutronDecayN2Ana_module.cc
//
// Generated at Sun Mar 24 09:05:02 2013 by Tingjun Yang using artmod
// from art v1_02_06.
//
//  Having a look at FD reconstruction quantities
//
//  Thomas Karl Warburton
//  k.warburton@sheffield.ac.uk
//
////////////////////////////////////////////////////////////////////////

// Framework includes
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h" 
#include "art/Framework/Principal/Event.h" 
#include "fhiclcpp/ParameterSet.h" 
#include "art/Framework/Principal/Handle.h" 
#include "canvas/Persistency/Common/Ptr.h" 
#include "canvas/Persistency/Common/PtrVector.h" 
#include "art/Framework/Services/Registry/ServiceHandle.h" 
#include "art/Framework/Services/Optional/TFileService.h" 
#include "art/Framework/Services/Optional/TFileDirectory.h" 
#include "messagefacility/MessageLogger/MessageLogger.h" 

// LArSoft includes
#include "larcore/Geometry/Geometry.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardataobj/RawData/ExternalTrigger.h"
#include "lardata/Utilities/AssociationUtil.h"
#include "larsim/MCCheater/BackTracker.h"
#include "lardataobj/AnalysisBase/Calorimetry.h"
#include "lardataobj/AnalysisBase/T0.h"
#include "lardataobj/AnalysisBase/ParticleID.h"

#include "nusimdata/SimulationBase/MCParticle.h"
#include "nusimdata/SimulationBase/MCTruth.h"

// ROOT includes
#include "TTree.h"
#include "TFile.h"

//standard library includes
#include <map>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <cmath>
#include <memory>
#include <limits> // std::numeric_limits<>

namespace NeutronDecayN2Ana {
  class NeutronDecayN2Ana;
}

class NeutronDecayN2Ana::NeutronDecayN2Ana : public art::EDAnalyzer {
public:
  explicit NeutronDecayN2Ana(fhicl::ParameterSet const & p);
  virtual ~NeutronDecayN2Ana();

  void analyze(art::Event const & e) override;

  void beginRun(art::Run& run);
  void beginJob() override;
  void endJob();
  void endRun();
  //void reconfigure(fhicl::ParameterSet const & p) override;
  
private:
  
  void ResetVars();
 
  // Handles
  art::ServiceHandle<geo::Geometry> geom;
  art::ServiceHandle<cheat::BackTracker> bktrk;

  // Parameter List
  std::string fHitsModuleLabel;
  std::string fTrackModuleLabel;
  std::string fMCTruthT0ModuleLabel;

  double ActiveBounds[6]; // Cryostat boundaries ( neg x, pos x, neg y, pos y, neg z, pos z )

  TTree* fRecoTree;
  int MCPdgCode;
  int MCTrackID;

 };
// ********************************** Begin Run *******************************************************
void NeutronDecayN2Ana::NeutronDecayN2Ana::beginRun(art::Run& run) {

}
// *********************************** Begin Job ********************************************************
void NeutronDecayN2Ana::NeutronDecayN2Ana::beginJob()
{
  // Build my Cryostat boundaries array...Taken from Tyler Alion in Geometry Core. Should still return the same values for uBoone.
  ActiveBounds[0] = ActiveBounds[2] = ActiveBounds[4] = DBL_MAX;
  ActiveBounds[1] = ActiveBounds[3] = ActiveBounds[5] = -DBL_MAX;

  std::cout << "The far detector has " << geom->Ncryostats() << " cryostats." << std::endl;

  // assume single cryostats
  auto const* geom = lar::providerFrom<geo::Geometry>();
  for (geo::TPCGeo const& TPC: geom->IterateTPCs()) {
    // get center in world coordinates
    double origin[3] = {0.};
    double center[3] = {0.};
    TPC.LocalToWorld(origin, center);
    double tpcDim[3] = {TPC.HalfWidth(), TPC.HalfHeight(), 0.5*TPC.Length() };
    
    if( center[0] - tpcDim[0] < ActiveBounds[0] ) ActiveBounds[0] = center[0] - tpcDim[0];
    if( center[0] + tpcDim[0] > ActiveBounds[1] ) ActiveBounds[1] = center[0] + tpcDim[0];
    if( center[1] - tpcDim[1] < ActiveBounds[2] ) ActiveBounds[2] = center[1] - tpcDim[1];
    if( center[1] + tpcDim[1] > ActiveBounds[3] ) ActiveBounds[3] = center[1] + tpcDim[1];
    if( center[2] - tpcDim[2] < ActiveBounds[4] ) ActiveBounds[4] = center[2] - tpcDim[2];
    if( center[2] + tpcDim[2] > ActiveBounds[5] ) ActiveBounds[5] = center[2] + tpcDim[2];
  } // for all TPC
  std::cout << "Active Boundaries: "
	    << "\n\tx: " << ActiveBounds[0] << " to " << ActiveBounds[1]
	    << "\n\ty: " << ActiveBounds[2] << " to " << ActiveBounds[3]
	    << "\n\tz: " << ActiveBounds[4] << " to " << ActiveBounds[5]
	    << std::endl;
  
  // Implementation of optional member function here.
  art::ServiceHandle<art::TFileService> tfs;
  fRecoTree = tfs->make<TTree>("ReconstructedTree","analysis tree");
  fRecoTree->Branch("MCPdgCode"        ,&MCPdgCode        ,"MCPdgCode/I"        );
  fRecoTree->Branch("MCTrackID"        ,&MCTrackID        ,"MCTrackID/I"        );
  
}
// ************************************ End Job *********************************************************
void NeutronDecayN2Ana::NeutronDecayN2Ana::endJob() {
  
}
// ************************************ End Run *********************************************************
void NeutronDecayN2Ana::NeutronDecayN2Ana::endRun() {
}

// ********************************** pset param *******************************************************
NeutronDecayN2Ana::NeutronDecayN2Ana::NeutronDecayN2Ana(fhicl::ParameterSet const & pset)
  : EDAnalyzer(pset)
  , fHitsModuleLabel         ( pset.get< std::string >("HitsModuleLabel"))
  , fTrackModuleLabel        ( pset.get< std::string >("TrackModuleLabel"))
  , fMCTruthT0ModuleLabel    ( pset.get< std::string >("MCTruthT0ModuleLabel"))
{

}
// ******************************************************************************************************
NeutronDecayN2Ana::NeutronDecayN2Ana::~NeutronDecayN2Ana()
{
  // Clean up dynamic memory and other resources here.

}
// ************************************ Analyse *********************************************************
void NeutronDecayN2Ana::NeutronDecayN2Ana::analyze(art::Event const & evt)
{
  //std::cout << "\n\n************* New Event / Module running *************\n\n" << std::endl;
  // Implementation of required member function here. 
  art::Handle< std::vector<recob::Track> > trackListHandle;
  std::vector<art::Ptr<recob::Track> > tracklist;
  if (evt.getByLabel(fTrackModuleLabel,trackListHandle))
    art::fill_ptr_vector(tracklist, trackListHandle);
  
  const sim::ParticleList& plist = bktrk->ParticleList();
  // Quickly get total number of MC Particles...
  for ( sim::ParticleList::const_iterator ipar = plist.begin(); ipar!=plist.end(); ++ipar){
    simb::MCParticle *particle = ipar->second;
    for ( int a=0; a<(int)particle->NumberTrajectoryPoints(); ++a ) {
    }
  }
  
  // ------------------ Section for if I want to use reconstructed tracks -------------------
  /*
  if ( trackListHandle.isValid() ) { // Check that trackListHandle is valid.....
    art::FindManyP<recob::Hit>        fmht   (trackListHandle, evt, fTrackModuleLabel);
    art::FindMany<anab::T0>           fmt0   (trackListHandle, evt, fMCTruthT0ModuleLabel);
    int ntracks_reco=tracklist.size();      
    for(int Track=0; Track < ntracks_reco; ++Track){
      if ( fmt0.isValid() ) {
	std::vector<const anab::T0*> T0s = fmt0.at(Track);
	for (size_t t0size =0; t0size < T0s.size(); t0size++) {
	  MCTrackID = T0s[t0size]->TriggerBits();
	} // T0 size
      } else std::cout << "fmt0 isn't valid" << std::endl;
      for ( sim::ParticleList::const_iterator ipar = plist.begin(); ipar!=plist.end(); ++ipar){
	simb::MCParticle *particle = ipar->second;
	if (particle->TrackId() == MCTrackID) {
	  MCPdgCode = particle->PdgCode();
	  std::cout << "Found my truth particle, it has pdg code " << MCPdgCode << std::endl;
	}
      }
      
      fRecoTree->Fill();
    } // Loop through tracks
  } // if trackListHandle.isValid()
  */
  // ------------------ Section for if I want to use reconstructed tracks -------------------
} // Analyse

// ******************************** Reset Variables *****************************************************
void NeutronDecayN2Ana::NeutronDecayN2Ana::ResetVars() {
}
// ******************************** Define Module *****************************************************
DEFINE_ART_MODULE(NeutronDecayN2Ana::NeutronDecayN2Ana)
