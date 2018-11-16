////////////////////////////////////////////////////////////////////////
//
// ProtoDUNEBeamlineFilter selects events based on reconstructed 
// beamline info
//
// 2018 Justin Hugon, jhugon@fnal.gov
//
////////////////////////////////////////////////////////////////////////

// Framework includes
#include "art/Framework/Core/EDFilter.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "fhiclcpp/ParameterSet.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "nusimdata/SimulationBase/MCTruth.h"

// dunetpc includes
#include "dunetpc/dune/Protodune/Analysis/ProtoDUNEDataUtils.h"

// ROOT includes
#include "TTree.h"
#include "TFile.h"
#include "TString.h"

// C++ Includes
#include <fstream>
#include <string>
#include <sstream>
#include <cmath>
#include <algorithm>
#include <iostream>
#include <vector>

const int NMAXBEAMPARTICLES = 500;

namespace protoana{
  class ProtoDUNEBeamlineFilter;
}

class protoana::ProtoDUNEBeamlineFilter : public art::EDFilter{
public:
  
  explicit ProtoDUNEBeamlineFilter(fhicl::ParameterSet const& pset);
  virtual ~ProtoDUNEBeamlineFilter();
  
  void beginJob();
  bool filter(art::Event& evt);
  void reconfigure(fhicl::ParameterSet const& pset);
  void endJob();
  
private:
  
  protoana::ProtoDUNEDataUtils fDataUtils;
  float fNominalBeamMomentum; // GeV/c
  bool fIsElectron;
  bool fIsMuon;
  bool fIsPion;
  bool fIsKaon;
  bool fIsProton;
  bool fAndParticles; // otherwise do or of the allowed particles
};
  
//-----------------------------------------------------------------------
protoana::ProtoDUNEBeamlineFilter::ProtoDUNEBeamlineFilter(fhicl::ParameterSet const& pset):
  fDataUtils(pset.get<fhicl::ParameterSet>("DataUtils"))
{

  this->reconfigure(pset);
  
}

//-----------------------------------------------------------------------
protoana::ProtoDUNEBeamlineFilter::~ProtoDUNEBeamlineFilter(){}

//-----------------------------------------------------------------------
void protoana::ProtoDUNEBeamlineFilter::beginJob() {
  
  //art::ServiceHandle<art::TFileService> tfs;
  
}

//-----------------------------------------------------------------------
void protoana::ProtoDUNEBeamlineFilter::reconfigure(fhicl::ParameterSet const& pset){
  fNominalBeamMomentum = pset.get<float>("NominalBeamMomentum"); // GeV/c
  fIsElectron = pset.get<bool>("IsElectron");
  fIsMuon = pset.get<bool>("IsMuon");
  fIsPion = pset.get<bool>("IsPion");
  fIsKaon = pset.get<bool>("IsKaon");
  fIsProton = pset.get<bool>("IsProton");
  fAndParticles = pset.get<bool>("AndParticles");
}

//-----------------------------------------------------------------------
bool protoana::ProtoDUNEBeamlineFilter::filter(art::Event& evt){

  const auto possibleParts = fDataUtils.GetBeamlineParticleID(evt,fNominalBeamMomentum);
  if(fAndParticles)
  {
    if(fIsElectron && !possibleParts.electron) return false;
    if(fIsMuon && !possibleParts.muon) return false;
    if(fIsPion && !possibleParts.pion) return false;
    if(fIsKaon && !possibleParts.kaon) return false;
    if(fIsProton && !possibleParts.proton) return false;
    return true;
  }
  else // or
  {
    if(fIsElectron && possibleParts.electron) return true;
    if(fIsMuon && possibleParts.muon) return true;
    if(fIsPion && possibleParts.pion) return true;
    if(fIsKaon && possibleParts.kaon) return true;
    if(fIsProton && possibleParts.proton) return true;
    return false;
  }
}

void protoana::ProtoDUNEBeamlineFilter::endJob() {}
 
DEFINE_ART_MODULE(protoana::ProtoDUNEBeamlineFilter)
