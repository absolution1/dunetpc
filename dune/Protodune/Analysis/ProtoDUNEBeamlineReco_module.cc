////////////////////////////////////////////////////////////////////////
//
// ProtoDUNEBeamlineReco does tracking and reconstruction from the beamline
// monitors
//
// 2018 Jake Calcutt, calcuttj@msu.edu 
//
////////////////////////////////////////////////////////////////////////

// Framework includes
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "fhiclcpp/ParameterSet.h"

// dunetpc includes
#include "dunetpc/dune/Protodune/Analysis/ProtoDUNEDataUtils.h"
#include "dunetpc/dune/Protodune/Analysis/ProtoDUNEBeamlineUtils.h"


namespace protoana{
  class ProtoDUNEBeamlineReco;
}

class protoana::ProtoDUNEBeamlineReco : public art::EDAnalyzer{
public:
  
  explicit ProtoDUNEBeamlineReco(fhicl::ParameterSet const& pset);
  virtual ~ProtoDUNEBeamlineReco();
  
  void beginJob();
  void analyze(art::Event const & evt) override;
  void reconfigure(fhicl::ParameterSet const& pset);
  void endJob();
  
private:
  
  protoana::ProtoDUNEBeamlineUtils fBeamlineUtils;

};
  
//-----------------------------------------------------------------------
protoana::ProtoDUNEBeamlineReco::ProtoDUNEBeamlineReco(fhicl::ParameterSet const& pset):
  EDAnalyzer(pset),
  fBeamlineUtils(pset.get<fhicl::ParameterSet>("BeamlineUtils"))
{
  this->reconfigure(pset);
}

//-----------------------------------------------------------------------
protoana::ProtoDUNEBeamlineReco::~ProtoDUNEBeamlineReco(){}

//-----------------------------------------------------------------------
void protoana::ProtoDUNEBeamlineReco::beginJob() {
  
  /*art::ServiceHandle<art::TFileService> tfs;
  fIsBeamTrigger = tfs->make<TH1F>("IsBeamTrigger", "Is the CTB trigger from the beamline?", 2,0,1);
  */
}

//-----------------------------------------------------------------------
void protoana::ProtoDUNEBeamlineReco::reconfigure(fhicl::ParameterSet const& pset){
/*
  fNominalBeamMomentum = pset.get<float>("NominalBeamMomentum"); // GeV/c
*/
}

//-----------------------------------------------------------------------
void protoana::ProtoDUNEBeamlineReco::analyze(art::Event const & evt){
  fBeamlineUtils.GetFibers( evt );
}

void protoana::ProtoDUNEBeamlineReco::endJob() {}
 
DEFINE_ART_MODULE(protoana::ProtoDUNEBeamlineReco)
