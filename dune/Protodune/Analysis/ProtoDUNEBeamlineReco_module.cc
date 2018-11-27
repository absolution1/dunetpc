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

#include "lardataobj/RecoBase/Track.h"


namespace protoana{
  class ProtoDUNEBeamlineReco;
}

class protoana::ProtoDUNEBeamlineReco : public art::EDAnalyzer{
public:
  
  explicit ProtoDUNEBeamlineReco(fhicl::ParameterSet const& pset);
  virtual ~ProtoDUNEBeamlineReco();
  
  virtual void beginJob() override;
  virtual void endJob() override;

  void analyze(art::Event const & evt) override;
  void reconfigure(fhicl::ParameterSet const& pset);
  
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

  std::vector< recob::Track > tracks = fBeamlineUtils.MakeTracks( evt );
  std::cout << "Got " << tracks.size() << " tracks" << std::endl;

  std::vector< double > momenta = fBeamlineUtils.MomentumSpec( evt );
  std::cout << "Got " << momenta.size() << " reconstructed momenta" << std::endl;
  for( size_t i = 0; i < momenta.size(); ++i ){
    std::cout << "\t" << momenta[i] << std::endl;
  }
  std::cout << std::endl;

}

void protoana::ProtoDUNEBeamlineReco::endJob() {}
 
DEFINE_ART_MODULE(protoana::ProtoDUNEBeamlineReco)
