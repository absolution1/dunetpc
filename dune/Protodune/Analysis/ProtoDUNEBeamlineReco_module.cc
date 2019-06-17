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
#include "art_root_io/TFileService.h"
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

  //std::cout << "Computed TOF "      << fBeamlineUtils.ComputeTOF( kProton, 2.0 ) << std::endl;
  //std::cout << "Computed Momentum " << fBeamlineUtils.ComputeMomentum( kProton, 130. ) << std::endl;



  //Access the Beam Event
  auto beamHandle = evt.getValidHandle<std::vector<beam::ProtoDUNEBeamEvent>>("beamevent");
  
  std::vector<art::Ptr<beam::ProtoDUNEBeamEvent>> beamVec;
  if( beamHandle.isValid()){
    art::fill_ptr_vector(beamVec, beamHandle);
  }

  const beam::ProtoDUNEBeamEvent & beamEvent = *(beamVec.at(0)); //Should just have one
  /////////////////////////////////////////////////////////////
  
  
  //Check the quality of the event
  std::cout << "Timing Trigger: " << beamEvent.GetTimingTrigger() << std::endl; 
  std::cout << "Is Matched: "     << beamEvent.CheckIsMatched() << std::endl << std::endl;

  if( !fBeamlineUtils.IsGoodBeamlineTrigger( evt ) ){
    std::cout << "Failed quality check" << std::endl;
    return;
  }

  std::cout << "Passed quality check!" << std::endl << std::endl;  
  /////////////////////////////////////////////////////////////
  

  //Access momentum
  const std::vector< double > & momenta = beamEvent.GetRecoBeamMomenta();
  std::cout << "Number of reconstructed momenta: " << momenta.size() << std::endl;

  if( momenta.size() > 0 ) 
    std::cout << "Measured Momentum: " << momenta.at(0) << std::endl;
  ///////////////////////////////////////////////////////////// 



  //Access time of flight
  const std::vector< double > & the_tofs  = beamEvent.GetTOFs();
  const std::vector< int    > & the_chans = beamEvent.GetTOFChans();

  std::cout << "Number of measured TOF: " << the_tofs.size() << std::endl;
  std::cout << "First TOF: "              << beamEvent.GetTOF()         << std::endl;
  std::cout << "First TOF Channel: "      << beamEvent.GetTOFChan()     << std::endl << std::endl;

  std::cout << "All (TOF, Channels): " << std::endl;
  for( size_t i = 0; i < the_tofs.size(); ++i ){
    std::cout << "\t(" << the_tofs[i] << ", " << the_chans[i] << ")" << std::endl;
  }
  std::cout << std::endl;
  /////////////////////////////////////////////////////////////
  

  //Access Cerenkov info
  std::cout << "Cerenkov status, pressure:" << std::endl;
  std::cout << "C0: " << beamEvent.GetCKov0Status() << ", " << beamEvent.GetCKov0Pressure() << std::endl;
  std::cout << "C1: " << beamEvent.GetCKov1Status() << ", " << beamEvent.GetCKov1Pressure() << std::endl << std::endl;
  ///////////////////////////////////////////////////////////// 



  //Access PID
  std::vector< int > pids = fBeamlineUtils.GetPID( beamEvent, 1. );

  std::cout << "Possible particles" << std::endl;

  for( size_t i = 0; i < pids.size(); ++i ){ 
    std::cout << pids[i] << std::endl;
  }
  std::cout << std::endl;

  PossibleParticleCands candidates = fBeamlineUtils.GetPIDCandidates( beamEvent, 1. );
  std::cout << std::left << std::setw(10) << "electron " << candidates.electron << std::endl;
  std::cout << std::left << std::setw(10) << "muon "     << candidates.muon     << std::endl;
  std::cout << std::left << std::setw(10) << "pion "     << candidates.pion     << std::endl;
  std::cout << std::left << std::setw(10) << "kaon "     << candidates.kaon     << std::endl;
  std::cout << std::left << std::setw(10) << "proton "   << candidates.proton   << std::endl << std::endl;

  std::string candidates_string = fBeamlineUtils.GetPIDCandidates( beamEvent, 1. );
  std::cout << candidates_string << std::endl;
  ///////////////////////////////////////////////////////////// 
  

}

void protoana::ProtoDUNEBeamlineReco::endJob() {}
 
DEFINE_ART_MODULE(protoana::ProtoDUNEBeamlineReco)
