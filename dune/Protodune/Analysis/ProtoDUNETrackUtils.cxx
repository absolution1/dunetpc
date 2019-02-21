#include "dune/Protodune/Analysis/ProtoDUNETrackUtils.h"

#include "larsim/MCCheater/BackTrackerService.h"
#include "larsim/MCCheater/ParticleInventoryService.h"
#include "lardataobj/RecoBase/Hit.h"
#include "art/Framework/Principal/Event.h"
#include "canvas/Persistency/Common/FindManyP.h"

#include "TFile.h"
#include "TH1F.h"

protoana::ProtoDUNETrackUtils::ProtoDUNETrackUtils(){

}

protoana::ProtoDUNETrackUtils::~ProtoDUNETrackUtils(){

}

std::vector<anab::CosmicTag> protoana::ProtoDUNETrackUtils::GetRecoTrackCosmicTag(const recob::Track &track, art::Event const &evt, std::string trackModule) const{

  auto recoTracks = evt.getValidHandle<std::vector<recob::Track> >(trackModule);

  unsigned int trackIndex = track.ID();

  // Convert to std::vector<anab::CosmicTag> from std::vector<art::Ptr<anab::CosmicTag>>
  std::vector<anab::CosmicTag> trackTags;

  try{
    const art::FindManyP<anab::CosmicTag> findCosmicTags(recoTracks,evt,trackModule);
    for(unsigned int t = 0; t < findCosmicTags.at(trackIndex).size(); ++t){
      trackTags.push_back((*(findCosmicTags.at(trackIndex)[t])));
    }
  }
  catch(...){
//    std::cerr << "Product not found - returning empty vector" << std::endl;
  }

  return trackTags;    
}

std::vector<anab::T0> protoana::ProtoDUNETrackUtils::GetRecoTrackT0(const recob::Track &track, art::Event const &evt, std::string trackModule) const{

  auto recoTracks = evt.getValidHandle<std::vector<recob::Track> >(trackModule);

  unsigned int trackIndex = track.ID();

  // Convert to std::vector<anab::T0> from std::vector<art::Ptr<anab::T0>>
  std::vector<anab::T0> trackT0s;
  
  try{
    const art::FindManyP<anab::T0> findTrackT0s(recoTracks,evt,trackModule);
    for(unsigned int t = 0; t < findTrackT0s.at(trackIndex).size(); ++t){
      trackT0s.push_back((*(findTrackT0s.at(trackIndex)[t])));
    }
  }
  catch(...){
//    std::cerr << "Product not found - returning empty vector" << std::endl;
  }
  
  return trackT0s;

}

// Get the Calorimetry(s) from a given reco track
std::vector<anab::Calorimetry> protoana::ProtoDUNETrackUtils::GetRecoTrackCalorimetry(const recob::Track &track, art::Event const &evt, const std::string trackModule, const std::string caloModule) const{

  auto recoTracks = evt.getValidHandle<std::vector<recob::Track> >(trackModule);
  std::vector<anab::Calorimetry> caloInfo;
  
  try{
    const art::FindManyP<anab::Calorimetry> findCalorimetry(recoTracks,evt,caloModule);
    std::vector<art::Ptr<anab::Calorimetry>> theseCalos = findCalorimetry.at(track.ID());

    for( auto calo : theseCalos){
      caloInfo.push_back(*calo);
    }
  }
  catch(...){
    std::cerr << "No calorimetry object found... returning empty vector" << std::endl;
  }

  return caloInfo;
}

// Get the hits from a given reco track
const std::vector<const recob::Hit*> protoana::ProtoDUNETrackUtils::GetRecoTrackHits(const recob::Track &track, art::Event const &evt, const std::string trackModule) const{

  auto recoTracks = evt.getValidHandle<std::vector<recob::Track> >(trackModule);
  art::FindManyP<recob::Hit> findHits(recoTracks,evt,trackModule);
  std::vector<art::Ptr<recob::Hit>> inputHits = findHits.at(track.ID());

  std::vector<const recob::Hit*> trackHits;

  for(const art::Ptr<recob::Hit> hit : inputHits){

    trackHits.push_back(hit.get());

  }

  return trackHits;  

}

const std::vector<const recob::Hit*> protoana::ProtoDUNETrackUtils::GetRecoTrackHitsFromPlane(const recob::Track &track, art::Event const &evt, const std::string trackModule, int planeID ) const{

  std::vector<const recob::Hit*> trackHits;
  if( planeID < 0 || planeID > 2 ){
    std::cout << "Please input plane 0, 1, or 2" << std::endl;
    return trackHits;
  }

  auto recoTracks = evt.getValidHandle<std::vector<recob::Track> >(trackModule);
  art::FindManyP<recob::Hit> findHits(recoTracks,evt,trackModule);
  std::vector<art::Ptr<recob::Hit>> inputHits = findHits.at(track.ID());

  for(const art::Ptr<recob::Hit> hit : inputHits){
    //Hacking this because idk how to get the plane id 
    std::string plane_string = hit.get()->WireID().asPlaneID().toString();
    size_t pos = plane_string.find( "P" );
    if ( pos == std::string::npos ) continue;
    std::string planeID_string = "";
    planeID_string += plane_string[ pos + 2 ];
    if( planeID_string != std::to_string( planeID ) ) continue;
       
    std::cout << "Found hit on " << planeID << std::endl;

    trackHits.push_back(hit.get());

  }

  return trackHits;  

}

std::vector< double >  protoana::ProtoDUNETrackUtils::CalibrateCalorimetry(  const recob::Track &track, art::Event const &evt, const std::string trackModule, const std::string caloModule, const fhicl::ParameterSet &ps ) {


  int planeID = ps.get< int >( "PlaneID" );
   
  double betap  = ps.get< double >( "betap"  );
  double Rho    = ps.get< double >( "Rho"    );
  double Efield = ps.get< double >( "Efield" );
  double Wion   = ps.get< double >( "Wion"   );
  double alpha  = ps.get< double >( "alpha"  );
  double norm_factor = ps.get< double >( "norm_factor" );
  double calib_factor = ps.get< double >( "calib_factor" );
  std::string X_correction_name = ps.get< std::string >( "X_correction" );
  TFile X_correction_file = TFile( X_correction_name.c_str(), "OPEN" );
  TH1F * X_correction_hist = (TH1F*)X_correction_file.Get( "dqdx_X_correction_hist" );


  std::vector< double > calibrated_dEdx;

  //Get the Calorimetry vector from the track
  std::vector< anab::Calorimetry > caloVector = GetRecoTrackCalorimetry( track, evt, trackModule, caloModule ); 
  
  size_t calo_position;
  bool found_plane = false;
  for( size_t i = 0; i < caloVector.size(); ++i ){
     //Hacking this because idk how to get the plane id 
     std::string plane_string = caloVector.at(i).PlaneID().toString();
     size_t pos = plane_string.find( "P" );
     if ( pos == std::string::npos ) continue;
     std::string planeID_string = "";
     planeID_string += plane_string[ pos + 2 ];
     if( planeID_string == std::to_string( planeID ) ){
       calo_position = i;
       found_plane = true;
       break;
     }
  }

  if( !found_plane ){
    std::cout << "Could not find the correct plane in the calorimetry vector" << std::endl;
    return calibrated_dEdx;
  }

  std::vector< float > dQdX = caloVector.at( calo_position).dQdx();
  auto theXYZPoints = caloVector.at( calo_position).XYZ();

  //Get the hits from the track from a specific plane
  const std::vector< const recob::Hit* > hits = GetRecoTrackHitsFromPlane( track, evt, trackModule, planeID ); 
  if( hits.size() == 0 ){
    std::cout << "Got empty hits vector" << std::endl;
    return calibrated_dEdx;
  }

  //Do Ajib's correction 
  for( size_t i = 0; i < dQdX.size(); ++i ){ 
    double hit_x = theXYZPoints[i].X();
    int X_bin = X_correction_hist->FindBin( hit_x );
    double X_correction = X_correction_hist->GetBinContent(X_bin);

    double corrected_dq_dx = dQdX[i] * X_correction * norm_factor;
    double scaled_corrected_dq_dx = corrected_dq_dx / calib_factor;
    double cal_de_dx = calc_dEdX( scaled_corrected_dq_dx,  betap,  Rho,  Efield,  Wion,  alpha );
 
    calibrated_dEdx.push_back( cal_de_dx );
  }


  return calibrated_dEdx;
}

double protoana::ProtoDUNETrackUtils::calc_dEdX(double dqdx, double betap, double Rho, double Efield, double Wion, double alpha){
  return (exp(dqdx*(betap/(Rho*Efield)*Wion))-alpha)/(betap/(Rho*Efield));  
}


// Get the hits from a given reco track
unsigned int protoana::ProtoDUNETrackUtils::GetNumberRecoTrackHits(const recob::Track &track, art::Event const &evt, const std::string trackModule) const{

  return GetRecoTrackHits(track,evt,trackModule).size();

}

// Get the PID from a given track
std::vector<anab::ParticleID> protoana::ProtoDUNETrackUtils::GetRecoTrackPID(const recob::Track &track, art::Event const &evt, const std::string trackModule, const std::string pidModule) const{

  auto recoTracks = evt.getValidHandle<std::vector<recob::Track> >(trackModule);
  std::vector<anab::ParticleID> pidvec;

  try{
    const art::FindManyP<anab::ParticleID> findPID(recoTracks,evt,pidModule);
    std::vector<art::Ptr<anab::ParticleID>> thePID = findPID.at(track.ID());
  
    for( auto pid : thePID){
      pidvec.push_back(*pid);
    }
  }
  catch(...){
    std::cerr << "No track PID object found... returning empty vector" << std::endl;
  }

  return pidvec;

}
