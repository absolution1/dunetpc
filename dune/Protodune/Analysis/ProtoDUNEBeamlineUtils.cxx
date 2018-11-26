#include "dunetpc/dune/Protodune/Analysis/ProtoDUNEBeamlineUtils.h"

#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"

#include "TVector3.h"

protoana::ProtoDUNEBeamlineUtils::ProtoDUNEBeamlineUtils(fhicl::ParameterSet const& p){
  this->reconfigure(p);
}

protoana::ProtoDUNEBeamlineUtils::~ProtoDUNEBeamlineUtils(){

}

void protoana::ProtoDUNEBeamlineUtils::GetFibers( art::Event const & evt ){
 
  auto beamHandle = evt.getValidHandle<std::vector<beam::ProtoDUNEBeamEvent>>(fBeamEventTag);
  
  std::vector<art::Ptr<beam::ProtoDUNEBeamEvent>> beamVec;
  if( beamHandle.isValid()){
    art::fill_ptr_vector(beamVec, beamHandle);
  }

  //Should just have one
  const beam::ProtoDUNEBeamEvent & beamEvent = *(beamVec.at(0));

  for( size_t i = 0; i < AllDevices.size(); ++i ){
    std::string monitor = AllDevices[i];
    ActiveFibers[ monitor ] = beamEvent.GetActiveFibers( monitor );

    
    if( ActiveFibers.at(monitor).size() > 0){
      std::cout << monitor << " has active fibers: " << std::endl;
      for( size_t j = 0; j < ActiveFibers.at(monitor).size(); ++j ){
        std::cout << ActiveFibers.at(monitor).at(j) << " ";    
      }
      std::cout << std::endl << std::endl;
    }
    else{ 
      std::cout << monitor << " has no active fibers" << std::endl;
    }
  }
}

std::vector< recob::Track > protoana::ProtoDUNEBeamlineUtils::MakeTracks( art::Event const & evt ){

  std::vector< recob::Track > tracks;
 
  //Load fibers 
  GetFibers( evt ); 

  std::vector<short> HorizUpstreamFibers   = ActiveFibers.at(HorizUpstream);
  std::vector<short> HorizDownstreamFibers = ActiveFibers.at(HorizDownstream);
  std::vector<short> VertUpstreamFibers    = ActiveFibers.at(VertUpstream);
  std::vector<short> VertDownstreamFibers  = ActiveFibers.at(VertDownstream);

  //Check that the monitors all have at least one hit each
  bool all_good ;
  all_good = ( ( HorizUpstreamFibers.size() > 0 )   && ( VertUpstreamFibers.size() > 0 ) 
            && ( HorizDownstreamFibers.size() > 0 ) && ( VertDownstreamFibers.size() > 0 ) );

  if( !all_good ){
    std::cout << "At least one empty monitor. Producing no track" << std::endl;
    return tracks;
  }

  //Convention: (Horiz, Vert)
  std::vector< std::pair<short, short> > UpstreamPairedFibers;
  std::vector< std::pair<short, short> > DownstreamPairedFibers;

  std::vector< TVector3 > UpstreamPositions;
  std::vector< TVector3 > DownstreamPositions;

  for(size_t iH = 0; iH < HorizUpstreamFibers.size(); ++iH){

    size_t HorizFiber = HorizUpstreamFibers[iH];

    for(size_t iV = 0; iV < VertUpstreamFibers.size(); ++iV){
      size_t VertFiber = VertUpstreamFibers[iV];

      //LOG_DEBUG("BeamEvent") << "Paired: " << HorizFiber << " " << VertFiber << "\n"; 
      std::cout << "Paired: " << HorizFiber << " " << VertFiber << std::endl; 
      UpstreamPairedFibers.push_back(std::make_pair(HorizFiber, VertFiber));

      //If there's 2 adjacent fibers. Skip the next. Could replace this with averaging the position
      //todo later
      if (iV < VertUpstreamFibers.size() - 1){
        if (VertUpstreamFibers[iV] == (VertUpstreamFibers[iV + 1] - 1)) ++iV;
      }    
    }    

    if (iH < HorizUpstreamFibers.size() - 1){
      if (HorizUpstreamFibers[iH] == (HorizUpstreamFibers[iH + 1] - 1)) ++iH;
    }    
  } 
  
  std::cout << "Upstream " << std::endl;

  for(size_t i = 0; i < UpstreamPairedFibers.size(); ++i){
    
    std::pair<short,short> thePair = UpstreamPairedFibers.at(i);
 
    double xPos = GetPosition(thePair.first);
    double yPos = GetPosition(thePair.second);
    
    std::cout << "normal " << xPos << " " << yPos << std::endl;
    TVector3 posInDet = ConvertMonitorCoordinates(xPos,yPos,0.,fFirstTrackingProfZ);
    std::cout << posInDet.X() << " " << posInDet.Y() << " " << posInDet.Z() << std::endl;
    UpstreamPositions.push_back( posInDet );
  }

  for(size_t iH = 0; iH < HorizDownstreamFibers.size(); ++iH){

    size_t HorizFiber = HorizDownstreamFibers[iH];

    for(size_t iV = 0; iV < VertDownstreamFibers.size(); ++iV){
      size_t VertFiber = VertDownstreamFibers[iV];

      //LOG_DEBUG("BeamEvent") << "Paired: " << HorizFiber << " " << VertFiber << "\n"; 
      std::cout << "Paired: " << HorizFiber << " " << VertFiber << std::endl; 
      DownstreamPairedFibers.push_back(std::make_pair(HorizFiber, VertFiber));

      //If there's 2 adjacent fibers. Skip the next. Could replace this with averaging the position
      //todo later
      if (iV < VertDownstreamFibers.size() - 1){
        if (VertDownstreamFibers[iV] == (VertDownstreamFibers[iV + 1] - 1)) ++iV;
      }    
    }    

    if (iH < HorizDownstreamFibers.size() - 1){
      if (HorizDownstreamFibers[iH] == (HorizDownstreamFibers[iH + 1] - 1)) ++iH;
    }    
  } 
  
  std::cout << "Downstream " << std::endl;

  for(size_t i = 0; i < DownstreamPairedFibers.size(); ++i){
    
    std::pair<short,short> thePair = DownstreamPairedFibers.at(i);
 
    double xPos = GetPosition(thePair.first);
    double yPos = GetPosition(thePair.second);
    
    std::cout << "normal " << xPos << " " << yPos << std::endl;
    TVector3 posInDet = ConvertMonitorCoordinates(xPos,yPos,0.,fFirstTrackingProfZ);
    std::cout << posInDet.X() << " " << posInDet.Y() << " " << posInDet.Z() << std::endl;
    DownstreamPositions.push_back( posInDet );
  }

  //Just for creating tracks
  std::vector< std::vector<double> > dummy;
  std::vector<double> mom(3, util::kBogusD);
  ///  

  for(size_t iU = 0; iU < UpstreamPositions.size(); ++iU){
    for(size_t iD = 0; iD < DownstreamPositions.size(); ++iD){
      std::vector<TVector3> thePoints;
      thePoints.push_back(UpstreamPositions.at(iU));
      thePoints.push_back(DownstreamPositions.at(iD));

      //Now project the last point to the TPC face
      thePoints.push_back( ProjectToTPC(thePoints[0],thePoints[1]) );    
      std::cout << "Projected: " << thePoints.back().X() << " " << thePoints.back().Y() << " " << thePoints.back().Z() << std::endl;

     
      std::vector<TVector3> theMomenta;
      //Just push back the unit vector for each point 
      theMomenta.push_back( ( DownstreamPositions.at(iD) - UpstreamPositions.at(iU) ).Unit() );
      theMomenta.push_back( ( DownstreamPositions.at(iD) - UpstreamPositions.at(iU) ).Unit() );
      theMomenta.push_back( ( DownstreamPositions.at(iD) - UpstreamPositions.at(iU) ).Unit() );

      recob::Track tempTrack(thePoints, theMomenta, dummy, mom, 1);     
      tracks.push_back( tempTrack );
    }    
  }
  
  return tracks;
}

void protoana::ProtoDUNEBeamlineUtils::reconfigure(fhicl::ParameterSet const& p){
  fBeamEventTag = p.get<art::InputTag>("BeamEventTag");

  fBeamX = p.get<double>("BeamX");
  fBeamY = p.get<double>("BeamY");
  fBeamZ = p.get<double>("BeamZ");

  fRotateMonitorXZ = p.get<double>("RotateMonitorXZ"); 
  fRotateMonitorYZ = p.get<double>("RotateMonitorYZ"); 

  fFirstTrackingProfZ  = p.get<double>("FirstTrackingProfZ");
  fSecondTrackingProfZ = p.get<double>("SecondTrackingProfZ");
  fNP04FrontZ          = p.get<double>("NP04FrontZ");

}

double protoana::ProtoDUNEBeamlineUtils::GetPosition( short theFiber ){
  //Fibers are 1mm in width
  return 1.*(96 - theFiber) - .5;
}

TVector3 protoana::ProtoDUNEBeamlineUtils::ConvertMonitorCoordinates(double x, double y, double z, double zOffset){

  if( !rotated ) BeamMonitorBasisVectors();

  double off = fNP04FrontZ - zOffset;

  TVector3 old(x,y,z);

  double newX = x*MonitorBasisX.X() + y*MonitorBasisY.X() + off*fabs(MonitorBasisZ.X());
  double newY = x*MonitorBasisX.Y() + y*MonitorBasisY.Y() + off*fabs(MonitorBasisZ.Y());
  double newZ = x*MonitorBasisX.Z() + y*MonitorBasisY.Z() - off*fabs(MonitorBasisZ.Z());

  newX += fBeamX*10.;
  newY += fBeamY*10.;
  newZ += fBeamZ*10.;

  TVector3 result(newX/10., newY/10., newZ/10.);
  return result;
}

void protoana::ProtoDUNEBeamlineUtils::BeamMonitorBasisVectors(){
  MonitorBasisX = TVector3(1.,0.,0.);
  MonitorBasisY = TVector3(0.,1.,0.);
  MonitorBasisZ = TVector3(0.,0.,1.);
  RotateMonitorVector(MonitorBasisX);
  RotateMonitorVector(MonitorBasisY);
  RotateMonitorVector(MonitorBasisZ);

  rotated = true;
}

void protoana::ProtoDUNEBeamlineUtils::RotateMonitorVector(TVector3 &vec){
  vec.RotateY(fRotateMonitorXZ * TMath::Pi()/180.);
  vec.RotateX(fRotateMonitorYZ * TMath::Pi()/180.);
}


TVector3 protoana::ProtoDUNEBeamlineUtils::ProjectToTPC(TVector3 firstPoint, TVector3 secondPoint){
  TVector3 dR = (secondPoint - firstPoint);

  double deltaZ = -1.*secondPoint.Z();
  double deltaX = deltaZ * (dR.X() / dR.Z());
  double deltaY = deltaZ * (dR.Y() / dR.Z());

  TVector3 lastPoint = secondPoint + TVector3(deltaX, deltaY, deltaZ);
  return lastPoint;
}

