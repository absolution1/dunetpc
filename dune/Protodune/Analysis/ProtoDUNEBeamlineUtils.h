#ifndef PROTODUNE_BEAMLINE_UTILS_H
#define PROTODUNE_BEAMLINE_UTILS_H

///////////////////////////////////////////////////////////////
// ProtoDUNEBeamlineUtils
//  - Class to create tracks and momentum from beamline 
//    information. Will be used when updated monitor positions
//    are available.
//   
// - Jake Calcutt (calcuttj@msu.edu) 
///////////////////////////////////////////////////////////////

#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Principal/Event.h"

#include "dunetpc/dune/DuneObj/ProtoDUNEBeamEvent.h"
#include "lardataobj/RecoBase/Track.h"

#include <set>
#include <map>
#include <string>

#include "TVector3.h"

namespace protoana {

  class ProtoDUNEBeamlineUtils {

  public:

    ProtoDUNEBeamlineUtils(fhicl::ParameterSet const &pset);
    ~ProtoDUNEBeamlineUtils();

    void reconfigure(fhicl::ParameterSet const &pset);

    void GetFibers( art::Event const & evt); 
    void GetCurrent( art::Event const & evt);

    std::vector< recob::Track> MakeTracks( art::Event const & evt);

    double   GetPosition( short );
    TVector3 ConvertMonitorCoordinates( double, double, double, double );
    void     BeamMonitorBasisVectors();
    void     RotateMonitorVector(TVector3&);
    TVector3 ProjectToTPC(TVector3, TVector3);

    std::vector< double > MomentumSpec( art::Event const & evt);
    double MomentumCosTheta( double, double, double );

  private:

    art::InputTag fBeamEventTag;

    std::map< std::string, std::vector< short > > ActiveFibers;

    double Current;
    
    //std::vector< recob::Track > Tracks;

    //Just to make it easier to go through all devices
    std::vector< std::string > AllDevices = { 
      "XBPF022707",
      "XBPF022708",
      "XBPF022716",
      "XBPF022717",
      "XBPF022697",
      "XBPF022701",
      "XBPF022702"    
    };

    //These will be static
    //
    //Tracking monitors
    std::string HorizUpstream   =  "XBPF022707";
    std::string VertUpstream    =  "XBPF022708";
    std::string HorizDownstream =  "XBPF022716";
    std::string VertDownstream  =  "XBPF022717";

    //Momentum monitors
    std::string BProf1 =  "XBPF022697";
    std::string BProf2 =  "XBPF022701";
    std::string BProf3 =  "XBPF022702";


    //Basis Vectors for transforming coords
    TVector3 MonitorBasisX;
    TVector3 MonitorBasisY;
    TVector3 MonitorBasisZ;

    bool rotated = false;

    double fRotateMonitorXZ;
    double fRotateMonitorYZ;
    double fRotateMonitorYX;

    double fFirstTrackingProfZ;
    double fSecondTrackingProfZ;
    double fNP04FrontZ;  
    double fBeamX, fBeamY, fBeamZ;


    double fBeamBend;
    double L1, L2, L3;

    //Hardware Parameters for magnetic field stuff
    double mag_P1 = 5.82044830e-3;
//    double mag_P2 = 0.;
    double mag_P3 = -4.68880000e-6;
    double mag_P4 = 324.573967;


  };

}

#endif

