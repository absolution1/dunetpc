#ifndef BEAMDATA_PROTODUNEBEAMSPILL_H
#define BEAMDATA_PROTODUNEBEAMSPILL_H

#include <vector>

namespace beamspill
{
  
  //Fiber Beam Monitor
  //
  struct FBM{

    //Bitmap for hit fibers in the monitor
    std::vector<short int> fibers = std::vector<short int>(24);
      
    long long int timeStamp;
    
    //ID (position in beamline?) of monitor
    int ID;
  };

  //Cerenkov Threshold Detector
  //
  struct CKov{          
    //Status at time of system trigger (on/off)
    bool trigger;
    
    long long int timeStamp;
  };



  class ProtoDUNEBeamSpill{
    public:
      ProtoDUNEBeamSpill();
     ~ProtoDUNEBeamSpill();
    
    private:

      //First time of anything in the spill
      //
      double t0;

      //Set of FBMs
      //
      std::vector< std::vector < FBM > > fiberMonitors;

      //Set of TOF detectors
      //
      std::vector< long long int > TOF1;
      std::vector< long long int > TOF2;

      //Set of Cerenkov detectors
      //
      std::vector< CKov > CKov1;
      std::vector< CKov > CKov2;
  
  };
}

#endif
