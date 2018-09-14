///////////////////////////////////////////////
//
// Class to store protoDUNE beam spill data 
// extracted from IFBeam database
//
// Jake Calcutt (calcuttj@msu.edu) - July 2018
// 
// ////////////////////////////////////////////


#include "dune/DuneObj/ProtoDUNEBeamEvent.h"

namespace beam
{
  ProtoDUNEBeamEvent::ProtoDUNEBeamEvent(){}
  void ProtoDUNEBeamEvent::InitFBMs(std::vector<std::string> monitors){
    nFBMs = monitors.size();
    std::cout << "Got " << nFBMs << " monitors and " << t0.size() << " good particles" << std::endl;
    FBM dummyFBM;
    dummyFBM.ID = -1;
    std::vector<FBM> dummyVec;
    for(size_t i = 0; i < t0.size(); ++i){
      dummyVec.push_back(dummyFBM);
    }
    for(size_t i = 0; i < nFBMs; ++i){
      std::cout << "Initing " << monitors[i] << std::endl;
      fiberMonitors[monitors[i]] = dummyVec;
    }
  }
  ProtoDUNEBeamEvent::~ProtoDUNEBeamEvent(){}
}     
