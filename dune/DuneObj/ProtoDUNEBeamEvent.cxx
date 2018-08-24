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
  void ProtoDUNEBeamEvent::InitFBMs(size_t nMonitors){
    nFBMs = nMonitors;
    FBM dummyFBM;
    dummyFBM.ID = -1;
    std::vector<FBM> dummyVec;
    dummyVec.push_back(dummyFBM);
    for(size_t i = 0; i < nFBMs; ++i){
      fiberMonitors.push_back(dummyVec);
    }
  }
  ProtoDUNEBeamEvent::~ProtoDUNEBeamEvent(){}
}     
