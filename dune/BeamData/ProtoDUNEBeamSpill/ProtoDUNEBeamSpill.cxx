///////////////////////////////////////////////
//
// Class to store protoDUNE beam spill data 
// extracted from the CERN database
//
// Jake Calcutt (calcuttj@msu.edu) - July 2018
// 
// ////////////////////////////////////////////


#include "dune/BeamData/ProtoDUNEBeamSpill/ProtoDUNEBeamSpill.h"

namespace beamspill 
{
  ProtoDUNEBeamSpill::ProtoDUNEBeamSpill(){}
  void ProtoDUNEBeamSpill::InitFBMs(size_t nMonitors){
    nFBMs = nMonitors;
    FBM dummyFBM;
    dummyFBM.ID = -1;
    std::vector<FBM> dummyVec;
    dummyVec.push_back(dummyFBM);
    for(size_t i = 0; i < nFBMs; ++i){
      fiberMonitors.push_back(dummyVec);
    }
  }
  ProtoDUNEBeamSpill::~ProtoDUNEBeamSpill(){}
}     
