///////////////////////////////////////////////
//
// Class to store protoDUNE beam spill data 
// extracted from IFBeam database
//
// Jake Calcutt (calcuttj@msu.edu) - July 2018
// 
// ////////////////////////////////////////////


#include "dune/DuneObj/ProtoDUNEBeamSpill.h"
#include <algorithm>
#include <iterator>

namespace beam
{
  ProtoDUNEBeamSpill::ProtoDUNEBeamSpill(){}
  void ProtoDUNEBeamSpill::InitFBMs(std::vector<std::string> monitors){
    nFBMs = monitors.size();
    std::cout << "Got " << nFBMs << " monitors and " << t0.size() << " good particles" << std::endl;
    FBM dummyFBM;
    dummyFBM.ID = -1;
    dummyFBM.fibers = {};
    std::uninitialized_fill( std::begin(dummyFBM.fiberData), std::end(dummyFBM.fiberData), 0. );
    std::uninitialized_fill( std::begin(dummyFBM.timeData), std::end(dummyFBM.timeData), 0. );
    dummyFBM.timeStamp = 0.;
    dummyFBM.decoded = false;
    dummyFBM.active = std::vector<short>();

    std::vector<FBM> dummyVec;
    for(size_t i = 0; i < t0.size(); ++i){
      dummyVec.push_back(dummyFBM);
    }
    for(size_t i = 0; i < nFBMs; ++i){
      std::cout << "Initing " << monitors[i] << std::endl;
      fiberMonitors[monitors[i]] = dummyVec;
    }
  }
  ProtoDUNEBeamSpill::~ProtoDUNEBeamSpill(){}
}     
