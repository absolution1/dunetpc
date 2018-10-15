#ifndef PROTODUNE_DATA_UTILS_H
#define PROTODUNE_DATA_UTILS_H

///////////////////////////////////////////////////////////////
// ProtoDUNEDataUtils
//  - Class to help analysers access useful beam data 
//    information
// 
// Leigh Whitehead - leigh.howard.whitehead@cern.ch
///////////////////////////////////////////////////////////////

#include "art/Framework/Principal/Event.h"

namespace protoana {

  class ProtoDUNEDataUtils {

  public:

    ProtoDUNEDataUtils();
    ~ProtoDUNEDataUtils();

    /// Access trigger information
    bool IsBeamTrigger(art::Event const & evt) const;

  private:




  };

}

#endif

