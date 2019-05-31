#ifndef PROTODUNE_DATA_UTILS_H
#define PROTODUNE_DATA_UTILS_H

///////////////////////////////////////////////////////////////
// ProtoDUNEDataUtils
//  - Class to help analysers access useful beam data 
//    information
// 
// Leigh Whitehead - leigh.howard.whitehead@cern.ch
///////////////////////////////////////////////////////////////

#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Principal/Event.h"

namespace protoana {

  class ProtoDUNEDataUtils {

  public:

    ProtoDUNEDataUtils(fhicl::ParameterSet const &pset);
    ~ProtoDUNEDataUtils();

    void reconfigure(fhicl::ParameterSet const &pset);

    /**
     * Returns true if the ProtoDUNE trigger says this is a beam trigger
     */
    bool IsBeamTrigger(art::Event const & evt) const;

    /// Get number of active fembs in an APA
    int GetNActiveFembsForAPA(art::Event const & evt, int apa) const;

  private:

    art::InputTag fTimingTag;

  };

}

#endif

