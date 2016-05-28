// ActivityTriggerService.h
//
// Michael Baird
// May 2016
//
// Basic trigger that looks for some minimal amount of activity in an APA.
// The trigger works by looking for "hits" above some threshold in the
// collection wires.
//

#ifndef ACTIVITYTRIGGERSERVICE_H
#define ACTIVITYTRIGGERSERVICE_H

#include "dune/DuneInterface/BaseTriggerService.h"

#include <memory>
#include <string>
#include <iostream>

namespace fhicl {
  class ParameterSet;
}
namespace art {
  class ActivityRegistry;
}

class ActivityTriggerService : public BaseTriggerService {
  
 public:
  
  // Simple constructor...
  ActivityTriggerService();
  
  // Constructor from fcl parameters...
  ActivityTriggerService(fhicl::ParameterSet const & pset, art::ActivityRegistry&);
  
  // Make a trigger decision based on TPC info only.
  //
  //   params:
  //   raw - A vector of raw digits, to be sliced and diced as the trigger algorithm
  //         sees fit...
  virtual bool TPCTrigger(std::vector<raw::RawDigit> const & raw);

  // Reset function used to reset any global variables
  void Reset();



 private:

  // Parameters.
  unsigned int fMinHits;     ///< minimum number of hits to make a trigger
  short        fThreshold;   ///< threshold for making a hit
  unsigned int fMinTicks;    ///< minimum number of ticks above threshold to be considered a hit
  unsigned int fNHits;       ///< total number of hits counted

};

DECLARE_ART_SERVICE_INTERFACE_IMPL(ActivityTriggerService, BaseTriggerService, LEGACY)

#endif
