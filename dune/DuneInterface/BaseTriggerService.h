// BaseTriggerService.h
//
// Michael Baird
// Feb. 2016
//
// Base trigger service class for other trigger services to inherit from.

#ifndef BaseTriggerService_H
#define BaseTriggerService_H

#include <memory>
#include "larcoreobj/SimpleTypesAndConstants/RawTypes.h"
#include "art/Framework/Services/Registry/ServiceMacros.h"
#include "art/Framework/Principal/Event.h"
#include "lardataobj/RawData/RawDigit.h"
#include "lardataobj/RawData/OpDetWaveform.h"

//jpd -- needed to allow for triggering on triggers
#include "dune/DAQTriggerSim/TriggerDataProducts/BasicTrigger.h"

namespace fhicl {
  class ParameterSet;
}

namespace art {
  class ActivityRegistry;
}


class BaseTriggerService {

public:

  // Destructor:
  virtual ~BaseTriggerService() = default;

  // Make your trigger decision.
  // Can be overridden by the inherited class.
  //
  //   params:
  //   evt - The ART event, to be sliced and diced as the trigger algorithm
  //         sees fit...
  virtual bool Trigger(art::Event const & evt)
  { return false; }

  // Make your trigger decision based on TPC info only.
  // Can be overridden by the inherited class.
  //
  //   params:
  //   rawTPC - A vector of raw digits, to be sliced and diced as the trigger algorithm
  //            sees fit...
  virtual bool TPCTrigger(std::vector<raw::RawDigit> const & rawTPC)
  { return false; }

  // Make your trigger decision based on optical (PD) info only.
  // Can be overridden by the inherited class.
  //
  //   params:
  //   rawPD - A vector of optical waveforms, to be sliced and diced as the trigger algorithm
  //           sees fit...
  virtual bool PDTrigger(std::vector<raw::OpDetWaveform> const & rawPD)
  { return false; }

  // Make your trigger decision based on a mixture of TPC and PD info.
  // Can be overridden by the inherited class.
  //
  //   params:
  //   rawTPC - A vector of raw digits, to be sliced and diced as the trigger algorithm
  //            sees fit...
  //   rawPD  - A vector of optical waveforms, to be sliced and diced as the trigger algorithm
  //            sees fit...
  virtual bool TPCPDTrigger(std::vector<raw::RawDigit>      const & rawTPC,
			    std::vector<raw::OpDetWaveform> const & rawPD)
  { return false; }


  // Make your trigger decision based on a vector of BasicTrigger objects
  // The idea is that this method triggers on more basic triggers that were run earlier in the chain
  // Should be overridden by the inherited class
  //
  //    params:
  //    triggerVec - A vector of BasicTrigger data products, which will have been produced by a lower 
  //                 level trigger. 
  virtual bool Trigger(std::vector<triggersim::BasicTrigger> const &triggerVec)
  { return false; }

  virtual std::string GetName() = 0;

};

DECLARE_ART_SERVICE_INTERFACE(BaseTriggerService, LEGACY)

#endif
