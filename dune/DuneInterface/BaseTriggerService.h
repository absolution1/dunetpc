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
#include "lardataobj/RawData/RawDigit.h"
#include "lardataobj/RawData/OpDetWaveform.h"

//jpd -- needed to allow for triggering on triggers
#include "dune/DAQTriggerSim/TriggerDataProducts/BasicTrigger.h"
//#ifndef __CLING__
#include "art/Framework/Principal/Event.h"
//#endif

namespace fhicl {
  class ParameterSet;
}

namespace art {
  class ActivityRegistry;
  class Event;
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
  virtual bool Trigger(const art::Event& evt)
  { return false; }

  // Make your trigger decision based on TPC info only.
  // Can be overridden by the inherited class.
  //
  //   params:
  //   rawTPC - A vector of raw digits, to be sliced and diced as the trigger algorithm
  //            sees fit...
  virtual bool TPCTrigger(const std::vector<raw::RawDigit>& rawTPC)
  { return false; }

  // Make your trigger decision based on optical (PD) info only.
  // Can be overridden by the inherited class.
  //
  //   params:
  //   rawPD - A vector of optical waveforms, to be sliced and diced as the trigger algorithm
  //           sees fit...
  virtual bool PDTrigger(const std::vector<raw::OpDetWaveform>& rawPD)
  { return false; }

  // Make your trigger decision based on a mixture of TPC and PD info.
  // Can be overridden by the inherited class.
  //
  //   params:
  //   rawTPC - A vector of raw digits, to be sliced and diced as the trigger algorithm
  //            sees fit...
  //   rawPD  - A vector of optical waveforms, to be sliced and diced as the trigger algorithm
  //            sees fit...
  virtual bool TPCPDTrigger(const std::vector<raw::RawDigit>&      rawTPC,
			    const std::vector<raw::OpDetWaveform>& rawPD)
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

#ifndef __CLING__
#include "art/Framework/Services/Registry/ServiceMacros.h"
DECLARE_ART_SERVICE_INTERFACE(BaseTriggerService, LEGACY)
#endif

#endif
