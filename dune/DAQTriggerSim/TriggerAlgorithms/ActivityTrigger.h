// ActivityTrigger.h
//
// Karl Warburton
// May 2017
//
// A activity trigger class for the DAQ trigger framework.
//

#ifndef ActivityTrigger_H
#define ActivityTrigger_H

// Framework includes
#include "lardataobj/RawData/raw.h"
#include "lardataobj/RawData/OpDetWaveform.h"
#include "lardataobj/RawData/RawDigit.h"

// Framework includes for the triggering framework.
#include "dune/DAQTriggerSim/TriggerDataProducts/BasicTrigger.h"
#include "dune/DAQTriggerSim/TriggerDataProducts/TriggerTypes.h"

// C++ includes
#include <vector>
#include <iostream>
#include <string>

namespace triggersim {
  class ActivityTrigger;
}

class triggersim::ActivityTrigger{
 public:

  // The constructor for the ActivityTrigger class.
  ActivityTrigger();

  // Destructor...
  ~ActivityTrigger(){};
  
  // Pass all of the fcl parameters to the class.
  void Configure( fhicl::ParameterSet const& pset );

  // The following is a list of example basic triggering functions which you may want to start from...

  // An example function for how you could trigger using information from the whole art event.
  //    Using this trigger means that you have the whole set of information which you would have
  //     in the art module which you called it from. However, this is useful to separate the 
  //     triggering algorithm from the rest of your code.
  bool TriggerOnWholeEvent( art::Event& event );
  
  // An example function for how you could trigger using information from the TPC Raw Digits.
  //    Using this trigger means that you ony have access to the Raw Digits.
  //      You can find the info on these here: http://nusoft.fnal.gov/larsoft/doxsvn/html/classraw_1_1RawDigit.html
  bool TriggerOnTPC(std::vector< raw::RawDigit> rawTPC);

  // An example function for how you could trigger using information from the Photon detector OpDetWaveforms
  //    Using this trigger means that you ony have access to the OpDetWaveforms
  //      You can find the info on these here: http://nusoft.fnal.gov/larsoft/doxsvn/html/classraw_1_1OpDetWaveform.html
  bool TriggerOnPD( std::vector< raw::OpDetWaveform > rawPD );

  // An example function for how to trigger using both RawDigits and OpDetWaveforms.
  //    Using this trigger will mean that you have access to both RawDigits and OpDetWaveforms.
  //      The methods to access stuff from these can be found in the documentation of the above functions...
  bool TriggerOnTPC_PD( std::vector< raw::RawDigit> rawTPC, std::vector< raw::OpDetWaveform> rawPD);

  // An example function for how to trigger using the output of other triggers.
  //    Using this trigger means that you only have access to what is stored in the triggersim::BasicTrigger data product.
  //      The methods to access stuff from these can be found in: dune/DAQTriggerSim/TriggerDataProducts/BasicTrigger.h
  bool TriggerOnTriggers( std::vector<triggersim::BasicTrigger> triggerVec);

 private:
  
  // Parameters:
  bool         fMakeTrig;        ///< Boolean which is passed as to whether to even attempt this trigger
  bool         fTrigDecision;    ///< did the trigger algorithm decide to keep this event?
  unsigned int fNumber;          ///< Some random number to change.

  std::string fMyString;
  std::string fRawDigLabel;
  std::string fOpWaveLabel;

};

#endif
