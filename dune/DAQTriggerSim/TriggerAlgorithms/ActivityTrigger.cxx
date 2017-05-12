// ActivityTrigger.cxx
//
// Karl Warburton --- karlwarb@iastate.edu
// May 2017
//
// A activity trigger class for the DAQ trigger framework.
//

#ifndef ActivityTrigger_CXX
#define ActivityTrigger_CXX

// Framework includes
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "larcore/Geometry/Geometry.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"
#include "larcoreobj/SimpleTypesAndConstants/RawTypes.h"

// DUNETPC specific includes
#include "dune/DAQTriggerSim/TriggerAlgorithms/ActivityTrigger.h"

// C++ includes
#include <iomanip>

//----------------------------------------------------------------------
//------------- The constructor for this trigger algorithm -------------
//----------------------------------------------------------------------
triggersim::ActivityTrigger::ActivityTrigger() {}

//----------------------------------------------------------------------
//------------ The configuration for this trigger algorithm ------------
//----------------------------------------------------------------------
void triggersim::ActivityTrigger::Configure( fhicl::ParameterSet const& pset ) {
  fMakeTrig    = true;
  fMyString    = pset.get<std::string> ("AString");
  fRawDigLabel = pset.get<std::string> ("RawDigitLabel");
  fOpWaveLabel = pset.get<std::string> ("OpDetWaveLabel");

  // --- We have got all of the fcl parameters here too, lets check that they are what we expect...
  std::cout << "\n------In my trigger class------\nThe fcl params have been set to :" 
	    << "\n  fMyString:    " << fMyString 
	    << "\n  fRawDigLabel: " << fRawDigLabel
	    << "\n  fOpWaveLabel: " << fOpWaveLabel
	    << "\n-------------------------------\n"
	    << std::endl;
} // Configure

//----------------------------------------------------------------------
//------------ The trigger algorithm on the whole art event ------------
//----------------------------------------------------------------------
bool triggersim::ActivityTrigger::TriggerOnWholeEvent( art::Event& event ) {
  // --- If for some reason you changed your mind about making this trigger...
  if (!fMakeTrig) return false;
  
  // --- Now do stuff...
  fNumber = event.event();
  std::cout << "  I have just got into my activity trigger for event " << event.event() << ", fMakeTrig is " << fMakeTrig << ", and fNumber is now " << fNumber << std::endl;
  
  // --- Lift out the TPC raw digits:
  auto digits = event.getValidHandle<std::vector<raw::RawDigit> >(fRawDigLabel);
  if ( digits.failedToGet() )
    mf::LogError("ActivityTrigger_Algorithm") << "The raw::RawDigit you gave me " << fRawDigLabel << " is not in the event..." << std::endl;
  
  // --- Lift out the Photon Detector OpDetWaveforms:
  auto waveforms = event.getValidHandle<std::vector<raw::OpDetWaveform> >(fOpWaveLabel);
  if ( waveforms.failedToGet() )
    mf::LogError("ActivityTrigger_Algorithm") << "The raw::OpDetWaveform you gave me " << fOpWaveLabel << " is not in the event..." << std::endl;
  
  // --- We have got the data products here too, look at the functions below to find out how to use them...
  std::cout << "  The size of digits is " << digits->size() << ", and the size of wave is " << waveforms->size() << std::endl;

  // --- Make a decision based on the event number...
  if (event.event() < 3) {
    fTrigDecision = true;
  } else {
    fTrigDecision = false;
  }

  // --- Return the result of the trigger.
  return fTrigDecision;
} // TriggerOnWholeEvent

//----------------------------------------------------------------------
//----------- The trigger algorithm on just the TPC RawDigits ----------
//----------------------------------------------------------------------
bool triggersim::ActivityTrigger::TriggerOnTPC( std::vector< raw::RawDigit> rawTPC ) {
  // --- If for some reason you changed your mind about making this trigger...
  if (!fMakeTrig) return false;

  // --- It is often helpful to make a geometry service handle.
  art::ServiceHandle<geo::Geometry> geom;
  
  // --- Now do stuff...
  fNumber = rawTPC.size();
  for (unsigned int Dig=0; Dig < rawTPC.size(); ++Dig) {
    raw::RawDigit ThisDig = rawTPC[Dig]; // Can also use rawTPC.at(Dig);
    int Chan = rawTPC[Dig].Channel();    // Alternatively can access stuff by doing rawTPC[Dig]->Channel()
    if (Dig < 5) {
      std::cout << "    Looking at Dig " << Dig << " of " << rawTPC.size() << ", it has " << ThisDig.Samples() << " samples on Channel " << ThisDig.Channel() << " ("<<Chan<<")" 
		<< ", and " << ThisDig.NADC() << " ADCs with a pedestal of " << ThisDig.GetPedestal()
		<< std::endl;
      for (unsigned int ADCLoop=0; ADCLoop < ThisDig.NADC(); ++ADCLoop) {
	if (ADCLoop < 5) {
	  std::cout << "      Looking at ADC " << ADCLoop << " of " << ThisDig.NADC() << ", it was " << ThisDig.ADC(ADCLoop) << std::endl;
	} // If ADCLoop < 5
      } // Loop over ADCs
    } // If Dig < 5

    // --- It can be useful to select only specific types of channels, normally you just continue for one of them.
    if (geom->SignalType(ThisDig.Channel()) == geo::kCollection) {
      // --- Something for collection plane wires...
    } else if (geom->SignalType(ThisDig.Channel()) == geo::kInduction) {
      // --- Something for induction plane wires...
    }

  } // Loop over rawPD

  std::cout << "  I have just got into my activity TPC trigger , fMakeTrig is " << fMakeTrig << ", and (size) fNumber is now " << fNumber << std::endl;
  
  // --- For now just trigger on how raw::RawDigits there are in this event...
  if (fNumber > 3700) {
    fTrigDecision = true;
  } else {
    fTrigDecision = false;
  }

  // --- Return the result of the trigger.
  return fTrigDecision;
} // TriggerOnTPC

//----------------------------------------------------------------------
//-- The trigger algorithm on just the Photon Detector OpDetWaveforms --
//----------------------------------------------------------------------
bool triggersim::ActivityTrigger::TriggerOnPD( std::vector< raw::OpDetWaveform > rawPD) {
  // --- If for some reason you changed your mind about making this trigger...
  if (!fMakeTrig) return false;

  // --- Now do stuff...
  fNumber = rawPD.size();
  for (unsigned int Wave=0; Wave < rawPD.size(); ++Wave) {
    if (Wave < 5) {
      raw::OpDetWaveform ThisWaveform = rawPD[Wave]; // Also, rawPD.at(Wave);
      std::cout << "    Looking at Wave " << Wave << " of " << rawPD.size() << ", it was on channel " << ThisWaveform.ChannelNumber() << ", at time " << ThisWaveform.TimeStamp()
		<< ", there are " << ThisWaveform.Waveform().size() << " ADCs in this waveform " << std::endl;
      std::vector< short > WaveformVec = ThisWaveform.Waveform();
      for (unsigned int WaveformLoop=0; WaveformLoop < WaveformVec.size(); ++WaveformLoop) {
	if (WaveformLoop < 5) {
	  std::cout << "      Element " << WaveformLoop << " of " << WaveformVec.size() << " has ADC value " << WaveformVec.at(WaveformLoop) << std::endl;
	} // If WaveformLoop < 5
      } // Loop over the Waveform() data member
    } // If Wave < 5
  } // Loop over rawPD
  
    // --- If there are any OpDetWaveforms return true.
  if (fNumber) {
    fTrigDecision = true;
  } else {
    fTrigDecision = false;
  }

  // --- Return the result of the trigger.
  return fTrigDecision;
} // TriggerOnPD

//----------------------------------------------------------------------
//----- The trigger algorithm on the RawDigits and OpDetWaveforms ------
//----------------------------------------------------------------------
bool triggersim::ActivityTrigger::TriggerOnTPC_PD( std::vector< raw::RawDigit > rawTPC, std::vector< raw::OpDetWaveform > rawPD) {
  // --- If for some reason you changed your mind about making this trigger...
  if (!fMakeTrig) return false;

  // --- Now do stuff...
  std::cout << "    The size of the RawDigits is " << rawTPC.size() << ", and the size of the OpDetWaveforms is " << rawPD.size() << std::endl;

  // --- Look at the above two function about how to access information out of either of these data products...

  // --- Just default the trigger to true.
  fTrigDecision = true;

  // --- Return the result of the trigger.
  return fTrigDecision;
} // TriggerOnTPC_PD

//----------------------------------------------------------------------
//----- The trigger algorithm on the RawDigits and OpDetWaveforms ------
//----------------------------------------------------------------------
bool triggersim::ActivityTrigger::TriggerOnTriggers( std::vector<triggersim::BasicTrigger> triggerVec) {

  std::cout << "    Looking at TriggerOnTriggers...I currently have " << triggerVec.size() << " triggers." << std::endl;
  for (unsigned int TrigVecLoop=0; TrigVecLoop < triggerVec.size(); ++TrigVecLoop) {
    //triggersim::BasicTrigger ThisTrig = triggerVec.at(TrigVecLoop);
    //std::cout << "      Looking at Trigger " << TrigVecLoop << ", TrigDecision " << ThisTrig.TrigDecision() << ", TrigType " << ThisTrig.TrigType()
    //	      << ", TrigSubtype " << ThisTrig.TrigSubType()
    std::cout << "      Looking at Trigger " << TrigVecLoop << ", TrigDecision " << triggerVec[TrigVecLoop].TrigDecision() << ", TrigType " << triggerVec[TrigVecLoop].TrigType()
	      << ", TrigSubtype " << triggerVec[TrigVecLoop].TrigSubType()

	      << std::endl;
  }


  // --- Just default the trigger to true.
  fTrigDecision = true;

  // --- Return the result of the trigger.
  return fTrigDecision;
} // TriggerOnTriggers


//----------------------------------------------------------------------
#endif
//----------------------------------------------------------------------
