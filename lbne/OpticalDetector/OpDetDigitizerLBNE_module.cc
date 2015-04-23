//=========================================================
// OpDetDigitizerLBNE_module.cc
// This module produces digitized output 
// (creating OpDetWaveform)
// from photon detectors taking SimPhotonsLite as input.
//
// Gleb Sinev, Duke, 2015
// Based on OpMCDigi_module.cc
//=========================================================

#ifndef OpDetDigitizerLBNE_h
#define OpDetDigitizerLBNE_h 1

// Framework includes

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
//#include "art/Framework/Services/Optional/RandomNumberGenerator.h"
#include "fhiclcpp/ParameterSet.h"

// ART extensions
//#include "artextensions/SeedService/SeedService.hh"

// LArSoft includes

#include "Simulation/sim.h"
#include "Simulation/SimPhotons.h"
#include "Simulation/LArG4Parameters.h"
//#include "Geometry/Geometry.h"
#include "Utilities/TimeService.h"
#include "OpticalDetector/OpDetResponseInterface.h"
#include "RawData/OpDetWaveform.h"
//#include "RawData/OpDetPulse.h"
//#include "OpticalDetectorData/OpticalRawDigit.h"

// CLHEP includes

//#include "CLHEP/Random/RandBinomial.h"

// C++ includes

#include <vector>
#include <map>
//#include <iostream>
#include <cmath>

namespace opdet {

  class OpDetDigitizerLBNE : public art::EDProducer{

    public:
      
      OpDetDigitizerLBNE(fhicl::ParameterSet const&);
      virtual ~OpDetDigitizerLBNE();
      
      void produce(art::Event&);

//      void beginJob();

    private:

      // The parameters read from the .fcl file
      std::string fInputModule;            // Input tag for OpDet collection
      //double fQE;                          // Quantum efficiency
      double fSampleFreq;                  // Frequency in GHz (number of ticks in one ns)
      double fTimeBegin;                   // Beginning of sample in ns
      double fTimeEnd;                     // End of sample in ns
      double fVoltageToADC;                // Conversion factor mV to ADC counts

      // Random number engines
      //CLHEP::RandBinomial *fRandBinomial;

      // Function that adds n pulses to a waveform
      void AddPulse(int timeBin, int scale, std::vector< double >& waveform);

      // Functional response to one photoelectron (time in ns)
      double Pulse1PE(double time);

      // One photoelectron pulse parameters
      double fPulseLength;  // 1PE pulse length in ticks
      double fPeakTime;     // Time when the pulse reaches its maximum in ns
      double fMaxAmplitude; // Maximum amplitude of the pulse in mV
      double fFrontTime;    // Constant in the exponential function in ns
      double fBackTime;     // Constant in the exponential function in ns

      std::vector< double > fSinglePEWaveform;
      void CreateSinglePEWaveform();

  };

}

#endif

namespace opdet {

  DEFINE_ART_MODULE(OpDetDigitizerLBNE)

}

namespace opdet {
  
  //---------------------------------------------------------------------
  // Constructor
  OpDetDigitizerLBNE::OpDetDigitizerLBNE(fhicl::ParameterSet const& pset)
  {

    // This module produces (infrastructure piece)
    //produces< std::vector< raw::OpDetPulse > >();
    //produces< std::vector< optdata::OpticalRawDigit > >();
    produces< std::vector< raw::OpDetWaveform > >();

    // Read the fcl-file
    fInputModule  = pset.get< std:: string >("InputModule" );
    //fQE           = pset.get< double       >("QE"          );
    fVoltageToADC = pset.get< double       >("VoltageToADC");
    //fSampleFreq   = pset.get< double >("SampleFreq");
    //fTimeBegin    = pset.get< double >("TimeBegin");
    //fTimeEnd      = pset.get< double >("TimeEnd");

    art::ServiceHandle< util::TimeService > timeService;
    // Converting MHz into GHz and us into ns
    fSampleFreq = timeService->OpticalClock().Frequency()  /1000.0;
    fTimeBegin  = timeService->OpticalClock().Time()       *1000.0;
    fTimeEnd    = timeService->OpticalClock().FramePeriod()*1000.0;
/*
    // Initializing random number engines
    unsigned int seed = pset.get< unsigned int >("Seed", sim::GetRandomNumberSeed());
    createEngine(seed);

    art::ServiceHandle< art::RandomNumberGenerator > rng;
    CLHEP::HepRandomEngine &engine = rng->getEngine();
    fRandBinomial = new CLHEP::RandBinomial(engine);
*/
    // Creating single photoelectron waveform
    fPulseLength = 4000.0;
    fPeakTime = 260.0;
    fMaxAmplitude = 0.12;
    fFrontTime = 9.0;
    fBackTime  = 476.0;
    CreateSinglePEWaveform();

  }

  //---------------------------------------------------------------------
  // Destructor
  OpDetDigitizerLBNE::~OpDetDigitizerLBNE()
  {
  }
/*
  //---------------------------------------------------------------------
  void OpDetDigitizerLBNE::beginJob()
  {
  }
 */ 
  //---------------------------------------------------------------------
  void OpDetDigitizerLBNE::produce(art::Event& evt)
  {
    
    // A pointer that will store produced OpDetWaveforms
    //std::unique_ptr< std::vector< raw::OpDetPulse > > pulseVecPtr(new std::vector< raw::OpDetPulse >);
    //std::unique_ptr< std::vector< optdata::OpticalRawDigit > > 
    //                                         pulseVecPtr(new std::vector< optdata::OpticalRawDigit >);
    std::unique_ptr< std::vector< raw::OpDetWaveform > > 
                                             pulseVecPtr(new std::vector< raw::OpDetWaveform >);
    
    art::ServiceHandle< sim::LArG4Parameters > lgp;
    bool fUseLitePhotons = lgp->UseLitePhotons();

    // Total number of ticks in our readout
    int nSamples = (fTimeEnd - fTimeBegin)*fSampleFreq;

    //art::ServiceHandle< geo::Geometry > geom;
    //int nOpChannels = geom->NOpChannels();

    // Service for determining opdet responses
    art::ServiceHandle< opdet::OpDetResponseInterface > odResponse;
    int nOpChannels = odResponse->NOpChannels();

    // This vector stores waveforms we create for each optical detector
    std::vector< std::vector< double > > opDetWaveforms(nOpChannels, std::vector< double >(nSamples,0.0));

    if (fUseLitePhotons)
    {
      // Get SimPhotonsLite from the event
      art::Handle< std::vector< sim::SimPhotonsLite > > litePhotonHandle;
      evt.getByLabel(fInputModule, litePhotonHandle);

      // For every optical detector:
      for (auto const& opDet : (*litePhotonHandle))
      {
        int channel = opDet.OpChannel;
        int readoutCh;
        std::map< int, int > photonsMap = opDet.DetectedPhotons;

        // For every pair of (arrival time, number of photons) in the map:
        for (auto const& pulse : photonsMap)
        {
          for (int i = 0; i < pulse.second; i++)
          {
            // Sample a random subset according to QE
            if (odResponse->detectedLite(channel, readoutCh) && 
                (pulse.first > fTimeBegin) && (pulse.first < fTimeEnd))
            {
              // Convert the time of the pulse to ticks
              int timeBin = int((pulse.first - fTimeBegin)*fSampleFreq);
              /*
              // Take into account the quantum efficiency
              int detPhotons = pulse.second;
              if (fQE < 1.0) detPhotons = fRandBinomial->fire(pulse.second,fQE);
              */
              // Add the pulse to the waveform
              //AddPulse(timeBin, detPhotons, opDetWaveforms[channel]);
              AddPulse(timeBin, 1, opDetWaveforms[readoutCh]);
            }
          }
        }
      }
    }
    else std::cout << "Sorry, but for now only Lite Photons are implemented!\n";

    for (int channel = 0; channel != nOpChannels; channel++)
    {
      // Produce ADC pulse of integers rather than doubles
      //std::vector< short > adcVec;
      //optdata::OpticalRawDigit adcVec(0, channel, nSamples);
      raw::OpDetWaveform adcVec(0, channel, nSamples);

      for (auto& value : opDetWaveforms[channel])
      {
        // 0.5 here is to round the value correctly
        adcVec.push_back(int(value + 0.5));
      }

      //pulseVecPtr->push_back(raw::OpDetPulse(channel, adcVec, 0, fTimeBegin));
      pulseVecPtr->push_back(adcVec);

    }

    evt.put(std::move(pulseVecPtr));

  }

  //---------------------------------------------------------------------
  void OpDetDigitizerLBNE::AddPulse(int timeBin, int scale, std::vector< double >& waveform)
  {

    // How many bins will be changed
    size_t pulseLength = fSinglePEWaveform.size();
    if ((timeBin + fSinglePEWaveform.size()) > waveform.size()) 
      pulseLength = (waveform.size() - timeBin);

    // Adding a pulse to the waveform
    for (size_t tick = 0; tick != pulseLength; tick++)
    {
      waveform.at(timeBin + tick) += scale*fSinglePEWaveform.at(tick);
    }

  }

  //---------------------------------------------------------------------
  double OpDetDigitizerLBNE::Pulse1PE(double time)
  {

    if (time < fPeakTime) return (fVoltageToADC*fMaxAmplitude*std::exp((time - fPeakTime)/fFrontTime));
    else return (fVoltageToADC*fMaxAmplitude*std::exp(-(time - fPeakTime)/fBackTime));

  }

  void OpDetDigitizerLBNE::CreateSinglePEWaveform()
  {

    int length = fPulseLength*fSampleFreq;
    fSinglePEWaveform.resize(length);
    for (size_t tick = 0; tick != fSinglePEWaveform.size(); tick++)
    {
      fSinglePEWaveform.at(tick) = Pulse1PE(double(tick)/fSampleFreq);
    }

  }

}
