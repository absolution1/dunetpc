//==========================================================
// OpDetDigiAnaLBNE_module.cc
// This analysis module creates histograms
// with information from OpDetWaveforms
//
// Gleb Sinev, Duke, 2015
// Based on OpDigiAna_module.cc
//==========================================================

#ifndef OpDetDigiAnaLBNE_H
#define OpDetDigiAnaLBNE_H 1

// Framework includes

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Optional/TFileDirectory.h"
#include "art/Persistency/Common/Ptr.h"
#include "fhiclcpp/ParameterSet.h"

// LArSoft includes

#include "Utilities/TimeService.h"
#include "RawData/OpDetWaveform.h"

// ROOT includes

#include "TH1.h"

// C++ includes

#include <vector>
#include <cstring>

namespace opdet {

  class OpDetDigiAnaLBNE : public art::EDAnalyzer {

    public:

      // Standard constructor and destructor for an ART module
      OpDetDigiAnaLBNE(fhicl::ParameterSet const&);
      virtual ~OpDetDigiAnaLBNE();

      // The analyzer routine, called one per event
      void analyze(art::Event const&);
   
    private:

      // Parameters we'll read from the fcl-file
      std::string fInputModule; // Module used to create OpDetWaveforms
      std::string fInstanceName;// Input tag for OpDetWaveforms collection
      float fSampleFreq;        // Sampling frequency in MHz 
      float fTimeBegin;         // Beginning of sample in us
      float fTimeEnd;           // End of sample in us

  };

}

#endif 

namespace opdet {

  DEFINE_ART_MODULE(OpDetDigiAnaLBNE)

}

namespace opdet {

  //---------------------------------------------------------------------------
  // Constructor
  OpDetDigiAnaLBNE::OpDetDigiAnaLBNE(fhicl::ParameterSet const& pset)
    : EDAnalyzer(pset)
  {

    // Read the fcl-file
    fInputModule = pset.get< std::string >("InputModule");
    fInstanceName = pset.get<std::string>("InstanceName");

    // Obtaining parameters from TimeService
    art::ServiceHandle< util::TimeService > timeService;
    fSampleFreq = timeService->OpticalClock().Frequency();
    fTimeBegin  = 0.0; //timeService->OpticalClock().Time();
    fTimeEnd    = 1600.0; //timeService->OpticalClock().FramePeriod();

  }

  //---------------------------------------------------------------------------
  // Destructor
  OpDetDigiAnaLBNE::~OpDetDigiAnaLBNE()
  {
  }

  //---------------------------------------------------------------------------
  void OpDetDigiAnaLBNE::analyze(art::Event const& evt)
  {

    // Create a string for histogram names
    char histName[50];

    // Get OpDetWaveforms from the event
    art::Handle< std::vector< raw::OpDetWaveform > > waveformHandle;
    evt.getByLabel(fInputModule, fInstanceName, waveformHandle);

    // Access ART's TFileService, which will handle creating and writing
    // histograms for us
    art::ServiceHandle< art::TFileService > tfs;

    for (size_t i = 0; i < waveformHandle->size(); i++)
    {

      // This was probably required to overcome the "const" problem 
      // with OpDetPulse::Waveform()
      art::Ptr< raw::OpDetWaveform > waveformPtr(waveformHandle, i);
      raw::OpDetWaveform pulse = *waveformPtr;
      // Make a name for the histogram
      sprintf(histName, "event_%d_opdet_%i", evt.id().event(), 
                                        pulse.ChannelNumber());

      TH1D * waveformHist = nullptr;

      waveformHist = tfs->make< TH1D >(histName, ";t (us);",
                                       int((fTimeEnd - fTimeBegin)
                                                 *fSampleFreq - 1),
                                              fTimeBegin, fTimeEnd);

      for (size_t tick = 0; tick < pulse.size(); tick++)
        waveformHist->SetBinContent(tick, (float) pulse[tick]);

    }

  }

}
