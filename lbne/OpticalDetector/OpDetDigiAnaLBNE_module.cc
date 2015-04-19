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

//#include "RawData/OpDetPulse.h"
#include "Utilities/TimeService.h"
//#include "OpticalDetectorData/OpticalRawDigit.h"
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
      double fSampleFreq;       // Frequency in GHz (number of ticks in one ns)
      double fTimeBegin;        // Beginning of sample in ns
      double fTimeEnd;          // End of sample in ns

  };

}

#endif 

namespace opdet {

  DEFINE_ART_MODULE(OpDetDigiAnaLBNE)

}

namespace opdet {

  //----------------------------------------------------------------------------
  // Constructor
  OpDetDigiAnaLBNE::OpDetDigiAnaLBNE(fhicl::ParameterSet const& pset)
    : EDAnalyzer(pset)
  {

    // Read the fcl-file
    fInputModule = pset.get< std::string >("InputModule");
//    fSampleFreq  = pset.get< double >("SampleFreq");
//    fTimeBegin   = pset.get< double >("TimeBegin");
//    fTimeEnd     = pset.get< double >("TimeEnd");

    art::ServiceHandle< util::TimeService > timeService;
    // Converting MHz into GHz and us into ns
    fSampleFreq = timeService->OpticalClock().Frequency()/1000.0;
    fTimeBegin  = 0; //timeService->OpticalClock().Time()*1000.0;
    fTimeEnd    = 8000; //timeService->OpticalClock().FramePeriod()*1000.0;

  }

  //----------------------------------------------------------------------------
  // Destructor
  OpDetDigiAnaLBNE::~OpDetDigiAnaLBNE()
  {
  }

  //----------------------------------------------------------------------------
  void OpDetDigiAnaLBNE::analyze(art::Event const& evt)
  {

    // Create a string for histogram names
    char histName[50];

    // Get OpDetWaveforms from the event
    //art::Handle< std::vector< raw::OpDetPulse > > waveformHandle;
    //art::Handle< std::vector< optdata::OpticalRawDigit > > waveformHandle;
    art::Handle< std::vector< raw::OpDetWaveform > > waveformHandle;
    evt.getByLabel(fInputModule, waveformHandle);

    // Access ART's TFileService, which will handle creating and writing
    // histograms for us
    art::ServiceHandle< art::TFileService > tfs;

    for (unsigned int i = 0; i < waveformHandle->size(); i++)
    {
      // This is probably required to overcome the "const" problem 
      // with OpDetPulse::Waveform()
      //art::Ptr< raw::OpDetPulse > waveformPtr(waveformHandle, i);
      //art::Ptr< optdata::OpticalRawDigit > waveformPtr(waveformHandle, i);
      art::Ptr< raw::OpDetWaveform > waveformPtr(waveformHandle, i);
      //raw::OpDetPulse pulse = *waveformPtr;
      //optdata::OpticalRawDigit pulse = *waveformPtr;
      raw::OpDetWaveform pulse = *waveformPtr;
      // Make a name for the histogram
      //sprintf(histName, "event_%d_opdet_%i", evt.id().event(), pulse.OpChannel());
      sprintf(histName, "event_%d_opdet_%i", evt.id().event(), pulse.ChannelNumber());

      int total = 0;
      for (unsigned int tick = 0; tick < pulse.size(); tick++)
        total += pulse[tick];

      if (total < 1) continue;
      
      TH1D * waveformHist = nullptr;

      waveformHist = tfs->make< TH1D >(histName, ";t (ns);",
                                       int((fTimeEnd - fTimeBegin)*fSampleFreq),
                                                            fTimeBegin, fTimeEnd);

      //for (unsigned int tick = 0; tick < pulse.Waveform().size(); tick++)
      for (unsigned int tick = 0; tick < pulse.size(); tick++)
        //waveformHist->SetBinContent(tick, (double) pulse.Waveform()[tick]);
        waveformHist->SetBinContent(tick, (double) pulse.at(tick));

    }

  }

}
