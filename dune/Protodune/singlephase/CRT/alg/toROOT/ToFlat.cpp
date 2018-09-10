//File: ToFlat.cpp
//Brief: Converts Triggers to a flat TTree.  Useful when Sarah wants to access CRT output 
//       in python.  
//Usage: Provide one or more EBuilder output files on
// the command line.  Files should have been produced by EBuilder.
//
//       fromEBuilder/toROOTFlat <one EBuilder file name> [another EBuilder file name]...
//
//Author: Andrew Olivier aolivier@ur.rochester.edu

//local includes
#include "dunetpc/dune/Protodune/singlephase/CRT/alg/toROOT/ToFlat.h"

namespace CRT
{
  ToFlat::~ToFlat() {}

  void ToFlat::AnalyzeEvent(const std::vector<CRT::Trigger>& triggers)
  {
    ResetBranches();

    size_t trigIt = 0;
    fNTriggers = triggers.size();
    for(const auto& trigger: triggers)
    {
      fTimestamps[trigIt] = trigger.Timestamp();
      fModules[trigIt] = trigger.Channel();

      const auto& hits = trigger.Hits();
      for(const auto& hit: hits)
      {
        fADCs[trigIt][hit.Channel()] = hit.ADC();
      }
      ++trigIt;
    }

    fTree->Fill();
  }

  void ToFlat::ResetBranches()
  {
    for(size_t trigger = 0; trigger < MaxTriggers; ++trigger)
    {
      fNTriggers = 0;
      fTimestamps[trigger] = -1;
      fModules[trigger] = -1;
      for(size_t channel = 0; channel < NChannels; ++channel) fADCs[trigger][channel] = 0;
    }
  }
}
