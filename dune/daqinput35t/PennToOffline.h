////////////////////////////////////////////////////////////////////////
// File:   PennToOffline.cc
// Author: Karl Warburton (Oct 2015)
//
// Utility to provide methods for reformatting the raw online PENN data
// into a structure which can be used in the splitter. 
// Heavily uses lbne-artdaq/OnlineMonitoring/DataReformatter.cxx as a base
////////////////////////////////////////////////////////////////////////

#ifndef PennToOffline_h
#define PennToOffline_h

#include "art/Framework/Principal/Handle.h"
#include "lbne-raw-data/Overlays/PennMilliSliceFragment.hh"
#include "artdaq-core/Data/Fragments.hh"
#include "utilities/UnpackFragment.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "Geometry/Geometry.h"
#include "RawData/ExternalTrigger.h"

#include <iostream>
#include <vector>
#include <map>
#include <bitset>
#include <utility>
#include <iterator>
#include "TMath.h"

namespace DAQToOffline {

    namespace TypeSizes {
    static int const CounterWordSize = 128;
    static int const TriggerWordSize = 32;
  }

  // Unpack the given artdaq::Fragment objects, and create a vector of raw::ExternalTrigger objects. The
  // Fragments are expected to be carrying Penn Trigger board data; this is not checked.
  std::vector<raw::ExternalTrigger> 
    PennFragmentToExternalTrigger( artdaq::Fragments const& Fragments );

  void CollectCounterBits(uint8_t* payload, size_t payload_size, std::vector<std::bitset<TypeSizes::CounterWordSize> > &fCounterBits);
  void CollectMuonTrigger(uint8_t* payload, size_t payload_size, std::map<int,int> &fMuonTriggerRates, std::vector<std::bitset<TypeSizes::TriggerWordSize> > &fMuonTriggerBits,
			  std::vector<int> &fMuonTriggerTimes, lbne::PennMicroSlice::Payload_Header::short_nova_timestamp_t timestamp);
  void ReverseBits(std::bitset<TypeSizes::TriggerWordSize> &bits);
  void ReverseBits(std::bitset<TypeSizes::CounterWordSize> &bits);
}
#endif
