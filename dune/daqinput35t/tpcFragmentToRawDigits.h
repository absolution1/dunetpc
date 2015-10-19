#ifndef tpcFragmentToRawDigits_h
#define tpcFragmentToRawDigits_h

#include "artdaq-core/Data/Fragments.hh"
#include "messagefacility/MessageLogger/MessageLogger.h"


// From lardata
#include "RawData/RawDigit.h"

// from lbne-raw-data

#include "lbne-raw-data/Overlays/TpcMilliSliceFragment.hh"

#include <vector>
#include <map>
#include <fstream>

namespace DAQToOffline {

  // Unpack the given Fragments, and create a vector of RawDigits. The
  // Fragments are expected to be carrying TPC data; this is not
  // checked.

  std::vector<raw::RawDigit> tpcFragmentToRawDigits(artdaq::Fragments const& rawFragments,
						    lbne::TpcNanoSlice::Header::nova_timestamp_t& firstTimestamp,
						    std::map<int,int> const& channelMap,
						    bool debug,
						    raw::Compress_t compression,
						    unsigned int zeroThreshold);

  void BuildTPCChannelMap(std::string channelMapFile, std::map<int,int>& channelMap);

}

#endif
