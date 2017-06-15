#ifndef tpcFragmentToRawDigits_h
#define tpcFragmentToRawDigits_h

// framework
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// lbne-raw-data
#include "lbne-raw-data/Services/ChannelMap/ChannelMapService.h"
#include "lbne-raw-data/Overlays/TpcMilliSliceFragment.hh"
#include "artdaq-core/Data/Fragment.hh"

// From dunetpc
#include "utilities/UnpackFragment.h"

// From larcore
#include "larcore/Geometry/Geometry.h"

// From lardata
#include "lardataobj/RawData/raw.h"
#include "lardataobj/RawData/RawDigit.h"
#include "TTimeStamp.h"

#include <map>
#include <iostream>
#include <vector>
#include <fstream>

namespace DAQToOffline {

  // Unpack the given Fragments, and create a vector of RawDigits. The
  // Fragments are expected to be carrying TPC data; this is not
  // checked.

  std::vector<raw::RawDigit> tpcFragmentToRawDigits(artdaq::Fragments const& rawFragments,
						    std::vector<std::pair<std::pair<unsigned int,unsigned int>, lbne::TpcNanoSlice::Header::nova_timestamp_t> > &DigitsIndexList,
						    lbne::TpcNanoSlice::Header::nova_timestamp_t& firstTimestamp,
						    art::ServiceHandle<lbne::ChannelMapService> const& channelMap, bool useChannelMap,
						    bool debug,
						    raw::Compress_t compression,
						    unsigned int zeroThreshold);

  const lbne::TpcNanoSlice::Header::nova_timestamp_t nova_time_ticks_per_second = 64e6;

  // Convert nova time to DUNE timestamp.
  // Please call instead DuneTimeConverter::fromNova(novaTime)
  art::Timestamp make_art_timestamp_from_nova_timestamp(lbne::TpcNanoSlice::Header::nova_timestamp_t this_nova_timestamp);

  // For those who want the old convention where seconds are stored in the low word.
  art::Timestamp old_make_art_timestamp_from_nova_timestamp(lbne::TpcNanoSlice::Header::nova_timestamp_t this_nova_timestamp);

}

#endif
