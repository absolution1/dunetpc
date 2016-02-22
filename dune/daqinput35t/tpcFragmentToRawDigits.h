#ifndef tpcFragmentToRawDigits_h
#define tpcFragmentToRawDigits_h

#include "artdaq-core/Data/Fragments.hh"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "art/Framework/Principal/Event.h"

// From lardata
#include "lardata/RawData/RawDigit.h"

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

  struct TPCChannel{
    int fOnlineChannel;
    int fRCE;
    int fRCEChannel;
    int fAPA;
    int fPlane;
    int fOfflineChannel;
    TPCChannel(int online_channel, int rce, int rce_channel, int apa, int plane, int offline_channel):
      fOnlineChannel(online_channel),
      fRCE(rce),
      fRCEChannel(rce_channel),
      fAPA(apa),
      fPlane(plane),
      fOfflineChannel(offline_channel)
    {
      
    }
  };

  struct TPCChannelMapDetailed{
    std::map<int,TPCChannel> fOnlineChannelMap;
    std::map<int,TPCChannel> fOfflineChannelMap;
    TPCChannelMapDetailed(std::string channelMapFile);
  };

  const lbne::TpcNanoSlice::Header::nova_timestamp_t nova_time_ticks_per_second = 64e6;

  art::Timestamp make_art_timestamp_from_nova_timestamp(lbne::TpcNanoSlice::Header::nova_timestamp_t this_nova_timestamp);

}

#endif
