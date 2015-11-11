// -*- mode: c++; c-basic-offset: 2; -*-
#ifndef SSPFragmentToOpDetWaveform_h
#define SSPFragmentToOpDetWaveform_h

//artdaq includes
#include "artdaq-core/Data/Fragments.hh"


//larsoft includes
#include "RawData/raw.h"
#include "RawData/OpDetWaveform.h"

// lbne-raw-data includes
#include "lbne-raw-data/Overlays/SSPFragment.hh"
#include "lbne-raw-data/Overlays/anlTypes.hh"

//lbnecode/daqinput35t includes
#include "utilities/UnpackFragment.h"
#include <vector>
#include <map>

namespace DAQToOffline {

  // Unpack the given artdaq::Fragment objects, and create a vector of raw::OpDetWaveform objects. The
  // Fragments are expected to be carrying Optical Detector data; this is not
  // checked.

  std::vector<raw::OpDetWaveform>
    SSPFragmentToOpDetWaveform(artdaq::Fragments const& raw, const double NOvAClockFrequency, const std::map<int,int> theChannelMap);


  void BuildOpDetChannelMap(std::string fChannelMapFile, std::map<int,int> &theChannelMap);

  // Load the milislice
  //const SSPDAQ::MillisliceHeader* CheckMilisliceMetadata(const auto& frag, lbne::SSPFragment sspf);
  
  // Extract data from the header
  uint32_t GetPeaksum(const SSPDAQ::EventHeader* daqHeader);
  unsigned short GetOpChannel(const SSPDAQ::EventHeader* daqHeader, std::map<int,int> theChannelMap);
  unsigned long GetGlobalFirstSample(const SSPDAQ::EventHeader* daqHeader);
  unsigned long GetInternalFirstSample(const SSPDAQ::EventHeader *daqHeader);
  void PrintHeaderInfo(const SSPDAQ::EventHeader *daqHeader, const double NOvAClockFrequency = 64);

}
#endif
