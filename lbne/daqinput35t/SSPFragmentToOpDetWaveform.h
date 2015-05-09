#ifndef SSPFragmentToOpDetWaveform_h
#define SSPFragmentToOpDetWaveform_h

//artdaq includes
#include "artdaq-core/Data/Fragments.hh"


//larsoft includes
#include "RawData/raw.h"
#include "RawData/OpDetWaveform.h"

//lbnecode/daqinput35t includes
#include "utilities/UnpackFragment.h"
#include <vector>

namespace DAQToOffline {

  // Unpack the given artdaq::Fragment objects, and create a vector of raw::OpDetWaveform objects. The
  // Fragments are expected to be carrying Optical Detector data; this is not
  // checked.

  std::vector<raw::OpDetWaveform>
    SSPFragmentToOpDetWaveform(artdaq::Fragments const& raw, const double NOvAClockFrequency);

}
#endif
