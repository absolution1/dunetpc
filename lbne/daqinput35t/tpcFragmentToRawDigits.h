#ifndef tpcFragmentToRawDigits_h
#define tpcFragmentToRawDigits_h

#include "art/Persistency/Provenance/EventID.h"
#include "artdaq-core/Data/Fragment.hh"
#include "artdaq-core/Data/Fragments.hh"

// From lardata
#include "RawData/RawDigit.h"

#include <vector>

namespace DAQToOffline {

  // Unpack the given Fragments, and create a vector of RawDigits. The
  // Fragments are expected to be carrying TPC data; this is not
  // checked.

  std::vector<raw::RawDigit>
  tpcFragmentToRawDigits(art::EventID const& id,
			 artdaq::Fragments const& rawFragments,
			 bool debug,
			 raw::Compress_t compression,
			 unsigned int zeroThreshold);
}

#endif
