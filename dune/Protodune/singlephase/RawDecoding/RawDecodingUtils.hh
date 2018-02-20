#ifndef RAW_DECODING_UTILS__XXX
#define RAW_DECODING_UTILS__XXX

#include "artdaq-core/Data/Fragment.hh"

#include "art/Framework/Principal/Event.h"

#include <string>

namespace dune {

  artdaq::Fragments
  getByLabelChecked(const art::Event& evt, const std::string& label, 
		    const std::string& fragtype, const bool& inContainer);
}

#endif
