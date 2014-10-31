#ifndef lbne_artdaq_Overlays_MilliSliceFragment_hh
#define lbne_artdaq_Overlays_MilliSliceFragment_hh

#include "lbne-raw-data/Overlays/MilliSlice.hh"
#include "lbne-raw-data/Overlays/MicroSlice.hh"
#include "artdaq-core/Data/Fragment.hh"

namespace lbne {
  class MilliSliceFragment;
}

class lbne::MilliSliceFragment : public lbne::MilliSlice {

public:

  // This constructor accepts an artdaq::Fragment that contains an existing
  // MilliSlice and allows the the data inside the MilliSlice to be accessed
  MilliSliceFragment(artdaq::Fragment const& frag);

protected:

  // returns a pointer to the header
  Header const* header_() const;

  // returns a pointer to the requested MicroSlice
  uint8_t* data_(int index) const;

private:

  artdaq::Fragment const& artdaq_fragment_;
};

#endif /* lbne_artdaq_Overlays_MilliSliceFragment_hh */
