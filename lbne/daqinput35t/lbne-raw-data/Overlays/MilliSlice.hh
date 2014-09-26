#ifndef lbne_artdaq_Overlays_MilliSlice_hh
#define lbne_artdaq_Overlays_MilliSlice_hh

#include "lbne-raw-data/Overlays/MicroSlice.hh"
#include "artdaq-core/Data/Fragment.hh"

namespace lbne {
  class MilliSlice;
}

class lbne::MilliSlice {

public:

  struct Header {
    typedef uint32_t data_t;

    typedef uint16_t pattern_t;
    typedef uint16_t version_t;
    typedef uint32_t millislice_size_t;  
    typedef uint32_t microslice_count_t;  

    // this structure is deliberately incomplete - it should be
    // improved to match what we want for the LBNE 35t detector

    data_t fixed_pattern : 16;
    data_t version       : 16;

    data_t millislice_size : 32;   // total size, data & header

    data_t microslice_count : 32;
  };

  // This constructor accepts a memory buffer that contains an existing
  // MilliSlice and allows the the data inside it to be accessed
  MilliSlice(uint8_t* address);

  // Returns the size of the MilliSlice
  Header::millislice_size_t size() const;

  // Returns the number of MicroSlices in this MilliSlice
  Header::microslice_count_t microSliceCount() const;

  // Returns the requested MicroSlice if the requested slice was found,
  // otherwise returns an empty pointer
  std::unique_ptr<MicroSlice> microSlice(uint32_t index) const;

protected:

  // returns a pointer to the header
  Header const* header_() const;

  // returns a pointer to the requested MicroSlice
  uint8_t* data_(int index) const;

  uint8_t* buffer_;
};

#endif /* lbne_artdaq_Overlays_MilliSlice_hh */
