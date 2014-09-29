#ifndef lbne_artdaq_Overlays_MicroSlice_hh
#define lbne_artdaq_Overlays_MicroSlice_hh

#include "lbne-raw-data/Overlays/NanoSlice.hh"
#include <memory>

namespace lbne {
  class MicroSlice;
}

class lbne::MicroSlice {

public:

  struct Header {
    typedef uint32_t data_t;

    typedef uint8_t version_t;
    typedef uint32_t microslice_size_t;  
    typedef uint32_t nanoslice_count_t;  

    // this structure is deliberately incomplete - it should be
    // improved to match what we want for the LBNE 35t detector

    data_t version         :  8;
    data_t microslice_size : 24;   // total size, data & header

    data_t nanoslice_count : 32;
  };

  // This constructor accepts a memory buffer that contains an existing
  // MicroSlice and allows the the data inside it to be accessed
  MicroSlice(uint8_t* address);

  // Returns the size of the MicroSlice
  Header::microslice_size_t size() const;

  // Returns the number of NanoSlices in this MicroSlice
  Header::nanoslice_count_t nanoSliceCount() const;

  // Returns the requested NanoSlice if the requested slice was found,
  // otherwise returns an empty pointer
  std::unique_ptr<NanoSlice> nanoSlice(uint32_t index) const;

protected:

  // returns a pointer to the header
  Header const* header_() const;

  // returns a pointer to the requested NanoSlice
  uint8_t* data_(uint32_t index) const;

  uint8_t* buffer_;
};

#endif /* lbne_artdaq_Overlays_MicroSlice_hh */
