#ifndef lbne_artdaq_Overlays_NanoSliceWriter_hh
#define lbne_artdaq_Overlays_NanoSliceWriter_hh

#include "lbne-raw-data/Overlays/NanoSlice.hh"

namespace lbne {
  class NanoSliceWriter;
}

class lbne::NanoSliceWriter : public lbne::NanoSlice {

public:

  // This constructor creates an empty NanoSlice which can be filled
  // with the appropriate data
  NanoSliceWriter(uint8_t* address, uint32_t max_size_bytes,
                  uint16_t channel_number = 0);

  // Sets the channel number for this NanoSlice.
  void setChannelNumber(uint16_t channel);

  // Adds the specified data sample NanoSlice.  Returns true
  // if it was successfully added, false if not (for example, if the
  // additional sample would overflow the maximum size of the
  // NanoSlice).
  bool addSample(uint16_t value);

  // Finalizes the NanoSlice.  This takes care of updating the
  // NanoSlice so that all internal pointers and data are consistent
  // with the data that has been added.  No more samples can be added
  // once the NanoSlice has been finalized.
  // This method returns the number of bytes that were reclaimed
  // when the maximum size for this NanoSlice was reduced to match
  // its actual current size.
  int32_t finalize();

protected:

  // returns a pointer to the header
  Header * header_();

  // returns a pointer to the first sample value
  uint16_t * data_();

  uint32_t max_size_bytes_;
};

#endif /* lbne_artdaq_Overlays_NanoSliceWriter_hh */
