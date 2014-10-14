#ifndef lbne_artdaq_Overlays_MicroSliceWriter_hh
#define lbne_artdaq_Overlays_MicroSliceWriter_hh

#include "lbne-raw-data/Overlays/MicroSlice.hh"
#include "lbne-raw-data/Overlays/NanoSliceWriter.hh"

namespace lbne {
  class MicroSliceWriter;
}

class lbne::MicroSliceWriter : public lbne::MicroSlice {

public:

  // This constructor creates an empty MicroSlice which can be filled
  // with the appropriate data
  MicroSliceWriter(uint8_t* address, uint32_t max_size_bytes);

  // Reserves the next NanoSlice in memory within this MicroSlice.
  // The NanoSliceWriter that is returned is initialized to an empty
  // state and is ready to be populated with data.  If an additional
  // NanoSlice cannot be reserved (for example if it would overflow
  // the MicroSlice buffer), an empty shared_ptr is returned.
  std::shared_ptr<NanoSliceWriter> reserveNanoSlice(uint32_t ns_max_bytes);

  // Finalizes the MicroSlice.  This takes care of updating the
  // MicroSlice so that all internal pointers and data are consistent
  // with the NanoSlices that have been added.  No more NanoSlices
  // can be added once the MicroSlice has been finalized.
  // This method returns the number of bytes that were reclaimed
  // when the maximum size for this MicroSlice was reduced to match
  // its actual current size.
  int32_t finalize();

protected:

  // finalizes the NanoSlice that was most recently added and
  // updates our size to reflect any memory that was reclaimed
  // by finalizing the NanoSlice
  void finalizeLatestNanoSlice_();

  // returns a pointer to the header
  Header* header_();

  uint32_t max_size_bytes_;
  std::shared_ptr<NanoSliceWriter> latest_nanoslice_ptr_;
};

#endif /* lbne_artdaq_Overlays_MicroSliceWriter_hh */
