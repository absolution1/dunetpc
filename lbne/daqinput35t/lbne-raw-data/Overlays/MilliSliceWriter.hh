#ifndef lbne_artdaq_Overlays_MilliSliceWriter_hh
#define lbne_artdaq_Overlays_MilliSliceWriter_hh

#include "lbne-raw-data/Overlays/MilliSlice.hh"
#include "lbne-raw-data/Overlays/MicroSliceWriter.hh"

namespace lbne {
  class MilliSliceWriter;
}

class lbne::MilliSliceWriter : public lbne::MilliSlice {

public:

  // This constructor creates an empty MilliSlice which can be filled
  // with the appropriate data
  MilliSliceWriter(uint8_t* address, uint32_t max_size_bytes);

  // Reserves the next MicroSlice in memory within this MilliSlice.
  // The MicroSliceWriter that is returned is initialized to an empty
  // state and is ready to be populated with data.  This method returns
  // true if the MicroSlice was successfully added, false if not (for example,
  // if the additional MicroSlice would overflow the maximum size of the
  // MilliSlice).
  std::shared_ptr<MicroSliceWriter> reserveMicroSlice(uint32_t ms_max_bytes);

  // Finalizes the MilliSlice.  This takes care of updating the
  // MilliSlice so that all internal pointers and data are consistent
  // with the MicroSlices that have been added.  No more MicroSlices
  // can be added once the MilliSlice has been finalized.
  // This method returns the number of bytes that were reclaimed
  // when the maximum size for this MilliSlice was reduced to match
  // its actual current size.
  int32_t finalize();

protected:

  // finalizes the MicroSlice that was most recently added and
  // updates our size to reflect any memory that was reclaimed
  // by finalizing the MicroSlice
  void finalizeLatestMicroSlice_();

  // returns a pointer to the header
  Header* header_();

  // returns a pointer to the requested MicroSlice
  uint8_t* data_(int index);

  uint32_t max_size_bytes_;
  std::shared_ptr<MicroSliceWriter> latest_microslice_ptr_;
};

#endif /* lbne_artdaq_Overlays_MilliSliceWriter_hh */
