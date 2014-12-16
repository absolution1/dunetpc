#ifndef lbne_artdaq_Overlays_TpcMillisliceWriter_hh
#define lbne_artdaq_Overlays_TpcMilliSliceWriter_hh

#include "lbne-raw-data/Overlays/TpcMilliSlice.hh"
//#include "lbne-raw-data/Overlays/TpcMicroSliceWriter.hh" TODO

namespace lbne {
  class TpcMilliSliceWriter;
}

class lbne::TpcMilliSliceWriter : public lbne::TpcMilliSlice {

public:

  // This constructor creates an empty TpcMilliSlice which can be filled
  // with the appropriate data
  TpcMilliSliceWriter(uint8_t* address, uint32_t max_size_bytes);

  // Reserves the next MicroSlice in memory within this TpcMilliSlice.
  // The MicroSliceWriter that is returned is initialized to an empty
  // state and is ready to be populated with data.  This method returns
  // true if the MicroSlice was successfully added, false if not (for example,
  // if the additional MicroSlice would overflow the maximum size of the
  // TpcMilliSlice).
  // std::shared_ptr<MicroSliceWriter> reserveMicroSlice(uint32_t ms_max_bytes); TODO ??

  // Finalizes the TpcMilliSlice.  This takes care of updating the
  // TpcMilliSlice so that all internal pointers and data are consistent
  // with the MicroSlices that have been added.  No more MicroSlices
  // can be added once the TpcMilliSlice has been finalized.
  // This method returns the number of bytes that were reclaimed
  // when the maximum size for this TpcMilliSlice was reduced to match
  // its actual current size.
  int32_t finalize(uint32_t data_size_bytes, uint32_t microslice_count);

protected:

  // finalizes the MicroSlice that was most recently added and
  // updates our size to reflect any memory that was reclaimed
  // by finalizing the MicroSlice
  // void finalizeLatestMicroSlice_(); TODO

  // returns a pointer to the header
  Header* header_();

  // returns a pointer to the requested MicroSlice
  // uint8_t* data_(int index); TODO

  uint32_t max_size_bytes_;
  // std::shared_ptr<MicroSliceWriter> latest_microslice_ptr_; TODO
};

#endif /* lbne_artdaq_Overlays_TpcMilliSliceWriter_hh */
