#ifndef lbne_artdaq_Overlays_MilliSliceFragmentWriter_hh
#define lbne_artdaq_Overlays_MilliSliceFragmentWriter_hh

#include "lbne-raw-data/Overlays/MilliSliceFragment.hh"
#include "lbne-raw-data/Overlays/MicroSliceWriter.hh"

namespace lbne {
  class MilliSliceFragmentWriter;
}

class lbne::MilliSliceFragmentWriter : public lbne::MilliSliceFragment {

public:

  // This constructor accepts an empty artdaq::Fragment and creates
  // a MilliSlice inside it that can be filled with the appropriate data
  MilliSliceFragmentWriter(artdaq::Fragment& frag, uint32_t max_size_bytes);

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

private:

  // Note that this non-const reference hides the const reference in the base
  // class.
  // Also Note that we store a reference to the artdaq::Fragment and look up
  // the dataBeginBytes() each time that we need it so that we are insulated
  // from changes in the overall layout of the Fragment.  Of course, what we
  // have in the MilliSlice, MicroSlice, and NanoSlice writers isn't infinitely
  // flexible:  If the Fragment is modified in the middle of populating data at
  // one of these levels, then users should understand that their existing
  // MicroSliceWriter and NanoSliceWriter references will be invalid.
  artdaq::Fragment& artdaq_fragment_;
};

#endif /* lbne_artdaq_Overlays_MilliSliceFragmentWriter_hh */
