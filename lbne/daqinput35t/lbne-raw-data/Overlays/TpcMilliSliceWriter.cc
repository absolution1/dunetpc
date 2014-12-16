#include "lbne-raw-data/Overlays/TpcMilliSliceWriter.hh"

lbne::TpcMilliSliceWriter::
TpcMilliSliceWriter(uint8_t* address, uint32_t max_size_bytes) :
    TpcMilliSlice(address), max_size_bytes_(max_size_bytes)
{
  header_()->version = 1;
  header_()->millislice_size = sizeof(Header);
  header_()->microslice_count = 0;
}

//std::shared_ptr<lbne::MicroSliceWriter> TODO
//lbne::TpcMilliSliceWriter::reserveMicroSlice(uint32_t ms_max_bytes)
//{
//  // finalize the most recent MicroSlice, in case that hasn't
//  // already been done
//  finalizeLatestMicroSlice_();
//
//  // test if this new MicroSlice could overflow our maximum size
//  if ((size() + ms_max_bytes) > max_size_bytes_) {
//    return false;
//  }
//
//  // create the next MicroSlice in our buffer, and update our
//  // counters to include the new MicroSlice
//  uint8_t* ms_ptr = data_(header_()->microslice_count);
//  latest_microslice_ptr_.reset(new MicroSliceWriter(ms_ptr, ms_max_bytes));
//  ++(header_()->microslice_count);
//  header_()->millislice_size += ms_max_bytes;
//  return latest_microslice_ptr_;
//}

int32_t lbne::TpcMilliSliceWriter::finalize(uint32_t data_size_bytes, uint32_t microslice_count)
{
  // first, we need to finalize the last MicroSlice, in case that
  // hasn't already been done
 // finalizeLatestMicroSlice_();


  header_()->millislice_size = sizeof(Header) + data_size_bytes;
  header_()->microslice_count = microslice_count;

  // next, we update our maximum size so that no more MicroSlices
  // can be added
  int32_t size_diff = max_size_bytes_ - header_()->millislice_size;
  max_size_bytes_ = header_()->millislice_size;
  return size_diff;
}

//void lbne::TpcMilliSliceWriter::finalizeLatestMicroSlice_()
//{
//  if (header_()->microslice_count > 0 &&
//      latest_microslice_ptr_.get() != 0) {
//    int size_change = latest_microslice_ptr_->finalize();
//    header_()->millislice_size -= size_change;
//  }
//}

lbne::TpcMilliSliceWriter::Header* lbne::TpcMilliSliceWriter::header_()
{
  return reinterpret_cast<Header *>(buffer_);
}

//uint8_t* lbne::TpcMilliSliceWriter::data_(int index) TODO
//{
//  uint8_t* ms_ptr = buffer_ + sizeof(Header);
//  for (int idx = 0; idx < index; ++idx) {
//    MicroSlice tmp_ms(ms_ptr);
//    ms_ptr += tmp_ms.size();
//  }
//  return ms_ptr;
//}
