#include "lbne-raw-data/Overlays/MilliSlice.hh"

lbne::MilliSlice::MilliSlice(uint8_t* address) : buffer_(address)
{
}

lbne::MilliSlice::Header::millislice_size_t lbne::MilliSlice::size() const
{
  return header_()->millislice_size;
}

lbne::MilliSlice::Header::microslice_count_t lbne::MilliSlice::microSliceCount() const
{
  return header_()->microslice_count;
}

std::unique_ptr<lbne::MicroSlice> lbne::MilliSlice::microSlice(uint32_t index) const
{
  std::unique_ptr<MicroSlice> mslice_ptr;
  if (index < microSliceCount()) {
    uint8_t* ms_ptr = data_(index);
    mslice_ptr.reset(new MicroSlice(ms_ptr));
  }
  return mslice_ptr;
}

lbne::MilliSlice::Header const* lbne::MilliSlice::header_() const
{
  return reinterpret_cast<Header const*>(buffer_);
}

uint8_t* lbne::MilliSlice::data_(int index) const
{
  uint8_t* ms_ptr = buffer_ + sizeof(Header);
  for (int idx = 0; idx < index; ++idx) {
    MicroSlice tmp_ms(ms_ptr);
    ms_ptr += tmp_ms.size();
  }
  return ms_ptr;
}
