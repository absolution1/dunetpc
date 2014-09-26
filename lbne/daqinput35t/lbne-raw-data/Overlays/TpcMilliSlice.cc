#include "lbne-raw-data/Overlays/TpcMilliSlice.hh"

lbne::TpcMilliSlice::TpcMilliSlice(uint8_t* address) : buffer_(address)
{
}

lbne::TpcMilliSlice::Header::millislice_size_t lbne::TpcMilliSlice::size() const
{
  return header_()->millislice_size;
}

lbne::TpcMilliSlice::Header::microslice_count_t lbne::TpcMilliSlice::microSliceCount() const
{
  return header_()->microslice_count;
}

std::unique_ptr<lbne::TpcMicroSlice> lbne::TpcMilliSlice::microSlice(uint32_t index) const
{
  std::unique_ptr<TpcMicroSlice> mslice_ptr;
  if (index < microSliceCount()) {
    uint8_t* ms_ptr = data_(index);
    mslice_ptr.reset(new TpcMicroSlice(ms_ptr));
  }
  return mslice_ptr;
}

lbne::TpcMilliSlice::Header const* lbne::TpcMilliSlice::header_() const
{
  return reinterpret_cast<Header const*>(buffer_);
}

uint8_t* lbne::TpcMilliSlice::data_(int index) const
{
  uint8_t* ms_ptr = buffer_ + sizeof(Header);
  for (int idx = 0; idx < index; ++idx) {
    TpcMicroSlice tmp_ms(ms_ptr);
    ms_ptr += tmp_ms.size();
  }
  return ms_ptr;
}
