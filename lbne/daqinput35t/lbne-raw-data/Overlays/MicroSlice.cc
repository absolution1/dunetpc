#include "lbne-raw-data/Overlays/MicroSlice.hh"

lbne::MicroSlice::MicroSlice(uint8_t* address) : buffer_(address)
{
}

lbne::MicroSlice::Header::microslice_size_t lbne::MicroSlice::size() const
{
  return header_()->microslice_size;
}

lbne::MicroSlice::Header::nanoslice_count_t lbne::MicroSlice::nanoSliceCount() const
{
  return header_()->nanoslice_count;
}

std::unique_ptr<lbne::NanoSlice> lbne::MicroSlice::nanoSlice(uint32_t index) const
{
  std::unique_ptr<NanoSlice> nslice_ptr;
  if (index < nanoSliceCount()) {
    uint8_t* ns_ptr = data_(index);
    nslice_ptr.reset(new NanoSlice(ns_ptr));
  }
  return nslice_ptr;
}

lbne::MicroSlice::Header const* lbne::MicroSlice::header_() const
{
  return reinterpret_cast<Header const *>(buffer_);
}

uint8_t* lbne::MicroSlice::data_(uint32_t index) const
{
  uint8_t* ns_ptr = reinterpret_cast<uint8_t *>(buffer_ + sizeof(Header));
  for (uint32_t idx = 0; idx < index; ++idx) {
    NanoSlice tmp_ns(ns_ptr);
    ns_ptr += tmp_ns.size();
  }
  return ns_ptr;
}
