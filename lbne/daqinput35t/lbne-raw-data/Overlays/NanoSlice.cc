#include "lbne-raw-data/Overlays/NanoSlice.hh"

lbne::NanoSlice::NanoSlice(uint8_t* address) : buffer_(address) { }

lbne::NanoSlice::Header::nanoslice_size_t lbne::NanoSlice::size() const
{
  return header_()->nanoslice_size;
}

lbne::NanoSlice::Header::channel_number_t lbne::NanoSlice::channelNumber() const
{
  return header_()->channel_number;
}

lbne::NanoSlice::Header::sample_count_t lbne::NanoSlice::sampleCount() const
{
  return header_()->sample_count;
}

bool lbne::NanoSlice::sampleValue(uint32_t index, uint16_t& value) const
{
  if (index >= sampleCount()) {return false;}
  value = data_()[index];
  return true;
}

lbne::NanoSlice::Header const* lbne::NanoSlice::header_() const
{
  return reinterpret_cast<Header const*>(buffer_);
}

uint16_t const* lbne::NanoSlice::data_() const
{
  return reinterpret_cast<uint16_t const*>(buffer_ + sizeof(Header));
}
