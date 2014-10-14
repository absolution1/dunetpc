#include "lbne-raw-data/Overlays/TpcNanoSlice.hh"
#include <iostream>

lbne::TpcNanoSlice::TpcNanoSlice(uint8_t* address) : buffer_(address) { }

lbne::TpcNanoSlice::nanoslice_size_t lbne::TpcNanoSlice::size() const
{
  return (lbne::TpcNanoSlice::nanoslice_size_t)(num_nanoslice_words * sizeof(raw_data_word_t));
}

lbne::TpcNanoSlice::sample_count_t lbne::TpcNanoSlice::sampleCount() const
{
  return (lbne::TpcNanoSlice::sample_count_t)(num_nanoslice_adcs);
}

bool lbne::TpcNanoSlice::sampleValue(uint32_t index, uint16_t& value) const
{
  if (index >= sampleCount()) {return false;}

  unsigned int highNibbleOffset = (index / 2) + num_nanoslice_adcs;
  unsigned int highNibbleShift  = (index % 2) * 4;

  value =  static_cast<uint16_t>(buffer_[index]) |
		  (static_cast<uint16_t>(((buffer_[highNibbleOffset] >> highNibbleShift) & 0xF) << 8));

  return true;
}

uint16_t const* lbne::TpcNanoSlice::data_() const
{
  return reinterpret_cast<uint16_t const*>(buffer_);
}
