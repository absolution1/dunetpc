#include "lbne-raw-data/Overlays/NanoSliceWriter.hh"

lbne::NanoSliceWriter::
NanoSliceWriter(uint8_t* address, uint32_t max_size_bytes,
                uint16_t channel_number) :
    NanoSlice(address), max_size_bytes_(max_size_bytes)
{
  header_()->version = 1;
  header_()->nanoslice_size = sizeof(Header);
  header_()->channel_number = channel_number;
  header_()->sample_count = 0;

  if (max_size_bytes_ < size()) {
    max_size_bytes_ = size();
  }
}

void lbne::NanoSliceWriter::setChannelNumber(uint16_t channel)
{
  header_()->channel_number = channel;
}

bool lbne::NanoSliceWriter::addSample(uint16_t value)
{
  if ((size() + sizeof(value)) > max_size_bytes_) {
    return false;
  }
  data_()[header_()->sample_count] = value;
  ++(header_()->sample_count);
  header_()->nanoslice_size += sizeof(value);
  return true;
}

int32_t lbne::NanoSliceWriter::finalize()
{
  int32_t size_diff = max_size_bytes_ - header_()->nanoslice_size;
  max_size_bytes_ = header_()->nanoslice_size;
  return size_diff;
}

lbne::NanoSlice::Header * lbne::NanoSliceWriter::header_()
{
  return reinterpret_cast<Header *>(buffer_);
}

uint16_t * lbne::NanoSliceWriter::data_()
{
  return reinterpret_cast<uint16_t *>(buffer_ + sizeof(Header));
}
