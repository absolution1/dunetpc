#include "lbne-raw-data/Overlays/TpcMicroSlice.hh"
#include <cmath>
#include <iostream>

lbne::TpcMicroSlice::TpcMicroSlice(uint8_t* address) : buffer_(address)
{
}

lbne::TpcMicroSlice::Header::microslice_size_t lbne::TpcMicroSlice::size() const
{

  // Microslice size is in bottom 20 bits of header word 0, encoded as length in 32 bit words
  return (header_()->raw_header_data[0] & 0xFFFFF) * sizeof(uint32_t);
}

lbne::TpcMicroSlice::Header::record_type_t lbne::TpcMicroSlice::record_type() const
{
	uint8_t header_record_type = (header_()->raw_header_data[0] >> 24) & 0xF;
	return static_cast<lbne::TpcMicroSlice::Header::record_type_t>(header_record_type);
}


lbne::TpcMicroSlice::Header::record_format_t lbne::TpcMicroSlice::record_format() const
{
	uint8_t header_record_format = (header_()->raw_header_data[0] >> 28) & 0xF;
	return static_cast<lbne::TpcMicroSlice::Header::record_format_t>(header_record_format);
}

lbne::TpcMicroSlice::Header::timestamp_t lbne::TpcMicroSlice::timestamp() const
{
	// Timestamp originating from the frontend is in word 1 of header
	return reinterpret_cast<lbne::TpcMicroSlice::Header::timestamp_t>(header_()->raw_header_data[1]);
}

lbne::TpcMicroSlice::Header::data_subtype_t lbne::TpcMicroSlice::data_subtype() const
{
	// Subtype is in bits 3:0 of word 2 of header
	uint8_t header_subtype = (header_()->raw_header_data[2]) & 0xF;
	return static_cast<lbne::TpcMicroSlice::Header::data_subtype_t>(header_subtype);
}

lbne::TpcMicroSlice::Header::trigger_type_t lbne::TpcMicroSlice::trigger_type() const
{
	// Trigger type  is in bits 7:4 of word 2 of header
	uint8_t trigger_type = (header_()->raw_header_data[2] >> 4) & 0xF;
	return static_cast<lbne::TpcMicroSlice::Header::trigger_type_t>(trigger_type);
}

lbne::TpcMicroSlice::Header::data_format_t lbne::TpcMicroSlice::data_format() const
{
	// Data format is in bits 11:8 of word 2 of header
	uint8_t data_format = (header_()->raw_header_data[2] >> 8) & 0xF;
	return static_cast<lbne::TpcMicroSlice::Header::data_format_t>(data_format);
}

lbne::TpcMicroSlice::Header::group_count_t lbne::TpcMicroSlice::channelGroupCount() const
{
	// Channel group count is in nchans field (bits 19:16) of word 2 of header
	uint8_t header_nchans = (header_()->raw_header_data[2] >> 16) & 0xF;
	return static_cast<lbne::TpcMicroSlice::Header::group_count_t>(header_nchans);
}

lbne::TpcMicroSlice::Header::nanoslice_count_t lbne::TpcMicroSlice::nanoSliceCount() const
{
  // Number of nanoslices is encoded in ntimes field (bits 23:20) of word 2 of header, as 2**ntimes
  uint8_t header_ntimes = (header_()->raw_header_data[2] >> 20) & 0xF;
  return static_cast<lbne::TpcMicroSlice::Header::nanoslice_count_t>(pow(2, header_ntimes));
}

lbne::TpcMicroSlice::Header::group_physical_id_t lbne::TpcMicroSlice::groupPhysicalId(uint8_t group) const
{
	lbne::TpcMicroSlice::Header::group_physical_id_t* phys_id =
			reinterpret_cast<lbne::TpcMicroSlice::Header::group_physical_id_t*>((void*)&(header_()->raw_header_data[3]));
	return phys_id[group];
}

lbne::TpcMicroSlice::Header::sequence_id_t lbne::TpcMicroSlice::sequence_id() const
{
  // Sequence ID is word 5 of header
  return reinterpret_cast<lbne::TpcMicroSlice::Header::sequence_id_t>(header_()->raw_header_data[5]);
}

std::unique_ptr<lbne::TpcNanoSlice> lbne::TpcMicroSlice::nanoSlice(uint8_t group, uint32_t index) const
{
  uint32_t idx_offset = (index * lbne::TpcMicroSlice::Header::groups_per_microslice) + group;
  std::unique_ptr<TpcNanoSlice> nslice_ptr;
  if (index < nanoSliceCount()) {
    uint8_t* ns_ptr = data_(idx_offset);
    nslice_ptr.reset(new TpcNanoSlice(ns_ptr));
  }
  return nslice_ptr;
}

lbne::TpcMicroSlice::Header const* lbne::TpcMicroSlice::header_() const
{
  return reinterpret_cast<Header const *>(buffer_);
}

uint8_t* lbne::TpcMicroSlice::data_(uint32_t index) const
{
  uint8_t* ns_ptr = reinterpret_cast<uint8_t *>(buffer_ + sizeof(Header));
  for (uint32_t idx = 0; idx < index; ++idx) {
    TpcNanoSlice tmp_ns(ns_ptr);
    ns_ptr += tmp_ns.size();
  }
  return ns_ptr;
}
