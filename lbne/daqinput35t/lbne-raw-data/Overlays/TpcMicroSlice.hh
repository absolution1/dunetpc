#ifndef lbne_artdaq_Overlays_TpcMicroSlice_hh
#define lbne_artdaq_Overlays_TpcMicroSlice_hh

#include "lbne-raw-data/Overlays/TpcNanoSlice.hh"
#include <memory>

namespace lbne {
  class TpcMicroSlice;
}

class lbne::TpcMicroSlice {

public:

  struct Header {

	// Generic raw header word type
    typedef uint32_t data_t;

    // Raw header word 0 fields
    typedef uint32_t microslice_size_t;
    typedef uint8_t  record_type_t;
    typedef uint8_t  record_format_t;

    // Raw header word 1 fields
    typedef uint32_t timestamp_t;

    // Raw header word 2 fields
    typedef uint8_t  data_subtype_t;
    typedef uint8_t  trigger_type_t;
    typedef uint8_t  data_format_t;
    typedef uint8_t  nchans_t;
    typedef uint8_t  ntimes_t;

    // Raw header words 3-4 fields
    typedef uint16_t group_physical_id_t;

    // Raw header word 5 fields
    typedef uint32_t sequence_id_t;

    // Parameters derived from header field translation
    typedef uint8_t  group_count_t;
    typedef uint32_t nanoslice_count_t;

    static constexpr size_t raw_header_words = 6;
    static constexpr uint8_t groups_per_microslice = 4;

    // Raw header words
    data_t raw_header_data[raw_header_words];

  };

  // This constructor accepts a memory buffer that contains an existing
  // TpcMicroSlice and allows the the data inside it to be accessed
  TpcMicroSlice(uint8_t* address);

  // Returns the size of the TpcMicroSlice
  Header::microslice_size_t size() const;

  // Returns the record type of the TpcMicroSlice
  Header::record_type_t record_type() const;

  // Returns the record format of the TpcMicroSlice
  Header::record_format_t record_format() const;

  // Returns the timestamp of the TpcMicroSlice
  Header::timestamp_t timestamp() const;

  // Returns the data subtype in the TpcMicroSlice
  Header::data_subtype_t data_subtype() const;

  // Returns the trigger type in the TpcMicroSlice
  Header::trigger_type_t trigger_type() const;

  // Return the data format in the TpcMicroSlice
  Header::data_format_t data_format() const;

  // Returns the number of groups of channels in the TpcMicroSlice
  Header::group_count_t channelGroupCount() const;

  // Returns the number of NanoSlices in the TpcMicroSlice
  Header::nanoslice_count_t nanoSliceCount() const;

  // Returns the physical ID of the specified channel group
  Header::group_physical_id_t groupPhysicalId(uint8_t group) const;

  // Returns the sequence ID of the TpcMicroSlice
  Header::sequence_id_t sequence_id() const;

  // Returns the requested NanoSlice if the requested slice was found,
  // otherwise returns an empty pointer
  std::unique_ptr<TpcNanoSlice> nanoSlice(uint8_t group, uint32_t index) const;

protected:

  // returns a pointer to the header
  Header const* header_() const;

  // returns a pointer to the requested NanoSlice
  uint8_t* data_(uint32_t index) const;

  uint8_t* buffer_;
};

#endif /* lbne_artdaq_Overlays_TpcMicroSlice_hh */
