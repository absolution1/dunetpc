#ifndef lbne_artdaq_Overlays_TpcNanoSlice_hh
#define lbne_artdaq_Overlays_TpcNanoSlice_hh

#include <stdint.h>

namespace lbne {
  class TpcNanoSlice;
}

class lbne::TpcNanoSlice {

public:

//  struct Header {
//    typedef uint32_t data_t;
//
//    typedef uint16_t version_t;
//    typedef uint16_t nanoslice_size_t;
//    typedef uint16_t channel_number_t;
//    typedef uint16_t sample_count_t;
//
//    // this structure is deliberately incomplete - it should be
//    // improved to match what we want for the LBNE 35t detector
//
//    data_t version        : 16;
//    data_t nanoslice_size : 16;   // total size, data & header
//    data_t channel_number : 16;
//    data_t sample_count   : 16;
//  };

  typedef uint32_t raw_data_word_t;
  typedef uint16_t nanoslice_size_t;
  typedef uint16_t sample_count_t;
  static constexpr uint8_t num_nanoslice_adcs = 32;
  static constexpr uint8_t num_nanoslice_words = 12;

  // This constructor accepts a memory buffer that contains an existing
  // nanoSlice and allows the the data inside it to be accessed
  TpcNanoSlice(uint8_t* address);

  // Returns the size of the TpcNanoSlice
  lbne::TpcNanoSlice::nanoslice_size_t size() const;

  // Returns the channel number for the TpcNanoSlice
  //Header::channel_number_t channelNumber() const;

  // Returns the number of samples in this TpcNanoSlice
  lbne::TpcNanoSlice::sample_count_t sampleCount() const;

  // Fetches the value for the requested sample.  Returns true if
  // the requested value was found, false if not.
  bool sampleValue(uint32_t index, uint16_t& value) const;

protected:

  // returns a pointer to the header
  //Header const* header_() const;

  // returns a pointer to the first sample value
  uint16_t const* data_() const;

  uint8_t* buffer_;
};

#endif /* lbne_artdaq_Overlays_TpcNanoSlice_hh */
