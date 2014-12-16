#ifndef artdaq_lbne_Overlays_SSPFragmentWriter_hh
#define artdaq_lbne_Overlays_SSPFragmentWriter_hh

////////////////////////////////////////////////////////////////////////
// SSPFragmentWriter
//
// Class derived from SSPFragment which allows writes to the data (for
// simulation purposes). Note that for this reason it contains
// non-const members which hide the const members in its parent class,
// SSPFragment, including its reference to the artdaq::Fragment
// object, artdaq_Fragment_, as well as its functions pointing to the
// beginning and end of ADC values in the fragment, dataBegin() and
// dataEnd()
//
////////////////////////////////////////////////////////////////////////

#include "artdaq-core/Data/Fragment.hh"
#include "lbne-raw-data/Overlays/SSPFragment.hh"

#include <iostream>

namespace lbne {
  class SSPFragmentWriter;
}


class lbne::SSPFragmentWriter: public lbne::SSPFragment {
public:


  SSPFragmentWriter(artdaq::Fragment & f); 

  // These functions form overload sets with const functions from
  // lbne::SSPFragment

  adc_t * dataBegin();
  adc_t * dataEnd();

  // We'll need to hide the const version of header in SSPFragment in
  // order to be able to perform writes

  Header * header_() {
    assert(frag_.dataSizeBytes() >= sizeof(Header) );
    return reinterpret_cast<Header *>(artdaq_Fragment_.dataBeginBytes());
  }

  void set_hdr_run_number(Header::run_number_t run_number) { 
    header_()->run_number = run_number;
  }

  void resize(size_t nAdcs);

private:
  size_t calc_event_size_words_(size_t nAdcs);

  static size_t adcs_to_words_(size_t nAdcs);

  // Note that this non-const reference hides the const reference in the base class
  artdaq::Fragment & artdaq_Fragment_;
};

// The constructor will expect the artdaq::Fragment object it's been
// passed to contain the artdaq::Fragment header + the
// SSPFragment::Metadata object, otherwise it throws

lbne::SSPFragmentWriter::SSPFragmentWriter(artdaq::Fragment& f ) :
  SSPFragment(f), artdaq_Fragment_(f) {
   
  if ( ! f.hasMetadata() || f.dataSizeBytes() > 0 ) {
    throw cet::exception("Error in SSPFragmentWriter: Raw artdaq::Fragment object does not appear to consist of (and only of) its own header + the SSPFragment::Metadata object");
  }
 
  // Allocate space for the header
  artdaq_Fragment_.resizeBytes( sizeof(Header) );
}


inline lbne::SSPFragment::adc_t * lbne::SSPFragmentWriter::dataBegin() {
  // Make sure there's data past the SSPFragment header
  assert(frag_.dataSizeBytes() >= sizeof(Header) + sizeof(artdaq::Fragment::value_type) );
  return reinterpret_cast<adc_t *>(header_() + 1);
}

inline lbne::SSPFragment::adc_t * lbne::SSPFragmentWriter::dataEnd() {
  return dataBegin() + total_adc_values();
}


inline void lbne::SSPFragmentWriter::resize(size_t nAdcs) {

  artdaq_Fragment_.resizeBytes(sizeof(Header::data_t) * calc_event_size_words_(nAdcs) );
  header_()->event_size = calc_event_size_words_(nAdcs);
}

inline size_t lbne::SSPFragmentWriter::calc_event_size_words_(size_t nAdcs) {
  return adcs_to_words_(nAdcs) + hdr_size_words();
}

inline size_t lbne::SSPFragmentWriter::adcs_to_words_(size_t nAdcs) {
  auto mod(nAdcs % adcs_per_word_());
  return (mod == 0) ?
    nAdcs / adcs_per_word_() :
    nAdcs / adcs_per_word_() + 1;
}

#endif /* artdaq_lbne_Overlays_SSPFragmentWriter_hh */
