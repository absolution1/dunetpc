#ifndef lbne_artdaq_Overlays_ToyFragmentWriter_hh
#define lbne_artdaq_Overlays_ToyFragmentWriter_hh

////////////////////////////////////////////////////////////////////////
// ToyFragmentWriter
//
// Class derived from ToyFragment which allows writes to the data (for
// simulation purposes). Note that for this reason it contains
// non-const members which hide the const members in its parent class,
// ToyFragment, including its reference to the artdaq::Fragment
// object, artdaq_Fragment_, as well as its functions pointing to the
// beginning and end of ADC values in the fragment, dataBegin() and
// dataEnd()
//
////////////////////////////////////////////////////////////////////////

#include "artdaq-core/Data/Fragment.hh"
#include "lbne-raw-data/Overlays/ToyFragment.hh"

#include <iostream>

namespace lbne {
  class ToyFragmentWriter;
}


class lbne::ToyFragmentWriter: public lbne::ToyFragment {
public:


  ToyFragmentWriter(artdaq::Fragment & f); 

  // These functions form overload sets with const functions from
  // lbne::ToyFragment

  adc_t * dataBegin();
  adc_t * dataEnd();

  // We'll need to hide the const version of header in ToyFragment in
  // order to be able to perform writes

  Header * header_() {
    assert(frag_.dataSizeBytes() >= sizeof(Header) );
    return reinterpret_cast<Header *>( artdaq_Fragment_.dataBeginBytes());
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
// ToyFragment::Metadata object, otherwise it throws

lbne::ToyFragmentWriter::ToyFragmentWriter(artdaq::Fragment& f ) :
  ToyFragment(f), artdaq_Fragment_(f) {
   
    // If this assert doesn't hold, then can't call
    // "words_to_frag_words_" below, translating between the
    // ToyFragment's standard data type size and the
    // artdaq::Fragment's data type size, on the Metadata object

    assert( sizeof(Metadata::data_t) == sizeof(Header::data_t) );

 
    if (! artdaq_Fragment_.hasMetadata() ) {
      throw cet::exception("Error in ToyFragmentWriter: passed artdaq::Fragment object doesn't contain metadata");
    }

    if (artdaq_Fragment_.dataSizeBytes() > 0) {
      throw cet::exception("Error in ToyFragmentWriter: passed artdaq::Fragment object already contains payload");
    }
 
    // Allocate space for the header
    artdaq_Fragment_.resizeBytes( sizeof(Header) );
}


inline lbne::ToyFragment::adc_t * lbne::ToyFragmentWriter::dataBegin() {
  assert(frag_.dataSizeBytes() > sizeof(Header) );
  return reinterpret_cast<adc_t *>(header_() + 1);
}

inline lbne::ToyFragment::adc_t * lbne::ToyFragmentWriter::dataEnd() {
  return dataBegin() + total_adc_values();
}


inline void lbne::ToyFragmentWriter::resize(size_t nAdcs) {
  auto es(calc_event_size_words_(nAdcs));
  artdaq_Fragment_.resizeBytes( sizeof(Header::data_t) * es );
  header_()->event_size = es;
}

inline size_t lbne::ToyFragmentWriter::calc_event_size_words_(size_t nAdcs) {
  return adcs_to_words_(nAdcs) + hdr_size_words();
}

inline size_t lbne::ToyFragmentWriter::adcs_to_words_(size_t nAdcs) {
  auto mod(nAdcs % adcs_per_word_());
  return (mod == 0) ?
    nAdcs / adcs_per_word_() :
    nAdcs / adcs_per_word_() + 1;
}

#endif /* lbne_artdaq_Overlays_ToyFragmentWriter_hh */
