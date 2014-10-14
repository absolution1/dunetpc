#include "lbne-raw-data/Overlays/ToyFragment.hh"

#include "cetlib/exception.h"

// namespace {
//   unsigned int pop_count (unsigned int n) {
//     unsigned int c; 
//     for (c = 0; n; c++) n &= n - 1; 
//     return c;
//   }
// }

void lbne::ToyFragment::checkADCData(int daq_adc_bits) const {
  lbne::ToyFragment::adc_t const * adcPtr(findBadADC(daq_adc_bits));
  if (adcPtr != dataEnd()) {
    throw cet::exception("IllegalADCVal")
        << "Illegal value of ADC word #"
        << (adcPtr - dataBegin())
        << ": 0x"
        << std::hex
        << *adcPtr
        << ".";
  }
}

std::ostream & lbne::operator << (std::ostream & os, ToyFragment const & f) {
  os << "ToyFragment event size: "
     << f.hdr_event_size()
     << ", run number: "
     << f.hdr_run_number()
     << "\n";

  return os;
}

