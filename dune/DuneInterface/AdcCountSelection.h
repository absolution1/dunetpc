// AdcCountSelection.h

#ifndef AdcCountSelection_H
#define AdcCountSelection_H

// David Adams
// May 2016
//
// Class AdcSelection holds the channel number, pedestal, (floating) ADC count vector
// and filter vector for one ADC channel. The latter indicates which channels are
// to be retained.

#include "dune/DuneInterface/AdcTypes.h"

struct AdcCountSelection {

  // typedefs
  typedef unsigned int Channel;

  // public data
  const AdcCountVector& counts;
  Channel channel;
  AdcPedestal pedestal;
  AdcFilterVector filter;

  // Ctor from an ADC count vector. The filter is set to accept all.
  AdcCountSelection(const AdcCountVector& a_counts, Channel a_channel, AdcPedestal a_pedestal)
  : counts(a_counts), channel(a_channel), pedestal(a_pedestal), filter(counts.size(), true) { }

};

typedef std::vector<AdcCountSelection> AdcCountSelectionVector;

#endif
