// AdcChannelViewer.h

// David Adams
// August 2017
//
// Interface for tools that provide a view of the data for an ADC channel.
// E.g. dump information to text stream or create plots.

#ifndef AdcChannelViewer_H
#define AdcChannelViewer_H

#include "dune/DuneInterface/Data/DataMap.h"
#include "dune/DuneInterface/AdcChannelData.h"

class AdcChannelViewer {

public:

  virtual ~AdcChannelViewer() =default;

  // Tool process method.
  //   acd - Input ADC channel data
  virtual DataMap view(const AdcChannelData& acd) const =0;

};

#endif
