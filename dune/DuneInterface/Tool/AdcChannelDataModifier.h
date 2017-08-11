// AdcChannelDataModifier.h

// David Adams
// August 2017
//
// Interface for tools that modify ADC channel data
// independently for each channel.
//
// E.g. read or write LarSoft digits or wires, add or remove
// noise or convert raw to/from prepared.

#ifndef AdcChannelDataModifier_H
#define AdcChannelDataModifier_H

#include "dune/DuneInterface/AdcChannelData.h"

class AdcChannelDataModifier {

public:

  virtual ~AdcChannelDataModifier() =default;

  // Tool process method.
  //   acd - ADC channel data
  virtual int update(AdcChannelData& acd) const =0;

};

#endif
