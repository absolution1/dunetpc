// AdcChannelDataModifier.h

// David Adams
// August 2017
//
// Interface for tools that modify the data for an ADC channel.
//
// E.g. read or write LarSoft digits or wires, add or remove
// noise or convert raw to/from prepared.

#ifndef AdcChannelDataModifier_H
#define AdcChannelDataModifier_H

#include "dune/DuneInterface/Tool/AdcChannelViewer.h"
#include "dune/DuneInterface/AdcChannelData.h"

class AdcChannelDataModifier : public AdcChannelViewer {

public:

  virtual ~AdcChannelDataModifier() =default;

  // Inherited method.
  // Default implementation does nothing and returns success.
  DataMap view(const AdcChannelData&) const override {return DataMap(0);}

  // Tool process method.
  //   acd - ADC channel data
  virtual DataMap update(AdcChannelData& acd) const =0;

};

#endif
