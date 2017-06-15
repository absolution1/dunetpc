// AdcDataViewer.h

// David Adams
// June 2017
//
// Interface for tools that provide a view of ADC channel data.
// E.g. dump information to text stream or create plots.

#ifndef AdcDataViewer_H
#define AdcDataViewer_H

#include "dune/DuneInterface/AdcChannelData.h"

class AdcDataViewer {

public:

  virtual ~AdcDataViewer() =default;

  // Tool process method.
  //   acds - Input ADC data map
  //  label - Label to appear on plots
  //   fpat - String to be included in output file names
  virtual int view(const AdcChannelDataMap& acds,
                   std::string label ="",
                   std::string fpat ="") const =0;

};

#endif
