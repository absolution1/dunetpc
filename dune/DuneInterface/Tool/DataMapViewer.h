// DataMapViewer.h

// David Adams
// October 2017
//
// Interface for tools that provide a view of DataMap data, i.e.
// of the result from another tool.

#ifndef DataMapViewer_H
#define DataMapViewer_H

#include "dune/DuneInterface/Data/DataMap.h"

class DataMapViewer {

public:

  virtual ~DataMapViewer() =default;

  // Tool process method.
  //   acd - Input ADC channel data
  virtual DataMap view(const DataMap& dm) const =0;

};

#endif
