// AcdRoiBuilder.h

// Christoph Alt
// October 2017
//
// Interface for tools that find signals and build ROI's


#ifndef AcdRoiBuilder_H
#define AcdRoiBuilder_H

#include <iostream>
#include "dune/DuneInterface/AdcChannelData.h"

class AcdRoiBuilder {

public:

  virtual ~AcdRoiBuilder() =default;

  // Build ROIs and set data.rois correspondingly.
  virtual int build(AdcChannelData& data) const =0;

  // Print parameters.
  virtual std::ostream& print(std::ostream& out =std::cout, std::string prefix ="") const =0;

};

#endif
