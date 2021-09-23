// AdcRoiShifter.h

// David Adams
// March 2019
//
// Tool to shift ROIS.
//
// Configuration:
//   LogLevel - 0=silent, 1=init, 2=each event, >2=more
//   BinOffset - # ticks to shift ROIs
//

#ifndef AdcRoiShifter_H
#define AdcRoiShifter_H

#include "art/Utilities/ToolMacros.h"
#include "fhiclcpp/ParameterSet.h"
#include "dune/DuneInterface/Tool/TpcDataTool.h"
#include <string>
#include <vector>

class AdcRoiShifter
: public TpcDataTool {

public:

  AdcRoiShifter(fhicl::ParameterSet const& ps);

  DataMap view(const AdcChannelData& acd) const override;

  DataMap update(AdcChannelData& acd) const override;

private:

  using Name = std::string;

  // Configuration data.
  int m_LogLevel;
  int m_BinOffset;

};


#endif
