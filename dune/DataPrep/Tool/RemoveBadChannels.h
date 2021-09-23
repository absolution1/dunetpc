// RemoveBadChannels.h
//
// Tool to remove bad channels
//
// Configuration:
//   LogLevel: 0=silent, 1=init, 2=each event, >2=more

#ifndef RemoveBadChannels_H
#define RemoveBadChannels_H

#include "art/Utilities/ToolMacros.h"
#include "fhiclcpp/ParameterSet.h"
#include "dune/DuneInterface/Tool/TpcDataTool.h"

namespace lariov {
  class ChannelStatusProvider;
}

class RemoveBadChannels : TpcDataTool {

public:

  RemoveBadChannels(fhicl::ParameterSet const& ps);

  ~RemoveBadChannels() override =default;

  DataMap update(AdcChannelData& acd) const override;

private:

  // Configuration data.
  int                m_LogLevel;
  bool               m_RemoveBadChs;
  bool               m_RemoveNoisyChs;

  // Channel status provider.
  const lariov::ChannelStatusProvider* m_pChannelStatusProvider;

};


#endif
