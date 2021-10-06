// ChannelStatusServiceTool.h
//
// David Adams
// November 2019
//
// Tool to return the channel status reported by ChannelStatusService.
//
// Parameters:
//   LogLevel - Message logging level (0=none, 1=ctor, 2=each call, ...)

#ifndef ChannelStatusServiceTool_H
#define ChannelStatusServiceTool_H

#include "art/Utilities/ToolMacros.h"
#include "fhiclcpp/ParameterSet.h"
#include "dune/DuneInterface/Tool/IndexMapTool.h"
#include <map>

namespace lariov {
  class ChannelStatusProvider;
}

class ChannelStatusServiceTool : public IndexMapTool {

public:

  using Index = IndexMapTool::Index;

  // Ctor.
  ChannelStatusServiceTool(fhicl::ParameterSet const& ps);

  // Dtor.
  ~ChannelStatusServiceTool() override =default;

  // Return the channel status.
  Index get(Index icha) const override;

private:

  // Parameters.
  Index m_LogLevel;

  // Derived from configuration.
  const lariov::ChannelStatusProvider* m_pChannelStatusProvider;

};


#endif
