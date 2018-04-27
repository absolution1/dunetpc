// StandardAdcChannelStringTool.h

// David Adams
// April 20187
//
// Tool to construct a string from AdcChannelData and input pattern.
//
// The string is the pattern with the following replacements:
//   %RUN% --> acd.run
//   %SUBRUN% --> acd.subRun
//   %EVENT% --> acd.event
//   %CHAN% --> acd.channel
//   %CHAN% --> acd.channel
//   %COUNT% --> count passed in call to build
// where w can be absent or any unsigned int.
// If w is 0 or blank, the value is written with its natural width.
// Otherwise it is written with that with padded with leading zeroes.
//
// Tool configuration:
//  LogLevel - 1 to log from ctor
//             2 to log every call to build
//  RunWidth - width for run
//  SubRunWidth - width for subrun
//  EventWidth - width for event
//  ChannelWidth - width for channel

#ifndef StandardAdcChannelStringTool_H
#define StandardAdcChannelStringTool_H

#include "art/Utilities/ToolMacros.h"
#include "fhiclcpp/ParameterSet.h"
#include "dune/DuneInterface/Tool/AdcChannelStringTool.h"

class StandardAdcChannelStringTool
: public AdcChannelStringTool {

public:

  using Index = unsigned int;

  StandardAdcChannelStringTool(fhicl::ParameterSet const& ps);

  std::string build(const AdcChannelData& acd, std::string sfmt, Index count) const override;

private:

  // Configuration data.
  int m_LogLevel;
  Index m_RunWidth;
  Index m_SubRunWidth;
  Index m_EventWidth;
  Index m_ChannelWidth;
  Index m_CountWidth;

  static const Index m_nrep = 5;
  Index m_wids[m_nrep];
  std::string m_reps[m_nrep];
  std::string m_bads[m_nrep];

};

DEFINE_ART_CLASS_TOOL(StandardAdcChannelStringTool)

#endif
