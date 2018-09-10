// StandardAdcChannelStringTool.h

// David Adams
// May 2018
//
// Tool to construct a string from an input template taking data
// from AdcChannelData and DataMap objects.
//
// The string is the pattern with the following replacements:
//   %RUN% --> acd.run
//   %SUBRUN% --> acd.subRun
//   %EVENT% --> acd.event
//   %CHAN% --> acd.channel
//   %FEMB% --> acd.fembID
//   %SUNIT% --> "sunit" where sunit = acd.sampleUnit
//   % SUNIT% --> " sunit" or "" if sunit is empty
//   %(SUNIT)% --> "(sunit)" or "" if sunit is empty
//   % (SUNIT)% --> " (sunit)" or "" if sunit is empty
//   %((SUNIT))% --> "(sunit)" or sunit if there is no space in sunit
//   % ((SUNIT))% --> add preceding space to this if sunit is not blank
//   %[SUNIT]% --> "[sunit]" or "" if sunit is empty
//   % [SUNIT]% --> " [sunit]" or "" if sunit is empty
//   %COUNT% --> dm.getInt("count") passed in call to build
//   %CHAN1% --> dm.getInt("chan1") passed in call to build
//   %CHAN2% --> dm.getInt("chan2") passed in call to build
// where acd is the AdcChannelData object and dm is the DataMap object
// passed in the call to build.
//
// E.g. if acd.run = 123, then "%RUN%" is replaced with "123".
//
// In addition, any number can be prefixed with a digit to indicate that
// the number should be padded with leading zeros to a width equal to that
// digit if its natural width is less that that value. E.g. %6RUN% will
// produce "000123" if acd.run = 123.
//
// If a number is prefixed with 0, then the width used is specified by
// the configuration here. E.g. if RunWidth == 5, then "%0RUN" will be
// replaced with "00123" if acd.run = 123.
//
// Tool configuration:
//  LogLevel - 1 to log from ctor
//             2 to log every call to build
//  RunWidth - width for run
//  SubRunWidth - width for subrun
//  EventWidth - width for event
//  ChannelWidth - width for channel
//  FembWidth - width for FEMB

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

  std::string build(const AdcChannelData& acd, const DataMap& dm, std::string spat) const override;

private:

  // Configuration data.
  int m_LogLevel;
  Index m_RunWidth;
  Index m_SubRunWidth;
  Index m_EventWidth;
  Index m_ChannelWidth;
  Index m_CountWidth;
  Index m_FembWidth;

  static const Index m_nrep = 8;
  Index m_wids[m_nrep];
  std::string m_reps[m_nrep];
  std::string m_bads[m_nrep];

};

DEFINE_ART_CLASS_TOOL(StandardAdcChannelStringTool)

#endif
