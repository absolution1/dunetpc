// AdcChannelDumper.h

// David Adams
// August 2017
//
// Tool to dump information about an ADC channel.
//
// Configuration:
//   FileName:  Name of the output file. Blank for std out.
//   Prefix:    Prefix for each line.
//   NewFile:   If true, a new file is created for each dump.
//   MaxSample: Maximum # sample values to display. -1 for all.

#ifndef AdcChannelDumper_H
#define AdcChannelDumper_H

#include "art/Utilities/ToolMacros.h"
#include "fhiclcpp/ParameterSet.h"
#include "dune/DuneInterface/Tool/AdcChannelTool.h"
#include <iostream>

class AdcChannelDumper : AdcChannelTool {

public:

  AdcChannelDumper(fhicl::ParameterSet const& ps);

  ~AdcChannelDumper() override;

  DataMap view(const AdcChannelData& acd) const override;
  bool updateWithView() const override { return true; }

private:

  // Configuration data.
  std::string m_FileName;
  std::string m_Prefix;
  bool m_NewFile;
  int  m_MaxSample;

  // Output stream.
  std::ostream* m_pout;

};

DEFINE_ART_CLASS_TOOL(AdcChannelDumper)

#endif
