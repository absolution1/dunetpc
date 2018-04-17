// AdcDataDumper.h

// David Adams
// June 2017
//
// Tool to dump information about an ADC channel data map.
//
// Configuration:
//   FileName: Name of the output file. Blank for std out.
//             The following substitutions are made:
//               %PAT% --> the pattern passed by the caller
//               %CHAN1% --> the first channel number
//               %CHAN2% --> the last channel number
//   Prefix:   Prefix for each line.
//   NewFile:  If true, a new file is created for each event.
//   ShowChannelCount: If true, display channel count
//   ShowTickCounts: Show the numbers of raw, sample flag, siginal entries
//                   Also show the number of ROIs.
//   ShowRaw - Show the raw data.
//   ShowPrepared - Show the prepared (pedstal removed, ...) data.
//   ShowFirst - First tick number to display
//   ShowRebin - Rebinning factor for show.
//               Displayed values are averages over this number of ticks.
//   ShowMax - Maximum # values to display
//   ShowThreshold - Threshold for ShowOpt = 2
//   ShowOpt - if 2, then symbols are displayed instead of values
//             * for above threshold, - for below -threshold, . otherwise

#ifndef AdcDataDumper_H
#define AdcDataDumper_H

#include "art/Utilities/ToolMacros.h"
#include "fhiclcpp/ParameterSet.h"
#include "dune/DuneInterface/Tool/AdcChannelTool.h"
#include <iostream>

class AdcDataDumper : AdcChannelTool {

public:

  AdcDataDumper(fhicl::ParameterSet const& ps);

  ~AdcDataDumper() override;

  DataMap viewMap(const AdcChannelDataMap& acds) const override;
  bool updateWithView() const override { return true; }

private:

  // Configuration data.
  std::string m_FileName;
  std::string m_Prefix;
  bool m_NewFile;
  bool m_ShowChannelCount;
  bool m_ShowTickCounts;
  bool m_ShowRaw;
  bool m_ShowPrepared;
  unsigned int m_ShowFirst;
  unsigned int m_ShowRebin;
  unsigned int m_ShowMax;
  float m_ShowThreshold;
  unsigned int m_ShowOpt;

  // Output stream.
  std::ostream* m_pout;

};

DEFINE_ART_CLASS_TOOL(AdcDataDumper)

#endif
