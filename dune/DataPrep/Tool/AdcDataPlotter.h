// AdcDataPlotter.h

// David Adams
// July 2017
//
// Tool to make event displays of prepared data from an ADC channel data map.
//
// Configuration:
//   FileName: Name of the output file. Blank for std out.
//             The following substitutions are made:
//               %PAT% --> the pattern passed by the caller
//               %CHAN1% --> the first channel number
//               %CHAN2% --> the last channel number
//   FirstTick - First tick number to display
//   LastTick - Last+1 tick number to display
//   MaxSignal - Displayed signal range is (-MaxSignal, MaxSignal)

#ifndef AdcDataPlotter_H
#define AdcDataPlotter_H

#include "art/Utilities/ToolMacros.h"
#include "fhiclcpp/ParameterSet.h"
#include "dune/DuneInterface/Tool/AdcDataViewer.h"
#include <iostream>

class AdcDataPlotter : AdcDataViewer {

public:

  AdcDataPlotter(fhicl::ParameterSet const& ps);

  ~AdcDataPlotter() override =default;

  int view(const AdcChannelDataMap& acds,
           std::string label ="",
           std::string fpat ="") const override;

private:

  // Configuration data.
  std::string m_FileName;
  unsigned long m_FirstTick;
  unsigned long m_LastTick;
  double m_MaxSignal;

  // Output stream.
  std::ostream* m_pout;

};

DEFINE_ART_CLASS_TOOL(AdcDataPlotter)

#endif
