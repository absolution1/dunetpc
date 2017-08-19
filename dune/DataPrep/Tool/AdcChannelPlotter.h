// AdcChannelPlotter.h

// David Adams
// August 2017
//
// Tool to plot data from an ADC channel.
//
// Configuration:
//   LogLevel - 0=silent, 1=init, 2=each event, >2=more
//   HistTypes: Types of histograms to create:
//                raw = raw data
//                prepared = prepared data
//   HistName:  Name for the histogram.
//   HistTitle: Title for the histogram.
//   RootFileName: If non-blank, histograms are written to this file.
//                 File is opened in UPDATE mode.
// The following subsitutions are made in the names:
//    %RUN% - run number
//    %SUBRUN% - event number
//    %EVENT% - event number
//    %CHAN% - channel number
//    %TYPE% - histogram type (see HistTypes)

#ifndef AdcChannelPlotter_H
#define AdcChannelPlotter_H

#include "art/Utilities/ToolMacros.h"
#include "fhiclcpp/ParameterSet.h"
#include "dune/DuneInterface/Tool/AdcChannelViewer.h"
#include <string>
#include <vector>

class AdcChannelPlotter : AdcChannelViewer {

public:

  AdcChannelPlotter(fhicl::ParameterSet const& ps);

  int view(const AdcChannelData& acd) const override;

private:

  using Name = std::string;
  using NameVector = std::vector<Name>;

  // Configuration data.
  int m_LogLevel;
  NameVector m_HistTypes;
  Name m_HistName;
  Name m_HistTitle;
  Name m_RootFileName;

  // Output stream.
  std::ostream* m_pout;

  // Make replacements in a name.
  Name nameReplace(Name name, const AdcChannelData& acd, Name type) const;

};

DEFINE_ART_CLASS_TOOL(AdcChannelPlotter)

#endif
