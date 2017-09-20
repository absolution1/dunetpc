// AdcPedestalFitter.h

// David Adams
// August 2017
//
// Tool to fit ADC data and extract a pedestal.
//
// Configuration:
//   LogLevel - 0=silent, 1=init, 2=each event, >2=more
//   HistName:  Name for the histogram.
//   HistTitle: Title for the histogram.
//   HistManager: Name of the tool that manages the output histogram.
//                Right now this is the only way to retrive the fit result.
// The following subsitutions are made in the names:
//    %RUN% - run number
//    %SUBRUN% - event number
//    %EVENT% - event number
//    %CHAN% - channel number

#ifndef AdcPedestalFitter_H
#define AdcPedestalFitter_H

#include "art/Utilities/ToolMacros.h"
#include "fhiclcpp/ParameterSet.h"
#include "dune/DuneInterface/Tool/AdcChannelViewer.h"
#include "dune/DuneInterface/Tool/AdcChannelDataModifier.h"
#include <string>
#include <vector>

class HistogramManager;
class TH1;

class AdcPedestalFitter
: public AdcChannelViewer,
  public AdcChannelDataModifier {

public:

  AdcPedestalFitter(fhicl::ParameterSet const& ps);

  int view(const AdcChannelData& acd) const override;

  int update(AdcChannelData& acd) const override;

private:

  using Name = std::string;
  using NameVector = std::vector<Name>;

  // Configuration data.
  int m_LogLevel;
  Name m_HistName;
  Name m_HistTitle;
  Name m_HistManager;

  // Histogram manager.
  HistogramManager* m_phm;

  // Make replacements in a name.
  Name nameReplace(Name name, const AdcChannelData& acd) const;

  // Find and return pedestal.
  struct Result {
    Result(int rstat) : stat(rstat) { }
    int stat = -99;
    double pedestal = 0.0;
    TH1* ph = nullptr;
  };
  Result getPedestal(const AdcChannelData& acd) const;

};

DEFINE_ART_CLASS_TOOL(AdcPedestalFitter)

#endif
