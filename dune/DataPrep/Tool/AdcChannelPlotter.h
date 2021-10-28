// AdcChannelPlotter.h

// David Adams
// August 2017
//
// Tool to plot data from an ADC channel.
//
// Either waveforms (signal vs. tick) or distributions (# tick vs signal)
// may be produced.
//
// Configuration:
//   LogLevel - 0=silent, 1=init, 2=each event, >2=more
//   HistTypes: Types of histograms to create:
//                raw = raw data: ADC vs tick
//                rawdist = raw data dist: # ticks vs ADC
//                rawdistlog = raw data dist: log(# ticks) vs ADC
//                prepared = prepared data: signal vs. tick
//   HistName:  Name for the histogram.
//   HistTitle: Title for the histogram.
//   RootFileName: If non-blank, histograms are written to this file.
//                 File is opened in UPDATE mode.
//   PlotFileName: Name of the file to which plots should be saved.
//   PlotSamMin: Min tick for raw and prepared plots.
//   PlotSamMax: Max tick for raw and prepared plots.
//   PlotSigOpt: Option for setting the plotted signal range:
//                  full - Full range determined from all ticks. Expanded to PlotSigMin ticks.
//                  fixed - Fixed range specified (PlotSigMin, PlotSigMax)
//                  pedestal - Fixed range around pedestal (PlotSigMin+ped, PlotSigMax+ped)
//   PlotSigMin: - Min for signal range. See PlotSigOpt.
//   PlotSigMax: - Max for signal range. See PlotSigOpt.
//   PlotDistMin - Min for y-axis in rawdist plots
//   PlotDistMax - Max for y-axis in rawdist plots
//   ColorBad - If nonzero, color for channels flagged bad.
//   ColorNoisy - If nonzero, color for channels flagged noisy.
//   SkipFlags - Samples with these flags are excluded from dist plots
//   LabelSize - Size for x and y labels and titles. This is fraction of the pad size.
//               Default of 0 ==> Root's 0.035 which can be tiny for many vertical suppads.
// The following subsitutions are made in the names:
//    %RUN% - run number
//    %SUBRUN% - event number
//    %EVENT% - event number
//    %CHAN% - channel number
//    %TYPE% - histogram type (see HistTypes)
//
// Datamapp returned by call to view includes:
//     TH1   HHH  - The waveform histogram (HHH = HistName in configuration)
//   float   plotSigMin_HHH - min value for the histogram
//   float   plotSigMax_HHH - max value for the histogram

#ifndef AdcChannelPlotter_H
#define AdcChannelPlotter_H

#include "art/Utilities/ToolMacros.h"
#include "fhiclcpp/ParameterSet.h"
#include "dune/DuneInterface/Tool/TpcDataTool.h"
#include <string>
#include <vector>
#include <map>
#include <set>

class AdcChannelStringTool;
namespace lariov {
  class ChannelStatusProvider;
}
class TH1;

class AdcChannelPlotter : TpcDataTool {

public:

  AdcChannelPlotter(fhicl::ParameterSet const& ps);

  ~AdcChannelPlotter();

  DataMap view(const AdcChannelData& acd) const override;
  DataMap viewMap(const AdcChannelDataMap& acds) const override;
  bool updateWithView() const override { return true; }

private:

  using Name = std::string;
  using NameVector = std::vector<Name>;
  using Index = unsigned int;
  using IndexVector = std::vector<Index>;
  using IndexSet = std::set<Index>;

  // Configuration data.
  int m_LogLevel;
  NameVector m_HistTypes;
  Name m_HistName;
  Name m_HistTitle;
  Name m_RootFileName;
  Name m_PlotFileName;
  Index m_PlotSamMin;      // Tick range to plot.
  Index m_PlotSamMax;
  Name m_PlotSigOpt;
  float m_PlotSigMin;
  float m_PlotSigMax;
  float m_PlotDistMin;
  float m_PlotDistMax;
  Index m_ColorBad;
  Index m_ColorNoisy;
  float m_LabelSize;
  IndexVector m_SkipFlags;

  // ADC string tool.
  const AdcChannelStringTool* m_adcStringBuilder;

  // Channel status provider.
  const lariov::ChannelStatusProvider* m_pChannelStatusProvider;

  // Derived from config.
  IndexSet m_skipFlags;

  // Make replacements in a name.
  Name nameReplace(Name name, const AdcChannelData& acd, Name type) const;

  using HistMap = std::map<Name, TH1*>;

  class State {
  public:
    HistMap hists;
  };
  mutable State m_state;
  State& getState() const { return m_state; }

};


#endif
