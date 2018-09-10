// AdcPedestalFitter.h

// David Adams
// August 2017
//
// Tool to fit ADC data and extract a pedestal.
//
// If FitRmsMin < FitRmsMax, the the RMS is constrained to the range
// (FitRmsMin, FitRmsMax) in the fit.
//
// Configuration:
//   LogLevel - 0=silent, 1=init, 2=each event, >2=more
//   FitRmsMin: Lower limit for RMS fit range.
//   FitRmsMax: Upper limit for RMS fit range.
//   HistName:  Name for the histogram.
//   HistTitle: Title for the histogram.
//   HistManager: Name of the tool that manages the output histogram.
//                This is obsolete.
//   PlotFileName: If nonblank, histograms are displayed in this file.
//   RootFileName: If nonblank, histogram is copied to this file.
//   PlotSizeX, PlotSizeY: Size in pixels of the plot file.
//                         Root default (700x500?) is used if either is zero.
//   PlotShowFit: Flag indicating how fit should be displayed.
//                  >= 1: Show final fit
//                  >= 2: Show starting fit function
//   PlotSplitX: If this is nonzero, plots are created in updateMap (not update)
//               and the drawing canvas is split NY x NX where NX = PlotSplitX.
//               If PlotSplitX == 0, one canvas/plot is created in update.
//   PlotSplitY: If PlotSplitY > 0, then the above NY = PlotSplitY. Otherwise
//               NY = PlotSplitX. No effect if PlotSplitX == 0.
//               and up to that many plots are shown on the screen.
//
// Tools:
//   adcStringBuilder is used to make the following
// substitutions in the names and title:
//      %RUN%    --> run number
//      %SUBRUN% --> subrun number
//      %EVENT%  --> event number
//      %CHAN%   --> channel number
//
// The updating methods add metadata to the ADC channel data:
//
//   fitPedFractionLow     - fraction of samples below the fit range
//   fitPedFractionHigh    - fraction of samples above the fit range
//   fitPedestal           - fitted pedestal (same as the assigned value)
//   fitPedRms             - Fit sigma of the pedestal
//   fitPedChiSquare       - Chi-square of the fit
//   fitPedPeakBinFraction - Fraction of the pedestal distribution in the peak channel
//   fitPedPeakBinExcess   - Fraction above fit the peak channel
//   fitPedNBinsRemoved    - Number of sticky bins removed before fit
//
// The single-channel methods return a data map with the following:
//   pedestal           - pedestal histogram
//   fitFractionLow     - fraction of samples below the fit range
//   fitFractionHigh    - fraction of samples above the fit range
//   fitPedestal        - mean from the pedestal fit
//   fitPedestalRms     - sigma from the pedestal
//   fitChiSquare       - chi-square from the pedestal fit
//   fitPeakBinFraction - Fraction of the pedestal distribution in the peak channel
//   fitPeakBinExcess   - Fraction above fit the peak channel
//   fitChannel         - ADC channel number
//   fitNBinsRemoved    - # bins removed before the fit

// The AdcChannelTool methods all return a data map with the following:
//
//   nPedFitGood       - # channels with a successful pedestalfit
//   nPedFitFail       - # channels with a failed pedestal fit
//   fitStats          - fit status for each channel
//   fitPedestals      - fitted pecdestal for each channel
//   fitPedestalRmss   - fittend pedestal RMS (sigma) for each channel

#ifndef AdcPedestalFitter_H
#define AdcPedestalFitter_H

#include "art/Utilities/ToolMacros.h"
#include "fhiclcpp/ParameterSet.h"
#include "dune/DuneInterface/Tool/AdcChannelTool.h"
#include <string>
#include <vector>

class HistogramManager;
class AdcChannelStringTool;
class TH1;
class TPadManipulator;

class AdcPedestalFitter
: public AdcChannelTool {

public:

  using Index = unsigned int;

  AdcPedestalFitter(fhicl::ParameterSet const& ps);

  DataMap view(const AdcChannelData& acd) const override;

  DataMap update(AdcChannelData& acd) const override;

  DataMap updateMap(AdcChannelDataMap& acds) const override;

private:

  using Name = std::string;
  using NameVector = std::vector<Name>;

  // Configuration data.
  int m_LogLevel;
  float m_FitRmsMin;
  float m_FitRmsMax;
  Name m_HistName;
  Name m_HistTitle;
  Name m_HistManager;
  Name m_PlotFileName;
  Name m_RootFileName;
  Index m_PlotSizeX;
  Index m_PlotSizeY;
  Index m_PlotShowFit;
  Index m_PlotSplitX;
  Index m_PlotSplitY;

  // ADC string tool.
  const AdcChannelStringTool* m_adcStringBuilder;

  // Histogram manager.
  HistogramManager* m_phm;

  // Make replacements in a name.
  Name nameReplace(Name name, const AdcChannelData& acd, bool isTitle) const;

  // Find and return pedestal.
  DataMap getPedestal(const AdcChannelData& acd) const;

  // Fill the pad for a channel.
  // Histogram "pedestal" from dm is drawn.
  int fillChannelPad(DataMap& dm, TPadManipulator* pman) const;

};

DEFINE_ART_CLASS_TOOL(AdcPedestalFitter)

#endif
