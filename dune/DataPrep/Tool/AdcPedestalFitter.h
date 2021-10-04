// AdcPedestalFitter.h

// David Adams
// August 2017
//
// Tool to fit ADC data and extract a pedestal. The pedestal is the mean of a Gaussian fit to
// the ADC distributions.
//
// If FitRmsMin < FitRmsMax, the the RMS is constrained to the range
// (FitRmsMin, FitRmsMax) in the fit.
//
// The starting mean and center of the histogram is the peak bin. If RemoveStickyCode is true, then
// one apparent sticky code may be removed when evaluating this peak.
//
// Configuration:
//   LogLevel - 0=silent, 1=init, 2=each event, >2=more
//   AdcRange - ADC values must less than this value, e.g. 4096 for 12 bits.
//   FitOpt - 0: Use histogram mean (no fit)
//            1: Chi-square fit
//            2: Likelihood fit
//            3: Chi-square fit fillowed by likelihood fit if chi-square fit fails
//   FitPrecision - Only used for chi-square fit? (TFitter::SetPrecision)
//   SkipFlags - Samples with these flags are excluded
//   AdcFitRange: Width [ADC counts] of the fit range.
//   FitRmsMin: Lower limit for RMS fit range.
//   FitRmsMax: Upper limit for RMS fit range.
//   RemoveStickyCode: If true, an apparent sticky code may be removed.
//   HistName:  Name for the histogram.
//   HistTitle: Title for the histogram.
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
//   fitPedFractionLow      - fraction of samples below the fit range
//   fitPedFractionHigh     - fraction of samples above the fit range
//   fitPedestal            - fitted pedestal (same as the assigned value)
//   fitPedRms              - Fit sigma of the pedestal
//   fitPedChiSquare        - Chi-square of the fit
//   fitPedReducedChiSquare - An estimate of the chi-square/DOF using only the central
//                            (3-sigma) part of the pedestal distribution
//   fitPedPeakBinFraction  - Fraction of the pedestal distribution in the peak channel
//   fitPedPeakBinExcess    - Fraction above fit the peak channel
//   fitPedNBinsRemoved     - Number of sticky bins removed before fit
//
// The single-channel methods return a data map with the following:
//   pedestal            - pedestal histogram
//   fitFractionLow      - fraction of samples below the fit range
//   fitFractionHigh     - fraction of samples above the fit range
//   fitPedestal         - mean from the pedestal fit
//   fitPedestalRms      - sigma from the pedestal
//   fitChiSquare        - chi-square from the pedestal fit
//   fitReducedChiSquare  - extimate of central reduced chi-square
//   fitPeakBinFraction  - Fraction of the pedestal distribution in the peak channel
//   fitPeakBinExcess    - Fraction above fit the peak channel
//   fitChannel          - ADC channel number
//   fitNBinsRemoved     - # bins removed before the fit

// The TpcDataTool methods all return a data map with the following:
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
#include "dune/DuneInterface/Tool/TpcDataTool.h"
#include <string>
#include <vector>
#include <set>

class ParFormula;
class RunDataTool;
class AdcChannelStringTool;
class TH1;
class TF1;
class TPadManipulator;

class AdcPedestalFitter
: public TpcDataTool {

public:

  using Index = unsigned int;
  using IndexVector = std::vector<Index>;
  using IndexSet = std::set<Index>;

  AdcPedestalFitter(fhicl::ParameterSet const& ps);

  DataMap view(const AdcChannelData& acd) const override;

  DataMap update(AdcChannelData& acd) const override;

  DataMap updateMap(AdcChannelDataMap& acds) const override;

  DataMap beginEvent(const DuneEventInfo&) const override;

private:

  using Name = std::string;
  using NameVector = std::vector<Name>;
  using FormulaMap = std::map<Name, ParFormula*>;
  using NameSet = std::set<Name>;
  using NameSetMap = std::map<Name, NameSet>;

  // Configuration data.
  int m_LogLevel;
  Name m_AdcRange;
  Index m_FitOpt;
  float m_FitPrecision;
  IndexVector m_SkipFlags;
  Name m_AdcFitRange;
  Name m_FitRmsMin;
  Name m_FitRmsMax;
  bool m_RemoveStickyCode;
  Name m_HistName;
  Name m_HistTitle;
  Name m_PlotFileName;
  Name m_RootFileName;
  Index m_PlotSizeX;
  Index m_PlotSizeY;
  Index m_PlotShowFit;
  Index m_PlotSplitX;
  Index m_PlotSplitY;

  // ADC string tool.
  const AdcChannelStringTool* m_adcStringBuilder;

  // Derived from config.
  IndexSet m_skipFlags;
  NameVector m_fitOpts;
  FormulaMap m_tfs;         // Formulas for AdcRange, AdcFitRange, FitRmsMin, FitRmsMax
  NameSetMap m_tfpars;       // Formula parameters.
  bool m_haveFormulaParams;  // Do any formulas have parameters?
  const RunDataTool* m_prdtool; // Run data tool.

  // State.
  class State {
  public:
    Index adcRange =0;
    TF1* pfitter = nullptr;
    Index ncall = 0;
    Index npeakBinSuppressed = 0;
    Index run = -1;
    Index nevt = 0;
  };
  std::unique_ptr<State> m_pstate;
  State& state() const { return *m_pstate; }

  // Make replacements in a name.
  Name nameReplace(Name name, const AdcChannelData& acd, bool isTitle) const;

  // Find and return pedestal.
  DataMap getPedestal(const AdcChannelData& acd) const;

  // Fill the pad for a channel.
  // Histogram "pedestal" from dm is drawn.
  int fillChannelPad(DataMap& dm, const AdcChannelData& acd, TPadManipulator* pman) const;

};


#endif
