// ExpTailPedRemover.h
//
// David Adams
// September 2020
//
// Tool to remove an exponential tail and tim-dependent pedestal from
// AC-coupled front-end electronics. See DUNE-doc-20618.
//
// The measured charge vector is assumed to be a sum of signal and tail:
//
//   Qdat = Qsig + Qtai + Qped + Qnoise
//
// where <Qnoise> = 0.
//
// This tool replaces the sample vector (assumed to be Qdat) in
// AdcChannelData with signal (pus noise), i.e. it removes the tail and pedestal.
//
// Any signal contribution is assumed to be accompanied by an exponential tail
// whose integral cancels that of the signal. The decay time of the exponential
// is a parameter of this tool. If 0, there is no tail.
//
// The data is fit with this model and Np + 1 varied parameters: the initial
// charge Qtail[0] and the Np parameters in the pedestal model and the initial charge of the tail csi.
//
// The fit is iterative using samples not flagged as signal to evaluate these
// parameters. Each iteration, the samples are updated and those values are used
// to evaluate the signal for the next iteration. Iterations stop when the signal
// evalaution does not change or the maximum # iterations is reached.
//
// The pedestal model is a polynomial in (tick-Pedtick0) with degree PedDgree
// plus a series of fixed frequncy terms for the frequncies listed in PedFreqs.
//
// Configuration:
//   LogLevel - 0=silent, 1=init, 3=each event, >4=more
//   DecayTime: Exponential decay time in ticks. If <= 0, the tail is not used.
//   MaxTick: Maximum # ticks expected. Pedestal arrays are this size.
//   PedDegree: Degree of the pedestal polynomial: <0=none, 1=flat, 2=linear, >=3=quadratic
//   PedTick0: polynomial = SUM lambda_i (isam - PedTick0)^i
//   PedFreqs: Frequency terms (sin and cos) added for each of these frequencies [1/tick]
//   SignalFlag: Specifies how signal regions are identified in the samples.
//               0 - Use all samples.
//               1 - Use the input acd.signal. If not enough values, pad with
//                   false and log warning.
//               2 - Apply the signal finder on updated samples each iteration.
//               3 - Same as 2 plus the signal finder is applied on the final
//                   samples.
//   SignalIterationLimit - Maximum nimber of fit iterations.
//   SignalTool - Name of the tool used identify signals. The method called is
//                AdcChannelTool::update, i.e. single-channel signal finding.
//   NoWarnStatuses - Status flags for which warnings should not be logged.
//   IncludeChannelRanges - List of channel ranges to correct. If empty, all channels
//                          Use empty or "all" for all channels or "none" for no channels.
//                          Tool channelRanges is used to map these names to channels.
//   ExcludeChannelRanges - List of channel ranges whose channels are excluded from
//                          the list to process derived from IncludeChannelRanges.
//                          If "all", no channels are processed.
//
// The following metadata is added to each channel (i.e. to AdcChannelData):
//   uscPedestal - fitted pedestal (w.r.t. input samples)
//   uscTail     - Fitted tail for the first bin (tau[0])
//   uscNoise    - Noise = RMS of (sample - pedestal) in the no-signal region
//
// The following are returned in the call toi view:
//   uscPedestal   - As above
//       uscTail   - As above.
//    uscNsamFit   - # non-signal samples
//   uscNiteration - # fit iterations

#ifndef ExpTailPedRemover_H
#define ExpTailPedRemover_H

#include "art/Utilities/ToolMacros.h"
#include "fhiclcpp/ParameterSet.h"
#include "dune/DuneInterface/Tool/TpcDataTool.h"
#include <vector>
#include <set>

class ExpTailPedRemover : TpcDataTool {

public:

  using Index = unsigned int;
  using IndexVector = std::vector<Index>;
  using IndexSet = std::set<Index>;
  using Vector = std::vector<double>;
  using Name = std::string;
  using NameVector = std::vector<Name>;
  using FloatVector = std::vector<float>;
  using FloatVectorVector = std::vector<FloatVector>;

  ExpTailPedRemover(fhicl::ParameterSet const& ps);

  ~ExpTailPedRemover() override =default;

  DataMap update(AdcChannelData& acd) const override;

private:

  // Configuration data.
  int                m_LogLevel;
  Index              m_SignalFlag;
  Index              m_SignalIterationLimit;
  Name               m_SignalTool;
  double             m_DecayTime;
  Index              m_MaxTick;
  int                m_PedDegree;
  Index              m_PedTick0;
  FloatVector        m_PedFreqs;
  IndexVector        m_NoWarnStatuses;
  NameVector         m_IncludeChannelRanges;
  NameVector         m_ExcludeChannelRanges;

  // Derived from configuration.
  bool m_useChannelRanges;  // If true, only channels in checkChannels are processed.
  IndexSet m_nowarnStatuses;
  std::vector<bool> m_checkChannels;
  AdcChannelTool* m_pSignalTool;
  FloatVectorVector m_pedVectors;
  NameVector m_fitNames;  // Names for the fitted params: {Tail, Pedestal, Slope, Curvature, Cos, Sin, ...)
  Name m_fitNamesString;

  // private methods

  // Return if tail is used.
  bool useTail() const { return m_DecayTime > 0.0; }
  bool usePolynomial() const { return m_PedDegree >= 0; }

  // Remove the pedestal and tail from a data vactor.
  void getSignal(const Vector& qdats, double ped, double csi, Vector& qsigs) const;

  // Use the samples qdats to estimate the pedestal ped and initial tail charge csi.
  // nsamFit is the # samples used in the fit, i.e. those not considered signal
  void estimatepars(const Vector& qdats, double& ped, double& csi, Index& nsamFit) const;

  // Fit to line.
  void wlinfit(const Vector& x, const Vector& y, const Vector& e,
               double& slope, double& intercept) const;

};


#endif
