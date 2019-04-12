////////////////////////////////////////////////////////////////////////
// ExpTailRemover.h
//
// Tool to remove an exponential tail from AC-coupled front-end electronics.
//
// The measured charge vector is assumed to be a sum of signal and tail:
//
//   Qdat = Qsig + Qtai
//
// and this tool replaces the sample vector (assumed to be the sum) in
// AdcChannelData with signal, i.e. it removes the tail.
//
// Any signal contribution is assumed to be accompanied by an exponential tail
// whose integral cancels that of the signal. The decay time of the exponential
// is a parameter of this tool.
//
// The data is fit with this model and two varied parameters: the (residual)
// pedestal ped and the initial charge of the tail csi.
//
// The fit is iterative using samples not flagged as signal to evaluate these
// parameters. Each iteration, the samples are updated and those values are used
// to evaluate the signal for the next iteration. Iterations stop when the signal
// evalaution does not change or the maximum # iterations is reached.
//
// The fit is done in time-domain, i.e there is no deconvolution following the
// the model in tool UndershootCorr and DUNE-doc-11662.
// Details of the fit procedure are in DUNE-doc-XXXXX.
//
// Configuration:
//   LogLevel - 0=silent, 1=init, 2=each event, >2=more
//   DecayTime: Exponential decay time in ticks.
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
//   CorrectFlag: switch per plane (induction, u, v) to turn the correction on and off
//                If empty, all planes corrected.
//
/////////////////////////////////////////////////////////////////////////
#ifndef ExpTailRemover_H
#define ExpTailRemover_H

#include "art/Utilities/ToolMacros.h"
#include "fhiclcpp/ParameterSet.h"
#include "dune/DuneInterface/Tool/AdcChannelTool.h"
#include <vector>

class ExpTailRemover : AdcChannelTool {

public:

  using Index = unsigned int;
  using Vector = std::vector<double>;
  using Name = std::string;

  ExpTailRemover(fhicl::ParameterSet const& ps);

  ~ExpTailRemover() override =default;

  DataMap update(AdcChannelData& acd) const override;

private:

  // Configuration data.
  int                m_LogLevel;
  Index              m_SignalFlag;
  Index              m_SignalIterationLimit;
  Name               m_SignalTool;
  double             m_DecayTime;
  std::vector<bool>  m_CorrectFlag;   // 0: U, 1: V, 2: Collection

  // Derived from configuration.
  AdcChannelTool* m_pSignalTool;

  // private methods

  // Remove the pedestal and tail from a data vactor.
  void getSignal(const Vector& qdats, double ped, double csi, Vector& qsigs) const;

  // Use the samples qdats to estimate the pedestal ped and initial tail charge csi.
  // nsamFit is the # samples used in the fit, i.e. those not considered signal
  void estimatepars(const Vector& qdats, double& ped, double& csi, Index& nsamFit) const;

  // Fit to line.
  void wlinfit(const Vector& x, const Vector& y, const Vector& e,
               double& slope, double& intercept) const;

};

DEFINE_ART_CLASS_TOOL(ExpTailRemover)

#endif
