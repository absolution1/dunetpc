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
// The fit is done in time-domain, i.e ther is no deconvolution following the
// the model in tool UndershootCorr and DUNE-doc-11662                .
//
// Configuration:
//   LogLevel - 0=silent, 1=init, 2=each event, >2=more
//   DecayTime: Exponential decay time in ticks.
//   SignalThreshold: Ticks with signals greated than this are deweighted in the fit
//   SignalUnit: If non-blank, then this is compared with sampleUnit and a warning
//               broadcast if they differ.
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

  ExpTailRemover(fhicl::ParameterSet const& ps);

  ~ExpTailRemover() override =default;

  DataMap update(AdcChannelData& acd) const override;

private:

  // Configuration data.
  int                m_LogLevel;
  std::string        m_SignalUnit;
  double             m_DecayTime;
  double             m_SignalThreshold;
  std::vector<bool>  m_CorrectFlag;   // 0: U, 1: V, 2: Collection

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
