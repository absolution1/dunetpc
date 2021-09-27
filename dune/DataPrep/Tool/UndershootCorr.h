////////////////////////////////////////////////////////////////////////
// UndershootCorr.h
//
// Tool to correct undershoot induced by AC-coupled front-end electronics
//  It implements a time-domain removal with a fit for the pedestal and charge
//  deposited before the event starts.
//
// Configuration:
//   LogLevel - 0=silent, 1=init, 2=each event, >2=more
//   CorrectFlag: switch per plane to turn the correction on and off
//   TDecayConst: constant used to decay accumulated charge sum on each tick.
//                csum *= TDecayConst on each tick and then
//                new charge is accumulated after pedestal subtraction
//   FSubConst: fraction of the accumulated charge to subtract from the next tick's data
//   RestoreBaseline: If true (for each plane), the median signal is added back after correction
//                    to restore the input baseline.
//   SignalThreshold: Ticks with signals greated than this are deweighted in the fit
//   SignalUnit: If non-blank, then this is compared with sampleUnit and a warning
//               broadcast if they differ.
//   LCA, LCB, LCC, LCD: Linear combination coefficients to translate from slope and intercept of a linear fit
//                 to the underhsoot-corrected waveform with zero initial charge and pedesal offset to get the
//                 fit initial charge and pedestal offset.
//
// See DUNE-doc-11662                
//
/////////////////////////////////////////////////////////////////////////
#ifndef UndershootCorr_H
#define UndershootCorr_H

#include "art/Utilities/ToolMacros.h"
#include "fhiclcpp/ParameterSet.h"
#include "dune/DuneInterface/Tool/TpcDataTool.h"
#include <vector>

class UndershootCorr : TpcDataTool {

public:

  UndershootCorr(fhicl::ParameterSet const& ps);

  ~UndershootCorr() override =default;

  DataMap update(AdcChannelData& acd) const override;

private:

  // Configuration data.
  int            m_LogLevel;
  std::vector<bool>           m_CorrectFlag;   // 0: U, 1: V, 2: Collection
  std::vector<double>         m_TDecayConst;
  std::vector<double>         m_FSubConst;
  std::vector<bool>           m_RestoreBaseline;
  std::vector<double>         m_SignalThreshold;
  std::string                 m_SignalUnit;
  std::vector<double>         m_LCA;
  std::vector<double>         m_LCB;
  std::vector<double>         m_LCC;
  std::vector<double>         m_LCD;

  // private methods

  void crc(std::vector<double> &orig, std::vector<double> &corr, double pedi, double csi, size_t plane) const;
  void estimatepars(std::vector<double> &x, std::vector<double> &y, double &pedi, double &csi, size_t plane) const;
  void wlinfit(std::vector<double> x, std::vector<double> &y, std::vector<double> &e, double &slope, double &intercept) const;
};


#endif
