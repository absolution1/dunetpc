////////////////////////////////////////////////////////////////////////
// UndershootCorr.h
//
// Tool to preform subtract baseline using linear interpolation between 
// regions defined by the datasize and fBaseSampleBins
//
// Configuration:
//   LogLevel - 0=silent, 1=init, 2=each event, >2=more
//   BaseSampleBins - 
//   BaseVarCut -
//
/////////////////////////////////////////////////////////////////////////
#ifndef UndershootCorr_H
#define UndershootCorr_H

#include "art/Utilities/ToolMacros.h"
#include "fhiclcpp/ParameterSet.h"
#include "dune/DuneInterface/Tool/AdcChannelTool.h"
#include <vector>

class UndershootCorr : AdcChannelTool {

public:

  UndershootCorr(fhicl::ParameterSet const& ps);

  ~UndershootCorr() override =default;

  DataMap update(AdcChannelData& acd) const override;

private:

  // Configuration data.
  int            m_LogLevel;
  std::vector<bool>          m_CorrectFlag;   // 0: U, 1: V, 2: Collection
  std::vector<double>         m_TDecayConst;
  std::vector<double>         m_FSubConst;
  std::vector<double>         m_LCA;
  std::vector<double>         m_LCB;
  std::vector<double>         m_LCC;
  std::vector<double>         m_LCD;

  // private methods

  void crc(std::vector<double> &orig, std::vector<double> &corr, double pedi, double csi, size_t plane) const;
  void estimatepars(std::vector<double> &x, std::vector<double> &y, double &pedi, double &csi, size_t plane) const;
  void wlinfit(std::vector<double> x, std::vector<double> &y, std::vector<double> &e, double &slope, double &intercept) const;
};

DEFINE_ART_CLASS_TOOL(UndershootCorr)

#endif
