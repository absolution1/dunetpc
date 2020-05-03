////////////////////////////////////////////////////////////////////////////////////
//
//  Remove trends in baseline fluctuations. This is basically a very low pass 
//  filter removing oscillations at ~1 kHz. The smoothed version of the pedestal
//  is contructed using LOWESS smoothing algorithm implemented in TSmoothGraph 
//  Parametres:
//     LogLevel: controls level of printout
//
//     Threshold: to bypass signals during smoothing
//
//     Pad ticks: number of ticks around signal candidates
//
//     MinFrac:   min number of pedestal samples to build the smoothed representation
//
//     Span:      span parameters between 0 ... 1. The larger it is the smoother the 
//                result. The value that works the best (follows baseline variations 
//                at ~1kHz rate) is around 0.04.
//
//
//  Written: V. Galymov
//           April 2020
//
////////////////////////////////////////////////////////////////////////////////////

#ifndef __BASELINEDETREND_H__
#define __BASELINEDETREND_H__

// dune interfaces
#include "art/Utilities/ToolMacros.h"
#include "fhiclcpp/ParameterSet.h"
#include "dune/DuneInterface/Tool/AdcChannelTool.h"

#include <vector>

class TGraphSmooth;
class TGraph;

//typedef std::vector<float> AdcSignalVector;

class BaselineDetrend : public AdcChannelTool
{
  
 public:
  BaselineDetrend( fhicl::ParameterSet const& ps );
  ~BaselineDetrend() override;

  DataMap update(AdcChannelData& acd) const override;
  
 private:
  int      m_LogLevel;
  float    m_Thresh;
  unsigned m_Pad;
  float    m_MinFrac;
  float    m_Span;

  //
  //AdcSignalVector m_Trend;
  
  // smoother
  TGraphSmooth *m_GS;

  //
  AdcSignalVector Smoother( const AdcSignalVector &data ) const;
};

DEFINE_ART_CLASS_TOOL(BaselineDetrend)

#endif
