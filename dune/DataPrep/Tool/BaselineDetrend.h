////////////////////////////////////////////////////////////////////////////////////
//
//  Remove trends in baseline fluctuations. This is basically a very low pass 
//  filter removing oscillations at ~1 kHz. The smoothed version of the pedestal
//  is contructed using LOWESS smoothing algorithm implemented in TSmoothGraph 
//  See docdb DUNE-doc-21111
//
//  Parametres:
//     LogLevel: controls level of printout
//
//     UseBasicROI: use basic ROI finder (para: Threshold and Pad ticks) 
//                  to bypass the signal regions
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
#include "dune/DuneInterface/Tool/TpcDataTool.h"

#include <vector>

class TGraphSmooth;
class TGraph;

//typedef std::vector<float> AdcSignalVector;

class BaselineDetrend : public TpcDataTool
{
  
 public:
  BaselineDetrend( fhicl::ParameterSet const& ps );
  ~BaselineDetrend() override;

  DataMap update(AdcChannelData& acd) const override;
  
 private:
  int      m_LogLevel;
  bool     m_UseBasicROI;
  float    m_Thresh;
  unsigned m_Pad;
  float    m_MinFrac;
  float    m_Span;

  //
  // smoother
  TGraphSmooth *m_GS;

  //
  AdcSignalVector Smoother( const AdcSignalVector &data,
			    const std::vector<unsigned> &pedidx ) const;
};


#endif
