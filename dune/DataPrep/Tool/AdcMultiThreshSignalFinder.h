// AdcMultiThreshSignalFinder.h
//
// Tool to find signals and build ROI's with respect to pedestal
// subtracted baseline. This is a more basic version of the algorithm 
// developed by Christoph Alt for the 3x1x1 dual-phase TPC detector
//
//
// Configuration:
//   LogLevel     - 0=silent, 1=init, 2=each event, >2=more
//   UseStandardDevidation - true/false to define thresholds wrt pedestal RMS
//   Threshold1   - 1st threshold to seed ROI (in ADC or units of ped RMS)
//   SamplesAboveThreshold1 - min number of continuous samples to be above 1st threshold
//   Threshold2   - 2nd threshold to refine ROI search (in ADC or units of ped RMS)
//   SamplesAboveThreshold2 - min number of continuous samples to be above 2nd threshold
//   ThresholdMax - threshold for which max value in seeded ROI
//   ThresholdMin - minimum value to which expand found ROIs
//   BinsBefore   - number of time bins to pad found ROI on the left
//   BinsAfter    - number of time bins to pad found ROI on the right
//


#ifndef AdcMultiThreshSignalFinder_H
#define AdcMultiThreshSignalFinder_H

#include "art/Utilities/ToolMacros.h"
#include "fhiclcpp/ParameterSet.h"
#include "dune/DuneInterface/Tool/TpcDataTool.h"
#include "dune/DuneInterface/Data/AdcTypes.h"

class AdcMultiThreshSignalFinder : TpcDataTool {

public:

  AdcMultiThreshSignalFinder(fhicl::ParameterSet const& ps);

  ~AdcMultiThreshSignalFinder() override;

  // Build ROIs
  DataMap update(AdcChannelData& acd) const override;

private:

  // Parameters.
  int       m_LogLevel;
  bool      m_UseStd;
  AdcSignal m_Thresh1;
  AdcIndex  m_NsaAbove1;
  AdcSignal m_Thresh2;
  AdcIndex  m_NsaAbove2;
  AdcSignal m_ThreshMax;
  AdcSignal m_ThreshMin;
  AdcIndex  m_BinsBefore;
  AdcIndex  m_BinsAfter;

  //
  typedef struct 
  {
    bool      isRoi;
    AdcIndex  StartRoi;
    AdcIndex  EndRoi;
    AdcIndex  NsaTh1;
    AdcIndex  NsaTh2;
    AdcIndex  NsaTh3;
    AdcSignal MaxValue;
    
    void init()
    {
      isRoi    = false;
      StartRoi = 0;
      EndRoi   = 0;
      NsaTh1   = 0;
      NsaTh2   = 0;
      NsaTh3   = 0;
      MaxValue = 0;
    }

  } ROICandidate_t;
};


#endif
