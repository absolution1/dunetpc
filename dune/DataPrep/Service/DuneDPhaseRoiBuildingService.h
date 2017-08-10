// DuneDPhaseRoiBuildingService.h
//
// David Adams
// June 2016
//
// Implementation of service to build ROIs in AdcChannelData using
// the same algorithm as in the original DUNE 35t module. An unpadded
// ROI starts for any signal above NSigmaStart sigma above the noise
// and continues until the level falls below NSigmaEnd sigma. The
// noise level sigma = sss.GetDeconNoise() where sss is the DUNE
// signal shaping service. Note that NSigmaEnd was fixed to 1.0 in
// the original code.
//
// The ROIs are then padded to include PadLow channels below and PadHigh
// channels above the unpadded region. Overlapping ROIs are then merged.
//
// The original code defines the upper value for an unpadded ROI to be one
// tick beyond the last tick above threshold unless it extends to the
// end of the sample array. Here the upper value is always set to be the
// last tick above threshold. one should be able to reproduce the old
// behavior by adding one to NSigmaEnd.
//
// The ROIs are recorded in in data.rois.
//
// Configuration:
//   LogLevel    - usual log level
//   NSigmaStart - Level in sigma at which an unpadded signal starts.
//   NSigmaEnd   - Level in sigma at which an unpadded signal ends.
//   PadLow      - Number of ticks to retain before signal above threshold.
//   PadHigh     - Number of ticks to retain after signal above threshold.
//
// R.Sulej:
//  - copied from DuneRoiBuildingService.h to substitute signal shaping
//    service for double-phase;
//  - added option to look for ROI's in the filtered signal, so one can
//    remove slow fluctuations of pedestal; this is identical filter as
//    in the noise removal and one can easyli inspect the filter output
//    in the event display if it is applied at noise removal step.
//
#ifndef DuneDPhaseRoiBuildingService_H
#define DuneDPhaseRoiBuildingService_H

#include "dune/DuneInterface/AdcRoiBuildingService.h"

class AdcSuppressService;

class DuneDPhaseRoiBuildingService : public AdcRoiBuildingService {

public:

  DuneDPhaseRoiBuildingService(fhicl::ParameterSet const& pset, art::ActivityRegistry&);

  int build(AdcChannelData& data) const;

  std::ostream& print(std::ostream& out =std::cout, std::string prefix ="") const;

private:
  AdcSignalVector getLowFreqFiltered(const AdcSignalVector& adc) const;

  // Parameters.
  bool m_UseFilter;
  std::vector< float > m_FltCoeffs;

  int m_LogLevel;
  AdcSignal m_NSigmaStart;
  AdcSignal m_NSigmaEnd;
  AdcIndex m_PadLow;
  AdcIndex m_PadHigh;

};

DECLARE_ART_SERVICE_INTERFACE_IMPL(DuneDPhaseRoiBuildingService, AdcRoiBuildingService, LEGACY)

#endif
