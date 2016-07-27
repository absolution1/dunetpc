// InterpolatingAdcMitigationService.h
//
// David Adams
// June 2016
//
// Implementation of service that replaces flagged samples in AdcChannelData
// with values interpolated from neighboring unflagged samples.
// This is a way to deal with the stuck bits in DUNE 35-ton data.
//
// Up to MaxConsecutiveSamples consecutive samples are replaced with interpolated values.
//
// Code is very similar to that in CalData/UnstickADCCodes_module.
//
// Consecutive flagged samples at the beginning or end of the data are treated in the
// same way as thos in regions with more than MaxConsecutiveSamples samples.
//
// Configuration:
//   LogLevel - message logging level: 0=none, 1=initialization, 2+=every event
//   SkipUnderflows - flag indicating not to update underflows
//   SkipOverflows - flag indicating not to update overflows
//   MaxConsecutiveSamples - maxumum number of samples to interpolate over (negative for no limit)
//   MaxConsecutiveFlag - Flag specifying what to do when there are too many consecutive samples.
//                        0 - Make no change.
//                        1 - Replace values with 0.0.

#ifndef InterpolatingAdcMitigationService_H
#define InterpolatingAdcMitigationService_H

#include "dune/DuneInterface/AdcMitigationService.h"

namespace lariov {
  class DetPedestalProvider;
}

class InterpolatingAdcMitigationService : public AdcMitigationService {

public:

  InterpolatingAdcMitigationService(fhicl::ParameterSet const& pset, art::ActivityRegistry&);

  int update(AdcChannelData& data) const;

  std::ostream& print(std::ostream& out =std::cout, std::string prefix ="") const;

private:

  // Configuration parameters.
  int  m_LogLevel;
  bool m_SkipUnderflows;
  bool m_SkipOverflows;
  int  m_MaxConsecutiveSamples;
  int  m_MaxConsecutiveFlag;

};

DECLARE_ART_SERVICE_INTERFACE_IMPL(InterpolatingAdcMitigationService, AdcMitigationService, LEGACY)

#endif
