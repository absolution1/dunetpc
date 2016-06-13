// Dune35tNoiseRemovalService.h
//
// David Adams
// June 2016
//
// Implementation of service that remove coherent noise from data
// produced with the DUNE 35-ton detector.
//
// This implementation is specific to the 35-ton detector because it uses
// the LBNE channel map service. This could be "fixed" by generalizing
// the channel map interface.
//
// The code is copied from dunetpc/dune/daqinput35t/FilterWF_module.cc but the
// following chnages are made.
//   1. Corrections are applied to floating signals after pedestal subtraction.
//      There might be change doe to int-to-float conversion.
//   2. The old code did not correct samples with stuck codes. Here it is an option.
//   3. The old code excluded samples with count < 10 or pedestal < 10 from the
//      correction. Such samples are retained here.
//
// Configuration:
//          LogLevel - message logging level: 0=none, 1=initialization, 2+=every event
//      GroupingFlag - 0=By regulator (128 chan), 1=By ASIC (32 channels)
//    SkipStuckCodes - Samples with stuck bits are not used to evaluate the correction.
// CorrectStuckCodes - Samples with stuck bits are corrected iff this is true.
//        ShowGroups - Display channel groups.
//

#ifndef Dune35tNoiseRemovalService_H
#define Dune35tNoiseRemovalService_H

#include "dune/DuneInterface/AdcNoiseRemovalService.h"
#include "dune/DuneInterface/AdcTypes.h"

namespace geo {
  class Geometry;
}

namespace lbne {
  class ChannelMapService;
}

class Dune35tNoiseRemovalService : public AdcNoiseRemovalService {

public:

  Dune35tNoiseRemovalService(fhicl::ParameterSet const& pset, art::ActivityRegistry&);

  int update(AdcChannelDataMap& datamap) const;

  std::ostream& print(std::ostream& out =std::cout, std::string prefix ="") const;

private:

  // Configuration parameters.
  int  m_LogLevel;
  bool m_GroupingFlag;
  bool m_SkipStuckCodes;
  bool m_CorrectStuckCodes;
  bool m_ShowGroups;

  // Services.
  const geo::Geometry* m_pGeometry;
  const lbne::ChannelMapService* m_pChannelMap;

  // List of channels for each orientation and group.
  std::vector<std::vector<AdcChannelVector>> m_GroupChannels;

};

DECLARE_ART_SERVICE_INTERFACE_IMPL(Dune35tNoiseRemovalService, AdcNoiseRemovalService, LEGACY)

#endif
