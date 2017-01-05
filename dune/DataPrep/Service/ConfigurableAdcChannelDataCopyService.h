// ConfigurableAdcChannelDataCopyService.h
//
// David Adams
// September 2016
//
// Simple implementation of a service that copies AdcChannelData.
// There is a flag for each field.
//
// Configuration:
//   LogLevel       - message logging level: 0=none, 1=initialization, 2+=every copy
//   CopyChannel    - copy the channel number [false]
//   CopyPedestal   - copy the pedestal [false]
//   CopyRaw        - copy the raw data [false]
//   CopySamples    - copy the samples [false]
//   CopyFlags      - copy the flags [false]
//   CopyRois       - copy the rois [false]
//   CopyDigit      - copy the digit pointer [false]
//   CopyWire       - copy the wire pointer [false]
//   CopyDigitIndex - copy the digit index [false]
//   CopyWireIndex  - copy the wire index [false]

#ifndef ConfigurableAdcChannelDataCopyService_H
#define ConfigurableAdcChannelDataCopyService_H

#include "dune/DuneInterface/AdcChannelDataCopyService.h"

class ConfigurableAdcChannelDataCopyService : public AdcChannelDataCopyService {

public:

  ConfigurableAdcChannelDataCopyService(fhicl::ParameterSet const& pset, art::ActivityRegistry&);

  int copy(const AdcChannelData& acdin, AdcChannelData& acdout) const;

  std::ostream& print(std::ostream& out =std::cout, std::string prefix ="") const;

private:

  // Configuration parameters.
  int  m_LogLevel       =1;
  bool m_CopyChannel    =false;
  bool m_CopyPedestal   =false;
  bool m_CopyRaw        =false;
  bool m_CopySamples    =false;
  bool m_CopyFlags      =false;
  bool m_CopySignal     =false;
  bool m_CopyRois       =false;
  bool m_CopyDigit      =false;
  bool m_CopyWire       =false;
  bool m_CopyDigitIndex =false;
  bool m_CopyWireIndex  =false;

};

DECLARE_ART_SERVICE_INTERFACE_IMPL(ConfigurableAdcChannelDataCopyService, AdcChannelDataCopyService, LEGACY)

#endif
