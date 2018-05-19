// VintageDeconvoluter.h
//
// David Adams
// May 2018
//
// Tool to preform deconvolution of data from an ADC.
//
// It uses the service with interface util::SignalShapingServiceDUNE.
//
// Configuration:
//   LogLevel - 0=silent, 1=init, 2=each event, >2=more

#ifndef VintageDeconvoluter_H
#define VintageDeconvoluter_H

#include "art/Utilities/ToolMacros.h"
#include "fhiclcpp/ParameterSet.h"
#include "dune/DuneInterface/Tool/AdcChannelTool.h"
#include <vector>

class VintageDeconvoluter : AdcChannelTool {

public:

  VintageDeconvoluter(fhicl::ParameterSet const& ps);

  ~VintageDeconvoluter() override =default;

  DataMap update(AdcChannelData& acd) const override;

private:

  // Configuration data.
  int            m_LogLevel;

};

DEFINE_ART_CLASS_TOOL(VintageDeconvoluter)

#endif
