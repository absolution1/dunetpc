// AcdWireReader.h

// David Adams
// August 2017
//
// Tool to fill AdcChannelData from recob::Wire.
//
// Caller must fill in acd.wire before calling.
//
// Configuration:
//   LogLevel - 0=silent, 1=init, 2=each event, >2=more

#ifndef AcdWireReader_H
#define AcdWireReader_H

#include "art/Utilities/ToolMacros.h"
#include "fhiclcpp/ParameterSet.h"
#include "dune/DuneInterface/Tool/AdcChannelTool.h"

class AcdWireReader : AdcChannelTool {

public:

  AcdWireReader(fhicl::ParameterSet const& ps);

  ~AcdWireReader() override =default;

  DataMap update(AdcChannelData& acd) const override;

private:

  // Configuration data.
  int   m_LogLevel;

};

DEFINE_ART_CLASS_TOOL(AcdWireReader)

#endif
