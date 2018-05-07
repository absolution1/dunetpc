// AcdDigitReader.h

// David Adams
// August 2017
//
// Tool to fill AdcChannelData from raw::RawDigit.
//
// Configuration:
//   LogLevel - 0=silent, 1=init, 2=each event, >2=more

#ifndef AcdDigitReader_H
#define AcdDigitReader_H

#include "art/Utilities/ToolMacros.h"
#include "fhiclcpp/ParameterSet.h"
#include "dune/DuneInterface/Tool/AdcChannelTool.h"

class AcdDigitReader : AdcChannelTool {

public:

  AcdDigitReader(fhicl::ParameterSet const& ps);

  ~AcdDigitReader() override =default;

  DataMap update(AdcChannelData& acd) const override;

private:

  // Configuration data.
  int   m_LogLevel;

};

DEFINE_ART_CLASS_TOOL(AcdDigitReader)

#endif
