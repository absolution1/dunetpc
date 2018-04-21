// FembMappedAdcModifier.h

// David Adams
// November, 2017
//
// AdcChannelTool that calls another modifier based on FEMB ID.
// The name of the called modifier is ToolBaseFFF where FFF is the FembID.
// If the FEMB ID is not set, the tool ToolBaseDefault is used.
//
// If the string DirName is not empty, then tool definition is taken from
// the fcl file DirName/ToolBaseFFF.fcl found on the usual fcl path.
// Otherwise, the tool definition must already exist.
//
// Parameters:
//   LogLevel - 0=silent, 1=init, 2=errors each event, 3=call each event
//   ToolBase - Prefix of the tool name.
//   DirName - Location for the tool fcl file.
//
// Output (appended to that from the called tool):
//   fembID - the FEMB ID

#ifndef FembMappedAdcModifier_H
#define FembMappedAdcModifier_H

#include "art/Utilities/ToolMacros.h"
#include "fhiclcpp/ParameterSet.h"
#include "dune/DuneInterface/Tool/AdcChannelTool.h"

class FembMappedAdcModifier
: public AdcChannelTool {

public:

  FembMappedAdcModifier(fhicl::ParameterSet const& ps);

  DataMap view(const AdcChannelData& acd) const override;

  DataMap update(AdcChannelData& acd) const override;

private:

  using Name = std::string;

  // Parameters.
  int m_LogLevel;
  Name m_ToolBase;
  Name m_DirName;

};

DEFINE_ART_CLASS_TOOL(FembMappedAdcModifier)

#endif
