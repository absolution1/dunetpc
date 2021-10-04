// ChannelStatusServiceTool_tool.cc

#include "ChannelStatusServiceTool.h"
#include "dune/DuneInterface/Data/AdcTypes.h"
#include "fhiclcpp/intermediate_table.h"
#include "fhiclcpp/make_ParameterSet.h"
#include "larevt/CalibrationDBI/Interface/ChannelStatusService.h"
#include "larevt/CalibrationDBI/Interface/ChannelStatusProvider.h"

#include <iostream>
#include <iomanip>

using std::cout;
using std::endl;
using std::setw;
using Name = std::string;
using Index = ChannelStatusServiceTool::Index;

//**********************************************************************

ChannelStatusServiceTool::ChannelStatusServiceTool(fhicl::ParameterSet const& ps)
: m_LogLevel(ps.get<Index>("LogLevel")) {
  const Name myname = "ChannelStatusServiceTool::ctor: ";
  m_pChannelStatusProvider = &art::ServiceHandle<lariov::ChannelStatusService>()->GetProvider();
  if ( m_pChannelStatusProvider == nullptr ) {
    cout << myname << "WARNING: Channel status provider not found." << endl;
  }
  if ( m_LogLevel >= 1 ) {
    cout << myname << "     LogLevel: " << m_LogLevel << endl;
  }
}

//**********************************************************************

Index ChannelStatusServiceTool::get(Index icha) const {
  const Name myname = "ChannelStatusServiceTool::get: ";
  Index chanStat = AdcChannelStatusUnknown;
  if ( m_pChannelStatusProvider != nullptr ) {
    if      ( m_pChannelStatusProvider->IsBad(icha)   ) chanStat = AdcChannelStatusBad;
    else if ( m_pChannelStatusProvider->IsNoisy(icha) ) chanStat = AdcChannelStatusNoisy;
    else                                                chanStat = AdcChannelStatusGood;
  }
  if ( m_LogLevel >= 2 ) {
    if ( m_LogLevel >=3 || chanStat ) {
      cout << myname << "Status for channel " << icha << ": " << chanStat << endl;
    }
  }
  return chanStat;
}

//**********************************************************************

DEFINE_ART_CLASS_TOOL(ChannelStatusServiceTool)
