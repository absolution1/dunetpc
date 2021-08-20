// RemoveBadChannels_tool.cc

#include "RemoveBadChannels.h"
#include "larevt/CalibrationDBI/Interface/ChannelStatusService.h"
#include "larevt/CalibrationDBI/Interface/ChannelStatusProvider.h"
#include <iostream>

using std::string;
using std::cout;
using std::endl;

using Index = unsigned int;

//**********************************************************************
// Class methods.
//**********************************************************************

RemoveBadChannels::RemoveBadChannels(fhicl::ParameterSet const& ps)
: m_LogLevel(ps.get<int>("LogLevel"))
, m_RemoveBadChs(ps.get<bool>("RemoveBadChs"))
, m_RemoveNoisyChs(ps.get<bool>("RemoveNoisyChs")) {
  const string myname = "RemoveBadChannels::ctor: ";
  if ( m_LogLevel >= 1 ) {
    cout << myname << "Parameters:" << endl;
    cout << myname << "         LogLevel: " << m_LogLevel << endl;
    cout << myname << "         RemoveBadChs: " << m_RemoveBadChs << endl;
    cout << myname << "         RemoveNoisyChs: " << m_RemoveNoisyChs << endl;
  }
  m_pChannelStatusProvider = &art::ServiceHandle<lariov::ChannelStatusService>()->GetProvider();
  if ( m_pChannelStatusProvider == nullptr ) {
    cout << myname << "WARNING: Channel status provider not found." << endl;
  }
}

//**********************************************************************

DataMap RemoveBadChannels::update(AdcChannelData& acd) const {
  const string myname = "RemoveBadChannels::update: ";
  if ( m_LogLevel >= 2 ) cout << "Processing run " << acd.run() << " event " << acd.event()
                              << " channel " << acd.channel() << endl;
  DataMap ret;

  //Set channel output to 0 for bad channels
  if ((m_pChannelStatusProvider->IsBad(acd.channel()) && m_RemoveBadChs) ||
       (m_pChannelStatusProvider->IsNoisy(acd.channel()) && m_RemoveNoisyChs)){
    if ( m_LogLevel >= 2 ) cout << "Zeroing " << acd.run() << " event " << acd.event() << " channel " << acd.channel() << endl;
    for ( float& sam : acd.samples ) sam = 0;
  }

  return ret;
}

//**********************************************************************

DEFINE_ART_CLASS_TOOL(RemoveBadChannels)
