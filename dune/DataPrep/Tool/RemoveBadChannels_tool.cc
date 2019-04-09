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
: m_LogLevel(ps.get<int>("LogLevel")) {
  const string myname = "RemoveBadChannels::ctor: ";
  if ( m_LogLevel >= 1 ) {
    cout << myname << "Parameters:" << endl;
    cout << myname << "         LogLevel: " << m_LogLevel << endl;
  }
  m_pChannelStatusProvider = &art::ServiceHandle<lariov::ChannelStatusService>()->GetProvider();
  if ( m_pChannelStatusProvider == nullptr ) {
    cout << myname << "WARNING: Channel status provider not found." << endl;
  }
}

//**********************************************************************

DataMap RemoveBadChannels::update(AdcChannelData& acd) const {
  const string myname = "RemoveBadChannels::update: ";
  if ( m_LogLevel >= 2 ) cout << "Processing run " << acd.run << " event " << acd.event
                              << " channel " << acd.channel << endl;
  DataMap ret;

  //Set channel output to 0 for bad channels
  if (m_pChannelStatusProvider->IsBad(acd.channel) ||
      m_pChannelStatusProvider->IsNoisy(acd.channel)){
    for ( float& sam : acd.samples ) sam = 0;
  }

  return ret;
}

//**********************************************************************
