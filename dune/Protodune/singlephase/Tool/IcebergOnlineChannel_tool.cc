// IcebergOnlineChannel.cxx

#include "IcebergOnlineChannel.h"
#include "dune-raw-data/Services/ChannelMap/IcebergChannelMapService.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include <iostream>

using std::string;
using std::cout;
using std::endl;
using Index = IcebergOnlineChannel::Index;

//**********************************************************************

IcebergOnlineChannel::IcebergOnlineChannel(const fhicl::ParameterSet& pset) 
: m_LogLevel(pset.get<Index>("LogLevel")),
  m_Ordering(pset.get<Name>("Ordering")),
  m_orderByWib(m_Ordering == "WIB"),
  m_orderByConnector(m_Ordering == "connector"),
  m_orderByFemb(m_Ordering == "FEMB") {
  const string myname = "IcebergOnlineChannel::ctor: ";
  cout << myname << "  LogLevel: " << m_LogLevel << endl;
  cout << myname << "  Ordering: " << m_Ordering << endl;
  if ( !m_orderByWib && !m_orderByConnector && !m_orderByFemb ) {
    cout << myname << "ERROR: Invalid ordering: " << m_Ordering << endl;
  }
}
  
//**********************************************************************

Index IcebergOnlineChannel::get(Index chanOff) const {
  const string myname = "IcebergOnlineChannel::get: ";
  if ( chanOff >= 1280 ) {
    if ( m_LogLevel > 1 ) cout << myname << "Invalid offline channel: " << chanOff << endl;
    return badIndex();
  }
  art::ServiceHandle<dune::IcebergChannelMapService> pms;
  // Fetch APA index.
  Index iapa = pms->APAFromOfflineChannel(chanOff);
  if ( iapa > 1 ) {
    if ( m_LogLevel > 1 ) cout << myname << "Invalid APA index: " << iapa << endl;
    return badIndex();
  }
  Index kapa = iapa;
  // Fetch WIB index.
  Index iwib = pms->WIBFromOfflineChannel(chanOff);
  if ( iwib > 4 ) {
    if ( m_LogLevel > 1 ) cout << myname << "Invalid WIB index: " << iwib << endl;
    return badIndex();
  }
  Index kwib = iwib;
  // Fetch connector index.
  Index icon = pms->FEMBFromOfflineChannel(chanOff);
  if ( icon < 1 || icon > 4 ) {
    if ( m_LogLevel > 1 ) cout << myname << "Invalid connector index: " << icon << endl;
    return badIndex();
  }
  Index kcon = icon - 1;
  // Fetch FEMB channel.
  Index ichf = pms->FEMBChannelFromOfflineChannel(chanOff);
  if ( ichf > 127 ) {
    if ( m_LogLevel > 1 ) cout << myname << "Invalid FEMB channel index: " << ichf << endl;
    return badIndex();
  }
  Index kchf = ichf;
  // Evauate online index.
  Index chanOn = 0;
  Index kfmb = 99;
  if ( m_orderByWib ) {
    cout << myname << "WIB ordering is not yest supported." << endl;
    return badIndex();
  } else if ( m_orderByConnector ) {
    cout << myname << "Connector ordering is not yest supported." << endl;
  } else if ( m_orderByFemb ) {
    // FEMB index mapping.
    Index ifmb = 99;
    if        ( chanOff <  200 ) {
      ifmb = ( 199-chanOff)/40 + 0;
    } else if ( chanOff <  400 ) {
      ifmb = ( 399-chanOff)/40 + 5;
    } else if ( chanOff <  800 ) {
      ifmb = (chanOff- 400)/40 + 0;
    } else if ( chanOff < 1040 ) {
      ifmb = (chanOff- 800)/48 + 0;
    } else if ( chanOff < 1280 ) {
      ifmb = (1279-chanOff)/48 + 5;
    } else {
      cout << "Invalid offline channel: : " << chanOff << endl;
      return badIndex();
    }
    kfmb = ifmb + 1;
    chanOn = 1280*kapa + 128*(ifmb) + kchf;
  }
  if ( m_LogLevel >= 4 ) {
    cout << myname << "ChanOff --> (kapa kwib kcon kchf) --> chanON: "
         << chanOff << " --> (" << kapa << " " << kwib << " " << kcon << " " << kchf
         << ") -->  " << chanOn
         << " --> " << 100*(kapa+1) + kfmb << "-" << kchf << endl;
  }
  return chanOn;
}

//**********************************************************************

DEFINE_ART_CLASS_TOOL(IcebergOnlineChannel)
