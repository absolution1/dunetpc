// PdspOnlineChannel.cxx

#include "PdspOnlineChannel.h"
#include "dune-raw-data/Services/ChannelMap/PdspChannelMapService.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include <iostream>

using std::string;
using std::cout;
using std::endl;
using Index = PdspOnlineChannel::Index;

//**********************************************************************

PdspOnlineChannel::PdspOnlineChannel(const fhicl::ParameterSet& pset) 
: m_LogLevel(pset.get<Index>("LogLevel")),
  m_Ordering(pset.get<Name>("Ordering")),
  m_orderByWib(m_Ordering == "WIB"),
  m_orderByConnector(m_Ordering == "connector"),
  m_orderByFemb(m_Ordering == "FEMB") {
  const string myname = "PdspOnlineChannel::ctor: ";
  cout << myname << "  LogLevel: " << m_LogLevel << endl;
  cout << myname << "  Ordering: " << m_Ordering << endl;
  if ( !m_orderByWib && !m_orderByConnector && !m_orderByFemb ) {
    cout << myname << "ERROR: Invalid ordering: " << m_Ordering << endl;
  }
}
  
//**********************************************************************

Index PdspOnlineChannel::get(Index chanOff) const {
  const string myname = "PdspOnlineChannel::get: ";
  if ( chanOff >= 15360 ) {
    if ( m_LogLevel > 1 ) cout << myname << "Invalid offline channel: " << chanOff << endl;
    return badIndex();
  }
  art::ServiceHandle<dune::PdspChannelMapService> pms;
  // Fetch APA index.
  Index iapa = pms->APAFromOfflineChannel(chanOff);
  if ( iapa > 5 ) {
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
  if ( m_orderByWib ) {
    chanOn = 2560*kapa + 512*kwib + 128*kcon + kchf;
  } else if ( m_orderByWib ) {
    chanOn = 2560*kapa + 640*kcon + 128*kwib + kchf;
  } else if ( m_orderByFemb ) {
    // FEMB index mapping.
    static Index ifmb[20] = {10,  9,  8,  7,  6,
                              5,  4,  3,  2,  1,
                             20, 19, 18, 17, 16,
                             15, 14, 13, 12, 11};
    Index jfmb = 5*kcon + kwib;
    Index kfmb = ifmb[jfmb] - 1;
    chanOn = 2560*kapa + 128*kfmb + kchf;
  }
  return chanOn;
}

//**********************************************************************
