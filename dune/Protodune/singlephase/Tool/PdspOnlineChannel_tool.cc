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

PdspOnlineChannel::PdspOnlineChannel(const fhicl::ParameterSet&) 
: m_LogLevel(2) { }
  
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
  Index chanOn = 2560*kapa + 512*kwib + 128*kcon + kchf;
  return chanOn;
}

//**********************************************************************
