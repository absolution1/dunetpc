// PdspOnlineChannel.cxx

#include "PdspOnlineChannel.h"
#include "dune-raw-data/Services/ChannelMap/PdspChannelMapService.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"

using Index = PdspOnlineChannel::Index;

//**********************************************************************

PdspOnlineChannel::PdspOnlineChannel(const fhicl::ParameterSet&) { }
  
//**********************************************************************

Index PdspOnlineChannel::get(Index chanOff) const {
  if ( chanOff >= 15360 ) return badIndex();
  art::ServiceHandle<dune::PdspChannelMapService> pms;
  Index chanOn = 2560*pms->APAFromOfflineChannel(chanOff);
  chanOn += 512*pms->WIBFromOfflineChannel(chanOff);
  chanOn += 128*pms->FEMBFromOfflineChannel(chanOff);
  chanOn += pms->FEMBChannelFromOfflineChannel(chanOff);
  return chanOn;
}

//**********************************************************************
