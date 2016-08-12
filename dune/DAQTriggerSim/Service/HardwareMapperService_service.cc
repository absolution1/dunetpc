// HardwareMapperService_service.cc
//
// Jonathan Davies j.p.davies@sussex.ac.uk
// August 2016
//
// Description

#include "dune/DAQTriggerSim/Service/HardwareMapperService.h"
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Principal/Event.h"


#include <set>

//......................................................
HardwareMapperService::HardwareMapperService(fhicl::ParameterSet const& pset, art::ActivityRegistry& reg){
  INFO_FILE_FUNCTION << std::endl;
}

//......................................................
int HardwareMapperService::getHardwareElement(int const& element, std::vector<int>& channelVec){
  int returnValue = -1;
  INFO_FILE_FUNCTION << std::endl;
  return returnValue;
}

//......................................................
void HardwareMapperService::printHardwareElement(int const& element){
  INFO_FILE_FUNCTION << std::endl;
  return;
}

//......................................................
void HardwareMapperService::printGeometryInfo(){
  INFO_FILE_FUNCTION << std::endl;
  INFO << std::endl;
  INFO << "DetectorName: " << geometryService->DetectorName() << std::endl;
  INFO << "TotalMass: " << geometryService->TotalMass() << std::endl;
  unsigned int Ncryostats  = geometryService->Ncryostats();
  unsigned int TotalNTPC   = geometryService->TotalNTPC();
  unsigned int Nchannels   = geometryService->Nchannels();
  unsigned int NOpChannels = geometryService->NOpChannels();
  INFO << "Ncryostats:   " << Ncryostats   << std::endl;
  INFO << "TotalNTPC:    " << TotalNTPC    << std::endl;
  INFO << "Nchannels:    " << Nchannels    << std::endl;
  INFO << "NOpChannels:  " << NOpChannels  << std::endl;
  INFO << std::endl;

  // geo::ChannelMapAlg const* channelMap =  geometryService->ChannelMap();
  // INFO << "Nchannels from channel map: " << channelMap->Nchannels() << std::endl;
  // INFO << std::endl;
    

  INFO << "Looping through cryostats" << std::endl;
  for(unsigned int cryostat = 0; cryostat < Ncryostats; cryostat++){
    unsigned int NTPC = geometryService->NTPC(cryostat);
    INFO << "  cryostat: " << cryostat << " has NTPC: " << NTPC << std::endl;
    INFO << "  Looping through tpcs" << std::endl;
    for(unsigned int tpc = 0; tpc < NTPC; tpc++){
      unsigned int Nplanes = geometryService->Nplanes(tpc, cryostat);
      geo::TPCGeo tpc_geo = geometryService->TPC(tpc, cryostat);
      INFO << "    tpc: " << tpc << " has Nplanes: " << Nplanes << std::endl;
    }
  }

  INFO << std::endl;
  INFO << "Looping through channels" << std::endl;
  for(raw::ChannelID_t channel=0; channel<Nchannels ;channel++){
    std::vector<geo::WireID> const Wires = geometryService->ChannelToWire(channel);
    //    INFO << "channel: " << channel " -> "
      
    
    /*
    SigType_t SignalType(raw::ChannelID_t const channel) const;
    View_t View(raw::ChannelID_t const channel) const;
    std::vector<geo::WireID> ChannelToWire
      (raw::ChannelID_t const channel) const;
    */
  }


  return;
}

//......................................................
void HardwareMapperService::printGeometryHelperInfo(){
  INFO_FILE_FUNCTION << std::endl;
  return;
}

//......................................................
void HardwareMapperService::printChannelMap(){
  INFO_FILE_FUNCTION << std::endl;
  
  std::set<unsigned int> offline_set;
  std::set<unsigned int> online_set;
  std::set<unsigned int> apa_set;
  std::set<unsigned int> plane_set;
  std::set<unsigned int> rce_set;
  std::set<unsigned int> rceChannel_set;
  std::set<unsigned int> regulator_set;
  std::set<unsigned int> asic_set;
  std::set<unsigned int> asicChannel_set;

  const unsigned int numChans = 2048;
  for(unsigned int online=0;online<numChans;online++){
    unsigned int offline       = channelMapService->Offline(online);
    unsigned int apa           = channelMapService->APAFromOnlineChannel(online);
    unsigned int plane         = channelMapService->PlaneFromOnlineChannel(online);
    unsigned int rce           = channelMapService->RCEFromOnlineChannel(online);
    unsigned int rceChannel    = channelMapService->RCEChannelFromOnlineChannel(online);
    unsigned int regulator     = channelMapService->RegulatorFromOnlineChannel(online);
    unsigned int asic          = channelMapService->ASICFromOnlineChannel(online);
    unsigned int asicChannel   = channelMapService->ASICChannelFromOnlineChannel(online);

    offline_set       .insert(offline    );
    apa_set           .insert(apa        );
    plane_set         .insert(plane      );
    rce_set           .insert(rce        );
    rceChannel_set    .insert(rceChannel );
    regulator_set     .insert(regulator  );
    asic_set          .insert(asic       );
    asicChannel_set   .insert(asicChannel);


    ASIC this_asic(asic,rce,0,apa);
    INFO << this_asic << std::endl;


    INFO << "online: " << online
         << " offline: " << offline
         << " apa: " << apa
         << " plane: " << plane
         << " rce: " << rce
         << " rceChannel: " << rceChannel
         << " regulator: " << regulator
         << " asic: " << asic
         << " asicChannel: " << asicChannel 
         << std::endl;
  }

  INFO << std::endl;
  INFO << "offline_set.size():     " << offline_set.size()     << std::endl;
  INFO << "apa_set.size():         " << apa_set.size()         << std::endl;
  INFO << "plane_set.size():       " << plane_set.size()       << std::endl;
  INFO << "rce_set.size():         " << rce_set.size()         << std::endl;
  INFO << "rceChannel_set.size():  " << rceChannel_set.size()  << std::endl;
  INFO << "regulator_set.size():   " << regulator_set.size()   << std::endl;
  INFO << "asic_set.size():        " << asic_set.size()        << std::endl;
  INFO << "asicChannel_set.size(): " << asicChannel_set.size() << std::endl;
  INFO << std::endl;
}

//......................................................
DEFINE_ART_SERVICE(HardwareMapperService)

