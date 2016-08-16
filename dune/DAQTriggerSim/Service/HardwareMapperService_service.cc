// HardwareMapperService_service.cc
//
// Jonathan Davies j.p.davies@sussex.ac.uk
// August 2016
//
// Description

#include "dune/DAQTriggerSim/Service/HardwareMapperService.h"
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Principal/Event.h"

#include <sstream>
#include <set>
#include <memory>

//......................................................
HardwareMapperService::HardwareMapperService(fhicl::ParameterSet const& pset, art::ActivityRegistry& reg){
  INFO_FILE_FUNCTION << std::endl;
  std::string detectorName = fGeometryService->DetectorName();
  //FIXME set the number of channels per ASIC and Board here

  fillHardwareMaps();
}

//......................................................
void HardwareMapperService::fillTPCMap(){
  INFO_FILE_FUNCTION << std::endl;
  unsigned int Nchannels   = fGeometryService->Nchannels();
  INFO << "Filling TPC Map" << std::endl;
  INFO << "Looping through channels" << std::endl;

  for(raw::ChannelID_t channel=0; channel<Nchannels ;channel++){
    std::vector<geo::WireID> const Wires = fGeometryService->ChannelToWire(channel);
    for(auto wire : Wires){
      auto tpc_id = wire.TPC;
      //See if we have already created a TPC object for this tpc_id in our map
      auto find_result = fTPCMap.find(tpc_id);
      std::shared_ptr<Hardware::TPC> this_tpc;
      if(find_result != fTPCMap.end()){
        //We already have one, add this channel
        this_tpc = (*find_result).second;
        //        INFO << "Found: " << *this_tpc << std::endl;
      }
      else{
        //We don't have one. Create a new one and add this channel
        this_tpc = std::shared_ptr<Hardware::TPC>(new Hardware::TPC(tpc_id));
        //        INFO << "Created: " << *this_tpc << std::endl;
      }
      this_tpc->addChannel(channel);
      fTPCMap[tpc_id] = this_tpc;
    }//loop over wires
  }//loop over channels
  INFO << "Finished filling TPC Map" << std::endl;
  INFO << std::endl;
}

//......................................................
void HardwareMapperService::fillAPAMap(){
  INFO_FILE_FUNCTION << std::endl;
  unsigned int Nchannels   = fGeometryService->Nchannels();
  INFO << "Looping through channels" << std::endl;
  INFO << "Filling APA Map" << std::endl;
  for(raw::ChannelID_t channel=0; channel<Nchannels ;channel++){
    std::vector<geo::WireID> const Wires = fGeometryService->ChannelToWire(channel);
    for(auto wire : Wires){
      auto apa_id = wire.TPC / 2;
      //See if we have already created a APA object for this apa_id in our map
      auto find_result = fAPAMap.find(apa_id);
      std::shared_ptr<Hardware::APA> this_apa;
      if(find_result != fAPAMap.end()){
        //We already have one, add this channel
        this_apa = (*find_result).second;
        //        INFO << "Found: " << *this_apa << std::endl;
      }
      else{
        //We don't have one. Create a new one and add this channel
        this_apa = std::shared_ptr<Hardware::APA>(new Hardware::APA(apa_id));
        //        INFO << "Created: " << *this_apa << std::endl;
      }
      this_apa->addChannel(channel);
      fAPAMap[apa_id] = this_apa;
    }//loop over wires
  }//loop over channels
  INFO << "Finished filling TPC Map" << std::endl;
  INFO << std::endl;
}

//......................................................
void HardwareMapperService::fillHardwareMaps(){
  fillTPCMap();
  fillAPAMap();
}

//......................................................
void HardwareMapperService::printTPCMap(unsigned int num_tpcs_to_print){
  INFO_FILE_FUNCTION << std::endl;
  INFO << "Printing the first: " << num_tpcs_to_print << " TPCs" << std::endl;
  unsigned int count = 0;
  for(auto this_pair : fTPCMap ){
    if(count++ >= num_tpcs_to_print) break;
    std::shared_ptr<Hardware::TPC> this_tpc = this_pair.second;
    INFO << *this_tpc << std::endl;
  }
  INFO << std::endl;
}

//......................................................
void HardwareMapperService::printAPAMap(unsigned int num_apas_to_print){
  INFO_FILE_FUNCTION << std::endl;
  INFO << "Printing the first: " << num_apas_to_print << " APAs" << std::endl;
  unsigned int count = 0;
  for(auto this_pair : fAPAMap ){
    if(count++ >= num_apas_to_print) break;
    std::shared_ptr<Hardware::APA> this_apa = this_pair.second;
    INFO << *this_apa << std::endl;
  }
  INFO << std::endl;
}

//......................................................
void HardwareMapperService::printGeometryInfo(){
  INFO_FILE_FUNCTION << std::endl;
  INFO << std::endl;
  INFO << "DetectorName: " << fGeometryService->DetectorName() << std::endl;
  INFO << "TotalMass: " << fGeometryService->TotalMass() << std::endl;
  unsigned int Ncryostats  = fGeometryService->Ncryostats();
  unsigned int TotalNTPC   = fGeometryService->TotalNTPC();
  unsigned int Nchannels   = fGeometryService->Nchannels();
  unsigned int NOpChannels = fGeometryService->NOpChannels();
  INFO << "Ncryostats:   " << Ncryostats   << std::endl;
  INFO << "TotalNTPC:    " << TotalNTPC    << std::endl;
  INFO << "Nchannels:    " << Nchannels    << std::endl;
  INFO << "NOpChannels:  " << NOpChannels  << std::endl;
  INFO << std::endl;

  printTPCMap();
  printAPAMap();

  std::set<raw::ChannelID_t> apa_channel_set = getAPAChannels(0);
  apa_channel_set = getAPAChannels(1e6);

  std::set<raw::ChannelID_t> tpc_channel_set = getTPCChannels(0);
  tpc_channel_set = getTPCChannels(1e6);


  return;
}

//......................................................
std::set<raw::ChannelID_t> const& HardwareMapperService::getTPCChannels(Hardware::ID tpc_id){
  INFO_FILE_FUNCTION << std::endl;
  INFO << "Finding channels for TPC: " << tpc_id << std::endl;

  auto find_result = fTPCMap.find(tpc_id);
  std::shared_ptr<Hardware::TPC> this_tpc;
  
  if(find_result == fTPCMap.end()){
    ERROR << "Failed to find this TPC: " << tpc_id << std::endl;
    ERROR << "Returning an empty set of channels" << std::endl;
    static std::set<raw::ChannelID_t> emptySet;
    return emptySet;
  }
  this_tpc = (*find_result).second;
  INFO << "Found " << *this_tpc << std::endl;
  INFO << std::endl;
  return this_tpc->getChannels();  
}

//......................................................
std::set<raw::ChannelID_t> const& HardwareMapperService::getAPAChannels(Hardware::ID apa_id){
  INFO_FILE_FUNCTION << std::endl;
  INFO << "Finding channels for APA: " << apa_id << std::endl;

  auto find_result = fAPAMap.find(apa_id);
  std::shared_ptr<Hardware::APA> this_apa;
  
  if(find_result == fAPAMap.end()){
    ERROR << "Failed to find this APA: " << apa_id << std::endl;
    ERROR << "Returning an empty set of channels" << std::endl;
    static std::set<raw::ChannelID_t> emptySet;
    return emptySet;
  }
  this_apa = (*find_result).second;
  INFO << "Found " << *this_apa << std::endl;
  INFO << std::endl;
  return this_apa->getChannels();  
}

//......................................................
void HardwareMapperService::setNumChannelsPerBoard(unsigned int N, bool refillMap){
  INFO_FILE_FUNCTION << std::endl;
  fNumChannelsPerBoard = N;
}
//......................................................
void HardwareMapperService::setNumChannelsPerASIC(unsigned int N, bool refillMap){
  INFO_FILE_FUNCTION << std::endl;
  fNumChannelsPerASIC = N;
}

//......................................................
DEFINE_ART_SERVICE(HardwareMapperService)

