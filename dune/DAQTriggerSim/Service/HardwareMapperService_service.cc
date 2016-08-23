// HardwareMapperService_service.cc
//
// Jonathan Davies j.p.davies@sussex.ac.uk
// August 2016
//
// Description: Dish out vectors of raw::ChannelID_t to the user for requested hardware element

#include "dune/DAQTriggerSim/Service/HardwareMapperService.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include <sstream>
#include <set>
#include <memory>

//......................................................
HardwareMapperService::HardwareMapperService(fhicl::ParameterSet const& pset, art::ActivityRegistry& reg)
  : fLogLevel(1), fServiceName("HardwareMapperService") {
  pset.get_if_present<unsigned int>("LogLevel", fLogLevel);
  const std::string func_name = "HardwareMapperService";
  if(fLogLevel>1) mf::LogInfo(fServiceName) << "In Function: " << func_name;


  //jpd -- This is how to register a callback function to the activity registry
  //    -- This means that upon the preBeginRun (which I assume is immediately before
  //    -- begin run) this function (HardwareMapperService::preBeginRun) will be called

  reg.sPreBeginRun.watch(this, &HardwareMapperService::checkGeomVsFileDetectorName);
  fillHardwareMaps();
}

//......................................................
void HardwareMapperService::checkGeomVsFileDetectorName(art::Run const& run){
  const std::string func_name = "checkGeomVsFileDetectorName";
  if(fLogLevel>1) mf::LogInfo(fServiceName) << "In Function: " << func_name;

  //jpd -- Check that the Detector Name matches between the Geometry and the File
  //Inspired by Geometry_service.cc
  mf::LogInfo loginfo(fServiceName);

  std::vector< art::Handle<sumdata::RunData> > rdcol;
  run.getManyByType(rdcol);
  if (rdcol.empty()) {
    mf::LogError(fServiceName) << "Can't find sumdata::RunData for geometry vs. file detector name checks\n";
    return;
  }

  fDetectorNameFromFile = rdcol.front()->DetName();

  loginfo << "This file has the following detector name:   " <<  fDetectorNameFromFile     << "\n";
  loginfo << "The geometry service has this detector name: " <<  fDetectorNameFromGeometry << "\n";
  
  if(fDetectorNameFromGeometry != fDetectorNameFromFile){
    throw cet::exception("HardwareMapperService") << "Detector Name missmatch!\n"
                                                  << "Detector Name from Geometry Service: " << fDetectorNameFromGeometry << "\n"
                                                  << "Detector Name from File: " << fDetectorNameFromFile << "\n"
                                                  << "This means that the HardwareMapperService is filled with information for the wrong Geometry\n"
                                                  << "Throwing an exception as the user should make sure the Geometry in their fcl and that used to produce their file match\n";
  }
}


//......................................................
void HardwareMapperService::fillTPCMap(){
  const std::string func_name = "fillTPCMap";
  if(fLogLevel>1) mf::LogInfo(fServiceName) << "In Function: " << func_name;
  unsigned int Nchannels   = fGeometryService->Nchannels();
  mf::LogInfo(fServiceName)  << "Filling TPC Map\n";

  for(raw::ChannelID_t channel=0; channel<Nchannels ;channel++){
    std::vector<geo::WireID> const Wires = fGeometryService->ChannelToWire(channel);
    for(auto wire : Wires){
      auto tpc_id = wire.TPC;
      //jpd -- See if we have already created a TPC object for this tpc_id in our map
      auto find_result = fTPCMap.find(tpc_id);
      std::shared_ptr<Hardware::TPC> this_tpc;
      if(find_result != fTPCMap.end()){
        //jpd -- We already have one, add this channel
        this_tpc = (*find_result).second;
      }
      else{
        //jpd -- We don't have one. Create a new one and add this channel
        this_tpc = std::shared_ptr<Hardware::TPC>(new Hardware::TPC(tpc_id));
      }
      this_tpc->addChannel(channel);
      fTPCMap[tpc_id] = this_tpc;
    }//loop over wires
  }//loop over channels
  mf::LogInfo(fServiceName)  << "Finished filling TPC Map\n";
}

//......................................................
void HardwareMapperService::fillAPAMap(){
  const std::string func_name = "fillAPAMap";
  if(fLogLevel>1) mf::LogInfo(fServiceName) << "In Function: " << func_name;
  unsigned int Nchannels   = fGeometryService->Nchannels();
  mf::LogInfo(fServiceName)  << "Filling APA Map\n";
  for(raw::ChannelID_t channel=0; channel<Nchannels ;channel++){
    std::vector<geo::WireID> const Wires = fGeometryService->ChannelToWire(channel);
    for(auto wire : Wires){
      auto apa_id = wire.TPC / 2;
      //jpd -- See if we have already created a APA object for this apa_id in our map
      auto find_result = fAPAMap.find(apa_id);
      std::shared_ptr<Hardware::APA> this_apa;
      if(find_result != fAPAMap.end()){
        //jpd -- We already have one, add this channel
        this_apa = (*find_result).second;
      }
      else{
        //jpd -- We don't have one. Create a new one and add this channel
        this_apa = std::shared_ptr<Hardware::APA>(new Hardware::APA(apa_id));
      }
      this_apa->addChannel(channel);
      fAPAMap[apa_id] = this_apa;
    }//loop over wires
  }//loop over channels
  mf::LogInfo(fServiceName)  << "Finished filling APA Map\n";
}

//......................................................
void HardwareMapperService::fillHardwareMaps(){
  const std::string func_name = "fillHardwareMaps";
  if(fLogLevel>1) mf::LogInfo(fServiceName) << "In Function: " << func_name;
  fDetectorNameFromGeometry = fGeometryService->DetectorName();

  {
    mf::LogInfo loginfo(fServiceName);
    
    loginfo << "HardwareMapperService\n";
    loginfo << "Filling Hardware Maps using geometry service\n";
    loginfo << "DetectorName: " << fGeometryService->DetectorName() << "\n";
    loginfo << "Ncryostats:   " << fGeometryService->Ncryostats()   << "\n";
    loginfo << "TotalNTPC:    " << fGeometryService->TotalNTPC()    << "\n";
    loginfo << "Nchannels:    " << fGeometryService->Nchannels()    << "\n";
    loginfo << "NOpChannels:  " << fGeometryService->NOpChannels()  << "\n";
  }//using annonymous namespace to force mf::LogInfo destructor call to flush output
  fillTPCMap();
  fillAPAMap();
}

//......................................................
void HardwareMapperService::printTPCMap(unsigned int num_tpcs_to_print){
  const std::string func_name = "printTPCMap";
  if(fLogLevel>1) mf::LogInfo(fServiceName) << "In Function: " << func_name;
  mf::LogInfo loginfo(fServiceName);
  loginfo << "Printing the first: " << num_tpcs_to_print << " TPCs\n";
  unsigned int total_channels = 0;
  unsigned int count = 0;
  for(auto this_pair : fTPCMap ){
    if(count++ >= num_tpcs_to_print) break;
    std::shared_ptr<Hardware::TPC> this_tpc = this_pair.second;
    total_channels += this_tpc->getNChannels();
    loginfo << *this_tpc << "\n";
  }
  loginfo << "\n";
  loginfo << "Total Channels: " << total_channels << "\n";
  loginfo << "\n";
}

//......................................................
void HardwareMapperService::printAPAMap(unsigned int num_apas_to_print){
  const std::string func_name = "printAPAMap";
  if(fLogLevel>1) mf::LogInfo(fServiceName) << "In Function: " << func_name;
  mf::LogInfo loginfo(fServiceName);

  loginfo << "Printing the first: " << num_apas_to_print << " APAs\n";
  unsigned int total_channels = 0;
  unsigned int count = 0;
  for(auto this_pair : fAPAMap ){
    if(count++ >= num_apas_to_print) break;
    std::shared_ptr<Hardware::APA> this_apa = this_pair.second;
    total_channels += this_apa->getNChannels();
    loginfo << *this_apa << "\n";
  }
  loginfo << "\n";
  loginfo << "Total Channels: " << total_channels << "\n";
  loginfo << "\n";
}

//......................................................
void HardwareMapperService::printHardwareMaps(){
  const std::string func_name = "printHardwareMaps";
  if(fLogLevel>1) mf::LogInfo(fServiceName) << "In Function: " << func_name;
  mf::LogInfo(fServiceName)  
    << "fDetectorName: " << fDetectorNameFromGeometry << " taken from Geometry Service at time of fillHardwareMaps()" 
    << "\n"
    << "N.B. the geometry service can reload a new geometry on beginRun(), which would be bad"
    << "\n";

  printTPCMap(getNTPCs());
  printAPAMap(getNAPAs());
}

//......................................................
void HardwareMapperService::printGeometryInfo(){
  const std::string func_name = "printGeometryInfo";
  if(fLogLevel>1) mf::LogInfo(fServiceName) << "In Function: " << func_name;
  {
    mf::LogInfo loginfo(fServiceName);
    loginfo << "Geometry Information with current Geometry Service\n";
    loginfo << "DetectorName: " << fGeometryService->DetectorName() << "\n";
    loginfo << "TotalMass: " << fGeometryService->TotalMass() << "\n";
    unsigned int Ncryostats  = fGeometryService->Ncryostats();
    unsigned int TotalNTPC   = fGeometryService->TotalNTPC();
    unsigned int Nchannels   = fGeometryService->Nchannels();
    unsigned int NOpChannels = fGeometryService->NOpChannels();
    loginfo << "Ncryostats:   " << Ncryostats   << "\n";
    loginfo << "TotalNTPC:    " << TotalNTPC    << "\n";
    loginfo << "Nchannels:    " << Nchannels    << "\n";
    loginfo << "NOpChannels:  " << NOpChannels  << "\n";
    loginfo << "\n";
  }//using annonymous namespace to force mf::LogInfo destructor (which makes flushes the output)

  printHardwareMaps();  

  return;
}

//......................................................
std::vector<raw::ChannelID_t> const& HardwareMapperService::getTPCChannels(Hardware::ID tpc_id){
  const std::string func_name = "getTPCChannels";
  mf::LogInfo loginfo(fServiceName);

  if(fLogLevel>1) loginfo << "In Function: " << func_name << "\n";
  if(fLogLevel>1) loginfo << "Finding channels for TPC: " << tpc_id << "\n";

  auto find_result = fTPCMap.find(tpc_id);
  std::shared_ptr<Hardware::TPC> this_tpc;
  
  if(find_result == fTPCMap.end()){
    mf::LogError(fServiceName) << "Failed to find this TPC: " << tpc_id << "\n";
    mf::LogError(fServiceName) << "Returning an empty vector of channels\n";
    static std::vector<raw::ChannelID_t> emptyVector;
    return emptyVector;
  }
  this_tpc = (*find_result).second;
  if(fLogLevel>1) loginfo << "Found " << *this_tpc << "\n";
  if(fLogLevel>1) loginfo << "\n";
  return this_tpc->getChannels();  
}

//......................................................
std::vector<raw::ChannelID_t> const& HardwareMapperService::getAPAChannels(Hardware::ID apa_id){
  const std::string func_name = "getAPAChannels";
  if(fLogLevel>1) mf::LogInfo(fServiceName) << "In Function: " << func_name;
  mf::LogInfo loginfo(fServiceName);
  if(fLogLevel>1) loginfo << "Finding channels for APA: " << apa_id << "\n";

  auto find_result = fAPAMap.find(apa_id);
  std::shared_ptr<Hardware::APA> this_apa;
  
  if(find_result == fAPAMap.end()){
    mf::LogError(fServiceName) << "Failed to find this APA: " << apa_id << "\n";
    mf::LogError(fServiceName) << "Returning an empty vector of channels\n";
    static std::vector<raw::ChannelID_t> emptyVector;
    return emptyVector;
  }
  this_apa = (*find_result).second;
  if(fLogLevel>1) loginfo << "Found " << *this_apa << "\n";
  if(fLogLevel>1) loginfo << "\n";
  return this_apa->getChannels();  
}

//......................................................
std::set<raw::ChannelID_t> const& HardwareMapperService::getTPCChannelsSet(Hardware::ID tpc_id){
  const std::string func_name = "getTPCChannelsSet";
  if(fLogLevel>1) mf::LogInfo(fServiceName) << "In Function: " << func_name;
  mf::LogInfo loginfo(fServiceName);
  if(fLogLevel>1) loginfo << "Finding channels for TPC: " << tpc_id << "\n";

  auto find_result = fTPCMap.find(tpc_id);
  std::shared_ptr<Hardware::TPC> this_tpc;
  
  if(find_result == fTPCMap.end()){
    mf::LogError(fServiceName) << "Failed to find this TPC: " << tpc_id << "\n";
    mf::LogError(fServiceName) << "Returning an empty set of channels\n";
    static std::set<raw::ChannelID_t> emptySet;
    return emptySet;
  }
  this_tpc = (*find_result).second;
  if(fLogLevel>1) loginfo << "Found " << *this_tpc << "\n";
  if(fLogLevel>1) loginfo << "\n";
  return this_tpc->getChannelsSet();  
}

//......................................................
std::set<raw::ChannelID_t> const& HardwareMapperService::getAPAChannelsSet(Hardware::ID apa_id){
  const std::string func_name = "getAPAChannelsSet";
  if(fLogLevel>1) mf::LogInfo(fServiceName) << "In Function: " << func_name;
  mf::LogInfo loginfo(fServiceName);
  if(fLogLevel>1) loginfo << "Finding channels for APA: " << apa_id << "\n";

  auto find_result = fAPAMap.find(apa_id);
  std::shared_ptr<Hardware::APA> this_apa;
  
  if(find_result == fAPAMap.end()){
    mf::LogError(fServiceName) << "Failed to find this APA: " << apa_id << "\n";
    mf::LogError(fServiceName) << "Returning an empty set of channels\n";
    static std::set<raw::ChannelID_t> emptySet;
    return emptySet;
  }
  this_apa = (*find_result).second;
  if(fLogLevel>1) loginfo << "Found " << *this_apa << "\n";
  if(fLogLevel>1) loginfo << "\n";
  return this_apa->getChannelsSet();  
}

//......................................................
DEFINE_ART_SERVICE(HardwareMapperService)

