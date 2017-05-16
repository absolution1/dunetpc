// HardwareMapperService.h
//
// Jonathan Davies j.p.davies@sussex.ac.uk
// August 2016
//
// Description: Dish out vectors of raw::ChannelID_t to the user for requested hardware element

#ifndef HARDWAREMAPPERSERVICE_H
#define HARDWAREMAPPERSERVICE_H

//The service declaration macros are defined in here
#include "art/Framework/Services/Registry/ServiceMacros.h"
//Needed to get a service handle
#include "art/Framework/Services/Registry/ServiceHandle.h"
//We need the geometry service for the hardwaremapper to work -- how else would it know what detector it was?!
#include "larcore/Geometry/Geometry.h"
//We need this for the preBeginRun thingy
#include "larcoreobj/SummaryData/RunData.h"

//We define some useful data structures (classes) in here
#include "dune/DAQTriggerSim/TriggerDataProducts/HardwareElements.h"

#include <vector>
#include <set>
#include <iostream>
#ifndef __GCCXML__
#include <iosfwd> // std::ostream
#endif

//Forward Class Declarations
namespace fhicl {
  class ParameterSet;
}

namespace art {
  class ActivityRegistry;
}

//......................................................
class HardwareMapperService{
 public:
  HardwareMapperService(fhicl::ParameterSet const& pset, art::ActivityRegistry& reg);
  
  void printASICMap(unsigned int num_asics_to_print=20);
  void printBoardMap(unsigned int num_boards_to_print=20);
  void printTPCMap(unsigned int num_tpcs_to_print=20);
  void printAPAMap(unsigned int num_apas_to_print=20);
  void printHardwareMaps(); //jpd -- prints both the tpc and apa maps

  void printGeometryInfo(); //jpd -- testing function to print geometry information once loaded

  unsigned int getNASICs() const { return fASICMap.size();}
  unsigned int getNBoards() const { return fBoardMap.size();}
  unsigned int getNAPAs() const { return fAPAMap.size();}
  unsigned int getNTPCs() const { return fTPCMap.size();}


  //jpd -- These are the main user accessible functions - dish out vectors of channel ids
  std::vector<raw::ChannelID_t> const& getTPCChannels(Hardware::ID tpc_id);
  std::vector<raw::ChannelID_t> const& getAPAChannels(Hardware::ID apa_id);

  //jpd -- For users that prefer a std::set
  std::set<raw::ChannelID_t> const& getTPCChannelsSet(Hardware::ID tpc_id);
  std::set<raw::ChannelID_t> const& getAPAChannelsSet(Hardware::ID apa_id);

  //jpd -- We register this such that it gets called just before we process a new run
  //    -- It double checks that the geometry we filled with is the same as that used 
  //    -- to produce this run
  void checkGeomVsFileDetectorName(art::Run const& run);  

 private:

  unsigned int fLogLevel;
  unsigned int fNASICsPerBoard;
  unsigned int fNBoardsPerAPA;
  std::string fServiceName;

  art::ServiceHandle<geo::Geometry> fGeometryService;

  std::string fDetectorNameFromGeometry;
  std::string fDetectorNameFromFile;

  //jpd -- these read in the geometry information and fill internal maps of ID->Hardware::Element
  void fillASICMap();
  void fillBoardMap();
  void fillTPCMap();
  void fillAPAMap();
  void fillHardwareMaps();

  Hardware::ASICMap fASICMap;
  Hardware::BoardMap fBoardMap;
  Hardware::TPCMap fTPCMap;
  Hardware::APAMap fAPAMap;
};

//......................................................
DECLARE_ART_SERVICE(HardwareMapperService, LEGACY)

#endif
