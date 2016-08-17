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

#include "larcore/Geometry/Geometry.h"

#include <vector>
#include <set>
#include <iostream>
#include <iosfwd>

//Forward Class Declarations
namespace fhicl {
  class ParameterSet;
}

namespace art {
  class ActivityRegistry;
}


//jpd -- these are useful for now. Should be replaced by message facility type streams later
#define INFO  std::cerr << "INFO   :   "
#define ERROR std::cerr << "ERROR  :   "
#define INFO_FUNCTION std::cerr << "INFO FN: " << __PRETTY_FUNCTION__ << " "

//......................................................
namespace Hardware{
  typedef unsigned int ID;

  class Element
  {
  public:
    Element(ID id, std::string this_type) : fID(id), fType(this_type) {}
    ID const& getID() const{ return fID;}
    std::string const& getType() const{ return fType; }
    std::vector<raw::ChannelID_t> const& getChannels() const{ return fChannelIDs;}
    std::set<raw::ChannelID_t> const& getChannelsSet() const{ return fChannelIDsSet;}
    size_t getNChannels() const{ return fChannelIDs.size();}
    size_t getNChannelsSet() const{ return fChannelIDsSet.size();}

    void addChannel(raw::ChannelID_t channel){ 
      //jpd -- Only add channel to the vector if it is not already in the set
      //addChannelToSet returns true if channel was no already in the set and false otherwise
      if(addChannelToSet(channel)) fChannelIDs.push_back(channel);
    }

    bool operator<( const Element& rhs ) const {
      return this->getID() < rhs.getID();
    }
    friend std::ostream & operator << (std::ostream &os,  Element &rhs){
      os << rhs.getType() << ": " << rhs.getID() << " has " << rhs.getNChannels() << " channels of which " << rhs.getNChannelsSet() << " are unique";
      std::set<raw::ChannelID_t> channels = rhs.getChannelsSet();
      unsigned int max_num_channels = 10;
      unsigned int this_channel_num = 0;
      for(auto channel : channels){
        if(this_channel_num==0)  os << ":";
        if(this_channel_num++ >= max_num_channels) break;
        os << " " << channel;
      }
      return os;
    }

  private:
    ID fID;
    std::string fType;
    std::vector<raw::ChannelID_t> fChannelIDs;
    std::set<raw::ChannelID_t> fChannelIDsSet;

    //jpd -- Returns true if channel was not already in set and was inserted
    //        false if channel was already in set / there was a problem
    bool addChannelToSet(raw::ChannelID_t this_channel){
      auto result = fChannelIDsSet.insert(this_channel);
      if(result.first != fChannelIDsSet.end()) return result.second;
      else return false;
    }
  };
  
  class ASIC      : public Element{ public: ASIC    (ID id) : Element(id, "ASIC"     ) {} };
  class Board     : public Element{ public: Board   (ID id) : Element(id, "Board"    ) {} };
  class TPC       : public Element{ public: TPC     (ID id) : Element(id, "TPC"      ) {} };
  class APA       : public Element{ public: APA     (ID id) : Element(id, "APA"      ) {} };
  class APAGroup  : public Element{ public: APAGroup(ID id) : Element(id, "APAGroup" ) {} };
  class Cryostat  : public Element{ public: Cryostat(ID id) : Element(id, "Cryostat" ) {} };

  using ASICMap  =   std::map<ID, std::shared_ptr<ASIC >>;
  using BoardMap =   std::map<ID, std::shared_ptr<Board>>;
  using TPCMap   =   std::map<ID, std::shared_ptr<TPC  >>;
  using APAMap   =   std::map<ID, std::shared_ptr<APA  >>;

  
}

//......................................................
class HardwareMapperService{
 public:
  HardwareMapperService(fhicl::ParameterSet const& pset, art::ActivityRegistry& reg);
  
  void printTPCMap(unsigned int num_tpcs_to_print=10);
  void printAPAMap(unsigned int num_apas_to_print=10);
  void printHardwareMaps(); //jpd -- prints both the tpc and apa maps

  void printGeometryInfo(); //jpd -- testing function to print geometry information once loaded

  unsigned int getNAPAs() const { return fAPAMap.size();}
  unsigned int getNTPCs() const { return fTPCMap.size();}

  //jpd -- These are the main user accessible functions - dish out vectors of channel ids
  std::vector<raw::ChannelID_t> const& getTPCChannels(Hardware::ID tpc_id);
  std::vector<raw::ChannelID_t> const& getAPAChannels(Hardware::ID apa_id);

  //jpd -- For users that prefer a std::set
  std::set<raw::ChannelID_t> const& getTPCChannelsSet(Hardware::ID tpc_id);
  std::set<raw::ChannelID_t> const& getAPAChannelsSet(Hardware::ID apa_id);

 private:

  unsigned int fLogLevel;

  art::ServiceHandle<geo::Geometry> fGeometryService;

  //jpd -- these read in the geometry information and fill internal maps of ID->Hardware::Element
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
