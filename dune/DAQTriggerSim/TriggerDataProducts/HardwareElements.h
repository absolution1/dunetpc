// HarwareElements.h
//
// Jonathan Davies j.p.davies@sussex.ac.uk
// August 2016
//
// Description: Define some data structures (classes) that will be used to define pieces of Hardware in the detector

#ifndef HARDWAREELEMENTS_H
#define HARDWAREELEMENTS_H

#ifndef __GCCXML__
#include <memory>
#include <map>
#include <vector>
#include <set>
#include <iostream>
#include <iosfwd> // std::ostream
#endif //__GCCXML__

#include "larcoreobj/SimpleTypesAndConstants/RawTypes.h"//jpd -- needed for raw::ChannelID_t
//......................................................
namespace Hardware{
  typedef unsigned int ID;

  //jpd -- A Hardware::HardwareID is a class that has a type ("APA", "TPC", "Board"), an ID (i.e. this is "APA" #1) 

  class HardwareID{
  public: 
    HardwareID(): fID(0), fType("Unknown") {}
    HardwareID(ID id, std::string this_type) : fID(id), fType(this_type) {}
    ID const& getID() const{ return fID;}
    std::string const& getType() const{ return fType; }

#ifndef __GCCXML__
    friend bool operator<( const HardwareID& lhs, const HardwareID& rhs ) {
      if( lhs.getTypeOrder() == rhs.getTypeOrder())
        return lhs.getID() < rhs.getID();
      else return lhs.getTypeOrder() < rhs.getTypeOrder();
    }

    friend std::ostream & operator << (std::ostream &os, HardwareID const &rhs){
      os << rhs.getType() << ": " << rhs.getID();
      return os;
    }

  protected:
    //jpd -- Helper for operator< (less than). The ordering is APA->manyTPC->manyBoard->manyASIC
    int getTypeOrder() const{
      if(fType == "APA") return 0;
      if(fType == "TPC") return 1;
      if(fType == "Board") return 2;
      if(fType == "ASIC") return 3;
      else return -1;
    }


#endif //__GCCXML__


  private:
    ID fID;
    std::string fType;
  };

  static const HardwareID UnknownID = HardwareID(0, "UnknownID");

  //jpd -- I don't think ROOT's dictionaries need to know about these derived classes for now
#ifndef __GCCXML__
//jpd -- A Hardware::Element is the same as a Hardware::HardwareID but also contains a list of channels 

  class Element : public HardwareID
  {
  public:
    Element(ID id, std::string this_type) : HardwareID(id, this_type) {}
    std::vector<raw::ChannelID_t> const& getChannels() const{ return fChannelIDs;}
    std::set<raw::ChannelID_t> const& getChannelsSet() const{ return fChannelIDsSet;}
    size_t getNChannels() const{ return fChannelIDs.size();}
    size_t getNChannelsSet() const{ return fChannelIDsSet.size();}

    std::vector<HardwareID> const& getHardwareIDs() const{ return fHardwareIDs;}
    std::set<HardwareID> const& getHardwareIDsSet() const{ return fHardwareIDsSet;}
    size_t getNHardwareIDs() const{ return fHardwareIDs.size();}
    size_t getNHardwareIDsSet() const{ return fHardwareIDsSet.size();}

    void addChannel(raw::ChannelID_t channel){ 
      //jpd -- Only add channel to the vector if it is not already in the set
      //addChannelToSet returns true if channel was no already in the set and false otherwise
      if(addChannelToSet(channel)) fChannelIDs.push_back(channel);
    }

    void addHardwareID(HardwareID id){
      //jpd -- Only add hardwareID to the vector if it is not already in the set
      //addHardwareIDToSet returns true if hardwareID was no already in the set and false otherwise
      if(addHardwareIDToSet(id)) fHardwareIDs.push_back(id);
    }

    friend std::ostream & operator << (std::ostream &os,  Element const &rhs){
      HardwareID const& base = rhs;
      os << base << " has "<< rhs.getNChannels() << " channels";

      std::set<raw::ChannelID_t> channels = rhs.getChannelsSet();
      unsigned int max_num_channels = 16;
      unsigned int this_channel_num = 0;
      for(auto channel : channels){
        if(this_channel_num==0)  os << ":";
        if(this_channel_num++ >= max_num_channels) { 
          os << " ...";
          break;
        }
        os << " " << channel;
      }

      os << "\n";
      os << "Contains: " << rhs.getNHardwareIDs() << " pieces of hardware\n";
      for(auto hardwareid : rhs.getHardwareIDsSet()){
        os << hardwareid << "\n";
      }

      return os;
    }

  private:
    std::vector<raw::ChannelID_t> fChannelIDs;
    std::set<raw::ChannelID_t> fChannelIDsSet;
    std::set<HardwareID> fHardwareIDsSet;
    std::vector<HardwareID> fHardwareIDs;


    //jpd -- Returns true if channel was not already in set and was inserted
    //        false if channel was already in set / there was a problem
    bool addChannelToSet(raw::ChannelID_t this_channel){
      auto result = fChannelIDsSet.insert(this_channel);
      if(result.first != fChannelIDsSet.end()) return result.second;
      else return false;
    }

    //jpd -- Returns true if channel was not already in set and was inserted
    //        false if channel was already in set / there was a problem
    bool addHardwareIDToSet(HardwareID this_hardwareid){
      auto result = fHardwareIDsSet.insert(this_hardwareid);
      if(result.first != fHardwareIDsSet.end()) return result.second;
      else return false;
    }
  };


  //jpd -- These derived classes are specialisations of a Hardware::Element
  class ASIC      : public Element{ public: ASIC    (ID id) : Element(id, "ASIC"     ) {} };
  class Board     : public Element{ public: Board   (ID id) : Element(id, "Board"    ) {} };
  class TPC       : public Element{ public: TPC     (ID id) : Element(id, "TPC"      ) {} };
  class APA       : public Element{ public: APA     (ID id) : Element(id, "APA"      ) {} };
  class APAGroup  : public Element{ public: APAGroup(ID id) : Element(id, "APAGroup" ) {} };
  class Cryostat  : public Element{ public: Cryostat(ID id) : Element(id, "Cryostat" ) {} };

  //jpd -- This just defines a bunch of shorthands for maps of Hardware::ID to fancy pointer to specialised Hardware::Element
  //jpd -- Used std::shared_ptr so we don't have to do the memory management
  using ASICMap  =   std::map<ID, std::shared_ptr<ASIC >>;
  using BoardMap =   std::map<ID, std::shared_ptr<Board>>;
  using TPCMap   =   std::map<ID, std::shared_ptr<TPC  >>;
  using APAMap   =   std::map<ID, std::shared_ptr<APA  >>;

#endif //#ifndef __GCCXML__
}

#endif //#ifndef HARDWAREELEMENTS_H
