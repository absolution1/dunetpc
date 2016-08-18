// HarwareElements.h
//
// Jonathan Davies j.p.davies@sussex.ac.uk
// August 2016
//
// Description: Define some data structures (classes) that will be used to define pieces of Hardware in the detector

#ifndef HARDWAREELEMENTS_H
#define HARDWAREELEMENTS_H
#include <vector>
#include <set>
#include <iostream>
#ifndef __GCCXML__
#include <iosfwd> // std::ostream
#endif

//jpd -- Needed to define things like raw::ChannelID_t
#include "larcore/Geometry/Geometry.h"
//......................................................
namespace Hardware{
  typedef unsigned int ID;

  //jpd -- A Hardware::HardwareID is a class that has a type ("APA", "TPC", "Board"), an ID (i.e. this is "APA" #1) 

  class HardwareID{
  public: 
    //    HardwareID(HardwareID old) : fID(old.getID()), fType(old.getType()) {}
    HardwareID(ID id, std::string this_type) : fID(id), fType(this_type) {}
    ID const& getID() const{ return fID;}
    std::string const& getType() const{ return fType; }

    bool operator<( const HardwareID& rhs ) const {
      return this->getID() < rhs.getID();
    }

    friend std::ostream & operator << (std::ostream &os, HardwareID const &rhs){
      os << rhs.getType() << ": " << rhs.getID();
      return os;
    }

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

    void addChannel(raw::ChannelID_t channel){ 
      //jpd -- Only add channel to the vector if it is not already in the set
      //addChannelToSet returns true if channel was no already in the set and false otherwise
      if(addChannelToSet(channel)) fChannelIDs.push_back(channel);
    }

    friend std::ostream & operator << (std::ostream &os,  Element const &rhs){
      HardwareID const& base = rhs;
      os << base << " has "<< rhs.getNChannels() << " channels of which " << rhs.getNChannelsSet() << " are unique";

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
