// HardwareMapperService.h
//
// Jonathan Davies j.p.davies@sussex.ac.uk
// August 2016
//
// Description

#ifndef HARDWAREMAPPERSERVICE_H
#define HARDWAREMAPPERSERVICE_H

//The service declaration macros are defined in here
#include "art/Framework/Services/Registry/ServiceMacros.h"

#include "art/Framework/Services/Registry/ServiceHandle.h"

//Testing
#include "lbne-raw-data/Services/ChannelMap/ChannelMapService.h"
#include "larcore/Geometry/Geometry.h"

#include <vector>
#include <set>
#include <iostream>

//Forward Class Declarations
namespace fhicl {
  class ParameterSet;
}

namespace art {
  class ActivityRegistry;
}

#define INFO  std::cerr << "INFO   : "
#define ERROR std::cerr << "ERROR  : "
#define INFO_FILE_FUNCTION  std::cerr << "INFO   : " << __FILE__ << " : " << __FUNCTION__


//......................................................
struct ASIC{
public:
  ASIC(unsigned int ID=fUnset, unsigned int RCE=fUnset, unsigned int COB=fUnset, unsigned int APA=fUnset, unsigned int Volume=fUnset)
  : fID(ID), fRCE(RCE), fCOB(COB), fAPA(APA), fVolume(Volume){ ;}

  bool addChannelID(unsigned int channel){
    auto result = fChannelIDSet.insert(channel);
    if( result.first != fChannelIDSet.end() && result.second) return true;
    else return false;
  }

  std::set<unsigned int> const& getChannelIDSet(){ return fChannelIDSet;}

  //FIXME -- need to think about this
  const bool operator < ( ASIC const & rhs ) const
  {
    if( fID != rhs.fID ) return fID < rhs.fID;
    if( fRCE != rhs.fRCE ) return fRCE < rhs.fRCE;
    return false;
  }

  friend std::ostream& operator<< (std::ostream &os, ASIC const &rhs){
    os << "ASIC"
       << " ID: " << rhs.fID
       << " RCE: " << rhs.fRCE
       << " COB: " << rhs.fCOB
       << " APA: " << rhs.fAPA
       << " Volume: " << rhs.fVolume
       << " Channels: " << rhs.fChannelIDSet.size();
    return os;
  }

private:
  unsigned int fID;
  unsigned int fRCE;
  unsigned int fCOB;
  unsigned int fAPA;
  unsigned int fVolume;
  static const unsigned int fUnset = std::numeric_limits<unsigned int>::max();
  std::set<unsigned int> fChannelIDSet;  
};

//......................................................
class HardwareMapperService{
 public:
  HardwareMapperService(fhicl::ParameterSet const& pset, art::ActivityRegistry& reg);
  //  ~HardwareMapperService() = default;
  
  int getHardwareElement(int const& element, std::vector<int>& channelVec);
  void printHardwareElement(int const& element);
  void printChannelMap();
  void printGeometryInfo();
  void printGeometryHelperInfo();

 private:
  art::ServiceHandle<lbne::ChannelMapService> channelMapService; //FIXME -- just for testing
  art::ServiceHandle<geo::Geometry> geometryService; //FIXME -- just for testing
  //Geometry Helper service handle -- FIXME What the duck does this actually do?

};

//......................................................
DECLARE_ART_SERVICE(HardwareMapperService, LEGACY)

#endif
