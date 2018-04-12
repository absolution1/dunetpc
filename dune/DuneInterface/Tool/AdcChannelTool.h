// AdcChannelTool.h

// David Adams
// April 2018
//
// Interface for tools that view or modify single objects o
// maps of ADC channel data.
//
// All methods return a DataMap object including a status that should be set
// to zero to indicate success.
//
// Modification might read or write LarSoft digits or wires, add or remove
// noise or convert raw to/from prepared.
//
// Viewing might return histograms, graphs or calculated data in its DataMap.
// Or it might write png or root files
//
// All methods have a default implementation. A typcial tool implemention
// will override one and, where appropriate, the defaults will use that
// ovewrride. If no methods are overriden or an inapppropriate call is made,
// all methods return status interfaceNotImplemented();
//
// Discussion of this interface may be found at
//   https://cdcvs.fnal.gov/redmine/issues/19661

#ifndef AdcChannelTool_H
#define AdcChannelTool_H

#include "dune/DuneInterface/AdcChannelData.h"
#include "dune/DuneInterface/Data/DataMap.h"

class AdcChannelTool {

public:

  virtual ~AdcChannelTool() =default;

  static int interfaceNotImplemented() { return -999; }

  // Modify data for one channel.
  // Default returns error.
  virtual DataMap modify(AdcChannelData&) const {
    return DataMap(interfaceNotImplemented());
  }

  // View data for one channel.
  // Default calls modify with a copy of the data.
  virtual DataMap view(const AdcChannelData& acd) const {
    AdcChannelData adctmp(acd);
    return modify(adctmp);
    //return modify(AdcChannelData(acd));
  }

  // Modify data for multiple channels.
  // Default loops over channels and appends the result for each success.
  // Returned status is the number of failures.
  virtual DataMap modifyMap(AdcChannelDataMap& acds) const {
    DataMap ret;
    int nfail = 0;
    for ( AdcChannelDataMap::value_type& iacd : acds ) {
      DataMap dm = modify(iacd.second);
      if ( dm.status() == interfaceNotImplemented() ) return DataMap(interfaceNotImplemented());
      else if ( dm.status() ) ++nfail;
      else ret += dm;
    }
    ret.setStatus(nfail);
    return ret;
  }

  // View data for multiple channels.
  // Default loops over channels and appends the result for each success.
  // Returned status is the number of failures.
  virtual DataMap viewMap(AdcChannelDataMap& acds) const;

};

//**********************************************************************

inline
DataMap AdcChannelTool::viewMap(AdcChannelDataMap& acds) const {
  DataMap ret;
  int nfail = 0;
  for ( AdcChannelDataMap::value_type& iacd : acds ) {
    DataMap dm = view(iacd.second);
    if ( dm.status() == interfaceNotImplemented() ) return DataMap(interfaceNotImplemented());
    else if ( dm.status() ) ++nfail;
    else ret += dm;
  }
  ret.setStatus(nfail);
  return ret;
}

//**********************************************************************

#endif
