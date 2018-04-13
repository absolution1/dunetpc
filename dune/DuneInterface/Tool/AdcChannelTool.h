// AdcChannelTool.h

// David Adams
// April 2018
//
// Interface for tools that view or modify single objects or
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
// All methods have a default implementation. A typical tool implemention
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
  virtual DataMap modify(AdcChannelData&) const;

  // View data for one channel.
  // Default calls modify with a copy of the data and returns the
  // DataMap from that call.
  virtual DataMap view(const AdcChannelData& acd) const;

  // Modify data for multiple channels.
  // Default calls modify for each channel and appends the result for each success.
  // The status is set to the number of failures.
  virtual DataMap modifyMap(AdcChannelDataMap& acds) const;

  // View data for multiple channels.
  // Default call view for each channel and appends the result for each success.
  // The status is set to the number of failures.
  virtual DataMap viewMap(AdcChannelDataMap& acds) const;

};

//**********************************************************************
// Definitions for the above declarations.
//**********************************************************************

inline
DataMap AdcChannelTool::modify(AdcChannelData&) const {
  return DataMap(interfaceNotImplemented());
}

//**********************************************************************

inline
DataMap AdcChannelTool::view(const AdcChannelData& acd) const {
  AdcChannelData adctmp(acd);
  return modify(adctmp);
}

//**********************************************************************

inline
DataMap AdcChannelTool::modifyMap(AdcChannelDataMap& acds) const {
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
