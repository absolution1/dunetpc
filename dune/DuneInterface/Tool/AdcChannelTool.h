// AdcChannelTool.h

// David Adams
// April 2018
//
// Interface for tools that view or update single objects or
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
// All methods have a default implementation so implementers to not have
// implement all.
//
// A typical viewer, i.e. a tool that does not modify the received data, will
// override view(acd) with the data handling and override updateWithView() to
// return true so the same handling is obtained with a call to moodify(acds).
//
// A typical updater, i.e. a tool that does modify its received data, will
// override update(acd) with that modification. If it overrides viewWithModify()
// to return true, then the default view(acd) will call modify(acd) with a
// temporary copy of the received data.
//
// By default, a call to viewMap(acds) loops over objects it its map, passes
// each in turn to view(acd) and concatenates the DataMap results to form the
// combined result (DataMap). If updateWithView() is true, the the call is
// instead forwarded to updateMap(acds).
//
// The default for updateMap(acds) does the same using update(acd) unless
// viewWithModify() is true in which case the call is forwarded to viewMap(acds).
//
// Implementers which handle multiple channels, e.g. correlated noise removal
// or event display, can instead override viewMap(acds) or updateMap(acds).
// The single-object methods view(acd) and update(acd) will return status
// interfaceNotImplemented().
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

  // Methods to update or view data. Implementer should override
  // at least one of these.

  // Modify data for one channel.
  // Default returns error.
  virtual DataMap update(AdcChannelData&) const;

  // View data for one channel.
  // Default calls update with a copy of the data and returns the
  // DataMap from that call.
  virtual DataMap view(const AdcChannelData& acd) const;

  // Modify data for multiple channels.
  // Default calls update for each channel and appends the result for each success.
  // The status is set to the number of failures.
  virtual DataMap updateMap(AdcChannelDataMap& acds) const;

  // View data for multiple channels.
  // Default call view for each channel and appends the result for each success.
  // The status is set to the number of failures.
  virtual DataMap viewMap(const AdcChannelDataMap& acds) const;

  // Methods indicating if calls to update should be forwarded to view
  // or vice versa. Implementers should override neither or one  to
  // return true.

  // If this is true, default calls to update are forwarded to view
  // and default calls to updateMap are forwarded to viewMap.
  virtual bool updateWithView() const { return false; }

  // If this is true, default calls to view are forwarded to view
  // and default calls to viewMap are forwarded to updateMap.
  // In both cases, the passed dat is first copied.
  virtual bool viewWithUpdate() const { return false; }

};

//**********************************************************************
// Definitions for the above declarations.
//**********************************************************************

inline
DataMap AdcChannelTool::update(AdcChannelData& acd) const {
  if ( updateWithView() ) return view(acd);
  return DataMap(interfaceNotImplemented());
}

//**********************************************************************

inline
DataMap AdcChannelTool::view(const AdcChannelData& acd) const {
  if ( viewWithUpdate() ) {
    AdcChannelData adctmp(acd);
    return update(adctmp);
  }
  return DataMap(interfaceNotImplemented());
}

//**********************************************************************

inline
DataMap AdcChannelTool::updateMap(AdcChannelDataMap& acds) const {
  if ( updateWithView() ) return viewMap(acds);
  DataMap ret;
  int nfail = 0;
  for ( AdcChannelDataMap::value_type& iacd : acds ) {
    DataMap dm = update(iacd.second);
    if ( dm.status() == interfaceNotImplemented() ) return DataMap(interfaceNotImplemented());
    else if ( dm.status() ) ++nfail;
    else ret += dm;
  }
  ret.setStatus(nfail);
  return ret;
}

//**********************************************************************

inline
DataMap AdcChannelTool::viewMap(const AdcChannelDataMap& acds) const {
  if ( viewWithUpdate() ) {
    AdcChannelDataMap adcstmp(acds);
    return updateMap(adcstmp);
  }
  DataMap ret;
  int nfail = 0;
  for ( const AdcChannelDataMap::value_type& iacd : acds ) {
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
