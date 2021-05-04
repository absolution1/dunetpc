// TpcDataTool.h

// David Adams
// April 2018
//
// Interface for tools that view or update TpcData objects.
//
// It inherits from AdcChannelTool and the default implementation here
// calls the ADC channel map methods of that class.

#ifndef TpcDataTool_H
#define TpcDataTool_H

#include "dune/DuneInterface/Tool/AdcChannelTool.h"
#include "dune/DuneInterface/Data/TpcData.h"

class TpcDataTool : public AdcChannelTool {

public:

  // Modify TpcData.
  virtual DataMap updateTpcData(TpcData&) const;

  // View TpcData.
  virtual DataMap viewTpcData(const TpcData&) const;

  // Method indicating how unimplemented calls to TpcData should be forwarded.
  //   0: No forward.
  //   1: AdcChannelDataMap
  virtual int forwardTpcData() const { return 1; }

};

//**********************************************************************
// Definitions for the above declarations.
//**********************************************************************

inline
DataMap TpcDataTool::updateTpcData(TpcData& tpd) const {
  if ( forwardTpcData() == 1 ) {
    DataMap dm;
    for ( TpcData::AdcDataPtr padc : tpd.getAdcData() ) {
      if ( padc ) dm = updateMap(*padc);
    }
    return dm;
  } else if ( updateWithView() ) {
    return viewTpcData(tpd);
  }
  return DataMap(interfaceNotImplemented());
}

//**********************************************************************

inline
DataMap TpcDataTool::viewTpcData(const TpcData& tpd) const {
  if ( forwardTpcData() == 1 ) {
    DataMap dm;
    for ( TpcData::AdcDataPtr padc : tpd.getAdcData() ) {
      if ( padc ) dm = viewMap(*padc);
    }
    return dm;
  }
  // We do not copy TpcData.
  return DataMap(interfaceNotImplemented());
}

//**********************************************************************

#endif
