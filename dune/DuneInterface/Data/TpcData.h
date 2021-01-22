// TpcData.h
//
// David Adams
// January 2021
//
// Class to hold TPC data in a variety of formats:
//   TpcData
//   AdcChannelDataMap
//   Tpc2dRoiVector
// Each collection is indexed by name.
//
// Names in slash-separated path format refer to objects in constituent
// TpcData, e.g. "apa12/apa12x/mydat" referes to the objec(s) named "mydat"
// in object "apa12x" of object "apa12".

#ifndef TpcData_H
#define TpcData_H

#include "dune/DuneInterface/Data/AdcChannelData.h"
#include "dune/DuneInterface/Data/Tpc2dRoi.h"
#include <string>
#include <map>

class TpcData {

public:

  using Name = std::string;
  using Tpc2dRoiVector = std::vector<Tpc2dRoi>;
  using TpcDataMap = std::map<Name, TpcData>;
  using AdcData = AdcChannelDataMap;
  using AdcDataMap = std::map<Name, AdcChannelDataMap>;
  using Tpc2dRoiMap = std::map<Name, Tpc2dRoiVector>;

  // Add a named constituent and return its pointer.
  // Fails and returns null if the name is already used or the
  // target TpcData object does not already exist..
  TpcData*           addTpcData(Name nam);
  AdcChannelDataMap* addAdcData(Name nam);
  Tpc2dRoiVector*    addTpcRois(Name nam);

  // Return a non-const constituent TPC data object.
  // The name "" or "." returns this object.
  TpcData* getTpcData(Name nam);
  AdcData* getAdcData(Name nam);

private:

  TpcDataMap  m_dat;
  AdcDataMap  m_acm;
  Tpc2dRoiMap m_r2d;

};

#endif
