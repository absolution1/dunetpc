// TpcData.h
//
// David Adams
// January 2021
//
// Class to hold TPC data in two formats:
//   AdcData: AdcChannelData objects held in an AdcChannelDataMap
//   TPC 2D ROIs: Tpc2dRoi objects held in a vector.
// Other may be added later.
//
// The object may also hold named objects of the same type so that data can
// can be split by detector region or processing stage.
// Branches below the first level, are referenced in slash-separated path format
// e.g. "apa3/z" refers to TPC data "z" in TPC data "apa3".

#ifndef TpcData_H
#define TpcData_H

#include "dune/DuneInterface/Data/AdcChannelData.h"
#include "dune/DuneInterface/Data/Tpc2dRoi.h"
#include <string>
#include <map>
#include <set>

class TpcData {

public:

  using Name = std::string;
  using Tpc2dRoiVector = std::vector<Tpc2dRoi>;
  using TpcDataMap = std::map<Name, TpcData>;
  using AdcData = AdcChannelDataMap;
  using AdcDataPtr = std::unique_ptr<AdcChannelDataMap>;
  using TpcDataSet = std::set<TpcData>;

  // Ctor.
  TpcData();

  // Copy ctor.
  TpcData(const TpcData& rhs) =delete;

  // Accessors.
  TpcData* getParent() { return m_parent; }
  const TpcData* getParent() const { return m_parent; }
  TpcDataMap& getData() { return m_dat; }
  const TpcDataMap& getData() const { return m_dat; }
  AdcData* getAdcData() { return m_padcs.get(); }
  const AdcData* getAdcData() const { return m_padcs.get(); }
  Tpc2dRoiVector& get2dRois() { return m_2drois; }
  const Tpc2dRoiVector& get2dRois() const { return m_2drois; }

  // Add a named TPC data constituent and return its pointer.
  // If the name is target/subname, constituent subname is inserted
  // in constituent target.
  // Fails and returns null if the name is already used or the
  // target does not already exist.
  TpcData* addTpcData(Name nam);

  // Return a constituent TPC data object by name.
  // The name "" or "." returns this object.
  // If the name is target/subname, constituent subname from
  // constituent target is returned.
  TpcData* getTpcData(Name nam);
  const TpcData* getTpcData(Name nam) const;

  // Create (empty) ADC data.
  // Fails and returns nullptr if the data is already present.
  AdcData* createAdcData();

  // Delete the ADC data.
  void clearAdcData();

  // Return the ADC data searching this object and then its ancestors.
  AdcData* findAdcData();
  const AdcData* findAdcData() const;

private:

  TpcData*       m_parent;
  TpcDataMap     m_dat;
  AdcDataPtr     m_padcs;
  Tpc2dRoiVector m_2drois;

};

#endif
