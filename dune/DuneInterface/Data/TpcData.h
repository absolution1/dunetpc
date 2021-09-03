// TpcData.h
//
// David Adams
// January 2021
//
// Class to hold TPC data in two formats:
//   AdcData: Vector of shared pointers to maps of AdcChannelData.
//   TPC 2D ROIs: Tpc2dRoi objects held in a vector.
// Other may be added later.
//
// The object may also hold named objects of the same type so that data can
// can be split by detector region (e.g. TPC plane) or processing stage.
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

  using Index = unsigned int;
  using Name = std::string;
  using Tpc2dRoiVector = std::vector<Tpc2dRoi>;
  using TpcDataMap = std::map<Name, TpcData>;
  using TpcDataVector = std::vector<TpcData*>;
  using AdcData = AdcChannelDataMap;
  using AdcDataPtr = std::shared_ptr<AdcChannelDataMap>;
  using AdcDataVector = std::vector<AdcDataPtr>;

  // Default ctor.
  TpcData();

  // Ctor to create npla AdcChannelData maps.
  TpcData(Index npla);

  // Ctor from ADC data.
  TpcData(const AdcDataVector& adcs);

  // Copy ctor.
  TpcData(const TpcData& rhs) =delete;

  // Accessors.
  TpcData* getParent() { return m_parent; }
  const TpcData* getParent() const { return m_parent; }
  TpcDataMap& getData() { return m_dat; }
  const TpcDataMap& getData() const { return m_dat; }
  AdcDataVector& getAdcData() { return m_adcs; }
  const AdcDataVector& getAdcData() const { return m_adcs; }
  Tpc2dRoiVector& get2dRois() { return m_2drois; }
  const Tpc2dRoiVector& get2dRois() const { return m_2drois; }

  // Add a named TPC data constituent and return its pointer.
  // If the name is target/subname, constituent subname is inserted
  // in constituent target.
  // Fails and returns null if the name is already used or the
  // target does not already exist.
  // If copyAdcData is true, the new object copies the AdcData (pointers).
  TpcData* addTpcData(Name nam, bool copyAdcData =true);

  // Return a constituent TPC data object by name.
  // The name "" or "." returns this object.
  // If the name is target/subname, constituent subname from
  // constituent target is returned.
  TpcData* getTpcData(Name nam);
  const TpcData* getTpcData(Name nam) const;

  // Appends named constituent TPC data objects to a vector.
  // The name "" or "." appends this object.
  // If the name is a/b/.../z all matching constituents are returned
  // with "*" used to select all at any level.
  // Returns non-zero if any errors arise.
  Index getTpcData(Name nam, TpcDataVector& out);

  // Add ADC data.
  // If updateParent is true, the same object is added to all ancestors.
  AdcDataPtr createAdcData(bool updateParent =true);
  AdcDataPtr addAdcData(AdcDataPtr padc, bool updateParent =true);

  // Delete the ADC data. References in acestor and descendants are not affected.
  void clearAdcData();

  // Print a brief description of thid object:
  // # ADC maps and channels for each.
  // # 2D ROIS
  // # constituents
  // Same for each constituent to the indicated depth.
  std::ostream& print(Name prefix, Index depth =10) const;

private:

  TpcData*       m_parent;
  TpcDataMap     m_dat;
  AdcDataVector  m_adcs;
  Tpc2dRoiVector m_2drois;

};

#endif
