// ApaChannelGroups.h
//
// David Adams
// July 2018
//
// Tool to return channel range groups for a detector with NAPA APAs
// assigned to TPC sets 0, 1, ..., NAPA-1.
//      Name    Label          Meaning
//      tpss    APAs           All TPC sets (set = pair of TPCs for one APA)
//      apas    APAs           All APAs
//     tppPs    TPC P planes   All TPC planes for view P
//     apaPs    APA P planes   All APA planes for view P
//   fembAFF    FEMB AFF       All channels in the FEMB (i.e. u, v and x)
//
// In the above
//    S is the TPC set (offline APA) number in range [0, NAPA-1]
//    A is the APA number in range [1,6]
//    P is the wire plane: u, v, z (TPC-side collection) or c (cryostat-side collection)
//      or the wire plane pair x for (z, c) or i for (u, v)
//   FF is the FEMB number in an APA in range [01,20]
//    V is a wire orientation: u, v or x (collection)
//
// Where relevant, the ranges are assigned a second label location as one of:
//   US-RaS, US-DaS, MS-RaS, MS-DaS, DS-RaS, DS-DaS
// and a third label indication the APA, e.g. APA3.
//
// Configuration parameters:
//   LogLevel: Logging level (0=none, 1=ctor, 2=every call)
//   ApaNumbers: Mapping from TPS to APA numbers. Size is NAPA.
//   IndexRangeTool: Tool that maps names here to channel ranges.
//
// Updates May 2018:
//   If ApaNumbers is empty, the detector consists of one unnumbered APA,
//   i.e. expected ranges are apau, apaz, etc.
//   Add P = z1, z2.

#ifndef ApaChannelGroups_H
#define ApaChannelGroups_H

#include "art/Utilities/ToolMacros.h"
#include "fhiclcpp/ParameterSet.h"
#include "dune/DuneInterface/Tool/IndexRangeGroupTool.h"
#include <map>

class IndexRangeTool;

class ApaChannelGroups : public IndexRangeGroupTool {

public:

  using Name = std::string;
  using NameVector = std::vector<Name>;
  using Index = IndexRangeGroup::Index;
  using IndexVector = std::vector<Index>;
  using GroupMap = std::map<Name, NameVector>;

  // Ctor.
  ApaChannelGroups(fhicl::ParameterSet const& ps);

  // Dtor.
  ~ApaChannelGroups() override =default;

  // Return a range.
  IndexRangeGroup get(Name nam) const override;

private:

  // Configuration parameters.
  Index m_LogLevel;
  IndexVector m_ApaNumbers;
  Name m_IndexRangeTool;

  // Derived from configuration.
  const IndexRangeTool* m_pIndexRangeTool =nullptr;
  GroupMap m_groups;
  GroupMap m_labels;

};


#endif
