// ProtoDuneChannelRanges.h
//
// David Adams
// July 2018
//
// Tool to return channel ranges for protoDUNE. The following are created:
//       all - all channels
//      tpsS - for TPC set S, e.g. tps0
//     tppSP - TPC plane or plane pair P in TPC set S, e.g. tps0z
//      apaA - APA, e.g. apa3
//  fembAFFV - FEMB AFF orientation x, eg femb302x
//
// In the above
//    S is the TPC set (offline APA) number in range [0,5]
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
// Note there is the option to append another index range tool which can
// oveerride or extend the above set of ranges.
//
// The TPC set numbering follows the convention:
//
//  --->   TPS1  TPS3  TPS5
//  beam   TPS0  TPS2  TPS4
//
// and APA numbering is
//
//  --->   APA5  APA6  APA4
//  beam   APA3  APA2  APA1
//
// both viewed from above.
//
// Parameters:
//   LogLevel - Message logging level (0=none, 1=ctor, 2=each call, ...)
//   ExtraRanges - Name of tool with additional ranges. Blank for none.

#ifndef ProtoDuneChannelRanges_H
#define ProtoDuneChannelRanges_H

#include "art/Utilities/ToolMacros.h"
#include "fhiclcpp/ParameterSet.h"
#include "dune/DuneInterface/Tool/IndexRangeTool.h"
#include <map>

class ProtoDuneChannelRanges : public IndexRangeTool {

public:

  using Name = std::string;
  using Index = IndexRange::Index;
  using IndexRangeMap = std::map<Name, IndexRange>;

  // Ctor.
  ProtoDuneChannelRanges(fhicl::ParameterSet const& ps);

  // Dtor.
  ~ProtoDuneChannelRanges() override =default;

  // Return a range.
  IndexRange get(Name nam) const override;

private:

  // Configuration parameters.
  Index m_LogLevel;
  Name m_ExtraRanges;

  IndexRangeMap m_Ranges;
  const IndexRangeTool* m_pExtraRanges =nullptr;

  // Add an entry to the range map.
  void insertLen(Name nam, Index begin, Index len, Name lab, Name lab1 ="", Name lab2 ="");

};


#endif
