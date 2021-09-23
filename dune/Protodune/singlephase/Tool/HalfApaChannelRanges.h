// HalfApaChannelRanges.h
//
// David Adams
// March 2020
//
// Tool to return channel ranges for one or more "half" APAs mapped to TPC sets.
// A half APA has 1280 channels: 400 u, 400 v and 240 z1 and 240 z2.
// Iceberg has one such APA.
//
// The TPC sets are numbered 0, 1, ..., NAPA-1.
// For inly one APA, the numbers are omitted.
// The following are created:
//       all - all channels
//      tpsS - for TPC set S, e.g. tps0
//     tppSP - TPC plane or plane pair P in TPC set S, e.g. tps0z
//      apaA - APA A -- one-to-one mapping to tps
//  fembAFFV - FEMB AFF orientation x, eg femb302x
//
// In the above
//  S=A is the TPC set (APA) number
//    A is the APA number
//    P is the wire plane: u, v, z (TPC-side collection) or c (cryostat-side collection)
//      or the wire plane pair x for (z, c) or i for (u, v)
//   FF is the FEMB number in an APA in range [01,20]
//    V is a wire orientation: u, v or x (collection)
//
// Note there is the option to append another index range tool which can
// oveerride or extend the above set of ranges.
//
// Parameters:
//   LogLevel - Message logging level (0=none, 1=ctor, 2=each call, ...)
//   ApaNumbers - Apa number for each TPC set. # entries is the # APAs.
//   ExtraRanges - Name of tool with additional ranges. Blank for none.

#ifndef ApaChannelRanges_H
#define ApaChannelRanges_H

#include "art/Utilities/ToolMacros.h"
#include "fhiclcpp/ParameterSet.h"
#include "dune/DuneInterface/Tool/IndexRangeTool.h"
#include <map>

class HalfApaChannelRanges : public IndexRangeTool {

public:

  using Name = std::string;
  using NameVector = std::vector<Name>;
  using Index = IndexRange::Index;
  using IndexVector = std::vector<Index>;
  using IndexRangeMap = std::map<Name, IndexRange>;

  // Ctor.
  HalfApaChannelRanges(fhicl::ParameterSet const& ps);

  // Dtor.
  ~HalfApaChannelRanges() override =default;

  // Return a range.
  IndexRange get(Name nam) const override;

private:

  // Configuration parameters.
  Index m_LogLevel;
  IndexVector m_ApaNumbers;
  Name m_ExtraRanges;

  IndexRangeMap m_Ranges;
  const IndexRangeTool* m_pExtraRanges =nullptr;

  // Add an entry to the range map.
  void insertLen(Name nam, Index begin, Index len, Name lab, Name lab1 ="", Name lab2 ="");

};


#endif
