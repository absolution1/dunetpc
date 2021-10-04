// ApaChannelRanges.h
//
// David Adams
// July 2018
//
// Tool to return channel ranges for one or more APAs mapped to TPC sets.
// The TPC sets are numbered 0, 1, ..., NAPA-1.
// The following are created:
//       all - all channels
//      tpsS - for TPC set S, e.g. tps0
//     tppSP - TPC plane or plane pair P in TPC set S, e.g. tps0z
//      apaA - APA, e.g. apa3
//  fembAFFV - FEMB AFF orientation x, eg femb302x
//
// In the above
//    S is the TPC set (offline APA) number
//    A is the APA number
//    P is the wire plane: u, v, z (TPC-side collection) or c (cryostat-side collection)
//      or the wire plane pair x for (z, c) or i for (u, v)
//   FF is the FEMB number in an APA in range [01,20]
//    V is a wire orientation: u, v or x (collection)
//
// Where relevant, the ranges are assigned a second label location
// taken from ApaLocatioNames, e.g US-RaS, ColdBox, FD-UL08,
// and a third label indication the APA, e.g. APA3.
//
// Note there is the option to append another index range tool which can
// oveerride or extend the above set of ranges.
//
// Parameters:
//   LogLevel - Message logging level (0=none, 1=ctor, 2=each call, ...)
//   ApaNumbers - Apa number for each TPC set.
//   ApaLocatioNames - Location name for each TPC set. Padded with "".
//   ExtraRanges - Name of tool with additional ranges. Blank for none.

#ifndef ApaChannelRanges_H
#define ApaChannelRanges_H

#include "art/Utilities/ToolMacros.h"
#include "fhiclcpp/ParameterSet.h"
#include "dune/DuneInterface/Tool/IndexRangeTool.h"
#include <map>

class ApaChannelRanges : public IndexRangeTool {

public:

  using Name = std::string;
  using NameVector = std::vector<Name>;
  using Index = IndexRange::Index;
  using IndexVector = std::vector<Index>;
  using IndexRangeMap = std::map<Name, IndexRange>;

  // Ctor.
  ApaChannelRanges(fhicl::ParameterSet const& ps);

  // Dtor.
  ~ApaChannelRanges() override =default;

  // Return a range.
  IndexRange get(Name nam) const override;

private:

  // Configuration parameters.
  Index m_LogLevel;
  IndexVector m_ApaNumbers;
  NameVector m_ApaLocationNames;
  Name m_ExtraRanges;

  IndexRangeMap m_Ranges;
  const IndexRangeTool* m_pExtraRanges =nullptr;

  // Add an entry to the range map.
  void insertLen(Name nam, Index begin, Index len, Name lab, Name lab1 ="", Name lab2 ="");

};


#endif
