// ProtoDuneChannelRanges_tool.cc

#include "ProtoDuneChannelRanges.h"
#include "dune/ArtSupport/DuneToolManager.h"
#include <iostream>
#include <iomanip>

using std::cout;
using std::endl;
using std::string;
using std::setw;
using Name = ProtoDuneChannelRanges::Name;
using Index = ProtoDuneChannelRanges::Index;
using IndexVector = std::vector<Index>;

//**********************************************************************

ProtoDuneChannelRanges::ProtoDuneChannelRanges(fhicl::ParameterSet const& ps)
: m_LogLevel(ps.get<Index>("LogLevel")),
  m_ExtraRanges(ps.get<Name>("ExtraRanges")) {
  const Name myname = "ProtoDuneChannelRanges::ctor: ";
  const Index ntps = 6;
  Index nwirApa = 2560;
  Index nwiru = 800;
  Index nwirv = 800;
  Index nwirz = 480;
  for ( Index itps=0; itps<ntps; ++itps ) {
    Index ch0 = itps*nwirApa;
    string sitps = std::to_string(itps);
    string stps = "tps" + sitps;
    insertLen(      stps, ch0, nwirApa, "TPC set " + sitps);
    string stpp = "tpp" + sitps;
    insertLen(stpp + "u", ch0, nwiru, "TPC plane " + sitps + "u");
    ch0 += nwiru;
    insertLen(stpp + "v", ch0, nwirv, "TPC plane " + sitps + "v");
    ch0 += nwirv;
    insertLen(stpp + "z", ch0, nwirz, "TPC plane " + sitps + "z");
  }
  if ( m_ExtraRanges.size() ) {
    DuneToolManager* ptm = DuneToolManager::instance();
    m_pExtraRanges = ptm->getShared<IndexRangeTool>(m_ExtraRanges);
    if ( m_pExtraRanges == nullptr ) {
      cout << myname << "WARNING: Extra range tool not found: " << m_ExtraRanges << endl;
    }
  }
  if ( m_LogLevel >= 1 ) {
    cout << myname << "          LogLevel: " << m_LogLevel << endl;
    cout << myname << "  Extra range tool: " << m_ExtraRanges << endl;
  }
}

//**********************************************************************

IndexRange ProtoDuneChannelRanges::get(Name nam) const {
  const Name myname = "ProtoDuneChannelRanges::runData: ";
  if ( m_pExtraRanges != nullptr ) {
    IndexRange rout = m_pExtraRanges->get(nam);
    if ( rout.isValid() ) return rout;
  }
  IndexRangeMap::const_iterator iran = m_Ranges.find(nam);
  if ( iran == m_Ranges.end() ) return IndexRange();
  return iran->second;
}

//**********************************************************************

void ProtoDuneChannelRanges::
insertLen(Name nam, Index begin, Index len, Name lab) {
  m_Ranges[nam] = IndexRange(nam, begin, begin+len, lab);
}

//**********************************************************************
