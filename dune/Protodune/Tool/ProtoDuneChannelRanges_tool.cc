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
  Index nchaApa = 2560;
  Index nchau = 800;
  Index nchav = 800;
  Index nchaz = 480;
  bool isEven = true;
  Index apaIdx[6] = { 3, 5, 2, 6, 1, 4 };  // Installation order.
  insertLen("all", 0, ntps*nchaApa, "All");
  for ( Index itps=0; itps<ntps; ++itps ) {
    string siapa = std::to_string(apaIdx[itps]);
    Index ch0 = itps*nchaApa;
    string sitps = std::to_string(itps);
    string stps = "tps" + sitps;
    string sapa = "apa" + siapa;
    insertLen(      stps, ch0, nchaApa, "TPC set " + sitps);
    insertLen(      sapa, ch0, nchaApa, "APA " + siapa);
    string stpp = "tpp" + sitps;
    insertLen(stpp + "u", ch0, nchau, "TPC plane " + sitps + "u");
    insertLen(sapa + "u", ch0, nchau, "APA plane " + siapa + "u");
    ch0 += nchau;
    insertLen(stpp + "v", ch0, nchav, "TPC plane " + sitps + "v");
    insertLen(sapa + "v", ch0, nchav, "APA plane " + siapa + "v");
    ch0 += nchav;
    Index chx1 = ch0;
    ch0 += nchaz;
    Index chx2 = ch0;
    Index chz = isEven ? chx2 : chx1;
    Index chc = isEven ? chx1 : chx2;
    insertLen(stpp + "c", chc, nchaz, "TPC plane " + sitps + "c");
    insertLen(sapa + "c", chc, nchaz, "APA plane " + siapa + "c");
    insertLen(stpp + "z", chz, nchaz, "TPC plane " + sitps + "z");
    insertLen(sapa + "z", chz, nchaz, "APA plane " + siapa + "z");
    isEven = ! isEven;
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
