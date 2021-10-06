// ProtoDuneChannelGroups_tool.cc

#include "ProtoDuneChannelGroups.h"
#include "dune/ArtSupport/DuneToolManager.h"
#include "dune/DuneInterface/Tool/IndexRangeTool.h"
#include <iostream>
#include <sstream>
#include <iomanip>

using std::cout;
using std::endl;
using std::string;
using std::ostringstream;
using std::istringstream;
using std::setw;
using std::setfill;
using Name = ProtoDuneChannelGroups::Name;
using NameVector = std::vector<Name>;
using Index = ProtoDuneChannelGroups::Index;
using IndexVector = std::vector<Index>;

//**********************************************************************

ProtoDuneChannelGroups::ProtoDuneChannelGroups(fhicl::ParameterSet const& ps)
: m_LogLevel(ps.get<Index>("LogLevel")),
  m_IndexRangeTool(ps.get<Name>("IndexRangeTool")) {
  const Name myname = "ProtoDuneChannelGroups::ctor: ";
  if ( m_IndexRangeTool.size() ) {
    DuneToolManager* ptm = DuneToolManager::instance();
    m_pIndexRangeTool = ptm->getShared<IndexRangeTool>(m_IndexRangeTool);
    if ( m_pIndexRangeTool == nullptr ) {
      cout << myname << "WARNING: Index range tool not found: " << m_IndexRangeTool << endl;
    }
  } else {
    cout << myname << "WARNING: No Index range tool name." << endl;
  }
  Index napa = 6;
  m_labels["tpss"].push_back("TPC sets");
  m_labels["apas"].push_back("APAs");
  NameVector soris = {"z", "c", "x", "u", "v", "i"};
  for ( Name sori : soris ) {
    m_labels["tpp" + sori + "s"].push_back(sori + " planes");
    m_labels["apa" + sori + "s"].push_back(sori + " planes");
  }
  for ( Index itps=0; itps<napa; ++itps ) {
    ostringstream sstps;
    sstps << itps;
    Name stps = sstps.str();
    Index iapa = itps + 1;
    ostringstream ssapa;
    ssapa << iapa;
    Name sapa = ssapa.str();
    m_groups["tpss"].push_back("tps" + stps);
    m_groups["apas"].push_back("apa" + sapa);
    for ( Name sori : soris ) {
      m_groups["tpp" + sori + "s"].push_back("tpp" + stps + sori);
      m_groups["apa" + sori + "s"].push_back("apa" + sapa + sori);
    }
  }
  Index nfmb = 20;
  for ( Index iapa=1; iapa<=napa; ++iapa ) {
    for ( Index ifmb=1; ifmb<=nfmb; ++ifmb ) {
      ostringstream ssgrp;
      ssgrp << "femb" << iapa;
      if ( ifmb < 10 ) ssgrp << "0";
      ssgrp << ifmb;
      Name sgrp = ssgrp.str();
      m_groups[sgrp].push_back(sgrp + "u");
      m_groups[sgrp].push_back(sgrp + "v");
      m_groups[sgrp].push_back(sgrp + "x");
      Name slab = "FEMB " + sgrp.substr(4);
      m_labels[sgrp].push_back(slab);
    }
  }
  if ( m_LogLevel >= 1 ) {
    cout << myname << "          LogLevel: " << m_LogLevel << endl;
  }
}

//**********************************************************************

IndexRangeGroup ProtoDuneChannelGroups::get(Name nam) const {
  const Name myname = "ProtoDuneChannelGroups::get: ";
  if ( m_pIndexRangeTool == nullptr ) {
    if ( m_LogLevel >= 2 ) cout << myname << "No IndexRangeTool." << endl;
    return IndexRangeGroup();
  }
  GroupMap::const_iterator igrp = m_groups.find(nam);
  if ( igrp == m_groups.end() ) {
    if ( m_LogLevel >= 2 ) cout << myname << "Invalid group name: " << nam << endl;
    return IndexRangeGroup();
  }
  IndexRangeGroup::RangeVector rans;
  for ( Name rnam : igrp->second ) {
    rans.push_back(m_pIndexRangeTool->get(rnam));
  }
  NameVector labs;
  GroupMap::const_iterator ilabs = m_labels.find(nam);
  if ( ilabs == m_labels.end() ) return IndexRangeGroup(nam, rans);
  return IndexRangeGroup(nam, ilabs->second, rans);
}

//**********************************************************************

DEFINE_ART_CLASS_TOOL(ProtoDuneChannelGroups)
