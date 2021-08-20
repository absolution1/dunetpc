// HalfApaChannelRanges_tool.cc

#include "HalfApaChannelRanges.h"
#include "dune/ArtSupport/DuneToolManager.h"
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
using Name = HalfApaChannelRanges::Name;
using Index = HalfApaChannelRanges::Index;
using IndexVector = std::vector<Index>;

//**********************************************************************

HalfApaChannelRanges::HalfApaChannelRanges(fhicl::ParameterSet const& ps)
: m_LogLevel(ps.get<Index>("LogLevel")),
  m_ApaNumbers(ps.get<IndexVector>("ApaNumbers")),
  m_ExtraRanges(ps.get<Name>("ExtraRanges")) {
  const Name myname = "HalfApaChannelRanges::ctor: ";
  const Index ntps = m_ApaNumbers.size();
  Index nchaApa = 1280;
  Index nfchau = 40;
  Index nfchav = 40;
  Index nfchax = 48;
  Index nchau = 400;
  Index nchav = 400;
  Index nchaz = 240;
  Index nchax = 2*nchaz;
  Index nchai = nchau + nchav;
  const IndexVector& apaIdx = m_ApaNumbers;
  insertLen("all", 0, ntps*nchaApa, "All", "");
  for ( Index itps=0; itps<ntps; ++itps ) {
    string sitps;
    string siapa;
    if ( ntps > 1 ) {
      sitps = std::to_string(itps);
      siapa = std::to_string(apaIdx[itps]);
    }
    string labTps = "TPC set " + sitps;
    string labApa = "APA " + siapa;
    string sloc = "";
    Index ch0 = itps*nchaApa;
    string stps = "tps" + sitps;
    string sapa = "apa" + siapa;
    insertLen(      stps, ch0, nchaApa, labTps, labApa);
    insertLen(      sapa, ch0, nchaApa, labApa);
    string stpp = "tpp" + sitps;
    Index chu0 = ch0;
    Index chv0 = chu0 + nchau;
    Index chx10 = chv0 + nchav;
    Index chx20 = chx10 + nchaz;
    Index chx0 = chx10;
    bool beamRight = 2*(itps/2) == itps;
    Index chz0 = beamRight ? chx20 : chx10;
    Index chc0 = beamRight ? chx10 : chx20;
    insertLen(stpp + "u", chu0, nchau, "TPC plane " + sitps + "u", labApa);
    insertLen(sapa + "u", chu0, nchau, "APA plane " + siapa + "u");
    insertLen(stpp + "v", chv0, nchav, "TPC plane " + sitps + "v", labApa);
    insertLen(sapa + "v", chv0, nchav, "APA plane " + siapa + "v");
    insertLen(stpp + "z1", chc0, nchaz, "TPC plane " + sitps + "z1", labApa);
    insertLen(sapa + "z1", chc0, nchaz, "APA plane " + siapa + "z1");
    insertLen(stpp + "z2", chz0, nchaz, "TPC plane " + sitps + "z2", labApa);
    insertLen(sapa + "z2", chz0, nchaz, "APA plane " + siapa + "z2");
    insertLen(stpp + "i", chu0, nchai, "TPC planes " + sitps + "i", labApa);
    insertLen(sapa + "i", chu0, nchai, "APA planes " + siapa + "i");
    insertLen(stpp + "x", chx0, nchax, "TPC planes " + sitps + "x", labApa);
    insertLen(sapa + "x", chx0, nchax, "APA planes " + siapa + "x");
    insertLen(stpp + "z", chx0, nchax, "collection planes " + sitps, labApa);
    insertLen(sapa + "z", chx0, nchax, "collection planes " + siapa);
    Index fchu0 = chu0;
    Index fchv0 = chv0;
    Index fchx0 = chx10;
    Index ifmbu =  6;
    Index ifmbv = 10;
    Index ifmbx = 10;
    // Loop over FEMBS in offline order.
    for ( Index ifmbOff=0; ifmbOff<20; ++ifmbOff ) {
      ostringstream ssnamu;
      ssnamu << siapa << setfill('0') << setw(2) << ifmbu << "u";
      string namu = ssnamu.str();
      insertLen("femb" + namu, fchu0, nfchau, "FEMB block " + namu);
      ostringstream ssnamv;
      ssnamv << siapa << setfill('0') << setw(2) << ifmbv << "v";
      string namv = ssnamv.str();
      insertLen("femb" + namv, fchv0, nfchav, "FEMB block " + namv);
      ostringstream ssnamx;
      ssnamx << siapa << setfill('0') << setw(2) << ifmbx << "x";
      string namx = ssnamx.str();
      insertLen("femb" + namx, fchx0, nfchax, "FEMB block " + namx);
      fchu0 += nfchau;
      fchv0 += nfchav;
      fchx0 += nfchax;
      ifmbu += 1; 
      if ( ifmbu > 10 ) ifmbu = 1;
      ifmbv -= 1;
      if ( ifmbv == 0 ) ifmbv = 10;
      if ( ifmbOff < 4 )       ifmbx -= 1;
      else if ( ifmbOff == 4 ) ifmbx = 1;
      else                     ifmbx += 1;
      
    }
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

IndexRange HalfApaChannelRanges::get(Name nam) const {
  const Name myname = "HalfApaChannelRanges::runData: ";
  if ( m_pExtraRanges != nullptr ) {
    IndexRange rout = m_pExtraRanges->get(nam);
    if ( rout.isValid() ) return rout;
  }
  // Special handling for online wire specifier fembAFFVCC
  // We look up FEMB block fembAFFV and pick the channel corresponding to CC,
  if ( nam.size() == 10 && nam.substr(0,4) == "femb" ) {
    IndexRange fbran = get(nam.substr(0,8));
    Index ifchmax = 0;
    if ( fbran.isValid() ) {
      char cpla = nam[7];
      bool dirSame = true;   // Is channel numbering dir same in FEMB and offline?
      if ( cpla == 'u' ) {
        ifchmax = 40;
        dirSame = false;
      } else if ( cpla == 'v' ) {
        ifchmax = 40;
      } else if ( cpla == 'x' || cpla == 'w') {
        ifchmax = 48;
        istringstream ssapa(nam.substr(4,1));
        int iapa = 0;
        ssapa >> iapa;
        istringstream ssfmb(nam.substr(5,2));
        int ifmb = 0;
        ssfmb >> ifmb;
        bool beamL = ifmb >=  1 && ifmb <= 10;
        bool beamR = ifmb >= 11 && ifmb <= 20;
        if ( beamL )      dirSame = false;
        else if ( beamR ) dirSame = true;
        else              ifchmax = 0;
      }
      if ( ifchmax != fbran.size() ) {
        cout << myname << "WARNING: FEMB block has unexpected size: "
             << ifchmax << " != " << fbran.size() << endl;
        ifchmax = 0;
      }
      string sch = nam.substr(8,2);
      istringstream ssfch(sch);
      Index ifch;
      ssfch >> ifch;
      if ( ifch > 0 && ifch <= ifchmax ) {
        if ( dirSame ) --ifch;
        else ifch = ifchmax - ifch;
        Index ich = fbran.begin + ifch;
        return IndexRange(nam, ich, ich+1, fbran.label(0) + sch, fbran.label(1));
      } 
    }
  }
  IndexRangeMap::const_iterator iran = m_Ranges.find(nam);
  if ( iran == m_Ranges.end() ) return IndexRange();
  return iran->second;
}

//**********************************************************************

void HalfApaChannelRanges::
insertLen(Name nam, Index begin, Index len, Name lab, Name lab1, Name lab2) {
  m_Ranges[nam] = IndexRange(nam, begin, begin+len, lab, lab1, lab2);
}

//**********************************************************************

DEFINE_ART_CLASS_TOOL(HalfApaChannelRanges)
