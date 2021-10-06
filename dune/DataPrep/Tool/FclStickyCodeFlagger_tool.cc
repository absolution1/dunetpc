// FclStickyCodeFlagger_tool.cc

#include "FclStickyCodeFlagger.h"
#include <iostream>
#include <sstream>
#include <set>

using std::string;
using std::cout;
using std::endl;
using std::istringstream;
using Index = FclStickyCodeFlagger::Index;
using IndexVector = FclStickyCodeFlagger::IndexVector;
using IndexSet = std::set<Index>;
using Name = FclStickyCodeFlagger::Name;
using NameVector = std::vector<Name>;
using IndexVectorMultiMap = std::multimap<Index, IndexVector>;
using IndexPair = FclStickyCodeFlagger::IndexPair;
using IndexPairMultiMap = FclStickyCodeFlagger::IndexPairMap;

//**********************************************************************
// Local definitions.
//**********************************************************************

// Read from a fcl map of arbitrary fcl-supported type indexed by unsigned int 
// encoded in fcl name into a multimap with the same index and value type T.
// Leading zeroes are stripped from the fcl index.
// Fcl:  preIII: TT, ...
// to map<string,T> m[III] = TT, ...
template<class T>
int readFclIndexMap(string pre, string sdesc, int logLevel, const fhicl::ParameterSet& pstab, std::multimap<Index,T>& dat) {
  using MMap = std::multimap<Index,T>;
  const string myname = "readFclIndexMap: ";
  Index lpre = pre.size();
  Index nerr = 0;
  NameVector chnams = pstab.get_names();
  if ( logLevel >= 3 ) cout << myname << sdesc << " fcl entry count: " << chnams.size() << endl;
  for ( Name chnam : chnams ) {
    if ( logLevel >= 3 ) cout << myname << sdesc << " fetching for " << chnam << endl;
    bool badnam = true;
    string::size_type ipos = lpre;
    string::size_type jpos = chnam.size() - 1;
    Name msg = "ERROR: " + sdesc + " invalid index specifier: ";
    if ( chnam.substr(0,ipos) == pre ) {
      // Remove leading zeroes to be extra careful. Not needed.
      while ( chnam[ipos] == '0' && ipos+1 < chnam.size() ) ++ipos;
      while ( jpos > ipos && chnam[jpos] == 'x' ) --jpos;
      const Index badIndex = -1;
      string scha = chnam.substr(ipos, jpos+1-ipos);
      Index icha = badIndex;
      istringstream sscha(scha);
      sscha >> icha;
      if ( dat.count(icha) ) {
        cout << myname << "WARNING: " + sdesc + " replacing data for index " << icha << endl;
      }
      T val = pstab.get<T>(chnam);
      dat.insert(typename MMap::value_type(icha, val));
      badnam = false;
    } else {
      ++nerr;
    }
    if ( badnam ) cout << myname << msg << chnam << endl;
  }
  return nerr;
}

//**********************************************************************
// Class methods.
//**********************************************************************

FclStickyCodeFlagger::FclStickyCodeFlagger(fhicl::ParameterSet const& ps)
: m_LogLevel(ps.get<int>("LogLevel")),
  m_StickyCode(ps.get<AdcFlag>("StickyCode")) {
  const string myname = "FclStickyCodeFlagger::ctor: ";
  // Display configuration.
  if ( m_LogLevel ) {
    cout << myname << "Configuration: " << endl;
    cout << myname << "          LogLevel: " << m_LogLevel << endl;
    cout << myname << "        StickyCode: " << m_StickyCode << endl;
  }
  // Check the sticky code.
  if ( m_StickyCode < AdcStuck || m_StickyCode >= AdcMitigated ) {
    cout << myname << "WARNING: Flag for sticky code has an unexpected value: "
         << m_StickyCode << endl;
  }
  // Build the sticky code map.
  fhicl::ParameterSet pstabsc = ps.get<fhicl::ParameterSet>("StickyCodes");
  IndexVectorMultiMap mvmap;
  int nerr = readFclIndexMap("chan", "Sticky codes", m_LogLevel, pstabsc, mvmap);
  for ( const IndexVectorMultiMap::value_type& imma : mvmap ) {
    Index icha = imma.first;
    const IndexVector& newcodes = imma.second;
    IndexVector& codes = m_stickyCodes[icha];
    codes.insert(codes.end(), newcodes.begin(), newcodes.end());
  }
  Index ncha = m_stickyCodes.size();
  Index ncod = 0;
  for ( IndexVectorMap::value_type& ient : m_stickyCodes ) ncod += ient.second.size();
  cout << myname << "Found " << ncod << " sticky code" << (ncod == 1 ? "" : "s") << " for "
       << ncha << " channel" << (ncha == 1 ? "" : "s") << "." << endl;
  if ( nerr ) {
    cout << myname << "WARNING: Found " << nerr << " error" << (nerr == 1 ? "" : "s")
                   << " while parsing the sticky code map." << endl;
  }
  // Build the sticky code range map.
  fhicl::ParameterSet pstabsr = ps.get<fhicl::ParameterSet>("StickyRanges");
  nerr = readFclIndexMap("chan", "Sticky ranges", m_LogLevel, pstabsr, m_stickyRanges);
  IndexSet chset;
  Index nran = m_stickyRanges.size();
  for ( IndexPairMap::value_type& ient : m_stickyRanges ) chset.insert(ient.first);
  ncha = chset.size();
  cout << myname << "Found " << nran << " sticky range" << (nran == 1 ? "" : "s") << " for "
       << ncha << " channel" << (ncha == 1 ? "" : "s") << "." << endl;
  if ( nerr ) {
    cout << myname << "WARNING: Found " << nerr << " error" << (nerr == 1 ? "" : "s")
                   << " while parsing the sticky range map." << endl;
  }
}

//**********************************************************************

DataMap FclStickyCodeFlagger::update(AdcChannelData& acd) const {
  const string myname = "FclStickyCodeFlagger::view: ";
  DataMap ret;
  // Make sure flags is as long as raw.
  if ( acd.flags.size() < acd.raw.size() ) {
    cout << myname << "WARNING: Increasing size of the flags vector for channel " << acd.channel() << endl;
    acd.flags.resize(acd.raw.size(), AdcGood);
  }
  IndexSet samplesToFlag;
  // Find samples with sticky codes.
  IndexVectorMap::const_iterator icod = m_stickyCodes.find(acd.channel());
  if ( icod != m_stickyCodes.end() ) {
    const IndexVector& codes = icod->second;
    for ( Index isam=0; isam<acd.raw.size(); ++isam ) {
      if ( std::find(codes.begin(), codes.end(), acd.raw[isam]) != codes.end() ) samplesToFlag.insert(isam);
    }
  }
  // Find samples with in sticky code ranges.
  using IndexPairRange = std::pair<IndexPairMap::const_iterator, IndexPairMap::const_iterator>;
  IndexPairRange irans = m_stickyRanges.equal_range(acd.channel());
  for ( IndexPairMap::const_iterator iran=irans.first; iran!=irans.second; ++iran ) {
    IndexPair ran = iran->second;
    for ( Index isam=0; isam<acd.raw.size(); ++isam ) {
      AdcIndex adc = acd.raw[isam];
      if ( adc >= ran.first && adc <= ran.second ) samplesToFlag.insert(isam);
    }
  }
  // Set flags.
  for ( AdcIndex isam : samplesToFlag ) acd.flags[isam] = m_StickyCode;
  // Fill output data map and return.
  ret.setInt("stickyChannel", acd.channel());
  ret.setInt("stickyCodeCount", samplesToFlag.size());
  return ret;
}

//**********************************************************************

DEFINE_ART_CLASS_TOOL(FclStickyCodeFlagger)
