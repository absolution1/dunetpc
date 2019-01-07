// FclStickyCodeFlagger_tool.cc

#include "FclStickyCodeFlagger.h"
#include <iostream>
#include <sstream>

using std::string;
using std::cout;
using std::endl;
using std::istringstream;
using Index = FclStickyCodeFlagger::Index;
using IndexVector = FclStickyCodeFlagger::IndexVector;
using Name = FclStickyCodeFlagger::Name;
using NameVector = std::vector<Name>;

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
  // Parse the sticky code map.
  fhicl::ParameterSet pstab = ps.get<fhicl::ParameterSet>("StickyCodes");
  NameVector chnams = pstab.get_names();
  Index ncha = 0;
  Index ncod = 0;
  if ( m_LogLevel >= 3 ) cout << myname << "Sticky code channel count: " << chnams.size() << endl;
  if ( chnams.size() == 0 ) cout << myname << "WARNING: Sticky code channel map is empty." << endl;
  for ( Name chnam : chnams ) {
    if ( m_LogLevel >= 3 ) cout << myname << "Fetching codes for " << chnam << endl;
    bool badnam = true;
    Name msg = "ERROR: Invalid channel specifier: ";
    if ( chnam.substr(0,4) == "chan" ) {
      istringstream sscha(chnam.substr(4));
      const Index badIndex = -1;
      Index icha = badIndex;
      sscha >> icha;
      if ( icha == badIndex ) {
        msg = "WARNING: Ignoring duplicate channel specifier: ";
      } else {
        m_stickyCodes[icha] = pstab.get<IndexVector>(chnam);
        badnam = false;
        ++ncha;
        ncod += m_stickyCodes[icha].size();
      }
    }
    if ( badnam ) cout << myname << msg << chnam << endl;
    else if ( m_LogLevel >= 2 ) cout << myname << "Read bad codes for " << chnam << endl;
  }
  cout << myname << "Found " << ncod << " sticky code" << (ncod == 1 ? "" : "s") << " for "
       << ncha << " channel" << (ncha == 1 ? "" : "s") << "." << endl;
}

//**********************************************************************

DataMap FclStickyCodeFlagger::update(AdcChannelData& acd) const {
  const string myname = "FclStickyCodeFlagger::view: ";
  DataMap ret;
  IndexMap::const_iterator icod = m_stickyCodes.find(acd.channel);
  if ( icod == m_stickyCodes.end() ) return ret;
  const IndexVector& codes = icod->second;
  if ( acd.flags.size() < acd.raw.size() ) {
    cout << myname << "WARNING: Increasing size of the flags vector for channel " << acd.channel << endl;
    acd.flags.resize(acd.raw.size(), AdcGood);
  }
  Index stickyCodeCount = 0;
  for ( Index isam=0; isam<acd.raw.size(); ++isam ) {
    if ( std::find(codes.begin(), codes.end(), acd.raw[isam]) != codes.end() ) {
      acd.flags[isam] = m_StickyCode;
      ++stickyCodeCount;
    }
  }
  ret.setInt("stickyChannel", acd.channel);
  DataMap::IntVector badCodes;
  for ( Index badCode : codes ) badCodes.push_back(badCode);
  ret.setIntVector("stickyCodes", badCodes);
  ret.setInt("stickyCodeCount", stickyCodeCount);
  return ret;
}

//**********************************************************************
