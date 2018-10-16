// AdcCodeMitigator_tool.cc

#include "AdcCodeMitigator.h"
#include <iostream>
#include <sstream>

using std::string;
using std::cout;
using std::endl;
using std::istringstream;
using Index = AdcCodeMitigator::Index;
using IndexVector = AdcCodeMitigator::IndexVector;
using Name = AdcCodeMitigator::Name;
using NameVector = std::vector<Name>;

//**********************************************************************
// Class methods.
//**********************************************************************

AdcCodeMitigator::AdcCodeMitigator(fhicl::ParameterSet const& ps)
: m_LogLevel(ps.get<int>("LogLevel")),
  m_FixFlags(ps.get<IndexVector>("FixFlags")),
  m_InterpolateFlags(ps.get<IndexVector>("InterpolateFlags")),
  m_SkipFlags(ps.get<IndexVector>("SkipFlags")) {
  const string myname = "AdcCodeMitigator::ctor: ";
  // Display configuration.
  if ( m_LogLevel ) {
    cout << myname << "Configuration: " << endl;
    cout << myname << "          LogLevel: " << m_LogLevel << endl;
    cout << myname << "          FixFlags: [";
    bool first = true;
    for ( Index iflg : m_FixFlags ) {
      if ( ! first ) cout << ", ";
      cout << iflg;
    }
    cout << "]" << endl;
    cout << myname << "  InterpolateFlags: [";
    first = true;
    for ( Index iflg : m_InterpolateFlags ) {
      if ( ! first ) cout << ", ";
      cout << iflg;
    }
    cout << "]" << endl;
    cout << myname << "         SkipFlags: [";
    first = true;
    for ( Index iflg : m_SkipFlags ) {
      if ( ! first ) cout << ", ";
      cout << iflg;
    }
    cout << "]" << endl;
  }
  for ( Index iflg : m_FixFlags ) m_fixSet.insert(iflg);
  for ( Index iflg : m_InterpolateFlags ) {
    m_interpolateSet.insert(iflg);
    m_skipSet.insert(iflg);
  }
  for ( Index iflg : m_SkipFlags ) m_skipSet.insert(iflg);
}

//**********************************************************************

DataMap AdcCodeMitigator::update(AdcChannelData& acd) const {
  const string myname = "AdcCodeMitigator::view: ";
  DataMap ret;
  Index nsam = acd.samples.size();
  if ( acd.flags.size() != acd.samples.size() ) {
    cout << myname << "ERROR: Flag and sample sizes disagree: "
         << acd.flags.size() << " != " << nsam << endl;
    return ret.setStatus(1);
  }
  Index mitCount = 0;
  // Fixed mitigation.
  for ( Index isam=0; isam<nsam; ++isam ) {
    AdcFlag iflg = acd.flags[isam];
    if ( m_fixSet.find(iflg) == m_fixSet.end() ) continue;
    ++mitCount;
    acd.samples[isam] = 0.0;
    acd.flags[isam] = AdcSetFixed;
  }
  // Interpolation mitigation.
  Index isamLo = 0;  // Sample used for high-side interpolation.
  bool haveLo = false;
  Index isamHi = 0;  // Sample used for high-side interpolation.
  for ( Index isam=0; isam<nsam; ++isam ) {
    AdcFlag iflg = acd.flags[isam];
    // Record if this sample can be used for low-side interpolation.
    if ( m_skipSet.find(iflg) == m_skipSet.end() ) {
      isamLo = isam;
      haveLo = true;
      continue;
    }
    // Skip sample if it is not to be interpolated.
    if ( m_interpolateSet.find(iflg) == m_interpolateSet.end() ) continue;
    // Check if we have a valid high-side sample.
    bool haveHi = isamHi > isam && isamHi < nsam;
    if ( ! haveHi && isamHi < nsam ) {
      isamHi = isam;
      while ( !haveHi && ++isamHi < nsam ) haveHi = m_skipSet.find(acd.flags[isamHi]) == m_skipSet.end();
    }
    double y = 0.0;
    if ( haveLo && haveHi ) {
      double ylo = acd.samples[isamLo];
      double yhi = acd.samples[isamHi];
      double x01 = isam - isamLo;
      double x12 = isamHi - isamLo;
      y = ylo + (yhi - ylo)*x01/x12;
      iflg = AdcInterpolated;
      if ( m_LogLevel >= 3 ) cout << "Interpolating isam=" << isam
                                  << ": [" << isamLo << "], [" << isamHi << "] ==> " << y << endl;
    } else if ( haveLo ) {
      y = acd.samples[isamLo];
      iflg = AdcExtrapolated;
      if ( m_LogLevel >= 3 ) cout << "Setting low extrapolation for isam=" << isam
                                  << ": [" << isamLo << "] = " << y << endl;
    } else if ( haveHi ) {
      y = acd.samples[isamHi];
      iflg = AdcExtrapolated;
      if ( m_LogLevel >= 3 ) cout << "Setting high extrapolation for isam=" << isam
                                  << ": [" << isamHi << "] = " << y << endl;
    } else {
      iflg = AdcSetFixed;
    }
    acd.samples[isam] = y;
    acd.flags[isam] = iflg;
    ++mitCount;
  }
  ret.setInt("mitCount", mitCount);
  return ret;
}

//**********************************************************************
