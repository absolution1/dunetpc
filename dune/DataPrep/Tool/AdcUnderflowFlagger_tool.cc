// AdcUnderflowFlagger_tool.cc

#include "AdcUnderflowFlagger.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>

using std::string;
using std::cout;
using std::endl;
using std::ofstream;
using std::ostream;
using std::ostringstream;
using std::setw;
using std::fixed;
using std::setprecision;
using std::vector;

using Index = unsigned int;

//**********************************************************************
// Class methods.
//**********************************************************************

AdcUnderflowFlagger::AdcUnderflowFlagger(fhicl::ParameterSet const& ps)
: m_LogLevel(ps.get<int>("LogLevel")),
  m_AdcThresholds(ps.get<AdcCountVector>("AdcThresholds")),
  m_DefaultThreshold(ps.get<int>("DefaultThreshold")),
  m_modifyCount(0) {
  const string myname = "AdcUnderflowFlagger::ctor: ";
  if ( m_LogLevel >= 1 ) {
    cout << myname << "Configuration parameters:" << endl;
    cout << myname << "      LogLevel: " << m_LogLevel << endl;
    cout << myname << " AdcThresholds: [";
    Index maxcha = 10;
    Index ncha = m_AdcThresholds.size();
    bool haveMore = ncha > maxcha;
    if ( ! haveMore ) maxcha = ncha;
    for ( Index icha=0; icha<maxcha; ++icha ) {
      if ( icha ) cout << ",";
      cout << " " << m_AdcThresholds[icha];
    }
    if ( haveMore ) cout << ",...";
    cout << "]" << endl;
  }
}

//**********************************************************************

DataMap AdcUnderflowFlagger::view(const AdcChannelData& acd) const {
  const string myname = "AdcUnderflowFlagger::view: ";
  AdcChannelData acdtmp;
  acdtmp.channel = acd.channel;
  acdtmp.raw = acd.raw;
  return update(acdtmp);
}

//**********************************************************************

DataMap AdcUnderflowFlagger::update(AdcChannelData& acd) const {
  const string myname = "AdcUnderflowFlagger::update: ";
  if ( m_LogLevel >= 3 ) cout << myname << "Calling " << endl;
  DataMap res(0);
  AdcChannel icha = acd.channel;
  res.setInt("channel", icha);
  res.setInt("nModify", m_modifyCount);
  if ( icha == AdcChannelData::badChannel ) {
    if ( m_LogLevel >= 2 ) cout << myname << "Data does not have an assigned channel." << endl;
    return res.setStatus(1);
  }
  IntVector::value_type sthr = -1;
  Index ncha = m_AdcThresholds.size();
  if ( icha < ncha ) sthr = m_AdcThresholds[icha];
  else sthr = m_DefaultThreshold;
  if ( sthr < 0 ) {
    if ( m_LogLevel >= 2 ) cout << myname << "No threhold set for channel " << icha << endl;
    return res.setStatus(2 + icha >= ncha);
  }
  AdcCount thr = sthr;
  Index nraw = acd.raw.size();
  Index nflg = acd.flags.size();
  if ( nflg < nraw ) acd.flags.resize(nraw, AdcGood);
  Index nunder = 0;
  for ( Index isam=0; isam<nraw; ++isam ) {
    if ( acd.raw[isam] <= thr ) {
      ++nunder;
      acd.flags[isam] = AdcUnderflow;
    }
  }
  res.setInt("nUnderflow", nunder);
  return res;
}

//**********************************************************************
