// AdcChannelTrimmer_tool.cc

#include "AdcChannelTrimmer.h"
#include <iostream>
#include <vector>

using std::string;
using std::cout;
using std::endl;

//**********************************************************************
// Class methods.
//**********************************************************************

AdcChannelTrimmer::AdcChannelTrimmer(fhicl::ParameterSet const& ps)
: m_LogLevel(ps.get<int>("LogLevel")), 
  m_Length(ps.get<Index>("Length")), 
  m_MaxTrim(ps.get<Index>("MaxTrim")) {
  const string myname = "AdcChannelTrimmer::ctor: ";
  // Display the configuration.
  if ( m_LogLevel >= 1 ) {
    cout << myname << "Configuration: " << endl;
    cout << myname << "         LogLevel: " << m_LogLevel << endl;
    cout << myname << "           Length: " << m_Length << endl;
    cout << myname << "          MaxTrim: " << m_MaxTrim << endl;
  }
}

//**********************************************************************

DataMap AdcChannelTrimmer::view(const AdcChannelData& acd) const {
  const string myname = "AdcChannelTrimmer::view: ";
  DataMap ret;
  Index nsam = acd.samples.size();
  if ( nsam == m_Length ) {
    ret.setInt("trimAction", 0);
    if ( m_LogLevel >= 3 ) cout << "Samples have requested length of "
                                << m_Length << "." << endl;
    return ret;
  }
  if ( nsam == 0 ) {
    ret.setInt("trimAction", 3);
    if ( m_LogLevel >= 3 ) cout << "WARNING: No samples found." << endl;
    return ret;
  }
  int idsam = nsam < m_Length ? -int(m_Length - nsam) : int(nsam - m_Length);
  Index udsam = idsam < 0 ? -idsam : idsam;
  ret.setInt("trimLength", idsam);
  if ( m_MaxTrim && udsam > m_MaxTrim ) {
    ret.setInt("trimAction", 2);
    cout << myname << "ERROR: Channel " << acd.channel() << " sample length " << nsam
         << " is outside of trim range " << m_MaxTrim  << "." << endl;
  } else {
    ret.setInt("trimAction", -1);
    if ( m_LogLevel >= 2 ) cout << myname << "Trimming channel " << acd.channel() << " from "
                                << nsam << " to " << m_Length << " samples." << endl;
  }
  return ret;
}

//**********************************************************************

DataMap AdcChannelTrimmer::update(AdcChannelData& acd) const {
  const string myname = "AdcChannelTrimmer::update: ";
  DataMap ret = view(acd);
  int action = ret.getInt("trimAction");
  if ( action != -1 ) return ret;
  Index nsam = acd.samples.size();
  acd.samples.resize(m_Length);
  Index jsam = 0;
  for ( Index isam=nsam; isam<m_Length; ++isam ) {
    acd.samples[isam] = acd.samples[jsam++];
  }
  ret.setInt("trimAction", 1);
  return ret;
}

//**********************************************************************

DEFINE_ART_CLASS_TOOL(AdcChannelTrimmer)
