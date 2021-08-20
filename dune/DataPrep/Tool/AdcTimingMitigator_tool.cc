// AdcTimingMitigator_tool.cc

#include "AdcTimingMitigator.h"
#include <iostream>

using std::string;
using std::cout;
using std::endl;

using Index = unsigned int;

//**********************************************************************
// Class methods.
//**********************************************************************

AdcTimingMitigator::AdcTimingMitigator(fhicl::ParameterSet const& ps)
: m_LogLevel(ps.get<int>("LogLevel")),
  m_SamplingRatio(ps.get<float>("SamplingRatio")),
  m_SamplingOffset(ps.get<float>("SamplingOffset")),
  m_FEMBs(ps.get<IndexVector>("FEMBs")) {
  const string myname = "AdcTimingMitigator::ctor: ";
  for ( Index ifmb : m_FEMBs ) m_fembSet.insert(ifmb);
  if ( m_LogLevel >= 1 ) {
    cout << myname << "Configuration parameters:" << endl;
    cout << myname << "        LogLevel: " << m_LogLevel << endl;
    cout << myname << "   SamplingRatio: " << m_SamplingRatio << endl;
    cout << myname << "  SamplingOffset: " << m_SamplingOffset << endl;
    cout << myname << "           FEMBs: [";
    bool first = true;
    for ( Index ifmb : m_fembSet ) {
      if ( first ) first = false;
      else cout << ", ";
      cout << ifmb;
    }
    cout << "]" << endl;
  }
}

//**********************************************************************

DataMap AdcTimingMitigator::update(AdcChannelData& acd) const {
  const string myname = "AdcTimingMitigator::update: ";
  DataMap ret;
  if ( m_fembSet.count(acd.fembID()) == 0 ) {
    if ( m_LogLevel >= 2 ) {
      cout << myname << "Skipping channel " << acd.channel() << " in FEMB " << acd.fembID() << endl;
    }
    return ret;
  }
  Index nsam = acd.samples.size();
  if ( nsam == 0 ) {
    if ( m_LogLevel >= 2 ) {
      cout << myname << "Skipping channel " << acd.channel() << " with no samples." << endl;
    }
    return ret;
  }
  if ( m_SamplingRatio < 0.01 ) return ret;
  AdcSignalVector newsams;
  newsams.reserve(nsam/m_SamplingRatio+1);
  // Loop over indices for new sample.
  if ( m_LogLevel >= 3 ) {
    cout << myname << "Updating samples for channel " << acd.channel() << "." << endl;
  }
  for ( Index j=0; ; ++j ) {
    double xj = j*m_SamplingRatio - m_SamplingOffset;
    double yj = acd.samples[0];
    if ( xj > 0.0 ) {
      Index i1 = xj;
      Index i2 = i1 + 1;
      if ( i2 >= acd.samples.size() ) break;
      double y1 = acd.samples[i1];
      double y2 = acd.samples[i2];
      yj = (y2 - y1)*(xj - i1) + y1;
      if ( m_LogLevel >= 3 ) {
        cout << myname << "  " << j << ": " << y1 << " , " << y2 << " ==> " << yj << endl;
      }
    } else {
      if ( m_LogLevel >= 3 ) {
        cout << myname << "  " << j << ": " << " Keeping first sample: " << yj << endl;
      }
    }
    newsams.push_back(yj);
  }
  if ( m_LogLevel >= 2 ) {
    cout << myname << "Updating samples for channel " << acd.channel() << ": size "
         << nsam << " --> " << newsams.size() << endl;
  }
  acd.samples = newsams;
  acd.signal.clear();
  acd.rois.clear();
  return ret;
}

//**********************************************************************

DEFINE_ART_CLASS_TOOL(AdcTimingMitigator)
