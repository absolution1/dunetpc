// AdcRangeSampleScaler_tool.cc

#include "AdcRangeSampleScaler.h"
#include <iostream>

using std::string;
using std::cout;
using std::endl;

using Index = unsigned int;

//**********************************************************************
// Class methods.
//**********************************************************************

AdcRangeSampleScaler::AdcRangeSampleScaler(fhicl::ParameterSet const& ps)
: m_LogLevel(ps.get<int>("LogLevel")),
  m_RangeLimits(ps.get<IndexVector>("RangeLimits")),
  m_RangeModulus(ps.get<Index>("RangeModulus")),
  m_ScaleFactors(ps.get<FloatVector>("ScaleFactors")) {
  const string myname = "AdcRangeSampleScaler::ctor: ";
  if ( m_LogLevel >= 1 ) {
    cout << myname << "Parameters:" << endl;
    cout << myname << "         LogLevel: " << m_LogLevel << endl;
    Index imin = 0;
    cout << myname << "          Modulus: " << m_RangeModulus << endl;
    cout << myname << "     Range: Scale: " << endl;
    Index imod = m_RangeModulus;
    for ( Index iran=0; ; ++iran ) {
      // No more ranges.
      // Unless protected with a modulus, we need one more scale factor.
      Index imax = 0;
      bool haveMax = true;
      if ( iran >= m_RangeLimits.size() )  {
        if ( imod && imin >= imod ) {
          break;
        }
        haveMax = false;
      } else {
        imax = m_RangeLimits[iran];
        if ( imod > 0 && imax >= imod ) {
          cout << myname << "WARNING: Index exceeds modulus: " << imax << " >= " << imod << endl;
          break;
        }
      }
      if ( m_ScaleFactors.size() < iran + 1 ) {
        cout << myname << "ERROR: Too few scale factors." << endl;
        break;
      }
      float sfac = m_ScaleFactors[iran];
      cout << myname << "       [" << imin << ", ";
      if ( haveMax ) cout << imax;
      else cout << "...";
      cout << "): " << sfac << endl;
      if ( haveMax ) imin = imax;
      else break;
    }
      
  }
}

//**********************************************************************

DataMap AdcRangeSampleScaler::update(AdcChannelData& acd) const {
  const string myname = "AdcRangeSampleScaler::update: ";
  if ( m_LogLevel >= 2 ) cout << "Processing run " << acd.run() << " event " << acd.event()
                              << " channel " << acd.channel() << endl;
  DataMap ret;

  // Find the index for the scale factor.
  if ( acd.channel() == AdcChannelData::badIndex() ) {
     cout << myname << "ERROR: Invalid channel." << endl;
    return ret.setStatus(1);
  }
  Index icha = acd.channel();
  Index ichaMod = m_RangeModulus ? icha%m_RangeModulus : icha;
  Index iran = 0;
  for ( iran=0; iran<m_RangeLimits.size(); ++iran ) {
    Index imax = m_RangeLimits[iran];
    if ( ichaMod < imax ) {
      break;
    }
  }

  if ( iran >= m_ScaleFactors.size() ) {
    cout << myname << "ERROR: No scale factor for channel " << icha;
    if ( ichaMod != icha ) cout << " (" << ichaMod << ")";
    cout << endl;
    return ret.setStatus(2);
  }
  float sfac = m_ScaleFactors[iran];

  // Scale samples.
  for ( float& sam : acd.samples ) sam *= sfac;

  ret.setFloat("arssScaleFactor", sfac);

  return ret;
}

//**********************************************************************

DEFINE_ART_CLASS_TOOL(AdcRangeSampleScaler)
