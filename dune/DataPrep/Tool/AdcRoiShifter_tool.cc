// AdcRoiShifter_tool.cc

#include "AdcRoiShifter.h"
#include <iostream>

using std::string;
using std::cout;
using std::endl;

//**********************************************************************
// Class methods.
//**********************************************************************

AdcRoiShifter::AdcRoiShifter(fhicl::ParameterSet const& ps)
: m_LogLevel(ps.get<int>("LogLevel")),
  m_BinOffset(ps.get<int>("BinOffset")) {
  const string myname = "AdcRoiShifter::ctor: ";
  if ( m_LogLevel >= 1 ) {
    cout << myname << "Configuration parameters:" << endl;
    cout << myname << "      LogLevel: " << m_LogLevel << endl;
    cout << myname << "     BinOffset: " << m_BinOffset<< endl;
  }
}

//**********************************************************************

DataMap AdcRoiShifter::update(AdcChannelData& acd) const {
  const string myname = "AdcRoiShifter::update: ";
  DataMap ret;
  AdcIndex nsam = acd.signal.size();
  if ( nsam ==0 || m_BinOffset == 0 ) {
    return ret;
  } else if ( m_BinOffset > 0 ) {
    AdcIndex shift = m_BinOffset;
    for ( AdcIndex isam=nsam-1; isam>=shift; --isam ) {
      acd.signal[isam] = acd.signal[isam-shift];
    }
    for ( AdcIndex isam=0; isam<shift; ++isam ) {
      acd.signal[isam] = false;
    }
  } else {
    AdcIndex shift = -m_BinOffset;
    for ( AdcIndex isam=0; isam<nsam-shift; ++isam ) {
      acd.signal[isam] = acd.signal[isam+shift];
    }
    for ( AdcIndex isam=nsam-shift; isam<nsam; ++isam ) {
      acd.signal[isam] = false;
    }
  }
  acd.roisFromSignal();
  if ( m_LogLevel >= 3 ) {
    cout << myname << "  # ROI: " << acd.rois.size() << endl;
  }
  ret.setInt("nroi", acd.rois.size());
  return ret;
}

//**********************************************************************

DataMap AdcRoiShifter::view(const AdcChannelData& acd) const {
  AdcChannelData acdtmp;
  acdtmp.signal = acd.signal;
  return update(acdtmp);
}

//**********************************************************************

DEFINE_ART_CLASS_TOOL(AdcRoiShifter)
