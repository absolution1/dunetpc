// AdcRegularSignalFinder_tool.cc

#include "AdcRegularSignalFinder.h"
#include <iostream>

using std::string;
using std::cout;
using std::endl;

using Index = AdcIndex;

//**********************************************************************
// Class methods.
//**********************************************************************

AdcRegularSignalFinder::
AdcRegularSignalFinder(AdcIndex per, AdcIndex len, int lev)
: m_LogLevel(lev),
  m_Period(per),
  m_Length(len) {
  const string myname = "AdcRegularSignalFinder::ctor: ";
  if ( m_LogLevel >= 1 ) {
    cout << myname << "Configuration parameters:" << endl;
    cout << myname << "  LogLevel: " << m_LogLevel << endl;
    cout << myname << "    Period: " << m_Period <<  endl;
    cout << myname << "    Length: " << m_Length << endl;
  }
}

//**********************************************************************

AdcRegularSignalFinder::
AdcRegularSignalFinder(fhicl::ParameterSet const& ps)
: AdcRegularSignalFinder(
   ps.get<unsigned int>("Period"),
   ps.get<unsigned int>("Length"),
   ps.get<int>("LogLevel")
  ) { }

//**********************************************************************

DataMap AdcRegularSignalFinder::update(AdcChannelData& acd) const {
  const string myname = "AdcRegularSignalFinder::update: ";
  Index nsam = acd.samples.size();
  Index nper = m_Period;
  Index nlen = m_Length == 0 ? m_Period : m_Length;
  Index nroi = 0;
  if ( nsam > 0 && nper > 0 ) {
    acd.signal.resize(nsam, true);
    Index nrem = nsam % nper;
    nroi = nsam/nper + (nrem>0);
cout << myname << "nsam=" << nsam << ", nper=" << nper << ", nroi=" << nroi << endl;
    for ( Index iroi=0; iroi<nroi; ++iroi ) {
      Index isam1 = nper*iroi;
      Index isam2 = isam1 + nlen;
      Index isam3 = isam1 + nper;
      if ( isam2 > nsam ) isam2 = nsam;
      if ( isam3 > nsam ) isam2 = nsam;
      for ( Index isam=isam2; isam<isam3; ++isam ) {
        acd.signal[isam] = false;
      }
      acd.rois.emplace_back(isam1, isam2-1);
    }
  } else {
    acd.signal.resize(nsam, false);
    acd.rois.clear();
  }
  DataMap res(0);
  res.setInt("roiPeriod", nper);
  res.setInt("roiLength", nlen);
  res.setInt("roiCount", nroi);
  return res;
}

//**********************************************************************

DataMap AdcRegularSignalFinder::view(const AdcChannelData& acd) const {
  AdcChannelData acdtmp;
  acdtmp.samples = acd.samples;
  return update(acdtmp);
}

//**********************************************************************
