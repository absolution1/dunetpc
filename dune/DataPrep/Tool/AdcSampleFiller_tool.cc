// AdcSampleFiller_tool.cc

#include "AdcSampleFiller.h"
#include <iostream>

using std::string;
using std::cout;
using std::endl;

//**********************************************************************
// Class methods.
//**********************************************************************

AdcSampleFiller::AdcSampleFiller(fhicl::ParameterSet const& ps)
: m_LogLevel(ps.get<int>("LogLevel")),
  m_AdcUnderflow(ps.get<unsigned int>("AdcUnderflow")),
  m_AdcOverflow(ps.get<unsigned int>("AdcOverflow")) {
  const string myname = "AdcSampleFiller::ctor: ";
  if ( m_LogLevel >= 1 ) {
    cout << myname << "Configuration parameters:" << endl;
    cout << myname << "      LogLevel: " << m_LogLevel << endl;
    cout << myname << "  AdcUnderflow: " << m_AdcUnderflow << endl;
    cout << myname << "   AdcOverflow: " << m_AdcOverflow << endl;
  }
}

//**********************************************************************

DataMap AdcSampleFiller::update(AdcChannelData& acd) const {
  const string myname = "AdcSampleFiller::update: ";
  AdcIndex nsam = acd.raw.size();
  acd.samples.resize(nsam, AdcChannelData::badSignal);
  acd.flags.resize(nsam, AdcGood);
  AdcPedestal ped = acd.pedestal;
  if ( m_LogLevel >= 2 ) {
    cout << myname << "# samples: " << nsam << endl;
    cout << myname << " Pedestal: " << ped << endl;
  }
  if ( ped == AdcChannelData::badSignal ) {
    cout << myname << "Pedestal is not set." << endl;
    return DataMap(1);
  }
  AdcIndex nlo = 0;
  AdcIndex nhi = 0;
  AdcIndex nbad = 0;
  for ( AdcIndex isam=0; isam<nsam; ++isam ) {
    AdcIndex adc = acd.raw[isam];
    acd.samples[isam] = adc - ped;
    if ( adc <= m_AdcUnderflow ) {
      acd.flags[isam] = AdcUnderflow;
      ++nlo;
    } else if ( adc >= m_AdcOverflow ) {
      acd.flags[isam] = AdcOverflow;
      ++nhi;
    } else {
      acd.flags[isam] = AdcGood;
    }
    if ( adc < m_AdcUnderflow || adc > m_AdcOverflow ) ++nbad;
  }
  acd.sampleUnit = "ADC counts";
  DataMap res(0);
  res.setInt("nSample",  nsam);
  res.setInt("nUnderflow",  nlo);
  res.setInt("nOverflow",   nhi);
  res.setInt("nOutOfRange", nbad);
  return res;
}

//**********************************************************************

DataMap AdcSampleFiller::view(const AdcChannelData& acd) const {
  AdcChannelData acdtmp;
  acdtmp.samples = acd.samples;
  return update(acdtmp);
}

//**********************************************************************
