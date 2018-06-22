// StickyCodeMetrics.cxx

#include "StickyCodeMetrics.h"
#include <map>
#include <iostream>
#include <sstream>

#include "TH1.h"

using std::map;
using std::string;
using std::cout;
using std::endl;
using std::ostringstream;

namespace {
using Index = unsigned int;
}

//**********************************************************************

StickyCodeMetrics::StickyCodeMetrics(const AdcCountVector& adcs) {
  for ( AdcCount adc : adcs ) {
    if ( m_counts.find(adc) == m_counts.end() ) m_counts[adc] = 0;
    ++m_counts[adc];
  }
  evaluateMetrics();
}

//**********************************************************************

StickyCodeMetrics::StickyCodeMetrics(const BinCounter& counts):
m_counts(counts) {
  evaluateMetrics();
}

//**********************************************************************

StickyCodeMetrics::StickyCodeMetrics(const TH1* pha) {
  const string myname = "StickyCodeMetrics::ctor: ";
  Index nbin = pha->GetNbinsX();
  Index iadc0 = pha->GetBinLowEdge(1);
  Index iadcLast = pha->GetBinLowEdge(nbin);
  Index nhadc = iadcLast - iadc0 + 1;
  if ( nhadc != nbin ) {
    cout << myname << "Histogram has inconssitent binning." << endl;
    return;
  }
  for ( Index ibin=0; ibin<nbin; ++ibin ) {
    m_counts[iadc0+ibin] = pha->GetBinContent(ibin+1);
  }
  evaluateMetrics();
}

//**********************************************************************

int StickyCodeMetrics::evaluateMetrics() {
  m_nsample = 0;
  Index badadc = 999999;
  Index maxadc = badadc;
  Index maxCount = 0;
  Index maxCount2 = 0;
  Index nmod0 = 0;
  Index nmod1 = 0;
  Index nmod63 = 0;
  double adcSum = 0.0;
  double countSum = 0.0;
  for ( BinCounter::value_type ibco : m_counts ) {
    Index iadc = ibco.first;
    double count = ibco.second;
    m_nsample += count;
    if ( count > maxCount ) {
      maxCount = count;
      maxadc = iadc;
    }
    AdcCount adcmod = iadc%64;
    if ( adcmod ==  0 ) nmod0 += count;
    if ( adcmod ==  1 ) nmod1 += count;
    if ( adcmod == 63 ) nmod63 += count;
    countSum += count;
    adcSum += count*iadc;
  }
  double countSum2 = 0;
  double adcSum2 = 0;
  Index maxadc2 = badadc;
  for ( BinCounter::value_type ibco : m_counts ) {
    Index iadc = ibco.first;
    double count = ibco.second;
    if ( iadc != maxadc ) {
      countSum2 += count;
      adcSum2 += count*iadc;
      if ( count > maxCount2 ) {
        maxCount2 = count;
        maxadc2 = iadc;
      }
    }
  }
  m_maxAdc = maxadc;
  m_maxAdc2 = maxadc2;
  m_meanAdc = countSum>0 ? adcSum/countSum : -1.0;
  m_meanAdc2 = countSum2>0 ? adcSum2/countSum2 : -1.0;
  m_maxFraction = maxCount/countSum;
  m_zeroFraction = nmod0/countSum;
  m_oneFraction = nmod1/countSum;
  m_highFraction = nmod63/countSum;
  return 0;
}

//**********************************************************************

DataMap StickyCodeMetrics::getMetrics(string prefix) const {
  DataMap res;
  res.setInt(  prefix + "MaxAdc", maxAdc());
  res.setInt(  prefix + "MaxAdc2", maxAdc2());
  res.setFloat(prefix + "MeanAdc", meanAdc());
  res.setFloat(prefix + "MeanAdc2", meanAdc2());
  res.setFloat(prefix + "MaxFraction", maxFraction());
  res.setFloat(prefix + "ZeroFraction", zeroFraction());
  res.setFloat(prefix + "OneFraction", oneFraction());
  res.setFloat(prefix + "HighFraction", highFraction());
  res.setFloat(prefix + "ClassicFraction", classicFraction());
  return res;
}

//**********************************************************************

void StickyCodeMetrics::print(string prefix) const {
  ostringstream sout;
  sout.precision(2);
  sout << prefix << "             # samples: " << nsample();
  sout << "\n";
  sout << prefix << "  Most common ADC code: " << maxAdc();
  sout << "\n";
  sout << prefix << " Next most common code: " << maxAdc2();
  sout << "\n";
  sout << prefix << "              Mean ADC: " << std::fixed << meanAdc();
  sout << "\n";
  sout << prefix << "      Mean ADC w/o max: " << std::fixed << meanAdc2();
  sout << "\n";
  sout.precision(3);
  sout << prefix << "       Frac in max bin: " << maxFraction();
  sout << "\n";
  sout << prefix << "            Frac LSB=0: " << zeroFraction();
  sout << "\n";
  sout << prefix << "            Frac LSB=1: " << oneFraction();
  sout << "\n";
  sout << prefix << "           Frac LSB=64: " << highFraction();
  sout << "\n";
  sout << prefix << "         Frac LSB=0,64: " << classicFraction();
  cout << sout.str() << endl;
}

//**********************************************************************
