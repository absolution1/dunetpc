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

StickyCodeMetrics::StickyCodeMetrics(const AdcCountVector& adcs, int) {
  m_nsample = adcs.size();
  map<AdcCount, Index> counts;
  int maxadc = -1;
  Index maxCount = 0;
  Index nmod0 = 0;
  Index nmod1 = 0;
  Index nmod63 = 0;
  Index adcsum = 0;
  for ( AdcCount adc : adcs ) {
    if ( counts.find(adc) == counts.end() ) counts[adc] = 0;
    ++counts[adc];
    if ( counts[adc] > maxCount ) {
      maxCount = counts[adc];
      maxadc = adc;
    }
    AdcCount adcmod = adc%64;
    if ( adcmod ==  0 ) ++nmod0;
    if ( adcmod ==  1 ) ++nmod1;
    if ( adcmod == 63 ) ++nmod63;
    adcsum += adc;
  }
  Index adcCount2 = 0;
  Index adcsum2 = 0;
  for ( auto ent : counts ) {
    AdcCount adc = ent.first;
    Index count = ent.second;
    if ( adc != maxadc ) {
      adcCount2 += count;
      adcsum2 += count*adc;
    }
  }
  double count = adcs.size();
  double count2 = adcCount2;
  m_maxAdc = maxadc;
  m_meanAdc = count>0 ? adcsum/count : -1.0;
  m_meanAdc2 = count2>0 ? adcsum2/count2 : -1.0;
  m_maxFraction = maxCount/count;
  m_zeroFraction = nmod0/count;
  m_oneFraction = nmod1/count;
  m_highFraction = nmod63/count;
}

//**********************************************************************

int StickyCodeMetrics::evaluateMetrics() {
  m_nsample = 0;
  Index badadc = 999999;
  Index maxadc = badadc;
  Index maxCount = 0;
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
  for ( BinCounter::value_type ibco : m_counts ) {
    Index iadc = ibco.first;
    double count = ibco.second;
    if ( iadc != maxadc ) {
      countSum2 += count;
      adcSum2 += count*iadc;
    }
  }
  m_maxAdc = maxadc;
  m_meanAdc = countSum>0 ? adcSum/countSum : -1.0;
  m_meanAdc2 = countSum2>0 ? adcSum2/countSum2 : -1.0;
  m_maxFraction = maxCount/countSum;
  m_zeroFraction = nmod0/countSum;
  m_oneFraction = nmod1/countSum;
  m_highFraction = nmod63/countSum;
  return 0;
}

//**********************************************************************

void StickyCodeMetrics::print(string prefix) const {
  ostringstream sout;
  sout.precision(2);
  sout << prefix << "        # samples: " << nsample();
  sout << "\n";
  sout << prefix << "          Max ADC: " << maxAdc();
  sout << "\n";
  sout << prefix << "         Mean ADC: " << std::fixed << meanAdc();
  sout << "\n";
  sout << prefix << " Mean ADC w/o max: " << std::fixed << meanAdc2();
  sout << "\n";
  sout.precision(3);
  sout << prefix << "  Frac in max bin: " << maxFraction();
  sout << "\n";
  sout << prefix << "       Frac LSB=0: " << zeroFraction();
  sout << "\n";
  sout << prefix << "       Frac LSB=1: " << oneFraction();
  sout << "\n";
  sout << prefix << "      Frac LSB=64: " << highFraction();
  sout << "\n";
  sout << prefix << "    Frac LSB=0,64: " << classicFraction();
  cout << sout.str() << endl;
}

//**********************************************************************
