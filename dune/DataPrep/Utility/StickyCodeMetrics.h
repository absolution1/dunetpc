// StickyCodeMetrics.h
//
// David Adams
// December 2017
//
// Class to evaluate sticky code metrics.
//
//  maxAdc = most common ADC code 
//  meanAdc = mean ADC
//  meanAdc2 = mean ADC without the most common code
//  maxFraction (s1) = fraction of ticks with the most common code
//  zeroFraction = fraction of ticks with ADC%64=0
//  oneFraction = fraction of ticks with ADC%64=1
//  highFraction = fraction of ticks with ADC%64=63
//  classicFraction (s2) = fraction of ticks with ADC%64 = 0 or 63
//
// Classic sticky codes are code%64 = 0 or 63

#ifndef StickyCodeMetrics_H
#define StickyCodeMetrics_H

#include "dune/DuneInterface/AdcTypes.h"
#include <map>

class TH1;

class StickyCodeMetrics {

public:

  using Index = unsigned int;
  using BinCounter = std::map<Index, double>;

  // Ctor from a count map.
  StickyCodeMetrics(const BinCounter& counts);

  // Ctor from a vector of ADC codes.
  // These typically are for a narrow range of input signals.
  StickyCodeMetrics(const AdcCountVector& adcs);

  // Ctor from a histogram.
  // Each histogram bin must correspond to one ADC bin.
  StickyCodeMetrics(const TH1* pha);

  // Old ctor for testing.
  StickyCodeMetrics(const AdcCountVector& adcs, int);

  // Metrics.
  Index nsample() const { return m_nsample; } 
  AdcIndex maxAdc() const { return m_maxAdc; }
  double meanAdc() const { return m_meanAdc; }
  double meanAdc2() const { return m_meanAdc2; }
  double maxFraction() const { return m_maxFraction; }
  double zeroFraction() const { return m_zeroFraction; }
  double oneFraction() const { return m_oneFraction; }
  double highFraction() const { return m_highFraction; }
  double classicFraction() const { return m_zeroFraction + m_highFraction; }

  // Display results.
  void print(std::string prefix ="") const;

private:

  Index m_nsample;
  AdcIndex m_maxAdc;
  double m_meanAdc;
  double m_meanAdc2;
  double m_maxFraction;
  double m_zeroFraction;
  double m_oneFraction;
  double m_highFraction;

  BinCounter m_counts;

  int evaluateMetrics();

};

#endif
