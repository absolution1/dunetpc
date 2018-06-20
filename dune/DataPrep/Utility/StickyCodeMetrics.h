// StickyCodeMetrics.h
//
// David Adams
// June 2018
//
// Class to evaluate sticky code metrics.
//
//  maxAdc = most common ADC code 
//  maxAdc2 = 2nd most common ADC code 
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
#include "dune/DuneInterface/Data/DataMap.h"
#include <map>

class TH1;

class StickyCodeMetrics {

public:

  using Index = unsigned int;
  using BinCounter = std::map<Index, double>;

  // Ctor from a count map.
  StickyCodeMetrics(const BinCounter& counts);

  // Ctor from a vector of ADC codes.
  StickyCodeMetrics(const AdcCountVector& adcs);

  // Ctor from a histogram.
  // Each histogram bin must correspond to one ADC bin.
  StickyCodeMetrics(const TH1* pha);

  // Metrics.
  Index nsample() const { return m_nsample; } 
  AdcIndex maxAdc() const { return m_maxAdc; }
  AdcIndex maxAdc2() const { return m_maxAdc2; }
  double meanAdc() const { return m_meanAdc; }
  double meanAdc2() const { return m_meanAdc2; }
  double maxFraction() const { return m_maxFraction; }
  double zeroFraction() const { return m_zeroFraction; }
  double oneFraction() const { return m_oneFraction; }
  double highFraction() const { return m_highFraction; }
  double classicFraction() const { return m_zeroFraction + m_highFraction; }

  // Return the ADC bin count map.
  const BinCounter& getBinCounts() { return m_counts; }

  // Return the metrics in a data map (typed name-value pairs).
  // Sets int or float values.
  // Name is prefix + capitalized metric name, e.g. scmMaxAdc
  DataMap getMetrics(std::string prefix ="scm") const;

  // Display results.
  void print(std::string prefix ="") const;

private:

  Index m_nsample;
  AdcIndex m_maxAdc;
  AdcIndex m_maxAdc2;
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
