// StickyCodeMetrics.h
//
// David Adams
// June 2018
//
// Class to evaluate sticky code metrics from a distribution
// of ADC codes presumed to correpond to the same input signal,
// e.g. a pedestal or tickmod distrobution.
//
// The metrics are:
//
//   maxAdc = most common ADC code 
//   maxAdc2 = 2nd most common ADC code 
//   meanAdc = mean ADC
//   meanAdc2 = mean ADC without the most common code
//   maxFraction (s1) = fraction of ticks with the most common code
//   zeroFraction = fraction of ticks with ADC%64=0
//   oneFraction = fraction of ticks with ADC%64=1
//   highFraction = fraction of ticks with ADC%64=63
//   classicFraction (s2) = fraction of ticks with ADC%64 = 0 or 63
//   fitStatus = status code from fit (0 for OK)
//   fitMean = mean ADC from fit
//   fitSigma = ADC sigma from fit
//   fitExcess = fraction with the most common code above fit
//
// If a histogram is created, the following are also defined:
//   hist - pointer to histogram accessed through getHist() or getSharedHist()
//
// Classic sticky codes are code%64 = 0 or 63
//
// Caller passes a distribution of ADC codes in one of three forms:
//   Bin counter - Map with count for each code
//   Code vector - Vector of codes (order is not used)
//   Histogram - TH1 with count or likelihood for each code
//
// The histogram must have one bin per ADC code.

#ifndef StickyCodeMetrics_H
#define StickyCodeMetrics_H

#include "dune/DuneInterface/AdcTypes.h"
#include "dune/DuneInterface/Data/DataMap.h"
#include <map>
#include <memory>

class TH1;

class StickyCodeMetrics {

public:

  using Index = unsigned int;
  using BinCounter = std::map<Index, double>;
  using HistPtr = std::shared_ptr<TH1>;
  using Name = std::string;

  // Ctor. No histogram is connfigured.
  StickyCodeMetrics() = default;

  // Ctor to configure for histogram with nbin bins with low edge
  // a multiple of lowbin;
  StickyCodeMetrics(Name hnam, Name httl, Index nbin, Index lowbin, float sigmaMin, float sigmaMax);

  // Evaluate a count map.
  int evaluate(const BinCounter& counts);

  // Evaluate a vector of ADC codes.
  int evaluate(const AdcCountVector& adcs);

  // Evaluate a histogram.
  // Each histogram bin must correspond to one ADC bin.
  int evaluate(const TH1* pha);

  // Clear the data and metrics.
  void clear();

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
  int fitStatus() const { return m_fitStatus; }
  double fitMean() const { return m_fitMean; }
  double fitSigma() const { return m_fitSigma; }
  double fitExcess() const { return m_fitExcess; }

  // Return the ADC bin count map.
  const BinCounter& getBinCounts() { return m_counts; }

  // Return the metrics in a data map (typed name-value pairs).
  // Sets int or float values.
  // Name is prefix + capitalized metric name, e.g. scmMaxAdc
  DataMap getMetrics(std::string prefix ="scm") const;

  // Return the histogram.
  TH1* getHist() { return getSharedHist().get(); };
  HistPtr getSharedHist() { return m_ph; }

  // Display results.
  void print(std::string prefix ="") const;

private:

  // Configuration.
  Name m_hnam;
  Name m_httl;
  int m_nbin =0;
  Index m_lowbin = 1;
  float m_sigmaMin = 0.1;
  float m_sigmaMax = 100.0;

  // Metrics.
  Index m_nsample;
  AdcIndex m_maxAdc;
  AdcIndex m_maxAdc2;
  double m_meanAdc;
  double m_meanAdc2;
  double m_maxFraction;
  double m_zeroFraction;
  double m_oneFraction;
  double m_highFraction;
  int m_fitStatus;
  double m_fitMean;
  double m_fitSigma;
  double m_fitExcess;

  BinCounter m_counts;
  HistPtr m_ph;

  int evaluateMetrics();

};

#endif
