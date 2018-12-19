// GausRmsFitter.h

#ifndef GausRmsFitter_H
#define GausRmsFitter_H

// David Adams
// December 2018
//
// Utility to evaluate the mean and RMS for a peak in a histogram and
// return those and the selected range in an attached gaus function with
// fitted height.
//
// There are two configuration parameters: sigma0 and nsigma.
// Calls to fit a histogram add a starting mean, mean0.
//
// If both sigma0 and nsigma are valid (i.e. are > 0), then the RMS is
// evaluated iteratively over the range (mean-nsigma*sigma, mean+nsigma*sigma)
// until it coverges on a value for sigma. The range is restricted to that of
// the input, i.e. its default or the value set with axis SetRange or SetRangeUser.
// The starting values for the iteration are mean0 and sigma0 and the iteration
// stops when the evaluated value for sigma is less or close to that used to
// define the range.
//
// If sigma0 is valid and nsigma is not, then those are used for the 
// corresponding function values.
//
// If both are invalid, mean0 is ignored and man and RMS for the full range
// are used to construct the function.
//
// Any other configuration results in a fit error.

#include <string>
class TH1;

class GausRmsFitter {

public:

  using Name = std::string;
  using Index = unsigned int;

  // Ctor for fit of full histogram.
  GausRmsFitter(Name fnam);

  // Ctor to fit range of n-sigma starting with sigma0.
  GausRmsFitter(double sigma0, double nsigma, Name fnam);

  // Fit a histogram with the provided starting mean.
  // Returns 0 for success.
  // If successful, a gaus function is attached to the histogram.
  // This function has the assigned mean, sigma and fit range and the fitted height.
  int fit(TH1* ph, double mean0) const;

  // Set the log level.
  void setLogLevel(Index lev) { m_LogLevel = lev; }

private:

  double m_sigma0;
  double m_nsigma;
  Name m_fnam;
  Index m_LogLevel;

};

#endif
