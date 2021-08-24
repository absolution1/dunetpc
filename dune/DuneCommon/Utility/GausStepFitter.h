// GausStepFitter.h

#ifndef GausStepFitter_H
#define GausStepFitter_H

// David Adams
// December 2018
//
// Utility that fits a histogram with a Gaussian stating from specified
// position, width and fit window.

#include <string>
class TH1;

class GausStepFitter {

public:

  using Name = std::string;
  using Index = unsigned int;

  // Ctor.
  //   pos - starting position
  //   width - starting sigma
  //   xwin - us d to construct the fit range (pos-xwin, pos+xwin)
  //   fopt - Root fitting option.
  GausStepFitter(double pos, double sigma, double height, Name fnam, Name fopt ="WWS");

  // Fit a histogram.
  // Returns 0 for success.
  // If successful, a gaus function is attahed to the histogram.
  int fit(TH1* ph) const;

private:

  double m_pos;
  double m_sigma;
  double m_height;
  Name m_fnam;
  Name m_fopt;
  Index m_LogLevel;

};

#endif
