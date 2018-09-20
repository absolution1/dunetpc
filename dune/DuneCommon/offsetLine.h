// offsetLine.h
//
// David Adams
// June 2018
//
// The offsetLine function is two line segments with equal slopes
// and equal-magnitude, opposite-sign x-offsets on either side of
// the origin. The function is continuous for offset >= 0 and has
// a discontinuity at x=0 for offset < 0.

#ifndef offsetLine_H
#define offsetLine_H

#include <string>
#include "TF1.h"

inline
double offsetLine(double* px, double* ppar) {
  double xoff = ppar[0];
  double gain = ppar[1];
  double x = px[0];
  if ( xoff > 0.0 ) {
    if ( x > xoff ) return gain*(x-xoff);
    if ( x < -xoff ) return gain*(x+xoff);
  } else {
    double y0 = -gain*xoff;
    if ( x > 0.0 ) return y0 + gain*x;
    else if ( x < 0.0 ) return -y0 + gain*x;
  }
  return 0.0;
}

inline
TF1* offsetLineTF1(double off =0.0, double slope =1.0,
                   double xmin =-10.0, double xmax =10.0,
                   std::string fname ="offsetLine") {
  TF1* pf = new TF1(fname.c_str(), offsetLine, xmin, xmax, 2);
  pf->SetParName(0, "Offset");
  pf->SetParName(1, "Slope");
  pf->SetParameter(0, off);
  pf->SetParameter(1, slope);
  return pf;
}

#endif
