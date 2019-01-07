// offsetLine.h
//
// David Adams
// June 2018
//
// The offsetLine function is two line segments with equal slopes
// and equal-magnitude, opposite-sign x-offsets on either side of
// the origin. The function is continuous for offset >= 0 and has
// a discontinuity at x=0 for offset < 0.
//
// The offsetLinePed function adds a pedestal to offsetLIne.

#ifndef offsetLine_H
#define offsetLine_H

#include <string>
#include "TF1.h"

inline
double offsetLineFull(double* px, double* ppar) {
  double xoff = ppar[0];
  double scal = ppar[1];
  double yped = ppar[2];
  double sneg = ppar[3];
  double x = px[0];
  double gain = scal;
  if ( x < 0.0 ) gain *= sneg;
  if ( xoff > 0.0 ) {
    if ( x > xoff ) return yped + gain*(x-xoff);
    if ( x < -xoff ) return yped + gain*(x+xoff);
  } else {
    double y0 = -gain*xoff;
    if ( x > 0.0 ) return yped + y0 + gain*x;
    else if ( x < 0.0 ) return yped - y0 + gain*x;
  }
  return 0.0;
}

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
double offsetLinePed(double* px, double* ppar) {
  double yped = ppar[2];
  return yped + offsetLine(px, ppar);
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

inline
TF1* offsetLinePedTF1(double off =0.0, double slope =1.0, double ped =0.0,
                      double xmin =-10.0, double xmax =10.0,
                      std::string fname ="offsetLinePed") {
  TF1* pf = new TF1(fname.c_str(), offsetLinePed, xmin, xmax, 3);
  pf->SetParName(0, "Offset");
  pf->SetParName(1, "Slope");
  pf->SetParName(2, "Pedestal");
  pf->SetParameter(0, off);
  pf->SetParameter(1, slope);
  pf->SetParameter(2, ped);
  return pf;
}
  
inline
TF1* offsetLineFullTF1(double off =0.0, double slope =1.0,
                       double ped =0.0, double sneg =1.0,
                       double xmin =-10.0, double xmax =10.0,
                       std::string fname ="offsetLinePed") {
  TF1* pf = new TF1(fname.c_str(), offsetLineFull, xmin, xmax, 4);
  pf->SetParName(0, "Offset");
  pf->SetParName(1, "Slope");
  pf->SetParName(2, "Pedestal");
  pf->SetParName(3, "NegScale");
  pf->SetParameter(0, off);
  pf->SetParameter(1, slope);
  pf->SetParameter(2, ped);
  pf->SetParameter(3, sneg);
  return pf;
}
  
#endif
