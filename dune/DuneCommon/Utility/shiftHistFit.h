// shiftHistFit.h
//
// David Adams
// May 2021
//
// Fit a histogram with a TF1. The fitting is done after shifting the histogram
// and function by a specified distance with results presented at the original
// position. The goal is to avoid numerical problame that can arise when the
// position coordinate, e.g. TPC ticks, reaches large values.
//
// Call TF1::Fit after suppressing Root logging.
//   ph - histogram to be fit
//   pf - fit function
//   fopt - Root fitting options
//   spar or ipar - name or index of the position parameter in the fit function
//   xshift - fit h'(x') = h(x) where x = x' + xshift
//   Return value is fit status. Include "S" in fopt to get IsValid
//   from the fit result pointer.

#ifndef shiftHistFit_H
#define shiftHistFit_H

#include <string>

class TH1;
class TF1;

int shiftHistFit(TH1* ph, TF1* pf, std::string fopt,         int ipar, double xshift);
int shiftHistFit(TH1* ph, TF1* pf, std::string fopt, std::string spar, double xshift);

#endif
