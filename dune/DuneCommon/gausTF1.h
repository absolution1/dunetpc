// gausTF1.h
//
// David Adams
// May 2018
//
// Create TF1 for Gaussian distribution.
//
// Normalized so that height = pf->Eval(0).

#ifndef gausTF1_H
#define gausTF1_H

#include "TF1.h"
#include <string>

using std::string;

double gausFunction(double* x, double* pars);

TF1* gausTF1(double heightIn =0.0, double meanIn =0.0, double sigmaIn =0.0, std::string fname ="dunegaus");

#endif
