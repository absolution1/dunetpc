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

double gausFunction(double* x, double* pars) {
  double pos = x[0];
  double height = pars[0];
  double mean = pars[1];
  double sigma = pars[2];
  if ( sigma <= 0.0 ) return 0;
  double rat = (pos-mean)/sigma;
  return height*exp(-0.5*rat*rat);
}

TF1* gausTF1(double heightIn =0.0, double meanIn =0.0, double sigmaIn =0.0, std::string fname ="gaus") {
  double height = heightIn != 0.0 ? heightIn : 1.0;
  double sigma = sigmaIn > 0.0 ? sigmaIn : 1.0;
  double mean = meanIn;
  TF1* pf = new TF1(fname.c_str(), gausFunction, -10, 10, 3);
  pf->SetParName(0, "Height");
  pf->SetParName(1, "Mean");
  pf->SetParName(2, "Sigma");
  pf->SetParameter(0, height);
  pf->SetParameter(1, mean);
  pf->SetParameter(2, sigma);
  double sigmaMin = 0.0;
  double sigmaMax = 1.e10;
  if ( sigmaIn > 0.0 ) {
    sigmaMin =  0.1*sigmaIn;
    sigmaMax = 10.0*sigmaIn;
  }
  pf->SetParLimits(2, sigmaMin, sigmaMax);
  if ( sigmaIn > 0.0 ) pf->SetParLimits(2, 0.1*sigma, 10*sigma);
  return pf;
}


#endif
