// quietHistFit.h
//
// David Adams
// May 2018
//
// Call TF1::Fit after suppressing Root logging.

#ifndef quietHistFit_H
#define quietHistFit_H

#include <string>

class TH1;
class TF1;

// Return ph->Fit(fname, fopt)
int quietHistFit(TH1* ph, std::string fname, std::string fopt);

// Return ph->Fit(pf, fopt)
int quietHistFit(TH1* ph, TF1* pf, std::string fopt);

#endif
