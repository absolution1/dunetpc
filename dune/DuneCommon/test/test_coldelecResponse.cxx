// test_coldelecResponse.cxx

// David Adams
// May 2018
//
// This is a test and demonstration for coldelecResponse.

#undef NDEBUG

#include "dune/DuneCommon/coldelecResponse.h"
#include <string>
#include <iostream>
#include <iomanip>
#include <vector>
#include <cassert>

using std::string;
using std::cout;
using std::endl;
using std::setw;
using std::vector;

using Index = unsigned int;

//**********************************************************************

int test_coldelecResponse() {
  const string myname = "test_coldelecResponse: ";
  cout << myname << "Starting test" << endl;
#ifdef NDEBUG
  cout << myname << "NDEBUG must be off." << endl;
  abort();
#endif
  string line = "-----------------------------";
  string scfg;

  cout << myname << line << endl;
  cout << "Testing TF1." << endl;
  double gain = 100.0;
  double shap = 4.0;
  double t0 = 10.0;
  TF1* pfs[4];
  pfs[0] = coldelecResponseTF1(gain, shap, t0);
  pfs[1] = coldelecResponseTF1(gain, shap, t0+0.2);
  pfs[2] = coldelecResponseTF1(gain, shap, t0+0.4);
  pfs[3] = coldelecResponseTF1(gain, shap, t0+0.6);
  cout.precision(4);
  for ( double t=0.0; t<30.0; t+=1.0 ) {
    cout << setw(10) << t << ": ";
    for ( unsigned int ifun=0; ifun<4; ++ifun ) {
      double wf = pfs[ifun]->Eval(t);
      cout << setw(10) << wf;
    }
    cout << endl;
  }

  cout << myname << line << endl;
  cout << myname << "Done." << endl;
  return 0;
}

//**********************************************************************

int main() {
  return test_coldelecResponse();
}

//**********************************************************************
