// test_offsetLine.cxx

// David Adams
// May 2018
//
// This is a test and demonstration for offsetLine.

#undef NDEBUG

#include "dune/DuneCommon/offsetLine.h"
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

int test_offsetLine() {
  const string myname = "test_offsetLine: ";
  cout << myname << "Starting test" << endl;
#ifdef NDEBUG
  cout << myname << "NDEBUG must be off." << endl;
  abort();
#endif
  string line = "-----------------------------";
  string scfg;

  cout << myname << line << endl;
  cout << "Testing TF1." << endl;

  cout << myname << line << endl;
  cout << "Build waveforms." << endl;
  vector<double> offs = {-3.5, -1.2, 4.4};
  vector<double> slps = {-1.0, 1.0, 2.0};
  vector<TF1*> pfs;
  for ( double slp : slps ) {
    for ( double off : offs ) {
      pfs.push_back(offsetLineTF1(off, slp));
    }
  }

  cout << myname << line << endl;
  cout << "Display waveforms." << endl;
  cout.precision(2);
  int w = 7;
  TF1* pf0 = pfs[0];
  for ( int ipar=0; ipar<pf0->GetNpar(); ++ipar ) {
    string parName = pf0->GetParName(ipar);
    cout << setw(10) << parName << ":";
    for ( TF1* pf : pfs ) cout << setw(w) << std::fixed << pf->GetParameter(ipar);
    cout << endl;
  }
  for ( double t=-8.0; t<8.0; t+=1.0 ) {
    cout << setw(10) << t << ":";
    for ( TF1* pf : pfs ) {
      double wf = pf->Eval(t);
      cout << setw(w) << std::fixed << wf;
    }
    cout << endl;
  }

  cout << myname << line << endl;
  cout << myname << "Done." << endl;
  return 0;
}

//**********************************************************************

int main() {
  return test_offsetLine();
}

//**********************************************************************
