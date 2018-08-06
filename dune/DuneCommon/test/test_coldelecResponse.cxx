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

  cout << myname << line << endl;
  cout << "Build waveforms." << endl;
  double gain0 = 100.0;
  double shap0 = 5.0;
  double toff0 = 10.0;
  vector<double> dtoffs = {0.0, 0.2, 0.4, 0.6, 0.8};
  Index noff = dtoffs.size();
  vector<TF1*> pfs;
  for ( double fac : {1.0, 0.8, 0.5} ) {
    double gain = gain0/fac;
    double shap = fac*shap0;
    for ( Index ioff=0; ioff<noff; ++ioff ) {
      double toff = toff0 + dtoffs[ioff];
      pfs.push_back(coldelecResponseTF1(gain, shap, toff));
    }
  }
  vector<double> areas;
  vector<double> arats;
  for ( TF1* pf : pfs ) {
    double a = pf->Integral(0,30);
    areas.push_back(a);
    double arat = a/pf->GetParameter(0)/pf->GetParameter(1);
    arats.push_back(arat);
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
  for ( double t=0.0; t<30.0; t+=1.0 ) {
    cout << setw(10) << t << ":";
    for ( TF1* pf : pfs ) {
      double wf = pf->Eval(t);
      cout << setw(w) << std::fixed << wf;
    }
    cout << endl;
  }
  cout << setw(10) << "area" << ":";
  for ( double area : areas ) cout << setw(w) << std::fixed << area;
  cout << endl;
  cout.precision(4);
  cout << setw(10) << "A/h/s" << ":";
  for ( double arat : arats ) cout << setw(w) << std::fixed << arat;
  cout << endl;

  cout << myname << line << endl;
  cout << myname << "Done." << endl;
  return 0;
}

//**********************************************************************

int main() {
  return test_coldelecResponse();
}

//**********************************************************************
