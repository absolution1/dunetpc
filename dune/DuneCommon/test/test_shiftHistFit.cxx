// test_shiftHistFit.cxx

// David Adams
// May 2021
//
// This is a test and demonstration for shiftHistFit.

#undef NDEBUG

#include "dune/DuneCommon/shiftHistFit.h"
#include "dune/DuneCommon/gausTF1.h"
#include "dune/DuneCommon/TPadManipulator.h"
#include <string>
#include <iostream>
#include <sstream>
#include <iomanip>
#include <vector>
#include <cassert>
#include "TH1F.h"
#include "TF1.h"
#include "TRandom.h"

using std::string;
using std::cout;
using std::endl;
using std::setw;
using std::vector;
using std::istringstream;

using Index = unsigned int;

//**********************************************************************

int test_shiftHistFit(double xshift) {
  const string myname = "test_shiftHistFit: ";
  cout << myname << "Starting test" << endl;
  cout << myname << "Shift is " << xshift << endl;
#ifdef NDEBUG
  cout << myname << "NDEBUG must be off." << endl;
  abort();
#endif
  string line = "-----------------------------";
  string scfg;

  cout << myname << line << endl;
  cout << "Testing shifHistFit." << endl;

  cout << myname << line << endl;
  cout << "Build histogram" << endl;
  // Truth parameters before shift.
  double x0 = 50.0;
  double sig0 = 10.0;
  // Histogram range.
  double hx1 = x0 + xshift - 50.0;
  double hx2 = x0 + xshift + 50.0;
  // Fit range.
  double fx1 = x0 + xshift - 25.0;
  double fx2 = x0 + xshift + 25.0;
  // Fill histgram with data distributed around shifted truth.
  TH1* ph = new TH1F("hgau", "My gaus", 100, hx1, hx2);
  for ( int ient=0; ient<1000; ++ient ) {
    ph->Fill(gRandom->Gaus(x0 + xshift, sig0));
  }

  cout << myname << line << endl;
  cout << "Build function." << endl;
  TF1* pf = gausTF1(10.0, x0 + xshift + 0.3*sig0, 0.8*sig0, "mygaus");
  pf->SetRange(fx1, fx2);
  pf->Print();

  //ph->Fit(pf);
  shiftHistFit(ph, pf, "R", "Mean", xshift);
  pf->Print();

  TPadManipulator man;
  man.add(ph, "hist");
  man.add(pf, "same");
  man.addAxis();
  string fnout = "test_shiftHistFit.png";
  cout << myname << "Printing " << fnout << endl;
  man.print(fnout);

  return 0;
}

//**********************************************************************

int main(int narg, const char* argc[]) {
  double xshift = 1000.0;
  if ( narg > 1 ) {
    istringstream ssin(argc[1]);
    ssin >> xshift;
  }
  return test_shiftHistFit(xshift);
}

//**********************************************************************
