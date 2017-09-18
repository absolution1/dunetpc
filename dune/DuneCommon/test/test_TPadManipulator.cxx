// test_TPadManipulator.cxx

// David Adams
// September 2017
//
// This is a test and demonstration for TPadManipulator.

#undef NDEBUG

#include "../TPadManipulator.h"
#include <string>
#include <iostream>
#include <iomanip>
#include <cassert>
#include "TH1F.h"
#include "TCanvas.h"

using std::string;
using std::cout;
using std::endl;
using std::setw;

//**********************************************************************

int test_TPadManipulator() {
  const string myname = "test_TPadManipulator: ";
  cout << myname << "Starting test" << endl;
#ifdef NDEBUG
  cout << myname << "NDEBUG must be off." << endl;
  abort();
#endif
  string line = "-----------------------------";
  string scfg;

  cout << myname << line << endl;
  cout << myname << "Test data" << endl;

  cout << myname << line << endl;
  cout << myname << "Create histogram." << endl;
  TH1* ph = new TH1F("h1", "My histo", 100, 0, 100);
  for ( unsigned int ibin=1; ibin<=100; ++ ibin ) {
    ph->SetBinContent(ibin, ibin-1);
  }

  cout << myname << line << endl;
  cout << myname << "Create canvas and manipulator" << endl;
  TCanvas* pcan = new TCanvas;
  ph->Draw();
  TPadManipulator man(pcan);
  assert( man.pad() == pcan );
  assert( man.hist() == ph );

  cout << myname << line << endl;
  cout << myname << "Add axis." << endl;
  man.addAxis();

  cout << myname << line << endl;
  cout << myname << "Add vertical mod lines." << endl;
  man.addVerticalModLines(20, 10);

  pcan->Print("test_TPadManipulator.png");

  cout << myname << line << endl;
  cout << myname << "Done." << endl;
  return 0;
}

//**********************************************************************

int main() {
  return test_TPadManipulator();
}

//**********************************************************************
