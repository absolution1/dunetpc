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
#include "TROOT.h"

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
  double xmax = 10000;
  TH1* ph = new TH1F("h1", "My histo; X values; Y value [units]", 100, 0, xmax);
  ph->SetStats(0);
  ph->SetFillColor(2);
  for ( unsigned int ibin=1; ibin<=100; ++ ibin ) {
    ph->SetBinContent(ibin, ibin-1);
  }

  cout << myname << line << endl;
  cout << myname << "Create manipulator" << endl;
  TPadManipulator* pman = new TPadManipulator;
  TPadManipulator& man = *pman;
  assert( man.pad() == nullptr );
  assert( man.hist() == nullptr );

  cout << myname << line << endl;
  cout << myname << "Add histogram." << endl;
  man.add(ph, "H");
  assert( man.hist() != nullptr );

  cout << myname << line << endl;
  cout << myname << "Add axis." << endl;
  man.addAxis();

  cout << myname << line << endl;
  cout << myname << "Add vertical mod lines." << endl;
  man.addVerticalModLines(0.2*xmax, 0.1*xmax);

  cout << myname << line << endl;
  cout << myname << "Draw." << endl;
  assert( man.pad() == nullptr );
  man.draw();
  assert( man.pad() != nullptr );

  man.pad()->Print("test_TPadManipulator.png");

  cout << myname << line << endl;
  cout << myname << "Root canvas count: " << gROOT->GetListOfCanvases()->GetEntries() << endl;
  cout << myname << "Deleting manipulator." << endl;
  delete pman;
  cout << myname << "Root canvas count: " << gROOT->GetListOfCanvases()->GetEntries() << endl;

  cout << myname << line << endl;
  cout << myname << "Done." << endl;
  return 0;
}

//**********************************************************************

int main() {
  return test_TPadManipulator();
}

//**********************************************************************
