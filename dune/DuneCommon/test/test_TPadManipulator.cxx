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
#include "TGraph.h"
#include "TTimeStamp.h"
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
  TPadManipulator* pmantop = new TPadManipulator(700, 1000, 1, 2);

  TPadManipulator* pman = pmantop->man(0);
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
  cout << myname << "Add sloped lines." << endl;
  man.addSlopedLine(0.01, 5, 2);
  man.addSlopedLine(0.01, -5, 3);

  cout << myname << line << endl;
  cout << myname << "Add label." << endl;
  assert( ! man.getLabel().size() );
  man.setLabel("This is my #beta data for E_{tot} = m#timesC^{2}.");
  assert( man.getLabel().size() );

  cout << myname << line << endl;
  cout << myname << "Copy the pad and draw copy with time format." << endl;
  *pmantop->man(1) = *pmantop->man(0);
  pmantop->man(1)->setRangeX(0, 8000);
  pmantop->man(1)->setTitle("Zoom of above");
  pmantop->man(1)->setRangeX(2000, 8000);
  pmantop->man(1)->setRangeY(10, 90);
  pmantop->man(1)->setTimeFormatX("%H:%M");

  cout << myname << line << endl;
  cout << myname << "Draw." << endl;
  assert( man.pad() == nullptr );
  pmantop->draw();
  assert( man.pad() != nullptr );

  cout << myname << line << endl;
  cout << myname << "Draw." << endl;
  pmantop->print("test_TPadManipulator.png");

  cout << myname << line << endl;
  cout << myname << "Copy the full plot." << endl;
  TPadManipulator man2(*pmantop);

  cout << myname << line << endl;
  cout << myname << "Draw copy." << endl;
  man2.print("test_TPadManipulator2.png");

  cout << myname << line << endl;
  cout << myname << "Create time graph." << endl;
  TGraph* pg = new TGraph;
  pg->SetTitle("Clock");
  pg->GetXaxis()->SetTitle("March 6 daytime minutes");
  pg->GetYaxis()->SetTitle("Minute");
  pg->SetMarkerStyle(2);
  srand(12345);
  TTimeStamp ts0(2019, 3, 6, 6, 0, 0);
  for ( int i=0; i<12*60; ++i ) {
    double tsec = ts0 + 60*i;
    double min = fmod(i, 60);
    pg->SetPoint(i, tsec, min);
  }
  TPadManipulator mant(1400, 500);
  mant.add(pg, "P");
  mant.addAxis();
  double t0 = ts0.AsDouble();
  mant.setRangeX(t0 - 7200, t0+14*3600);
  mant.setRangeY(0, 61);
  mant.setTimeFormatX("%H:%M%F2019-06-06 00:00:00");
  mant.print("test_TPadManipulator-time.png");

  cout << myname << line << endl;
  cout << myname << "Root canvas count: " << gROOT->GetListOfCanvases()->GetEntries() << endl;
  cout << myname << "Deleting manipulator." << endl;
  delete pmantop;
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
