// test_LineColors.cxx

// David Adams
// Feruary 2018
//
// This is a test and demonstration for LineColors.

#undef NDEBUG

#include "../LineColors.h"
#include <string>
#include <iostream>
#include <iomanip>
#include <cassert>
#include <vector>
#include "TH1F.h"
#include "TCanvas.h"

using std::string;
using std::cout;
using std::endl;
using std::setw;
using std::vector;

//**********************************************************************

int test_LineColors() {
  const string myname = "test_LineColors: ";
  cout << myname << "Starting test" << endl;
#ifdef NDEBUG
  cout << myname << "NDEBUG must be off." << endl;
  abort();
#endif
  string line = "-----------------------------";
  string scfg;

  cout << myname << line << endl;
  cout << myname << "Create histos." << endl;
  TH1* ph0 = new TH1F("hrate", "Rate; x; y", 50, 0, 10);
  ph0->SetStats(0);
  ph0->SetLineWidth(2);
  for ( int ibin=0; ibin<50; ++ibin ) {
    ph0->SetBinContent(ibin+1, 0.02*ibin);
  }
  vector<TH1*> hsts;
  LineColors lc;
  LineColors::Index ncol = 12;
  for ( LineColors::Index iidx=0; iidx<ncol; ++iidx ) {
    LineColors::ColorType icol = lc.color(iidx, ncol);
    cout << myname << iidx << ": " << icol << endl;
    TH1* ph = dynamic_cast<TH1*>(ph0->Clone());
    ph->Scale(iidx+1);
    ph->SetLineColor(icol);
    hsts.push_back(ph);
  }

  cout << myname << line << endl;
  cout << myname << "Draw histos." << endl;
  TCanvas* pcan = new TCanvas;
  ph0->SetMaximum(ncol + 0.2);
  ph0->Draw("axis");
  for ( TH1* ph : hsts ) ph->Draw("hist same");
  pcan->RedrawAxis();
  pcan->Print("test_LineColors.png");

  return 0;
}

//**********************************************************************

int main() {
  return test_LineColors();
}

//**********************************************************************
