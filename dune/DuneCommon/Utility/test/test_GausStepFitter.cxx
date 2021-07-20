// test_GausStepFitter.cxx

// David Adams
// December 2018.
//
// This is a test of GausStepFitter.

#undef NDEBUG

#include "../GausStepFitter.h"
#include "../gausTF1.h"
#include "TH1F.h"
#include "TF1.h"
#include "TCanvas.h"

#include <string>
#include <iostream>
#include <sstream>
#include <iomanip>
#include <cassert>
#include <vector>

using std::string;
using std::cout;
using std::endl;
using std::setw;
using std::ostringstream;
using std::istringstream;
using DoubleVector = std::vector<double>;

//**********************************************************************

int test_GausStepFitter(double satfrac) {
  const string myname = "test_GausStepFitter: ";
  cout << myname << "Starting test" << endl;
#ifdef NDEBUG
  cout << myname << "NDEBUG must be off." << endl;
  abort();
#endif
  string line = "-----------------------------";
  string scfg;

  cout << myname << line << endl;
  cout << myname << "Test data" << endl;
  ostringstream ssttl;
  ssttl << "Test data with frac " << satfrac << endl;
  TH1* ph = new TH1F("hdat", ssttl.str().c_str(), 100, 0, 100);
  ph->SetDirectory(nullptr);
  TF1* pf1 = gausTF1(1.0, 50.0, 4.0, "f1");
  cout << myname << "Created function " << pf1->GetName() << endl;
  TF1* pf2 = gausTF1(1.0, 30.0, 4.0, "f2");
  cout << myname << "Created function " << pf2->GetName() << endl;
  int npt = 1000;
  ph->FillRandom("f1", npt);
  ph->FillRandom("f2", satfrac*npt);
  ph->GetListOfFunctions()->Clear();
  cout << myname << "  Hist mean: " << ph->GetMean() << endl;
  cout << myname << "   Hist RMS: " << ph->GetRMS() << endl;
  delete pf1;
  delete pf2;

  cout << myname << line << endl;
  GausStepFitter gsf(50, 3.0, 20, "ff");
  assert( gsf.fit(ph) == 0 );
  TF1* pffit = ph->GetFunction("ff");
  if ( pffit == nullptr ) {
    cout << myname << "  Fit function not found!" << endl;
    ph->GetListOfFunctions()->Print();
  } else {
    cout << myname << "   Fit mean: " << pffit->GetParameter("Mean") << endl;
    cout << myname << "  Fit sigma: " << pffit->GetParameter("Sigma") << endl;
  }

  cout << myname << line << endl;
  TCanvas* pcan = new TCanvas;
  ph->Draw();
  ostringstream ssfnam;
  static int ncan = 0;
  ssfnam << "gausstep" << ncan << ".png";
  ++ncan;
  string fnam = ssfnam.str();
  cout << myname << "Writing " << fnam << endl;
  pcan->Print(fnam.c_str());
  delete pcan;

  return 0;
}

//**********************************************************************

int main(int narg, const char* argc[]) {
  DoubleVector satfracs = {0.0, 0.1, 0.2, 0.5};
  if ( narg > 1 ) {
    string sarg = argc[1];
    if ( sarg == "-h" ) {
      cout << "Usage: " << argc[0] << " [F1 F2 ...]" << endl;
      return 0;
    }
    satfracs.clear();
    for ( int iarg=1; iarg<narg; ++iarg ) {
      double satfrac = 0.0;
      istringstream ssin(argc[iarg]);
      ssin >> satfrac;
      satfracs.push_back(satfrac);
    }
  }
  int err = 0;
  for ( double satfrac : satfracs ) err += test_GausStepFitter(satfrac);
  return err;
}

//**********************************************************************
