// test_GausRmsFitter.cxx

// David Adams
// December 2018.
//
// This is a test of GausRmsFitter.

#undef NDEBUG

#include "../GausRmsFitter.h"
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

int test_GausRmsFitter(int dbg, double satfrac) {
  const string myname = "test_GausRmsFitter: ";
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
  GausRmsFitter gsf(4.0, 4.0, "ff");
  if ( dbg > 0 ) gsf.setLogLevel(dbg);
  assert( gsf.fit(ph, 50) == 0 );
  TF1* pffit = ph->GetFunction("ff");
  if ( pffit == nullptr ) {
    cout << myname << "  Fit function not found!" << endl;
    ph->GetListOfFunctions()->Print();
  } else {
    double x1, x2;
    pffit->GetRange(x1, x2);
    double sum = pffit->Integral(x1, x2);
    cout << myname << "      Fit name: " << pffit->GetName() << endl;
    cout << myname << "      Fit mean: " << pffit->GetParameter("Mean") << endl;
    cout << myname << "     Fit sigma: " << pffit->GetParameter("Sigma") << endl;
    cout << myname << "     Fit range: (" << x1 << ", " << x2 << ")" << endl;
    cout << myname << "  Fit integral: " << sum << endl;
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
  int dbg = 0;
  DoubleVector satfracs = {0.0, 0.1, 0.2, 0.5};
  if ( narg > 1 ) {
    string sarg = argc[1];
    if ( sarg == "-h" ) {
      cout << "Usage: " << argc[0] << " [DBG] [F1 F2 ...]" << endl;
      return 0;
    }
    istringstream ssin(sarg);
    ssin >> dbg;
  }
  if ( narg > 2 ) {
    string sarg = argc[2];
    satfracs.clear();
    for ( int iarg=2; iarg<narg; ++iarg ) {
      double satfrac = 0.0;
      istringstream ssin(argc[iarg]);
      ssin >> satfrac;
      satfracs.push_back(satfrac);
    }
  }
  int err = 0;
  for ( double satfrac : satfracs ) err += test_GausRmsFitter(dbg, satfrac);
  return err;
}

//**********************************************************************
