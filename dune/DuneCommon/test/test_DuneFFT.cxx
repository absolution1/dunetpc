// test_DuneFFT.cxx
//
// David Adams
// April 2019
//
// Test DuneFFT.

#include "dune/DuneCommon/DuneFFT.h"
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <iomanip>
#include <cmath>

#undef NDEBUG
#include <cassert>

using std::string;
using std::cout;
using std::endl;
using std::ostringstream;
using std::ofstream;
using std::vector;
using std::setw;
using std::fixed;

using Index = unsigned int;
using FloatVector = DuneFFT::FloatVector;

//**********************************************************************

int test_DuneFFT(bool useExistingFcl, Index len) {
  const string myname = "test_DuneFFT: ";
#ifdef NDEBUG
  cout << myname << "NDEBUG must be off." << endl;
  abort();
#endif
  string line = "-----------------------------";

  Index normOpt = 1;
  Index logLevel = 1;

  cout << myname << line << endl;
  cout << myname << "Create data." << endl;
  vector<float> sams = {   3.0,  6.0,  16.1,  28.6,  30.2,  27.7,  16.3,   9.6,  4.2, -1.0,
                          -2.3,  -4.2,  -9.2, -18.6, -21.9, -29.0, -24.3, -14.2, -5.0, -3.0};
  if ( len > 0 ) sams.resize(len, 0.0);
  Index nsam = sams.size();
  assert( sams.size() == nsam );
  cout << myname << "Sample length: " << nsam << endl;
  float samsum = 0.0;
  for ( float sam : sams ) samsum += sam;
  cout << myname << "Sample mean: " << samsum/nsam << endl;
  Index nmag = (nsam+2)/2;
  Index npha = (nsam+1)/2;
  cout << myname << "        # samples: " << nsam << endl;
  cout << myname << "Expected mag size: " << nmag << endl;
  cout << myname << "Expected mag size: " << nmag << endl;
  cout << myname << "Expected pha size: " << npha << endl;

  cout << myname << line << endl;
  cout << myname << "Transform forward." << endl;
  FloatVector xmgs, xphs, xres, xims;
  assert( DuneFFT::fftForward(normOpt, sams, xres, xims, xmgs, xphs, logLevel) == 0 );
  assert( sams.size() == nsam );
  assert( xmgs.size() == nmag );
  assert( xphs.size() == npha );
  assert( xres.size() == nsam );
  assert( xims.size() == nsam );

  cout << myname << line << endl;
  cout << myname << "Check power." << endl;
  float pwr1 = 0.0;
  for ( float sam : sams ) pwr1 += sam*sam;
  float pwr2 = 0.0;
  for ( float mag : xmgs ) pwr2 += mag*mag;
  float pwr3 = 0.0;
  for ( float x : xres ) pwr3 += x*x;
  for ( float x : xims ) pwr3 += x*x;
  cout << myname << "Tick power: " << pwr1 << endl;
  cout << myname << "Freq power: " << pwr2 << endl;
  cout << myname << "Frq2 power: " << pwr3 << endl;
  assert( fabs(pwr2 - pwr1) < 1.e-5*pwr1 );
  assert( fabs(pwr3 - pwr1) < 1.e-5*pwr1 );

  cout << myname << line << endl;
  cout << myname << "Call inverse." << endl;
  FloatVector sams2, xres2, xims2;
  int rstat = DuneFFT::fftInverse(normOpt, xmgs, xphs, xres, xims, sams2, logLevel);
  if ( rstat != 0 ) {
    cout << myname << "FFT invert returned " << rstat << endl;
    assert( false );
  }

  assert( sams2.size() == nsam );
  for ( Index isam=0; isam<nsam; ++isam ) {
    cout << setw(4) << isam << ":" << setw(10) << fixed << sams[isam]
         << " ?= " << setw(10) << fixed << sams2[isam] << endl;
    assert( fabs(sams2[isam] - sams[isam]) < 1.e-4 );
  }

  cout << myname << line << endl;
  cout << myname << "Done." << endl;
  return 0;
}

//**********************************************************************

int main(int argc, char* argv[]) {
  bool useExistingFcl = false;
  Index len = 0;
  if ( argc > 1 ) {
    string sarg(argv[1]);
    if ( sarg == "-h" ) {
      cout << "Usage: " << argv[0] << " [ARG] [LEN]" << endl;
      cout << "  If ARG = true, existing FCL file is used." << endl;
      return 0;
    }
    useExistingFcl = sarg == "true" || sarg == "1";
  }
  if ( argc > 2 ) {
    string sarg(argv[2]);
    len = std::stoi(sarg);
  }
  return test_DuneFFT(useExistingFcl, len);
}

//**********************************************************************
