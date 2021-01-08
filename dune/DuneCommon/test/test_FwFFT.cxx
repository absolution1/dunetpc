// test_FwFFT.cxx
//
// David Adams
// April 2019
//
// Test FwFFT.

#include "dune/DuneCommon/FwFFT.h"
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
using FloatVector = FwFFT::FloatVector;
using DFT = FwFFT::DFT;

//**********************************************************************

int test_FwFFT(Index ignorm, Index itnorm, int loglev, Index len) {
  const string myname = "test_FwFFT: ";
#ifdef NDEBUG
  cout << myname << "NDEBUG must be off." << endl;
  abort();
#endif
  string line = "-----------------------------";

  RealDftNormalization norm(ignorm, itnorm);
  cout << line << endl;
  cout << myname << " Global norm: " << ignorm << endl;
  cout << myname << "   Term norm: " << itnorm << endl;

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
  Index namp = (nsam+2)/2;
  Index npha = (nsam+1)/2;
  cout << myname << "        # samples: " << nsam << endl;
  cout << myname << "Expected amp size: " << namp << endl;
  cout << myname << "Expected pha size: " << npha << endl;

  cout << myname << line << endl;
  cout << myname << "Create empty DFT." << endl;
  DFT dft(norm);
  assert( dft.size() == 0 );

  cout << myname << line << endl;
  cout << myname << "Transform forward." << endl;
  FwFFT xf(50, 0);

  cout << myname << line << endl;
  cout << myname << "Transform forward." << endl;
  assert( xf.fftForward(sams, dft, loglev) == 0 );
  assert( sams.size() == nsam );
  cout << myname << "DFT size: " << dft.size() << endl;
  assert( dft.size() == nsam );
  assert( dft.nCompact() == namp );
  assert( dft.nAmplitude() == namp );
  assert( dft.nPhase() == npha );
  cout << "Frequency components:" << endl;
  for ( Index ifrq=0; ifrq<namp; ++ifrq ) {
    cout << myname << setw(4) << ifrq << ": " << dft.amplitude(ifrq);
    if ( ifrq < npha ) cout << " @ " << dft.phase(ifrq);
    cout << endl;
  }

  cout << myname << line << endl;
  cout << myname << "Check power." << endl;
  float pwr1 = 0.0;
  for ( float sam : sams ) pwr1 += sam*sam;
  float pwr2 = dft.power();
  float pwr3 = 0.0;
  float pwr4 = 0.0;
  for ( Index ifrq=0; ifrq<nsam; ++ifrq ) {
    float amp = dft.amplitude(ifrq);
    float xre = dft.real(ifrq);
    float xim = dft.imag(ifrq);
    if ( !dft.normalization().isPower() || ifrq < dft.nCompact() ) {
      pwr3 += amp*amp;
      pwr4 += xre*xre + xim*xim;
    }
  }
  double pfac = 1.0;
  if ( dft.normalization().isStandard() ) pfac = 1.0/nsam;
  if ( dft.normalization().isBin() ) pfac = nsam;
  pwr3 *= pfac;
  pwr4 *= pfac;
  cout << myname << "Tick power: " << pwr1 << endl;
  cout << myname << " DFT power: " << pwr2 << endl;
  cout << myname << " Amp power: " << pwr3 << endl;
  cout << myname << " ReI power: " << pwr4 << endl;
  assert( fabs(pwr2 - pwr1) < 1.e-5*pwr1 );

  cout << myname << line << endl;
  cout << myname << "Call inverse." << endl;
  FloatVector sams2, xres2, xims2;
  int rstat = xf.fftInverse(dft, sams2, loglev);
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
  Index ignorm = 1;
  Index itnorm = 1;
  int loglev = 0;
  Index len = 0;
  if ( argc > 1 ) {
    string sarg(argv[1]);
    if ( sarg == "-h" ) {
      cout << "Usage: " << argv[0] << " [ignorm [itnorm [loglev [LEN]]]]" << endl;
      return 0;
    }
    ignorm = std::stoi(sarg);
  }
  if ( argc > 2 ) {
    string sarg(argv[2]);
    itnorm = std::stoi(sarg);
  }
  if ( argc > 3 ) {
    string sarg(argv[3]);
    loglev = std::stoi(sarg);
  }
  if ( argc > 4 ) {
    string sarg(argv[4]);
    len = std::stoi(sarg);
  }
  return test_FwFFT(ignorm, itnorm, loglev, len);
}

//**********************************************************************
