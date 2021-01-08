// test_CompactRealDftData.cxx
//
// David Adams
// April 2019
//
// Test CompactRealDftData.

#include "dune/DuneCommon/CompactRealDftData.h"
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

using Dft = CompactRealDftData<float>;
using Index = Dft::Index;
using FloatVector = Dft::FloatVector;

//**********************************************************************

int test_CompactRealDftData(bool useExistingFcl) {
  const string myname = "test_RealDftData: ";
#ifdef NDEBUG
  cout << myname << "NDEBUG must be off." << endl;
  abort();
#endif
  string line = "-----------------------------";
  const Dft::Norm mynorm(RealDftNormalization::Standard, RealDftNormalization::Unit);

  cout << myname << line << endl;
  cout << myname << "Create DFT." << endl;
  FloatVector amps = {   5.0, 4.0, 3.0, 2.0, 1.0 };
  FloatVector phas = {   0.0, 0.5, 1.0, 2.0, 4.0 };
  Index namp = amps.size();
  Index npha = phas.size();
  Index nsam = namp + npha - 1;
  assert( amps.size() == namp );
  assert( phas.size() == npha );
  cout << myname << "Sample length: " << nsam << endl;

  cout << myname << line << endl;
  cout << myname << "Create invalid DFT objects." << endl;
  FloatVector badvals = {11.1, 22.2};
  Dft dftbad1(mynorm, amps, badvals);
  assert ( ! dftbad1.size() );
  Dft dftbad2(mynorm, badvals, phas);
  assert ( ! dftbad2.size() );

  cout << myname << line << endl;
  cout << myname << "Create DFT object with data." << endl;
  Dft dft(mynorm, amps, phas);
  cout << myname << "  # samples: " << dft.size() << endl;
  cout << myname << "      Power: " << dft.power() << endl;
  assert( dft.normalization().isValid() );
  assert( dft.isValid() );
  assert( dft.size() );
  assert( dft.nCompact() == namp );
  assert( dft.nAmplitude() == namp );
  assert( dft.nPhase() == npha );
  assert( dft.nSample() == nsam );
  assert( dft.normalization().isStandard() );
  assert( ! dft.normalization().isConsistent() );
  assert( ! dft.normalization().isBin() );
  assert( dft.normalization().isUnit() );
  assert( ! dft.normalization().isPower() );

  // Copy data out.
  cout << myname << line << endl;
  cout << myname << "Copy data out." << endl;
  FloatVector amps1, phas1;
  dft.copyOut(amps1, phas1);
  assert( dft.isValid() );
  assert( dft.size() );
  assert( dft.nCompact() == namp );
  assert( dft.nAmplitude() == namp );
  assert( dft.nPhase() == npha );
  assert( dft.nSample() == nsam );
  assert( amps1 == amps );
  assert( phas1 == phas );

  // Move data out.
  cout << myname << line << endl;
  cout << myname << "Move data out." << endl;
  FloatVector amps2, phas2;
  dft.moveOut(amps2, phas2);
  assert( ! dft.isValid() );
  assert( dft.normalization().isValid() );
  assert( dft.size() == 0 );
  assert( dft.nCompact() == 0 );
  assert( dft.nAmplitude() == 0 );
  assert( dft.nPhase() == 0 );
  assert( dft.nSample() == 0 );
  assert( amps2 == amps );
  assert( phas2 == phas );

  cout << myname << line << endl;
  cout << myname << "Done." << endl;
  return 0;
}

//**********************************************************************

int main(int argc, char* argv[]) {
  bool useExistingFcl = false;
  if ( argc > 1 ) {
    string sarg(argv[1]);
    if ( sarg == "-h" ) {
      cout << "Usage: " << argv[0] << " [ARG]" << endl;
      cout << "  If ARG = true, existing FCL file is used." << endl;
      return 0;
    }
    useExistingFcl = sarg == "true" || sarg == "1";
  }
  return test_CompactRealDftData(useExistingFcl);
}

//**********************************************************************
