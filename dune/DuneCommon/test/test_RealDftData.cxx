// test_RealDftData.cxx
//
// David Adams
// April 2019
//
// Test RealDftData.

#include "dune/DuneCommon/RealDftDataCompact.h"
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

using DftSU = RealDftDataCompact<float, RealDftNormalization::Standard, RealDftNormalization::Unit>;
using DftCP = RealDftDataCompact<float, RealDftNormalization::Consistent, RealDftNormalization::Power>;

using Index = DftSU::Index;
using FloatVector = DftSU::FloatVector;

//**********************************************************************

int test_RealDftData(bool useExistingFcl) {
  const string myname = "test_RealDftData: ";
#ifdef NDEBUG
  cout << myname << "NDEBUG must be off." << endl;
  abort();
#endif
  string line = "-----------------------------";

  cout << myname << line << endl;
  cout << myname << "Create DFT." << endl;
  vector<float> amps = {   5.0, 4.0, 3.0, 2.0, 1.0 };
  vector<float> phas = {   0.0, 0.5, 1.0, 2.0, 4.0 };
  Index namp = amps.size();
  Index npha = phas.size();
  Index nsam = namp + npha - 1;
  assert( amps.size() == namp );
  assert( phas.size() == npha );
  cout << myname << "Sample length: " << nsam << endl;

  cout << myname << line << endl;
  cout << myname << "Create invalid DFT objects." << endl;
  FloatVector badvals = {11.1, 22.2};
  DftSU dftbad1(amps, badvals);
  assert ( ! dftbad1.size() );
  DftSU dftbad2(badvals, phas);
  assert ( ! dftbad2.size() );

  cout << myname << line << endl;
  cout << myname << "Create DFT object with data." << endl;
  DftSU dft(amps, phas);
  cout << myname << "  # samples: " << dft.size() << endl;
  cout << myname << "      Power: " << dft.power() << endl;
  assert( dft.hasValidNormalization() );
  assert( dft.isValid() );
  assert( dft.size() );
  assert( dft.nCompact() == namp );
  assert( dft.nAmplitude() == namp );
  assert( dft.nPhase() == npha );
  assert( dft.nSample() == nsam );
  assert( dft.isStandard() );
  assert( ! dft.isConsistent() );
  assert( ! dft.isBin() );
  assert( dft.isUnit() );
  assert( ! dft.isPower() );

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
  return test_RealDftData(useExistingFcl);
}

//**********************************************************************
