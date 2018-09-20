// test_RunData.cxx
//
// David Adams
// May 2018
//
// Test RunData.

#include <string>
#include <iostream>
#include <fstream>
#include "dune/DuneInterface/Data/RunData.h"

#undef NDEBUG
#include <cassert>

using std::string;
using std::cout;
using std::endl;

//**********************************************************************

int test_RunData(bool useExistingFcl =false) {
  const string myname = "test_RunData: ";
#ifdef NDEBUG
  cout << myname << "NDEBUG must be off." << endl;
  abort();
#endif
  string line = "-----------------------------";

  cout << myname << line << endl;
  cout << myname << "Creating bad data." << endl;
  RunData dat;
  cout << dat << endl;
  assert( ! dat.haveRun() );
  assert( ! dat.haveCryostat() );
  assert( ! dat.haveApas() );

  cout << myname << line << endl;
  cout << myname << "Make data good." << endl;
  dat.setRun(123);
  assert( dat.haveRun() );
  assert( dat.run() == 123 );
  cout << myname << "Add cryostat." << endl;
  assert( ! dat.haveCryostat() );
  dat.setCryostat("coldbox");
  assert( dat.haveCryostat() );
  assert( dat.cryostat() == "coldbox" );
  cout << myname << "Add APAs." << endl;
  assert( ! dat.haveApas() );
  RunData::IndexVector myapas = {1, 2};
  dat.setApas(myapas);
  assert( dat.haveApas() );
  assert( dat.apas().size() == 2 );
  assert( dat.apas() == myapas );
  cout << myname << "Add gain." << endl;
  assert( ! dat.haveGain() );
  dat.setGain(14);
  assert( dat.haveGain() );
  assert( dat.gain() == float(14.0) );
  cout << myname << "Add shaping." << endl;
  assert( ! dat.haveShaping() );
  dat.setShaping(2.0);
  assert( dat.haveShaping() );
  assert( dat.shaping() == float(2.0) );
  cout << myname << "Add leakage." << endl;
  assert( ! dat.haveLeakage() );
  dat.setLeakage(2.0);
  assert( dat.haveLeakage() );
  assert( dat.leakage() == float(2.0) );
  cout << myname << "Add HV frac." << endl;
  assert( ! dat.haveHvfrac() );
  dat.setHvfrac(0.5);
  assert( dat.haveHvfrac() );
  assert( dat.hvfrac() == float(0.5) );
  cout << myname << "Add pulser amplitude" << endl;
  assert( ! dat.havePulserAmplitude() );
  dat.setPulserAmplitude(16);
  assert( dat.havePulserAmplitude() );
  assert( dat.pulserAmplitude() == 16 );
  cout << myname << "Add pulser period" << endl;
  assert( ! dat.havePulserPeriod() );
  dat.setPulserPeriod(1234);
  assert( dat.havePulserPeriod() );
  assert( dat.pulserPeriod() == 1234 );

  cout << dat << endl;

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
  return test_RunData(useExistingFcl);
}

//**********************************************************************
