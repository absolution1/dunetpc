// test_DuneTimeConverter.cxx

// David Adams
// September 2015
//
// This test demonstrates how to use art::Timestamp.

#undef NDEBUG

#include "../DuneTimeConverter.h"
#include <string>
#include <iostream>
#include <iomanip>
#include <cassert>

using std::string;
using std::cout;
using std::endl;
using std::setw;
using std::uint32_t;
using std::uint64_t;

//**********************************************************************

int test_DuneTimeConverter() {
  const string myname = "test_DuneTimeConverter: ";
  cout << myname << "Starting test" << endl;
#ifdef NDEBUG
  cout << myname << "NDEBUG must be off." << endl;
  abort();
#endif
  string line = "-----------------------------";
  string scfg;

  cout << myname << line << endl;
  cout << myname << "Test data" << endl;
  assert( DuneTimeConverter::novaT0Sec() > 0 );
  assert( DuneTimeConverter::novaTicksPerUsec() > 0 );
  assert( DuneTimeConverter::novaTicksPerSec() > 0 );
  assert( DuneTimeConverter::novaTicksPerSec() == 1000000*DuneTimeConverter::novaTicksPerUsec() );

  cout << myname << line << endl;
  cout << myname << "Create test time" << endl;
  //uint32_t tsec = 2208988800;   // 2040
  uint32_t tsec = 1577836800;   // 2020
  //uint32_t tsec = 1483228800;   // 2017
  uint32_t trem = 123456789;

  cout << myname << line << endl;
  cout << myname << "Create timestamp." << endl;
  art::Timestamp ts1 = DuneTimeConverter::makeTimestamp(tsec, trem);
  cout << myname << DuneTimeConverter::toString(ts1) << " sec " << endl;
  assert( ts1.timeLow() == trem );
  assert( ts1.timeHigh() == tsec );

  cout << myname << line << endl;
  cout << myname << "Create nova time." << endl;
  uint64_t tnova = DuneTimeConverter::toNova(ts1);
  art::Timestamp ts2 = DuneTimeConverter::fromNova(tnova);
  cout << myname << "Nova time: " <<tnova << endl;
  cout << myname << "Timestamp: " << setw(9) << DuneTimeConverter::toString(ts2) << " sec" << endl;
  assert( ts2.timeHigh() == tsec );
  int64_t lodiff = ts2.timeLow();
  lodiff -= trem;
  int64_t tol = 1000/DuneTimeConverter::novaTicksPerUsec();;
  cout << myname << "Low tolerance: " << tol << endl;
  assert( abs(lodiff) < tol );

  cout << myname << line << endl;
  cout << myname << "Test creation from string." << endl;
  string stime = "1577836800.123456789";
  art::Timestamp ts5 = DuneTimeConverter::fromString(stime);
  cout << myname << "   String: " << stime << endl;
  cout << myname << "Timestamp: " << DuneTimeConverter::toString(ts5) << endl;
  assert( ts5.timeLow() == ts1.timeLow() );
  assert( ts5.timeHigh() == ts1.timeHigh() );
  assert( DuneTimeConverter::toString(ts5) == stime );

  cout << myname << line << endl;
  cout << myname << "Done." << endl;
  return 0;
}

//**********************************************************************

int main() {
  return test_DuneTimeConverter();
}

//**********************************************************************
