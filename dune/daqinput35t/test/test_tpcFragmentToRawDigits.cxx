// test_tpcFragmentToRawDigits.cxx

// David Adams
// September 2015
//
// This test demonstrates how to use art::Timestamp.

#undef NDEBUG

#include "../tpcFragmentToRawDigits.h"
#include "dune/DuneCommon/DuneTimeConverter.h"
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

int test_tpcFragmentToRawDigits() {
  const string myname = "test_tpcFragmentToRawDigits: ";
  cout << myname << "Starting test" << endl;
#ifdef NDEBUG
  cout << myname << "NDEBUG must be off." << endl;
  abort();
#endif
  string line = "-----------------------------";
  string scfg;

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
  cout << myname << "Test make_art_timestamp_from_nova_timestamp." << endl;
  art::Timestamp ts3 = DAQToOffline::make_art_timestamp_from_nova_timestamp(tnova);
  cout << myname << "Timestamp: " << setw(9) << ts3.timeHigh() << "."
       << setw(9) << ts3.timeLow() << " sec" << endl;
  assert( ts3.timeHigh() == ts2.timeHigh() );
  assert( ts3.timeLow() == ts2.timeLow() );

  cout << myname << line << endl;
  cout << myname << "Test old_make_art_timestamp_from_nova_timestamp." << endl;
  art::Timestamp ts4 = DAQToOffline::old_make_art_timestamp_from_nova_timestamp(tnova);
  cout << myname << "Timestamp: " << DuneTimeConverter::toString(ts4) << endl;
  assert( ts4.timeHigh() == 0 );
  assert( ts4.timeLow() == ts2.timeHigh() );

  cout << myname << line << endl;
  cout << myname << "Done." << endl;
  return 0;
}

//**********************************************************************

int main() {
  return test_tpcFragmentToRawDigits();
}

//**********************************************************************
