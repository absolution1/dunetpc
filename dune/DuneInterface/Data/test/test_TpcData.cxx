// test_TpcData.cxx
//
// David Adams
// January 2021
//
// Test TpcData.

#include <string>
#include <iostream>
#include <fstream>
#include <vector>
#include "dune/DuneInterface/Data/TpcData.h"

#undef NDEBUG
#include <cassert>

using std::string;
using std::cout;
using std::endl;
using std::vector;

//**********************************************************************

int test_TpcData() {
  const string myname = "test_TpcData: ";
#ifdef NDEBUG
  cout << myname << "NDEBUG must be off." << endl;
  abort();
#endif
  string line = myname + "-----------------------------";

  cout << myname << "Checking empty data." << endl;
  TpcData tpd;
  assert( tpd.getTpcData("") == &tpd );
  assert( tpd.getTpcData(".") == &tpd );
  assert( tpd.getTpcData("sub1") == nullptr );
  assert( tpd.getParent() == nullptr );

  cout << myname << "Add sub1" << endl;
  TpcData* pdat1 = tpd.addTpcData("sub1");
  assert( pdat1 != nullptr );
  assert( tpd.addTpcData("sub1") == nullptr );
  assert( tpd.getTpcData("sub1") == pdat1 );
  assert( pdat1->getParent() == &tpd );
  const TpcData* pdat = &tpd;
  assert( pdat->getTpcData("sub1") == pdat1 );
  assert( pdat->getParent() == nullptr );
  const TpcData* pdat1c = pdat1;
  assert( pdat1c->getParent() == &tpd );

  cout << myname << "Add sub11" << endl;
  assert( tpd.getTpcData("sub1/sub11") == nullptr );
  TpcData* pdat11 = tpd.addTpcData("sub1/sub11");
  assert( pdat11 != nullptr );
  assert( tpd.addTpcData("sub1") == nullptr );
  assert( tpd.addTpcData("sub1/sub11") == nullptr );
  assert( tpd.getTpcData("sub1") == pdat1 );
  assert( tpd.getTpcData("sub1/sub11") == pdat11 );
  assert( pdat11->getParent() == pdat1 );
  const TpcData* pdat11c = pdat11;
  assert( pdat11c->getParent() == pdat1 );

  cout << myname << "Add sub12" << endl;
  assert( tpd.getTpcData("sub1/sub12") == nullptr );
  TpcData* pdat12 = tpd.addTpcData("sub1/sub12");
  assert( pdat12 != nullptr );
  assert( pdat12 != pdat11 );
  assert( tpd.addTpcData("sub1") == nullptr );
  assert( tpd.addTpcData("sub1/sub12") == nullptr );
  assert( tpd.getTpcData("sub1") == pdat1 );
  assert( tpd.getTpcData("sub1/sub12") == pdat12 );
  assert( pdat12->getParent() == pdat1 );
  const TpcData* pdat12c = pdat12;
  assert( pdat12c->getParent() == pdat1 );

  cout << line << endl;
  cout << myname << "Add ADC data." << endl;
  assert( tpd.getAdcData().size() == 0 );
  assert( pdat1->getAdcData().size() == 0 );
  assert( pdat11->getAdcData().size() == 0 );
  TpcData::AdcDataPtr pacm0 = tpd.createAdcData();
  assert( pacm0 );
  assert( tpd.getAdcData().size() == 1 );
  assert( pdat1->getAdcData().size() == 0 );
  TpcData::AdcDataPtr pacm1 = tpd.createAdcData();
  assert( pacm1 );
  assert( pacm1 != pacm0 );
  assert( tpd.getAdcData().size() == 2 );
  assert( tpd.getAdcData()[0] == pacm0 );
  assert( pdat1->getAdcData().size() == 0 );
  assert( pdat11->getAdcData().size() == 0 );
  assert( pdat1c->getAdcData().size() == 0 );
  assert( pdat11c->getAdcData().size() == 0 );

  cout << line << endl;
  cout << myname << "Check TPC data retrieval." << endl;
  assert(tpd.getTpcData("") == &tpd );
  assert(tpd.getTpcData(".") == &tpd );
  assert(tpd.getTpcData("sub3") == nullptr );
  assert(tpd.getTpcData("sub1") == pdat1 );
  assert(tpd.getTpcData("sub1/sub13") == nullptr );
  assert(tpd.getTpcData("sub1/sub11") == pdat11 );
  assert(tpd.getTpcData("sub1/sub12") == pdat12 );
  TpcData::TpcDataVector dats1;
  cout << myname << "...sub1" << endl;
  assert( tpd.getTpcData("sub1", dats1) == 0 );
  assert( dats1.size() == 1 );
  assert( dats1[0] == pdat1 );
  cout << myname << "...sub1" << endl;
  TpcData::TpcDataVector dats11;
  assert( tpd.getTpcData("sub1/sub11", dats11) == 0 );
  assert( dats11.size() == 1 );
  assert( dats11[0] == pdat11 );
  cout << myname << "...sub1/sub12" << endl;
  TpcData::TpcDataVector dats12;
  assert( tpd.getTpcData("sub1/sub12", dats12) == 0 );
  assert( dats12.size() == 1 );
  assert( dats12[0] == pdat12 );
  cout << myname << "...sub1/*" << endl;
  TpcData::TpcDataVector dats1x;
  assert( tpd.getTpcData("sub1/*", dats1x) == 0 );
  assert( dats1x.size() == 2 );
  assert( dats1x[0] == pdat11 );
  assert( dats1x[1] == pdat12 );

  cout << line << endl;
  cout << myname << "All tests passed." << endl;
  return 0;
}

//**********************************************************************

int main(int argc, char* argv[]) {
  if ( argc > 1 ) {
    cout << "Usage: " << argv[0] << " [ARG]" << endl;
    return 0;
  }
  return test_TpcData();
}

//**********************************************************************
