// test_DataMap.cxx
//
// David Adams
// September 2017
//
// Test DataMap.

#include <string>
#include <iostream>
#include <fstream>
#include "dune/DuneInterface/Data/DataMap.h"

#undef NDEBUG
#include <cassert>

using std::string;
using std::cout;
using std::endl;

//**********************************************************************

int test_DataMap(bool useExistingFcl =false) {
  const string myname = "test_DataMap: ";
#ifdef NDEBUG
  cout << myname << "NDEBUG must be off." << endl;
  abort();
#endif
  string line = "-----------------------------";

  cout << myname << line << endl;
  cout << myname << "Creating failed result." << endl;
  DataMap resf(3);
  assert( resf );
  assert( resf != 0 );
  assert( resf == 3 );

  cout << myname << line << endl;
  cout << myname << "Creating success data." << endl;
  DataMap res;
  assert( ! res );

  cout << myname << line << endl;
  cout << myname << "Add int." << endl;
  assert( res.getIntMap().size() == 0 );
  assert( ! res.haveInt("myint") );
  res.setInt("myint", 123);
  assert( res.getIntMap().size() == 1 );
  assert( ! res.haveInt("nosuch") );
  assert( res.haveInt("myint") );
  assert( ! res.haveFloat("myint") );
  assert( res.getInt("nosuch") == 0 );
  assert( res.getInt("myint") == 123 );
  res.setInt("myint", 456);
  assert( res.getIntMap().size() == 1 );
  assert( res.getFloatMap().size() == 0 );
  assert( res.getHistMap().size() == 0 );
  assert( res.getInt("myint") == 456 );

  cout << myname << line << endl;
  cout << myname << "Add float." << endl;
  assert( res.getFloatMap().size() == 0 );
  assert( ! res.haveFloat("myflt") );
  res.setFloat("myflt", 1.23);
  assert( res.haveFloat("myflt") );
  assert( res.getFloat("myflt") == 1.23 );
  assert( res.getIntMap().size() == 1 );
  assert( res.getFloatMap().size() == 1 );
  assert( res.getHistMap().size() == 0 );

  cout << myname << line << endl;
  cout << myname << "Add histogram." << endl;
  TH1* ph = (TH1*) 0xdeadbeef;
  assert( res.getHistMap().size() == 0 );
  assert( ! res.haveHist("myhst") );
  res.setHist("myhst", ph);
  assert( res.haveHist("myhst") );
  assert( res.getHist("myhst") == ph );
  assert( res.getIntMap().size() == 1 );
  assert( res.getFloatMap().size() == 1 );
  assert( res.getHistMap().size() == 1 );

  cout << myname << line << endl;
  cout << myname << "Extend data map." << endl;
  DataMap res2;
  res2.setInt("myint", 135);
  res2.setInt("myint2", 578);
  res2.setFloat("myflt2", 5.78);
  res2.setHist("myhst2", ph);
  res.extend(res2);
  assert( res == 0 );
  assert( res.getIntMap().size() == 2 );
  assert( res.getFloatMap().size() == 2 );
  assert( res.getHistMap().size() == 2 );
  assert( res.getInt("myint") == 135 );
  assert( res.getInt("myint2") == 578 );
  assert( res.getFloat("myflt") == 1.23 );
  assert( res.getFloat("myflt2") == 5.78 );
  assert( res.getHist("myhst") == ph );
  assert( res.getHist("myhst2") == ph );


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
  return test_DataMap(useExistingFcl);
}

//**********************************************************************
