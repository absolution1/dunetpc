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
  assert( res.getFloatVectorMap().size() == 0 );
  assert( res.getHistMap().size() == 0 );
  assert( res.getInt("myint") == 456 );
  res.print();

  cout << myname << line << endl;
  cout << myname << "Add int vector." << endl;
  assert( res.getIntVectorMap().size() == 0 );
  assert( ! res.haveIntVector("myintvec") );
  DataMap::IntVector myintvec;
  myintvec.push_back(1);
  myintvec.push_back(2);
  myintvec.push_back(3);
  res.setIntVector("myintvec", myintvec);
  assert( res.getIntVectorMap().size() == 1 );
  assert( ! res.haveIntVector("nosuch") );
  assert( res.haveIntVector("myintvec") );
  assert( ! res.haveFloat("myint") );
  res.setInt("myint", 456);
  assert( res.getIntMap().size() == 1 );
  assert( res.getFloatMap().size() == 0 );
  assert( res.getFloatVectorMap().size() == 0 );
  assert( res.getHistMap().size() == 0 );
  assert( res.getIntVector("myintvec") == myintvec );
  res.print();

  cout << myname << line << endl;
  cout << myname << "Add float." << endl;
  assert( res.getFloatMap().size() == 0 );
  assert( ! res.haveFloat("myflt") );
  float fltval = 1.23;
  res.setFloat("myflt", fltval);
  res.print();
  assert( res.haveFloat("myflt") );
  assert( res.getFloat("myflt") == fltval );
  assert( res.getIntMap().size() == 1 );
  assert( res.getFloatMap().size() == 1 );
  assert( res.getFloatVectorMap().size() == 0 );
  assert( res.getHistMap().size() == 0 );
  assert( res.getHistVectorMap().size() == 0 );

  cout << myname << line << endl;
  cout << myname << line << endl;
  cout << myname << "Add float vector." << endl;
  std::vector<float> flts;
  flts.push_back(1.234);
  flts.push_back(12.34);
  flts.push_back(123.4);
  assert( res.getFloatVectorMap().size() == 0 );
  assert( ! res.haveFloatVector("myflt") );
  res.setFloatVector("myfltvec", flts);
  assert( res.haveFloatVector("myfltvec") );
  assert( res.getFloatVector("myfltvec") == flts );
  assert( res.getIntMap().size() == 1 );
  assert( res.getFloatMap().size() == 1 );
  assert( res.getFloatVectorMap().size() == 1 );
  assert( res.getHistMap().size() == 0 );
  assert( res.getHistVectorMap().size() == 0 );
  res.print();

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
  assert( res.getFloatVectorMap().size() == 1 );
  assert( res.getHistMap().size() == 1 );
  assert( res.getHistVectorMap().size() == 0 );
  res.print();

  cout << myname << line << endl;
  cout << myname << "Add histogram vector." << endl;
  std::vector<TH1*> hsts = { (TH1*) 0xdeadbee1, (TH1*) 0xdeadbee2, (TH1*) 0xdeadbee3 };
  assert( res.getHistVectorMap().size() == 0 );
  assert( ! res.haveHistVector("myhstvec") );
  res.setHistVector("myhstvec", hsts);
  assert( res.haveHistVector("myhstvec") );
  assert( res.getHistVector("myhstvec") == hsts );
  assert( ! res.haveHist("nosuch") );
  assert( res.getHistVector("nosuch").size() == 0 );
  assert( res.getIntMap().size() == 1 );
  assert( res.getFloatMap().size() == 1 );
  assert( res.getFloatVectorMap().size() == 1 );
  assert( res.getHistMap().size() == 1 );
  assert( res.getHistVectorMap().size() == 1 );
  assert( res.getGraphMap().size() == 0 );
  res.print();

  cout << myname << line << endl;
  cout << myname << "Add graph." << endl;
  TGraph* pg = new TGraph;
  assert( ! res.haveGraph("mygrf") );
  res.setGraph("mygrf", pg);
  assert( res.haveGraph("mygrf") );
  assert( res.getGraph("mygrf") == pg );
  assert( res.getIntMap().size() == 1 );
  assert( res.getFloatMap().size() == 1 );
  assert( res.getFloatVectorMap().size() == 1 );
  assert( res.getHistMap().size() == 1 );
  assert( res.getHistVectorMap().size() == 1 );
  assert( res.getGraphMap().size() == 1 );
  res.print();

  cout << myname << line << endl;
  cout << myname << line << endl;
  cout << myname << "Extend data map." << endl;
  DataMap res2;
  TH1* ph2 = (TH1*) 0xdeadbee5;
  std::vector<TH1*> hsts2 = { (TH1*) 0xdeadbee6, (TH1*) 0xdeadbee7 };
  res2.setInt("myint", 135);
  res2.setInt("myint2", 578);
  float fltval2 = 5.78;
  res2.setFloat("myflt2", fltval2);
  res2.setHist("myhst2", ph2);
  res2.setHistVector("myhstvec2", hsts2);
  res.extend(res2);
  assert( res == 0 );
  assert( res.getIntMap().size() == 2 );
  assert( res.getFloatMap().size() == 2 );
  assert( res.getFloatVectorMap().size() == 1 );
  assert( res.getHistMap().size() == 2 );
  assert( res.getInt("myint") == 135 );
  assert( res.getInt("myint2") == 578 );
  assert( res.getFloat("myflt") == fltval );
  assert( res.getFloat("myflt2") == fltval2 );
  assert( res.getHist("myhst") == ph );
  assert( res.getHist("myhst2") == ph2 );
  assert( res.getHistVector("myhstvec") == hsts );
  assert( res.getHistVector("myhstvec2") == hsts2 );
  res2.print();

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
