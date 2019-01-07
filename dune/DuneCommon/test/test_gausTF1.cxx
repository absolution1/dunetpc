// test_gausTF1.cxx

// David Adams
// November 2018
//
// This is a test and demonstration for gausTF1.

#undef NDEBUG

#include "dune/DuneCommon/gausTF1.h"
#include <string>
#include <iostream>
#include <iomanip>
#include <vector>
#include <cassert>

using std::string;
using std::cout;
using std::endl;
using std::setw;
using std::vector;

using Index = unsigned int;

//**********************************************************************

int test_gausTF1() {
  const string myname = "test_gausTF1: ";
  cout << myname << "Starting test" << endl;
#ifdef NDEBUG
  cout << myname << "NDEBUG must be off." << endl;
  abort();
#endif
  string line = "-----------------------------";
  string scfg;

  cout << myname << line << endl;
  cout << "Testing TF1." << endl;

  cout << myname << line << endl;
  cout << "Build function." << endl;
  TF1* pf = gausTF1(20, 50, 10, "mygaus");
  pf->Print();
  return 0;
}

//**********************************************************************

int main() {
  return test_gausTF1();
}

//**********************************************************************
