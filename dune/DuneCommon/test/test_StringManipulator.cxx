// test_StringManipulator.cxx

// David Adams
// April 2018
//
// This is a test and demonstration for StringManipulator.

#undef NDEBUG

#include "../StringManipulator.h"
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

int test_StringManipulator() {
  const string myname = "test_StringManipulator: ";
  cout << myname << "Starting test" << endl;
#ifdef NDEBUG
  cout << myname << "NDEBUG must be off." << endl;
  abort();
#endif
  string line = "-----------------------------";
  string scfg;

  char cfill = 'x';
  cout << myname << "Checking fill char." << endl;
  int vint = 123;
  cfill = StringManipulator::getFill<int>(vint);
  cout << myname << "  int " << vint << ": " << cfill << endl;
  assert( cfill == '0' );
  vint = -123;
  cfill = StringManipulator::getFill<int>(vint);
  cout << myname << "  int " << vint << ": " << cfill << endl;
  assert( cfill == '-' );
  float vflo = 1.23;
  cfill = StringManipulator::getFill<float>(vflo);
  cout << myname << "  float " << vflo << ": " << cfill << endl;
  assert( cfill == '_' );

  // Build raw and expected strings.
  vector<string> raws = {
    "abc*SUB*ghi", 
    "abc*POS*ghi", 
    "abc*NEG*ghi",
    "abc*SUB*ghi_*NEG*_*POS*.txt"
  };
  vector<string> exps = {
    "abcdefghi",
    "abc123ghi",
    "abc-123ghi",
    "abcdefghi_-123_123.txt"
  };
  vector<string> expws = {
    "abc__defghi",
    "abc00123ghi",
    "abc--123ghi",
    "abc__defghi_--123_00123.txt"
  };

  cout << myname << "Testing no width with " << raws.size() << " strings." << endl;
  Index w = 12;
  for ( Index istr=0; istr<raws.size(); ++istr ) {
    string raw = raws[istr];
    string exp = exps[istr];
    string val = raw;
    StringManipulator sman(val);
    sman.replace("*SUB*", "def");
    sman.replace("*POS*", 123);
    sman.replace("*NEG*", -123);
    cout << myname << "  "
         << setw(w) << raw << ": "
         << setw(w) << val << " ?= "
         << setw(w) << exp << endl;
    assert( raw != exp );
    assert( val == exp );
  }

  cout << myname << "Testing width 5 with " << raws.size() << " strings." << endl;
  for ( Index istr=0; istr<raws.size(); ++istr ) {
    string raw = raws[istr];
    string exp = expws[istr];
    string val = raw;
    StringManipulator sman(val);
    sman.replaceFixedWidth("*SUB*", "def", 5);
    sman.replaceFixedWidth("*POS*", 123, 5);
    sman.replaceFixedWidth("*NEG*", -123, 5);
    cout << myname << "  "
         << setw(w) << raw << ": "
         << setw(w) << val << " ?= "
         << setw(w) << exp << endl;
    assert( raw != exp );
    assert( val == exp );
  }

  cout << myname << "Done." << endl;
  return 0;
}

//**********************************************************************

int main() {
  return test_StringManipulator();
}

//**********************************************************************
