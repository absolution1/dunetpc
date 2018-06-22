// test_StickyCodeMetrics.cxx
//
// David Adams
// April 2017
//
// Test StickyCodeMetrics.

#include <string>
#include <iostream>
#include <fstream>
#include <vector>
#include "dune/DataPrep/Utility/StickyCodeMetrics.h"
#include "TH1F.h"

#undef NDEBUG
#include <cassert>

using std::string;
using std::cout;
using std::endl;
using std::ofstream;
using std::vector;

using Index = unsigned int;
using BinCounter = StickyCodeMetrics::BinCounter;

//**********************************************************************

int test_StickyCodeMetrics() {
  const string myname = "test_StickyCodeMetrics: ";
#ifdef NDEBUG
  cout << myname << "NDEBUG must be off." << endl;
  abort();
#endif
  string line = "-----------------------------";

  cout << myname << line << endl;

  cout << myname << line << endl;
  cout << myname << "Create data." << endl;
  StickyCodeMetrics::BinCounter counts;
  counts[111] =  3;
  counts[112] =  6;
  counts[113] = 20;
  counts[114] = 72;
  counts[115] = 85;
  counts[116] = 90;
  counts[117] = 86;
  counts[118] = 77;
  counts[119] = 35;
  counts[120] = 18;
  counts[121] =  7;
  counts[122] =  4;
  counts[123] =  1;
  counts[127] = 20;
  counts[128] = 40;
  counts[129] =  4;
  AdcCountVector vals;
  TH1F* ph = new TH1F("hraw", "Input ADC distribution; ADC count; # entries", 4096, 0, 4096);
  TH1F* ph2 = new TH1F("hraw2", "Input ADC distribution; ADC count; # entries", 30, 100, 130);
  Index nadc = 0;
  for ( BinCounter::value_type icnt : counts ) {
    AdcCount iadc = icnt.first;
    Index nval = icnt.second;
    nadc += nval;
    for ( Index ival=0; ival<nval; ++ival ) {
      vals.push_back(iadc);
      ph->Fill(iadc);
      ph2->Fill(iadc);
    }
  }
  cout << myname << "  Bin counter count: " << nadc << endl;
  cout << myname << "       Vector count: " << vals.size() << endl;
  cout << myname << "    Histogram count: " << ph->GetEntries() << endl;
  assert( vals.size() == nadc );
  assert( ph->GetEntries() == nadc );

  cout << myname << line << endl;
  cout << myname << "Bin Counter evaluation." << endl;
  StickyCodeMetrics scmbco(counts);
  scmbco.print();

  cout << myname << line << endl;
  cout << myname << "ADC samples evaluation." << endl;
  StickyCodeMetrics scmsam(vals);
  scmsam.print();

  cout << myname << line << endl;
  cout << myname << "Histogram evaluation." << endl;
  StickyCodeMetrics scmhst(ph);
  scmhst.print();

  cout << myname << line << endl;
  cout << myname << "Limited-range histogram evaluation." << endl;
  StickyCodeMetrics scmhs2(ph2);
  scmhs2.print();

  cout << myname << line << endl;
  cout << myname << "Bin Counter data map." << endl;
  scmbco.getMetrics().print();

  cout << myname << line << endl;
  cout << myname << "Done." << endl;
  return 0;
}

//**********************************************************************

int main(int argc, char* argv[]) {
  if ( argc > 1 ) {
    string sarg(argv[1]);
    if ( sarg == "-h" ) {
      cout << "Usage: " << argv[0] << endl;
      return 0;
    }
  }
  return test_StickyCodeMetrics();
}

//**********************************************************************
