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
#include "TCanvas.h"

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

  StickyCodeMetrics scm;

  cout << myname << line << endl;
  cout << myname << "Bin Counter evaluation." << endl;
  scm.evaluate(counts);
  scm.print();
  assert( scm.getHist() == nullptr );
  assert( ! scm.getSharedHist() );

  cout << myname << line << endl;
  cout << myname << "ADC samples evaluation." << endl;
  scm.evaluate(vals);
  scm.print();

  cout << myname << line << endl;
  cout << myname << "Histogram evaluation." << endl;
  scm.evaluate(ph);
  scm.print();

  cout << myname << line << endl;
  cout << myname << "Limited-range histogram evaluation." << endl;
  scm.evaluate(ph2);
  scm.print();

  cout << myname << line << endl;
  cout << myname << "Bin Counter data map." << endl;
  scm.getMetrics().print();

  TCanvas* pcan = new TCanvas;
  pcan->SetGridx();

  cout << myname << line << endl;
  cout << myname << "Evaluation creating histogram." << endl;
  StickyCodeMetrics scmh("hadctest", "ADC spectrum for test", 50, 10, 1, 20);
  assert( scmh.evaluate(counts) == 0 );
  scmh.print();
  scmh.getMetrics().print();
  assert( scmh.getHist() != nullptr );
  scmh.getHist()->Print();
  cout << myname << "Histogram integral: " << scmh.getHist()->Integral() << endl;
  scmh.getHist()->Draw();
  pcan->Print("test.png");
  assert( int(scmh.getHist()->Integral()+0.1) == int(nadc) );

  cout << myname << line << endl;
  cout << myname << "Evaluation creating wide histogram." << endl;
  StickyCodeMetrics scmhw("hadctest", "ADC spectrum for test", 100, 10, 1, 20);
  assert( scmhw.evaluate(counts) == 0 );
  scmhw.print();
  scmhw.getMetrics().print();
  assert( scmhw.getHist() != nullptr );
  scmhw.getHist()->Print();
  cout << myname << "Histogram integral: " << scmhw.getHist()->Integral() << endl;
  scmhw.getHist()->Draw();
  pcan->Print("testw.png");
  assert( int(scmhw.getHist()->Integral()+0.1) == int(nadc) );

  cout << myname << line << endl;
  cout << myname << "Evaluation creating narrow histogram." << endl;
  StickyCodeMetrics scmh2("hadctest", "ADC spectrum for test", 15, 5, 1, 20);
  assert( scmh2.evaluate(counts) == 0 );
  scmh2.print();
  scmh2.getMetrics().print();
  assert( scmh2.getHist() != nullptr );
  scmh2.getHist()->Print();
  cout << myname << "Histogram integral: " << scmh2.getHist()->Integral() << endl;
  cout << myname << "Histogram undrflow: " << scmh2.getHist()->GetBinContent(16) << endl;
  cout << myname << "Histogram overflow: " << scmh2.getHist()->GetBinContent(0) << endl;
  scmh2.getHist()->Draw();
  pcan->Print("test2.png");
  Index nadc2 = 0;
  Index nadc2u = 0;
  Index nadc2o = 0;
  for ( int iadc=110; iadc<125; ++iadc ) {
    if ( counts.find(iadc) != counts.end() ) nadc2 += counts[iadc];
  }
  for ( int iadc=100; iadc<110; ++iadc ) {
    if ( counts.find(iadc) != counts.end() ) nadc2u += counts[iadc];
  }
  for ( int iadc=125; iadc<135; ++iadc ) {
    if ( counts.find(iadc) != counts.end() ) nadc2o += counts[iadc];
  }
  cout << myname << "Expect integral: " << nadc2 << endl;
  cout << myname << "Expect undrflow: " << nadc2u << endl;
  cout << myname << "Expect overflow: " << nadc2o << endl;
  assert( int(scmh2.getHist()->Integral()+0.1) == int(nadc2) );

  delete pcan;

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
