// test_SampleTailer.cxx

// David Adams
// April 2019
//
// This is a test and demonstration for SampleTailer.

#undef NDEBUG

#include "../SampleTailer.h"
#include <string>
#include <iostream>
#include <iomanip>
#include <cassert>
#include <vector>
#include "TH1F.h"
#include "TCanvas.h"
#include "TLine.h"

using std::string;
using std::cout;
using std::endl;
using std::setw;
using std::vector;

using FloatVector = SampleTailer::FloatVector;
using Index = unsigned int;
using IndexVector = std::vector<Index>;

//**********************************************************************

namespace {

string line() {
  return "-----------------------------";
}

void drawResults(const SampleTailer& sta, string httl, string fnam) {
  const string myname = "drawResults: ";
  cout << myname << line() << endl;
  cout << myname << "Draw data." << endl;
  TH1* phd = sta.dataHist();
  TH1* phs = sta.signalHist();
  TH1* pht = sta.tailHist();
  assert( phd != nullptr );
  assert( phs != nullptr );
  assert( pht != nullptr );
  TCanvas* pcan = new TCanvas("mycan", "mycan", 700, 500);
  pcan->SetGrid();
  phd->SetLineColor(1);
  phs->SetLineColor(2);
  pht->SetLineColor(9);
  TH1* pha = dynamic_cast<TH1*>(phd->Clone("hax"));
  pha->SetDirectory(nullptr);
  pha->GetYaxis()->SetRangeUser(-50, 220);
  pha->Clear();
  pha->SetTitle(httl.c_str());
  pha->Draw("");
  TLine* plin = new TLine(0, 0, sta.size(), 0);
  plin->Draw();
  phd->Draw("SAME");
  pht->Draw("SAME");
  phs->Draw("SAME");
  pcan->Print(fnam.c_str());
}

void printResults(const SampleTailer& sta) {
  const string myname = "printResults: ";
  cout << myname << line() << endl;
  cout << myname << "Display results." << endl;
  cout << myname << "    i       data    signal      tail" << endl;
  for ( Index isam=0; isam<sta.size(); ++isam ) {
    cout << myname << setw(5) << isam << ":";
    cout << setw(10) << std::fixed << std::setprecision(2) << sta.data(isam);
    cout << setw(10) << std::fixed << std::setprecision(2) << sta.signal(isam);
    cout << setw(10) << std::fixed << std::setprecision(2) << sta.tail(isam);
    cout << endl;
  }
}

}  // end unnamed namespace

//**********************************************************************

int test_SampleTailer() {
  const string myname = "test_SampleTailer: ";
  cout << myname << "Starting test" << endl;
#ifdef NDEBUG
  cout << myname << "NDEBUG must be off." << endl;
  abort();
#endif
  string scfg;

  cout << myname << "Create signals." << endl;
  Index nsam = 300;
  FloatVector pulse = {  0.1,  4.5, 15.2, 66.4, 94.3, 100.0, 96.5, 88.4, 72.6, 58.4,
                        42.3, 35.1, 26.0, 18.6, 12.6,   8.8,  6.9,  4.4,  2.0, 0.3 };
  Index npul = pulse.size();
  FloatVector sigs1(nsam, 0.0);
  IndexVector peakPoss = {10, 100, 115, 230};
  FloatVector peakAmps = {0.5, 2.0, 0.7, 1.0};
  Index npea = peakPoss.size();
  for ( Index ipea=0; ipea<npea; ++ipea ) {
    Index iposPeak = peakPoss[ipea];
    float norm = peakAmps[ipea];
    for ( Index ipul=0; ipul<npul; ++ipul ) {
      Index isam = iposPeak + ipul;
      if ( isam >= nsam ) break;
      sigs1[isam] += norm*pulse[ipul];
    }
  }

  cout << myname << "Create utility." << endl;
  SampleTailer sta(100.0);
  sta.setPedestal(5.0);
  sta.setTail0(-15.0);
  sta.setUnit("ADC count");
  cout << myname << "  decayTime: " << sta.decayTime() << endl;
  cout << myname << "       beta: " << sta.beta() << endl;
  cout << myname << "      alpha: " << sta.alpha() << endl;
  cout << myname << "   pedestal: " << sta.pedestal() << endl;
  cout << myname << "      tail0: " << sta.tail0() << endl;

  cout << myname << line() << endl;
  cout << myname << "Create data from signal." << endl;
  assert( sta.setSignal(sigs1) == 0 );
  assert( sta.data().size() == nsam );
  assert( sta.signal().size() == nsam );
  assert( sta.tail().size() == nsam );
  FloatVector dats1 = sta.data();
  sta.makeHists();
  printResults(sta);
  drawResults(sta, "Data, signal and tail from signal", "SampleTrailerFromSignal.png");

  cout << myname << line() << endl;
  cout << myname << "Create signal from data." << endl;
  sta.clear();
  assert( sta.data().size() == 0 );
  assert( sta.signal().size() == 0 );
  assert( sta.tail().size() == 0 );
  assert( sta.setData(dats1) == 0 );
  assert( sta.data().size() == nsam );
  assert( sta.signal().size() == nsam );
  assert( sta.tail().size() == nsam );
  FloatVector sigs2 = sta.signal();
  sta.makeHists();
  printResults(sta);
  drawResults(sta, "Data, signal and tail from data", "SampleTrailerFromData.png");

  cout << myname << line() << endl;
  cout << myname << "Compare new signal with original.";
  cout << myname << "    i       sig1      sig2      diff" << endl;
  assert( sigs1.size() == sigs2.size() );
  Index nbad = 0;
  for ( Index isam=0; isam<nsam; ++isam ) {
    float dif = sigs2[isam] - sigs1[isam];
    float adif = fabs(dif);
    cout << myname << setw(5) << isam << ":";
    cout << setw(10) << std::fixed << std::setprecision(2) << sigs1[isam];
    cout << setw(10) << std::fixed << std::setprecision(2) << sigs2[isam];
    cout << setw(10) << std::fixed << std::setprecision(2) << dif;
    if ( adif > 1.e-3 ) ++nbad;
    cout << endl;
  }
  assert( nbad == 0 );

  return 0;
}

//**********************************************************************

int main() {
  return test_SampleTailer();
}

//**********************************************************************
