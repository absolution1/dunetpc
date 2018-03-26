// test_RootPalette.cxx

// David Adams
// March 2018
//
// This is a test and demonstration for RootPalette.

#undef NDEBUG

#include "dune/DuneCommon/RootPalette.h"
#include <string>
#include <iostream>
#include <iomanip>
#include <cassert>
#include <vector>
#include "TH2F.h"
#include "TCanvas.h"
#include "TColor.h"  // To get Root palette indices
#include "TStyle.h"
#include "TRandom.h"

using std::string;
using std::cout;
using std::endl;
using std::setw;
using std::vector;

//**********************************************************************

int test_RootPalette() {
  const string myname = "test_RootPalette: ";
  cout << myname << "Starting test" << endl;
#ifdef NDEBUG
  cout << myname << "NDEBUG must be off." << endl;
  abort();
#endif
  string line = "-----------------------------";
  string scfg;

  cout << myname << line << endl;
  cout << myname << "Test data" << endl;

  vector<double> peak1 = { 1.5, 6.0, 22.3, 50.3, 74.1, 92.3, 100.0, 95.1, 89.2, 64.3, 40.2, 18.0, 4.8, 0.6 };
  vector<double> peak2 = { 1.5, 6.0, 22.3, 50.3, 74.1, 92.3, 100.0, 95.1, 89.2, 64.3, 40.2, 18.0, 4.8,
                          -12.0, -26.1, -51.2, -93.7, -95.8, -80.6, -64.3, -49.3, -35.0, -22.1, -15.3,
                          -8.6, -5.1, -2.8};
  vector<double> wf1(100, 0.0);
  vector<double> wf2(100, 0.0);
  for ( unsigned int ipk=0; ipk<peak1.size(); ++ipk ) {
    wf1[ipk] = 0.80*peak1[ipk];
  }
  for ( unsigned int ipk=0; ipk<peak2.size(); ++ipk ) {
    wf2[ipk] = 0.80*peak2[ipk];
  }

  cout << myname << line << endl;
  cout << myname << "Create histogram." << endl;
  TH1* ph1 = new TH2F("h2", "My histo; Tick; Channel", 100, 0, 100, 50, 0, 50);
  TH1* ph2 = new TH2F("h2", "My histo; Tick; Channel", 100, 0, 100, 50, 0, 50);
  for ( TH1* ph : {ph1, ph2} ) {
    ph->SetStats(0);
    ph->SetFillColor(2);
    ph->SetContour(40);
    ph->SetMinimum(-100);
    ph->SetMaximum(100);
  }
  ph1->SetMinimum(0.0);
  ph1->SetContour(20);
  for ( unsigned int jbin=1; jbin<=100; ++jbin ) {
    for ( unsigned int ibin=1; ibin<=100; ++ibin ) {
      int iwf = ( 190 + ibin - jbin)%100;
      double z1 = wf1[iwf];
      z1 += gRandom->Gaus(0.0, 5.0);
      ph1->SetBinContent(ibin, jbin, z1);
      double z2 = wf2[iwf];
      z2 += gRandom->Gaus(0.0, 5.0);
      ph2->SetBinContent(ibin, jbin, z2);
    }
  }

  TCanvas* pcan = new TCanvas;
  pcan->SetGridx();
  pcan->SetGridy();
  int ncol = 0;
  TH1* ph = ph2;

  cout << myname << line << endl;
  cout << myname << "Draw with default palette." << endl;
  ph->SetTitle("Default palette");
  ncol = gStyle->GetNumberOfColors();
  cout << myname << "  Color count: " << ncol << endl;
  ph->Draw("colz");
  pcan->Print("test_RootPalette_default.png");

  RootPalette oldPalette;

  cout << myname << line << endl;
  cout << myname << "Draw with palette 1010." << endl;
  ph = ph1;
  assert( RootPalette::set(1010) );
  ncol = gStyle->GetNumberOfColors();
  cout << myname << "  Color count: " << ncol << endl;
  ph->SetTitle("Palette 1010");
  ph->Draw("colz");
  pcan->Print("test_RootPalette_1010.png");

  cout << myname << line << endl;
  cout << myname << "Draw with palette 1016." << endl;
  ph = ph1;
  assert( RootPalette::set(1016) );
  ncol = gStyle->GetNumberOfColors();
  cout << myname << "  Color count: " << ncol << endl;
  ph->SetTitle("Palette 1016");
  ph->Draw("colz");
  pcan->Print("test_RootPalette_1016.png");

  cout << myname << line << endl;
  cout << myname << "Draw with palette 1017." << endl;
  ph = ph1;
  assert( RootPalette::set(1017) );
  ncol = gStyle->GetNumberOfColors();
  cout << myname << "  Color count: " << ncol << endl;
  ph->SetTitle("Palette 1017");
  ph->Draw("colz");
  pcan->Print("test_RootPalette_1017.png");

  cout << myname << line << endl;
  cout << myname << "Draw with palette 1019." << endl;
  ph = ph1;
  assert( RootPalette::set(1019) );
  ncol = gStyle->GetNumberOfColors();
  cout << myname << "  Color count: " << ncol << endl;
  ph->SetTitle("Palette 1019");
  ph->Draw("colz");
  pcan->Print("test_RootPalette_1019.png");

  cout << myname << line << endl;
  cout << myname << "Draw with palette 1020." << endl;
  ph = ph2;
  assert( RootPalette::set(1020) );
  ncol = gStyle->GetNumberOfColors();
  cout << myname << "  Color count: " << ncol << endl;
  ph->SetTitle("Palette 1020");
  ph->Draw("colz");
  pcan->Print("test_RootPalette_1020.png");

  cout << myname << line << endl;
  cout << myname << "Draw with palette 1026." << endl;
  ph = ph2;
  assert( RootPalette::set(1026) );
  ncol = gStyle->GetNumberOfColors();
  cout << myname << "  Color count: " << ncol << endl;
  ph->SetTitle("Palette 1026");
  ph->Draw("colz");
  pcan->Print("test_RootPalette_1026.png");

  cout << myname << line << endl;
  cout << myname << "Draw with palette 1027." << endl;
  ph = ph2;
  assert( RootPalette::set(1027) );
  ncol = gStyle->GetNumberOfColors();
  cout << myname << "  Color count: " << ncol << endl;
  ph->SetTitle("Palette 1027");
  ph->Draw("colz");
  pcan->Print("test_RootPalette_1027.png");

  cout << myname << line << endl;
  cout << myname << "Draw with palette 1029." << endl;
  ph = ph2;
  assert( RootPalette::set(1029) );
  ncol = gStyle->GetNumberOfColors();
  cout << myname << "  Color count: " << ncol << endl;
  ph->SetTitle("Palette 1029");
  ph->Draw("colz");
  pcan->Print("test_RootPalette_1029.png");

  cout << myname << line << endl;
  cout << myname << "Draw with rainbow palette." << endl;
  ph = ph2;
  assert( RootPalette::set(kRainBow) );
  ncol = gStyle->GetNumberOfColors();
  cout << myname << "  Color count: " << ncol << endl;
  ph->SetTitle("Palette kRainBow");
  ph->Draw("colz");
  pcan->Print("test_RootPalette_rainbow.png");

  cout << myname << line << endl;
  cout << myname << "Draw with palette zero." << endl;
  assert( RootPalette::set(0) );
  ncol = gStyle->GetNumberOfColors();
  cout << myname << "  Color count: " << ncol << endl;
  ph->SetTitle("Palette zero (default)");
  ph->Draw("colz");
  pcan->Print("test_RootPalette_zero.png");

  cout << myname << line << endl;
  cout << myname << "Draw with invalid palette 720." << endl;
  assert( RootPalette::set(720) );
  ncol = gStyle->GetNumberOfColors();
  cout << myname << "  Color count: " << ncol << endl;
  ph->SetTitle("Palette 720 (invalid)");
  ph->Draw("colz");
  pcan->Print("test_RootPalette_0720.png");

  cout << myname << line << endl;
  cout << myname << "Draw with invalid palette 12." << endl;
  ph = ph1;
  assert( RootPalette::set(12) );
  ncol = gStyle->GetNumberOfColors();
  cout << myname << "  Color count: " << ncol << endl;
  ph->SetTitle("Palette 12 (invalid)");
  ph->Draw("colz");
  pcan->Print("test_RootPalette_0012.png");

  cout << myname << line << endl;
  cout << myname << "Draw with kGistEarthlette inverted." << endl;
  ph = ph1;
  assert( RootPalette::set(-kGistEarth) );
  ncol = gStyle->GetNumberOfColors();
  cout << myname << "  Color count: " << ncol << endl;
  ph->SetTitle("Palette kGistEarth inverted)");
  ph->Draw("colz");
  pcan->Print("test_RootPalette_gistearth.png");

  cout << myname << line << endl;
  cout << myname << "Draw with back body palette." << endl;
  ph = ph2;
  assert( RootPalette::set(kBlackBody) );
  ncol = gStyle->GetNumberOfColors();
  cout << myname << "  Color count: " << ncol << endl;
  ph->SetTitle("Palette kBlackBody");
  ph->Draw("colz");
  pcan->Print("test_RootPalette_blackbody.png");

  cout << myname << line << endl;
  cout << myname << "Draw with inverted deep sea palette." << endl;
  ph = ph1;
  assert( RootPalette::set(kDeepSea) );
  TColor::InvertPalette();
  ncol = gStyle->GetNumberOfColors();
  cout << myname << "  Color count: " << ncol << endl;
  ph->SetTitle("Palette kDeepSea inverted");
  ph->Draw("colz");
  pcan->Print("test_RootPalette_deepsea_inv.png");

  cout << myname << line << endl;
  cout << myname << "Draw with inverted cherry palette." << endl;
  ph = ph1;
  assert( RootPalette::set(-kCherry) );
  ncol = gStyle->GetNumberOfColors();
  cout << myname << "  Color count: " << ncol << endl;
  ph->SetTitle("Palette kCherry inverted");
  ph->Draw("colz");
  pcan->Print("test_RootPalette_cherry_inv.png");

  cout << myname << line << endl;
  cout << myname << "Draw with default palette positive." << endl;
  ph = ph1;
  assert( RootPalette::defaultPalette()->setRootPalette() == 0 );
  ncol = gStyle->GetNumberOfColors();
  cout << myname << "  Color count: " << ncol << endl;
  ph->SetTitle("Palette default pos");
  ph->Draw("colz");
  pcan->Print("test_RootPalette_default_pos.png");

  cout << myname << line << endl;
  cout << myname << "Done." << endl;
  return 0;
}

//**********************************************************************

int main() {
  return test_RootPalette();
}

//**********************************************************************
