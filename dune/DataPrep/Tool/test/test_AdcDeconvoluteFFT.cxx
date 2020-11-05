// test_AdcDeconvoluteFFT.cxx
//
// David Adams
// April 2019
//
// Test AdcDeconvoluteFFT.

#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <iomanip>
#include "dune/DuneInterface/Tool/AdcChannelTool.h"
#include "dune/ArtSupport/DuneToolManager.h"
#include "dune/DuneCommon/TPadManipulator.h"
#include <TRandom.h>
#include <TH1F.h>
#include <TCanvas.h>
#include "TLegend.h"

#undef NDEBUG
#include <cassert>

using std::string;
using std::cout;
using std::endl;
using std::ostringstream;
using std::istringstream;
using std::ofstream;
using fhicl::ParameterSet;
using std::vector;
using std::setw;
using std::fixed;

using Index = unsigned int;
using Name = std::string;

//**********************************************************************
namespace {

class Plot {
public:
  TPadManipulator man;
  TLegend* pleg;
  Index nhst = 0;
  Plot() : man(1000, 500) {
    pleg = man.addLegend(0.70, 0.60, 0.90, 0.85);
  }
  int addHist(const AdcSignalVector& vals, Name nam, int col, int lwid, int lsty) {
    Index nsam = vals.size();
    Name hnam = "h" + nam;
    Name httl = "Samples for " + nam + "; Tick; Signal";
    TH1* ph = new TH1F(hnam.c_str(), httl.c_str(), nsam, 0, nsam);
    ph->SetStats(0);
    ph->SetLineWidth(lwid);
    ph->SetLineStyle(lsty);
    ph->SetLineColor(col);
    for ( Index isam=0; isam<nsam; ++isam ) ph->Fill(isam+1, vals[isam]);
    Name dopt = nhst ? "HIST SAME" : "HIST";
    man.add(ph, dopt);
    if ( nhst == 0 ) {
      man.setGridY();
    }
    pleg->AddEntry(ph, nam.c_str(), "l");
    ++nhst;
    //pcan->cd();
    //ph->Draw(dopt.c_str());
    return 0;
  }
  int print(Name fname, float ymin =0.0, float ymax =0.0) {
    if ( ymax > ymin ) man.setRangeY(ymin, ymax);
    man.draw();
    man.print(fname.c_str());
    return 0;
  }
};

}  // end unnamed namepsace
//**********************************************************************

int test_AdcDeconvoluteFFT(bool useExistingFcl, float noiseLev, float sigmaFilter, Index nsam) {
  const string mypre = "test_AdcDeconvoluteFFT";
  const string myname = mypre + ": ";
#ifdef NDEBUG
  cout << myname << "NDEBUG must be off." << endl;
  abort();
#endif
  string line = "-----------------------------";

  cout << myname << line << endl;
  string fclfile = "test_AdcDeconvoluteFFT.fcl";
  AdcSignalVector res = {0.03, 0.31, 0.77, 0.97, 1.00, 0.96, 0.92, 0.84, 0.77, 0.66,
                         0.53, 0.39, 0.31, 0.26, 0.21, 0.17, 0.14, 0.115, 0.08, 0.06,
                         0.043, 0.029, 0.022, 0.015, 0.009, 0.005, 0.003, 0.001};
  float sum = 0.0;
  for ( float val : res ) sum += val;
  for ( float& val : res ) val /= sum;
  if ( ! useExistingFcl ) {
    cout << myname << "Creating top-level FCL." << endl;
    ofstream fout(fclfile.c_str());
    fout << "tools: {" << endl;
    fout << "  mydco: {" << endl;
    fout << "           tool_type: AdcDeconvoluteFFT" << endl;
    fout << "            LogLevel: 2" << endl;
    fout << "     ResponseVectors: [[";
    for ( Index ival=0; ival<res.size(); ++ival ) {
      if ( ival ) fout << ", ";
      fout << res[ival];
    }
    fout << "]]" << endl;
    fout << "              Action: 1" << endl;
    fout << "   GausFilterSigmas: [" << sigmaFilter << "]" << endl;
    fout << "       IndexMapTool: \"\"" << endl;
    fout << "  }" << endl;
    fout << "}" << endl;
    fout << "tools.mycondir: @local::tools.mydco" << endl;
    fout << "tools.mycondir.Action: 4" << endl;
    fout << "tools.myconfft: @local::tools.mydco" << endl;
    fout << "tools.myconfft.Action: 2" << endl;
    fout.close();
  } else {
    cout << myname << "Using existing top-level FCL." << endl;
  }

  cout << myname << line << endl;
  cout << myname << "Fetching tool manager." << endl;
  DuneToolManager* ptm = DuneToolManager::instance(fclfile);
  assert ( ptm != nullptr );
  DuneToolManager& tm = *ptm;
  tm.print();
  assert( tm.toolNames().size() == 3 );

  cout << myname << line << endl;
  cout << myname << "Fetching tool." << endl;
  auto pcondir = tm.getPrivate<AdcChannelTool>("mycondir");
  assert( pcondir != nullptr );
  auto pconfft = tm.getPrivate<AdcChannelTool>("myconfft");
  assert( pconfft != nullptr );
  auto pdco = tm.getPrivate<AdcChannelTool>("mydco");
  assert( pdco != nullptr );

  cout << myname << line << endl;
  cout << myname << "Create true samples." << endl;
  AdcSignalVector samsTru(nsam, 0.0);
  samsTru[20] = 1000.0;

  cout << myname << line << endl;
  cout << myname << "Create data by direct convolution with response." << endl;
  AdcChannelData acd;
  acd.run = 123;
  acd.event = 456;
  acd.channel = 789;
  acd.samples = samsTru;
  assert( acd.samples.size() == nsam );
  DataMap ret = pcondir->update(acd);
  assert( ret == 0 );
  cout << "      # samples: " << acd.samples.size() << endl;
  cout << "     # DFT amps: " << acd.dftmags.size() << endl;
  cout << "   # DFT phases: " << acd.dftphases.size() << endl;
  assert( acd.samples.size() == nsam );
  AdcSignalVector samsDatNoNoise = acd.samples;

  cout << myname << line << endl;
  cout << myname << "Create data by FFT convolution with response." << endl;
  AdcChannelData acdchk;
  acdchk.run = 123;
  acdchk.event = 456;
  acdchk.channel = 789;
  acdchk.samples = samsTru;
  DataMap retchk = pconfft->update(acdchk);
  assert( retchk == 0 );
  cout << "      # samples: " << acdchk.samples.size() << endl;
  cout << "     # DFT amps: " << acdchk.dftmags.size() << endl;
  cout << "   # DFT phases: " << acdchk.dftphases.size() << endl;
  assert( acdchk.samples.size() == nsam );
  AdcSignalVector samsDatCheck = acdchk.samples;

  cout << myname << line << endl;
  cout << myname << "Compare convolutions." << endl;
  Index nbad = 0;
  ostringstream eout;
  eout << "           Direct       FFT" << endl;
  for ( Index isam=0; isam<nsam; ++isam ) {
    float qdir = samsDatNoNoise[isam];
    float qfft = samsDatCheck[isam];
    eout << setw(5) << isam <<  ": " << fixed << setw(10) << qdir << fixed << setw(10) << qfft;
    if ( fabs(qfft - qdir) > 1.e-5 ) {
      eout << " *";
      ++nbad;
    }
    eout << endl;
  }
  if ( nbad != 0 ) {
    cout << myname << "Convolutions disagree:" << endl;
    cout << eout.str() << endl;
  }
  
  assert( nbad == 0 );

  cout << myname << line << endl;
  cout << myname << "Add noise." << endl;
  if ( noiseLev > 0.0 ) {
    for ( float& val : acd.samples ) val += gRandom->Gaus(0.0, noiseLev);
  }
  AdcSignalVector samsDat = acd.samples;

  cout << myname << line << endl;
  cout << myname << "Deconvolute." << endl;
  ret = pdco->update(acd);
  AdcSignalVector samsDco = acd.samples;
  assert( ret == 0 );

  cout << myname << line << endl;
  cout << myname << "Check integrals." << endl;
  float ares = 0.0;
  for ( float sig : res ) ares += sig;
  float aTru = 0.0;
  for ( float sig : samsTru ) aTru += sig;
  float aDatNoNoise = 0.0;
  for ( float sig : samsDatNoNoise ) aDatNoNoise += sig;
  float aDat = 0.0;
  for ( float sig : samsDat ) aDat += sig;
  float aDco = 0.0;
  for ( float sig : samsDco ) aDco += sig;
  cout << myname << "  Response function: " << ares << endl;
  cout << myname << "       Input signal: " << aTru << endl;
  cout << myname << " Data without noise: " << aDatNoNoise << endl;
  cout << myname << "    Data with noise: " << aDat << endl;
  cout << myname << "      Deconvolution: " << aDco << endl;

  cout << myname << line << endl;
  cout << myname << "Plot spectra." << endl;
  Plot timPlot;
  timPlot.addHist(samsTru, "True", 3, 2, 1);
  //timPlot.addHist(res, "Response", 5, 2, 1);
  timPlot.addHist(samsDatNoNoise, "DataNoNoise", 38, 4, 1);
  timPlot.addHist(samsDat, "Data", 1, 2, 1);
  timPlot.addHist(samsDco, "Deconvoluted", 2, 2, 1);
  timPlot.print(mypre + "_time.png", -50, 200);

  cout << myname << line << endl;
  cout << myname << "Done." << endl;
  return 0;
}

//**********************************************************************

int main(int argc, char* argv[]) {
  bool useExistingFcl = false;
  float noiseLev = 4.0;
  float sigmaFilter = 3.0;
  Index nsam = 200;
  if ( argc > 1 ) {
    string sarg(argv[1]);
    if ( sarg == "-h" ) {
      cout << "Usage: " << argv[0] << " [KEEPFCL [NOISE [SIGMAFIL [NSAM [SETSEED]]]]]" << endl;
      cout << "  KEEPFCL = if 1 (or true), existing FCL file is used [false]" << endl;
      cout << "  NOISE = noise sigma [2.0]" << endl;
      cout << "  SIGMAFIL = Filter sigma [Tick] [2.0]" << endl;
      cout << "  NSAM = data length [200]" << endl;
      cout << "  SETSEED = if 1 (or true), a new sed is used for the noise" << endl;
      return 0;
    }
    useExistingFcl = sarg == "true" || sarg == "1";
  }
  if ( argc > 2 ) {
    string sarg(argv[2]);
    istringstream ssarg(sarg);
    ssarg >> noiseLev;
  }
  if ( argc > 3 ) {
    string sarg(argv[3]);
    istringstream ssarg(sarg);
    ssarg >> sigmaFilter;
  }
  if ( argc > 4 ) {
    string sarg(argv[4]);
    nsam = std::stoi(sarg);
  }
  if ( argc > 5 ) {
    string sarg(argv[5]);
    if ( sarg == "true" || sarg == "1" ) gRandom->SetSeed();
  }
  return test_AdcDeconvoluteFFT(useExistingFcl, noiseLev, sigmaFilter, nsam);
}

//**********************************************************************
