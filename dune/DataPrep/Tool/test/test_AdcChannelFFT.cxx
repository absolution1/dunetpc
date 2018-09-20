// test_AdcChannelFFT.cxx
//
// David Adams
// April 2017
//
// Test AdcChannelFFT.

#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <iomanip>
#include "dune/DuneInterface/Tool/AdcChannelTool.h"
#include "dune/ArtSupport/DuneToolManager.h"
#include <TRandom.h>

#undef NDEBUG
#include <cassert>

using std::string;
using std::cout;
using std::endl;
using std::ostringstream;
using std::ofstream;
using fhicl::ParameterSet;
using std::vector;
using std::setw;
using std::fixed;

using Index = unsigned int;

//**********************************************************************

int test_AdcChannelFFT(bool useExistingFcl, Index len) {
  const string myname = "test_AdcChannelFFT: ";
#ifdef NDEBUG
  cout << myname << "NDEBUG must be off." << endl;
  abort();
#endif
  string line = "-----------------------------";

  cout << myname << line << endl;
  string fclfile = "test_AdcChannelFFT.fcl";
  if ( ! useExistingFcl ) {
    cout << myname << "Creating top-level FCL." << endl;
    ofstream fout(fclfile.c_str());
    fout << "tools: {" << endl;
    fout << "  mytool: {" << endl;
    fout << "           tool_type: AdcChannelFFT" << endl;
    fout << "            LogLevel: 2" << endl;
    fout << "           FirstTick: 0" << endl;
    fout << "               NTick: 0" << endl;
    fout << "             NormOpt: 1" << endl;
    fout << "              Action: 3" << endl;
    fout << "           ReturnOpt: 13" << endl;
    fout << "  }" << endl;
    fout << "  mytoolinv: {" << endl;
    fout << "           tool_type: AdcChannelFFT" << endl;
    fout << "            LogLevel: 3" << endl;
    fout << "           FirstTick: 0" << endl;
    fout << "               NTick: 0" << endl;
    fout << "             NormOpt: 1" << endl;
    fout << "              Action: 13" << endl;
    fout << "           ReturnOpt: 13" << endl;
    fout << "  }" << endl;
    fout << "}" << endl;
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
  assert( tm.toolNames().size() == 2 );

  cout << myname << line << endl;
  cout << myname << "Fetching tool." << endl;
  auto pfft = tm.getPrivate<AdcChannelTool>("mytool");
  assert( pfft != nullptr );

  cout << myname << line << endl;
  cout << myname << "Create data." << endl;
  AdcChannelData acd;
  acd.channel = 123;
  vector<float> sams = {   3.0,  6.0,  16.1,  28.6,  30.2,  27.7,  16.3,   9.6,  4.2, -1.0,
                          -2.3,  -4.2,  -9.2, -18.6, -21.9, -29.0, -24.3, -14.2, -5.0, -3.0};
  if ( len > 0 ) sams.resize(len, 0.0);
  Index nsam = sams.size();
  assert( sams.size() == nsam );
  cout << myname << "Sample length: " << nsam << endl;
  float samsum = 0.0;
  for ( float sam : sams ) samsum += sam;
  cout << myname << "Sample mean: " << samsum/nsam << endl;
  acd.samples = sams;
  Index nmag = (nsam+2)/2;
  Index npha = (nsam+1)/2;
  cout << myname << "        # samples: " << nsam << endl;
  cout << myname << "Expected mag size: " << nmag << endl;
  cout << myname << "Expected mag size: " << nmag << endl;
  cout << myname << "Expected pha size: " << npha << endl;

  cout << myname << line << endl;
  cout << myname << "Call tool view." << endl;
  DataMap dm = pfft->view(acd);
  dm.print();
  assert( dm == 0 );
  assert( Index(dm.getInt("fftNTick")) == nsam );
  assert( Index(dm.getInt("fftTick0")) == 0 );
  assert( Index(dm.getInt("fftNMag")) == nmag );
  assert( Index(dm.getInt("fftNPhase")) == npha );
  assert( dm.haveFloatVector("fftMags") );
  assert( dm.getFloatVector("fftMags").size() == nmag );
  assert( dm.haveFloatVector("fftPhases") );
  assert( dm.getFloatVector("fftPhases").size() == npha );
  assert( dm.haveFloatVector("fftReals") );
  assert( dm.getFloatVector("fftReals").size() == nsam );
  assert( dm.haveFloatVector("fftImags") );
  assert( dm.getFloatVector("fftImags").size() == nsam );

  cout << myname << line << endl;
  cout << myname << "Check power." << endl;
  float pwr1 = 0.0;
  for ( float sam : sams ) pwr1 += sam*sam;
  float pwr2 = 0.0;
  for ( float mag : dm.getFloatVector("fftMags") ) pwr2 += mag*mag;
  float pwr3 = 0.0;
  for ( float x : dm.getFloatVector("fftReals") ) pwr3 += x*x;
  for ( float x : dm.getFloatVector("fftImags") ) pwr3 += x*x;
  cout << myname << "Tick power: " << pwr1 << endl;
  cout << myname << "Freq power: " << pwr2 << endl;
  cout << myname << "Frq2 power: " << pwr3 << endl;
  assert( fabs(pwr2 - pwr1) < 1.e-5*pwr1 );
  assert( fabs(pwr3 - pwr1) < 1.e-5*pwr1 );

  cout << myname << line << endl;
  cout << myname << "Call tool update." << endl;
  dm = pfft->update(acd);
  dm.print();
  assert( acd.dftmags.size() == nmag );
  assert( acd.dftphases.size() == npha );

  cout << myname << line << endl;
  cout << myname << "Fetching inverse tool." << endl;
  auto pffi = tm.getPrivate<AdcChannelTool>("mytoolinv");
  assert( pffi != nullptr );

  cout << myname << line << endl;
  cout << myname << "Call inverse tool." << endl;
  acd.samples.clear();
  assert( acd.samples.size() == 0 );
  assert( acd.dftmags.size() == nmag );
  assert( acd.dftphases.size() == npha );
  dm = pffi->update(acd);
  dm.print();
  assert( dm == 0 );
  assert( acd.samples.size() == nsam );
  for ( Index isam=0; isam<nsam; ++isam ) {
    cout << setw(4) << isam << ":" << setw(10) << fixed << acd.samples[isam]
         << " ?= " << setw(10) << fixed << sams[isam] << endl;
    assert( fabs(acd.samples[isam] - sams[isam]) < 1.e-4 );
  }

  cout << myname << line << endl;
  cout << myname << "Done." << endl;
  return 0;
}

//**********************************************************************

int main(int argc, char* argv[]) {
  bool useExistingFcl = false;
  Index len = 0;
  if ( argc > 1 ) {
    string sarg(argv[1]);
    if ( sarg == "-h" ) {
      cout << "Usage: " << argv[0] << " [ARG] [LEN]" << endl;
      cout << "  If ARG = true, existing FCL file is used." << endl;
      return 0;
    }
    useExistingFcl = sarg == "true" || sarg == "1";
  }
  if ( argc > 2 ) {
    string sarg(argv[2]);
    len = std::stoi(sarg);
  }
  return test_AdcChannelFFT(useExistingFcl, len);
}

//**********************************************************************
