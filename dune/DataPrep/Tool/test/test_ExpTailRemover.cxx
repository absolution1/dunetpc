// test_ExpTailRemover.cxx
//
// David Adams
// April 2019
//
// Test ExpTailRemover.

#include <string>
#include <iostream>
#include <fstream>
#include <iomanip>
#include "dune/DuneInterface/Tool/AdcChannelTool.h"
#include "dune/DuneCommon/SampleTailer.h"
#include "dune/ArtSupport/DuneToolManager.h"
#include "TRandom.h"

#undef NDEBUG
#include <cassert>

using std::string;
using std::cout;
using std::endl;
using std::ofstream;
using std::setw;
using fhicl::ParameterSet;

using Index = unsigned int;
using IndexVector = std::vector<Index>;
using FloatVector = AdcSignalVector;

//**********************************************************************

int test_ExpTailRemover(bool useExistingFcl, Index flag) {
  const string myname = "test_ExpTailRemover: ";
#ifdef NDEBUG
  cout << myname << "NDEBUG must be off." << endl;
  abort();
#endif
  string line = "-----------------------------";

  cout << myname << line << endl;
  string fclfile = "test_ExpTailRemover.fcl";
  float decayTime = 100.0;
  float ped = 5.0;
  float tail0 = -15.0;
  if ( ! useExistingFcl ) {
    cout << myname << "Creating top-level FCL." << endl;
    ofstream fout(fclfile.c_str());
    fout << "tools: {" << endl;
    fout << "  sigfind: {" << endl;
    fout << "    tool_type: AdcThresholdSignalFinder" << endl;
    fout << "    LogLevel: 1" << endl;
    fout << "    Threshold: 10" << endl;
    fout << "    BinsAfter: 10" << endl;
    fout << "    BinsBefore: 5" << endl;
    fout << "    FlagPositive: true" << endl;
    fout << "    FlagNegative: true" << endl;
    fout << "  }" << endl;
    fout << "  mytool: {" << endl;
    fout << "               tool_type: ExpTailRemover" << endl;
    fout << "                LogLevel: 3" << endl;
    fout << "              SignalFlag: " << flag << endl;
    fout << "    SignalIterationLimit: 10" << endl;
    fout << "              SignalTool: \"\"" << endl;
    fout << "               DecayTime: " << decayTime << endl;
    fout << "             CorrectFlag: []" << endl;
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
  auto ptoo = tm.getPrivate<AdcChannelTool>("mytool");
  assert( ptoo != nullptr );

  cout << myname << "Create signals." << endl;
  Index nsam = 300;
  FloatVector pulse = {  0.1,  4.5, 15.2, 66.4, 94.3, 100.0, 96.5, 88.4, 72.6, 58.4,
                        42.3, 35.1, 26.0, 18.6, 12.6,   8.8,  6.9,  4.4,  2.0, 0.3 };
  Index npul = pulse.size();
  FloatVector sigs1(nsam, 0.0);
  AdcFilterVector isSignal(nsam, false);
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
      isSignal[isam] = true;
    }
  }

  cout << myname << "Create sample tailer." << endl;
  SampleTailer sta(decayTime);
  sta.setPedestal(ped);
  sta.setTail0(tail0);
  sta.setUnit("ADC count");
  cout << myname << "  decayTime: " << sta.decayTime() << endl;
  cout << myname << "       beta: " << sta.beta() << endl;
  cout << myname << "      alpha: " << sta.alpha() << endl;
  cout << myname << "   pedestal: " << sta.pedestal() << endl;
  cout << myname << "      tail0: " << sta.tail0() << endl;

  cout << myname << line << endl;
  cout << myname << "Create data from signal." << endl;
  assert( sta.setSignal(sigs1) == 0 );

  cout << myname << line << endl;
  cout << myname << "Create data." << endl;
  AdcSignalVector dats1 = sta.data();
  float noiseLevel = 2.0;
  for ( float& dat : dats1 ) dat += gRandom->Gaus(noiseLevel, 0.0);

  cout << myname << line << endl;
  cout << myname << "Create channel data." << endl;
  AdcChannelData acd;
  acd.run = 123;
  acd.event = 456;
  acd.channel = 789;
  acd.pedestal = 1000.0;
  acd.samples = sta.data();
  acd.signal = isSignal;

  cout << myname << line << endl;
  cout << myname << "Use tool to remove tail from data." << endl;
  DataMap res = ptoo->update(acd);
  res.print();
  assert ( res == 0 );

  cout << myname << "Done." << endl;
  return 0;
}

//**********************************************************************

int main(int argc, char* argv[]) {
  bool useExistingFcl = false;
  Index flag = 1;
  if ( argc > 1 ) {
    string sarg(argv[1]);
    if ( sarg == "-h" ) {
      cout << "Usage: " << argv[0] << " [ARG]" << endl;
      cout << "  If ARG = true, existing FCL file is used." << endl;
      return 0;
    }
    useExistingFcl = sarg == "true" || sarg == "1";
  }
  if ( argc > 2 ) {
    string sarg(argv[2]);
    flag = std::stoi(sarg);
  }
  return test_ExpTailRemover(useExistingFcl, flag);
}

//**********************************************************************
