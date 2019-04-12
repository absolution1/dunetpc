// test_ExpTailRemover.cxx
//
// David Adams
// April 2019
//
// Test ExpTailRemover.

#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
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
using std::istringstream;
using std::setw;
using fhicl::ParameterSet;

using Index = unsigned int;
using IndexVector = std::vector<Index>;
using FloatVector = AdcSignalVector;

//**********************************************************************

namespace {

struct SigStats {
  Index count = 0;
  float mean = 0.0;
  float rms = 0.0;
  SigStats(const AdcChannelData& acd, string lab);
};

SigStats::SigStats(const AdcChannelData& acd, string lab) {
  Index nSig = 0;
  double sumSig = 0.0;
  double sumSig2 = 0.0;
  Index nsam = acd.samples.size();
  if ( acd.signal.size() < nsam ) return;
  bool dbg = false;
  for ( Index isam=0; isam<nsam; ++isam ) {
    if ( acd.signal[isam] ) continue;
    ++nSig;
    double sig = acd.samples[isam];
    sumSig += sig;
    sumSig2 += sig*sig;
    if ( dbg ) cout << lab << setw(5) << isam << ":" << setw(10) << sig << endl;
  }
  if ( nSig == 0 ) return;
  count = nSig;
  mean = sumSig/nSig;
  rms = sqrt(sumSig2/nSig - mean*mean);
  if ( lab.size() ) {
    cout << lab << "  # signal: " << count << endl;
    cout << lab << "      mean: " << mean << endl;
    cout << lab << "       RMS: " << rms << endl;
  }
}

}  // end unnamed namespace

//**********************************************************************

int test_ExpTailRemover(bool useExistingFcl, Index flag, float noiseSigma, bool setSeed) {
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
    fout << "              SignalTool: \"sigfind\"" << endl;
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

  cout << myname << line << endl;
  cout << myname << "Add noise to the signal." << endl;
  if ( setSeed ) gRandom->SetSeed();
  for ( float& sig : sigs1 ) sig += gRandom->Gaus(0.0, noiseSigma);

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
  cout << myname << "Add pedestal and tail to the signal to create data." << endl;
  assert( sta.setSignal(sigs1) == 0 );

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
  cout << myname << "Check input noise." << endl;
  SigStats inStats(acd, myname);
  assert( inStats.count > 0 );

  cout << myname << line << endl;
  cout << myname << "Use tool to remove tail from data." << endl;
  DataMap res = ptoo->update(acd);
  res.print();
  assert ( res == 0 );

  cout << myname << line << endl;
  cout << myname << "Check output noise." << endl;
  SigStats ouStats(acd, myname);
  assert( ouStats.count > 0 );

  cout << myname << "Done." << endl;
  return 0;
}

//**********************************************************************

int main(int argc, char* argv[]) {
  bool useExistingFcl = false;
  Index flag = 2;
  float noiseSigma = 2.0;
  bool setSeed = false;
  if ( argc > 1 ) {
    string sarg(argv[1]);
    if ( sarg == "-h" ) {
      cout << "Usage: " << argv[0] << " [ARG [OPT [noise [setSeed]]]]" << endl;
      cout << "  If ARG = true, existing FCL file is used." << endl;
      cout << "  OPT [2] is SignalOption = 0, 1, 2, or 3" << endl;
      cout << "  noise [2.0]  is the sigma of the noise added to the data" << endl;
      cout << "  setSeed nonzero means a random random seed for the noise" << endl;
      return 0;
    }
    useExistingFcl = sarg == "true" || sarg == "1";
  }
  if ( argc > 2 ) {
    string sarg(argv[2]);
    flag = std::stoi(sarg);
  }
  if ( argc > 3 ) {
    string sarg(argv[3]);
    istringstream ssarg(sarg);
    ssarg >> noiseSigma;
  }
  if ( argc > 4 ) {
    string sarg(argv[3]);
    setSeed = sarg != "0";
  }
  return test_ExpTailRemover(useExistingFcl, flag, noiseSigma, setSeed);
}

//**********************************************************************
