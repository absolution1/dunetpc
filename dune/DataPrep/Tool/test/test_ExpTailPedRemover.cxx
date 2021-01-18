// test_ExpTailPedRemover.cxx
//
// David Adams
// April 2019
//
// Test ExpTailPedRemover.

#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include "dune/DuneInterface/Tool/AdcChannelTool.h"
#include "dune/DuneCommon/SampleTailer.h"
#include "dune/DuneCommon/LineColors.h"
#include "dune/DuneCommon/TPadManipulator.h"
#include "dune/ArtSupport/DuneToolManager.h"
#include "TRandom.h"
#include "TH1F.h"
#include "TGraphErrors.h"

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
  float rms = 0.0;
  Index bgcount = 0;
  float bgmean = 0.0;
  float bgrms = 0.0;
  SigStats(const AdcChannelData& acd, const FloatVector& sigTrue, string lab);
};

SigStats::SigStats(const AdcChannelData& acd, const FloatVector& sigTrue, string lab) {
  Index nSig = 0;
  double sumSig2 = 0.0;
  Index nNsg = 0;
  double sumNsg = 0.0;
  double sumNsg2 = 0.0;
  Index nsam = acd.samples.size();
  if ( acd.signal.size() < nsam ) return;
  for ( Index isam=0; isam<nsam; ++isam ) {
    if ( acd.signal[isam] ) {
      ++nSig;
      double sig = acd.samples[isam] - sigTrue[isam];
      sumSig2 += sig*sig;
    } else {
      ++nNsg;
      double sig = acd.samples[isam];
      sumNsg += sig;
      sumNsg2 += sig*sig;
    }
  }
  count = nSig;
  rms = sqrt(sumSig2/nSig);
  bgcount = nNsg;
  bgmean = sumNsg/nNsg;
  bgrms = sqrt(sumNsg2/nNsg - bgmean*bgmean);
  
  if ( lab.size() ) {
    cout << lab << "   # signal: " << count << endl;
    cout << lab << "       # BG: " << bgcount << endl;
    cout << lab << "    BG mean: " << bgmean << endl;
    cout << lab << "     BG RMS: " << bgrms << endl;
    cout << lab << " Signal RMS: " << rms << endl;
  }
}

void drawResults(const FloatVector sigin, const FloatVector& sigraw, const FloatVector& sigout) {
  const string myname = "drawResults: ";
  string line = "-----------------------------";
  cout << myname << line << endl;
  Index nsam = 300;
  cout << myname << "Draw data." << endl;
  LineColors lc;
  TH1* phi = new TH1F("hin", "ExpTailPedRemover fit;Tick;signal", nsam, 0, nsam);
  phi->SetStats(0);
  phi->SetLineWidth(2);
  phi->SetLineColor(lc.blue());
  TGraphErrors *pgr = new TGraphErrors(nsam);
  pgr->SetMarkerStyle(2);
  pgr->SetMarkerColor(lc.black());
  TGraphErrors *pgo = new TGraphErrors(nsam);
  pgo->SetMarkerStyle(24);
  pgo->SetMarkerColor(lc.red());
  for ( Index isam=0; isam<nsam; ++ isam ) {
    phi->SetBinContent(isam+1, sigin[isam]);
    pgr->SetPoint(isam, isam, sigraw[isam]);
    pgo->SetPoint(isam, isam, sigout[isam]);
  }
  TPadManipulator man(1000, 500);
  man.add(phi, "HIST");
  man.add(pgr, "P");
  man.add(pgo, "P");
  string fnam = "test_ExpTailPedRemover.png";
  man.setRangeY(-40, 240);
  man.addHorizontalLine(0.0);
  man.addAxis();
  man.print(fnam);
  cout << myname << "Plot saved as " << fnam << endl;
}

}  // end unnamed namespace

//**********************************************************************

int test_ExpTailPedRemover(bool useExistingFcl, Index flag, float ped, float slope, float noiseSigma, bool setSeed) {
  const string myname = "test_ExpTailPedRemover: ";
#ifdef NDEBUG
  cout << myname << "NDEBUG must be off." << endl;
  abort();
#endif
  string line = "-----------------------------";

  cout << myname << line << endl;
  string fclfile = "test_ExpTailPedRemover.fcl";
  float decayTime = 100.0;
  float tail0 = -15.0;
  int tick0 = 150;
  int pedDegree = slope == 0.0 ? 0 : 1;
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
    fout << "               tool_type: ExpTailPedRemover" << endl;
    fout << "                LogLevel: 3" << endl;
    fout << "              SignalFlag: " << flag << endl;
    fout << "    SignalIterationLimit: 10" << endl;
    fout << "              SignalTool: \"sigfind\"" << endl;
    fout << "               DecayTime: " << decayTime << endl;
    fout << "                 MaxTick: 400" << endl;
    fout << "               PedDegree: " << pedDegree << endl;
    fout << "                PedTick0: " << tick0 << endl;
    fout << "                PedFreqs: []" << endl;
    fout << "    NoWarnStatuses: []" << endl;
    fout << "    IncludeChannelRanges: [\"all\"]" << endl;
    fout << "    ExcludeChannelRanges: []" << endl;
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
  IndexVector peakPoss = {70, 100, 115, 230};
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
  FloatVector sigs0 = sigs1;

  cout << myname << line << endl;
  cout << myname << "Add noise to the signal." << endl;
  if ( setSeed ) gRandom->SetSeed();
  for ( float& sig : sigs1 ) sig += gRandom->Gaus(0.0, noiseSigma);

  cout << myname << "Create sample tailer." << endl;
  SampleTailer sta(decayTime);
  FloatVector peds;
  float curv = 0.0;
  if ( slope == 0.0 ) {
    sta.setPedestal(ped);
  } else {
    peds.resize(nsam);
    for ( Index isam=0; isam<nsam; ++isam ) {
      float x = float(isam) - tick0;
      peds[isam] = ped + slope*x + curv*x*x;
    }
    sta.setPedestalVector(&peds);
  }
  sta.setTail0(tail0);
  sta.setUnit("ADC count");
  cout << myname << "  decayTime: " << sta.decayTime() << endl;
  cout << myname << "       beta: " << sta.beta() << endl;
  cout << myname << "      alpha: " << sta.alpha() << endl;
  cout << myname << "   pedestal: " << ped << endl;
  cout << myname << "      slope: " << slope << endl;
  cout << myname << "  curvature: " << curv << endl;
  cout << myname << "      tail0: " << sta.tail0() << endl;

  cout << myname << line << endl;
  cout << myname << "Add pedestal and tail to the signal to create data." << endl;
  assert( sta.setSignal(sigs1) == 0 );

  cout << myname << line << endl;
  cout << myname << "Create channel data." << endl;
  AdcChannelData acd;
  acd.setEventInfo(123, 456);
  acd.channel = 789;
  acd.pedestal = 1000.0;
  acd.samples = sta.data();
  acd.signal = isSignal;

  cout << myname << line << endl;
  cout << myname << "Check input noise." << endl;
  SigStats inStats(acd, sigs0, myname);
  assert( inStats.count > 0 );

  cout << myname << line << endl;
  cout << myname << "Use tool to remove tail from data." << endl;
  DataMap res = ptoo->update(acd);
  res.print();
  assert ( res == 0 );

  cout << myname << line << endl;
  cout << myname << "Check output noise." << endl;
  SigStats ouStats(acd, sigs0, myname);
  assert( ouStats.count > 0 );

  drawResults(sigs0, sta.data(), acd.samples);

  cout << myname << "Done." << endl;
  return 0;
}

//**********************************************************************

int main(int argc, char* argv[]) {
  bool useExistingFcl = false;
  Index flag = 2;
  float noiseSigma = 2.0;
  float ped = 5.0;
  float slope = 0.02;
  bool setSeed = false;
  if ( argc > 1 ) {
    string sarg(argv[1]);
    if ( sarg == "-h" ) {
      cout << "Usage: " << argv[0] << " [ARG [SLOPE [OPT [noise [setSeed]]]]]" << endl;
      cout << "  If ARG = true, existing FCL file is used." << endl;
      cout << "  SLOPE [0.02] is the pedestal slope" << endl;
      cout << "  OPT [2] is SignalOption = 0, 1, 2, or 3" << endl;
      cout << "  noise [2.0]  is the sigma of the noise added to the data" << endl;
      cout << "  setSeed nonzero means a random random seed for the noise" << endl;
      return 0;
    }
    useExistingFcl = sarg == "true" || sarg == "1";
  }
  if ( argc > 2 ) {
    string sarg(argv[2]);
    istringstream ssarg(sarg);
    ssarg >> slope;
  }
  if ( argc > 3 ) {
    string sarg(argv[3]);
    flag = std::stoi(sarg);
  }
  if ( argc > 4 ) {
    string sarg(argv[4]);
    istringstream ssarg(sarg);
    ssarg >> noiseSigma;
  }
  if ( argc > 5 ) {
    string sarg(argv[5]);
    setSeed = sarg != "0";
  }
  return test_ExpTailPedRemover(useExistingFcl, flag, ped, slope, noiseSigma, setSeed);
}

//**********************************************************************
