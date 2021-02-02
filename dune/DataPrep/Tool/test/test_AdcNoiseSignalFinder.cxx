// test_AdcNoiseSignalFinder.cxx
//
// David Adams
// April 2017
//
// Test AdcNoiseSignalFinder.

#include <string>
#include <iostream>
#include <fstream>
#include "dune/DuneInterface/Tool/TpcDataTool.h"
#include "dune/ArtSupport/DuneToolManager.h"

#undef NDEBUG
#include <cassert>

using std::string;
using std::cout;
using std::endl;
using std::ofstream;
using fhicl::ParameterSet;

//**********************************************************************

int test_AdcNoiseSignalFinder(bool useExistingFcl =false) {
  const string myname = "test_AdcNoiseSignalFinder: ";
#ifdef NDEBUG
  cout << myname << "NDEBUG must be off." << endl;
  abort();
#endif
  string line = "-----------------------------";

  cout << myname << line << endl;
  string fclfile = "test_AdcNoiseSignalFinder.fcl";
  if ( ! useExistingFcl ) {
    cout << myname << "Creating top-level FCL." << endl;
    ofstream fout(fclfile.c_str());
    fout << "tools: {" << endl;
    fout << "  mytool: {" << endl;
    fout << "    tool_type: AdcNoiseSignalFinder" << endl;
    fout << "    LogLevel: 1" << endl;
    fout << "    SigFracMax: 0.2" << endl;
    fout << "    ThresholdMin: 20" << endl;
    fout << "    ThresholdRatio: 5" << endl;
    fout << "    ThresholdRatioTol: 0.1" << endl;
    fout << "    MaxLoop: 10" << endl;
    fout << "    BinsAfter: 10" << endl;
    fout << "    BinsBefore: 5" << endl;
    fout << "    FlagPositive: true" << endl;
    fout << "    FlagNegative: true" << endl;
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
  assert( tm.toolNames().size() == 1 );

  cout << myname << line << endl;
  cout << myname << "Fetching tool." << endl;
  auto psgf = tm.getPrivate<TpcDataTool>("mytool");
  assert( psgf != nullptr );
  auto psgfmod = tm.getPrivate<TpcDataTool>("mytool");
  assert( psgfmod != nullptr );

  cout << myname << line << endl;
  cout << myname << "Create data and call tool." << endl;
  AdcChannelData data;
  int wnoi = 20;
  float signoi = wnoi/sqrt(12.0);
  srand(12345);
  for ( AdcIndex itic=0; itic<100; ++itic ) {
    float xadc = rand()%wnoi - 0.5*wnoi;
    data.samples.push_back(xadc);
  }
  data.samples[30] = 150.0;
  assert( data.signal.size() == 0 );
  assert( data.rois.size() == 0 );
  assert( data.samples[30] = 150 );

  cout << myname << line << endl;
  cout << myname << "Running tool." << endl;
  DataMap resmod = psgfmod->update(data);
  resmod.print();

  cout << myname << line << endl;
  cout << myname << "Checking results." << endl;
  assert( resmod == 0 );
  assert( resmod.getInt("nsfLoopCount") > 1 );
  assert( resmod.getInt("nsfRoiCount") == 1 );
  assert( data.signal.size() == 100 );
  assert( data.rois.size() == 1 );
  assert( data.rois[0].first == 25 );
  assert( data.rois[0].second == 40 );
  float mdNoise = data.getMetadata("nsfNoise");
  assert( mdNoise > 0.9*signoi );
  assert( mdNoise < 1.1*signoi );
  float mdSigFrac = data.getMetadata("nsfSigFrac");
  int mdNsig = int(100*mdSigFrac + 0.1);
  assert( mdNsig == 16 );
  float mdThresh = data.getMetadata("nsfThreshold");
  float expThresh = 5.0*signoi;
  assert( mdThresh > 0.9*expThresh );
  assert( mdThresh < 1.1*expThresh );

  cout << myname << line << endl;
  cout << myname << "Done." << endl;
  return 0;
}

//**********************************************************************

int main(int argc, char* argv[]) {
  bool useExistingFcl = false;
  if ( argc > 1 ) {
    string sarg(argv[1]);
    if ( sarg == "-h" ) {
      cout << "Usage: " << argv[0] << " [ARG]" << endl;
      cout << "  If ARG = true, existing FCL file is used." << endl;
      return 0;
    }
    useExistingFcl = sarg == "true" || sarg == "1";
  }
  return test_AdcNoiseSignalFinder(useExistingFcl);
}

//**********************************************************************
