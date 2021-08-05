// test_AdcThresholdSignalFinder.cxx
//
// David Adams
// April 2017
//
// Test AdcThresholdSignalFinder.

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

int test_AdcThresholdSignalFinder(bool useExistingFcl =false) {
  const string myname = "test_AdcThresholdSignalFinder: ";
#ifdef NDEBUG
  cout << myname << "NDEBUG must be off." << endl;
  abort();
#endif
  string line = "-----------------------------";

  cout << myname << line << endl;
  string fclfile = "test_AdcThresholdSignalFinder.fcl";
  if ( ! useExistingFcl ) {
    cout << myname << "Creating top-level FCL." << endl;
    ofstream fout(fclfile.c_str());
    fout << "tools: {" << endl;
    fout << "  mytool: {" << endl;
    fout << "    tool_type: AdcThresholdSignalFinder" << endl;
    fout << "    LogLevel: 1" << endl;
    fout << "    Threshold: \"100*[gain]/14.0\"" << endl;
    fout << "    BinsAfter: 10" << endl;
    fout << "    BinsBefore: 5" << endl;
    fout << "    FlagPositive: true" << endl;
    fout << "    FlagNegative: true" << endl;
    fout << "  }" << endl;
    fout << "  runDataTool: {" << endl;
    fout << "    tool_type: FclRunDataTool" << endl;
    fout << "    LogLevel: 1" << endl;
    fout << "    FileNames: [\"rundata.fcl\"]" << endl;
    fout << "  }" << endl;
    fout << "}" << endl;
    fout.close();
    ofstream fout2("rundata.fcl");
    fout2 << "run: 123" << endl;
    fout2 << "gain: 14.0" << endl;
    fout2 << "shaping: 2.0" << endl;
    fout2.close();
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
  auto psgf = tm.getPrivate<TpcDataTool>("mytool");
  assert( psgf != nullptr );
  auto psgfmod = tm.getPrivate<TpcDataTool>("mytool");
  assert( psgfmod != nullptr );

  cout << myname << line << endl;
  cout << myname << "Create data and call tool." << endl;
  AdcChannelData data;
  for ( AdcIndex itic=0; itic<100; ++itic ) {
    float xadc = rand()%20 - 10.0;
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
  assert( fabs(resmod.getFloat("threshold") - 100) < 0.1 );
  assert( resmod.getInt("nThresholdBins") == 1 );
  assert( resmod.getInt("nroi") == 1 );
  assert( data.signal.size() == 100 );
  assert( data.rois.size() == 1 );
  assert( data.rois[0].first == 25 );
  assert( data.rois[0].second == 40 );

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
  return test_AdcThresholdSignalFinder(useExistingFcl);
}

//**********************************************************************
