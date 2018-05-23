// test_DuneAdcSignalFinder.cxx
//
// David Adams
// April 2017
//
// Test DuneAdcSignalFinder.

#include <string>
#include <iostream>
#include <fstream>
#include "dune/DuneInterface/Tool/AdcChannelTool.h"
#include "dune/ArtSupport/DuneToolManager.h"

#undef NDEBUG
#include <cassert>

using std::string;
using std::cout;
using std::endl;
using std::ofstream;
using fhicl::ParameterSet;

//**********************************************************************

int test_DuneAdcSignalFinder(bool useExistingFcl =false) {
  const string myname = "test_DuneAdcSignalFinder: ";
#ifdef NDEBUG
  cout << myname << "NDEBUG must be off." << endl;
  abort();
#endif
  string line = "-----------------------------";

  cout << myname << line << endl;
  string fclfile = "test_DuneAdcSignalFinder.fcl";
  if ( ! useExistingFcl ) {
    cout << myname << "Creating top-level FCL." << endl;
    ofstream fout(fclfile.c_str());
    fout << "tools: {" << endl;
    fout << "  mytool: {" << endl;
    fout << "    tool_type: DuneAdcSignalFinder" << endl;
    fout << "    LogLevel: 1" << endl;
    fout << "    NoiseSigma: 0.0" << endl;
    fout << "    NSigmaStart: 3.0" << endl;
    fout << "    NSigmaEnd: 1.0" << endl;
    fout << "    TicksAfter: 10" << endl;
    fout << "    TicksBefore: 5" << endl;
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
  auto psgf = tm.getPrivate<AdcChannelTool>("mytool");
  assert( psgf != nullptr );
  auto psgfmod = tm.getPrivate<AdcChannelTool>("mytool");
  assert( psgfmod != nullptr );

  cout << myname << line << endl;
  cout << myname << "Create data and call tool." << endl;
  AdcChannelData data;
  for ( AdcIndex itic=0; itic<100; ++itic ) {
    float xadc = rand()%20 - 10.0;
    data.samples.push_back(xadc);
  }
  data.sampleNoise = 40.0;
  data.samples[30] = 150.0;
  assert( data.signal.size() == 0 );
  assert( data.rois.size() == 0 );
  assert( data.samples[30] = 150 );
  assert( data.samples[31] = 130 );
  assert( data.samples[32] =  90 );
  assert( data.samples[33] =  70 );
  assert( data.samples[34] =  45 );
  assert( data.samples[35] =  30 );
  assert( data.samples[36] =  20 );
  assert( data.samples[37] =  15 );

  cout << myname << line << endl;
  cout << myname << "Running tool." << endl;
  DataMap resmod = psgfmod->update(data);
  resmod.print();

  cout << myname << line << endl;
  cout << myname << "Checking results." << endl;
  assert( resmod == 0 );
  assert( resmod.getInt("nroi") == 1 );
  assert( data.signal.size() == 100 );
  assert( data.rois.size() == 1 );
  assert( data.rois[0].first == 25 );
  assert( data.rois[0].second == 44 );

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
  return test_DuneAdcSignalFinder(useExistingFcl);
}

//**********************************************************************
