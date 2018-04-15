// test_AdcSampleFiller.cxx
//
// David Adams
// April 2017
//
// Test AdcSampleFiller.

#include <string>
#include <iostream>
#include <fstream>
#include "dune/DuneInterface/Tool/AdcChannelTool.h"
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

int test_AdcSampleFiller(bool useExistingFcl =false) {
  const string myname = "test_AdcSampleFiller: ";
#ifdef NDEBUG
  cout << myname << "NDEBUG must be off." << endl;
  abort();
#endif
  string line = "-----------------------------";

  cout << myname << line << endl;
  string fclfile = "test_AdcSampleFiller.fcl";
  if ( ! useExistingFcl ) {
    cout << myname << "Creating top-level FCL." << endl;
    ofstream fout(fclfile.c_str());
    fout << "tools: {" << endl;
    fout << "  mytool: {" << endl;
    fout << "    tool_type: AdcSampleFiller" << endl;
    fout << "    LogLevel: 1" << endl;
    fout << "    AdcUnderflow: 0" << endl;
    fout << "    AdcOverflow: 4095" << endl;
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
  auto pasf = tm.getPrivate<AdcChannelTool>("mytool");
  assert( pasf != nullptr );
  auto pasfmod = tm.getPrivate<AdcChannelTool>("mytool");
  assert( pasfmod != nullptr );

  cout << myname << line << endl;
  cout << myname << "Create data and call tool." << endl;
  AdcChannelData data;
  for ( AdcIndex itic=0; itic<60; ++itic ) {
    AdcIndex adc = 100*itic;
    data.raw.push_back(adc);
  }
  assert( data.raw.size() == 60 );
  assert( data.samples.size() == 0 );
  assert( data.signal.size() == 0 );
  assert( data.rois.size() == 0 );

  cout << myname << line << endl;
  cout << myname << "Running tool without pedestal." << endl;
  DataMap resmodBad = pasfmod->update(data);
  resmodBad.print();
  assert( resmodBad != 0 );

  cout << myname << line << endl;
  cout << myname << "Running tool." << endl;
  data.pedestal = 1000.0;
  DataMap resmod = pasfmod->update(data);
  resmod.print();

  cout << myname << line << endl;
  cout << myname << "Checking results." << endl;
  assert( resmod == 0 );
  assert( !resmod );
  assert( resmod.getInt("nUnderflow") == 1 );
  assert( resmod.getInt("nOverflow") == 19 );
  assert( resmod.getInt("nOutOfRange") == 19 );
  assert( data.raw.size() == 60 );
  assert( data.samples.size() == 60 );
  assert( data.signal.size() == 0 );

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
  return test_AdcSampleFiller(useExistingFcl);
}

//**********************************************************************
