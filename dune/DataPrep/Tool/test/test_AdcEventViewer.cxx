// test_AdcEventViewer.cxx
//
// David Adams
// October 2017
//
// Test AdcEventViewer.

#include <string>
#include <iostream>
#include <fstream>
#include "dune/DuneInterface/Tool/AdcChannelTool.h"
#include "dune/ArtSupport/DuneToolManager.h"
#include "TRandom.h"

#undef NDEBUG
#include <cassert>

using std::string;
using std::cout;
using std::endl;
using std::ofstream;
using std::vector;
using fhicl::ParameterSet;

using Index = unsigned int;

//**********************************************************************

int test_AdcEventViewer(bool useExistingFcl =false) {
  const string myname = "test_AdcEventViewer: ";
#ifdef NDEBUG
  cout << myname << "NDEBUG must be off." << endl;
  abort();
#endif
  string line = "-----------------------------";

  cout << myname << line << endl;
  string fclfile = "test_AdcEventViewer.fcl";
  if ( ! useExistingFcl ) {
    cout << myname << "Creating top-level FCL." << endl;
    ofstream fout(fclfile.c_str());
    fout << "#include \"dataprep_tools.fcl\"" << endl;
    fout << "tools.mytool: {" << endl;
    fout << "  tool_type: AdcEventViewer" << endl;
    fout << "  LogLevel: 1" << endl;
    fout << "  EventHists: []" << endl;
    fout << "  EventGraphs: []" << endl;
    fout << "  ChannelRanges: []" << endl;
    fout << "  ChannelRangeLabel: \"\"" << endl;
    fout << "  ClockUnit: Mtick" << endl;
    fout << "  ClockRate: 50000000" << endl;
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
  assert( tm.toolNames().size() >= 2 );

  cout << myname << line << endl;
  cout << myname << "Fetching tool." << endl;
  auto ptoo = tm.getPrivate<AdcChannelTool>("mytool");
  assert( ptoo != nullptr );

  cout << myname << line << endl;
  cout << myname << "Create test data." << endl;
  AdcChannelDataMap acds;
  vector<float> pulse = { 2.0, -3.0, 0.0, 5.0, 24.0, 56.0, 123.0, 71.0, 52.1, 26.3,
                         12.5,  8.1, 4.5, 2.0, -1.0,  3.2,   1.1, -2.2,  0.1, -0.1};
  Index npul = 100;
  for ( Index icha=250; icha<270; ++icha ) {
    AdcChannelData& acd = acds[icha];
    acd.run = 111;
    acd.event = 123;
    acd.channel = icha;
    acd.samples.resize(npul*pulse.size());
    acd.sampleUnit = "ADC count";
    double sigma = 5.0;
    Index itck = 0;
    for ( Index ipul=0; ipul<npul; ++ipul ) {
      Index itck0 = itck;
      for ( Index ismp=0; ismp<pulse.size(); ++ismp, ++itck ) {
        acd.samples[itck] = pulse[ismp] + gRandom->Gaus(0.0, sigma);
      }
      acd.rois.emplace_back( itck0+2, itck-2);
    }
    assert( acd.rois.size() == npul );
  }

  cout << myname << line << endl;

  cout << myname << line << endl;
  cout << myname << "Call tool for channel map." << endl;
  DataMap res;
  res = ptoo->viewMap(acds);
  res = ptoo->viewMap(acds);
  res = ptoo->viewMap(acds);

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
  return test_AdcEventViewer(useExistingFcl);
}

//**********************************************************************
