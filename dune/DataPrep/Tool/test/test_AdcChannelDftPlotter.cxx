// test_AdcChannelDftPlotter.cxx
//
// David Adams
// April 2017
//
// Test AdcChannelDftPlotter.

#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
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

using Index = unsigned int;

//**********************************************************************

int test_AdcChannelDftPlotter(bool useExistingFcl =false) {
  const string myname = "test_AdcChannelDftPlotter: ";
#ifdef NDEBUG
  cout << myname << "NDEBUG must be off." << endl;
  abort();
#endif
  string line = "-----------------------------";

  cout << myname << line << endl;
  string fclfile = "test_AdcChannelDftPlotter.fcl";
  if ( ! useExistingFcl ) {
    cout << myname << "Creating top-level FCL." << endl;
    ofstream fout(fclfile.c_str());
    fout << "#include \"dataprep_tools.fcl\"" << endl;
    fout << "tools.myphases: {" << endl;
    fout << "    tool_type: AdcChannelDftPlotter" << endl;
    fout << "     LogLevel: 2" << endl;
    fout << "     Variable: phase" << endl;
    fout << "   SampleFreq: 0" << endl;
    fout << "     HistName: \"hdftphase_run%RUN%_evt%EVENT%\"" << endl;
    fout << "    HistTitle: \"DFT phases for run %RUN% event %EVENT%\"" << endl;
    fout << "     PlotName: \"dftphase_run%RUN%_evt%EVENT%.png\"" << endl;
    fout << "}" << endl;
    fout << "tools.mymags: {" << endl;
    fout << "    tool_type: AdcChannelDftPlotter" << endl;
    fout << "     LogLevel: 2" << endl;
    fout << "     Variable: magnitude" << endl;
    fout << "   SampleFreq: 0" << endl;
    fout << "         YMax: 0.0" << endl;
    fout << "     HistName: \"hdftmags_run%RUN%_evt%EVENT%\"" << endl;
    fout << "    HistTitle: \"DFT amplitudes for run %RUN% event %EVENT%\"" << endl;
    fout << "     PlotName: \"dftmag_run%RUN%_evt%EVENT%.png\"" << endl;
    fout << "}" << endl;
    fout << "tools.mypwr: {" << endl;
    fout << "    tool_type: AdcChannelDftPlotter" << endl;
    fout << "     LogLevel: 2" << endl;
    fout << "     Variable: power" << endl;
    fout << "   SampleFreq: 2000" << endl;
    fout << "         YMax: 50.0" << endl;
    fout << "        NBinX: 5" << endl;
    fout << "     HistName: \"hdftpower_run%RUN%_evt%EVENT%\"" << endl;
    fout << "    HistTitle: \"DFT power for run %RUN% event %EVENT%\"" << endl;
    fout << "     PlotName: \"dftpower_run%RUN%_evt%EVENT%.png\"" << endl;
    fout << "}" << endl;
    fout << "tools.mypwt: {" << endl;
    fout << "    tool_type: AdcChannelDftPlotter" << endl;
    fout << "     LogLevel: 2" << endl;
    fout << "     Variable: \"power/tick\"" << endl;
    fout << "   SampleFreq: 2000" << endl;
    fout << "         YMax: 5.0" << endl;
    fout << "        NBinX: 5" << endl;
    fout << "     HistName: \"hdftpowt_run%RUN%_evt%EVENT%\"" << endl;
    fout << "    HistTitle: \"DFT power for run %RUN% event %EVENT%\"" << endl;
    fout << "     PlotName: \"dftpowt_run%RUN%_evt%EVENT%.png\"" << endl;
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
  assert( tm.toolNames().size() >= 1 );

  cout << myname << line << endl;
  cout << myname << "Fetching tools." << endl;
  auto ppha = tm.getPrivate<AdcChannelTool>("myphases");
  auto pmag = tm.getPrivate<AdcChannelTool>("mymags");
  auto ppwr = tm.getPrivate<AdcChannelTool>("mypwr");
  auto ppwt = tm.getPrivate<AdcChannelTool>("mypwt");
  assert( ppha != nullptr );
  assert( pmag != nullptr );
  assert( ppwr != nullptr );
  assert( ppwt != nullptr );

  cout << myname << line << endl;
  cout << myname << "Create data." << endl;
  AdcChannelData acd;
  acd.channel = 100123;
  acd.run = 123;
  acd.subRun = 45;
  acd.event = 2468;
  acd.sampleUnit = "fC";
  vector<float> mags = {  1.0, 2.0, 3.0, 4.0, 5.0, 4.0, 2.0,  1.0,  1.0,  1.0,  1.0 };
  vector<float> phas = {  0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, -2.0, -1.5, -1.0 };
  acd.dftmags = mags;
  acd.dftphases = phas;
  cout << myname << "  Mag size: " << mags.size() << endl;
  cout << myname << "Phase size: " << phas.size() << endl;

  cout << myname << line << endl;
  cout << myname << "Call phase tool." << endl;
  DataMap dm = ppha->view(acd);
  dm.print();
  assert( dm == 0 );

  cout << myname << line << endl;
  cout << myname << "Call mag tool." << endl;
  dm = pmag->view(acd);
  dm.print();
  assert( dm == 0 );

  cout << myname << line << endl;
  cout << myname << "Call power tool." << endl;
  dm = ppwr->view(acd);
  dm.print();
  assert( dm == 0 );

  cout << myname << line << endl;
  cout << myname << "Call tick power tool." << endl;
  dm = ppwt->view(acd);
  dm.print();
  assert( dm == 0 );

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
  return test_AdcChannelDftPlotter(useExistingFcl);
}

//**********************************************************************
