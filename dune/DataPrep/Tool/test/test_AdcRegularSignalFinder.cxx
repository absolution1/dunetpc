// test_AdcRegularSignalFinder.cxx
//
// David Adams
// April 2017
//
// Test AdcRegularSignalFinder.

#include <string>
#include <iostream>
#include <fstream>
#include <iomanip>
#include "dune/DuneInterface/Tool/AdcChannelTool.h"
#include "dune/DuneInterface/Tool/AdcChannelTool.h"
#include "dune/ArtSupport/DuneToolManager.h"

#undef NDEBUG
#include <cassert>

using std::string;
using std::cout;
using std::endl;
using std::ofstream;
using std::setw;
using fhicl::ParameterSet;
using Index = AdcIndex;

//**********************************************************************

int test_AdcRegularSignalFinder(bool useExistingFcl =false) {
  const string myname = "test_AdcRegularSignalFinder: ";
#ifdef NDEBUG
  cout << myname << "NDEBUG must be off." << endl;
  abort();
#endif
  string line = "-----------------------------";

  cout << myname << line << endl;
  string fclfile = "test_AdcRegularSignalFinder.fcl";
  if ( ! useExistingFcl ) {
    cout << myname << "Creating top-level FCL." << endl;
    ofstream fout(fclfile.c_str());
    fout << "tools: {" << endl;
    fout << "  mytool: {" << endl;
    fout << "    tool_type: AdcRegularSignalFinder" << endl;
    fout << "    LogLevel: 1" << endl;
    fout << "    Period: 30" << endl;
    fout << "    Length: 20" << endl;
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
  Index nsam = 100;
  for ( AdcIndex itic=0; itic<nsam; ++itic ) {
    float xadc = rand()%20 - 10.0;
    data.samples.push_back(xadc);
  }
  data.samples[30] = 150.0;
  assert( data.signal.size() == 0 );
  assert( data.rois.size() == 0 );
  assert( data.samples.size() == nsam );
  assert( data.samples[30] = 150 );

  cout << myname << line << endl;
  cout << myname << "Running tool." << endl;
  DataMap resmod = psgfmod->update(data);
  resmod.print();

  cout << myname << line << endl;
  cout << myname << "Checking results." << endl;
  Index per = 30;
  Index len = 20;
  Index nroi = 4;
  
  assert( resmod == 0 );
  assert( resmod.getInt("roiPeriod") == int(per) );
  assert( resmod.getInt("roiLength") == int(len) );
  assert( resmod.getInt("roiCount") == int(nroi) );
  assert( data.signal.size() == nsam );
  assert( data.rois.size() == nroi );
  for ( Index iroi=0; iroi<nroi; ++iroi ) {
    AdcRoi roi = data.rois[iroi];
    cout << myname << "ROI " << iroi << ": ["
         << roi.first << ", "
         << roi.second << "]" << endl;
    Index isam1 = iroi*per;
    Index isam2 = isam1 + len;
    if ( isam2 > nsam ) isam2 = nsam;
    Index isam3 = isam1 + per;
    if ( isam3 > nsam ) isam3 = nsam;
    for ( Index isam=isam1; isam<isam3; ++isam ) {
      bool sig = data.signal[isam];
      string ssig = sig ? "true" : "false";
      cout << myname << setw(4) << isam << ": " << ssig << endl;
      if ( isam < isam2 ) assert( data.signal[isam]);
      else assert( ! data.signal[isam]);
    }
    assert( roi.first == isam1 );
    assert( roi.second == isam2 - 1 );
  }
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
  return test_AdcRegularSignalFinder(useExistingFcl);
}

//**********************************************************************
