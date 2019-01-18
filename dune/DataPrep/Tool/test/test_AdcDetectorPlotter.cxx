// test_AdcDetectorPlotter.cxx
//
// David Adams
// April 2018
//
// Test AdcDetectorPlotter.

#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include "dune/DuneInterface/Tool/AdcChannelTool.h"
#include "dune/ArtSupport/DuneToolManager.h"
#include "dune/ArtSupport/ArtServiceHelper.h"
#include "dune/Geometry/AdcChannelDataTester.h"
#include "larcore/Geometry/Geometry.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "TRandom.h"

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

int test_AdcDetectorPlotter(bool useExistingFcl =false) {
  const string myname = "test_AdcDetectorPlotter: ";
#ifdef NDEBUG
  cout << myname << "NDEBUG must be off." << endl;
  abort();
#endif
  string line = "-----------------------------";

  cout << myname << line << endl;
  string fclfile = "test_AdcDetectorPlotter.fcl";
  string gname = "protodune_geo";
  if (useExistingFcl) {
    cout << myname << "Using existing top-level FCL." << endl;
  } else {
    cout << myname << "Creating top-level FCL." << endl;
    std::ofstream config{fclfile};
    config << "#include \"geometry_dune.fcl\"" << endl;
    config << "services.Geometry:                   @local::" + gname << endl;
    config << "services.ExptGeoHelperInterface:     @local::dune_geometry_helper" << endl;
    config << "#include \"dataprep_tools.fcl\"" << endl;  // Need adcStringBuilder
    config << "tools.mytool: {" << endl;
    config << "    tool_type: AdcDetectorPlotter" << endl;
    config << "    LogLevel: 2" << endl;
    config << "   WireAngle: 0.0" << endl;
    //config << "   WireAngle: 0.623" << endl;
    config << "    DataType: 0" << endl;
    config << "       Tick0: 0" << endl;
    config << "  DriftSpeed: 0.01" << endl;
    config << "    XMin:  600.0" << endl;
    config << "    XMax: -600.0" << endl;
    config << "    ZMin:  -20.0" << endl;
    config << "    ZMax:  720.0" << endl;
    config << "    SignalThreshold: 10" << endl;
    config << "    SkipBadChannels: false" << endl;
    config << "    ShowAllTicks: false" << endl;
    config << "    FirstTick: 0" << endl;
    config << "     LastTick: 0" << endl;
    config << "    ShowWires: true" << endl;
    config << "    ShowCathode: true" << endl;
    config << "    ShowTpcSets: []" << endl;
    config << "    ShowGrid: true" << endl;
    config << "    Title: \"Prepared ADC run %RUN% event %EVENT%\"" << endl;
    config << "    PlotTitle: \"%UTCTIME2\"" << endl;
    config << "    FileName: \"test_AdcDetectorPlotter-run%0RUN%-evt%0EVENT%.png\"" << endl;
    config << "}" << endl;
  }
  ArtServiceHelper::load_services(fclfile, ArtServiceHelper::FileOnPath);

  cout << myname << line << endl;
  cout << myname << "Get Geometry service." << endl;
  art::ServiceHandle<geo::Geometry> pgeo;
  cout << myname << "Geometry name: " << pgeo->DetectorName() << endl;

  cout << myname << line << endl;
  cout << myname << "Fetching tool manager." << endl;
  DuneToolManager* ptm = DuneToolManager::instance(fclfile);
  assert ( ptm != nullptr );
  DuneToolManager& tm = *ptm;
  tm.print();
  assert( tm.toolNames().size() >= 1 );

  cout << myname << line << endl;
  cout << myname << "Fetching tool." << endl;
  auto padv = tm.getPrivate<AdcChannelTool>("mytool");
  assert( padv != nullptr );

  cout << myname << line << endl;
  cout << myname << "Create data." << endl;
  Index nevt = 2;
  vector<vector<AdcChannelDataMap>> evtdatamaps(nevt);
  Index icry = 0;
  const geo::CryostatGeo& gcry = pgeo->Cryostat(icry);
  Index ntpc = gcry.NTPC();
  Index ncha = pgeo->Nchannels();
  // Loop over events.
  AdcChannelDataTester tester;
  tester.run = 123;
  tester.event = 0;
  tester.nsam = 20000;
  bool strumWires = false;
  // Loop over events.
  for ( vector<AdcChannelDataMap>& datamaps : evtdatamaps ) {
    if ( strumWires ) {
      datamaps.resize(ntpc);
      // Loop over TPCs.  No, no, no...
      for ( Index itpc=0; itpc<ntpc; ++itpc ) {
        tester.strumTpcWires(datamaps[itpc], itpc, icry);
      }
    } else {
      datamaps.resize(1);
      AdcChannelDataMap& datamap = datamaps[0];
      for ( Index icha=0; icha<ncha; ++ icha ) datamap[icha];
      tester.strumChannels(datamap);
    }
    Index ngrp = datamaps.size();
    cout << myname << "  Event " << tester.event << endl;
    cout << myname << "    Group count is " << ngrp << endl;
    for ( Index igrp=0; igrp<ngrp; ++ igrp ) {
      cout << myname << "      Group " << igrp << " channel count: " << datamaps[igrp].size() << endl;
    }
    ++tester.event;
  }

  cout << myname << line << endl;
  cout << myname << "Call tool." << endl;
  assert( nevt == 2 );
  // Loop over events.
  for ( AdcIndex ievt=0; ievt<nevt; ++ievt ) {
    cout << myname << "Event " << ievt << endl;
    const vector<AdcChannelDataMap>& datamaps = evtdatamaps[ievt];
    Index ngrp = datamaps.size();
    // Loop over ADC channel groups (1 group per TPC).
    for ( Index igrp=0; igrp<ngrp; ++igrp ) {
      const AdcChannelDataMap& datamap = datamaps[igrp];
      assert( datamap.size() > 0 );
      assert( datamap.begin()->second.event == ievt );
      assert( datamap.begin()->second.run == tester.run );
      cout << myname << "Calling tool for event " << ievt << ", group " << igrp
           << ", Nchan = " << datamap.size() << endl;
      assert( padv->viewMap(datamap) == 0 );
    }
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
  return test_AdcDetectorPlotter(useExistingFcl);
}

//**********************************************************************
