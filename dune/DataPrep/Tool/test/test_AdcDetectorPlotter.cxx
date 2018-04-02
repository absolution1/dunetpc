// test_AdcDetectorPlotter.cxx
//
// David Adams
// April 2017
//
// Test AdcDetectorPlotter.

#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include "dune/DuneInterface/Tool/AdcDataViewer.h"
#include "dune/ArtSupport/DuneToolManager.h"
#include "dune/ArtSupport/ArtServiceHelper.h"
#include "larcore/Geometry/Geometry.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
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
  if ( ! useExistingFcl ) {
    cout << myname << "Creating top-level FCL." << endl;
    ofstream fout(fclfile.c_str());
    fout << "#include \"geometry_dune.fcl\"" << endl;
    fout << "services.Geometry:                   @local::" + gname << endl;
    fout << "services.ExptGeoHelperInterface:     @local::dune_geometry_helper" << endl;
    fout << "tools: {" << endl;
    fout << "  mytool: {" << endl;
    fout << "    tool_type: AdcDetectorPlotter" << endl;
    fout << "    LogLevel: 2" << endl;
    fout << "   WireAngle: 0.0" << endl;
    fout << "    DataType: 0" << endl;
    fout << "    XMin:    0.0" << endl;
    fout << "    XMax:  900.0" << endl;
    fout << "    ZMin: -400.0" << endl;
    fout << "    ZMax:  400.0" << endl;
    fout << "    SignalThreshold: 10" << endl;
    fout << "    Title: \"Prepared ADC run %RUN% event %EVENT%\"" << endl;
    fout << "    FileName: \"test_AdcDetectorPlotter-run%RUN%-evt%EVENT%.png\"" << endl;
    fout << "  }" << endl;
    fout << "}" << endl;
    fout.close();
  } else {
    cout << myname << "Using existing top-level FCL." << endl;
  }

  cout << myname << line << endl;
  cout << myname << "Fetch art service helper." << endl;
  ArtServiceHelper& myash = ArtServiceHelper::instance();
  myash.setLogLevel(3);

  cout << myname << line << endl;
  cout << myname << "Add services from " << fclfile << endl;
  assert( myash.addServices(fclfile, true) == 0 );

  cout << myname << line << endl;
  cout << myname << "Display services" << endl;
  myash.print();

  cout << myname << line << endl;
  cout << myname << "Load the services." << endl;
  assert( myash.loadServices() == 1 );
  myash.print();

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
  assert( tm.toolNames().size() == 1 );

  cout << myname << line << endl;
  cout << myname << "Fetching tool." << endl;
  auto padv = tm.getPrivate<AdcDataViewer>("mytool");
  assert( padv != nullptr );

  cout << myname << line << endl;
  cout << myname << "Create data." << endl;
  AdcIndex nevt = 2;
  float peds[20] = {701, 711, 733, 690, 688, 703, 720, 720, 695, 702,
                    410, 404, 388, 389, 400, 401, 410, 404, 395, 396};
  vector<double> wf = {5.0, 20.1, 53.2, 80.6, 130.2, 160.1, 150.4, 125.7, 72.5, 41.3, 18.4,
                       -6.5, -34.9, -56.6, -88.9, -132.6, -170.8, -172.9, -144.6, -112.6,
                        -79.4, -44.9, -22.1, -12.6, -4.7};
  vector<vector<AdcChannelDataMap>> evtdatamaps(2);
  for ( vector<AdcChannelDataMap>& datamaps : evtdatamaps ) datamaps.resize(4);
  for ( AdcIndex ievt=0; ievt<nevt; ++ievt ) {
    cout << myname << "Event " << ievt << endl;
    assert( ievt < evtdatamaps.size() );
    vector<AdcChannelDataMap>& datamaps = evtdatamaps[ievt];
    AdcIndex ncha = 20;
    AdcIndex icha1 = 10000;
    for ( AdcIndex kcha=0; kcha<ncha; ++kcha ) {
      AdcIndex icha = icha1 + kcha;
      Index igrp = kcha/5;
      cout << myname << "  Channel " << icha << ", group " << igrp << endl;
      assert( igrp < datamaps.size() );
      AdcChannelDataMap& datamap = datamaps[igrp];
      std::pair<AdcChannelDataMap::iterator, bool> kdat = datamap.emplace(icha, AdcChannelData());
      assert(kdat.second);
      AdcChannelDataMap::iterator idat = kdat.first;
      AdcChannelData& data = idat->second;
      data.run = 123;
      data.event = ievt;
      float ped = peds[icha-icha1];
      data.channel = icha;
      data.pedestal = ped;
      for ( AdcIndex itic=0; itic<100; ++itic ) {
        float xadc = ped + gRandom->Gaus(0.0, 10.0);
        AdcCount iadc = xadc;
        data.raw.push_back(iadc);
        data.samples.push_back(iadc - ped);
      }
      AdcIndex tp = 10*ievt + 60 - 2.3*(icha-icha1);
      for ( unsigned int iwf=0; iwf<wf.size(); ++iwf ) {
        unsigned int isam = tp+iwf;
        if ( isam < data.samples.size() ) {
          data.raw[isam] += wf[iwf];
          data.samples[isam] += wf[iwf];
        }
      }
      for ( unsigned int isam=0; isam<data.samples.size(); ++isam ) {
        data.samples[isam] *= 0.04;
      }
      data.sampleUnit = "ke";
    }
    ostringstream sslab;
    sslab << "event " << ievt << " plane 3u";
    string slab = sslab.str();
    string fpat = slab;
    for ( char& ch : fpat ) if ( ch == ' ' ) ch = '-';
  }

  cout << myname << line << endl;
  cout << myname << "Call tool." << endl;
  string slab = "my data";
  string fpat = "mydata";
  assert( nevt == 2 );
  for ( AdcIndex ievt=0; ievt<nevt; ++ievt ) {
    cout << myname << "Event " << ievt << endl;
    const vector<AdcChannelDataMap>& datamaps = evtdatamaps[ievt];
    Index ngrp = datamaps.size();
    assert( ngrp == 4 );
    for ( Index igrp=0; igrp<ngrp; ++igrp ) {
      const AdcChannelDataMap& datamap = datamaps[igrp];
      assert(datamap.size() == 5 );
      assert( padv->view(datamap, slab, fpat) == 0 );
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
