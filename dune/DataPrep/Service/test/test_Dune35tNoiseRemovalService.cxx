// test_Dune35tNoiseRemovalService.cxx
//
// David Adams
// May 2016
//
// Test Dune35tNoiseRemovalService.

#include <string>
#include <iostream>
#include <sstream>
#include <fstream>
#include <iomanip>
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "dune/ArtSupport/ArtServiceHelper.h"
#include "lbne-raw-data/Services/ChannelMap/ChannelMapService.h"
#include "dune/DuneInterface/AdcNoiseRemovalService.h"

#undef NDEBUG
#include <cassert>

using std::string;
using std::vector;
using std::cout;
using std::endl;
using std::istringstream;
using std::ofstream;
using std::setw;
using std::setprecision;
using std::fixed;
using art::ServiceHandle;

//**********************************************************************

int test_Dune35tNoiseRemovalService(bool useExistingFcl) {
  const string myname = "test_Dune35tNoiseRemovalService: ";
#ifdef NDEBUG
  cout << myname << "NDEBUG must be off." << endl;
  abort();
#endif
  string line = "-----------------------------";

  cout << myname << line << endl;
  if (useExistingFcl) {
    cout << myname << "Using existing top-level FCL." << endl;
    ArtServiceHelper::load_services("test_Dune35tNoiseRemovalService.fcl",
                                    ArtServiceHelper::FileOnPath);
  } else {
    cout << myname << "Creating top-level FCL." << endl;
    string gname = "dune35t4apa_v6";
    std::stringstream config;
    config << "services.Geometry: {" << endl;
    config << "  DisableWiresInG4: true" << endl;
    config << "  GDML: \"" << gname << ".gdml\"" << endl;
    config << "  Name: \"" << gname << "\"" << endl;
    config << "  ROOT: \"" << gname << "\"" << endl;
    config << "  SortingParameters: { DetectorVersion: \"" << gname << "\" ChannelsPerOpDet: 12} " << endl;
    config << "  SurfaceY: 0" << endl;
    config << "}" << endl;
    config << "services.ExptGeoHelperInterface: {" << endl;
    config << "  service_provider: \"DUNEGeometryHelper\"" << endl;
    config << "}" << endl;
    config << "services.ChannelMapService: {" << endl;
    config << "  LogLevel: 1" << endl;
    config << "  FileName: \"35tTPCChannelMap_v6.txt\"" << endl;
    config << "}" << endl;
    config << "services.AdcNoiseRemovalService: {" << endl;
    config << "  service_provider: Dune35tNoiseRemovalService" << endl;
    config << "           LogLevel: 1" << endl;
    config << "       GroupingFlag: 1" << endl;
    config << "     SkipStuckCodes: false" << endl;
    config << "        SkipSignals: false" << endl;
    config << "  CorrectStuckCodes: true" << endl;
    config << "         ShowGroups: 2" << endl;
    config << "  ShowGroupsChannel: 4" << endl;
    config << "}" << endl;
    ArtServiceHelper::load_services(config);
  }

  cout << myname << line << endl;
  cout << myname << "Fetch channel map." << endl;
  ServiceHandle<lbne::ChannelMapService> hlcm;

  cout << myname << line << endl;
  cout << myname << "Fetch noise removal." << endl;
  ServiceHandle<AdcNoiseRemovalService> hanr;
  hanr->print();

  // Each channel is filled with a bipolar signal offset by one from preceding channel.
  cout << myname << line << endl;
  cout << myname << "Build offline channel list." << endl;
  cout << myname << "Create data." << endl;
  unsigned int nchan = 64;
  unsigned int nsig = 200;
  float fac = 25.0;
  AdcChannelDataMap datamap;
  for ( AdcChannel chan=0; chan<nchan; ++chan ) {
    cout << "  Online Channel: " << chan << endl;
    AdcChannel chanoff = hlcm->Offline(chan);
    cout << " Offline Channel: " << chanoff << endl;
    assert( hlcm->Online(chanoff) == chan );
    AdcChannelData& data = datamap[chanoff];
    data.channel = chanoff;
    AdcSignalVector& sigs = data.samples;
    //AdcFlagVector& flags = data.flags;
    unsigned int isig1 = 10 + chan;
    for ( unsigned int isig=0; isig<isig1; ++isig ) sigs.push_back(0);
    for ( unsigned int i=0; i<10; ++i )             sigs.push_back(fac*i);
    for ( unsigned int i=10; i<1000; --i )          sigs.push_back(fac*i);
    for ( unsigned int i=19; i<1000; --i )          sigs.push_back(-sigs[i+isig1]);
    for ( unsigned int isig=sigs.size(); isig<nsig; ++isig ) sigs.push_back(0);
  }
  // Add noise.
  for ( auto& chdata : datamap ) {
    AdcChannel chanoff = chdata.first;
    AdcChannelData& data = chdata.second;
    int reg = hlcm->RegulatorFromOfflineChannel(chanoff);
    int ori = hlcm->PlaneFromOfflineChannel(chanoff);
    for ( unsigned int isig=0; isig<nsig; ++isig ) {
      AdcSignal& sig = data.samples[isig];
      AdcSignal noise = 50.0*reg + 10*ori + isig%5;
      sig += noise;
    }
  }

  cout << myname << line << endl;
  cout << myname << "Data before noise removal:" << endl;
  unsigned int wcha = 5;
  unsigned int wsig = 5;
  for ( AdcChannelDataMap::value_type& idata : datamap ) {
    AdcChannel chan = idata.first;
    AdcChannelData& data = idata.second;
    cout << setw(wcha) << chan << ": ";
    for ( unsigned int isig=0; isig<40; ++isig ) {
      cout << fixed << setprecision(0) << setw(wsig) << data.samples[isig];
    }
    cout << endl;
  }

  cout << myname << line << endl;
  cout << myname << "Remove noise." << endl;
  assert( hanr->update(datamap) == 0 );

  cout << myname << line << endl;
  cout << myname << "Data after noise removal:" << endl;
  for ( AdcChannelDataMap::value_type& idata : datamap ) {
    AdcChannel chan = idata.first;
    AdcChannelData& data = idata.second;
    cout << setw(wcha) << chan << ": ";
    for ( unsigned int isig=0; isig<40; ++isig ) {
      cout << fixed << setprecision(0) << setw(wsig) << data.samples[isig];
    }
    cout << endl;
  }

  cout << myname << line << endl;
  cout << myname << "Check unpolluted data after noise removal:" << endl;
  for ( auto& chdata : datamap ) {
    AdcChannelData& data = chdata.second;
    for ( unsigned int isig=0; isig<10; ++isig ) {
      assert(data.samples[isig] == 0.0);
    }
  }

  cout << myname << line << endl;
  cout << myname << "Done." << endl;
  return 0;
}

//**********************************************************************

int main(int argc, char* argv[]) {
  bool useExistingFcl = false;
  bool skipTest = false;
  if ( argc > 1 ) {
    string sarg(argv[1]);
    if ( sarg == "-h" ) {
      cout << "Usage: " << argv[0] << " [ARG]" << endl;
      cout << "  If ARG = true, existing FCL file is used." << endl;
      return 0;
    }
    useExistingFcl = sarg == "true" || sarg == "1";
  } else if ( skipTest ) {
    cout << "Skipping test_Dune35tNoiseRemovalService to avoid hang." << endl;
    cout << "See https://cdcvs.fnal.gov/redmine/issues/19206" << endl;
    return 0;
  }
  return test_Dune35tNoiseRemovalService(useExistingFcl);
}

//**********************************************************************
