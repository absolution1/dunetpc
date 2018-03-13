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
  string fclfile = "test_Dune35tNoiseRemovalService.fcl";
  if ( ! useExistingFcl ) {
    cout << myname << "Creating top-level FCL." << endl;
    string gname = "dune35t4apa_v6";
    ofstream fout(fclfile.c_str());
    fout << "services.Geometry: {" << endl;
    fout << "  DisableWiresInG4: true" << endl;
    fout << "  GDML: \"" << gname << ".gdml\"" << endl;
    fout << "  Name: \"" << gname << "\"" << endl;
    fout << "  ROOT: \"" << gname << "\"" << endl;
    fout << "  SortingParameters: { DetectorVersion: \"" << gname << "\" ChannelsPerOpDet: 12} " << endl;
    fout << "  SurfaceY: 0" << endl;
    fout << "}" << endl;
    fout << "services.ExptGeoHelperInterface: {" << endl;
    fout << "  service_provider: \"DUNEGeometryHelper\"" << endl;
    fout << "}" << endl;
    fout << "services.ChannelMapService: {" << endl;
    fout << "  LogLevel: 1" << endl;
    fout << "  FileName: \"35tTPCChannelMap_v6.txt\"" << endl;
    fout << "}" << endl;
    fout << "services.AdcNoiseRemovalService: {" << endl;
    fout << "  service_provider: Dune35tNoiseRemovalService" << endl;
    fout << "           LogLevel: 1" << endl;
    fout << "       GroupingFlag: 1" << endl;
    fout << "     SkipStuckCodes: false" << endl;
    fout << "        SkipSignals: false" << endl;
    fout << "  CorrectStuckCodes: true" << endl;
    fout << "         ShowGroups: 2" << endl;
    fout << "  ShowGroupsChannel: 4" << endl;
    fout << "}" << endl;
    fout.close();
  } else {
    cout << myname << "Using existing top-level FCL." << endl;
  }

  cout << myname << "Fetch art service helper." << endl;
  ArtServiceHelper& ash = ArtServiceHelper::instance();
  ash.print();

  if ( ash.serviceStatus() == 0 ) {

    cout << myname << line << endl;
    cout << myname << "Add the Geometry services." << endl;
    assert( ash.addService("Geometry", fclfile, true) == 0 );
    assert( ash.addService("ExptGeoHelperInterface", fclfile, true) == 0 );

    cout << myname << line << endl;
    cout << myname << "Add LBNE channel map service." << endl;
    assert( ash.addService("ChannelMapService", fclfile, true) == 0 );

    cout << myname << line << endl;
    cout << myname << "Add noise removal service." << endl;
    assert( ash.addService("AdcNoiseRemovalService", fclfile, true) == 0 );

    cout << myname << line << endl;
    cout << myname << "Load services." << endl;
    assert( ash.loadServices() == 1 );
    ash.print();
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
  cout << myname << "Close service helper." << endl;
  ArtServiceHelper::close();

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
