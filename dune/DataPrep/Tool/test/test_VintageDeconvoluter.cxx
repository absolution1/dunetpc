// test_VintageDeconvoluter.cxx
//
// David Adams
// May 2018
//
// Test VintageDeconvoluter.

#include <string>
#include <iostream>
#include <fstream>
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "dune/DuneInterface/Tool/TpcDataTool.h"
#include "dune/ArtSupport/DuneToolManager.h"
#include "dune/ArtSupport/ArtServiceHelper.h"

#undef NDEBUG
#include <cassert>

using std::string;
using std::cout;
using std::endl;
using std::ofstream;
using fhicl::ParameterSet;
using Index = unsigned int;

//**********************************************************************

int test_VintageDeconvoluter(bool useExistingFcl =false) {
  const string myname = "test_VintageDeconvoluter: ";
#ifdef NDEBUG
  cout << myname << "NDEBUG must be off." << endl;
  abort();
#endif
  string line = "-----------------------------";

  cout << myname << line << endl;
  string fclfile = "test_VintageDeconvoluter.fcl";
  if (useExistingFcl) {
    cout << myname << "Using existing top-level FCL." << endl;
  } else {
    cout << myname << "Creating top-level FCL." << endl;
    std::ofstream config{fclfile};
    config << "#include \"services_dune.fcl\"" << endl;
    //config << "services: @local::protodune_reco_services" << endl;
    config << "services: @local::dune35t_services_legacy" << endl;
    config << "tools: {" << endl;
    config << "  mytool: {" << endl;
    config << "    tool_type: VintageDeconvoluter" << endl;
    config << "    LogLevel: 1" << endl;
    config << "  }" << endl;
    config << "}" << endl;
  }
  ArtServiceHelper::load_services(fclfile, ArtServiceHelper::FileOnPath);

  cout << myname << line << endl;
  cout << myname << "Fetching tool manager." << endl;
  DuneToolManager* ptm = DuneToolManager::instance(fclfile);
  assert ( ptm != nullptr );
  DuneToolManager& tm = *ptm;
  tm.print();
  assert( tm.toolNames().size() == 1 );

  cout << myname << line << endl;
  cout << myname << "Fetching tool." << endl;
  auto psgf = tm.getPrivate<TpcDataTool>("mytool");
  assert( psgf != nullptr );
  auto psgfmod = tm.getPrivate<TpcDataTool>("mytool");
  assert( psgfmod != nullptr );

  cout << myname << line << endl;
  cout << myname << "Create data and call tool." << endl;
  AdcChannelData data;
  data.setChannelInfo(123);
  for ( AdcIndex itic=0; itic<100; ++itic ) {
    float xadc = rand()%20 - 10.0;
    data.samples.push_back(xadc);
  }
  data.samples[30] = 150.0;
  assert( data.signal.size() == 0 );
  assert( data.rois.size() == 0 );
  assert( data.samples[30] = 150 );
  AdcChannelData data0 = data;
  data0.signal  = data.signal;
  data0.rois    = data.rois;
  data0.samples = data.samples;

  cout << myname << line << endl;
  cout << myname << "Running tool." << endl;
  DataMap resmod = psgfmod->update(data);
  resmod.print();

  cout << myname << line << endl;
  cout << myname << "Checking results." << endl;
  assert( data.samples.size() == data0.samples.size() );
  for ( Index isam=0; isam<data0.samples.size(); ++isam ) {
    cout << isam << ": " << data0.samples[isam] << "  " << data.samples[isam] << endl;
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
  return test_VintageDeconvoluter(useExistingFcl);
}

//**********************************************************************
