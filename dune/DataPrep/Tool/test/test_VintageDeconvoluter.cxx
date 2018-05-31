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
#include "dune/DuneInterface/Tool/AdcChannelTool.h"
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
  if ( ! useExistingFcl ) {
    cout << myname << "Creating top-level FCL." << endl;
    ofstream fout(fclfile.c_str());
    fout << "#include \"services_dune.fcl\"" << endl;
    //fout << "services: @local::protodune_reco_services" << endl;
    fout << "services: @local::dune35t_services" << endl;
    fout << "tools: {" << endl;
    fout << "  mytool: {" << endl;
    fout << "    tool_type: VintageDeconvoluter" << endl;
    fout << "    LogLevel: 1" << endl;
    fout << "  }" << endl;
    fout << "}" << endl;
    fout.close();
  } else {
    cout << myname << "Using existing top-level FCL." << endl;
  }

  cout << myname << line << endl;
  cout << myname << "Fetch art service helper." << endl;
  ArtServiceHelper& ash = ArtServiceHelper::instance();
  assert( ash.addServices(fclfile, true) == 0 );
  assert( ash.loadServices() == 1 );
  ash.print();

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
  data.channel = 123;
  for ( AdcIndex itic=0; itic<100; ++itic ) {
    float xadc = rand()%20 - 10.0;
    data.samples.push_back(xadc);
  }
  data.samples[30] = 150.0;
  assert( data.signal.size() == 0 );
  assert( data.rois.size() == 0 );
  assert( data.samples[30] = 150 );
  AdcChannelData data0 = data;

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
  cout << "Closing service helper." << endl;
  ArtServiceHelper::close();    // See https://cdcvs.fnal.gov/redmine/issues/10618

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
