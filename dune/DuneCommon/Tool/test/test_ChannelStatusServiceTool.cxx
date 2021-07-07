// test_ChannelStatusServiceTool.cxx
//
// David Adams
// Novenber 2019
//
// Test ChannelStatusServiceTool.

#include "dune/DuneInterface/Data/AdcTypes.h"
#include "dune/ArtSupport/ArtServiceHelper.h"
#include "dune/ArtSupport/DuneToolManager.h"
#include "dune/DuneInterface/Tool/IndexMapTool.h"
#include "TH1F.h"
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>

#undef NDEBUG
#include <cassert>

using std::string;
using std::cout;
using std::endl;
using std::ofstream;
using std::istringstream;
using std::setw;
using std::vector;
using fhicl::ParameterSet;
using Index = unsigned int;

//**********************************************************************

int test_ChannelStatusServiceTool(bool useExistingFcl =false, string extrafcl ="") {
  const string myname = "test_ChannelStatusServiceTool: ";
#ifdef NDEBUG
  cout << myname << "NDEBUG must be off." << endl;
  abort();
#endif
  string line = "-----------------------------";

  cout << myname << line << endl;
  string fclfile = "test_ChannelStatusServiceTool.fcl";
  if ( ! useExistingFcl ) {
    cout << myname << "Creating top-level FCL." << endl;
    ofstream fout(fclfile.c_str());
    fout << "#include \"services_dune.fcl\"" << endl;
    fout << "services.Geometry: @local::protodune_geo" << endl;
    fout << "services.ExptGeoHelperInterface: @local::dune_geometry_helper" << endl;
    fout << "services.ChannelStatusService: {" << endl;
    fout << "  service_provider: SimpleChannelStatusService" << endl;
    fout << "  BadChannels: [ 1 ]" << endl;
    fout << "  NoisyChannels: [3, 4]" << endl;
    fout << "}" << endl;
    fout << "tools: {" << endl;
    fout << "  mytool: {" << endl;
    fout << "    tool_type: ChannelStatusServiceTool" << endl;
    fout << "    LogLevel: 1" << endl;
    fout << "  }" << endl;
    fout << "}" << endl;
    if ( extrafcl.size() ) {
      fout << "#include \"" + extrafcl + "\"" << endl;
    }
    fout.close();
  } else {
    cout << myname << "Using existing top-level FCL." << endl;
  }

  cout << myname << line << endl;
  cout << myname << "Load services." << endl;
  ArtServiceHelper::load(fclfile);

  cout << myname << line << endl;
  cout << myname << "Fetching tool manager." << endl;
  DuneToolManager* ptm = DuneToolManager::instance(fclfile);
  assert ( ptm != nullptr );
  DuneToolManager& tm = *ptm;
  tm.print();
  assert( tm.toolNames().size() >= 1 );

  cout << myname << line << endl;
  cout << myname << "Fetching tool." << endl;
  auto pcs = tm.getPrivate<IndexMapTool>("mytool");
  assert( pcs != nullptr );

  cout << myname << line << endl;
  cout << "Check some channels." << endl;
  vector<Index> chks(20, 0);
  chks[1] = 1;
  chks[3] = 2;
  chks[4] = 2;
  for ( Index icha=0; icha<20; ++icha ) {
    Index chanStat = pcs->get(icha);
    cout << setw(8) << icha << ": " << setw(2) << chanStat << endl;
    assert( chanStat == chks[icha] );
  }

  cout << myname << line << endl;
  cout << myname << "Done." << endl;
  return 0;
}

//**********************************************************************

int main(int argc, char* argv[]) {
  bool useExistingFcl = false;
  string extrafcl;
  if ( argc > 1 ) {
    string sarg(argv[1]);
    if ( sarg == "-h" ) {
      cout << "Usage: " << argv[0] << " [keepFCL] [extrafcl]" << endl;
      cout << "  If keepFCL = true, existing FCL file is used." << endl;
      cout << "  If extrafcl is give, that file is included in the config" << endl;
      return 0;
    }
    useExistingFcl = sarg == "true" || sarg == "1";
  }
  if ( argc > 2 ) {
    extrafcl = argv[2];
  }
  return test_ChannelStatusServiceTool(useExistingFcl, extrafcl);
}

//**********************************************************************
