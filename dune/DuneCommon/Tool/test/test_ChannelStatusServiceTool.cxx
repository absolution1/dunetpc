// test_ChannelStatusServiceTool.cxx
//
// David Adams
// Novenber 2019
//
// Test ChannelStatusServiceTool.

#include "dune/DuneInterface/AdcTypes.h"
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
using fhicl::ParameterSet;
using Index = unsigned int;

//**********************************************************************

int test_ChannelStatusServiceTool(bool useExistingFcl =false, bool useDuneServices = false) {
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
    fout << "tools: {" << endl;
    fout << "  mytool: {" << endl;
    fout << "    tool_type: ChannelStatusServiceTool" << endl;
    fout << "    LogLevel: 2" << endl;
    fout << "  }" << endl;
    fout << "}" << endl;
    if ( useDuneServices ) {
      fout << "#include \"dune_services.fcl\"" << endl;
    }
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
  cout << myname << "Fetching tool." << endl;
  auto pcs = tm.getPrivate<IndexMapTool>("mytool");
  assert( pcs != nullptr );

  cout << myname << line << endl;
  cout << "Check some channels." << endl;
  for ( Index icha=0; icha<20; ++icha ) {
    Index chanStat = pcs->get(icha);
    cout << setw(8) << icha << ": " << setw(2) << chanStat << endl;
    assert( chanStat == AdcChannelStatusUnknown );
  }

  cout << myname << line << endl;
  cout << myname << "Done." << endl;
  return 0;
}

//**********************************************************************

int main(int argc, char* argv[]) {
  bool useExistingFcl = false;
  bool useDuneServices = false;
  Index run = 0;
  if ( argc > 1 ) {
    string sarg(argv[1]);
    if ( sarg == "-h" ) {
      cout << "Usage: " << argv[0] << " [keepFCL] [useDuneServices]" << endl;
      cout << "  If keepFCL = true, existing FCL file is used." << endl;
      cout << "  If useDuneServices, include dune services." << endl;
      return 0;
    }
    useExistingFcl = sarg == "true" || sarg == "1";
  }
  if ( argc > 2 ) {
    string sarg(argv[2]);
    istringstream ssarg(argv[2]);
    ssarg >> run;
  }
  return test_ChannelStatusServiceTool(useExistingFcl, useDuneServices);
}

//**********************************************************************
