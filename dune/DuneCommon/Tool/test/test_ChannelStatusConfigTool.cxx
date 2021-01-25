// test_ChannelStatusConfigTool.cxx
//
// David Adams
// April 2020
//
// Test ChannelStatusConfigTool.

#include "dune/DuneInterface/Data/AdcTypes.h"
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

int test_ChannelStatusConfigTool(bool useExistingFcl =false) {
  const string myname = "test_ChannelStatusConfigTool: ";
#ifdef NDEBUG
  cout << myname << "NDEBUG must be off." << endl;
  abort();
#endif
  string line = "-----------------------------";

  cout << myname << line << endl;
  string fclfile = "test_ChannelStatusConfigTool.fcl";
  if ( ! useExistingFcl ) {
    cout << myname << "Creating top-level FCL." << endl;
    ofstream fout(fclfile.c_str());
    fout << "tools: {" << endl;
    fout << "  mytool: {" << endl;
    fout << "    tool_type: ChannelStatusConfigTool" << endl;
    fout << "    LogLevel: 1" << endl;
    fout << "    DefaultIndex: 99" << endl;
    fout << "    IndexVectors: [" << endl;
    fout << "      [  0,  1,  2]," << endl;
    fout << "      [ 10, 11, 12]," << endl;
    fout << "      [ 20, 21, 22]" << endl;
    fout << "    ]" << endl;
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
  assert( tm.toolNames().size() >= 1 );

  cout << myname << line << endl;
  cout << myname << "Fetching tool." << endl;
  auto pcs = tm.getPrivate<IndexMapTool>("mytool");
  assert( pcs != nullptr );

  cout << myname << line << endl;
  cout << "Check some channels." << endl;
  vector<Index> chks(30, 99);
  chks[ 0] =  0;
  chks[ 1] =  0;
  chks[ 2] =  0;
  chks[10] =  1;
  chks[11] =  1;
  chks[12] =  1;
  chks[20] =  2;
  chks[21] =  2;
  chks[22] =  2;
  for ( Index icha=0; icha<chks.size(); ++icha ) {
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
  if ( argc > 1 ) {
    string sarg(argv[1]);
    if ( sarg == "-h" ) {
      cout << "Usage: " << argv[0] << " [keepFCL]" << endl;
      cout << "  If keepFCL = true, existing FCL file is used." << endl;
      return 0;
    }
    useExistingFcl = sarg == "true" || sarg == "1";
  }
  return test_ChannelStatusConfigTool(useExistingFcl);
}

//**********************************************************************
