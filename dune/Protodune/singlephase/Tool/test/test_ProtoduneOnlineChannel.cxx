// test_ProtoduneOnlineChannel.cxx
//
// David Adams
// May 2018
//
// Test ProtoduneOnlineChannel.

#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include "dune/ArtSupport/DuneToolManager.h"
#include "dune/DuneInterface/Tool/IndexMapTool.h"
#include "dune/Protodune/singlephase/Utility/ProtoduneChannelHelper.h"
#include "TH1F.h"

#undef NDEBUG
#include <cassert>

using std::string;
using std::cout;
using std::endl;
using std::ofstream;
using std::istringstream;
using std::setw;
using fhicl::ParameterSet;
using Index = IndexMapTool::Index;
using IndexVector = std::vector<Index>;

//**********************************************************************

int test_ProtoduneOnlineChannel(bool useExistingFcl =false, Index nshow =64) {
  const string myname = "test_ProtoduneOnlineChannel: ";
#ifdef NDEBUG
  cout << myname << "NDEBUG must be off." << endl;
  abort();
#endif
  string line = "-----------------------------";

  cout << myname << line << endl;
  string fclfile = "test_ProtoduneOnlineChannel.fcl";
  if ( ! useExistingFcl ) {
    cout << myname << "Creating top-level FCL." << endl;
    ofstream fout(fclfile.c_str());
    fout << "tools: {" << endl;
    fout << "  mytool: {" << endl;
    fout << "    tool_type: ProtoduneOnlineChannel" << endl;
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
  auto cma = tm.getPrivate<IndexMapTool>("mytool");
  assert( cma != nullptr );

  Index badIndex = IndexMapTool::badIndex();

  cout << myname << line << endl;
  cout << myname << "Check some good values." << endl;
  for ( Index ichaOff : { 0, 102, 1234, 1407, 3967, 15359 } ) {
    Index ichaOn = cma->get(ichaOff);
    cout << myname << ichaOff << " --> " << ichaOn << endl;
    assert( ichaOff != badIndex );
  }

  cout << myname << line << endl;
  cout << myname << "Check some bad values." << endl;
  for ( Index ichaOff : { -1, 15360, 20000 } ) {
    Index ichaOn = cma->get(ichaOff);
    cout << myname << ichaOff << " --> " << ichaOn << endl;
    assert( ichaOn == badIndex );
  }

  cout << myname << line << endl;
  cout << myname << "Check each online index appears exactly once." << endl;
  const Index ncha = 15360;
  IndexVector onlineCounts(ncha);
  IndexVector offlineChannel(ncha, badIndex);
  ProtoduneChannelHelper chh(false);
  for ( Index ichaOff=0; ichaOff<ncha; ++ichaOff ) {
    Index ichaOn = cma->get(ichaOff);
    Index irem = ichaOn;
    Index itps = irem/2560;
    irem = irem%2560;
    Index ifmb = irem/128 + 1;
    irem = irem%128;
    Index iasc = irem/16 + 1;
    Index iach = irem%16;
    if ( nshow*(ichaOff/nshow) == ichaOff || ichaOn >= ncha ) {
      cout <<  myname << " " << setw(4) << ichaOff << " --> " << setw(4) << ichaOn
           << " (" << itps << ", " << setw(2) << ifmb << ", "
           << iasc << ", " << setw(2) << iach << ")" << endl;
    }
    assert( itps == chh.tpcSet(ichaOn) );
    assert( ifmb == chh.femb(ichaOn) );
    assert( iasc == chh.asic(ichaOn) );
    assert( iach == chh.asicChannel(ichaOn) );
    assert( ichaOn < ncha );
    if ( offlineChannel[ichaOn] != badIndex ) {
      cout << myname << "ERROR: Online channel " << ichaOn
           << " is mapped to two offline channels:" << endl;
      cout << "  " << offlineChannel[ichaOn] << endl;
      cout << "  " << ichaOff << endl;
      assert( false );
    }
    assert( onlineCounts[ichaOn] == 0 );
    onlineCounts[ichaOn] += 1;
    offlineChannel[ichaOn] = ichaOff;
  }
  for ( Index ichaOn=0; ichaOn<ncha; ++ichaOn ) {
    assert( onlineCounts[ichaOn] == 1 );
    assert( offlineChannel[ichaOn] != badIndex );
  }

  cout << myname << line << endl;
  cout << myname << "Done." << endl;
  return 0;
}

//**********************************************************************

int main(int argc, char* argv[]) {
  bool useExistingFcl = false;
  Index nshow = 64;
  if ( argc > 1 ) {
    string sarg(argv[1]);
    if ( sarg == "-h" ) {
      cout << "Usage: " << argv[0] << " [keepFCL] [NSHOW]" << endl;
      cout << "  keepFCL [false]: If true, existing FCL file is used." << endl;
      cout << "  NSHOW [64]: Every nshow'th channels will be displayed in log." << endl;
      return 0;
    }
    useExistingFcl = sarg == "true" || sarg == "1";
  }
  if ( argc > 2 ) {
    string sarg(argv[2]);
    nshow = std::stoi(sarg);
  }
  return test_ProtoduneOnlineChannel(useExistingFcl, nshow);
}

//**********************************************************************
