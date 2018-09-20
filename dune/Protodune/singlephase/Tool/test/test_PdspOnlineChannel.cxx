// test_PdspOnlineChannel.cxx
//
// David Adams
// May 2018
//
// Test PdspOnlineChannel.

#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "dune/ArtSupport/DuneToolManager.h"
#include "dune/DuneInterface/Tool/IndexMapTool.h"
#include "dune/ArtSupport/ArtServiceHelper.h"
#include "TH1F.h"

#undef NDEBUG
#include <cassert>

using std::string;
using std::cout;
using std::endl;
using std::ofstream;
using std::istringstream;
using fhicl::ParameterSet;
using Index = IndexMapTool::Index;
using IndexVector = std::vector<Index>;

//**********************************************************************

int test_PdspOnlineChannel(bool useExistingFcl =false) {
  const string myname = "test_PdspOnlineChannel: ";
#ifdef NDEBUG
  cout << myname << "NDEBUG must be off." << endl;
  abort();
#endif
  string line = "-----------------------------";

  cout << myname << line << endl;
  string fclfile = "test_PdspOnlineChannel.fcl";
  if ( ! useExistingFcl ) {
    cout << myname << "Creating top-level FCL." << endl;
    ofstream fout(fclfile.c_str());
    fout << "#include \"PdspChannelMapService.fcl\"" << endl;
    fout << "services: { PdspChannelMapService: @local::pdspchannelmap }" << endl;
    fout << "tools: {" << endl;
    fout << "  mytool: {" << endl;
    fout << "    tool_type: PdspOnlineChannel" << endl;
    fout << "     LogLevel: 1" << endl;
    fout << "     Ordering: \"FEMB\"" << endl;
    fout << "  }" << endl;
    fout << "  reftool: {" << endl;
    fout << "    tool_type: ProtoduneOnlineChannel" << endl;
    fout << "     LogLevel: 1" << endl;
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

  // Load services with the service helper.
  cout << myname << line << endl;
  cout << myname << "Fetch art service helper." << endl;
  ArtServiceHelper& ash = ArtServiceHelper::instance();
  assert( ash.addServices(fclfile, true) == 0 );
  assert( ash.loadServices() == 1 );
  ash.print();

  cout << myname << line << endl;
  cout << myname << "Fetching tool." << endl;
  auto cma = tm.getPrivate<IndexMapTool>("mytool");
  assert( cma != nullptr );

  Index badIndex = IndexMapTool::badIndex();

  cout << myname << line << endl;
  cout << myname << "Check some good values." << endl;
  for ( Index ichaOff : { 0, 102, 1234, 15359 } ) {
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
  Index nshow = 64;
  for ( Index ichaOff=0; ichaOff<ncha; ++ichaOff ) {
    Index ichaOn = cma->get(ichaOff);
    if ( nshow*(ichaOff/nshow) == ichaOff || ichaOn >= ncha )
      cout <<  myname << "  "  << ichaOff << " --> " << ichaOn << endl;
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
  cout << myname << "Compare with ProtoduneChannelmap." << endl;
  auto ref = tm.getPrivate<IndexMapTool>("reftool");
  assert( ref != nullptr );
  bool skipDiv1 = false;  // Set this false to check every channel.
  for ( Index idiv : {2560, 128, 1} ) {
    if ( skipDiv1 && idiv == 1 ) {
      cout << myname << "WARNING: Skipping div 1 test" << endl;
      continue;
    }
    cout << myname << "...checking div " << idiv << endl;
    for ( Index ichaOff=0; ichaOff<ncha; ++ichaOff ) {
      Index ichaOnl = cma->get(ichaOff);
      Index ichaRef = ref->get(ichaOff);
      Index ichaOnlDiv = ichaOnl/idiv;
      Index ichaRefDiv = ichaRef/idiv;
      if ( ichaOnlDiv != ichaRefDiv ) {
        cout << myname << "Maps disagree:" << endl;
        cout << myname << "  Offline: " << ichaOff << endl;
        cout << myname << "   Online: " << ichaOnl << " ("  << ichaOnlDiv << ")" << endl;
        cout << myname << "      Ref: " << ichaRef << " ("  << ichaRefDiv << ")" << endl;
        assert(false);
      }
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
      cout << "Usage: " << argv[0] << " [keepFCL] [RUN]" << endl;
      cout << "  If keepFCL = true, existing FCL file is used." << endl;
      cout << "  If RUN is nonzero, the data for that run are displayed." << endl;
      return 0;
    }
    useExistingFcl = sarg == "true" || sarg == "1";
  }
  return test_PdspOnlineChannel(useExistingFcl);
}

//**********************************************************************
