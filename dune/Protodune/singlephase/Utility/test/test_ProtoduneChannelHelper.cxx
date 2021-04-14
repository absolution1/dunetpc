// test_ProtoduneChannelHelper.cxx
//
// David Adams
// May 2018
//
// Test ProtoduneChannelHelper.

#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
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
using Index = ProtoduneChannelHelper::Index;
using IndexVector = std::vector<Index>;

//**********************************************************************

int test_ProtoduneChannelHelper(Index nshow =64) {
  const string myname = "test_ProtoduneChannelHelper: ";
#ifdef NDEBUG
  cout << myname << "NDEBUG must be off." << endl;
  abort();
#endif
  string line = "-----------------------------";

  cout << myname << line << endl;

  ProtoduneChannelHelper chh(false);
  Index badIndex = chh.badIndex();

  cout << myname << "Check some good values." << endl;
  for ( Index ichaOff : { 0, 102, 1234, 1407, 3967, 15359 } ) {
    Index ichaOn = chh.onlineChannel(ichaOff);
    cout << myname << ichaOff << " --> " << ichaOn << endl;
    assert( ichaOff != badIndex );
  }

  cout << myname << line << endl;
  cout << myname << "Check some bad values." << endl;
  for ( Index ichaOff : { -1, 15360, 20000 } ) {
    Index ichaOn = chh.onlineChannel(ichaOff);
    cout << myname << ichaOff << " --> " << ichaOn << endl;
    assert( ichaOn == badIndex );
  }

  cout << myname << line << endl;
  cout << myname << "Check each online index appears exactly once." << endl;
  const Index ncha = 15360;
  IndexVector onlineCounts(ncha);
  IndexVector offlineChannel(ncha, badIndex);
  for ( Index ichaOff=0; ichaOff<ncha; ++ichaOff ) {
    Index ichaOn = chh.onlineChannel(ichaOff);
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
  Index nshow = 64;
  if ( argc > 1 ) {
    string sarg(argv[1]);
    if ( sarg == "-h" ) {
      cout << "Usage: " << argv[0] << "[NSHOW]" << endl;
      cout << "  NSHOW [64]: Every nshow'th channels will be displayed in log." << endl;
      return 0;
    }
    nshow = std::stoi(sarg);
  }
  return test_ProtoduneChannelHelper(nshow);
}

//**********************************************************************
