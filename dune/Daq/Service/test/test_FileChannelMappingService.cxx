// test_FileChannelMappingService.cxx

// David Adams
// February 2015
//
// Test FileChannelMappingService.

#include "../FileChannelMappingService.h"
#include <string>
#include <iostream>
#include <sstream>
#include <fstream>
#include <iomanip>
//#include "dune/ArtSupport/ArtServiceHelper.h"

using std::string;
using std::cout;
using std::endl;
using std::setw;
typedef ChannelMappingService::Channel Channel;

#undef NDEBUG
#include <cassert>

int test_FileChannelMappingService(string fname) {
  const string myname = "test_FileChannelMappingService: ";
#ifdef NDEBUG
  cout << myname << "NDEBUG must be off." << endl;
  abort();
#endif
  const string line = "-----------------------------";

  cout << myname << line << endl;
  cout << myname << "Create channel mapping service." << endl;
  fhicl::ParameterSet pset;
  pset.put("FileName", fname);
  FileChannelMappingService mapsvc(pset);

  cout << myname << line << endl;
  cout << myname << "Fetch offline channels." << endl;
  Channel maxch = 20;
  unsigned int w = 12;
  cout << myname << setw(w) << "online" << setw(w) << "offline" << endl;
  for ( unsigned int chin=0; chin<maxch; ++chin ) {
    cout << myname << setw(w) << chin << setw(w) << mapsvc.offline(chin) << endl;
  }

  cout << myname << line << endl;
  cout << myname << "Fetch online channels." << endl;
  cout << myname << setw(w) << "offline" << setw(w) << "online" << endl;
  for ( unsigned int chin=0; chin<maxch; ++chin ) {
    cout << myname << setw(w) << chin << setw(w) << mapsvc.online(chin) << endl;
  }

  cout << myname << line << endl;
  cout << myname << "Check mapping failures raise exceptions." << endl;
  Channel nch = 2048;
  try {
    mapsvc.offline(nch);
    cout << myname << "ERROR: Invalid online channel did not raise exception!" << endl;
    return 1;
  } catch(...) {
    cout << myname << "Invalid online channel raised exception." << endl;
  }
  try {
    mapsvc.online(nch);
    cout << myname << "ERROR: Invalid offline channel did not raise exception!" << endl;
    return 1;
  } catch(...) {
    cout << myname << "Invalid offline channel raised exception." << endl;
  }

  cout << myname << line << endl;
  cout << myname << "Check map consistency." << endl;
  for ( Channel chin=0; chin<nch; ++chin ) {
    Channel choff = mapsvc.offline(chin);
    Channel chon = mapsvc.online(choff);
    assert( chon == chin );
  }

  cout << myname << line << endl;
  cout << myname << "Done." << endl;
  return 0;
}

int main(int argc, char* argv[]) {
  test_FileChannelMappingService("rce_channel_map_dune35t.txt");
  return 0;
}
