// test_DuneToolManager.cxx

//#include "AXService/DuneToolManager.h"
#include "dune/ArtSupport/DuneToolManager.h"

#include <string>
#include <iostream>

using std::string;
using std::cout;
using std::endl;

#undef NDEBUG
#include <cassert>

int test_DuneToolManager() {
  const string myname = "test_DuneToolManager: ";
  cout << myname << "Starting test" << endl;
#ifdef NDEBUG
  cout << myname << "NDEBUG must be off." << endl;
  abort();
#endif
  string line = "-----------------------------";
  string scfg;

  cout << line << endl;
  cout << myname << "Retrieve tool manager." << endl;
  DuneToolManager* pdtm = DuneToolManager::instance();
  assert( pdtm != nullptr );
  DuneToolManager& dtm = *pdtm;

  cout << line << endl;
  cout << myname << "Display tools." << endl;
  dtm.print();

  cout << line << endl;
  cout << myname << "Check tool list." << endl;
  assert( dtm.toolNames().size() > 0 );

  cout << myname << line << endl;
  cout << myname << "Done." << endl;
  return 0;
}

int main() {
  int rstat = test_DuneToolManager();
  return rstat;
}
