// test_DuneToolManager.cxx

//#include "AXService/DuneToolManager.h"
#include "dune/ArtSupport/DuneToolManager.h"
#include "dune/ArtSupport/Tool/TestTool.h"

#include <string>
#include <iostream>
#include <fstream>

using std::string;
using std::cout;
using std::endl;
using std::ofstream;

#undef NDEBUG
#include <cassert>

int test_DuneToolManager(bool useExistingFcl =false) {
  const string myname = "test_DuneToolManager: ";
  cout << myname << "Starting test" << endl;
#ifdef NDEBUG
  cout << myname << "NDEBUG must be off." << endl;
  abort();
#endif
  string line = "-----------------------------";
  string scfg;

  cout << myname << line << endl;
  string fclfile = "test_DuneToolManager.fcl";
  if ( ! useExistingFcl ) {
    cout << myname << "Creating top-level FCL." << endl;
    ofstream fout(fclfile.c_str());
    fout << "tools: {" << endl;
    fout << "  mytool1: {" << endl;
    fout << "    tool_type: TestTool" << endl;
    fout << "    Label: \"Tool 1\"" << endl;
    fout << "  }" << endl;
    fout << "  mytool2: {" << endl;
    fout << "    tool_type: TestTool" << endl;
    fout << "    Label: \"Tool 2\"" << endl;
    fout << "  }" << endl;
    fout << "}" << endl;
    fout.close();
  } else {
    cout << myname << "Using existing top-level FCL." << endl;
  }

  cout << line << endl;
  cout << myname << "Retrieve tool manager." << endl;
  DuneToolManager* ptm = DuneToolManager::instance(fclfile);
  assert( ptm != nullptr );
  DuneToolManager& dtm = *ptm;

  cout << myname << line << endl;
  cout << myname << "Display tools." << endl;
  dtm.print();

  cout << myname << line << endl;
  cout << myname << "Check tool list." << endl;
  assert( dtm.toolNames().size() > 0 );

  cout << myname << line << endl;
  cout << myname << "Fetch a private tool." << endl;
  auto ptp1 = ptm->getPrivate<TestTool>("mytool1");
  assert( ptp1 != nullptr );
  assert( ptp1->label() == "Tool 1" );

  cout << myname << line << endl;
  cout << myname << "Fetch a private tool." << endl;
  auto ptp2 = ptm->getPrivate<TestTool>("mytool1");
  assert( ptp2 != nullptr );
  assert( ptp2 != ptp1 );
  assert( ptp2->label() == "Tool 1" );

  cout << myname << line << endl;
  cout << myname << "Fetch a shared tool." << endl;
  auto pts1 = ptm->getShared<TestTool>("mytool1");
  assert( pts1 != nullptr );
  assert( pts1 != ptp1.get() );
  assert( pts1->label() == "Tool 1" );

  cout << myname << line << endl;
  cout << myname << "Fetch shared tool again." << endl;
  auto pts2 = ptm->getShared<TestTool>("mytool1");
  assert( pts2 != nullptr );
  assert( pts2 == pts1 );
  assert( pts2->label() == "Tool 1" );
  
  cout << myname << line << endl;
  cout << myname << "Fetch another shared tool." << endl;
  auto pts3 = ptm->getShared<TestTool>("mytool2");
  assert( pts3 != nullptr );
  assert( pts3 != pts2 );
  assert( pts3->label() == "Tool 2" );

  cout << myname << line << endl;
  cout << myname << "Fetch a private tool with config." << endl;
  string scfgTool = "{ tool_type:TestTool Label:\"Tool C\" }";
  cout << myname << "Config: " << scfgTool << endl;
  auto ptpc = ptm->getPrivate<TestTool>(scfgTool);
  assert( ptpc != nullptr );
  assert( ptpc != ptp1 );
  assert( ptpc->label() == "Tool C" );

  cout << myname << line << endl;
  cout << myname << "Done." << endl;
  return 0;
}

int main() {
  int rstat = test_DuneToolManager();
  return rstat;
}
