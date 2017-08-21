// test_SimpleHistogramManager.cxx
//
// David Adams
// April 2017
//
// Test SimpleHistogramManager.

#include <string>
#include <iostream>
#include <fstream>
#include "dune/DuneInterface/Tool/HistogramManager.h"
#include "dune/ArtSupport/DuneToolManager.h"
#include "TH1F.h"

#undef NDEBUG
#include <cassert>

using std::string;
using std::cout;
using std::endl;
using std::ofstream;
using fhicl::ParameterSet;

//**********************************************************************

int test_SimpleHistogramManager(bool useExistingFcl =false) {
  const string myname = "test_SimpleHistogramManager: ";
#ifdef NDEBUG
  cout << myname << "NDEBUG must be off." << endl;
  abort();
#endif
  string line = "-----------------------------";

  cout << myname << line << endl;
  string fclfile = "test_SimpleHistogramManager.fcl";
  if ( ! useExistingFcl ) {
    cout << myname << "Creating top-level FCL." << endl;
    ofstream fout(fclfile.c_str());
    fout << "tools: {" << endl;
    fout << "  mytool: {" << endl;
    fout << "    tool_type: SimpleHistogramManager" << endl;
    fout << "    LogLevel: 1" << endl;
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
  assert( tm.toolNames().size() == 1 );

  cout << myname << line << endl;
  cout << myname << "Fetching tool." << endl;
  auto phm = tm.getPrivate<HistogramManager>("mytool");
  assert( phm != nullptr );

  cout << myname << line << endl;
  cout << myname << "Creating histogram." << endl;
  TH1* ph = new TH1F("h1", "h1", 10, 0, 10);
  assert(ph->GetDirectory() != nullptr);

  cout << myname << line << endl;
  cout << myname << "Managing histogram." << endl;
  assert(phm->manage(ph) == 0);
  assert(ph->GetDirectory() == nullptr);

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
  return test_SimpleHistogramManager(useExistingFcl);
}

//**********************************************************************
