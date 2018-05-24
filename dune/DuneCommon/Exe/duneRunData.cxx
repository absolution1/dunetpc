// duneRunData.cxx
//
// David Adams
// May 2018
//
// Executable that displays the run data for a run.

#include <string>
#include <iostream>
#include <vector>
#include <sstream>
#include <memory>
#include "dune/ArtSupport/DuneToolManager.h"
#include "dune/DuneInterface/Tool/RunDataTool.h"

using std::string;
using std::cout;
using std::endl;
using std::istringstream;

int main(int argc, char** argv) {
  const string myname = "duneRunData: ";
  bool help = false;
  unsigned int irun = 0;
  bool verbose = false;
  string fname;
  if ( argc > 1 ) {
    for ( int iarg=1; iarg< argc; ++iarg ) {
      string sarg = argv[iarg];
      if ( sarg == "-h" ) {
        help = true;
        break;
      }
      if ( sarg == "-v" ) verbose = true;
    }
    if ( ! help ) {
      string sarg = argv[argc-1];
      istringstream ssarg(sarg);
      ssarg >> irun;
      if ( irun == 0 ) {
        cout << myname << "ERROR: Last field does not appear to be a number: \""
             << sarg << "\"" << endl;
        return 1;
      }
    }
  } else {
    help = true;
  } 
  if ( help ) {
    cout << "Usage: duneRunData [-h] [-v] [-c FCLFILE] RUN" << endl;
    cout << "  Displays the run info for a run." << endl;
    cout << "  -h - Display this message." << endl;
    cout << "  -v - Verbose display including info about fcl file and run data tool." << endl;
    return 0;
  }
  string tname = "protoduneRunDataTool";
  if ( ! verbose ) cout.setstate(std::ios_base::failbit);   // Turn off cout
  cout << endl;
  cout << myname << "Fetching tool manager." << endl;
  DuneToolManager* ptm = DuneToolManager::instance();
  cout << endl;
  cout << myname << "Fetching run data tool." << endl;
  const RunDataTool* prdt = ptm->getShared<RunDataTool>(tname);
  cout << endl;
  if ( ! verbose ) cout.clear();  // Turn cout back on.
  if ( prdt == nullptr ) {
    cout << myname << "ERROR: Unable to find RunDataTool " << tname << endl;
    return 1;
  }
  RunData rdat = prdt->runData(irun);
  if ( rdat.isValid() ) {
    cout << rdat << endl;
  } else {
    cout << myname << "Unable to find run " << irun << endl;
    return 2;
  }
  return 0;
}
