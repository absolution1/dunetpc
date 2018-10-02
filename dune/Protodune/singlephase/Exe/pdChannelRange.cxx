// pdChannelRange.cxx
//
// David Adams
// October 2018
//
// Display the range of channels for a named channel range.

#include "dune/ArtSupport/DuneToolManager.h"
#include "dune/DuneInterface/Tool/IndexRangeTool.h"
#include <string>
#include <iostream>
#include <vector>
#include <sstream>

using std::string;
using std::cout;
using std::endl;
using std::ostream;
using std::istringstream;

//**********************************************************************

int main(int argc, char** argv) {
  const string myname = "pdChannelRange: ";
  bool help = argc == 1;
  string crname;
  string fclfile = "protodune_tools_dune.fcl";
  string tname = "channelRanges";
  if ( argc > 1 ) {
    string arg = argv[1];
    if ( arg == "-h" ) help = true;
    else crname = argv[1];
  }
  if ( help ) {
    cout << "Usage: " << argv[0] << " NAME" << endl;
    cout << "  Displays the range for channel range NAME." << endl;
    cout << "  Channel range names include:" << endl;
    cout << "    apa1 - apa6: APA with installation numbering." << endl;
    cout << "    tps0 - tps5: APA (TPC set) with offline numbering." << endl;
    cout << "    tpp0u, tpp0v, tpp0z, tpp0c, ..., tpp5c: APA plane with offline numbering." << endl;
    cout << "    (u and v are induction, z is TPC-side collection, c is cryostat side." << endl;
    return 0;
  }

  DuneToolManager* ptm = DuneToolManager::instance(fclfile, 0);
  assert ( ptm != nullptr );
  DuneToolManager& tm = *ptm;
  //tm.print();
  assert( tm.toolNames().size() >= 1 );

  IndexRangeTool* pcrt = tm.getShared<IndexRangeTool>(tname);
  IndexRange cr = pcrt->get(crname);
  cout << cr << endl;

  return 0;
}
