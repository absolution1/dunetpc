// pdChannelRange.cxx
//
// David Adams
// October 2018
//
// Display the range of channels for a named channel range.

#include "dune/ArtSupport/DuneToolManager.h"
#include "dune/DuneInterface/Tool/IndexRangeTool.h"
#include "dune/DuneInterface/Tool/IndexRangeGroupTool.h"
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
  string tnameRange = "channelRanges";
  string tnameGroup = "channelGroups";
  if ( argc > 1 ) {
    string arg = argv[1];
    if ( arg == "-h" ) help = true;
    else crname = argv[1];
  }
  if ( help ) {
    cout << "Usage: " << argv[0] << " NAME" << endl;
    cout << "  Displays the group for channel group NAME" << endl;
    cout << "  or the range for channel range NAME." << endl;
    cout << "  Channel group names include:" << endl;
    cout << "    apas - APAs." << endl;
    cout << "    tpss - TPC sets." << endl;
    cout << "    apaus, apavs, apazs, apacs - APA planes with install naming." << endl;
    cout << "    tppus, tppvs, tppzs, tppcs - APA planes with offline naming." << endl;
    cout << "  Channel range names include:" << endl;
    cout << "    apa1 - apa6: APA with install numbering." << endl;
    cout << "    tps0 - tps5: APA (TPC set) with offline numbering." << endl;
    cout << "    apa0u, apa0v, apa0z, apa0c, ..., apa5c: APA plane with install numbering." << endl;
    cout << "    tpp0u, tpp0v, tpp0z, tpp0c, ..., tpp5c: APA plane with offline numbering." << endl;
    cout << "    (u and v are induction, z is TPC-side collection, c is cryostat side." << endl;
    cout << "    femb101u, femb102u, ..., femb620u: FEMB u channels (same for v and x)." << endl;
    return 0;
  }

  DuneToolManager* ptm = DuneToolManager::instance(fclfile, 0);
  assert ( ptm != nullptr );
  DuneToolManager& tm = *ptm;
  //tm.print();
  assert( tm.toolNames().size() >= 1 );

  IndexRangeGroupTool* pcgt = tm.getShared<IndexRangeGroupTool>(tnameGroup);
  if ( pcgt == nullptr ) {
    cout << "WARNING: Channel group tool not found: " << tnameGroup << endl;
  } else {
    IndexRangeGroup cg = pcgt->get(crname);
    if ( cg.isValid() ) {
      cout << cg << endl;
      return 0;
    }
  }
  IndexRangeTool* pcrt = tm.getShared<IndexRangeTool>(tnameRange);
  if ( pcrt == nullptr ) {
    cout << "WARNING: Channel range tool not found: " << tnameRange << endl;
  } else {
    IndexRange cr = pcrt->get(crname);
    if ( cr.isValid() ) {
      cout << cr << endl;
      return 0;
    }
  }
  cout << crname << " is not a valid channel group or range." << endl;
  return 0;
}
