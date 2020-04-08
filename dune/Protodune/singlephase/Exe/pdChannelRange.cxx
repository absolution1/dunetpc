// pdChannelRange.cxx
//
// David Adams
// October 2018
//
// Display the range of channels for a named channel range or group.

#include "dune/ArtSupport/DuneToolManager.h"
#include "dune/DuneInterface/Tool/IndexRangeTool.h"
#include "dune/DuneInterface/Tool/IndexRangeGroupTool.h"
#include <string>
#include <iostream>
#include <vector>
#include <sstream>

using std::string;
using std::cout;
using std::cin;
using std::endl;
using std::ostream;
using std::istringstream;

//**********************************************************************

int main(int argc, char** argv) {
  const string myname = "pdChannelRange: ";
  bool help = argc == 1;
  string crname;
  string tnameRange = "channelRanges";
  string tnameGroup = "channelGroups";
  int iarg = 1;
  string sdet = "protodune";
  while ( iarg < argc ) {
    string arg = argv[iarg++];
    if ( arg[0] == '-' ) {
      if ( arg == "-h" ) {
        help = true;
      } else if ( arg == "-d" ) {
        if ( iarg < argc ) sdet = argv[iarg++];
      } else {
        cout << "ERROR: Invalid flag: " << arg << endl;
        help = true;
      }
    } else {
      crname = arg;
    }
  }
  string fclfile = sdet + "_tools_dune.fcl";
  if ( help ) {
    cout << "Usage: " << argv[0] << " [-d DET] NAME" << endl;
    cout << "  Displays the range names group for channel group NAME" << endl;
    cout << "  or the range for channel range NAME." << endl;
    cout << "  DET = protodune (default) or apa7 or any value for which\n"
         << "        DET_tools_dune.fcl provides definitions for\n"
         << "        channelRanges and channelGroups." << endl;
    if ( sdet == "protodune" ) {
      cout << "  Protodune channel group names include:" << endl;
      cout << "    apas - APAs." << endl;
      cout << "    tpss - TPC sets." << endl;
      cout << "    apaus, apavs, apazs, apacs - APA planes with install naming." << endl;
      cout << "    tppus, tppvs, tppzs, tppcs - APA planes with offline naming." << endl;
      cout << "  Protodune channel range names include:" << endl;
      cout << "    apa1 - apa6: APA with APA numbering." << endl;
      cout << "    tps0 - tps5: APA (TPC set) with offline numbering." << endl;
      cout << "    apa0u, apa0v, apa0z, apa0c, ..., apa5c: APA plane with install numbering." << endl;
      cout << "    tpp0u, tpp0v, tpp0z, tpp0c, ..., tpp5c: APA plane with offline numbering." << endl;
      cout << "    (u and v are induction, z is TPC-side collection, c is cryostat side." << endl;
      cout << "    femb101u, femb102u, ..., femb620u: FEMB u channels (same for v and x)." << endl;
    } else if ( sdet == "apa7" ) {
      cout << "  CERN coldbox channel range names include:" << endl;
      cout << "    apa7: same as all." << endl;
      cout << "    apa7u, apa7v, apa7z, apa7c: APA plane." << endl;
      cout << "    tpp7u, tpp7v, tpp7z, tpp07: APA plane." << endl;
      cout << "    (u and v are induction, z is TPC-side collection, c is cryostat side." << endl;
      cout << "    femb701u, femb702u, ..., femb720u: FEMB u channels (same for v and x)." << endl;
    } else if ( sdet == "iceberg" ) {
      cout << "  Iceberg channel range names include:" << endl;
      cout << "    apa: same as all." << endl;
      cout << "    apau, apav, apaz1, apaz2: APA plane." << endl;
      cout << "    tppu, tppv, tppz1, apaz2: APA plane." << endl;
      cout << "    (u and v are induction, z is TPC-side collection, c is cryostat side." << endl;
      cout << "    femb701u, femb702u, ..., femb720u: FEMB u channels (same for v and x)." << endl;
      cout << "    apai: Induction channels." << endl;
      cout << "    apaz or apax: Collection channels." << endl;
    }
    return 0;
  }

  DuneToolManager* ptm = DuneToolManager::instance(fclfile, 0);
  assert ( ptm != nullptr );
  DuneToolManager& tm = *ptm;
  //tm.print();
  assert( tm.toolNames().size() >= 1 );

  IndexRangeGroupTool* pcgt = tm.getShared<IndexRangeGroupTool>(tnameGroup);
  IndexRangeTool* pcrt = tm.getShared<IndexRangeTool>(tnameRange);

  int warnGroup = 0;
  int warnRange = 0;
  // If user provides name, exit after displaying result.
  bool exitLoop = crname.size();
  while ( true ) {
    if ( crname.size() == 0 ) {
      cout << "Group or range name> ";
      cout.flush();
      std::getline(cin, crname);
    }
    if ( crname.size() == 0  || crname == "q" ) return 0;
    IndexRangeGroup cg;
    if ( pcgt == nullptr ) {
      if ( ! warnGroup++ ) cout << "WARNING: Channel group tool not found: " << tnameGroup << endl;
    } else {
      cg = pcgt->get(crname);
    }
    if ( cg.isValid() ) {
      cout << cg << endl;
    } else if ( pcrt == nullptr ) {
      if ( ! warnRange ) cout << "WARNING: Channel range tool not found: " << tnameRange << endl;
    } else {
      IndexRange cr = pcrt->get(crname);
      if ( cr.isValid() ) {
        cout << cr << endl;
      } else {
        cout << "ERROR: " << crname << " is not a valid channel group or range." << endl;
      }
    }
    if ( warnGroup && warnRange ) return 1;
    if ( exitLoop ) return 0;
    crname = "";
  }
  return 0;
}
