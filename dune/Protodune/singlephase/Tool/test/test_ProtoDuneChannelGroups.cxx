// test_ProtoDuneChannelGroups.cxx
//
// David Adams
// July 2018
//
// Test ProtoDuneChannelGroups.

#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include "dune/ArtSupport/DuneToolManager.h"
#include "dune/DuneInterface/Tool/IndexRangeGroupTool.h"
#include "TH1F.h"

#undef NDEBUG
#include <cassert>

using std::string;
using std::cout;
using std::endl;
using std::ofstream;
using std::istringstream;
using std::ostringstream;
using std::setw;
using std::setfill;
using std::vector;
using fhicl::ParameterSet;
using Index = unsigned int;
using IndexVector = std::vector<Index>;
using Name = std::string;

//**********************************************************************

int test_ProtoDuneChannelGroups(bool useExistingFcl =false, int show =1) {
  const string myname = "test_ProtoDuneChannelGroups: ";
#ifdef NDEBUG
  cout << myname << "NDEBUG must be off." << endl;
  abort();
#endif
  string line = "-----------------------------";

  cout << myname << line << endl;
  string fclfile = "test_ProtoDuneChannelGroups.fcl";
  if ( ! useExistingFcl ) {
    cout << myname << "Creating top-level FCL." << endl;
    ofstream fout(fclfile.c_str());
    fout << "tools: {" << endl;
    fout << "  crtool: {" << endl;
    fout << "    tool_type: ProtoDuneChannelRanges" << endl;
    fout << "    LogLevel: 2" << endl;
    fout << "    ExtraRanges: \"\"" << endl;
    fout << "  }" << endl;
    fout << "  mytool: {" << endl;
    fout << "    tool_type: ProtoDuneChannelGroups" << endl;
    fout << "    LogLevel: 2" << endl;
    fout << "    IndexRangeTool: crtool" << endl;
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

  cout << myname << line << endl;
  cout << myname << "Fetching tool." << endl;
  auto pcgt = tm.getPrivate<IndexRangeGroupTool>("mytool");
  assert( pcgt != nullptr );

  cout << myname << line << endl;
  cout << myname << "Fetch invalid." << endl;
  IndexRangeGroup chg = pcgt->get("bad");
  cout << chg << endl;
  assert( ! chg.isValid() );

  vector<string> expApas = {"apa1", "apa2", "apa3", "apa4", "apa5", "apa6"};
  vector<string> expTpss = {"tps0", "tps1", "tps2", "tps3", "tps4", "tps5"};

  cout << myname << line << endl;
  cout << myname << "Check apas." << endl;
  chg = pcgt->get("apas");
  cout << myname << chg << endl;
  assert( chg.isValid() );
  assert( chg.name == "apas" );
  Index nran = chg.ranges.size();
  assert( nran == expApas.size() );
  for ( Index iran=0; iran<nran; ++iran ) {
    IndexRange ran = chg.range(iran);
    cout << myname << setw(5) << iran << ": " << ran << endl;
    assert( ran.name == expApas[iran] );
  }

  cout << myname << line << endl;
  cout << myname << "Check tpss." << endl;
  chg = pcgt->get("tpss");
  cout << myname << chg << endl;
  assert( chg.isValid() );
  assert( chg.name == "tpss" );
  nran = chg.ranges.size();
  assert( nran == expTpss.size() );
  for ( Index iran=0; iran<nran; ++iran ) {
    IndexRange ran = chg.range(iran);
    cout << myname << setw(5) << iran << ": " << ran << endl;
    assert( ran.name == expTpss[iran] );
  }

  // Check plane groups: tpszs, apaus, ...
  for ( Name sapt : {"apa", "tps"} ) {
    Name rbas = sapt;
    if ( sapt == "tps" ) rbas = "tpp";
    for ( Name sori : {"z", "c", "u", "v"} ) {
      vector<string> rnams;
      if ( sapt == "apa" ) {
        for ( Name sapa : { "1", "2", "3", "4", "5", "6" } ) rnams.push_back(rbas + sapa + sori);
      } else {
        for ( Name stps : { "0", "1", "2", "3", "4", "5" } ) rnams.push_back(rbas + stps + sori);
      }
      Name gnam = rbas + sori + "s";
      cout << myname << line << endl;
      cout << myname << "Check " << gnam << "." << endl;
      chg = pcgt->get(gnam);
      cout << myname << chg << endl;
      assert( chg.isValid() );
      assert( chg.name == gnam );
      Index nran = chg.ranges.size();
      assert( nran == rnams.size() );
      for ( Index iran=0; iran<nran; ++iran ) {
        IndexRange ran = chg.range(iran);
        cout << myname << setw(5) << iran << ": " << ran << endl;
        assert( ran.name == rnams[iran] );
      }
    }
  }

  // Check FEMBs.
  cout << myname << "Checking FEMBs" << endl;
  for ( Name sapa : { "1", "2", "3", "4", "5", "6" } ) {
    for ( Name sfmb : { "01", "02", "03", "04", "05", "06", "07", "08", "09", "10",
                        "11", "12", "13", "14", "15", "16", "17", "18", "19", "20" } ) {
      Name gnam = "femb" + sapa + sfmb;
      chg = pcgt->get(gnam);
      assert( chg.isValid() );
      assert( chg.name == gnam );
      Index nran = chg.ranges.size();
      assert( nran == 3 );
      cout << myname << "  " << gnam << ": ";
      for ( Index iran=0; iran<nran; ++iran ) {
        if ( iran ) cout << ",";
        cout << " " << chg.ranges[iran].name;
      }
      cout << endl;
    }
  }
  
  cout << myname << line << endl;
  cout << myname << "Done." << endl;
  return 0;
}

//**********************************************************************

int main(int argc, char* argv[]) {
  bool useExistingFcl = false;
  int show = 1;
  if ( argc > 1 ) {
    string sarg(argv[1]);
    if ( sarg == "-h" ) {
      cout << "Usage: " << argv[0] << " [keepFCL] [SHOW]" << endl;
      cout << "  If keepFCL = true, existing FCL file is used." << endl;
      cout << "  SHOW > 1 also shows the 480 FEMB blocks." << endl;
      return 0;
    }
    useExistingFcl = sarg == "true" || sarg == "1";
  }
  if ( argc > 2 ) {
    string sarg(argv[2]);
    show = std::stoi(sarg);
  }
  return test_ProtoDuneChannelGroups(useExistingFcl, show);
}

//**********************************************************************
