// test_ProtoDuneChannelRanges.cxx
//
// David Adams
// July 2018
//
// Test ProtoDuneChannelRanges.

#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include "dune/ArtSupport/DuneToolManager.h"
#include "dune/DuneInterface/Tool/IndexRangeTool.h"
#include "TH1F.h"

#undef NDEBUG
#include <cassert>

using std::string;
using std::cout;
using std::endl;
using std::ofstream;
using std::istringstream;
using std::setw;
using std::vector;
using fhicl::ParameterSet;
using Index = unsigned int;
using IndexVector = std::vector<Index>;

//**********************************************************************

int test_ProtoDuneChannelRanges(bool useExistingFcl =false) {
  const string myname = "test_ProtoDuneChannelRanges: ";
#ifdef NDEBUG
  cout << myname << "NDEBUG must be off." << endl;
  abort();
#endif
  string line = "-----------------------------";

  cout << myname << line << endl;
  string fclfile = "test_ProtoDuneChannelRanges.fcl";
  if ( ! useExistingFcl ) {
    cout << myname << "Creating top-level FCL." << endl;
    ofstream fout(fclfile.c_str());
    fout << "tools: {" << endl;
    fout << "  mytool: {" << endl;
    fout << "    tool_type: ProtoDuneChannelRanges" << endl;
    fout << "    LogLevel: 2" << endl;
    fout << "    ExtraRanges: \"\"" << endl;
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
  auto irt = tm.getPrivate<IndexRangeTool>("mytool");
  assert( irt != nullptr );

  cout << myname << line << endl;
  cout << myname << "Fetching TPC ranges." << endl;
  int nbad = 0;
  vector<string> namSufs = {   "",   "u",   "v",   "z",  "c"};
  vector<string> namPres = {"tps", "tpp", "tpp", "tpp", "tpp"};
  for ( Index inam=0; inam<namPres.size(); ++inam ) {
    for ( Index itps=0; itps<6; ++itps ) {
      string stps = std::to_string(itps);
      string nam = namPres[inam] + stps + namSufs[inam];
      IndexRange ir = irt->get(nam);
      if ( ! ir.isValid() ) {
        cout << myname << "Invalid range: " << nam << endl;
        ++nbad;
      } else {
        cout << myname << setw(10) << ir.name << setw(20) << ir.rangeString()
             << " " << ir.label();
        for ( Index ilab=1; ilab<ir.labels.size(); ++ilab ) cout << ", " << ir.label(ilab);
        cout << endl;
        assert( ir.name == nam );
      }
    }
  }
  assert( nbad == 0 );

  cout << myname << line << endl;
  cout << myname << "Fetching APA ranges." << endl;
  nbad = 0;
  vector<string> namPresApa = {"apa", "apa", "apa", "apa", "apa"};
  for ( Index inam=0; inam<namPres.size(); ++inam ) {
    for ( Index iapa=1; iapa<=6; ++iapa ) {
      string stps = std::to_string(iapa);
      string nam = namPresApa[inam] + stps + namSufs[inam];
      IndexRange ir = irt->get(nam);
      if ( ! ir.isValid() ) {
        cout << myname << "Invalid range: " << nam << endl;
        ++nbad;
      } else {
        cout << myname << setw(10) << ir.name << setw(20) << ir.rangeString()
             << " " << ir.label();
        for ( Index ilab=1; ilab<ir.labels.size(); ++ilab ) cout << ", " << ir.label(ilab);
        cout << endl;
        assert( ir.name == nam );
      }
    }
  }
  assert( nbad == 0 );

  cout << myname << line << endl;
  cout << "Fetch bad range" << endl;
  IndexRange irb = irt->get("rangebad");
  cout << irb.rangeString() << endl;
  assert( ! irb.isValid() );

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
      cout << "Usage: " << argv[0] << " [keepFCL] [RUN]" << endl;
      cout << "  If keepFCL = true, existing FCL file is used." << endl;
      cout << "  If RUN is nonzero, the data for that run are displayed." << endl;
      return 0;
    }
    useExistingFcl = sarg == "true" || sarg == "1";
  }
  return test_ProtoDuneChannelRanges(useExistingFcl);
}

//**********************************************************************
