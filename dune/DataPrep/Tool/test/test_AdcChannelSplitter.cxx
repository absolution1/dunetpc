// test_AdcChannelSplitter.cxx
//
// David Adams
// September 2019
//
// Test AdcChannelSplitter.

#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include "dune/DuneInterface/Tool/TpcDataTool.h"
#include "dune/ArtSupport/DuneToolManager.h"
#include <TRandom.h>

#undef NDEBUG
#include <cassert>

using std::string;
using std::cout;
using std::endl;
using std::ostringstream;
using std::ofstream;
using fhicl::ParameterSet;
using std::vector;

using Name = string;
using Index = unsigned int;

//**********************************************************************

int test_AdcChannelSplitter(bool useExistingFcl =false) {
  const string myname = "test_AdcChannelSplitter: ";
#ifdef NDEBUG
  cout << myname << "NDEBUG must be off." << endl;
  abort();
#endif
  string line = "-----------------------------";

  cout << myname << line << endl;
  string fclfile = "test_AdcChannelSplitter.fcl";
  if ( ! useExistingFcl ) {
    cout << myname << "Creating top-level FCL." << endl;
    ofstream fout(fclfile.c_str());
    fout << "#include \"dataprep_tools.fcl\"" << endl;
    fout << "tools.mytool: {" << endl;
    fout << "  tool_type: AdcChannelSplitter" << endl;
    fout << "  LogLevel: 2" << endl;
    fout << "  Length: 8" << endl;
    fout << "  DataPath: \"\"" << endl;
    fout << "  DataView: split" << endl;
    fout << "}" << endl;
    fout << "tools.mytool2: @local::tools.mytool" << endl;
    fout << "tools.mytool2.Length: 3 " << endl;
    fout << "tools.mytool2.DataPath: \"split\"" << endl;
    fout << "tools.mytool2.DataView: \"split2\"" << endl;
    fout << "" << endl;
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
  cout << myname << "Fetching tools." << endl;
  auto ptoo = tm.getPrivate<TpcDataTool>("mytool");
  auto ptoo2 = tm.getPrivate<TpcDataTool>("mytool2");
  assert(ptoo);
  assert(ptoo2);

  cout << myname << line << endl;
  cout << myname << "Create data." << endl;
  AdcChannelData acd;
  acd.setChannelInfo(100123);
  acd.setEventInfo(123, 2468, 45);
  acd.sampleUnit = "fC";
  acd.raw.resize(20);
  acd.samples.resize(20);
  for ( Index isam=0; isam<acd.samples.size(); ++isam ) {
    acd.raw[isam] = 10 + isam;
    acd.samples[isam] = 100 + isam;
  }
  assert( acd.samples.size() == 20 );
  
  cout << myname << line << endl;
  cout << myname << "Split the top data." << endl;
  DataMap dmv = ptoo->update(acd);
  dmv.print();
  assert( dmv == 0 );
  assert( dmv.getInt("splitInputCount") == 1 );
  assert( dmv.getInt("splitOutputCount") == 2);
  assert( dmv.getInt("splitRawCopyCount") == 16);
  assert( dmv.getInt("splitSampleCopyCount") == 16);
  cout << myname << "Data view names:";
  for ( Name vnam : acd.viewNames() ) cout << " " << vnam;
  cout << endl;
  assert( acd.viewSize("split") == 2 );
  assert( acd.viewSize("split2") == 0 );
  for ( Name vnam : acd.viewNames() ) {
    for ( const AdcChannelData& acv : acd.view(vnam) ) {
      assert( acv.viewParent != nullptr );
      assert( acv.viewParent == &acd );
    }
  }

  cout << myname << line << endl;
  cout << myname << "Split the split." << endl;
  dmv = ptoo2->update(acd);
  dmv.print();
  assert( dmv == 0 );
  assert( dmv.getInt("splitInputCount") == 2 );
  assert( dmv.getInt("splitOutputCount") == 4);
  assert( dmv.getInt("splitRawCopyCount") == 12);
  assert( dmv.getInt("splitSampleCopyCount") == 12);
  assert( acd.viewSize("split") == 2 );
  assert( acd.viewSize("split/split") == 0 );
  assert( acd.viewSize("split/split2") == 4 );

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
  return test_AdcChannelSplitter(useExistingFcl);
}

//**********************************************************************
