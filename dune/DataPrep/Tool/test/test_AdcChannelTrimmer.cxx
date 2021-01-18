// test_AdcChannelTrimmer.cxx
//
// David Adams
// April 2017
//
// Test AdcChannelTrimmer.

#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include "dune/DuneInterface/Tool/AdcChannelTool.h"
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

using Index = unsigned int;

//**********************************************************************

int test_AdcChannelTrimmer(bool useExistingFcl =false) {
  const string myname = "test_AdcChannelTrimmer: ";
#ifdef NDEBUG
  cout << myname << "NDEBUG must be off." << endl;
  abort();
#endif
  string line = "-----------------------------";

  cout << myname << line << endl;
  string fclfile = "test_AdcChannelTrimmer.fcl";
  if ( ! useExistingFcl ) {
    cout << myname << "Creating top-level FCL." << endl;
    ofstream fout(fclfile.c_str());
    fout << "#include \"dataprep_tools.fcl\"" << endl;
    fout << "tools.mytool: {" << endl;
    fout << "  tool_type: AdcChannelTrimmer" << endl;
    fout << "  LogLevel: 2" << endl;
    fout << "  Length: 20" << endl;
    fout << "  MaxTrim: 10" << endl;
    fout << "}" << endl;
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
  cout << myname << "Fetching tool." << endl;
  auto ptoo = tm.getPrivate<AdcChannelTool>("mytool");

  cout << myname << line << endl;
  cout << myname << "Create data." << endl;
  AdcChannelData acd;
  acd.channel = 100123;
  acd.setEventInfo(123, 2468, 45);
  acd.sampleUnit = "fC";
  acd.samples.resize(15);
  for ( Index isam=0; isam<acd.samples.size(); ++isam ) {
    acd.samples[isam] = 1.0 + 0.01*isam;
  }
  assert( acd.samples.size() == 15 );
  
  cout << myname << line << endl;
  cout << myname << "Call tool view." << endl;
  DataMap dmv = ptoo->view(acd);
  dmv.print();
  assert( dmv == 0 );
  assert( dmv.getInt("trimAction") == -1 );
  assert( dmv.getInt("trimLength") == -5 );
  assert( acd.samples.size() == 15 );

  cout << myname << line << endl;
  cout << myname << "Call tool update." << endl;
  Index len1 = acd.samples.size();
  DataMap dmu = ptoo->update(acd);
  Index len2 = acd.samples.size();
  dmu.print();
  assert( dmu == 0 );
  assert( dmu.getInt("trimAction") == 1 );
  assert( dmu.getInt("trimLength") == -5 );
  assert( len1 == 15 );
  assert( len2 == 20 );
  for ( Index isam=len1; isam<len2; ++isam ) {
    float qsam = acd.samples[isam];
    float qchk = acd.samples[isam-len1];
    assert( fabs(qsam-qchk) < 1.e-6 );
  }

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
  return test_AdcChannelTrimmer(useExistingFcl);
}

//**********************************************************************
