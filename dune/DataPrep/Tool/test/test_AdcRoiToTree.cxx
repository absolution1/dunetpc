// test_AdcRoiToTree.cxx
//
// David Adams
// October 2017
//
// Test AdcRoiToTree.

#include <string>
#include <iostream>
#include <fstream>
#include "dune/DuneInterface/Tool/TpcDataTool.h"
#include "dune/ArtSupport/DuneToolManager.h"
#include "TRandom.h"

#undef NDEBUG
#include <cassert>

using std::string;
using std::cout;
using std::endl;
using std::ofstream;
using std::vector;
using fhicl::ParameterSet;

using Index = unsigned int;

//**********************************************************************

int test_AdcRoiToTree(bool useExistingFcl =false) {
  const string myname = "test_AdcRoiToTree: ";
#ifdef NDEBUG
  cout << myname << "NDEBUG must be off." << endl;
  abort();
#endif
  string line = "-----------------------------";

  cout << myname << line << endl;
  string fclfile = "test_AdcRoiToTree.fcl";
  if ( ! useExistingFcl ) {
    cout << myname << "Creating top-level FCL." << endl;
    ofstream fout(fclfile.c_str());
    fout << "tools.mytool: {" << endl;
    fout << "  tool_type: AdcRoiToTree" << endl;
    fout << "  LogLevel: 4" << endl;
    fout << "  OutFile: \"adcrois.root\"" << endl;
    fout << "  MetadataFields: [\"runevt\", \"xcha\"]" << endl;
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
  auto ptoo = tm.getPrivate<TpcDataTool>("mytool");
  assert( ptoo != nullptr );

  cout << myname << line << endl;
  cout << myname << "Create test data." << endl;
  Index irun = 123;
  Index ievt = 456;
  float runevt = 1000*irun + ievt;
  AdcChannelDataMap acds;
  Index ncha = 4;
  Index icha1 = 100;
  Index icha2 = icha1 + ncha;
  for ( Index icha=icha1; icha<icha2; ++icha ) {
    AdcChannelData& acd = acds[icha];
    acd.setEventInfo(irun, ievt);
    acd.setMetadata("runevt", runevt);
    float xcha = 5*(icha%100);   // Wire posn in mm.
    acd.setMetadata("xcha", xcha);
    acd.setChannelInfo(icha);
    acd.samples.resize(1000, icha);
    Index nroi = icha - icha1 + 1;
    for ( Index iroi=0; iroi<nroi; ++iroi ) {
      Index isam1 = 100*iroi;
      Index isam2 = isam1 + iroi;
      acd.rois.emplace_back(isam1, isam2);
    }
  }

  cout << myname << line << endl;
  cout << myname << "Call tool." << endl;
  DataMap res = ptoo->viewMap(acds);
  res.print();

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
  return test_AdcRoiToTree(useExistingFcl);
}

//**********************************************************************
