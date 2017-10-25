// test_AdcRoiViewer.cxx
//
// David Adams
// October 2017
//
// Test AdcRoiViewer.

#include <string>
#include <iostream>
#include <fstream>
#include "dune/DuneInterface/Tool/AdcChannelViewer.h"
#include "dune/ArtSupport/DuneToolManager.h"

#undef NDEBUG
#include <cassert>

using std::string;
using std::cout;
using std::endl;
using std::ofstream;
using std::vector;
using fhicl::ParameterSet;

//**********************************************************************

int test_AdcRoiViewer(bool useExistingFcl =false) {
  const string myname = "test_AdcRoiViewer: ";
#ifdef NDEBUG
  cout << myname << "NDEBUG must be off." << endl;
  abort();
#endif
  string line = "-----------------------------";

  cout << myname << line << endl;
  string fclfile = "test_AdcRoiViewer.fcl";
  if ( ! useExistingFcl ) {
    cout << myname << "Creating top-level FCL." << endl;
    ofstream fout(fclfile.c_str());
    fout << "tools: {" << endl;
    fout << "  mytool: {" << endl;
    fout << "    tool_type: AdcRoiViewer" << endl;
    fout << "    LogLevel: 1" << endl;
    fout << "    HistOpt: 1" << endl;
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
  auto ptoo = tm.getPrivate<AdcChannelViewer>("mytool");
  assert( ptoo != nullptr );

  cout << myname << line << endl;
  cout << myname << "Create test data." << endl;
  AdcChannelData acd;
  unsigned int nsam = 80;
  vector<float> pulse = { 2.0, -3.0, 0.0, 5.0, 24.0, 56.0, 123.0, 71.0, 52.1, 26.3,
                         12.5,  8.1, 4.5, 2.0, -1.0,  3.2,   1.1, -2.2,  0.1, -0.1};
  acd.samples.resize(nsam);
  for ( unsigned int ismp=0; ismp<pulse.size(); ++ismp ) {
    float smp = pulse[ismp];
    acd.samples[ismp] = smp;
    acd.samples[20+ismp] = -smp;
    acd.samples[ismp] = 2*smp;
    acd.samples[ismp] = -2*smp;
  }
  acd.rois.emplace_back( 0, 17);
  acd.rois.emplace_back(20, 37);
  acd.rois.emplace_back(40, 57);
  acd.rois.emplace_back(60, 77);
  assert( acd.rois.size() == 4 );

  cout << myname << line << endl;
  cout << myname << "Call tool." << endl;
  DataMap res = ptoo->view(acd);
  res.print();
  cout << myname << "roiCount: " << res.getInt("roiCount") << endl;
  cout << myname << "roiHists:" << endl;
  for ( TH1* ph : res.getHistVector("roiHists") ) {
    cout << myname << "  " << ph->GetName() << endl;
  }
  assert( res == 0 );
  assert( res.getInt("roiCount") == 4 );
  assert( res.getHistVector("roiHists").size() == 4 );

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
  return test_AdcRoiViewer(useExistingFcl);
}

//**********************************************************************
