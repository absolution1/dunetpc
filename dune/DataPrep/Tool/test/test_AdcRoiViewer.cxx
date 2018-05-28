// test_AdcRoiViewer.cxx
//
// David Adams
// October 2017
//
// Test AdcRoiViewer.

#include <string>
#include <iostream>
#include <fstream>
#include "dune/DuneInterface/Tool/AdcChannelTool.h"
#include "dune/ArtSupport/DuneToolManager.h"

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
    fout << "#include \"dataprep_tools.fcl\"" << endl;
    fout << "tools.mytool1: {" << endl;
    fout << "  tool_type: AdcRoiViewer" << endl;
    fout << "  LogLevel: 3" << endl;
    fout << "  RoiHistOpt: 1" << endl;
    fout << "  FitOpt: 1" << endl;
    fout << "  SumHists: []" << endl;
    fout << "  RoiRootFileName: \"roi.root\"" << endl;
    fout << "  SumRootFileName: \"\"" << endl;
    fout << "}" << endl;
    fout << "tools.mytool2: @local::tools.mytool1" << endl;
    fout << "tools.mytool2.SumHists: [" << endl;
    fout << "    {var:\"fitHeight\" name:\"hfh_%0RUN%\""
         << " title:\"Fit height run %RUN%\""
         << " nbin:60 xmin:-300 xmax:300}, " << endl;
    fout << "    {var:\"fitWidth\" name:\"hfw_%0RUN%_chan%0CHAN%\""
         << " title:\"Fit width run %RUN%\""
         << " nbin:40 xmin:0 xmax:4.0}" << endl;
    fout << "]" << endl;
    fout << "tools.mytool2.RoiRootFileName: \"\"" << endl;
    fout << "tools.mytool2.SumRootFileName: \"roisum.root\"" << endl;
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
  assert( tm.toolNames().size() >= 2 );

  cout << myname << line << endl;
  cout << myname << "Fetching tool." << endl;
  auto ptoo1 = tm.getPrivate<AdcChannelTool>("mytool1");
  assert( ptoo1 != nullptr );
  auto ptoo2 = tm.getPrivate<AdcChannelTool>("mytool2");
  assert( ptoo2 != nullptr );

  cout << myname << line << endl;
  cout << myname << "Create test data." << endl;
  AdcChannelDataMap acds;
  for ( Index icha=250; icha<254; ++icha ) {
    AdcChannelData& acd = acds[icha];
    acd.run = 111;
    acd.event = 123;
    acd.channel = icha;
    unsigned int nsam = 80;
    vector<float> pulse = { 2.0, -3.0, 0.0, 5.0, 24.0, 56.0, 123.0, 71.0, 52.1, 26.3,
                           12.5,  8.1, 4.5, 2.0, -1.0,  3.2,   1.1, -2.2,  0.1, -0.1};
    acd.samples.resize(nsam);
    float chfac = float(icha)/250.0;
    for ( unsigned int ismp=0; ismp<pulse.size(); ++ismp ) {
      float smp = chfac*pulse[ismp];
      acd.samples[ismp] = smp;
      acd.samples[20+ismp] = -smp;
      acd.samples[40+ismp] = 2*smp;
      acd.samples[60+ismp] = -2*smp;
    }
    acd.rois.emplace_back( 0, 17);
    acd.rois.emplace_back(20, 37);
    acd.rois.emplace_back(40, 57);
    acd.rois.emplace_back(60, 77);
    assert( acd.rois.size() == 4 );
  }
  const AdcChannelData& acd = acds[250];

  cout << myname << line << endl;
  cout << myname << "Call tool for one channel." << endl;
  DataMap res = ptoo1->view(acd);
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
  cout << myname << "Call tool for channel map." << endl;
  res = ptoo2->viewMap(acds);
  res.print();
  cout << myname << "roiCount: " << res.getInt("roiCount") << endl;
  cout << myname << "roiHists:" << endl;
  assert( res == 0 );
  assert( res.getInt("roiCount") == 16 );

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
