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
    fout << "  SigThresh: 0.0" << endl;
    fout << "  TickBorder: 0" << endl;
    fout << "  RoiHistOpt: 1" << endl;
    fout << "  FitOpt: 1" << endl;
    fout << "  PulserStepCharge: 0.0" << endl;
    fout << "  PulserDacOffset: 0.0" << endl;
    fout << "  PulserChargeUnit: \"\"" << endl;
    fout << "  RunDataTool: \"\"" << endl;
    fout << "  TickOffsetTool: \"\"" << endl;
    fout << "  SumHists: []" << endl;
    fout << "  ChannelRanges: []" << endl;
    fout << "  ChanSumHists: []" << endl;
    fout << "  RoiRootFileName: \"roi.root\"" << endl;
    fout << "  SumRootFileName: \"\"" << endl;
    fout << "  ChanSumRootFileName: \"\"" << endl;
    fout << "}" << endl;
    fout << "tools.mytool2: @local::tools.mytool1" << endl;
    fout << "tools.mytool2.SumHists: [" << endl;
    fout << "    {var:\"fitHeight\" name:\"hfh_%0RUN%\""
         << " title:\"Combined fit height run %RUN%\""
         << " nbin:60 xmin:-300 xmax:300}, " << endl;
    fout << "    {var:\"fitHeight\" name:\"hfh_%0RUN%_chan%0CHAN%\""
         << " title:\"Fit height run %RUN% channel %CHAN%\""
         << " nbin:40 xmin:0 xmax:0 fit:gaus}," << endl;
    fout << "    {var:\"fitWidth\" name:\"hfw_%0RUN%_chan%0CHAN%\""
         << " title:\"Fit width run %RUN% channel %CHAN%\""
         << " nbin:40 xmin:0 xmax:4.0 fit:gaus}" << endl;
    fout << "]" << endl;
    fout << "tools.mytool2.ChannelRanges: [" << endl;
    fout << "  {name:apa1x label:APA1x begin:250 end:270}" << endl;
    fout << "]" << endl;
    fout << "tools.mytool2.ChanSumHists: [" << endl;
    fout << "  {name:\"hcsHeight_%CRNAME%\" title:\"Pulse heights for run %RUN% %CRLABEL%\" "
         <<    "valHist:\"hfh_%0RUN%_chan%0CHAN%\" valType:fitMean errType:fitSigma cr:apa1x}," << endl;
    fout << "  {name:\"hcsWidth_%CRNAME%\" title:\"Shaping times for run %RUN% %CRLABEL%\" "
         <<    "valHist:\"hfw_%0RUN%_chan%0CHAN%\" valType:fitMean errType:fitSigma cr:apa1x}" << endl;
    fout << "]" << endl;
    //fout << "  {val:\"hfw_%0RUN%_chan%0CHAN%:FitMean\" hist:\"hcsWidth_%0RUN%\"
    fout << "tools.mytool2.RoiRootFileName: \"\"" << endl;
    fout << "tools.mytool2.SumRootFileName: \"roisum.root\"" << endl;
    fout << "tools.mytool2.ChanSumRootFileName: \"roichan.root\"" << endl;
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
  vector<float> pulse = { 2.0, -3.0, 0.0, 5.0, 24.0, 56.0, 123.0, 71.0, 52.1, 26.3,
                         12.5,  8.1, 4.5, 2.0, -1.0,  3.2,   1.1, -2.2,  0.1, -0.1};
  Index npul = 100;
  for ( Index icha=250; icha<270; ++icha ) {
    AdcChannelData& acd = acds[icha];
    acd.run = 111;
    acd.event = 123;
    acd.channel = icha;
    acd.samples.resize(npul*pulse.size());
    acd.sampleUnit = "ADC count";
    double sigma = 5.0;
    Index itck = 0;
    for ( Index ipul=0; ipul<npul; ++ipul ) {
      Index itck0 = itck;
      for ( Index ismp=0; ismp<pulse.size(); ++ismp, ++itck ) {
        acd.samples[itck] = pulse[ismp] + gRandom->Gaus(0.0, sigma);
      }
      acd.rois.emplace_back( itck0+2, itck-2);
    }
    assert( acd.rois.size() == npul );
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
  Index nroi = res.getInt("roiCount");
  assert( nroi == npul );
  assert( res.getIntVector("roiTick0s").size() == nroi );
  assert( res.getIntVector("roiNTicks").size() == nroi );
  assert( res.getIntVector("roiNUnderflows").size() == nroi );
  assert( res.getIntVector("roiNOverflows").size() == nroi );
  assert( res.getIntVector("roiTickMins").size() == nroi );
  assert( res.getIntVector("roiTickMaxs").size() == nroi );
  assert( res.getFloatVector("roiSigMins").size() == nroi );
  assert( res.getFloatVector("roiSigMaxs").size() == nroi );
  assert( res.getFloatVector("roiSigAreas").size() == nroi );
  assert( res.getHistVector("roiHists").size() == nroi );
  assert( res.getFloatVector("roiFitHeights").size() == nroi );
  assert( res.getFloatVector("roiFitWidths").size() == nroi );
  assert( res.getFloatVector("roiFitPositions").size() == nroi );

  cout << myname << line << endl;
  cout << myname << "Call tool for channel map." << endl;
  res = ptoo2->viewMap(acds);
  res.print();
  cout << myname << "roiCount: " << res.getInt("roiCount") << endl;
  cout << myname << "roiHists:" << endl;
  assert( res == 0 );
  assert( res.getInt("roiCount") == int(npul*acds.size()) );

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
