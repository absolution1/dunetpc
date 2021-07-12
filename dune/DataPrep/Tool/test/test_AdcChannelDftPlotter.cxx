// test_AdcChannelDftPlotter.cxx
//
// David Adams
// April 2017
//
// Test AdcChannelDftPlotter.

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

using Index = unsigned int;
using Name = string;

//**********************************************************************

int test_AdcChannelDftPlotter(bool useExistingFcl =false) {
  const string myname = "test_AdcChannelDftPlotter: ";
#ifdef NDEBUG
  cout << myname << "NDEBUG must be off." << endl;
  abort();
#endif
  string line = "-----------------------------";

  cout << myname << line << endl;
  string fclfile = "test_AdcChannelDftPlotter.fcl";
  Name vnam = "myview";
  if ( ! useExistingFcl ) {
    cout << myname << "Creating top-level FCL." << endl;
    ofstream fout(fclfile.c_str());
    fout << "#include \"dataprep_tools.fcl\"" << endl;
    fout << "tools.mytemplate: {" << endl;
    fout << "  tool_type: AdcChannelDftPlotter" << endl;
    fout << "  LogLevel: 2" << endl;
    fout << "  ChannelStatusFlag: 0" << endl;
    fout << "  ChannelSelection: \"\"" << endl;
    fout << "  SampleFreq: 2000" << endl;
    fout << "  XMin: 0.0" << endl;
    fout << "  XMax: 0.0" << endl;
    fout << "  YMinLog: 0.0" << endl;
    fout << "  HistName: \"\"" << endl;
    fout << "  HistTitle: \"\"" << endl;
    fout << "  HistSummaryTitles: []" << endl;
    fout << "  PlotDataView: \"\"" << endl;
    fout << "  PlotName: \"\"" << endl;
    fout << "  PlotSummaryNames: []" << endl;
    fout << "  PlotChannelRanges: []" << endl;
    fout << "  PlotChannelGroups: []" << endl;
    fout << "  PlotOverlayGroups: 1" << endl;
    fout << "  PlotSizeX: 1000" << endl;
    fout << "  PlotSizeY:  700" << endl;
    fout << "  PlotSplitX: 2" << endl;
    fout << "  PlotSplitY: 0" << endl;
    fout << "}" << endl;
    fout << "" << endl;
    fout << "tools.myphases: @local::tools.mytemplate" << endl;
    fout << "tools.myphases.Variable: phase" << endl;
    fout << "tools.myphases.HistName: \"hdftphase_run%0RUN%_evt%0EVENT%\"" << endl;
    fout << "tools.myphases.HistTitle: \"DFT phases for run %RUN% event %EVENT% channel %CHAN%\"" << endl;
    fout << "tools.myphases.PlotName: \"dftphase_run%0RUN%_evt%0EVENT%.png\"" << endl;
    fout << "" << endl;
    fout << "tools.mymags: @local::tools.mytemplate" << endl;
    fout << "tools.mymags.Variable: magnitude" << endl;
    fout << "tools.mymags.YMax: 0.0" << endl;
    fout << "tools.mymags.HistName: \"hdftmags_run%0RUN%_evt%0EVENT%\"" << endl;
    fout << "tools.mymags.HistTitle: \"DFT amplitudes for run %RUN% event %EVENT% channel %CHAN%\"" << endl;
    fout << "tools.mymags.PlotName: \"dftmag_run%0RUN%_evt%0EVENT%.png\"" << endl;
    fout << "" << endl;
    fout << "tools.mypwr: @local::tools.mytemplate" << endl;
    fout << "tools.mypwr.Variable: power" << endl;
    fout << "tools.mypwr.YMax: 65.0" << endl;
    fout << "tools.mypwr.NBinX: 5" << endl;
    fout << "tools.mypwr.HistName: \"hdftpower_run%0RUN%_evt%0EVENT%\"" << endl;
    fout << "tools.mypwr.HistTitle: \"DFT power for run %RUN% event %EVENT% channel %CHAN%\"" << endl;
    fout << "tools.mypwr.PlotName: \"dftpower_run%0RUN%_evt%0EVENT%_ch%0CHAN%.png\"" << endl;
    fout << "" << endl;
    fout << "tools.mypwt: @local::tools.mytemplate" << endl;
    fout << "tools.mypwt.Variable: \"power/tick\"" << endl;
    fout << "tools.mypwt.YMax: 6.0" << endl;
    fout << "tools.mypwt.NBinX: 5" << endl;
    fout << "tools.mypwt.HistName: \"hdftpowt_run%0RUN%_evt%0EVENT%\"" << endl;
    fout << "tools.mypwt.HistTitle: \"DFT power for run %RUN% event %EVENT% channel %CHAN%\"" << endl;
    fout << "tools.mypwt.PlotName: \"dftpowt_run%0RUN%_evt%0EVENT%_ch%0CHAN%.png\"" << endl;
    fout << "tools.mvpwr: @local::tools.mypwr" << endl;
    fout << "tools.mvpwr.PlotDataView: " << vnam << endl;
    fout << "tools.mvpwr.HistTitle: \"DFT power for run %RUN% event %EVENT% channel %CHAN% view %VIEW%\"" << endl;
    fout << "tools.mvpwr.PlotName: \"dftpower_run%0RUN%_evt%0EVENT%_ch%0CHAN%_view%VIEW%.png\"" << endl;
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
  auto ppha = tm.getPrivate<TpcDataTool>("myphases");
  auto pmag = tm.getPrivate<TpcDataTool>("mymags");
  auto ppwr = tm.getPrivate<TpcDataTool>("mypwr");
  auto ppwt = tm.getPrivate<TpcDataTool>("mypwt");
  auto ppwrv = tm.getPrivate<TpcDataTool>("mvpwr");
  assert( ppha != nullptr );
  assert( pmag != nullptr );
  assert( ppwr != nullptr );
  assert( ppwt != nullptr );
  assert( ppwrv != nullptr );

  cout << myname << line << endl;
  cout << myname << "Create data." << endl;
  AdcChannelData acd;
  acd.setChannelInfo(123456);
  acd.setEventInfo(123, 2468, 45);
  acd.sampleUnit = "fC";
  float pi = acos(-1.0);
  vector<float> mags = {  1.0, 2.0, 3.0, 4.0, 5.0, 4.0, 2.0,  1.0,  1.0,  1.0,  1.0 };
  vector<float> phas = {  0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, -2.0, -1.5, -1.0 };
  acd.dftmags = mags;
  acd.dftphases = phas;
  cout << myname << "  Mag size: " << mags.size() << endl;
  cout << myname << "Phase size: " << phas.size() << endl;

  cout << myname << line << endl;
  cout << myname << "Create data map." << endl;
  AdcChannelDataMap acds;
  for ( Index icha=10000; icha<10010; ++icha ) {
    acds[icha].setEventInfo(acd.run(), acd.event(), acd.subRun());;
    acds[icha].setChannelInfo(icha);
    acds[icha].dftmags = mags;
    acds[icha].dftphases = phas;
    for ( float& mag : mags ) mag *= 1.02;
    for ( float& pha : phas ) { pha += 0.2; if ( pha > pi ) pha -= 2.0*pi; }
    // Build the view vnam.
    vector<float> vmags = mags;
    for ( float& mag : vmags ) mag *= sqrt(0.70);
    vector<float> vphas = phas;
    for ( float& pha : vphas ) { pha += 0.4; if ( pha > pi ) pha -= 2.0*pi; }
    acds[icha].updateView(vnam).push_back(acds[icha]);
    acds[icha].updateView(vnam).back().dftmags = vmags;
    acds[icha].updateView(vnam).back().dftphases = vphas;
    vmags = mags;
    for ( float& mag : vmags ) mag *= sqrt(0.30);
    vphas = phas;
    for ( float& pha : vphas ) { pha += 0.4; if ( pha > pi ) pha -= 2.0*pi; }
    acds[icha].updateView(vnam).push_back(acds[icha]);
    acds[icha].updateView(vnam).back().dftmags = vmags;
    acds[icha].updateView(vnam).back().dftphases = vphas;
  }

  cout << myname << line << endl;
  cout << myname << "Call phase tool." << endl;
  DataMap dm = ppha->view(acd);
  dm.print();
  assert( dm == 0 );

  cout << myname << line << endl;
  cout << myname << "Call mag tool." << endl;
  dm = pmag->view(acd);
  dm.print();
  assert( dm == 0 );

  cout << myname << line << endl;
  cout << myname << "Call power tool." << endl;
  dm = ppwr->view(acd);
  dm.print();
  assert( dm == 0 );

  cout << myname << line << endl;
  cout << myname << "Call tick power tool." << endl;
  dm = ppwt->view(acd);
  dm.print();
  assert( dm == 0 );

  cout << myname << line << endl;
  cout << myname << "Call power tool with map." << endl;
  dm = ppwr->viewMap(acds);
  dm.print();
  assert( dm == 0 );

  cout << myname << line << endl;
  cout << myname << "Call tick power tool with map." << endl;
  dm = ppwt->viewMap(acds);
  dm.print();
  assert( dm == 0 );

  cout << myname << line << endl;
  cout << myname << "Call tick power view tool with map." << endl;
  dm = ppwrv->viewMap(acds);
  dm.print();
  assert( dm == 0 );

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
  return test_AdcChannelDftPlotter(useExistingFcl);
}

//**********************************************************************
