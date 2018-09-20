// test_AdcDataPlotter.cxx
//
// David Adams
// April 2017
//
// Test AdcDataPlotter.

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

//**********************************************************************

int test_AdcDataPlotter(bool useExistingFcl =false) {
  const string myname = "test_AdcDataPlotter: ";
#ifdef NDEBUG
  cout << myname << "NDEBUG must be off." << endl;
  abort();
#endif
  string line = "-----------------------------";

  cout << myname << line << endl;
  string fclfile = "test_AdcDataPlotter.fcl";
  if ( ! useExistingFcl ) {
    cout << myname << "Creating top-level FCL." << endl;
    ofstream fout(fclfile.c_str());
    fout << "#include \"dataprep_tools.fcl\"" << endl;  // Need adcStringBuilder
    fout << "#include \"dunecommon_tools.fcl\"" << endl;  // Need tickRanges
    fout << "tools.tickRanges.myTicks: { begin:10 end:90 labels:[\"My ticks\"] }" << endl;
    fout << "tools.mytool: {" << endl;
    fout << "           tool_type: AdcDataPlotter" << endl;
    fout << "            DataType: 0" << endl;
    fout << "            LogLevel: 2" << endl;
    fout << "           TickRange: \"myTicks\"" << endl;
    fout << "       ChannelRanges: []" << endl;
    fout << "     FembTickOffsets: []" << endl;
    fout << "           MaxSignal: 10" << endl;
    fout << "  ChannelLineModulus:  4" << endl;
    fout << "  ChannelLinePattern:  [1]" << endl;
    fout << "             Palette: 1026" << endl;
    fout << "            HistName: \"hadc\"" << endl;
    fout << "           HistTitle: \"Prepared ADC run %RUN% event %EVENT%\"" << endl;
    fout << "        PlotFileName: \"myplot-run%0RUN%-evt%0EVENT%.png\"" << endl;
    fout << "           PlotSizeX: 0" << endl;
    fout << "           PlotSizeY: 0" << endl;
    fout << "        RootFileName: \"adc.root\"" << endl;
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
  auto padv = tm.getPrivate<AdcChannelTool>("mytool");
  assert( padv != nullptr );

  cout << myname << line << endl;
  cout << myname << "Create data and call tool." << endl;
  AdcIndex nevt = 2;
  float peds[20] = {701, 711, 733, 690, 688, 703, 720, 720, 695, 702,
                    410, 404, 388, 389, 400, 401, 410, 404, 395, 396};
  vector<double> wf = {5.0, 20.1, 53.2, 80.6, 130.2, 160.1, 150.4, 125.7, 72.5, 41.3, 18.4,
                       -6.5, -34.9, -56.6, -88.9, -132.6, -170.8, -172.9, -144.6, -112.6,
                        -79.4, -44.9, -22.1, -12.6, -4.7};
  for ( AdcIndex ievt=0; ievt<nevt; ++ievt ) {
    cout << myname << "Event " << ievt << endl;
    AdcChannelDataMap datamap;
    AdcIndex ncha = 20;
    AdcIndex icha1 = 10000;
    AdcIndex icha2 = icha1 + ncha;
    for ( AdcIndex icha=icha1; icha<icha2; ++icha ) {
      std::pair<AdcChannelDataMap::iterator, bool> kdat = datamap.emplace(icha, AdcChannelData());
      assert(kdat.second);
      AdcChannelDataMap::iterator idat = kdat.first;
      AdcChannelData& data = idat->second;
      data.run = 123;
      data.event = ievt;
      float ped = peds[icha-icha1];
      data.channel = icha;
      data.pedestal = ped;
      for ( AdcIndex itic=0; itic<100; ++itic ) {
        float xadc = ped + gRandom->Gaus(0.0, 10.0);
        AdcCount iadc = xadc;
        data.raw.push_back(iadc);
        data.samples.push_back(iadc - ped);
      }
      AdcIndex tp = 10*ievt + 60 - 2.3*(icha-icha1);
      for ( unsigned int iwf=0; iwf<wf.size(); ++iwf ) {
        unsigned int isam = tp+iwf;
        if ( isam < data.samples.size() ) {
          data.raw[isam] += wf[iwf];
          data.samples[isam] += wf[iwf];
        }
      }
      for ( unsigned int isam=0; isam<data.samples.size(); ++isam ) {
        data.samples[isam] *= 0.04;
      }
      data.sampleUnit = "ke";
    }
    assert( padv->viewMap(datamap) == 0 );
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
  return test_AdcDataPlotter(useExistingFcl);
}

//**********************************************************************
