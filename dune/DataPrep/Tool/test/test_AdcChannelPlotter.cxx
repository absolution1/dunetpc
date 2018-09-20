// test_AdcChannelPlotter.cxx
//
// David Adams
// April 2017
//
// Test AdcChannelPlotter.

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
using fhicl::ParameterSet;

//**********************************************************************

int test_AdcChannelPlotter(bool useExistingFcl =false) {
  const string myname = "test_AdcChannelPlotter: ";
#ifdef NDEBUG
  cout << myname << "NDEBUG must be off." << endl;
  abort();
#endif
  string line = "-----------------------------";

  cout << myname << line << endl;
  string fclfile = "test_AdcChannelPlotter.fcl";
  if ( ! useExistingFcl ) {
    cout << myname << "Creating top-level FCL." << endl;
    ofstream fout(fclfile.c_str());
    fout << "#include \"dataprep_tools.fcl\"" << endl;  // Need adcStringBuilder
    fout << "tools.mytool: {" << endl;
    fout << "  tool_type: AdcChannelPlotter" << endl;
    fout << "  LogLevel: 1" << endl;
    fout << "  HistTypes: [\"raw\", \"rawdist\", \"prepared\"]" << endl;
    fout << "  HistName: \"adc%TYPE%_%EVENT%_%CHAN%\"" << endl;
    fout << "  HistTitle: \"ADC %TYPE% event %EVENT% channel %CHAN%\"" << endl;
    fout << "  RootFileName: \"adcplot.root\"" << endl;
    fout << "  PlotFileName: \"adcsigs.png\"" << endl;
    fout << "  HistManager: \"\"" << endl;
    fout << "  MaxSample: 80" << endl;
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
  cout << myname << "Create data and call too." << endl;
  AdcIndex nevt = 2;
  string lab = "plane 3u";
  float peds[10] = {701.1, 711.2, 733.3, 690.4, 688.5, 703.6, 720.7, 720.8, 695.9, 702.0};
  for ( AdcIndex ievt=0; ievt<nevt; ++ievt ) {
    cout << myname << "Event " << ievt << endl;
    AdcChannelDataMap datamap;
    AdcIndex ncha = 10;
    for ( AdcIndex icha=0; icha<ncha; ++icha ) {
      std::pair<AdcChannelDataMap::iterator, bool> kdat = datamap.emplace(icha, AdcChannelData());
      assert(kdat.second);
      AdcChannelDataMap::iterator idat = kdat.first;
      AdcChannelData& data = idat->second;
      float ped = peds[icha];
      data.run = 101;
      data.subRun = 23;
      data.event = ievt;
      data.channel = icha;
      data.pedestal = ped;
      for ( AdcIndex itic=0; itic<100; ++itic ) {
        float xadc = ped + rand()%20 - 10.0;
        AdcIndex iticeff = itic - 3*icha;
        if ( iticeff > 20 && iticeff < 40 ) xadc +=600;
        AdcCount iadc = xadc;
        data.raw.push_back(iadc);
        data.flags.push_back(0);
        data.samples.push_back(iadc - ped);
        data.signal.push_back(xadc - ped > 300.0 );
      }
      AdcIndex tp = 10*ievt + 60 - 2.3*icha;
      AdcIndex tm = tp - 8;
      data.raw[tp] += 100;
      data.samples[tp] += 100;
      data.raw[tm] -= 100;
      data.samples[tm] -= 100;
      data.flags[tm] = 4;
      data.roisFromSignal();
      assert( padv->view(datamap[icha]) == 0 );
    }
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
  return test_AdcChannelPlotter(useExistingFcl);
}

//**********************************************************************
