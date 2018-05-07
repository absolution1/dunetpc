// test_AdcDataDumper.cxx
//
// David Adams
// April 2017
//
// Test AdcDataDumper.

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

int test_AdcDataDumper(bool useExistingFcl =false) {
  const string myname = "test_AdcDataDumper: ";
#ifdef NDEBUG
  cout << myname << "NDEBUG must be off." << endl;
  abort();
#endif
  string line = "-----------------------------";

  cout << myname << line << endl;
  string fclfile = "test_AdcDataDumper.fcl";
  if ( ! useExistingFcl ) {
    cout << myname << "Creating top-level FCL." << endl;
    ofstream fout(fclfile.c_str());
    fout << "tools: {" << endl;
    fout << "  mytool: {" << endl;
    fout << "    tool_type: AdcDataDumper" << endl;
    fout << "    FileName: \"\"" << endl;
    fout << "    Prefix: \"ADC dump for \"" << endl;
    fout << "    NewFile: false" << endl;
    fout << "    ShowChannelCount: true" << endl;
    fout << "    ShowTickCounts: true" << endl;
    fout << "    ShowRaw: true" << endl;
    fout << "    ShowPrepared: true" << endl;
    fout << "    ShowFirst: 10" << endl;
    fout << "    ShowRebin: 5" << endl;
    fout << "    ShowMax: 16" << endl;
    fout << "    ShowThreshold: 20" << endl;
    fout << "    ShowOpt: 1" << endl;
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
  auto padv = tm.getPrivate<AdcChannelTool>("mytool");
  assert( padv != nullptr );

  cout << myname << line << endl;
  cout << myname << "Create data and call too." << endl;
  AdcIndex nevt = 2;
  float peds[20] = {701, 711, 733, 690, 688, 703, 720, 720, 695, 702,
                    410, 404, 388, 389, 400, 401, 410, 404, 395, 396};
  for ( AdcIndex ievt=0; ievt<nevt; ++ievt ) {
    cout << myname << "Event " << ievt << endl;
    AdcChannelDataMap datamap;
    AdcIndex ncha = 20;
    for ( AdcIndex icha=0; icha<ncha; ++icha ) {
      std::pair<AdcChannelDataMap::iterator, bool> kdat = datamap.emplace(icha, AdcChannelData());
      assert(kdat.second);
      AdcChannelDataMap::iterator idat = kdat.first;
      AdcChannelData& data = idat->second;
      float ped = peds[icha];
      data.channel = icha;
      data.pedestal = ped;
      for ( AdcIndex itic=0; itic<100; ++itic ) {
        float xadc = ped + rand()%20 - 10.0;
        AdcCount iadc = xadc;
        data.raw.push_back(iadc);
        data.samples.push_back(iadc - ped);
      }
      AdcIndex tp = 10*ievt + 60 - 2.3*icha;
      AdcIndex tm = tp - 8;
      data.raw[tp] += 100;
      data.samples[tp] += 100;
      data.raw[tm] -= 100;
      data.samples[tm] -= 100;
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
  return test_AdcDataDumper(useExistingFcl);
}

//**********************************************************************
