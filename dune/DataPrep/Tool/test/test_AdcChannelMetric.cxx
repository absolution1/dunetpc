// test_AdcChannelMetric.cxx
//
// David Adams
// April 2017
//
// Test AdcChannelMetric.

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

int test_AdcChannelMetric(bool useExistingFcl =false) {
  const string myname = "test_AdcChannelMetric: ";
#ifdef NDEBUG
  cout << myname << "NDEBUG must be off." << endl;
  abort();
#endif
  string line = "-----------------------------";

  cout << myname << line << endl;
  string fclfile = "test_AdcChannelMetric.fcl";
  string hname = "hchped_all";
  if ( ! useExistingFcl ) {
    cout << myname << "Creating top-level FCL." << endl;
    ofstream fout(fclfile.c_str());
    fout << "#include \"dataprep_tools.fcl\"" << endl;   // Need adcNameManipulator
    fout << "#include \"protodune_dataprep_tools.fcl\"" << endl;   // Need adcNameManipulator
    fout << "tools.mytool: {" << endl;
    fout << "           tool_type: AdcChannelMetric" << endl;
    fout << "            LogLevel: 3" << endl;
    fout << "              Metric: \"pedestal\"" << endl;
    fout << "       ChannelRanges: [all, tpp3c, tpp3z]" << endl;
    fout << "           MetricMin: 0.0" << endl;
    fout << "           MetricMax: 2000.0" << endl;
    fout << "  ChannelLineModulus:  40" << endl;
    fout << "  ChannelLinePattern:  [10]" << endl;
    fout << "            HistName: \"hchped_%CRNAME%\"" << endl;
    fout << "           HistTitle: \"ADC pedestals for run %RUN% event %EVENT% %CRLABEL%\"" << endl;
    fout << "         MetricLabel: \"Pedestal [ADC counts]\"" << endl;
    fout << "           PlotSizeX: 0" << endl;
    fout << "           PlotSizeY: 0" << endl;
    fout << "        PlotFileName: \"mypeds-run%RUN%-evt%EVENT%_%CRNAME%.png\"" << endl;
    fout << "        RootFileName: \"\"" << endl;
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
  const int nped = 20;
  float peds[nped] = {701, 711, 733, 690, 688, 703, 720, 720, 695, 702,
                    410, 404, 388, 389, 400, 401, 410, 404, 395, 396};
  vector<double> wf = {5.0, 20.1, 53.2, 80.6, 130.2, 160.1, 150.4, 125.7, 72.5, 41.3, 18.4,
                       -6.5, -34.9, -56.6, -88.9, -132.6, -170.8, -172.9, -144.6, -112.6,
                        -79.4, -44.9, -22.1, -12.6, -4.7};
  for ( AdcIndex ievt=0; ievt<nevt; ++ievt ) {
    cout << myname << "Event " << ievt << endl;
    AdcChannelDataMap datamap;
    AdcIndex ncha = 200;
    AdcIndex icha1 = 10000;
    AdcIndex icha2 = icha1 + ncha;
    for ( AdcIndex icha=icha1; icha<icha2; ++icha ) {
      std::pair<AdcChannelDataMap::iterator, bool> kdat = datamap.emplace(icha, AdcChannelData());
      assert(kdat.second);
      AdcChannelDataMap::iterator idat = kdat.first;
      AdcChannelData& data = idat->second;
      data.run = 123;
      data.event = ievt;
      float ped = peds[(icha-icha1)%nped];
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
    DataMap ret = padv->viewMap(datamap);
    ret.print();
    assert( ret == 0 );
    cout << myname << "Checking histogram " << hname << endl;
    TH1* phout = ret.getHist(hname);
    assert( phout != nullptr );
    assert( phout->GetName() == hname );
    assert( phout->GetEntries() == ncha );
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
  return test_AdcChannelMetric(useExistingFcl);
}

//**********************************************************************
