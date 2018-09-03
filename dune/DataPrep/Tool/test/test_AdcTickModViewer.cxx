// test_AdcTickModViewer.cxx
//
// David Adams
// April 2017
//
// Test AdcTickModViewer.

#include <string>
#include <iostream>
#include <fstream>
#include <vector>
#include "dune/DuneInterface/Tool/AdcChannelTool.h"
#include "dune/ArtSupport/DuneToolManager.h"
#include "TError.h"

#undef NDEBUG
#include <cassert>

using std::string;
using std::cout;
using std::endl;
using std::ofstream;
using std::vector;
using fhicl::ParameterSet;

//**********************************************************************

void TestErrorHandler(Int_t, Bool_t, const char*, const char* msg) {
  const string myname = "TestErrorHandler: ";
  cout << myname << msg << endl;
  abort();
}

//**********************************************************************

int test_AdcTickModViewer(bool useExistingFcl, bool doUpdate, bool doUpdateMap) {
  const string myname = "test_AdcTickModViewer: ";
#ifdef NDEBUG
  cout << myname << "NDEBUG must be off." << endl;
  abort();
#endif
  string line = "-----------------------------";

  // Set a Root error handler and make sure it is not called when
  // the tool is used.
  SetErrorHandler(TestErrorHandler);

  cout << myname << line << endl;
  string fclfile = "test_AdcTickModViewer.fcl";
  if ( ! useExistingFcl ) {
    cout << myname << "Creating top-level FCL." << endl;
    ofstream fout(fclfile.c_str());
    fout << "#include \"dataprep_tools.fcl\"" << endl;   // Need adcNameManipulator
    fout << "tools.mytool: {" << endl;
    fout << "  tool_type: AdcTickModViewer" << endl;
    fout << "  LogLevel: 1" << endl;
    fout << "  TickModPeriod: 4" << endl;
    fout << "  TimeOffsetTool: myTimeOffsetTool" << endl;
    fout << "  FitSigmaMin: 1.0" << endl;
    fout << "  FitSigmaMax: 20.0" << endl;
    fout << "  HistName: \"adctm_ch%0CHAN%_tm%0TICKMOD%\"" << endl;
    fout << "  HistTitle: \"ADC spectrum for channel %CHAN% tickmod %TICKMOD%\"" << endl;
    fout << "  HistChannelCount: 100" << endl;
    fout << "  PlotChannels: []" << endl;
    fout << "  AllPlotFileName: \"adctm%TICKMOD%_ch%0CHAN%.png\"" << endl;
    fout << "  MinPlotFileName: \"adctmMin_ch%0CHAN%.png\"" << endl;
    fout << "  MaxPlotFileName: \"adctmMax_ch%0CHAN%.png\"" << endl;
    fout << "  PhaseGrouping: channel" << endl;
    fout << "  PhasePlotFileName: \"adcPhase_ch%0CHAN%.png\"" << endl;
    fout << "  PhaseVariable: phase" << endl;
    fout << "  RootFileName: \"adctm.root\"" << endl;
    fout << "  TreeFileName: \"tickmod.root\"" << endl;
    fout << "  TreeFileName: \"\"" << endl;
    fout << "  PlotSizeX: 1400" << endl;
    fout << "  PlotSizeY: 1000" << endl;
    fout << "  PlotShowFit:  0" << endl;
    fout << "  PlotSplitX:  2" << endl;
    fout << "  PlotSplitY:  2" << endl;
    fout << "  PlotFrequency: 0" << endl;
    fout << "  PhasePlotSizeX: 1400" << endl;
    fout << "  PhasePlotSizeY: 1000" << endl;
    fout << "  PhasePlotSplitX:  2" << endl;
    fout << "  PhasePlotSplitY:  2" << endl;
    fout << "}" << endl;
    fout << "tools.myTimeOffsetTool: {" << endl;
    fout << "  tool_type: FixedTimeOffsetTool" << endl;
    fout << "  LogLevel: 2" << endl;
    fout << "  Value: 1000" << endl;
    fout << "  Rem: 0.0" << endl;
    fout << "  Unit: tick" << endl;
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
  auto pvtm = tm.getPrivate<AdcChannelTool>("mytool");
  assert( pvtm != nullptr );
  if ( ! doUpdate ) pvtm = nullptr;

  cout << myname << line << endl;
  cout << myname << "Create data." << endl;
  AdcIndex nevt = 2;
  string lab = "plane 3u";
  float peds[10] = {701.1, 711.2, 733.3, 690.4, 688.5, 703.6, 720.7, 720.8, 695.9, 702.0};
  vector<AdcChannelDataMap> acms(nevt);
  AdcIndex ncha = 10;
  for ( AdcIndex ievt=0; ievt<nevt; ++ievt ) {
    cout << myname << "Event " << ievt << endl;
    AdcChannelDataMap& datamap = acms[ievt];
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
        if ( iticeff > 20 && iticeff < 45 ) xadc +=600;
        AdcCount iadc = xadc;
        data.raw.push_back(iadc);
        data.flags.push_back(0);
        data.samples.push_back(iadc - ped);
        data.signal.push_back(xadc - ped > 300.0 );
      }
      data.pedestal += 10.0;
      AdcIndex tp = 10*ievt + 60 - 2.3*icha;
      AdcIndex tm = tp - 8;
      data.raw[tp] += 100;
      data.samples[tp] += 100;
      data.raw[tm] -= 100;
      data.samples[tm] -= 100;
      data.flags[tm] = 4;
      data.roisFromSignal();
    }
  }

  cout << myname << line << endl;
  cout << myname << "Call tool." << endl;
  unsigned int expCount = 0;
  for ( const AdcChannelDataMap& acm : acms ) {
    cout << myname << "Event " << acm.begin()->second.event << endl;
    expCount += 25;
    for ( const AdcChannelDataMap::value_type& iacd : acm ) {
      const AdcChannelData& acd = iacd.second;
      cout << myname << "Event " << acd.event << ", channel " << acd.channel << endl;
      cout << myname << "ADC channel data size: " << acd.raw.size() << endl;
      DataMap dm = pvtm->view(acd);
      dm.print(myname);
      assert( dm.status() == 0 );
      assert( dm.haveHistVector("tmHists") );
      const DataMap::HistVector& phs = dm.getHistVector("tmHists");
      assert( phs.size() == 4 );
      for ( const TH1* ph : phs ) {
        cout << myname << "  " << ph->GetName() << ": " << ph->GetTitle() << endl;
        double countSum = ph->Integral();
        double countUnder = ph->GetBinContent(0);
        double countOver = ph->GetBinContent(ph->GetNbinsX()+1);
        cout << myname << "    Hist integral:" << countSum << endl;
        cout << myname << "    Hist undrflow:" << countUnder << endl;
        cout << myname << "    Hist overflow:" << countOver << endl;
        double countTot = countSum + countUnder + countOver;
        assert( countTot == expCount );
      }
    }
  }

  cout << myname << line << endl;
  cout << myname << "Done." << endl;
  return 0;
}

//**********************************************************************

int main(int argc, char* argv[]) {
  bool useExistingFcl = false;
  bool doUpdate = true;
  bool doUpdateMap = true;
  if ( argc > 1 ) {
    string sarg(argv[1]);
    if ( sarg == "-h" ) {
      cout << "Usage: " << argv[0] << " [ARG]" << endl;
      cout << "  If ARG = true, existing FCL file is used." << endl;
      return 0;
    }
    useExistingFcl = sarg == "true" || sarg == "1";
  }
  return test_AdcTickModViewer(useExistingFcl, doUpdate, doUpdateMap);
}

//**********************************************************************
