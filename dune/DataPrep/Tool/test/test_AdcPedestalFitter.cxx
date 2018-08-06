// test_AdcPedestalFitter.cxx
//
// David Adams
// April 2017
//
// Test AdcPedestalFitter.

#include <string>
#include <iostream>
#include <fstream>
#include "dune/DuneInterface/Tool/AdcChannelTool.h"
#include "dune/ArtSupport/DuneToolManager.h"
#include "TError.h"

#undef NDEBUG
#include <cassert>

using std::string;
using std::cout;
using std::endl;
using std::ofstream;
using fhicl::ParameterSet;

//**********************************************************************

void TestErrorHandler(Int_t, Bool_t, const char*, const char* msg) {
  const string myname = "TestErrorHandler: ";
  cout << myname << msg << endl;
  abort();
}

//**********************************************************************

int test_AdcPedestalFitter(bool useExistingFcl, bool doUpdate, bool doUpdateMap) {
  const string myname = "test_AdcPedestalFitter: ";
#ifdef NDEBUG
  cout << myname << "NDEBUG must be off." << endl;
  abort();
#endif
  string line = "-----------------------------";

  // Set a Root error handler and make sure it is not called when
  // the tool is used.
  SetErrorHandler(TestErrorHandler);

  cout << myname << line << endl;
  string fclfile = "test_AdcPedestalFitter.fcl";
  if ( ! useExistingFcl ) {
    cout << myname << "Creating top-level FCL." << endl;
    ofstream fout(fclfile.c_str());
    fout << "#include \"dataprep_tools.fcl\"" << endl;   // Need adcNameManipulator
    fout << "tools.mytool: {" << endl;
    fout << "  tool_type: AdcPedestalFitter" << endl;
    fout << "  LogLevel: 1" << endl;
    fout << "  FitRmsMin: 1.0" << endl;
    fout << "  FitRmsMax: 20.0" << endl;
    fout << "  HistName: \"adcped_%EVENT%_%CHAN%\"" << endl;
    fout << "  HistTitle: \"ADC pedestal for event %EVENT% channel %CHAN%\"" << endl;
    fout << "  PlotFileName: \"adcped_ev%EVENT%_chan%CHAN%.png\"" << endl;
    fout << "  RootFileName: \"adcped.root\"" << endl;
    fout << "  HistManager: \"\"" << endl;
    fout << "  PlotSizeX:  700" << endl;
    fout << "  PlotSizeY:  500" << endl;
    fout << "  PlotShowFit:  2" << endl;
    fout << "  PlotSplitX:  0" << endl;
    fout << "  PlotSplitY:  0" << endl;
    fout << "}" << endl;
    fout << "tools.mymaptool: @local::tools.mytool" << endl;
    fout << "tools.mymaptool.LogLevel: 2" << endl;
    fout << "tools.mymaptool.PlotFileName: \"adcpedmap_ev%EVENT%_chan%CHAN%.png\"" << endl;
    fout << "tools.mymaptool.PlotSplitX: 2" << endl;
    fout << "tools.mymaptool.RootFileName: \"\"" << endl;
    fout << "tools.mymaptool.PlotSizeX: 1400" << endl;
    fout << "tools.mymaptool.PlotSizeY: 1000" << endl;
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
  assert( tm.toolNames().size() > 1 );

  cout << myname << line << endl;
  cout << myname << "Fetching histogram manaager." << endl;
  auto phm = tm.getShared<AdcChannelTool>("mytool");
  assert( phm != nullptr );

  cout << myname << line << endl;
  cout << myname << "Fetching tool." << endl;
  auto padvNotUsed = tm.getPrivate<AdcChannelTool>("mytool");
  assert( padvNotUsed != nullptr );
  auto padvsin = tm.getPrivate<AdcChannelTool>("mytool");
  assert( padvsin != nullptr );
  auto padvmap = tm.getPrivate<AdcChannelTool>("mymaptool");
  assert( padvmap != nullptr );
  if ( ! doUpdate ) padvsin = nullptr;
  if ( ! doUpdateMap ) padvmap = nullptr;

  cout << myname << line << endl;
  cout << myname << "Create data and call tool." << endl;
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
      data.pedestal += 10.0;
      AdcIndex tp = 10*ievt + 60 - 2.3*icha;
      AdcIndex tm = tp - 8;
      data.raw[tp] += 100;
      data.samples[tp] += 100;
      data.raw[tm] -= 100;
      data.samples[tm] -= 100;
      data.flags[tm] = 4;
      data.roisFromSignal();
      if ( padvsin != nullptr ) {
        double ped0 = datamap[icha].pedestal;
        //assert( padv->view(datamap[icha]) == 0 );
        double ped1 = datamap[icha].pedestal;
        assert( ! datamap[icha].hasMetadata("fitPedPeakBinFraction") );
        assert( padvsin->update(datamap[icha]) == 0 );
        double ped2 = datamap[icha].pedestal;
        cout << "Old pedestal: " << ped0 << endl;
        cout << "New pedestal: " << ped2 << endl;
        assert( ped1 == ped0 );
        assert( ped2 != ped0 );
        assert( ped2 != 0.0 );
        assert( datamap[icha].hasMetadata("fitPedPeakBinFraction") );
        assert( datamap[icha].hasMetadata("fitPedFractionLow") );
        assert( datamap[icha].hasMetadata("fitPedFractionHigh") );
        //assert( fabs(ped2-ped) < 0.01 );
      }
    }
    if ( padvmap != nullptr ) {
      for ( AdcIndex icha=0; icha<ncha; ++icha ) datamap[icha].metadata.clear();
      assert( ! datamap[0].hasMetadata("fitPedPeakBinFraction") );
      assert( padvmap->updateMap(datamap) == 0 );
      string mname = "fitPedPeakBinFraction";
      for ( AdcIndex icha=0; icha<ncha; ++icha ) {
        cout << myname << "Checking channel " << icha << endl;
        cout << myname << "  Metadata size: " << datamap[icha].metadata.size() << endl;
        cout << myname << "  " << mname << " = " << datamap[icha].metadata[mname] << endl;
        assert( datamap[icha].hasMetadata(mname) );
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
  return test_AdcPedestalFitter(useExistingFcl, doUpdate, doUpdateMap);
}

//**********************************************************************
