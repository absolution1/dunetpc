// test_FclStickyCodeFlagger.cxx
//
// David Adams
// April 2017
//
// Test FclStickyCodeFlagger.

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
using Index = unsigned int;

//**********************************************************************

int test_FclStickyCodeFlagger(bool useExistingFcl =false) {
  const string myname = "test_FclStickyCodeFlagger: ";
#ifdef NDEBUG
  cout << myname << "NDEBUG must be off." << endl;
  abort();
#endif
  string line = "-----------------------------";

  cout << myname << line << endl;
  string fclfile = "test_FclStickyCodeFlagger.fcl";
  if ( ! useExistingFcl ) {
    cout << myname << "Creating top-level FCL." << endl;
    ofstream fout(fclfile.c_str());
    fout << "tools: {" << endl;
    fout << "  mytool: {" << endl;
    fout << "    tool_type: FclStickyCodeFlagger" << endl;
    fout << "    LogLevel: 3" << endl;
    fout << "    StickyCode: 8" << endl;
    fout << "    StickyCodes: {" << endl;
    fout << "      chan012: [105, 110]" << endl;
    fout << "      chan012x: [109]" << endl;
    fout << "    }" << endl;
    fout << "    StickyRanges: {" << endl;
    fout << "      chan014: [110, 111]" << endl;
    fout << "      chan014x: [119, 119]" << endl;
    fout << "    }" << endl;
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
  cout << myname << "Create data and call tool." << endl;
  AdcIndex nevt = 2;
  for ( AdcIndex ievt=0; ievt<nevt; ++ievt ) {
    cout << myname << "Event " << ievt << endl;
    AdcChannelDataMap datamap;
    AdcIndex ncha = 20;
    for ( AdcIndex icha=0; icha<ncha; ++icha ) {
      std::pair<AdcChannelDataMap::iterator, bool> kdat = datamap.emplace(icha, AdcChannelData());
      assert(kdat.second);
      AdcChannelDataMap::iterator idat = kdat.first;
      AdcChannelData& data = idat->second;
      float ped = 100.0;
      data.setEventInfo(new AdcChannelData::EventInfo(101, ievt, 23));
      data.channel = icha;
      data.pedestal = ped;
      for ( AdcIndex itic=0; itic<20; ++itic ) {
        AdcCount iadc = ped + itic;
        data.raw.push_back(iadc);
        data.flags.push_back(0);
      }
      DataMap dm = padv->update(datamap[icha]);
      assert( dm == 0 );
      dm.print();
      if ( data.channel == 12 ) {
        for ( Index isam=0; isam<data.raw.size(); ++isam ) {
          cout << myname << "    " << data.raw[isam] << "  " << data.flags[isam] << endl;
          if ( data.raw[isam] == 105 || data.raw[isam] == 109 || data.raw[isam] == 110 ) assert( data.flags[isam] == 8 );
          else assert( data.flags[isam] == 0 );
        }
      }
      if ( data.channel == 14 ) {
        for ( Index isam=0; isam<data.raw.size(); ++isam ) {
          cout << myname << "    " << data.raw[isam] << "  " << data.flags[isam] << endl;
          if ( data.raw[isam] == 110 || data.raw[isam] == 111 || data.raw[isam] == 119 ) assert( data.flags[isam] == 8 );
          else assert( data.flags[isam] == 0 );
        }
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
  if ( argc > 1 ) {
    string sarg(argv[1]);
    if ( sarg == "-h" ) {
      cout << "Usage: " << argv[0] << " [ARG]" << endl;
      cout << "  If ARG = true, existing FCL file is used." << endl;
      return 0;
    }
    useExistingFcl = sarg == "true" || sarg == "1";
  }
  return test_FclStickyCodeFlagger(useExistingFcl);
}

//**********************************************************************
