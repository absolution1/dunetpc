// test_FembLinearCalibration.cxx
//
// David Adams
// November 2017
//
// Test FembLinearCalibration.

#include <string>
#include <iostream>
#include <fstream>
#include <iomanip>
#include "dune/DuneInterface/Tool/AdcChannelTool.h"
#include "dune/ArtSupport/DuneToolManager.h"

#undef NDEBUG
#include <cassert>

using std::string;
using std::cout;
using std::endl;
using std::ofstream;
using std::setw;
using fhicl::ParameterSet;

using Index = unsigned int;

//**********************************************************************

int test_FembLinearCalibration(bool useExistingFcl =false) {
  const string myname = "test_FembLinearCalibration: ";
#ifdef NDEBUG
  cout << myname << "NDEBUG must be off." << endl;
  abort();
#endif
  string line = "-----------------------------";

  cout << myname << line << endl;
  string fclfile = "test_FembLinearCalibration.fcl";
  if ( ! useExistingFcl ) {
    cout << myname << "Creating top-level FCL." << endl;
    ofstream fout(fclfile.c_str());
    fout << "tools: {" << endl;
    fout << "  mytool: {" << endl;
    fout << "    tool_type: FembLinearCalibration" << endl;
    fout << "    LogLevel: 1" << endl;
    fout << "    Units: Coulombs" << endl;
    fout << "    Gains: [1.0, 2.0, 3.0, 4.0, 5.0]" << endl;
    fout << "    AdcMin: 0" << endl;
    fout << "    AdcMins: [1200, 1300, 1400, 1500, 1600]" << endl;
    fout << "    AdcMax: 1800" << endl;
    fout << "    AdcMaxs: []" << endl;
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
  auto pmod = tm.getPrivate<AdcChannelTool>("mytool");
  assert( pmod != nullptr );

  cout << myname << line << endl;
  cout << myname << "Create data." << endl;
  AdcChannelData acd0;
  acd0.pedestal = 1000.0;
  acd0.raw.push_back(1100);
  acd0.raw.push_back(1200);
  acd0.raw.push_back(1300);
  acd0.raw.push_back(1400);
  acd0.raw.push_back(1500);
  acd0.pedestal = 1000.0;
  AdcSignal dped = 0.0;
  AdcChannelDataMap acds;
  for ( AdcChannel icha=0; icha<5; ++icha ) {
    AdcChannelData& acd = acds[icha];
    acd.channel = icha;
    for ( AdcCount& adc : acd0.raw ) acd.raw.push_back(adc+dped);
    acd.pedestal = acd0.pedestal + dped;
    dped += 100.0;
  }
  int w = 8;
  Index nsam = acd0.raw.size();
  for ( auto ient : acds ) {
    AdcChannel icha = ient.first;
    AdcChannelData& acd = ient.second;
    AdcSignalVector sigchk(nsam);
    for ( Index isam=0; isam<nsam; ++isam ) {
      sigchk[isam] = (icha+1)*(acd0.raw[isam] - acd0.pedestal);
    }
    cout << myname << line << endl;
    assert( acd.sampleUnit == "" );
    cout << myname << "Channel " << icha << endl;
    DataMap res = pmod->update(acd);
    cout << myname << "Modify:" << endl;
    res.print();
    Index nsamRes = res.getInt("calibSampleCount");
    Index nunder = res.getInt("calibUnderflowCount");
    Index nover = res.getInt("calibOverflowCount");
    assert( res.status() == 0 );
    assert( nsamRes == nsam );
    cout << myname << "  Raw:";
    for ( AdcCount adc : acd.raw ) cout << setw(w) << adc;
    cout << endl;
    cout << myname << " Prep:";
    for ( AdcSignal sam : acd.samples ) cout << setw(w) << sam;
    cout << endl;
    cout << myname << "Check:";
    for ( AdcSignal sam : sigchk ) cout << setw(w) << sam;
    cout << endl;
    cout << myname << " Flag:";
    for ( AdcFlag flg : acd.flags ) cout << setw(w) << flg;
    cout << endl;
    cout << myname << "Check samples." << endl;
    for ( Index isam=0; isam<nsam; ++isam ) {
      assert( acd.samples[isam] == sigchk[isam] );
    }
    cout << myname << "Check under and overflows." << endl;
    Index nunderExp = 2;
    Index noverExp = icha > 2 ? icha - 2 : 0;
    assert ( nunder == nunderExp );
    assert ( nover == noverExp );
    for ( Index isam=0; isam<nsam; ++isam ) {
      AdcFlag flgExp = isam < nunderExp ? AdcUnderflow :
                       isam >= nsam-noverExp ? AdcOverflow :
                       AdcGood;
      if ( acd.flags[isam] != flgExp ) {
         cout << myname << "Sample " << isam << ": " << acd.flags[isam]
              << " != " << flgExp << endl;
      }
      assert( acd.flags[isam] == flgExp  );
    }
    assert( acd.sampleUnit == "Coulombs" );
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
  return test_FembLinearCalibration(useExistingFcl);
}

//**********************************************************************
