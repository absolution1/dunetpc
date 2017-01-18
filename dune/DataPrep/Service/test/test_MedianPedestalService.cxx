// test_MedianPedestalService.cxx
//
// David Adams
// May 2016
//
// Test MedianPedestalService.

#include <string>
#include <iostream>
#include <sstream>
#include <fstream>
#include <iomanip>
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "dune/ArtSupport/ArtServiceHelper.h"
#include "dune/DuneInterface/PedestalEvaluationService.h"

#undef NDEBUG
#include <cassert>

using std::string;
using std::vector;
using std::cout;
using std::endl;
using std::istringstream;
using std::ofstream;
using std::setw;
using std::setprecision;
using std::fixed;
using art::ServiceHandle;

//**********************************************************************

int test_MedianPedestalService(bool useExistingFcl) {
  const string myname = "test_MedianPedestalService: ";
#ifdef NDEBUG
  cout << myname << "NDEBUG must be off." << endl;
  abort();
#endif
  string line = "-----------------------------";

  cout << myname << line << endl;
  string fclfile = "test_MedianPedestalService.fcl";
  if ( ! useExistingFcl ) {
    cout << myname << "Creating top-level FCL." << endl;
    ofstream fout(fclfile.c_str());
    fout << "services.PedestalEvaluationService: {" << endl;
    fout << "  service_provider: MedianPedestalService" << endl;
    fout << "            LogLevel:     1" << endl;
    fout << "             UseMean: false" << endl;
    fout << "  SkipFlaggedSamples:  true" << endl;
    fout << "         SkipSignals:  true" << endl;
    fout << "}" << endl;
    fout.close();
  } else {
    cout << myname << "Using existing top-level FCL." << endl;
  }

  cout << myname << "Fetch art service helper." << endl;
  ArtServiceHelper& ash = ArtServiceHelper::instance();
  ash.print();

  if ( ash.serviceStatus() == 0 ) {

    cout << myname << line << endl;
    cout << myname << "Add pedestal evaluation service." << endl;
    assert( ash.addService("PedestalEvaluationService", fclfile, true) == 0 );

    cout << myname << line << endl;
    cout << myname << "Load services." << endl;
    assert( ash.loadServices() == 1 );
    ash.print();
  }

  cout << myname << line << endl;
  cout << myname << "Make data." << endl;
  AdcSignalVector oddsigs = { 1.0, 8.0, 2.0, 3.0, 4.0 };
  AdcSignalVector evnsigs = { 1.0, 8.0, 2.0, 4.0, 6.0, 9.0 };
  AdcSignalVector flgsigs = { 1.0, 8.0, 111.0, 2.0, 4.0, 6.0, 222.0, 9.0 };
  AdcFlagVector flags(8, AdcGood);
  flags[2] = AdcStuckOff;
  flags[6] = AdcOverflow;
  AdcChannelData oddData;
  oddData.samples = oddsigs;
  AdcChannelData evnData;
  oddData.samples = evnsigs;
  AdcChannelData flgData;
  flgData.samples = flgsigs;
  flgData.flags = flags;

  cout << myname << line << endl;
  cout << myname << "Fetch service." << endl;
  ServiceHandle<PedestalEvaluationService> hpev;
  hpev->print();

  cout << myname << line << endl;
  cout << myname << "Checking odd data." << endl;
  AdcSignal ped = 0.0;
  assert( hpev->evaluate(oddData, &ped) == 0 );
  assert( ped = 3.0 );

  cout << myname << line << endl;
  cout << myname << "Checking even data." << endl;
  ped = 0.0;
  assert( hpev->evaluate(evnData, &ped) == 0 );
  assert( ped = 3.0 );

  cout << myname << line << endl;
  cout << myname << "Checking flagged data." << endl;
  ped = 0.0;
  assert( hpev->evaluate(flgData, &ped) == 0 );
  assert( ped = 3.0 );

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
  return test_MedianPedestalService(useExistingFcl);
}

//**********************************************************************
