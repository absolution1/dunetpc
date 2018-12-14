// test_AdcSuppressSignalFindingService.cxx
//
// David Adams
// June 2016
//
// Test AdcSuppressSignalFindingService.

#include <string>
#include <iostream>
#include <sstream>
#include <fstream>
#include <iomanip>
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "lardataobj/RawData/RawDigit.h"
#include "dune/ArtSupport/ArtServiceHelper.h"
#include "dune/DuneInterface/AdcTypes.h"
#include "dune/DuneInterface/AdcSignalFindingService.h"

#undef NDEBUG
#include <cassert>

using std::string;
using std::cout;
using std::endl;
using std::istringstream;
using std::ofstream;
using std::setw;
using std::setprecision;
using std::fixed;
using art::ServiceHandle;

//**********************************************************************

int test_AdcSuppressSignalFindingService() {
  const string myname = "test_AdcSuppressSignalFindingService: ";
#ifdef NDEBUG
  cout << myname << "NDEBUG must be off." << endl;
  abort();
#endif
  string line = "-----------------------------";

  cout << myname << line << endl;
  cout << myname << "Create top-level FCL." << endl;
  std::ostringstream oss;
  oss << "services.AdcSuppressService: {" << endl;
  oss << "  service_provider: Legacy35tZeroSuppressService" << endl;
  oss << "  AdcThreshold: 30.0" << endl;
  oss << "  TickRange: 4" << endl;
  oss << "  MinTickGap: 2" << endl;
  oss << "  SuppressStickyBits: false" << endl;
  oss << "}" << endl;
  oss << "services.AdcSignalFindingService: {" << endl;
  oss << "  service_provider: AdcSuppressSignalFindingService" << endl;
  oss << "}" << endl;
  ArtServiceHelper::load_services(oss.str());

  cout << myname << line << endl;
  cout << myname << "Create data." << endl;
  AdcChannelData data;
  data.channel = 1234;
  data.pedestal = 1000.0;
  for ( unsigned int tic=0; tic<50; ++tic ) {
    AdcSignal sig = tic%10 - 5.0;
    if ( tic > 10 && tic < 20 ) sig += 100;
    AdcCount cnt = sig + data.pedestal;
    data.raw.push_back(cnt);
    data.samples.push_back(sig);
    data.flags.push_back(AdcGood);
  }
  unsigned int nsig = data.samples.size();

  cout << myname << line << endl;
  cout << myname << "Fetch signal finding service." << endl;
  ServiceHandle<AdcSignalFindingService> hsfs;
  hsfs->print();

  AdcFilterVector exp(50, false);
  for ( unsigned int isig=7; isig<24; ++isig ) exp[isig] = true;

  cout << myname << line << endl;
  cout << myname << "Find signal." << endl;
  assert( hsfs->find(data) == 0 );
  assert( data.raw.size() == nsig );
  assert( data.samples.size() == nsig );
  assert( data.flags.size() == nsig );
  assert( data.signal.size() == nsig );
  for ( unsigned int isig=0; isig<nsig; ++isig ) {
    cout << setw(4) << isig << ": " << setw(6) << data.raw[isig]
         << fixed << setprecision(1) << setw(8) << data.samples[isig]
         << " " << ( data.signal[isig] ? "*" : " ") << endl;
    assert( data.signal[isig] == exp[isig] );
  }

  return 0;
}

//**********************************************************************

int main(int argc, char* argv[]) {
  int logLevel = 1;
  if ( argc > 1 ) {
    istringstream ssarg(argv[1]);
    ssarg >> logLevel;
  }
  return test_AdcSuppressSignalFindingService();
}

//**********************************************************************
