// test_MultiChannelNoiseRemovalService.cxx
//
// David Adams
// May 2016
//
// Test MultiChannelNoiseRemovalService.
//

#include <string>
#include <iostream>
#include <fstream>
#include <iomanip>
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "dune/ArtSupport/ArtServiceHelper.h"
#include "dune/DuneInterface/AdcNoiseRemovalService.h"

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

int test_MultiChannelNoiseRemovalService(int a_LogLevel =1) {
  const string myname = "test_MultiChannelNoiseRemovalService: ";
#ifdef NDEBUG
  cout << myname << "NDEBUG must be off." << endl;
  abort();
#endif
  string line = "-----------------------------";

  cout << myname << line << endl;
  cout << myname << "Create top-level FCL." << endl;
  string fclfile = "test_MultiChannelNoiseRemovalService.fcl";
  ofstream fout(fclfile.c_str());
  AdcSignal threshold = 4.0;
  fout << "services.AdcChannelNoiseRemovalService: {" << endl;
  fout << "  service_provider: ThresholdNoiseRemovalService" << endl;
  fout << "  Threshold:    " << threshold <<  endl;
  fout << "  LogLevel:       " << a_LogLevel << endl;
  fout << "}" << endl;
  fout << "services.AdcNoiseRemovalService: {" << endl;
  fout << "  service_provider: MultiChannelNoiseRemovalService" << endl;
  fout << "  LogLevel:       " << a_LogLevel << endl;
  fout << "}" << endl;
  fout.close();

  cout << myname << "Fetch art service helper." << endl;
  ArtServiceHelper& ash = ArtServiceHelper::instance();
  ash.print();

  if ( ash.serviceStatus() == 0 ) {
    cout << myname << line << endl;
    cout << myname << "Add noise removal service." << endl;
    assert( ash.addService("AdcChannelNoiseRemovalService", fclfile, true) == 0 );
    assert( ash.addService("AdcNoiseRemovalService", fclfile, true) == 0 );
    ash.print();

    cout << myname << line << endl;
    cout << myname << "Load services." << endl;
    assert( ash.loadServices() == 1 );
    ash.print();
  }

  cout << myname << line << endl;
  cout << myname << "Create samples." << endl;
  const unsigned int nchan = 5;
  double f[nchan] = {0.11, 0.22, 0.33, 0.44, 0.55};
  const unsigned int nsig = 100;
  AdcChannelDataMap acdmap;
  vector<AdcSignalVector> sigsin(nchan);     // Record the stae before noise removal.
  vector<AdcSignalVector> sigsexp(nchan);    // Record the expected state.
  for ( unsigned int icha=0; icha<nchan; ++icha ) {
    AdcChannel chan = 100 + icha;
    acdmap[chan] = AdcChannelData();
    AdcChannelData& acd = acdmap[chan];
    AdcSignalVector& exp = sigsexp[icha];
    AdcSignalVector& in = sigsin[icha];
    acd.channel = chan;
    for ( unsigned int isig=0; isig<nsig; ++isig ) {
      acd.samples.push_back(20*sin(f[icha]*isig));
      acd.flags.push_back(AdcGood);
    }
    AdcSignalVector& sigs = acd.samples;
    AdcSignalVector sigsin = sigs;
    in = sigs;
    exp = sigs;
    for ( AdcSignal& sig : exp ) if ( fabs(sig) < threshold ) sig = 0.0;
  }

  cout << myname << line << endl;
  cout << myname << "Fetch the noise removal." << endl;
  ServiceHandle<AdcNoiseRemovalService> hanr;
  hanr->print(cout, myname);

  cout << myname << line << endl;
  cout << myname << "Remove noise." << endl;
  assert( hanr->update(acdmap) == 0 );
  for ( unsigned int icha=0; icha<nchan; ++icha ) {
    AdcChannel chan = 100 + icha;
    AdcChannelData& acd = acdmap[chan];
    AdcSignalVector& sigs = acd.samples;
    cout << myname << "-------" << endl;
    cout << myname << "Channel " << chan << endl;
    for ( unsigned int isig=0; isig<nsig; ++isig ) {
      cout << myname << setw(4) << isig << ": "
           << fixed << setprecision(1) << setw(10) << sigsin[icha][isig]
           << fixed << setprecision(1) << setw(10) << sigs[isig]
           << endl;
      assert( sigs[isig] == sigsexp[icha][isig] );
    }
  }
  
  cout << myname << line << endl;
  cout << myname << "Done." << endl;
  return 0;
}

//**********************************************************************

int main(int argc, char* argv[]) {
  int a_LogLevel = -1;
  if ( argc > 1 ) {
    istringstream ssarg(argv[1]);
    ssarg >> a_LogLevel;
  }
  return test_MultiChannelNoiseRemovalService(a_LogLevel);
}

//**********************************************************************
