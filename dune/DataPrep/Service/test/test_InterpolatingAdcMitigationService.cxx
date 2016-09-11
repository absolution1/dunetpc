// test_InterpolatingAdcMitigationService.cxx
//
// David Adams
// May 2016
//
// Test InterpolatingAdcMitigationService.

#include <string>
#include <iostream>
#include <sstream>
#include <fstream>
#include <iomanip>
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "dune/ArtSupport/ArtServiceHelper.h"
#include "dune/DuneInterface/AdcMitigationService.h"

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

typedef vector<unsigned int> IndexVector;

bool sigequal(AdcSignal sig1, AdcSignal sig2) {
  AdcSignal sigdiff = sig2 - sig1;
  if ( sigdiff < -0.5 ) return false;
  if ( sigdiff >  0.5 ) return false;
  return true;
}

//**********************************************************************

int test_InterpolatingAdcMitigationService(int a_LogLevel =1, int a_MaxConsecutiveFlag =0) {
  const string myname = "test_InterpolatingAdcMitigationService: ";
#ifdef NDEBUG
  cout << myname << "NDEBUG must be off." << endl;
  abort();
#endif
  string line = "-----------------------------";

  cout << myname << line << endl;
  cout << myname << "Create top-level FCL." << endl;
  string fclfile = "test_InterpolatingAdcMitigationService.fcl";
  ofstream fout(fclfile.c_str());
  fout << "services.AdcMitigationService: {" << endl;
  fout << "  service_provider: InterpolatingAdcMitigationService" << endl;
  fout << "  LogLevel:              " << a_LogLevel << endl;
  fout << "  SkipUnderflows:     true" << endl;
  fout << "  SkipOverflows:     false" << endl;
  fout << "  MaxConsecutiveSamples: 3 " << endl;
  fout << "  MaxConsecutiveFlag:    " << a_MaxConsecutiveFlag << endl;
  fout << "}" << endl;
  fout.close();

  cout << myname << "Fetch art service helper." << endl;
  ArtServiceHelper& ash = ArtServiceHelper::instance();
  ash.print();

  if ( ash.serviceStatus() == 0 ) {
    cout << myname << line << endl;
    cout << myname << "Add ADC mitigation service." << endl;
    assert( ash.addService("AdcMitigationService", fclfile, true) == 0 );
    ash.print();

    cout << myname << line << endl;
    cout << myname << "Load services." << endl;
    assert( ash.loadServices() == 1 );
    ash.print();
  }

  const unsigned int nsig = 50;
  AdcChannelData acd;
  for ( unsigned int isig=0; isig<50; ++isig ) {
    acd.samples.push_back(10+isig);
    acd.flags.push_back(AdcGood);
  }
  AdcSignalVector& sigs = acd.samples;
  AdcFlagVector& flags = acd.flags;

  // Expected results.
  AdcSignalVector sigsexp = acd.samples;
  AdcFlagVector flagsexp = acd.flags;

  // Add overflows. These are mitigated.
  IndexVector isigOvers = {5, 6, 9};
  for ( unsigned int isig : isigOvers ) {
    sigs[isig] = 9999;
    flags[isig] = AdcOverflow;
    flagsexp[isig] = AdcInterpolated;
  }

  // Add underflows. These are not mitigated.
  IndexVector isigUnders = {15, 16, 19};
  for ( unsigned int isig : isigUnders ) {
    sigs[isig] = -9999;
    flags[isig] = AdcUnderflow;
    sigsexp[isig] = sigs[isig];
    flagsexp[isig] = flags[isig];
  }

  // Add samples with stuck bits on.
  // These are over short ranges so they are mitigated.
  IndexVector isigOns = {21, 23, 24, 28, 29, 30};
  for ( unsigned int isig : isigOns ) {
    sigs[isig] = 1111;
    flags[isig] = AdcStuckOn;
    flagsexp[isig] = AdcInterpolated;
  }

  // Add samples with struck bits off.
  // These are at edges or over long ranges so they are zeroed.
  IndexVector isigOffs = {0, 1, 33, 34, 35, 36, 41, 42, 43, 44, 45, 46, nsig-1};
  for ( unsigned int isig : isigOffs ) {
    sigs[isig] = 1100;
    flags[isig] = AdcStuckOff;
    if ( a_MaxConsecutiveFlag == 1 ) {
      sigsexp[isig] = 0.0;
      flagsexp[isig] = AdcSetFixed;
    } else {
      assert( a_MaxConsecutiveFlag == 0 );
      sigsexp[isig] = sigs[isig];
      flagsexp[isig] = flags[isig];
    }
  }

  // Data before mitigation.
  AdcSignalVector sigsin = acd.samples;
  AdcFlagVector flagsin = acd.flags;

  cout << myname << "Fetch ADC mitigation service." << endl;
  ServiceHandle<AdcMitigationService> hams;
  hams->print();
  ash.print();

  cout << myname << line << endl;
  cout << myname << "Mitigate." << endl;
  assert( hams->update(acd) == 0 );
  cout << myname << "Output vector size: " << sigs.size() << endl;
  cout << myname << " Output flags size: " << flags.size() << endl;
  assert( sigs.size() == nsig );
  assert( flags.size() == nsig );
  for ( unsigned int isig=0; isig<nsig; ++isig ) {
    cout << setw(4) << isig << ": "
         << fixed << setprecision(1) << setw(8) << sigsin[isig]
         << " [" << flagsin[isig] << "]"
         << fixed << setprecision(1) << setw(10) << sigs[isig]
         << " [" << flags[isig] << "]"
         << " [" << flagsexp[isig] << "]"
         << endl;
    assert( flags[isig] == flagsexp[isig] );
    assert( sigs[isig] == sigsexp[isig] );
  }

  return 0;
}

//**********************************************************************

int main(int argc, char* argv[]) {
  int a_LogLevel = 1;
  int a_MaxConsecutiveFlag = 0;
  if ( argc > 1 ) {
    istringstream ssarg(argv[1]);
    ssarg >> a_LogLevel;
  }
  if ( argc > 2 ) {
    istringstream ssarg(argv[1]);
    ssarg >> a_MaxConsecutiveFlag;
  }
  return test_InterpolatingAdcMitigationService(a_LogLevel, a_MaxConsecutiveFlag);
}

//**********************************************************************
