// test_StandardRawDigitExtractService.cxx
//
// David Adams
// May 2016
//
// Test StandardRawDigitExtractService.

#include <string>
#include <iostream>
#include <sstream>
#include <fstream>
#include <iomanip>
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "lardataobj/RawData/RawDigit.h"
#include "dune/ArtSupport/ArtServiceHelper.h"
#include "dune/DuneInterface/AdcTypes.h"
#include "dune/DuneInterface/RawDigitExtractService.h"

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

bool sigequal(AdcSignal sig1, AdcSignal sig2) {
  AdcSignal sigdiff = sig2 - sig1;
  if ( sigdiff < -0.5 ) return false;
  if ( sigdiff >  0.5 ) return false;
  return true;
}

//**********************************************************************

int test_StandardRawDigitExtractService() {
  const string myname = "test_StandardRawDigitExtractService: ";
#ifdef NDEBUG
  cout << myname << "NDEBUG must be off." << endl;
  abort();
#endif
  string line = "-----------------------------";

  cout << myname << line << endl;
  cout << myname << "Create top-level FCL." << endl;
  string fclfile = "test_StandardRawDigitExtractService.fcl";
  ofstream fout(fclfile.c_str());
  int loglevel = 1;
  fout << "#include \"dataprep_tools.fcl\"" << endl;
  fout << "services.RawDigitExtractService: {" << endl;
  fout << "  service_provider: StandardRawDigitExtractService" << endl;
  fout << "  LogLevel: " << loglevel << endl;
  fout << "  DigitReadTool: digitReader" << endl;
  fout << "  PedestalOption: 1" << endl;
  fout << "  FlagStuckOff: true" << endl;
  fout << "  FlagStuckOn: true" << endl;
  fout << "}" << endl;
  fout.close();

  cout << myname << "Fetch art service helper." << endl;
  ArtServiceHelper& ash = ArtServiceHelper::instance();
  ash.print();

  cout << myname << line << endl;
  cout << myname << "Add raw digit extract service." << endl;
  assert( ash.addService("RawDigitExtractService", fclfile, true) == 0 );
  ash.print();

  cout << myname << line << endl;
  cout << myname << "Load services." << endl;
  assert( ash.loadServices() == 1 );
  ash.print();

  AdcCount lowbits = 63;
  AdcCount highbits = 4095 - 63;

  cout << myname << line << endl;
  cout << myname << "Create a raw digit." << endl;
  AdcSignalVector sigsin;
  float fac = 250.0;
  float ped = 2000.2;
  for ( int i=0; i<10; ++i ) sigsin.push_back(fac*i);
  for ( int i=10; i>=0; --i ) sigsin.push_back(fac*i);
  for ( int i=19; i>=0; --i ) sigsin.push_back(-sigsin[i]);
  unsigned int nsig = 41;
  assert(sigsin.size() == nsig);
  AdcCountVector adcsin;
  unsigned int isig_stucklo = 5;
  unsigned int isig_stuckhi = 15;
  for ( unsigned int isig=0; isig<nsig; ++isig) {
    AdcSignal sig = sigsin[isig] + ped;
    AdcCount adc = 0.0;
    if ( sig > 0.0 ) adc = int(sig+0.5);
    if ( adc > 4095 ) adc = 4095;
    AdcCount adchigh = adc & highbits;
    if ( isig == isig_stucklo ) adc = adchigh;           // Stuck low bits to zero.
    if ( isig == isig_stuckhi ) adc = adchigh + lowbits; // Stuck low bits to one.
    adcsin.push_back(adc);
  }
  assert(adcsin.size() == nsig);
  unsigned int chan = 123;
  raw::RawDigit dig(chan, nsig, adcsin, raw::kNone);
  dig.SetPedestal(ped);
  cout << myname << "    Compressed size: " << dig.NADC() << endl;
  cout << myname << "  Uncompressed size: " << dig.Samples() << endl;
  cout << myname << "           Pedestal: " << dig.GetPedestal() << endl;
  cout << myname << "            Channel: " << dig.Channel() << endl;
  assert(dig.Samples() == nsig);
  assert(dig.Channel() == chan);
  assert(dig.GetPedestal() == ped);

  cout << myname << line << endl;
  cout << myname << "Create the expected flag vector." << endl;
  AdcFlagVector expflags(nsig, AdcGood);
  for ( unsigned int isig=0; isig<nsig; ++isig) {
    AdcCount adc = adcsin[isig];
    AdcCount adclow = adc & lowbits;
    if ( adc <= 0 ) expflags[isig] = AdcUnderflow;
    else if ( adc >= 4095 ) expflags[isig] = AdcOverflow;
    else if ( adclow == 0 ) expflags[isig] = AdcStuckOff;
    else if ( adclow == lowbits ) expflags[isig] = AdcStuckOn;
  }

  cout << myname << line << endl;
  cout << myname << "Fetch raw digit extract service." << endl;
  ServiceHandle<RawDigitExtractService> hrdx;
  hrdx->print();
  ash.print();

  cout << myname << line << endl;
  cout << myname << "Extract data from digit." << endl;
  AdcChannelData acd;
  acd.digit = &dig;
  assert( hrdx->extract(acd) == 0 );
  const AdcCountVector& raw = acd.raw;
  const AdcSignalVector& sigs = acd.samples;
  const AdcFlagVector& flags = acd.flags;
  AdcChannel chanout = acd.channel;
  AdcSignal pedout = acd.pedestal;
  cout << myname << "Output raw vector size: " << sigs.size() << endl;
  cout << myname << "Output prep vector size: " << sigs.size() << endl;
  cout << myname << " Output flags size: " << flags.size() << endl;
  cout << myname << "          Pedestal: " << pedout << endl;
  assert( raw.size() == nsig );
  assert( sigs.size() == nsig );
  assert( flags.size() == nsig );
  assert( chanout == chan );
  assert( pedout == ped );
  for ( unsigned int isig=0; isig<nsig; ++isig ) {
    cout << setw(4) << isig << ": " << setw(4) << adcsin[isig]
         << fixed << setprecision(1) << setw(8) << sigs[isig]
         << " [" << flags[isig] << "]" << endl;
    if ( flags[isig] == AdcGood ) assert( sigequal(sigs[isig], sigsin[isig]) );
    assert( flags[isig] == expflags[isig] );
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
  return test_StandardRawDigitExtractService();
}

//**********************************************************************
