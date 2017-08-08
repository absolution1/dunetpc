// test_AcdDigitReader.cxx
//
// David Adams
// April 2017
//
// Test AcdDigitReader.

#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include "dune/DuneInterface/Tool/AdcChannelDataModifier.h"
#include "dune/ArtSupport/DuneToolManager.h"
#include "lardataobj/RawData/RawDigit.h"

#undef NDEBUG
#include <cassert>

using std::string;
using std::cout;
using std::endl;
using std::ostringstream;
using std::ofstream;
using std::setw;
using std::setprecision;
using std::fixed;
using fhicl::ParameterSet;

//**********************************************************************

bool sigequal(AdcSignal sig1, AdcSignal sig2) {
  AdcSignal sigdiff = sig2 - sig1;
  if ( sigdiff < -0.5 ) return false;
  if ( sigdiff >  0.5 ) return false;
  return true;
}

//**********************************************************************

int test_AcdDigitReader(bool useExistingFcl =false) {
  const string myname = "test_AcdDigitReader: ";
#ifdef NDEBUG
  cout << myname << "NDEBUG must be off." << endl;
  abort();
#endif
  string line = "-----------------------------";

  cout << myname << line << endl;
  string fclfile = "test_AcdDigitReader.fcl";
  if ( ! useExistingFcl ) {
    cout << myname << "Creating top-level FCL." << endl;
    ofstream fout(fclfile.c_str());
    fout << "tools: {" << endl;
    fout << "  mytool: {" << endl;
    fout << "    tool_type: AcdDigitReader" << endl;
    fout << "    LogLevel: 2" << endl;
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
  auto prdr = tm.getPrivate<AdcChannelDataModifier>("mytool");
  assert( prdr != nullptr );

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
  cout << myname << "Extract data from digit." << endl;
  AdcChannelData acd;
  acd.digit = &dig;
  assert( prdr->update(acd) == 0 );
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
  //assert( sigs.size() == nsig );
  //assert( flags.size() == nsig );
  assert( chanout == chan );
  assert( pedout == ped );
  //for ( unsigned int isig=0; isig<nsig; ++isig ) {
  //  cout << setw(4) << isig << ": " << setw(4) << adcsin[isig]
  //       << fixed << setprecision(1) << setw(8) << sigs[isig]
  //       << " [" << flags[isig] << "]" << endl;
  //  if ( flags[isig] == AdcGood ) assert( sigequal(sigs[isig], sigsin[isig]) );
  //  assert( flags[isig] == expflags[isig] );
  //}

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
  return test_AcdDigitReader(useExistingFcl);
}

//**********************************************************************
