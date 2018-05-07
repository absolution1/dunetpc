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
#include "dune/DuneInterface/Tool/AdcChannelTool.h"
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

#include "TestDigit.h"

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
  auto prdr = tm.getPrivate<AdcChannelTool>("mytool");
  assert( prdr != nullptr );

  cout << myname << line << endl;
  cout << myname << "Construct test digit." << endl;
  TestDigit tdig(0);
  AdcIndex nsig = tdig.nsig;
  const raw::RawDigit* pdig = tdig.pdig;

  cout << myname << line << endl;
  cout << myname << "Extract data from digit." << endl;
  AdcChannelData acd;
  acd.digit = pdig;
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
  assert( chanout == pdig->Channel() );
  assert( pedout == pdig->GetPedestal() );
  assert( acd.raw.size() == pdig->Samples() );
  int dbg = 0;
  for ( unsigned int isig=0; isig<nsig; ++isig ) {
    AdcCount sigin = tdig.sigsin[isig] + tdig.pedestal + 0.1;
    AdcCount sigout = acd.raw[isig];
    AdcFlag flag = tdig.expflags[isig];
    if ( dbg ) cout << isig << ": " << sigin << " --> " << sigout 
                    << " (flag=" << flag << ")" << endl;
    if ( flag == AdcGood) assert( sigout == sigin );
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
  return test_AcdDigitReader(useExistingFcl);
}

//**********************************************************************
