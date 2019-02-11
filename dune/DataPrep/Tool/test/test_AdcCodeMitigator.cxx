// test_AdcCodeMitigator.cxx
//
// David Adams
// April 2017
//
// Test AdcCodeMitigator.

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
using std::setw;
using std::ofstream;
using fhicl::ParameterSet;
using Index = unsigned int;

bool isEqual(float x1, float x2) {
  if ( fabs(x2 - x1) < 1.e-5 ) return true;
  cout << "isEqual: " << x1 << " != " << x2 << endl;
  return false;
}

//**********************************************************************

int test_AdcCodeMitigator(bool useExistingFcl =false, bool interp =false) {
  const string myname = "test_AdcCodeMitigator: ";
#ifdef NDEBUG
  cout << myname << "NDEBUG must be off." << endl;
  abort();
#endif
  string line = "-----------------------------";

  cout << myname << line << endl;
  string fclfile = "test_AdcCodeMitigator.fcl";
  if ( ! useExistingFcl ) {
    cout << myname << "Creating top-level FCL." << endl;
    ofstream fout(fclfile.c_str());
    fout << "tools: {" << endl;
    fout << "  mytool: {" << endl;
    fout << "    tool_type: AdcCodeMitigator" << endl;
    fout << "    LogLevel: 1" << endl;
    if ( interp ) {
      fout << "    FixFlags: []" << endl;
      fout << "    InterpolateFlags: [8]" << endl;
    } else {
      fout << "    FixFlags: [8]" << endl;
      fout << "    InterpolateFlags: []" << endl;
    }
    fout << "    SkipFlags: []" << endl;
    fout << "  FixedCurvThresh: 0.0" << endl;
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
  auto pmit = tm.getPrivate<AdcChannelTool>("mytool");
  assert( pmit != nullptr );

  cout << myname << line << endl;
  cout << myname << "Create data and call tool." << endl;
  AdcSignalVector chkSigs = {   5,  5, 11.0, 11.5, 12.0, 12.5, 13, 14,   14};
  AdcSignalVector inpSigs = {1000,  5, 11.0, 1000, 1000, 1000, 13, 14, 1000};
  AdcFlagVector   inpFlgs = {   8,  0,    0,    8,    8,    8,  0,  0,    8};
  AdcFlagVector   chkFlgs = {  AdcExtrapolated,
                                    0,    0,   AdcInterpolated, AdcInterpolated, AdcInterpolated,
                                                                0,  0,   AdcExtrapolated};
  Index nsam = chkSigs.size();
  if ( ! interp ) {
    for ( AdcIndex isam=0; isam<nsam; ++isam ) {
      if ( inpFlgs[isam] != 0 ) {
        chkSigs[isam] = 0.0;
        chkFlgs[isam] = AdcSetFixed;
      }
    }
  }
  AdcChannelData acd;
  acd.samples = inpSigs;
  acd.flags = inpFlgs;
  DataMap ret = pmit->update(acd); 
  ret.print();
  assert( ret == 0 );
  assert( ret.getInt("mitCount") == 5 );
  for ( AdcIndex isam=0; isam<nsam; ++isam ) {
    cout << isam << ":" << setw(6) << inpSigs[isam] << setw(6) << inpFlgs[isam]
         << setw(6) << acd.samples[isam] << setw(6) << acd.flags[isam] << endl;
    assert( isEqual(acd.samples[isam], chkSigs[isam]) );
    assert( acd.flags[isam] == chkFlgs[isam] );
  }

  cout << myname << line << endl;
  cout << myname << "Done." << endl;
  return 0;
}

//**********************************************************************

int main(int argc, char* argv[]) {
  bool useExistingFcl = false;
  bool interp = true;
  if ( argc > 1 ) {
    string sarg(argv[1]);
    if ( sarg == "-h" ) {
      cout << "Usage: " << argv[0] << " [ARG] [FIX]" << endl;
      cout << "  If ARG = true, existing FCL file is used." << endl;
      cout << "  If FIX = true, fixed value is used instead of interpolation." << endl;
      return 0;
    }
    useExistingFcl = sarg == "true" || sarg == "1";
    if ( argc > 2 ) {
      interp = false;
      string sarg(argv[2]);
      interp = sarg == "true" || sarg == "1";
    }
  }
  return test_AdcCodeMitigator(useExistingFcl, interp);
}

//**********************************************************************
