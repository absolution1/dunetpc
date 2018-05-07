// test_AcdWireReader.cxx
//
// David Adams
// April 2017
//
// Test AcdWireReader.

#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "dune/DuneInterface/RawDigitPrepService.h"
#include "dune/DuneInterface/Tool/AdcChannelTool.h"
#include "dune/ArtSupport/ArtServiceHelper.h"
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
using art::ServiceHandle;

#include "TestDigit.h"

//**********************************************************************

bool sigequal(AdcSignal sig1, AdcSignal sig2) {
  AdcSignal sigdiff = sig2 - sig1;
  if ( sigdiff < -0.5 ) return false;
  if ( sigdiff >  0.5 ) return false;
  return true;
}

//**********************************************************************

int test_AcdWireReader(bool useExistingFcl =false) {
  const string myname = "test_AcdWireReader: ";
#ifdef NDEBUG
  cout << myname << "NDEBUG must be off." << endl;
  abort();
#endif
  string line = "-----------------------------";

  cout << myname << line << endl;
  string fclfile = "test_AcdWireReader.fcl";
  if ( ! useExistingFcl ) {
    cout << myname << "Creating top-level FCL." << endl;
    // To convert digit to wire, use the DUNE fcl for 35-ton reco.
    // Disable noise removal and deconvolution because these make it difficult
    // to predict the result.
    ofstream fout(fclfile.c_str());
    fout << "#include \"services_dune.fcl\"" << endl;
    fout << "services: @local::dune35tdata_reco_services" << endl;
    fout << "services.RawDigitPrepService.DoNoiseRemoval: false" << endl;
    fout << "services.RawDigitPrepService.DoDeconvolution: false" << endl;
    fout << "services.RawDigitPrepService.DoIntermediateStates: true" << endl;
    fout << "services.AdcChannelDataCopyService.CopyFlags: true" << endl;
    // Need the standard tool to read channel data from a digit.
    fout << "#include \"dataprep_tools.fcl\"" << endl;
    // Build local wire reader.
    fout << "tools.mytool: {" << endl;
    fout << "  tool_type: AcdWireReader" << endl;
    fout << "  LogLevel: 2" << endl;
    fout << "}" << endl;
    fout.close();
  } else {
    cout << myname << "Using existing top-level FCL." << endl;
  }

  cout << myname << line << endl;
  cout << myname << "Fetch art service helper." << endl;
  ArtServiceHelper& ash = ArtServiceHelper::instance();
  assert( ash.addServices(fclfile, true) == 0 );
  assert( ash.loadServices() == 1 );
  ash.print();

  cout << myname << line << endl;
  cout << myname << "Fetching tool manager." << endl;
  DuneToolManager* ptm = DuneToolManager::instance(fclfile);
  assert ( ptm != nullptr );
  DuneToolManager& tm = *ptm;
  tm.print();
  assert( tm.toolNames().size() );

  cout << myname << line << endl;
  cout << myname << "Fetch services." << endl;
  ServiceHandle<RawDigitPrepService> hrdp;

  cout << myname << line << endl;
  cout << myname << "Construct test digit." << endl;
  TestDigit tdig(2);
  AdcIndex nsig = tdig.nsig;
  const raw::RawDigit* pdig = tdig.pdig;

  cout << myname << line << endl;
  cout << myname << "Create channel dat map and set digits." << endl;
  AdcChannelDataMap acds;
  AdcChannelData& acd = acds[tdig.channel];
  acd.digit = pdig;
  assert( acd.raw.size() == 0 );

  cout << myname << line << endl;
  cout << myname << "Construct test wire." << endl;
  std::vector<recob::Wire> wires;
  wires.reserve(acds.size());
  assert( acd.samples.size() == 0 );
  assert( hrdp->prepare(acds, &wires) == 0 );
  assert( acd.raw.size() == nsig );
  assert( acd.samples.size() == nsig );

  cout << myname << line << endl;
  cout << myname << "Fetch wire read tool." << endl;
  auto prdr = tm.getPrivate<AdcChannelTool>("mytool");
  assert( prdr != nullptr );

  cout << myname << line << endl;
  cout << myname << "Read wire into channel data." << endl;
  AdcChannelData newacd;
  newacd.wire = &wires[0];
  assert( newacd.samples.size() == 0 );
  assert( prdr->update(newacd) == 0 );
  assert( newacd.samples.size() == nsig );

  cout << myname << line << endl;
  cout << myname << "Check samples." << endl;
  int dbg = 1;
  for ( unsigned int isig=0; isig<nsig; ++isig ) {
    if ( dbg ) cout << myname << isig << ": " << acd.samples[isig]
                    << " ?= " << newacd.samples[isig]
                    << " (" << acd.signal[isig] << " "
                    << newacd.signal[isig] << ")" << endl;
    assert(sigequal(acd.samples[isig], newacd.samples[isig]));
  }
  cout << myname << "NROI: " << wires[0].SignalROI().size() << " " 
                 << acd.rois.size() << " " << newacd.rois.size() << endl;
  assert( newacd.rois.size() == acd.rois.size() );

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
  return test_AcdWireReader(useExistingFcl);
}

//**********************************************************************
