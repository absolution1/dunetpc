// test_ConfigurableAdcChannelDataCopyService.cxx
//
// David Adams
// August 2016
//
// Test ConfigurableAdcChannelDataCopyService.

#include <string>
#include <iostream>
#include <sstream>
#include <fstream>
#include <iomanip>
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "dune/ArtSupport/ArtServiceHelper.h"
#include "dune/DuneInterface/AdcChannelDataCopyService.h"

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

//**********************************************************************

int test_ConfigurableAdcChannelDataCopyService(int a_LogLevel =1, int a_MaxConsecutiveFlag =0) {
  const string myname = "test_ConfigurableAdcChannelDataCopyService: ";
#ifdef NDEBUG
  cout << myname << "NDEBUG must be off." << endl;
  abort();
#endif
  string line = "-----------------------------";

  cout << myname << line << endl;
  cout << myname << "Create top-level FCL." << endl;
  string fclfile = "test_ConfigurableAdcChannelDataCopyService.fcl";
  ofstream fout(fclfile.c_str());
  fout << "services.AdcChannelDataCopyService: {" << endl;
  fout << "  service_provider: ConfigurableAdcChannelDataCopyService" << endl;
  fout << "  LogLevel:              " << a_LogLevel << endl;
  fout << "  CopyChannel: true" << endl;
  fout << "  CopyPedestal: true" << endl;
  fout << "  CopySamples: true" << endl;
  fout << "}" << endl;
  fout.close();

  cout << myname << "Fetch art service helper." << endl;
  ArtServiceHelper& ash = ArtServiceHelper::instance();
  ash.print();

  if ( ash.serviceStatus() == 0 ) {
    cout << myname << line << endl;
    cout << myname << "Add ADC mitigation service." << endl;
    assert( ash.addService("AdcChannelDataCopyService", fclfile, true) == 0 );
    ash.print();

    cout << myname << line << endl;
    cout << myname << "Load services." << endl;
    assert( ash.loadServices() == 1 );
    ash.print();
  }

  AdcChannel channel = 123;
  AdcSignal pedestal = 123.45;
  const unsigned int nsig = 50;
  AdcChannelData acd1;
  for ( unsigned int isig=0; isig<50; ++isig ) {
    acd1.samples.push_back(10+isig);
    acd1.flags.push_back(AdcGood);
  }
  AdcSignalVector& sigs = acd1.samples;
  acd1.channel = channel;
  acd1.pedestal = pedestal;

  cout << myname << line << endl;
  cout << myname << "Check original." << endl;
  assert( acd1.channel == channel );
  assert( acd1.pedestal == pedestal );
  assert( acd1.raw.size() == 0);
  assert( acd1.samples == sigs );
  assert( acd1.flags.size() == nsig );
  assert( acd1.raw.size() == 0 );
  assert( acd1.signal.size() == 0 );
  assert( acd1.rois.size() == 0 );
  assert( acd1.digit == nullptr );
  assert( acd1.wire == nullptr );
  assert( acd1.digitIndex == AdcChannelData::badIndex );
  assert( acd1.wireIndex == AdcChannelData::badIndex );

  cout << myname << line << endl;
  cout << myname << "Fetch service." << endl;
  ServiceHandle<AdcChannelDataCopyService> hcop;

  cout << myname << line << endl;
  cout << myname << "Copy." << endl;
  AdcChannelData acd2;
  hcop->copy(acd1, acd2);

  cout << myname << line << endl;
  cout << myname << "Check copy." << endl;
  assert( acd2.channel == channel );
  assert( acd2.pedestal == pedestal );
  assert( acd2.raw.size() == 0);
  assert( acd2.samples == sigs );
  assert( acd2.flags.size() == 0 );
  assert( acd2.raw.size() == 0 );
  assert( acd2.signal.size() == 0 );
  assert( acd2.rois.size() == 0 );
  assert( acd2.digit == nullptr );
  assert( acd2.wire == nullptr );
  assert( acd2.digitIndex == AdcChannelData::badIndex );
  assert( acd2.wireIndex == AdcChannelData::badIndex );

  cout << myname << line << endl;
  cout << myname << "Done." << endl;
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
  return test_ConfigurableAdcChannelDataCopyService(a_LogLevel, a_MaxConsecutiveFlag);
}

//**********************************************************************
