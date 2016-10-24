// test_DuneRoiBuildingService.cxx
//
// David Adams
// May 2016
//
// Test DuneRoiBuildingService.
//

#include <string>
#include <iostream>
#include <sstream>
#include <fstream>
#include <iomanip>
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "dune/ArtSupport/ArtServiceHelper.h"
#include "dune/Utilities/SignalShapingServiceDUNE.h"
#include "dune/DuneInterface/AdcRoiBuildingService.h"

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

int test_DuneRoiBuildingService(int a_LogLevel =1) {
  const string myname = "test_DuneRoiBuildingService: ";
#ifdef NDEBUG
  cout << myname << "NDEBUG must be off." << endl;
  abort();
#endif
  string line = "-----------------------------";

  cout << myname << line << endl;
  cout << myname << "Create top-level FCL." << endl;
  string fclfile = "test_DuneRoiBuildingService.fcl";
  ofstream fout(fclfile.c_str());
  fout << "#include \"services_dune.fcl\"" << endl;
  fout << "services:      @local::dune35t_services" << endl;
  fout << "services.AdcRoiBuildingService: {" << endl;
  fout << "  service_provider: DuneRoiBuildingService" << endl;
  fout << "  NSigmaStart:  4.0" << endl;
  fout << "  NSigmaEnd:    1.0" << endl;
  fout << "  PadLow:         5" << endl;
  fout << "  PadHigh:       10" << endl;
  fout << "  LogLevel:       " << a_LogLevel << endl;
  fout << "}" << endl;
  fout.close();

  cout << myname << "Fetch art service helper." << endl;
  ArtServiceHelper& ash = ArtServiceHelper::instance();
  ash.print();

  if ( ash.serviceStatus() == 0 ) {

    cout << myname << line << endl;
    cout << myname << "Add supporting services." << endl;
    assert( ash.addService("ExptGeoHelperInterface", fclfile, true) == 0 );
    assert( ash.addService("Geometry", fclfile, true) == 0 );
    assert( ash.addService("LArPropertiesService", fclfile, true) == 0 );
    assert( ash.addService("DetectorClocksService", fclfile, true) == 0 );
    assert( ash.addService("DetectorPropertiesService", fclfile, true) == 0 );
    assert( ash.addService("LArFFT", fclfile, true) == 0 );
    assert( ash.addService("SignalShapingServiceDUNE", fclfile, true) == 0 );
    ash.print();

    cout << myname << line << endl;
    cout << myname << "Add ROI building service." << endl;
    assert( ash.addService("AdcRoiBuildingService", fclfile, true) == 0 );
    ash.print();

    cout << myname << line << endl;
    cout << myname << "Load services." << endl;
    assert( ash.loadServices() == 1 );
    ash.print();
  }

  const unsigned int nsig = 100;
  AdcChannelData acd;
  acd.channel = 100;
  for ( unsigned int isig=0; isig<nsig; ++isig ) {
    acd.samples.push_back(0);
    acd.flags.push_back(AdcGood);
  }
  acd.samples[10] =  3.0;
  acd.samples[29] =  2.0;
  acd.samples[30] =  8.0;
  acd.samples[31] = 10.0;
  acd.samples[32] =  8.0;
  acd.samples[33] =  4.0;
  acd.samples[34] =  2.0;
  acd.samples[35] =  0.5;
  acd.samples[60] = 20.0;
  acd.samples[70] = 20.0;
  acd.samples[90] = 20.0;
  AdcSignalVector& sigs = acd.samples;
  AdcSignalVector sigsin = sigs;

  //cout << myname << "Fetch FFT service." << endl;
  //art::ServiceHandle<util::LArFFT> hfft;
  //unsigned int fftsize = hfft->FFTSize();
  //cout << myname << "Resizing signal vector to FFT size " << fftsize
       //<< " as required for convolution." << endl;
  //sigs.resize(fftsize, 0.0);
  //cout << myname << "Samples size: " << acd.samples.size() << endl;

  cout << myname << "Fetch shaping service." << endl;
  ServiceHandle<util::SignalShapingServiceDUNE> hsss;
  cout << myname << "Decon noise: " << hsss->GetDeconNoise(acd.channel) << endl;

  cout << myname << "Fetch ROI building service." << endl;
  ServiceHandle<AdcRoiBuildingService> hroi;
  hroi->print(cout, myname);

  cout << myname << "Build ROIs." << endl;
  hroi->build(acd);
  cout << myname << "Samples size: " << acd.samples.size() << endl;
  cout << myname << "Output ROIS size: " << acd.rois.size() << endl;
  assert( sigs.size() == nsig );
  for ( unsigned int isig=0; isig<nsig; ++isig ) {
    cout << myname << setw(4) << isig << ": "
         << fixed << setprecision(1) << setw(10) << sigsin[isig]
         << fixed << setprecision(1) << setw(10) << sigs[isig]
         << endl;
    assert( sigs[isig] == sigsin[isig] );
  }
  cout << myname << line << endl;
  cout << myname << "ROIS:" << endl;
  for ( AdcRoi roi : acd.rois ) {
    cout << setw(6) << roi.first << setw(6) << roi.second << endl;
  }

  cout << myname << line << endl;
  cout << myname << "Check ROIs." << endl;
  assert( acd.rois.size() == 3 );
  assert( acd.rois[0].first == 25 );
  assert( acd.rois[0].second == 44 );
  assert( acd.rois[1].first == 55 );
  assert( acd.rois[1].second == 80 );
  assert( acd.rois[2].first == 85 );
  assert( acd.rois[2].second == 99 );

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
  return test_DuneRoiBuildingService(a_LogLevel);
}

//**********************************************************************
