// test_DuneDeconvolutionService.cxx
//
// David Adams
// May 2016
//
// Test DuneDeconvolutionService.
//
// This test crashes in cleanup until this problem is resolved:
// https://cdcvs.fnal.gov/redmine/issues/10618

#include <string>
#include <iostream>
#include <sstream>
#include <fstream>
#include <iomanip>
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "dune/ArtSupport/ArtServiceHelper.h"
#include "dune/Utilities/SignalShapingServiceDUNE.h"
#include "dune/DuneInterface/AdcDeconvolutionService.h"

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

int test_DuneDeconvolutionService(int a_LogLevel =-1) {
  const string myname = "test_DuneDeconvolutionService: ";
#ifdef NDEBUG
  cout << myname << "NDEBUG must be off." << endl;
  abort();
#endif
  string line = "-----------------------------";

  if ( a_LogLevel < 0 ) {
    cout << myname << "Skipping test while we wat for resolution of larsoft issue 10618" << endl;
    return 0;
  }

  cout << myname << line << endl;
  cout << myname << "Create top-level FCL." << endl;
  string fclfile = "test_DuneDeconvolutionService.fcl";
  ofstream fout(fclfile.c_str());
  fout << "#include \"services_dune.fcl\"" << endl;
  fout << "services:      @local::dune35t_services" << endl;
  fout << "services.AdcDeconvolutionService: {" << endl;
  fout << "  service_provider: DuneDeconvolutionService" << endl;
  fout << "  LogLevel:              " << a_LogLevel << endl;
  fout << "}" << endl;
  fout.close();

  cout << myname << "Fetch art service helper." << endl;
  ArtServiceHelper& ash = ArtServiceHelper::instance();
  ash.print();

  if ( ash.serviceStatus() == 0 ) {

    cout << myname << line << endl;
    cout << myname << "Add services." << endl;
    assert( ash.addServices(fclfile, true) == 0 );

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
  acd.samples[50] = 100000.0;
  AdcSignalVector& sigs = acd.samples;
  AdcSignalVector sigsin = sigs;

  cout << myname << "Fetch FFT service." << endl;
  art::ServiceHandle<util::LArFFT> hfft;
  unsigned int fftsize = hfft->FFTSize();
  cout << myname << "Resizing signal vector to FFT size " << fftsize
       << " as required for convolution." << endl;
  sigs.resize(fftsize, 0.0);

  cout << myname << "Fetch shaping service." << endl;
  ServiceHandle<util::SignalShapingServiceDUNE> hsss;

  cout << myname << "Convolute." << endl;
  hsss->Convolute(acd.channel, sigs);
  AdcSignalVector sigsco = sigs;

  cout << myname << "Fetch deconvolution service." << endl;
  ServiceHandle<AdcDeconvolutionService> hdco;
  hdco->print();

  cout << myname << line << endl;
  cout << myname << "Deconvolute." << endl;
  assert( hdco->update(acd) == 0 );
  cout << myname << "Output vector size: " << sigs.size() << endl;
  assert( sigs.size() == sigsco.size() );
  for ( unsigned int isig=0; isig<nsig; ++isig ) {
    cout << setw(4) << isig << ": "
         << fixed << setprecision(1) << setw(8) << sigsin[isig]
         << fixed << setprecision(1) << setw(10) << sigsco[isig]
         << fixed << setprecision(1) << setw(10) << sigs[isig]
         << endl;
    //assert( sigs[isig] == sigsexp[isig] );
  }

  cout << myname << "Done." << endl;
  return 0;
}

//**********************************************************************

int main(int argc, char* argv[]) {
  const string myname = "main: ";
  int a_LogLevel = 1;
  if ( argc > 1 ) {
    istringstream ssarg(argv[1]);
    ssarg >> a_LogLevel;
  }
  int rstat = test_DuneDeconvolutionService(a_LogLevel);
  cout << myname << "Closing service helper." << endl;
  ArtServiceHelper::close();
  cout << myname << "Exiting." << endl;
  return rstat;
}

//**********************************************************************
