// DuneDeconvolutionService_service.cc

#include "DuneDeconvolutionService.h"
#include <iostream>
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "lardata/Utilities/LArFFT.h"
#include "dune/Utilities/SignalShapingServiceDUNE.h"

using std::string;
using std::cout;
using std::endl;

//**********************************************************************

DuneDeconvolutionService::
DuneDeconvolutionService(fhicl::ParameterSet const& pset, art::ActivityRegistry&)
: m_LogLevel(1) {
  const string myname = "DuneDeconvolutionService::ctor: ";
  pset.get_if_present<int>("LogLevel", m_LogLevel);
  print(cout, myname);
}

//**********************************************************************

int DuneDeconvolutionService::
update(AdcChannelData& data) const {
  const string myname = "DuneDeconvolutionService::update: ";
  AdcChannel chan = data.channel;
  AdcSignalVector& samples = data.samples;
  unsigned int nsam = samples.size();
  if ( m_LogLevel >= 2 ) cout << myname << "Deconvoluting channel " << chan
                             << " with " << nsam << " samples." << endl;
  // Fetch the FFT size.
  art::ServiceHandle<util::LArFFT> hFFT;
  unsigned int fftsize = hFFT->FFTSize();
  if ( nsam > fftsize ) {
    cout << myname << "ERROR: Data has too many ticks for FFT: "
         << nsam << " > " << fftsize << "." << endl;
    return 1;
  }
  // Pad the data to the FFT size.
  bool pad = fftsize > nsam;
  if ( pad ) {
    if ( m_LogLevel >= 3 ) cout << myname << "  Padding sample vector to " << fftsize << endl;
    samples.resize(fftsize);
    for ( unsigned int isam=nsam; isam<fftsize; ++isam ) {
      samples[isam] = samples[isam-nsam];
    }
  }
  // Deconvolute.
  art::ServiceHandle<util::SignalShapingServiceDUNE> hsss;
  if ( m_LogLevel >= 3 ) cout << myname << "  Deconvoluting." << endl;
  hsss->Deconvolute(chan, samples);
  if ( pad ) samples.resize(nsam);
  if ( m_LogLevel >= 3 ) cout << myname << "  Normalizing." << endl;
  float norm = 1.0/hsss->GetDeconNorm();
  if ( m_LogLevel >= 3 ) cout << myname << "  Scale factor: " << norm << endl;
  for ( float& sam : samples ) sam *= norm;
  // Done.
  return 0;
}

//**********************************************************************

std::ostream& DuneDeconvolutionService::
print(std::ostream& out, std::string prefix) const {
  out << prefix << "DuneDeconvolutionService:"                   << endl;
  out << prefix << "               LogLevel: " << m_LogLevel              << endl;
  return out;
}

//**********************************************************************

DEFINE_ART_SERVICE_INTERFACE_IMPL(DuneDeconvolutionService, AdcDeconvolutionService)

//**********************************************************************
