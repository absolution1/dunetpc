// DuneDPhaseDeconvolutionService_service.cc

#include "DuneDPhaseDeconvolutionService.h"
#include <iostream>
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "lardata/Utilities/LArFFT.h"
#include "dune/Utilities/SignalShapingServiceDUNEDPhase.h"

using std::string;
using std::cout;
using std::endl;

//**********************************************************************

DuneDPhaseDeconvolutionService::
DuneDPhaseDeconvolutionService(fhicl::ParameterSet const& pset, art::ActivityRegistry&)
: m_LogLevel(1) {
  const string myname = "DuneDPhaseDeconvolutionService::ctor: ";
  pset.get_if_present<int>("LogLevel", m_LogLevel);
  print(cout, myname);
}

//**********************************************************************

int DuneDPhaseDeconvolutionService::
update(AdcChannelData& data) const {
  const string myname = "DuneDPhaseDeconvolutionService::update: ";
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
  art::ServiceHandle<util::SignalShapingServiceDUNEDPhase> hsss;
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

std::ostream& DuneDPhaseDeconvolutionService::
print(std::ostream& out, std::string prefix) const {
  out << prefix << "DuneDPhaseDeconvolutionService:"                   << endl;
  out << prefix << "               LogLevel: " << m_LogLevel              << endl;
  return out;
}

//**********************************************************************

DEFINE_ART_SERVICE_INTERFACE_IMPL(DuneDPhaseDeconvolutionService, AdcDeconvolutionService)

//**********************************************************************
