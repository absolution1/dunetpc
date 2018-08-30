// VintageDeconvoluter_tool.cc

#include "VintageDeconvoluter.h"
#include <iostream>
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "lardata/Utilities/LArFFT.h"
#include "dune/Utilities/SignalShapingServiceDUNE.h"

using std::string;
using std::cout;
using std::endl;

//**********************************************************************
// Class methods.
//**********************************************************************

VintageDeconvoluter::VintageDeconvoluter(fhicl::ParameterSet const& ps)
: m_LogLevel(ps.get<int>("LogLevel")) {
  const string myname = "VintageDeconvoluter::ctor: ";
  if ( m_LogLevel >= 1 ) {
    cout << myname << "Parameters:" << endl;
    cout << myname << "  LogLevel: " << m_LogLevel << endl;
    art::ServiceHandle<util::SignalShapingServiceDUNE> hsss;
    float dnorm = hsss->GetDeconNorm();
    float norm = dnorm > 0.0 ? 1.0/dnorm : 1.0;
    cout << myname << "Normalization factor: " << norm << endl;
  }
}

//**********************************************************************

DataMap VintageDeconvoluter::update(AdcChannelData& acd) const {
  const string myname = "VintageDeconvoluter::view: ";
  DataMap ret;
  AdcChannel chan = acd.channel;
  AdcSignalVector& samples = acd.samples;
  unsigned int nsam = samples.size();
  if ( m_LogLevel >= 2 ) cout << myname << "Deconvoluting channel " << chan
                             << " with " << nsam << " samples." << endl;
  // Fetch the FFT size.
  art::ServiceHandle<util::LArFFT> hFFT;
  unsigned int fftsize = hFFT->FFTSize();
  if ( nsam > fftsize ) {
    cout << myname << "ERROR: Data has too many ticks for FFT: "
         << nsam << " > " << fftsize << "." << endl;
    return ret.setStatus(1);
  }
  // Pad the data to the FFT size.
  bool pad = fftsize > nsam;
  if ( pad ) {
    if ( m_LogLevel >= 3 ) cout << myname << "  Padding sample vector to " << fftsize << endl;
    samples.resize(fftsize);
    for ( unsigned int isam=nsam; isam<fftsize; ++isam ) {
      //samples[isam] = samples[isam-nsam];
      samples[isam] = 0;
    }
  }
  // Deconvolute.
  art::ServiceHandle<util::SignalShapingServiceDUNE> hsss;
  if ( m_LogLevel >= 3 ) cout << myname << "  Deconvoluting." << endl;
  hsss->Deconvolute(chan, samples);
  if ( pad ) samples.resize(nsam);
  if ( m_LogLevel >= 3 ) cout << myname << "  Normalizing." << endl;
  float dnorm = hsss->GetDeconNorm();
  float norm = dnorm > 0.0 ? 1.0/dnorm : 1.0;
  if ( m_LogLevel >= 3 ) cout << myname << "  Scale factor: " << norm << endl;
  for ( float& sam : samples ) sam *= norm;
  acd.sampleUnit = "Q_{dco}";
  acd.sampleNoise = hsss->GetDeconNoise(chan);
  // Done.
  return ret;
}

//**********************************************************************
