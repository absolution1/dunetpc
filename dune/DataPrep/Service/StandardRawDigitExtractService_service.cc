// StandardRawDigitExtractService_service.cc

#include "StandardRawDigitExtractService.h"
#include <iostream>
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "lardata/RawData/RawDigit.h"
#include "lardata/RawData/raw.h"
#include "larevt/CalibrationDBI/Interface/DetPedestalService.h"
#include "larevt/CalibrationDBI/Interface/DetPedestalProvider.h"

using std::string;
using std::cout;
using std::endl;

//**********************************************************************

StandardRawDigitExtractService::
StandardRawDigitExtractService(fhicl::ParameterSet const& pset, art::ActivityRegistry&)
: m_LogLevel(1),
  m_pPedProv(nullptr) {
  const string myname = "StandardRawDigitExtractService::ctor: ";
  pset.get_if_present<int>("LogLevel", m_LogLevel);
  m_PedestalOption = pset.get<int>("PedestalOption");
  m_FlagStuckOff   = pset.get<bool>("FlagStuckOff");
  m_FlagStuckOn    = pset.get<bool>("FlagStuckOn");
  if ( m_PedestalOption == 2 ) {
    if ( m_LogLevel ) cout << myname << "Fetching pedestal provider." << endl;
    m_pPedProv = &art::ServiceHandle<lariov::DetPedestalService>()->GetPedestalProvider();
    if ( m_LogLevel ) cout << myname << "  Pedestal provider: @" <<  m_pPedProv << endl;
  }
  print(cout, myname);
}

//**********************************************************************

int StandardRawDigitExtractService::
extract(const raw::RawDigit& dig, AdcChannel* pchan, AdcSignal* pped,
        AdcSignalVector* psigs_in, AdcFlagVector* pflgs_in) const {
  const string myname = "StandardRawDigitExtractService:extract: ";
  if ( m_LogLevel >= 2 ) {
    cout << myname << "Entering..." << endl;
    cout << myname << "Input vector size: " << dig.Samples() << endl;
    cout << myname << "Signal vector ";
    if ( psigs_in == nullptr ) cout << "not ";
    cout << "requested." << endl;
    cout << myname << "Flags ";
    if ( pflgs_in == nullptr ) cout << "not ";
    cout << "requested." << endl;
  }
  AdcChannel chan = dig.Channel();
  if ( pchan != nullptr ) *pchan = chan;
  unsigned int nsig = dig.Samples();
  // Initialize the output signal and flag vectors.
  AdcSignalVector* psigs_local = nullptr;
  if ( psigs_in == nullptr ) psigs_local = new AdcSignalVector;
  AdcSignalVector& sigs = psigs_in == nullptr ? *psigs_local : *psigs_in;
  sigs.resize(nsig, -999);     // See https://cdcvs.fnal.gov/redmine/issues/11572.
  AdcFlagVector* pflgs_local = nullptr;
  if ( pflgs_in == nullptr ) pflgs_local = new AdcFlagVector;
  AdcFlagVector& flgs = pflgs_in == nullptr ? *pflgs_local : *pflgs_in;
  flgs.resize(nsig, AdcGood);
  // Extract the signals from the digit.
  AdcCountVector adcs(nsig, -999);   // See https://cdcvs.fnal.gov/redmine/issues/11572.
  raw::Uncompress(dig.ADCs(), adcs, dig.GetPedestal(), dig.Compression());
  // Retrieve pedestal.
  AdcSignal ped = 0.0;
  if ( m_PedestalOption == 1 ) {
    ped = dig.GetPedestal();
  } else if ( m_PedestalOption == 2 ) {
    ped = m_pPedProv->PedMean(chan);
  }
  if ( pped != nullptr ) *pped = ped;
  // Convert int -> float, subtract pedestal and set conversion flag.
  const AdcCount lowbits = 0x3f;
  for ( unsigned int isig=0; isig<nsig; ++isig ) {
    AdcCount adc = adcs[isig];
    AdcCount adclow = adc & lowbits;
    sigs[isig] = adc - ped;
    if      ( adc == 0 )                            flgs[isig] = AdcUnderflow;
    else if ( adc >= 4095 )                         flgs[isig] = AdcOverflow;
    else if ( m_FlagStuckOff && adclow == 0 )       flgs[isig] = AdcStuckOff;
    else if ( m_FlagStuckOn  && adclow == lowbits ) flgs[isig] = AdcStuckOn;
  }

  return 0;
}

//**********************************************************************

std::ostream& StandardRawDigitExtractService::
print(std::ostream& out, std::string prefix) const {
  out << prefix << "StandardRawDigitExtractService:"        << endl;
  out << prefix << "        LogLevel: " << m_LogLevel       << endl;
  out << prefix << "  PedestalOption: " << m_PedestalOption << endl;
  out << prefix << "    FlagStuckOff: " << m_FlagStuckOff   << endl;
  out << prefix << "     FlagStuckOn: " << m_FlagStuckOn    << endl;
  return out;
}

//**********************************************************************

DEFINE_ART_SERVICE_INTERFACE_IMPL(StandardRawDigitExtractService, RawDigitExtractService)

//**********************************************************************

