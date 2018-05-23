// StandardRawDigitExtractService_service.cc

#include "StandardRawDigitExtractService.h"
#include <iostream>
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "lardataobj/RawData/RawDigit.h"
#include "larevt/CalibrationDBI/Interface/DetPedestalService.h"
#include "larevt/CalibrationDBI/Interface/DetPedestalProvider.h"
#include "dune/ArtSupport/DuneToolManager.h"
#include "dune/DuneInterface/Tool/AdcChannelTool.h"

using std::string;
using std::cout;
using std::endl;

//**********************************************************************

StandardRawDigitExtractService::
StandardRawDigitExtractService(fhicl::ParameterSet const& pset, art::ActivityRegistry&)
: m_LogLevel(1),
  m_pDigitReadTool(nullptr),
  m_pROIBuilderTool(nullptr),
  m_pPedProv(nullptr) {
  const string myname = "StandardRawDigitExtractService::ctor: ";
  pset.get_if_present<int>("LogLevel", m_LogLevel);
  pset.get_if_present<string>("ROIBuilderTool", m_ROIBuilderTool);
  m_DigitReadTool  = pset.get<string>("DigitReadTool");
  m_PedestalOption = pset.get<int>("PedestalOption");
  m_FlagStuckOff   = pset.get<bool>("FlagStuckOff");
  m_FlagStuckOn    = pset.get<bool>("FlagStuckOn");
  // Retrieve tools
  DuneToolManager* ptm = DuneToolManager::instance("");
  if ( ptm == nullptr ) {
    cout << myname << "ERROR: Unable to retrieve tool manaager." << endl;
  } else {
    m_pDigitReadTool = ptm->getPrivate<AdcChannelTool>(m_DigitReadTool);
    m_pROIBuilderTool = ptm->getPrivate<AdcChannelTool>(m_ROIBuilderTool);
    if ( m_pDigitReadTool == nullptr ) {
      cout << myname << "ERROR: Unable to retrieve digit reader " << m_DigitReadTool << endl;
    } else {
      if ( m_LogLevel ) cout << myname << "Retrieved digit read tool " << m_DigitReadTool << endl;
    }
  }
  // Retrieve pedestal provider.
  if ( m_PedestalOption == 2 ) {
    if ( m_LogLevel ) cout << myname << "Fetching pedestal provider." << endl;
    m_pPedProv = &art::ServiceHandle<lariov::DetPedestalService>()->GetPedestalProvider();
    if ( m_LogLevel ) cout << myname << "  Pedestal provider: @" <<  m_pPedProv << endl;
  }
  if ( m_PedestalOption == 3 ) {
    if ( m_LogLevel ) cout << myname << "Fetching pedestal evaluation service." << endl;
    m_PedestalEvaluationService = &*art::ServiceHandle<PedestalEvaluationService>();
    if ( m_LogLevel ) cout << myname << "  Pedestal evaluation service: @"
                           <<  m_PedestalEvaluationService << endl;
  }
  print(cout, myname);
}

//**********************************************************************

int StandardRawDigitExtractService::extract(AdcChannelData& acd) const {
  const string myname = "StandardRawDigitExtractService:extract: ";
  const raw::RawDigit* pdig = acd.digit;
  if ( pdig == nullptr ) {
    cout << myname << "ERROR: ADC channel does not have a larsoft digit." << endl;
    return 1;
  }
  const raw::RawDigit& dig = *pdig;
  if ( m_LogLevel >= 2 ) {
    cout << myname << "Entering..." << endl;
    cout << myname << "Input vector size: " << dig.Samples() << endl;
  }
  if ( acd.samples.size() ) {
    cout << myname << "ERROR: Channel has data." << endl;
    return 2;
  }
  if ( acd.flags.size() ) {
    cout << myname << "ERROR: Channel has flags." << endl;
    return 3;
  }
  if ( m_pDigitReadTool ) {
    m_pDigitReadTool->update(acd);
  } else {
    cout << myname << "ERROR: Digit read tool " << m_DigitReadTool << " was not found." << endl;
  }
  unsigned int nsig = acd.raw.size();
  acd.flags.resize(nsig, AdcGood);
  acd.samples.resize(nsig, -999);
  // Retrieve pedestal.
  AdcSignal ped = 0.0;
  AdcSignal pedrms = 0.0;
  if ( m_PedestalOption == 1 ) {
    ped = dig.GetPedestal();
    pedrms = dig.GetSigma();
  } else if ( m_PedestalOption == 2 ) {
    ped = m_pPedProv->PedMean(acd.channel);
  }
  acd.pedestal = ped;
  acd.pedestalRms = pedrms;
  // Convert int -> float, subtract pedestal and set conversion flag.
  const AdcCount lowbits = 0x3f;
  for ( unsigned int isig=0; isig<nsig; ++isig ) {
    AdcCount adc = acd.raw[isig];
    AdcCount adclow = adc & lowbits;
    acd.samples[isig] = adc - ped;
    if      ( adc == 0 )                            acd.flags[isig] = AdcUnderflow;
    else if ( adc >= 4095 )                         acd.flags[isig] = AdcOverflow;
    else if ( m_FlagStuckOff && adclow == 0 )       acd.flags[isig] = AdcStuckOff;
    else if ( m_FlagStuckOn  && adclow == lowbits ) acd.flags[isig] = AdcStuckOn;
  }
  if ( m_PedestalOption == 3 ) {

    //optional: build ROI and exclude ROI from pedestal calculation
    if( m_ROIBuilderTool.size() )
    {
      m_pROIBuilderTool->update(acd);
    }

    m_PedestalEvaluationService->evaluate(acd, &ped);
    for ( unsigned int isig=0; isig<nsig; ++isig ) {
      acd.samples[isig] -= ped;
    }
    acd.pedestal += ped;
  }
  return 0;
}

//**********************************************************************

std::ostream& StandardRawDigitExtractService::
print(std::ostream& out, std::string prefix) const {
  out << prefix << "StandardRawDigitExtractService:"        << endl;
  out << prefix << "        LogLevel: " << m_LogLevel       << endl;
  out << prefix << "   DigitReadTool: " << m_DigitReadTool  << endl;
  out << prefix << "  PedestalOption: " << m_PedestalOption << endl;
  out << prefix << "    FlagStuckOff: " << m_FlagStuckOff   << endl;
  out << prefix << "     FlagStuckOn: " << m_FlagStuckOn    << endl;
  return out;
}

//**********************************************************************

DEFINE_ART_SERVICE_INTERFACE_IMPL(StandardRawDigitExtractService, RawDigitExtractService)

//**********************************************************************

