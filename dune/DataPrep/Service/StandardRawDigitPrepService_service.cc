// StandardRawDigitPrepService_service.cc

#include "StandardRawDigitPrepService.h"
#include <iostream>
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "lardata/RawData/RawDigit.h"
#include "dune/DuneInterface/AdcChannelData.h"
#include "dune/DuneInterface/RawDigitExtractService.h"
#include "dune/DuneInterface/AdcMitigationService.h"

using std::string;
using std::cout;
using std::endl;
using std::vector;
using raw::RawDigit;

//**********************************************************************

StandardRawDigitPrepService::
StandardRawDigitPrepService(fhicl::ParameterSet const& pset, art::ActivityRegistry&)
: m_LogLevel(1),
  m_pExtractSvc(nullptr),
  m_pmitigateSvc(nullptr) {
  const string myname = "StandardRawDigitPrepService::ctor: ";
  pset.get_if_present<int>("LogLevel", m_LogLevel);
  m_DoMitigation = pset.get<bool>("DoMitigation");
  if ( m_LogLevel ) cout << myname << "Fetching extract service." << endl;
  m_pExtractSvc = &*art::ServiceHandle<RawDigitExtractService>();
  if ( m_LogLevel ) cout << myname << "  Extract service: @" <<  m_pExtractSvc << endl;
  if ( m_DoMitigation ) {
    if ( m_LogLevel ) cout << myname << "Fetching mitigation service." << endl;
    m_pmitigateSvc = &*art::ServiceHandle<AdcMitigationService>();
    if ( m_LogLevel ) cout << myname << "  Mitigation service: @" <<  m_pmitigateSvc << endl;
  }
  print(cout, myname);
}

//**********************************************************************

int StandardRawDigitPrepService::
prepare(const vector<RawDigit>& digs, AdcChannelDataMap& prepdigs) const {
  const string myname = "StandardRawDigitPrepService:prepare: ";
  if ( m_LogLevel >= 2 ) {
    cout << myname << "Entering..." << endl;
    cout << myname << "Input # input digits: " << digs.size() << endl;
    cout << myname << "Input # prepared digits: " << prepdigs.size() << endl;
  }
  // Extract digits.
  int nbad = 0;
  for ( const RawDigit& dig : digs ) {
    AdcChannelData data;
    AdcChannel& chan = data.channel;
    m_pExtractSvc->extract(dig, &chan, &data.samples, &data.flags);
    data.digit = &dig;
    AdcChannelDataMap::const_iterator idig = prepdigs.find(chan);
    if ( idig != prepdigs.end() ) {
      cout << myname << "WARNING: Data already exists for channel " << chan << ". Skipping." << endl;
      ++nbad;
      continue;
    }
    if ( m_DoMitigation ) {
      m_pmitigateSvc->update(data);
    }
    prepdigs[chan] = data;
  }
  return nbad;
}

//**********************************************************************

std::ostream& StandardRawDigitPrepService::
print(std::ostream& out, std::string prefix) const {
  out << prefix << "StandardRawDigitPrepService:"        << endl;
  out << prefix << "        LogLevel: " << m_LogLevel       << endl;
  return out;
}

//**********************************************************************

DEFINE_ART_SERVICE_INTERFACE_IMPL(StandardRawDigitPrepService, RawDigitPrepService)

//**********************************************************************

