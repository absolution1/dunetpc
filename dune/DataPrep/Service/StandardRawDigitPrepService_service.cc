// StandardRawDigitPrepService_service.cc

#include "StandardRawDigitPrepService.h"
#include <iostream>
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "lardataobj/RawData/RawDigit.h"
#include "dune/DuneInterface/AdcChannelData.h"
#include "dune/DuneInterface/RawDigitExtractService.h"
#include "dune/DuneInterface/AdcMitigationService.h"
#include "dune/DuneInterface/AdcSignalFindingService.h"
#include "dune/DuneInterface/AdcNoiseRemovalService.h"
#include "dune/DuneInterface/PedestalEvaluationService.h"
#include "dune/DuneInterface/AdcDeconvolutionService.h"
#include "dune/DuneInterface/AdcRoiBuildingService.h"
#include "dune/DuneInterface/AdcWireBuildingService.h"

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
  m_pmitigateSvc(nullptr),
  m_pAdcSignalFindingService(nullptr),
  m_pNoiseRemoval(nullptr),
  m_pPedestalEvaluation(nullptr),
  m_pDeconvolutionService(nullptr),
  m_pRoiBuildingService(nullptr),
  m_pWireBuildingService(nullptr) {
  const string myname = "StandardRawDigitPrepService::ctor: ";
  pset.get_if_present<int>("LogLevel", m_LogLevel);
  m_DoMitigation = pset.get<bool>("DoMitigation");
  m_DoEarlySignalFinding = pset.get<bool>("DoEarlySignalFinding");
  m_DoNoiseRemoval       = pset.get<bool>("DoNoiseRemoval");
  m_DoPedestalAdjustment = pset.get<bool>("DoPedestalAdjustment");
  m_DoDeconvolution      = pset.get<bool>("DoDeconvolution");
  m_DoROI                = pset.get<bool>("DoROI");
  m_DoWires              = pset.get<bool>("DoWires");
  if ( m_LogLevel ) cout << myname << "Fetching extract service." << endl;
  m_pExtractSvc = &*art::ServiceHandle<RawDigitExtractService>();
  if ( m_LogLevel ) cout << myname << "  Extract service: @" <<  m_pExtractSvc << endl;
  if ( m_DoMitigation ) {
    if ( m_LogLevel ) cout << myname << "Fetching mitigation service." << endl;
    m_pmitigateSvc = &*art::ServiceHandle<AdcMitigationService>();
    if ( m_LogLevel ) cout << myname << "  Mitigation service: @" <<  m_pmitigateSvc << endl;
  }
  if ( m_DoEarlySignalFinding ) {
    if ( m_LogLevel ) cout << myname << "Fetching signal finding service." << endl;
    m_pAdcSignalFindingService = &*art::ServiceHandle<AdcSignalFindingService>();
    if ( m_LogLevel ) cout << myname << "  Signal finding service: @" <<  m_pAdcSignalFindingService << endl;
  }
  if ( m_DoNoiseRemoval ) {
    if ( m_LogLevel ) cout << myname << "Fetching noise removal service." << endl;
    m_pNoiseRemoval = &*art::ServiceHandle<AdcNoiseRemovalService>();
    if ( m_LogLevel ) cout << myname << "  Noise removal service: @" <<  m_pNoiseRemoval << endl;
  }
  if ( m_DoPedestalAdjustment ) {
    if ( m_LogLevel ) cout << myname << "Fetching pedestal evaluation service." << endl;
    m_pPedestalEvaluation = &*art::ServiceHandle<PedestalEvaluationService>();
    if ( m_LogLevel ) cout << myname << "  Pedestal evalution service: @" <<  m_pPedestalEvaluation << endl;
  }
  if ( m_DoDeconvolution ) {
    if ( m_LogLevel ) cout << myname << "Fetching deconvolution service." << endl;
    m_pDeconvolutionService = &*art::ServiceHandle<AdcDeconvolutionService>();
    if ( m_LogLevel ) cout << myname << "  Deconvolution service: @" <<  m_pDeconvolutionService << endl;
  }
  if ( m_DoROI ) {
    if ( m_LogLevel ) cout << myname << "Fetching ROI building service." << endl;
    m_pRoiBuildingService = &*art::ServiceHandle<AdcRoiBuildingService>();
    if ( m_LogLevel ) cout << myname << "  ROI building service: @" <<  m_pRoiBuildingService << endl;
  }
  if ( m_DoWires ) {
    if ( m_LogLevel ) cout << myname << "Fetching wire building service." << endl;
    m_pWireBuildingService = &*art::ServiceHandle<AdcWireBuildingService>();
    if ( m_LogLevel ) cout << myname << "  Wire building service: @" <<  m_pWireBuildingService << endl;
  }
  print(cout, myname);
}

//**********************************************************************

int StandardRawDigitPrepService::
prepare(const vector<RawDigit>& digs, AdcChannelDataMap& datamap,
        std::vector<recob::Wire>* pwires) const {
  const string myname = "StandardRawDigitPrepService:prepare: ";
  if ( m_LogLevel >= 2 ) {
    cout << myname << "Entering..." << endl;
    cout << myname << "Input # input digits: " << digs.size() << endl;
    cout << myname << "Input # prepared digits: " << datamap.size() << endl;
  }
  // Extract digits.
  int nbad = 0;
  for ( const RawDigit& dig : digs ) {
    AdcChannelData data;
    AdcChannel& chan = data.channel;
    AdcSignal& ped = data.pedestal;
    m_pExtractSvc->extract(dig, &chan, &ped, &data.raw, &data.samples, &data.flags);
    data.digit = &dig;
    AdcChannelDataMap::const_iterator idig = datamap.find(chan);
    if ( idig != datamap.end() ) {
      cout << myname << "WARNING: Data already exists for channel " << chan << ". Skipping." << endl;
      ++nbad;
      continue;
    }
    if ( m_DoMitigation ) {
      m_pmitigateSvc->update(data);
    }
    if ( m_DoEarlySignalFinding ) {
      m_pAdcSignalFindingService->find(data);
    }
    datamap[chan] = data;
  }
  if ( m_DoNoiseRemoval ) {
    m_pNoiseRemoval->update(datamap);
  }
  if ( m_DoDeconvolution ) {
    for ( AdcChannelDataMap::value_type chdata : datamap ) {
      m_pDeconvolutionService->update(chdata.second);
    }
  }
  if ( m_DoPedestalAdjustment ) {
    for ( auto& chdata : datamap ) {
      AdcChannelData& data = chdata.second;
      AdcSignal ped = 0.0;
      m_pPedestalEvaluation->evaluate(data, &ped);
      for ( AdcSignal& sig : data.samples ) sig -= ped;
    }
  }
  if ( m_DoROI ) {
    for ( auto& chdata : datamap ) {
      AdcChannelData& acd = chdata.second;
      m_pRoiBuildingService->build(acd);
    }
  }
  if ( m_DoWires ) {
    for ( auto& chdata : datamap ) {
      AdcChannelData& acd = chdata.second;
      m_pWireBuildingService->build(acd, pwires);
    }
  }
  if ( m_LogLevel >=1 ) print(cout, myname);
  return nbad;
}

//**********************************************************************

std::ostream& StandardRawDigitPrepService::
print(std::ostream& out, std::string prefix) const {
  out << prefix << "StandardRawDigitPrepService:"                      << endl;
  out << prefix << "             LogLevel: " << m_LogLevel             << endl;
  out << prefix << "         DoMitigation: " << m_DoMitigation         << endl;
  out << prefix << " DoEarlySignalFinding: " << m_DoEarlySignalFinding << endl;
  out << prefix << "       DoNoiseRemoval: " << m_DoNoiseRemoval       << endl;
  out << prefix << "      DoDeconvolution: " << m_DoDeconvolution      << endl;
  out << prefix << " DoPedestalAdjustment: " << m_DoPedestalAdjustment << endl;
  out << prefix << "                DoROI: " << m_DoROI                << endl;
  out << prefix << "               DoWires " << m_DoWires              << endl;
  return out;
}

//**********************************************************************

DEFINE_ART_SERVICE_INTERFACE_IMPL(StandardRawDigitPrepService, RawDigitPrepService)

//**********************************************************************

