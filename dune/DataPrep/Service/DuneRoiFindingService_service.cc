// DuneRoiFindingService_service.cc

#include "DuneRoiFindingService.h"
#include <iostream>
#include <sstream>
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "dune/DuneInterface/AdcSuppressService.h"
#include "dune/Utilities/SignalShapingServiceDUNE.h"

using std::vector;
using std::string;
using std::ostream;
using std::cout;
using std::endl;
using art::ServiceHandle;

//**********************************************************************

DuneRoiFindingService::
DuneRoiFindingService(fhicl::ParameterSet const& pset, art::ActivityRegistry&)
: m_LogLevel(1) {
  const string myname = "DuneRoiFindingService::ctor: ";
  pset.get_if_present<int>("m_LogLevel", m_LogLevel);
  m_NSigmaStart = pset.get<AdcSignal>("m_NSigmaStart");
  m_NSigmaEnd   = pset.get<AdcSignal>("m_NSigmaEnd");
  m_PadLow      = pset.get<AdcSignal>("m_PadLow");
  m_PadHigh     = pset.get<AdcSignal>("m_PadHigh");
  if ( m_LogLevel > 0 ) print(cout, myname);
}


//**********************************************************************

int DuneRoiFindingService::find(AdcChannelData& data) const {
  const string myname = "DuneRoiFindingService:find: ";
  if ( m_LogLevel >= 2 ) cout << myname << "Finding ROIs for channel "
                              << data.channel << "." << endl;
  data.rois.clear();
  // Get signal shaping service.
  art::ServiceHandle<util::SignalShapingServiceDUNE> hsss;
  AdcSignal sigma = hsss->GetDeconNoise(data.channel);
  const AdcSignalVector& sigs = data.samples;
  // Find ROIS before merging.
  AdcRoiVector rois;
  bool inroi = false;
  AdcIndex isig0 = 0;
  AdcSignal siglow = m_NSigmaStart*sigma;
  AdcSignal sighigh = m_NSigmaEnd*sigma;
  AdcIndex nsig = sigs.size();
  for ( AdcIndex isig=0; isig<sigs.size(); ++isig ) {
    AdcSignal sig = sigs[isig];
    if ( inroi ) {
      if ( sig < siglow || isig == nsig-1 ) {
        rois.push_back(AdcRoi(isig0, isig));
        isig0 = 0;
        inroi = false;
      }
    } else {
      if ( sig > sighigh ) {
        isig0 = isig;
        inroi = true;
      }
    }
  }
  if ( m_LogLevel >= 2 ) cout << myname << "  ROI count before merge: "
                              << rois.size() << endl;
  if ( rois.size() == 0 ) return 0;
  // Merge ROIs.
  AdcIndex iroi = 0;
  AdcIndex lo1 = rois[iroi].first;
  if ( lo1 > m_PadLow ) lo1 =- m_PadLow;
  else lo1 = 0;
  AdcIndex hi1 = rois[iroi].second + m_PadHigh;
  if ( hi1 > nsig ) hi1 = nsig;
  while ( ++iroi < rois.size() ) {
    AdcIndex lo2 = rois[iroi].first;
    if ( lo2 > m_PadLow ) lo1 =- m_PadLow;
    else lo2 = 0;
    AdcIndex hi2 = rois[iroi].second + m_PadHigh;
    if ( hi2 > nsig ) hi2 = nsig;
    bool endroi = false;
    if ( lo2 <= hi1+1 ) {
      hi1 = hi2;
      endroi = iroi+1 == rois.size();
    } else {
      endroi = true;
    }
    if ( endroi ) {
      data.rois.push_back(AdcRoi(lo1, hi1));
      lo1 = lo2;
      hi1 = hi2;
    }
  }
  if ( m_LogLevel >= 2 ) cout << myname << "   ROI count after merge: "
                              << data.rois.size() << endl;
  return 0;
}

//**********************************************************************

ostream& DuneRoiFindingService::
print(ostream& out, string prefix) const {
  out << prefix << "DuneRoiFindingService:" << endl;
  out << "    LogLevel: " << m_LogLevel << endl;
  out << " NSigmaStart: " << m_NSigmaStart << endl;
  out << "   NSigmaEnd: " << m_NSigmaEnd << endl;
  out << "      PadLow: " << m_PadLow << endl;
  out << "     PadHigh: " << m_PadHigh << endl;
  return out;
}

//**********************************************************************

DEFINE_ART_SERVICE_INTERFACE_IMPL(DuneRoiFindingService, AdcRoiFindingService)

//**********************************************************************
