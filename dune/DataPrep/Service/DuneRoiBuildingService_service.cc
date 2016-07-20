// DuneRoiBuildingService_service.cc

#include "DuneRoiBuildingService.h"
#include <iostream>
#include <sstream>
#include <iomanip>
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "dune/DuneInterface/AdcSuppressService.h"
#include "dune/Utilities/SignalShapingServiceDUNE.h"

using std::vector;
using std::string;
using std::ostream;
using std::cout;
using std::endl;
using std::setw;
using art::ServiceHandle;

//**********************************************************************

DuneRoiBuildingService::
DuneRoiBuildingService(fhicl::ParameterSet const& pset, art::ActivityRegistry&)
: m_LogLevel(1) {
  const string myname = "DuneRoiBuildingService::ctor: ";
  pset.get_if_present<int>("LogLevel", m_LogLevel);
  m_NSigmaStart = pset.get<AdcSignal>("NSigmaStart");
  m_NSigmaEnd   = pset.get<AdcSignal>("NSigmaEnd");
  m_PadLow      = pset.get<AdcSignal>("PadLow");
  m_PadHigh     = pset.get<AdcSignal>("PadHigh");
  if ( m_LogLevel > 0 ) print(cout, myname);
}


//**********************************************************************

int DuneRoiBuildingService::build(AdcChannelData& data) const {
  const string myname = "DuneRoiBuildingService:build: ";
  if ( m_LogLevel >= 2 ) cout << myname << "Building ROIs for channel "
                              << data.channel << "." << endl;
  data.rois.clear();
  // Get signal shaping service.
  art::ServiceHandle<util::SignalShapingServiceDUNE> hsss;
  AdcSignal sigma = hsss->GetDeconNoise(data.channel);
  const AdcSignalVector& sigs = data.samples;
  // Build ROIS before merging.
  AdcRoiVector rois;
  bool inroi = false;
  AdcIndex isig0 = 0;
  AdcSignal siglow = m_NSigmaEnd*sigma;
  AdcSignal sighigh = m_NSigmaStart*sigma;
  AdcIndex nsig = sigs.size();
  for ( AdcIndex isig=0; isig<sigs.size(); ++isig ) {
    AdcSignal sig = sigs[isig];
    if ( inroi ) {
      if ( sig < siglow || isig == nsig-1 ) {
        rois.push_back(AdcRoi(isig0, isig-1));
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
  // Display ROIs before padding and merging.
  if ( m_LogLevel >= 3 ) {
    cout << myname << "  ROIs before merge (size = " << rois.size() << "):" << endl;
    for ( const AdcRoi& roi : rois ) {
      cout << myname << setw(8) << roi.first << " " << setw(8) << roi.second << endl;
    }
  } else if ( m_LogLevel >= 2 ) {
    cout << myname << "  ROI count before merge: " << rois.size() << endl;
  }
  if ( rois.size() == 0 ) return 0;
  // Pad and merge ROIs. The current ROI is (lo1, hi1).
  AdcIndex lo1 = 0;
  AdcIndex hi1 = 0;
  // Loop over unpadded and unmerged ROIs.
  for ( AdcIndex iroi=0; iroi < rois.size(); ++iroi ) {
    // Pad the new ROI and store it in (lo2, hi2).
    AdcIndex lo2 = rois[iroi].first;
    if ( lo2 > m_PadLow ) lo2 -= m_PadLow;
    else lo2 = 0;
    AdcIndex hi2 = rois[iroi].second + m_PadHigh;
    if ( hi2 >= nsig ) hi2 = nsig-1;
    // First ROI. Make it the current ROI.
    if ( iroi == 0 ) {
      lo1 = lo2;
      hi1 = hi2;
    // ROI overlaps the current ROI. Merge by extending the current ROI.
    } else if ( lo2 <= hi1 + 1 ) {
      hi1 = hi2;
    // New ROI. Save the current ROI and make the new ROI current.
    } else {
      data.rois.push_back(AdcRoi(lo1, hi1));
      lo1 = lo2;
      hi1 = hi2;
    }
  }
  // Save the last ROI.
  data.rois.push_back(AdcRoi(lo1, hi1));
  // Display final ROIs.
  if ( m_LogLevel >= 3 ) {
    cout << myname << "  ROIs after merge (size = " << data.rois.size() << "):" << endl;
    for ( const AdcRoi& roi : data.rois ) {
      cout << myname << setw(8) << roi.first << " " << setw(8) << roi.second << endl;
    }
  } else if ( m_LogLevel >= 2 ) {
    cout << myname << "  ROI count after merge: " << data.rois.size() << endl;
  }
  return 0;
}

//**********************************************************************

ostream& DuneRoiBuildingService::
print(ostream& out, string prefix) const {
  out << prefix << "DuneRoiBuildingService:" << endl;
  out << prefix << "    LogLevel: " << m_LogLevel << endl;
  out << prefix << " NSigmaStart: " << m_NSigmaStart << endl;
  out << prefix << "   NSigmaEnd: " << m_NSigmaEnd << endl;
  out << prefix << "      PadLow: " << m_PadLow << endl;
  out << prefix << "     PadHigh: " << m_PadHigh << endl;
  return out;
}

//**********************************************************************

DEFINE_ART_SERVICE_INTERFACE_IMPL(DuneRoiBuildingService, AdcRoiBuildingService)

//**********************************************************************
