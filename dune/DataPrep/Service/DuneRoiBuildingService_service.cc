// DuneRoiBuildingService_service.cc

#include "DuneRoiBuildingService.h"
#include <iostream>
#include <sstream>
#include <iomanip>
#include "art/Framework/Services/Registry/ServiceHandle.h"
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
  if ( m_LogLevel >= 3 ) {
    cout << myname << "  Noise level: " << sigma << " " << data.sampleUnit << endl;
  }
  const AdcSignalVector& sigs = data.samples;
  // Build ROIS before padding and merging.
  AdcFilterVector& signal = data.signal;
  AdcRoiVector& rois = data.rois;
  signal.clear();
  signal.resize(sigs.size(), false);
  bool inroi = false;
  AdcSignal siglow = m_NSigmaEnd*sigma;
  AdcSignal sighigh = m_NSigmaStart*sigma;
  AdcIndex nsig = sigs.size();
  if ( nsig < 1 ) {
    if ( m_LogLevel >= 2 ) cout << myname << "Channel " << data.channel
                                << " has no samples." << endl;
    return 0;
  }
  for ( AdcIndex isig=0; isig<sigs.size(); ++isig ) {
    AdcSignal sig = sigs[isig];
    if ( inroi ) {
      if ( sig > siglow ) {
        signal[isig] = true;
      } else  {
        inroi = false;
      }
    } else {
      if ( sig > sighigh ) {
        signal[isig] = true;
        inroi = true;
      }
    }
  }
  // Fill the unpadded ROIs.
  data.roisFromSignal();
  // Display ROIs before padding and merging.
  if ( m_LogLevel >= 3 ) {
    cout << myname << "  ROIs before merge (size = " << rois.size() << "):" << endl;
    for ( const AdcRoi& roi : rois ) {
      cout << myname << setw(8) << roi.first << " " << setw(8) << roi.second << endl;
    }
  } else if ( m_LogLevel >= 2 ) {
    cout << myname << "  ROI count before merge: " << data.rois.size() << endl;
  }
  if ( rois.size() == 0 ) return 0;
  // Pad ROIs.
  unsigned int isig1 = 0;
  unsigned int isig2 = 0;
  for ( AdcRoi roi : rois ) {
    isig2 = roi.first;
    isig1 = isig2 > m_PadLow ? isig2 - m_PadLow : 0;
    for ( unsigned int isig=isig1; isig<isig2; ++isig ) signal[isig] = true;
    isig1 = roi.second + 1;
    isig2 = isig1 + m_PadHigh;
    if ( isig2 > nsig ) isig2 = nsig;
    for ( unsigned int isig=isig1; isig<isig2; ++isig ) signal[isig] = true;
  }
  // Fill the final ROIs.
  data.roisFromSignal();
  // Display final ROIs.
  if ( m_LogLevel >= 3 ) {
    cout << myname << "  ROIs after merge (size = " << rois.size() << "):" << endl;
    for ( const AdcRoi& roi : rois ) {
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
