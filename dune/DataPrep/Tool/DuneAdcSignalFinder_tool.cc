// DuneAdcSignalFinder_tool.cc

#include "DuneAdcSignalFinder.h"
#include <iostream>
#include <iomanip>

using std::string;
using std::cout;
using std::endl;
using std::setw;

//**********************************************************************
// Class methods.
//**********************************************************************

DuneAdcSignalFinder::DuneAdcSignalFinder(fhicl::ParameterSet const& ps)
: m_LogLevel(ps.get<int>("LogLevel")),
  m_NoiseSigma(ps.get<float>("NoiseSigma")),
  m_NSigmaStart(ps.get<float>("NSigmaStart")),
  m_NSigmaEnd(ps.get<float>("NSigmaEnd")),
  m_TicksBefore(ps.get<AdcIndex>("TicksBefore")),
  m_TicksAfter(ps.get<AdcIndex>("TicksAfter"))
{
  const string myname = "DuneAdcSignalFinder::ctor: ";
  if ( m_LogLevel >= 1 ) {
    cout << myname << "Configuration parameters:" << endl;
    cout << myname << "      LogLevel: " << m_LogLevel << endl;
    cout << myname << "    NoiseSigma: " << m_NoiseSigma << endl;
    cout << myname << "   NSigmaStart: " << m_NSigmaStart << endl;
    cout << myname << "     NSigmaEnd: " << m_NSigmaEnd << endl;
    cout << myname << "   TicksBefore: " << m_TicksBefore << endl;
    cout << myname << "    TicksAfter: " << m_TicksAfter << endl;
  }
}

//**********************************************************************

DataMap DuneAdcSignalFinder::update(AdcChannelData& acd) const {
  const string myname = "DuneAdcSignalFinder::update: ";
  DataMap ret;
  if ( m_LogLevel >= 2 ) {
    cout << myname << "Finding ROIs for channel " << acd.channel << endl;
  }
  AdcIndex nsam = acd.samples.size();
  if ( nsam == 0 ) {
    cout << myname << "ERROR: No samples found in channel " << acd.channel << endl;
    acd.signal.clear();
    acd.rois.clear();
    return ret.setStatus(1);
  }
  AdcSignal sigma = m_NoiseSigma > 0 ? m_NoiseSigma : acd.sampleNoise;
  if ( sigma <= 0 ) {
    cout << myname << "ERROR: Unable to determine noise level." << endl;
    return ret.setStatus(2);
  }
  if ( m_LogLevel >= 3 ) cout << myname << "Noise level: " << sigma << endl;
  // Build ROIS before padding and merging.
  const AdcSignalVector& sams = acd.samples;
  AdcFilterVector& signal = acd.signal;
  AdcRoiVector& rois = acd.rois;
  signal.clear();
  signal.resize(sams.size(), false);
  bool inroi = false;
  AdcSignal siglow = m_NSigmaEnd*sigma;
  AdcSignal sighigh = m_NSigmaStart*sigma;
  for ( AdcIndex isam=0; isam<nsam; ++isam ) {
    AdcSignal sig = sams[isam];
    if ( inroi ) {
      if ( sig > siglow ) {
        signal[isam] = true;
      } else  {
        inroi = false;
      }
    } else {
      if ( sig > sighigh ) {
        signal[isam] = true;
        inroi = true;
      }
    }
  }
  // Fill the unpadded ROIs.
  acd.roisFromSignal();
  if ( rois.size() == 0 ) {
    if ( m_LogLevel >= 3 ) cout << myname << "  No ROIs found." << endl;
  } else {
    // Display ROIs before padding and merging.
    if ( m_LogLevel >= 3 ) {
      cout << myname << "  ROIs before merge (size = " << rois.size() << "):" << endl;
      for ( const AdcRoi& roi : rois ) {
        cout << myname << setw(8) << roi.first << " " << setw(8) << roi.second << endl;
      }
    } else if ( m_LogLevel >= 2 ) {
      cout << myname << "  ROI count before merge: " << acd.rois.size() << endl;
    }
    // Pad ROIs.
    unsigned int isam1 = 0;
    unsigned int isam2 = 0;
    for ( AdcRoi roi : rois ) {
      isam2 = roi.first;
      isam1 = isam2 > m_TicksBefore ? isam2 - m_TicksBefore : 0;
      for ( unsigned int isam=isam1; isam<isam2; ++isam ) signal[isam] = true;
      isam1 = roi.second + 1;
      isam2 = isam1 + m_TicksAfter;
      if ( isam2 > nsam ) isam2 = nsam;
      for ( unsigned int isam=isam1; isam<isam2; ++isam ) signal[isam] = true;
    }
    // Fill the final ROIs.
    acd.roisFromSignal();
    // Display final ROIs.
    if ( m_LogLevel >= 3 ) {
      cout << myname << "  ROIs after merge (size = " << rois.size() << "):" << endl;
      for ( const AdcRoi& roi : rois ) {
        cout << myname << setw(8) << roi.first << " " << setw(8) << roi.second << endl;
      }
    } else if ( m_LogLevel >= 2 ) {
      cout << myname << "  ROI count after merge: " << acd.rois.size() << endl;
    }
  }
  ret.setInt("nroi", acd.rois.size());
  return ret;
}

//**********************************************************************

DataMap DuneAdcSignalFinder::view(const AdcChannelData& acd) const {
  AdcChannelData acdtmp;
  acdtmp.samples = acd.samples;
  return update(acdtmp);
}

//**********************************************************************
