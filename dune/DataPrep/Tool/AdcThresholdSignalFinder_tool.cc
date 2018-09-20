// AdcThresholdSignalFinder_tool.cc

#include "AdcThresholdSignalFinder.h"
#include <iostream>

using std::string;
using std::cout;
using std::endl;

//**********************************************************************
// Local defintitions.
//**********************************************************************

namespace {

string boolToString(bool val) {
  return val ? "true" : "false";
}

}  // end unnamed namespace

//**********************************************************************
// Class methods.
//**********************************************************************

AdcThresholdSignalFinder::AdcThresholdSignalFinder(fhicl::ParameterSet const& ps)
: m_LogLevel(ps.get<int>("LogLevel")),
  m_Threshold(ps.get<float>("Threshold")),
  m_BinsBefore(ps.get<unsigned int>("BinsBefore")),
  m_BinsAfter(ps.get<unsigned int>("BinsAfter")),
  m_FlagPositive(ps.get<bool>("FlagPositive")),
  m_FlagNegative(ps.get<bool>("FlagNegative")) {
  const string myname = "AdcThresholdSignalFinder::ctor: ";
  if ( m_LogLevel >= 1 ) {
    cout << myname << "Configuration parameters:" << endl;
    cout << myname << "      LogLevel: " << m_LogLevel << endl;
    cout << myname << "     Threshold: " << m_Threshold<< endl;
    cout << myname << "    BinsBefore: " << m_BinsBefore << endl;
    cout << myname << "     BinsAfter: " << m_BinsAfter << endl;
    cout << myname << "  FlagPositive: " << boolToString(m_FlagPositive) << endl;
    cout << myname << "  FlagNegative: " << boolToString(m_FlagNegative) << endl;
  }
}

//**********************************************************************

DataMap AdcThresholdSignalFinder::update(AdcChannelData& acd) const {
  const string myname = "AdcThresholdSignalFinder::update: ";
  DataMap ret;
  AdcIndex nsam = acd.samples.size();
  if ( nsam == 0 ) {
    cout << myname << "ERROR: No samples found in channel " << acd.channel << endl;
    acd.signal.clear();
    acd.rois.clear();
    return ret.setStatus(1);
  }
  if ( m_LogLevel >= 2 ) cout << myname << "Finding ROIs for channel " << acd.channel << endl;
  AdcIndex nsamlo = m_BinsBefore;
  AdcIndex nsamhi = m_BinsAfter;
  acd.signal.resize(nsam, false);
  AdcIndex nbinAbove = 0;
  AdcIndex isamNotRoi = 0;  // This is the 1st sample after the last ROI.
  for ( AdcIndex isam=0; isam<nsam; ++isam ) {
    bool keep = ( m_FlagPositive && acd.samples[isam] >  m_Threshold ) ||
                ( m_FlagNegative && acd.samples[isam] < -m_Threshold );
    if ( keep ) { // Leave isam after the new ROI
      ++nbinAbove;
      AdcIndex jsam1 = isam > nsamlo ? isam - nsamlo : 0;
      if ( jsam1 < isamNotRoi ) jsam1 = isamNotRoi;
      AdcIndex jsam2 = isam + nsamhi + 1;
      if ( jsam2 > nsam ) jsam2 = nsam;
      for ( isam=jsam1; isam<jsam2; ++isam ) acd.signal[isam] = true;
      isamNotRoi = isam;
    } else {
      ++isam;
    }
  }
  acd.roisFromSignal();
  if ( m_LogLevel >= 3 ) {
    cout << myname << "  # ticks above threshold: " << nbinAbove << endl;
    cout << myname << "  # ROI: " << acd.rois.size() << endl;
  }
  ret.setInt("nThresholdBins", nbinAbove);
  ret.setInt("nroi", acd.rois.size());
  return ret;
}

//**********************************************************************

DataMap AdcThresholdSignalFinder::view(const AdcChannelData& acd) const {
  AdcChannelData acdtmp;
  acdtmp.samples = acd.samples;
  return update(acdtmp);
}

//**********************************************************************
