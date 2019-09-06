// AdcRoiSlicer_tool.cc

#include "AdcRoiSlicer.h"
#include <iostream>

using Name = AdcRoiSlicer::Name;
using Index = unsigned int;
using std::cout;
using std::endl;

//**********************************************************************
// Class methods.
//**********************************************************************

AdcRoiSlicer::AdcRoiSlicer(fhicl::ParameterSet const& ps)
: m_LogLevel(ps.get<int>("LogLevel")),
  m_OutViewName(ps.get<Name>("OutViewName")),
  m_SliceOpt(ps.get<int>("SliceOpt")),
  m_CopyRaw(ps.get<bool>("CopyRaw")) {
  const Name myname = "AdcRoiSlicer::ctor: ";
  if ( m_LogLevel >= 1 ) {
    cout << myname << "Configuration parameters:" << endl;
    cout << myname << "       LogLevel: " << m_LogLevel << endl;
    cout << myname << "    OutViewName: " << m_OutViewName << endl;
    cout << myname << "       SliceOpt: " << m_SliceOpt << endl;
    cout << myname << "        CopyRaw: " << (m_CopyRaw ? "true" : "false") << endl;
  }
}

//**********************************************************************

DataMap AdcRoiSlicer::update(AdcChannelData& acd) const {
  const Name myname = "AdcRoiSlicer::update: ";
  DataMap ret;
  bool keepRoi = m_SliceOpt == 1 || m_SliceOpt == 3;
  bool keepNot = m_SliceOpt == 2 || m_SliceOpt == 3;
  bool copyRaw = m_CopyRaw;
  if ( ! keepRoi && ! keepNot ) {
    cout << "ERROR: Invalid slice option: " << m_SliceOpt << endl;
    return ret.setStatus(1);
  }
  Index nsam = acd.samples.size();
  if ( acd.signal.size() != nsam ) {
    cout << "ERROR: Signal size does not match samples: " << acd.signal.size()
         << " != " << nsam << endl;
    return ret.setStatus(2);
  }
  if ( acd.hasView(m_OutViewName) ) {
    cout << "ERROR: Data for channel " << acd.channel << " already has view " << m_OutViewName << endl;
    return ret.setStatus(3);
  }
  if ( copyRaw && acd.raw.size() < nsam ) {
    cout << "ERROR: Insufficient raw data for channel " << acd.channel << ": "
         << acd.raw.size() << " < " << nsam << "." << endl;
    copyRaw = false;
  }
  AdcChannelData::View& view = acd.updateView(m_OutViewName);
  for ( Index isam=0; isam<nsam; ++isam ) {
    bool isRoi = acd.signal[isam];
    bool keep = (isRoi && keepRoi) || (!isRoi && keepNot);
    bool changeRoi = isam==0 || acd.signal[isam] != acd.signal[isam-1];
    bool startData = keep && changeRoi;
    if ( startData ) {
      view.push_back(acd);
      AdcChannelData& acdout = view.back();
      acdout.tick0 = isam;
    }
    if ( keep ) {
      AdcChannelData& acdout = view.back();
      if ( copyRaw ) acdout.raw.push_back(acd.raw[isam]);
      acdout.samples.push_back(acd.samples[isam]);
      acdout.signal.push_back(acd.signal[isam]);
    }
  }
  ret.setInt("nRoiView", view.size());
  return ret;
}

//**********************************************************************
