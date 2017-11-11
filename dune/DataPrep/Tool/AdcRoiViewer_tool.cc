// AdcRoiViewer_tool.cc

#include "AdcRoiViewer.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include "TH1F.h"

using std::string;
using std::cout;
using std::endl;
using std::ostringstream;

//**********************************************************************
// Class methods.
//**********************************************************************

AdcRoiViewer::AdcRoiViewer(fhicl::ParameterSet const& ps)
: m_LogLevel(ps.get<int>("LogLevel")),
  m_HistOpt(ps.get<int>("HistOpt")) {
  const string myname = "AdcRoiViewer::ctor: ";
  if ( m_LogLevel>= 1 ) {
    cout << myname << "  LogLevel: " << m_LogLevel << endl;
    cout << myname << "   HistOpt: " << m_HistOpt << endl;
  }
}

//**********************************************************************

DataMap AdcRoiViewer::view(const AdcChannelData& acd) const {
  const string myname = "AdcRoiViewer::view: ";
  DataMap res;
  unsigned int nraw = acd.raw.size();
  unsigned int nsam = acd.samples.size();
  unsigned int ntickChannel = nsam > nraw ? nsam : nraw;
  unsigned int nroi = acd.rois.size();
  bool doHist = m_HistOpt != 0;
  bool histRelativeTick = false;
  int histType = 0;
  if ( doHist ) {
    if        ( m_HistOpt ==  1 ) {
      histType = 1;
    } else if ( m_HistOpt ==  2 ) {
      histType = 2;
    } else if ( m_HistOpt == 11 ) {
      histType = 1;
      histRelativeTick = true;
    } else if ( m_HistOpt == 12 ) {
      histType = 2;
      histRelativeTick = true;
    } else {
      cout << myname << "Invalid value for HistOpt: " << m_HistOpt << endl;
      return res.setStatus(1);
    }
  }
  DataMap::HistVector roiHists;
  DataMap::FloatVector roiSigMins;
  DataMap::FloatVector roiSigMaxs;
  DataMap::FloatVector roiSigAreas;
  DataMap::IntVector nUnderflow(nroi, 0);
  DataMap::IntVector nOverflow(nroi, 0);
  DataMap::IntVector tick1(nroi, 0);
  DataMap::IntVector ntick(nroi, 0);
  for ( unsigned int iroi=0; iroi<nroi; ++iroi ) {
    AdcRoi roi = acd.rois[iroi];
    ostringstream sshnam;
    sshnam << "hroi";
    if ( iroi < 100 ) sshnam << "0";
    if ( iroi <  10 ) sshnam << "0";
    sshnam << iroi;
    string hnam = sshnam.str();
    ostringstream sshttl;
    sshttl << "ROI " << iroi;
    sshttl << " ;Tick ;Signal";
    string httl = sshttl.str();
    unsigned int isam1 = roi.first;
    unsigned int isam2 = roi.second + 1;
    tick1[iroi] = isam1;
    ntick[iroi] = isam2 - isam1;
    float x1 = histRelativeTick ? 0.0 : isam1;
    float x2 = histRelativeTick ? isam2 - isam1 : isam2;
    TH1* ph = new TH1F(hnam.c_str(), httl.c_str(), isam2-isam1, x1, x2);
    ph->SetDirectory(nullptr);
    ph->SetStats(0);
    ph->SetLineWidth(2);
    unsigned int ibin = 0;
    float sigmin = 0.0;
    float sigmax = 0.0;
    float sigarea = 0.0;
    for ( unsigned int isam=isam1; isam<isam2; ++isam ) {
      float sig = 0.0;
      if ( histType == 1 && isam<nsam ) sig = acd.samples[isam];
      if ( histType == 2 && isam<nraw ) sig = isam<nraw ? acd.raw[isam] : 0.0;
      if ( ibin == 0 ) {
        sigmin = sig;
        sigmax = sig;
      } else {
        if ( sig < sigmin ) sigmin = sig;
        if ( sig > sigmax ) sigmax = sig;
      }
      sigarea += sig;
      ph->SetBinContent(++ibin, sig);
      AdcFlag flag = acd.flags.size() > isam ? acd.flags[isam] : 0;
      if ( flag == AdcUnderflow ) ++nUnderflow[iroi];
      if ( flag == AdcOverflow ) ++nOverflow[iroi];
    }
    roiHists.push_back(ph);
    roiSigMins.push_back(sigmin);
    roiSigMaxs.push_back(sigmax);
    roiSigAreas.push_back(sigarea);
  }
  res.setInt("roiCount", nroi);
  res.setInt("roiNTickChannel", ntickChannel);
  res.setIntVector("roiTick0s", tick1);
  res.setIntVector("roiNTicks", ntick);
  res.setIntVector("roiNUnderflows", nUnderflow);
  res.setIntVector("roiNOverflows", nOverflow);
  res.setFloatVector("roiSigMins", roiSigMins);
  res.setFloatVector("roiSigMaxs", roiSigMaxs);
  res.setFloatVector("roiSigAreas", roiSigAreas);
  res.setHistVector("roiHists", roiHists, true);
  return res;
}

//**********************************************************************
