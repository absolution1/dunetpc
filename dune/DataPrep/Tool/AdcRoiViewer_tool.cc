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
: m_LogLevel(ps.get<int>("LogLevel")) {
  const string myname = "AdcRoiViewer::ctor: ";
  if ( m_LogLevel>= 1 ) {
    cout << myname << "  LogLevel: " << m_LogLevel << endl;
  }
}

//**********************************************************************

DataMap AdcRoiViewer::view(const AdcChannelData& acd) const {
  DataMap res;
  unsigned int nsam = acd.samples.size();
  unsigned int nroi = acd.rois.size();
  res.setInt("nROI", nroi);
  DataMap::HistVector roiHists;
  DataMap::FloatVector roiSigMins;
  DataMap::FloatVector roiSigMaxs;
  DataMap::FloatVector roiSigAreas;
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
    TH1* ph = new TH1F(hnam.c_str(), httl.c_str(), isam2-isam1, isam1, isam2);
    ph->SetDirectory(nullptr);
    ph->SetStats(0);
    ph->SetLineWidth(2);
    unsigned int ibin = 0;
    float sigmin = 0.0;
    float sigmax = 0.0;
    float sigarea = 0.0;
    for ( unsigned int isam=isam1; isam<isam2; ++isam ) {
      float sam = isam<nsam ? acd.samples[isam] : 0.0;
      if ( ibin == 0 ) {
        sigmin = sam;
        sigmax = sam;
      } else {
        if ( sam < sigmin ) sigmin = sam;
        if ( sam > sigmax ) sigmax = sam;
      }
      sigarea += sam;
      ph->SetBinContent(++ibin, sam);
    }
    roiHists.push_back(ph);
    roiSigMins.push_back(sigmin);
    roiSigMaxs.push_back(sigmax);
    roiSigAreas.push_back(sigarea);
  }
  res.setHistVector("roiHists", roiHists, true);
  res.setFloatVector("roiSigMins", roiSigMins);
  res.setFloatVector("roiSigMaxs", roiSigMaxs);
  res.setFloatVector("roiSigAreas", roiSigAreas);
  return res;
}

//**********************************************************************
