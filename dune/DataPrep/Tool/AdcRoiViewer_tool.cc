// AdcRoiViewer_tool.cc

#include "AdcRoiViewer.h"
#include "dune/DuneInterface/Tool/AdcChannelStringTool.h"
#include "dune/ArtSupport/DuneToolManager.h"
#include "dune/DuneCommon/coldelecResponse.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include "TH1F.h"
#include "TDirectory.h"
#include "TFile.h"
#include "TF1.h"

using std::string;
using std::cout;
using std::endl;
using std::ostringstream;

using Index = unsigned int;

//**********************************************************************
// Class methods.
//**********************************************************************

AdcRoiViewer::AdcRoiViewer(fhicl::ParameterSet const& ps)
: m_LogLevel(ps.get<int>("LogLevel")),
  m_HistOpt(ps.get<int>("HistOpt")),
  m_FitOpt(ps.get<int>("FitOpt")),
  m_RootFileName(ps.get<string>("RootFileName"))
{
  const string myname = "AdcRoiViewer::ctor: ";
  string stringBuilder = "adcStringBuilder";
  DuneToolManager* ptm = DuneToolManager::instance();
  m_adcStringBuilder = ptm->getShared<AdcChannelStringTool>(stringBuilder);
  if ( m_adcStringBuilder == nullptr ) {
    cout << myname << "WARNING: AdcChannelStringTool not found: " << stringBuilder << endl;
  }
  if ( m_LogLevel>= 1 ) {
    cout << myname << "      LogLevel: " << m_LogLevel << endl;
    cout << myname << "       HistOpt: " << m_HistOpt << endl;
    cout << myname << "        FitOpt: " << m_FitOpt << endl;
    cout << myname << "  RootFileName: " << m_RootFileName << endl;
  }
}

//**********************************************************************

DataMap AdcRoiViewer::view(const AdcChannelData& acd) const {
  DataMap res;
  doView(acd, m_LogLevel, res);
  return res;
}

//**********************************************************************

DataMap AdcRoiViewer::viewMap(const AdcChannelDataMap& acds) const {
  const string myname = "AdcRoiViewer::viewMap: ";
  DataMap ret;
  Index ncha = 0;
  Index nroi = 0;
  Index nfail = 0;
  DataMap::IntVector failedChannels;
  for ( AdcChannelDataMap::value_type iacd : acds ) {
    const AdcChannelData acd = iacd.second;
    if ( m_LogLevel >= 3 ) cout << myname << "Processing channel " << acd.channel << endl;
    DataMap dm;
    int dbg = m_LogLevel > 3 ? m_LogLevel - 2 : 0;
    doView(acd, dbg, dm);
    if ( dm.status() ) {
      ++nfail;
      failedChannels.push_back(acd.channel);
    }
    ++ncha;
     nroi += dm.getInt("roiCount");
  }
  ret.setInt("roiChannelCount", ncha);
  ret.setInt("roiFailedChannelCount", nfail);
  ret.setIntVector("roiFailedChannels", failedChannels);
  ret.setInt("roiCount", nroi);
  return ret;
}

//**********************************************************************

int AdcRoiViewer::doView(const AdcChannelData& acd, int dbg, DataMap& res) const {
  const string myname = "AdcRoiViewer::view: ";
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
      return res.setStatus(1).status();
    }
  }
  if ( dbg >=2 ) cout << myname << "Processing channel " << acd.channel << ". ROI count is " << nroi << endl;
  DataMap::HistVector roiHists;
  DataMap::FloatVector roiSigMins;
  DataMap::FloatVector roiSigMaxs;
  DataMap::FloatVector roiSigAreas;
  DataMap::FloatVector roiFitHeights;
  DataMap::FloatVector roiFitWidths;
  DataMap::FloatVector roiFitPositions;
  DataMap::IntVector roiTickMins;
  DataMap::IntVector roiTickMaxs;
  DataMap::IntVector nUnderflow(nroi, 0);
  DataMap::IntVector nOverflow(nroi, 0);
  DataMap::IntVector tick1(nroi, 0);
  DataMap::IntVector ntick(nroi, 0);
  for ( unsigned int iroi=0; iroi<nroi; ++iroi ) {
    AdcRoi roi = acd.rois[iroi];
    if ( dbg >=3 ) cout << myname << "  ROI " << iroi << ": ["
                        << roi.first << ", " << roi.second << "]" << endl;
    
    ostringstream sshnam;
    sshnam << "hroi_evt%EVENT%_chan%CHAN%_roi";
    if ( iroi < 100 ) sshnam << "0";
    if ( iroi <  10 ) sshnam << "0";
    sshnam << iroi;
    string hnam = AdcChannelStringTool::AdcChannelStringTool::build(m_adcStringBuilder, acd, sshnam.str());
    ostringstream sshttl;
    sshttl << "Run %RUN% event %EVENT% channel %CHAN% ROI " << iroi;
    sshttl << " ;Tick ;";
    if ( histType == 1 ) sshttl << "Signal% [SUNIT]%";
    if ( histType == 2 ) sshttl << "ADC count";
    string httl = AdcChannelStringTool::build(m_adcStringBuilder, acd, sshttl.str());
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
    int roiTickMin = 0;
    int roiTickMax = 0;
    for ( unsigned int isam=isam1; isam<isam2; ++isam ) {
      float sig = 0.0;
      if ( histType == 1 && isam<nsam ) sig = acd.samples[isam];
      if ( histType == 2 && isam<nraw ) sig = isam<nraw ? acd.raw[isam] : 0.0;
      if ( ibin == 0 ) {
        sigmin = sig;
        sigmax = sig;
      } else {
        if ( sig < sigmin ) {
          sigmin = sig;
          roiTickMin = ibin;
        }
        if ( sig > sigmax ) {
          sigmax = sig;
          roiTickMax = ibin;
        }
      }
      sigarea += sig;
      ph->SetBinContent(++ibin, sig);
      AdcFlag flag = acd.flags.size() > isam ? acd.flags[isam] : 0;
      if ( flag == AdcUnderflow ) ++nUnderflow[iroi];
      if ( flag == AdcOverflow ) ++nOverflow[iroi];
    }
    roiHists.push_back(ph);
    roiTickMins.push_back(roiTickMin);
    roiTickMaxs.push_back(roiTickMax);
    roiSigMins.push_back(sigmin);
    roiSigMaxs.push_back(sigmax);
    roiSigAreas.push_back(sigarea);
    if ( m_FitOpt == 1 ) {
      if ( dbg >= 3 ) cout << "  Fitting with coldelecResponse" << endl;
      bool isNeg = fabs(sigmin) > sigmax;
      double h = isNeg ? sigmin : sigmax;
      double shap = 2.5*ph->GetRMS();
      double t0 = x1 + (isNeg ? roiTickMin : roiTickMax) - shap;
      TF1* pf = coldelecResponseTF1(h, shap, t0, "coldlec");
      TF1* pfinit = dynamic_cast<TF1*>(pf->Clone("coldelec0"));
      pfinit->SetLineColor(3);
      pfinit->SetLineStyle(2);
      string fopt = "0";
      fopt = "WWB";
      if ( dbg < 3 ) fopt += "Q";
      // Block Root info message for new Canvas produced in fit.
      int levelSave = gErrorIgnoreLevel;
      gErrorIgnoreLevel = 1001;
      // Block non-default (e.g. art) from handling the Root "error".
      // We switch to the Root default handler while making the call to Print.
      ErrorHandlerFunc_t pehSave = nullptr;
      ErrorHandlerFunc_t pehDefault = DefaultErrorHandler;
      if ( GetErrorHandler() != pehDefault ) {
        pehSave = SetErrorHandler(pehDefault);
      }
      ph->Fit(pf, fopt.c_str());
      if ( pehSave != nullptr ) SetErrorHandler(pehSave);
      gErrorIgnoreLevel = levelSave;
      ph->GetListOfFunctions()->AddLast(pfinit, "0");
      ph->GetListOfFunctions()->Last()->SetBit(TF1::kNotDraw, true);
      roiFitHeights.push_back(pf->GetParameter(0));
      roiFitWidths.push_back(pf->GetParameter(1));
      roiFitPositions.push_back(pf->GetParameter(2));
      delete pf;
    }
  }
  res.setInt("roiCount", nroi);
  res.setInt("roiNTickChannel", ntickChannel);
  res.setIntVector("roiTick0s", tick1);
  res.setIntVector("roiNTicks", ntick);
  res.setIntVector("roiNUnderflows", nUnderflow);
  res.setIntVector("roiNOverflows", nOverflow);
  res.setIntVector("roiTickMins", roiTickMins);
  res.setIntVector("roiTickMaxs", roiTickMaxs);
  res.setFloatVector("roiSigMins", roiSigMins);
  res.setFloatVector("roiSigMaxs", roiSigMaxs);
  res.setFloatVector("roiSigAreas", roiSigAreas);
  res.setHistVector("roiHists", roiHists, true);
  if ( roiFitHeights.size() ) {
    res.setFloatVector("roiFitHeights", roiFitHeights);
    res.setFloatVector("roiFitWidths", roiFitWidths);
    res.setFloatVector("roiFitPositions", roiFitPositions);
  }
  if ( m_RootFileName.size () ) {
    TDirectory* savdir = gDirectory;
    string ofrname = AdcChannelStringTool::build(m_adcStringBuilder, acd, m_RootFileName);
    if ( dbg >= 2 ) cout << myname << "Writing histograms to " << ofrname << endl;
    TFile* pfile = TFile::Open(ofrname.c_str(), "UPDATE");
    for ( TH1* ph : roiHists ) {
      TH1* phnew = dynamic_cast<TH1*>(ph->Clone());
      phnew->Write();
      if ( dbg >= 3 ) cout << myname << "  Wrote " << phnew->GetName() << endl;
    }
    delete pfile;
    savdir->cd();
  }
  return res.status();
}

//**********************************************************************
