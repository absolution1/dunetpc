// AdcRoiViewer_tool.cc

#include "AdcRoiViewer.h"
#include "dune/DuneInterface/Tool/AdcChannelStringTool.h"
#include "dune/ArtSupport/DuneToolManager.h"
#include "dune/DuneCommon/gausTF1.h"
#include "dune/DuneCommon/coldelecResponse.h"
#include "dune/DuneCommon/quietHistFit.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <algorithm>
#include "TH1F.h"
#include "TDirectory.h"
#include "TFile.h"
#include "TF1.h"

using std::string;
using std::cout;
using std::endl;
using std::ostringstream;
using fhicl::ParameterSet;

using Index = AdcRoiViewer::Index;
using Name = AdcRoiViewer::Name;
using NameVector = std::vector<Name>;
using ParameterSetVector = std::vector<ParameterSet>;
using FloatVector = std::vector<float>;
using IntVector = std::vector<int>;

//**********************************************************************
// Subclass methods.
//**********************************************************************

AdcRoiViewer::State::~State() {
  for ( TH1*& ph : histTemplates ) {
    delete ph;
    ph = nullptr;
  }
  for ( HistMap::value_type ihst : hists ) {
    TH1*& ph = ihst.second;
    delete ph;
    ph = nullptr;
  }
}

//**********************************************************************

TH1* AdcRoiViewer::State::getHist(Name hname) {
  HistMap::iterator ihst = hists.find(hname);
  if ( ihst == hists.end() ) return nullptr;
  return ihst->second;
}

//**********************************************************************
// Class methods.
//**********************************************************************

AdcRoiViewer::AdcRoiViewer(fhicl::ParameterSet const& ps)
: m_LogLevel(ps.get<int>("LogLevel")),
  m_SigMinThresh(ps.get<float>("SigMinThresh")),
  m_SigMaxThresh(ps.get<float>("SigMaxThresh")),
  m_RoiHistOpt(ps.get<int>("RoiHistOpt")),
  m_FitOpt(ps.get<int>("FitOpt")),
  m_RoiRootFileName(ps.get<string>("RoiRootFileName")),
  m_SumRootFileName(ps.get<string>("SumRootFileName")),
  m_state(new AdcRoiViewer::State)
{
  const string myname = "AdcRoiViewer::ctor: ";
  string stringBuilder = "adcStringBuilder";
  DuneToolManager* ptm = DuneToolManager::instance();
  m_adcStringBuilder = ptm->getShared<AdcChannelStringTool>(stringBuilder);
  if ( m_adcStringBuilder == nullptr ) {
    cout << myname << "WARNING: AdcChannelStringTool not found: " << stringBuilder << endl;
  }
  ParameterSetVector pshists = ps.get<ParameterSetVector>("SumHists");
  for ( const ParameterSet& psh : pshists ) {
    Name hvar  = psh.get<Name>("var");
    Name hnam  = psh.get<Name>("name");
    Name httl  = psh.get<Name>("title");
    int nbin   = psh.get<int>("nbin");
    float xmin = psh.get<float>("xmin");
    float xmax = psh.get<float>("xmax");
    Name sfit;
    psh.get_if_present("fit", sfit);
    Name xlab = hvar;
    if      ( hvar == "fitHeight"    ) xlab = "Fit height% [SUNIT]%";
    else if ( hvar == "fitWidth"     ) xlab = "Fit width [Ticks]";
    else if ( hvar == "fitPosition"  ) xlab = "Fit position [Ticks]";
    else if ( hvar == "fitTickRem"   ) xlab = "Fit position remainder [Ticks]";
    else if ( hvar == "fitPeriodRem" )
      xlab = "Fit period " + std::to_string(m_TickPeriod) + " remainder [Ticks]";
    else if ( hvar == "fitChiSquare" ) xlab = "Fit #chi^{2}";
    else {
      cout << myname << "WARNING: Unknown summary variable: " << hvar << endl;
    }
    getState().vars.push_back(hvar);
    getState().fits.push_back(sfit);
    TH1* ph = new TH1F(hnam.c_str(), httl.c_str(), nbin, xmin, xmax);
    ph->SetDirectory(nullptr);
    ph->SetLineWidth(2);
    ph->GetXaxis()->SetTitle(xlab.c_str());
    ph->GetYaxis()->SetTitle("# ROI");
    // Add fit to template so it will be used for each child histogram.
    if ( sfit.size() ) {
      if ( m_LogLevel >= 1 ) cout << myname << "Adding fitter " << sfit
                                  << " to hist template " << hnam << endl;
      TF1* pf = new TF1(sfit.c_str(), sfit.c_str());
      ph->GetListOfFunctions()->AddLast(pf);
      ph->GetListOfFunctions()->SetOwner(kTRUE);
    }
    getState().histTemplates.push_back(ph);
  }
  if ( m_LogLevel>= 1 ) {
    cout << myname << "         LogLevel: " << m_LogLevel << endl;
    cout << myname << "       RoiHistOpt: " << m_RoiHistOpt << endl;
    cout << myname << "     SigMinThresh: " << m_SigMinThresh << endl;
    cout << myname << "     SigMaxThresh: " << m_SigMaxThresh << endl;
    cout << myname << "           FitOpt: " << m_FitOpt << endl;
    cout << myname << "  RoiRootFileName: " << m_RoiRootFileName << endl;
    cout << myname << "  SumRootFileName: " << m_SumRootFileName << endl;
    cout << myname << "         SumHists: [";
    for ( Index ihst=0; ihst<getState().vars.size(); ++ihst ) {
      if ( ihst ) cout << ", ";
      TH1* ph = getState().histTemplates[ihst];
      cout << ph->GetName() << "(" << getState().vars[ihst] << ")";
      if ( ph->GetListOfFunctions()->GetEntries() ) {
        cout << "-" << ph->GetListOfFunctions()->At(0)->GetName();
      }
    }
    cout << "]" << endl;
  }
}

//**********************************************************************

AdcRoiViewer::~AdcRoiViewer() {
  const string myname = "AdcRoiViewer::dtor: ";
  if ( m_LogLevel >= 1 ) cout << myname << "Exiting." << endl;
  fitSumHists();
  writeSumHists();
}

//**********************************************************************

DataMap AdcRoiViewer::view(const AdcChannelData& acd) const {
  DataMap res;
  doView(acd, m_LogLevel, res);
  if ( m_RoiRootFileName.size() ) {
    writeRoiHists(res, m_LogLevel);
  }
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
  bool save = m_RoiRootFileName.size();
  Index ndm = save ? acds.size() : 1;
  int dbg = m_LogLevel > 3 ? m_LogLevel - 2 : 0;
  Index nacd = acds.size();
  DataMapVector dms;  // Cache to hold results from doView
  dms.reserve(ndm);
  Index nroiLimit = save ? 1000 : 1000; // Clear cache after we get this many ROIs.
  Index nroiCached = 0;
  for ( AdcChannelDataMap::value_type iacd : acds ) {
    const AdcChannelData& acd = iacd.second;
    if ( m_LogLevel >= 3 ) {
      cout << myname << "Processing channel " << acd.channel
           << " (" << ncha << "/" << nacd << ")" << endl;
    }
    dms.emplace_back();
    DataMap& dm = dms.back();
    doView(acd, dbg, dm);
    if ( dm.status() ) {
      ++nfail;
      failedChannels.push_back(acd.channel);
    }
    ++ncha;
    nroi += dm.getInt("roiCount");
    nroiCached += dm.getInt("roiCount");
    if ( nroiCached > nroiLimit ) {
      if ( m_LogLevel >= 3 ) cout << myname << "  Clearing result cache." << endl;
      if ( save && dms.size() ) writeRoiHists(dms, dbg);
      dms.clear();
      nroiCached = 0;
    }
  }
  if ( save && dms.size() ) writeRoiHists(dms, dbg);
  ret.setInt("roiChannelCount", ncha);
  ret.setInt("roiFailedChannelCount", nfail);
  ret.setIntVector("roiFailedChannels", failedChannels);
  ret.setInt("roiCount", nroi);
  return ret;
}

//**********************************************************************

int AdcRoiViewer::doView(const AdcChannelData& acd, int dbg, DataMap& res) const {
  const string myname = "AdcRoiViewer::doView: ";
  unsigned int nraw = acd.raw.size();
  unsigned int nsam = acd.samples.size();
  unsigned int ntickChannel = nsam > nraw ? nsam : nraw;
  unsigned int nroiRaw = acd.rois.size();
  bool doHist = m_RoiHistOpt != 0;
  bool histRelativeTick = false;
  int histType = 0;
  if ( doHist ) {
    if        ( m_RoiHistOpt ==  1 ) {
      histType = 1;
    } else if ( m_RoiHistOpt ==  2 ) {
      histType = 2;
    } else if ( m_RoiHistOpt == 11 ) {
      histType = 1;
      histRelativeTick = true;
    } else if ( m_RoiHistOpt == 12 ) {
      histType = 2;
      histRelativeTick = true;
    } else {
      cout << myname << "Invalid value for RoiHistOpt: " << m_RoiHistOpt << endl;
      return res.setStatus(1).status();
    }
  }
  if ( dbg >=2 ) cout << myname << "Processing channel " << acd.channel << "."
                      << " Input ROI count is " << nroiRaw << endl;
  DataMap::HistVector roiHists;
  DataMap::FloatVector roiSigMins;
  DataMap::FloatVector roiSigMaxs;
  DataMap::FloatVector roiSigAreas;
  DataMap::FloatVector roiFitHeights;
  DataMap::FloatVector roiFitWidths;
  DataMap::FloatVector roiFitPositions;
  DataMap::FloatVector roiFitChiSquares;
  DataMap::IntVector roiTickMins;
  DataMap::IntVector roiTickMaxs;
  DataMap::IntVector roiNUnderflows;
  DataMap::IntVector roiNOverflows;
  DataMap::IntVector tick1;
  DataMap::IntVector ntick;
  DataMap::IntVector roiFitStats;
  Index nroi = 0;
  for ( unsigned int iroiRaw=0; iroiRaw<nroiRaw; ++iroiRaw ) {
    AdcRoi roi = acd.rois[iroiRaw];
    if ( dbg >=3 ) cout << myname << "  ROI " << nroi << "(raw " << iroiRaw << "): ["
                        << roi.first << ", " << roi.second << "]" << endl;
    ostringstream sshnam;
    sshnam << "hroi_evt%EVENT%_chan%CHAN%_roi";
    if ( nroi < 100 ) sshnam << "0";
    if ( nroi <  10 ) sshnam << "0";
    sshnam << nroi;
    string hnam = AdcChannelStringTool::AdcChannelStringTool::build(m_adcStringBuilder, acd, sshnam.str());
    ostringstream sshttl;
    sshttl << "Run %RUN% event %EVENT% channel %CHAN% ROI " << nroi;
    sshttl << " ;Tick ;";
    if ( histType == 1 ) sshttl << "Signal% [SUNIT]%";
    if ( histType == 2 ) sshttl << "ADC count";
    string httl = AdcChannelStringTool::build(m_adcStringBuilder, acd, sshttl.str());
    unsigned int isam1 = roi.first;
    unsigned int isam2 = roi.second + 1;
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
    Index nunder = 0;
    Index nover = 0;
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
      if ( flag == AdcUnderflow ) ++nunder;
      if ( flag == AdcOverflow ) ++nover;
    }
    // Check if this a ROI to keep.
    if ( m_SigMinThresh < 0.0 && sigmin > m_SigMinThresh ) continue;
    if ( m_SigMaxThresh > 0.0 && sigmax < m_SigMaxThresh ) continue;
    ++nroi;
    roiHists.push_back(ph);
    roiTickMins.push_back(roiTickMin);
    roiNUnderflows.push_back(nunder);
    roiNOverflows.push_back(nover);
    roiTickMaxs.push_back(roiTickMax);
    roiSigMins.push_back(sigmin);
    roiSigMaxs.push_back(sigmax);
    roiSigAreas.push_back(sigarea);
    tick1.push_back(isam1);
    ntick.push_back(isam2 - isam1);
    if ( m_FitOpt == 1 ) {
      if ( dbg >= 3 ) cout << "  Fitting with coldelecResponse" << endl;
      bool isNeg = fabs(sigmin) > sigmax;
      double h = isNeg ? sigmin : sigmax;
      //double shap = 2.5*ph->GetRMS();  // No! Negative entries break RMS calculation.
      double shap = 0.8*fabs(sigarea)/h;
      double t0 = x1 + (isNeg ? roiTickMin : roiTickMax) - shap;
      TF1* pf = coldelecResponseTF1(h, shap, t0, "coldlec");
      TF1* pfinit = coldelecResponseTF1(h, shap, t0, "coldlec");
      // The following was very slow when function was written to a file.
      //TF1* pfinit = dynamic_cast<TF1*>(pf->Clone("coldelec0"));
      pfinit->SetLineColor(3);
      pfinit->SetLineStyle(2);
      string fopt = "0";
      fopt = "WWB";
      if ( dbg < 3 ) fopt += "Q";
      int fstat = quietHistFit(ph, pf, fopt.c_str());
      ph->GetListOfFunctions()->AddLast(pfinit, "0");
      ph->GetListOfFunctions()->Last()->SetBit(TF1::kNotDraw, true);
      ph->GetListOfFunctions()->SetOwner(kTRUE);  // So the histogram owns pfinit
      roiFitHeights.push_back(pf->GetParameter(0));
      roiFitWidths.push_back(pf->GetParameter(1));
      roiFitPositions.push_back(pf->GetParameter(2));
      roiFitStats.push_back(fstat);
      roiFitChiSquares.push_back(pf->GetChisquare());
      delete pf;
      //delete pfinit;  This give error: list accessing deleted object
    }
  }
  res.setInt("roiEvent",   acd.event);
  res.setInt("roiRun",     acd.run);
  res.setInt("roiSubRun",  acd.subRun);
  res.setInt("roiChannel", acd.channel);
  res.setInt("roiCount", nroi);
  res.setInt("roiRawCount", nroiRaw);
  res.setInt("roiNTickChannel", ntickChannel);
  res.setIntVector("roiTick0s", tick1);
  res.setIntVector("roiNTicks", ntick);
  res.setIntVector("roiNUnderflows", roiNUnderflows);
  res.setIntVector("roiNOverflows", roiNOverflows);
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
    res.setIntVector("roiFitStats", roiFitStats);
    res.setFloatVector("roiFitChiSquares", roiFitChiSquares);
  }
  fillSumHists(acd, res);
  return res.status();
}

//**********************************************************************

void AdcRoiViewer::writeRoiHists(const DataMap& dm, int dbg) const {
  DataMapVector dms(1, dm);
  writeRoiHists(dms, dbg);
}

//**********************************************************************

void AdcRoiViewer::writeRoiHists(const DataMapVector& dms, int dbg) const {
  const string myname = "AdcRoiViewer::writeRoiHists: ";
  if ( m_RoiRootFileName.size() == 0 ) return;
  TDirectory* savdir = gDirectory;
  string ofrnameOld = "";
  TFile* pfile = nullptr;
  for ( const DataMap& dm : dms ) {
    AdcChannelData acd;
    acd.run     = dm.getInt("roiRun");
    acd.subRun  = dm.getInt("roiSubRun");
    acd.event   = dm.getInt("roiEvent");
    acd.channel = dm.getInt("roiChannel");
    string ofrname = AdcChannelStringTool::build(m_adcStringBuilder, acd, m_RoiRootFileName);
    if ( ofrname != ofrnameOld ) {
      if ( pfile != nullptr ) pfile->Close();
      delete pfile;
      if ( m_LogLevel >= 2 ) cout << myname << "Writing histograms to " << ofrname << endl;
      pfile = TFile::Open(ofrname.c_str(), "UPDATE");
      ofrnameOld = ofrname;
    }
    const DataMap::HistVector& roiHists = dm.getHistVector("roiHists");
    for ( TH1* ph : roiHists ) {
      TH1* phnew = dynamic_cast<TH1*>(ph->Clone());
      phnew->GetListOfFunctions()->SetOwner(kTRUE);  // So the histogram owns pfinit
      phnew->Write();
      if ( dbg >= 3 ) cout << myname << "  Wrote " << phnew->GetName() << endl;
    }
  }
  if ( pfile != nullptr ) pfile->Close();
  delete pfile;
  savdir->cd();
}

//**********************************************************************

void AdcRoiViewer::fillSumHists(const AdcChannelData acd, const DataMap& dm) const {
  const string myname = "AdcRoiViewer::fillSumHists: ";
  const NameVector vars = getState().vars;
  const HistVector histTemplates = getState().histTemplates;
  for ( Index ihst=0; ihst<vars.size(); ++ihst ) {
    Name var = vars[ihst];
    TH1* ph0 = histTemplates[ihst];
    FloatVector vals;
    IntVector ivals;
    if      ( var == "fitHeight"    ) vals = dm.getFloatVector("roiFitHeights");
    else if ( var == "fitWidth"     ) vals = dm.getFloatVector("roiFitWidths");
    else if ( var == "fitPosition"  ) vals = dm.getFloatVector("roiFitPositions");
    else if ( var == "fitTickRem"   ) vals = dm.getFloatVector("roiFitPositions");
    else if ( var == "fitPeriodRem" ) vals = dm.getFloatVector("roiFitPositions");
    else if ( var == "fitStat" )     ivals = dm.getIntVector("roiFitStats");
    else if ( var == "fitChiSquare" ) vals = dm.getFloatVector("roiFitChiSquares");
    else {
      if ( m_LogLevel >= 2 ) {
        cout << myname << "ERROR: Invalid variable name: " << var << endl;
        continue;
      }
    }
    Name hnam0 = ph0->GetName();
    Name hnam = AdcChannelStringTool::build(m_adcStringBuilder, acd, hnam0);
    Name xlab = AdcChannelStringTool::build(m_adcStringBuilder, acd, ph0->GetXaxis()->GetTitle());
    Name ylab = ph0->GetYaxis()->GetTitle();
    TH1* ph = getState().getHist(hnam);
    if ( ivals.size() && !vals.size() ) for ( int ival : ivals ) vals.push_back(ival);
    if ( var == "fitTickRem" ) for ( float& val : vals ) val = std::remainder(val,1);
    if ( var == "fitPeriodRem" ) for ( float& val : vals ) val = std::remainder(val,m_TickPeriod);
    if ( ph == nullptr && vals.size() ) {
      if ( m_LogLevel >= 2 ) cout << myname << "Creating histogram " << hnam << endl;
      Name httl0 = ph0->GetTitle();
      Name httl = AdcChannelStringTool::build(m_adcStringBuilder, acd, httl0);
      int nbin = ph0->GetNbinsX();
      float xmin = ph0->GetXaxis()->GetXmin();
      float xmax = ph0->GetXaxis()->GetXmax();
      // If xmin > xmax and xmin > 0, then we center histogram on median and use width = xmin.
      // If also xmax >0, then we round the first bin edge to that value.
      if ( xmin > xmax && xmin > 0.0 ) {
        FloatVector tmpvals = vals;
        std::sort(tmpvals.begin(), tmpvals.end());
        Index nval = tmpvals.size();
        if ( m_LogLevel >= 3 ) cout << myname << "  Centering histogram on median of "
                                    << nval << " value" << (nval==1 ? "" : "s") << endl;
        float xmed = 0.5*(tmpvals[(nval-1)/2] + tmpvals[nval/2]);
        float width = xmin;
        xmin = xmed - 0.5*width;
        bool roundXmin = xmax > 0.0;
        if ( roundXmin ) {
          //float rexp = log10(width/50.0);
          //rexp = rexp > 0 ? int(rexp) : int(rexp-1);
          //float rfac = pow(10, rexp);
          float rfac = xmax;
          xmin = rfac*int(xmin/rfac + (xmin > 0.0 ? 0.5 : -0.5));
        }
        xmax = xmin + width;
      }
      if ( m_LogLevel >= 3 ) {
        cout << myname << "   Name: " << hnam << endl;
        cout << myname << "  Title: " << httl << endl;
        cout << myname << "   nbin: " << nbin << endl;
        cout << myname << "   xmin: " << xmin << endl;
        cout << myname << "   xmax: " << xmax << endl;
      }
      ph = new TH1F(hnam.c_str(), httl.c_str(), nbin, xmin, xmax);
      ph->SetDirectory(nullptr);
      ph->SetStats(0);
      ph->SetLineWidth(2);
      ph->GetXaxis()->SetTitle(xlab.c_str());
      ph->GetYaxis()->SetTitle(ylab.c_str());
      if ( ph0->GetListOfFunctions()->GetEntries() ) {
        TF1* pf = dynamic_cast<TF1*>(ph0->GetListOfFunctions()->At(0)->Clone());
        ph->GetListOfFunctions()->AddLast(pf);
        ph->GetListOfFunctions()->SetOwner(kTRUE);
      }
      getState().hists[hnam] = ph;
    }
    if ( m_LogLevel >= 3 ) cout << myname << "Filling histogram " << hnam << endl;
    for ( float val : vals ) {
      ph->Fill(val);
    }
  }
}

//**********************************************************************

void AdcRoiViewer::fitSumHists() const {
  const string myname = "AdcRoiViewer::fitSumHists: ";
  if ( getState().hists.size() == 0 ) {
    cout << myname << "No summary histograms found." << endl;
    return;
  }
  if ( m_LogLevel >= 1 ) cout << myname << "Fitting summary histograms. Count is "
                              << getState().hists.size() << "." << endl;
  for ( HistMap::value_type ihst : getState().hists ) {
    TH1* ph = ihst.second;
    Index nfun = ph->GetListOfFunctions()->GetEntries();
    if ( nfun ) {
      Name fname = ph->GetListOfFunctions()->At(0)->GetName();
      if ( m_LogLevel >= 3 ) cout << myname << "Fitting hist " << ph->GetName() << " with " << fname << endl;
      TF1* pf = nullptr;
      if ( fname == "gaus" ) {
        double mean = ph->GetMean();
        double sigma = ph->GetRMS();
        double height = ph->GetMaximum();
        pf = gausTF1(height, mean, sigma);
      } else {
        pf = new TF1(fname.c_str(), fname.c_str());
      }
      if ( m_LogLevel >= 4 ) cout << myname << "  Created function " << pf->GetName() << " at " << std::hex << pf << endl;
      int fstat = quietHistFit(ph, pf, "WWQ");
      if ( fstat != 0 ) {
        cout << myname << "  WARNING: Fit " << pf->GetName() << " of " << ph->GetName() << " returned " << fstat << endl;
      } else {
        if ( m_LogLevel >=4 ) cout << myname << "  Fit succeeded." << endl;
      }
      ph->GetListOfFunctions()->SetOwner(kTRUE);  // So the histogram owns pf
      delete pf;
    }
  }
}

//**********************************************************************

void AdcRoiViewer::writeSumHists() const {
  const string myname = "AdcRoiViewer::writeSumHists: ";
  if ( m_SumRootFileName.size() == 0 ) return;
  if ( getState().hists.size() == 0 ) {
    cout << myname << "No summary histograms found." << endl;
    return;
  }
  AdcChannelData acd;
  Name ofrname = AdcChannelStringTool::build(m_adcStringBuilder, acd, m_SumRootFileName);
  TDirectory* savdir = gDirectory;
  TFile* pfile = TFile::Open(ofrname.c_str(), "UPDATE");
  if ( m_LogLevel >= 1 ) cout << myname << "Writing summary histograms. Count is "
                              << getState().hists.size() << "." << endl;
  for ( HistMap::value_type ihst : getState().hists ) {
    TH1* ph = ihst.second;
    TH1* phnew = dynamic_cast<TH1*>(ph->Clone());
    phnew->GetListOfFunctions()->SetOwner(kTRUE);
    phnew->Write();
    if ( m_LogLevel >= 2 ) cout << myname << "  Wrote " << phnew->GetName() << endl;
  }
  if ( pfile != nullptr ) pfile->Close();
  delete pfile;
  if ( m_LogLevel >= 1 ) cout << myname << "Closed summary histogram file " << ofrname << endl;
  savdir->cd();
}

//**********************************************************************
