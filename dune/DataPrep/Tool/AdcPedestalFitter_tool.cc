// AdcPedestalFitter_tool.cc

#include "AdcPedestalFitter.h"
#include "dune/DuneCommon/StringManipulator.h"
#include "dune/DuneInterface/Tool/HistogramManager.h"
#include "dune/ArtSupport/DuneToolManager.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include "TDirectory.h"
#include "TFile.h"
#include "TH1F.h"
#include "TF1.h"
#include "TROOT.h"

using std::string;
using std::cout;
using std::endl;
using std::ofstream;
using std::ostream;
using std::ostringstream;
using std::setw;
using std::fixed;
using std::setprecision;
using std::vector;

using Index = unsigned int;

//**********************************************************************
// Class methods.
//**********************************************************************

AdcPedestalFitter::AdcPedestalFitter(fhicl::ParameterSet const& ps)
: m_LogLevel(ps.get<int>("LogLevel")),
  m_HistName(ps.get<string>("HistName")),
  m_HistTitle(ps.get<string>("HistTitle")),
  m_HistManager(ps.get<string>("HistManager")),
  m_phm(nullptr) {
  const string myname = "AdcPedestalFitter::ctor: ";
  if ( m_HistManager.size() ) {
    DuneToolManager* ptm = DuneToolManager::instance();
    m_phm = ptm->getShared<HistogramManager>(m_HistManager);
    if ( m_phm == nullptr ) {
      cout << myname << "WARNING: Histogram manager not found: " << m_HistManager << endl;
    }
  }
  if ( m_LogLevel >= 1 ) {
    cout << myname << "Configuration parameters:" << endl;
    cout << myname << "    LogLevel: " << m_LogLevel << endl;
    cout << myname << "    HistName: " << m_HistName << endl;
    cout << myname << "   HistTitle: " << m_HistTitle << endl;
    cout << myname << " HistManager: " << m_HistManager << endl;
  }
}

//**********************************************************************

DataMap AdcPedestalFitter::view(const AdcChannelData& acd) const {
  const string myname = "AdcPedestalFitter::view: ";
  if ( m_LogLevel >= 3 ) cout << myname << "Calling " << endl;
  DataMap res = getPedestal(acd);
  TH1* phped = res.getHist("pedestal");
  if ( m_LogLevel >= 3 ) cout << myname << "Viewing " << endl;
  if ( res.status() != 0 ) {
    delete phped;
    return res;
  }
  if ( m_phm == nullptr ) {
    res.setHist("pedestal", phped, true);
  } else {
    int rstat = m_phm->manage(phped);
    if ( rstat != 0 ) {
      cout << myname << "WARNING: Attempt to manage histogram " << phped->GetName()
           << " returned error " << rstat << endl;
      delete phped;
      res.setHist("pedestal", nullptr);
      return res.setStatus(11);
    }
  }
  if ( m_LogLevel >= 3 ) cout << myname << "Exiting" << endl;
  return res;
}

//**********************************************************************

DataMap AdcPedestalFitter::update(AdcChannelData& acd) const {
  const string myname = "AdcPedestalFitter::update: ";
  if ( m_LogLevel >= 3 ) cout << myname << "Calling " << endl;
  DataMap res = getPedestal(acd);
  if ( m_LogLevel >= 3 ) cout << myname << "Called " << endl;
  TH1* phped = res.getHist("pedestal");
  delete phped;
  res.setHist("pedestal", nullptr);
  if ( res.status() != 0 ) return res;
  if ( m_LogLevel >= 3 ) cout << myname << "Old pedestal: " << acd.pedestal << endl;
  acd.pedestal = res.getFloat("fitPedestal");
  if ( m_LogLevel >= 3 ) cout << myname << "New pedestal: " << acd.pedestal << endl;
  return res;
}

//**********************************************************************

string AdcPedestalFitter::
nameReplace(string name, const AdcChannelData& acd) const {
  string nameout = name;
  StringManipulator sman(nameout);
  if ( acd.run != AdcChannelData::badIndex ) sman.replace("%RUN%", acd.run);
  else sman.replace("%RUN%", "RunNotFound");
  if ( acd.subRun != AdcChannelData::badIndex ) sman.replace("%SUBRUN%", acd.subRun);
  else sman.replace("%SUBRUN%", "SubRunNotFound");
  if ( acd.event != AdcChannelData::badIndex ) sman.replace("%EVENT%", acd.event);
  else sman.replace("%EVENT%", "EventNotFound");
  if ( acd.channel != AdcChannelData::badChannel ) sman.replace("%CHAN%", acd.channel);
  else sman.replace("%CHAN%", "ChannelNotFound");
  return nameout;
}

//**********************************************************************

DataMap
AdcPedestalFitter::getPedestal(const AdcChannelData& acd) const {
  const string myname = "AdcPedestalFitter::getPedestal: ";
  DataMap res;
  if ( m_LogLevel >= 3 ) cout << myname << "Entering..." << endl;
  string hnameBase = m_HistName;
  string htitlBase = m_HistTitle;
  Index nsam = acd.raw.size();
  if ( nsam == 0 ) {
    if ( m_LogLevel >= 2 ) cout << myname << "WARNING: Raw data is empty." << endl;
    return res.setStatus(1);
  }
  string hname = nameReplace(hnameBase, acd);
  string htitl = nameReplace(htitlBase, acd);
  htitl += "; ADC count; # samples";
  unsigned int nadc = 4096;
  unsigned int rebin = 10;
  unsigned int nbin = (nadc + rebin - 0.01)/rebin;
  double xmax = rebin*nbin;
  string hnamr = hname + "_rebin";
  TH1* phr = new TH1F(hname.c_str(), htitl.c_str(), nbin, 0, xmax);
  phr->SetDirectory(0);
  for ( Index isam=0; isam<nsam; ++isam ) {
    phr->Fill(acd.raw[isam]);
  }
  int rbinmax1 = phr->GetMaximumBin();
  double radcmax1 = phr->GetBinCenter(rbinmax1);
  double adcmean = phr->GetMean();
  // Max may just be a sticky code. Reduce it and find the next maximum.
  double tmpval = 0.5*(phr->GetBinContent(rbinmax1-1)+phr->GetBinContent(rbinmax1+1));
  phr->SetBinContent(rbinmax1, tmpval);
  int rbinmax2 = phr->GetMaximumBin();
  // Define the max to be the first value if the two maxima are close or the
  // average if they are far part.
  int rbinmax = rbinmax1;
  if ( abs(rbinmax2-rbinmax1) > 1 ) {
    rbinmax = (rbinmax1 + rbinmax2)/2;
    adcmean = phr->GetMean();
  }
  double adcmax = phr->GetBinCenter(rbinmax);
  delete phr;
  double wadc = 100.0;
  // Make sure the peak bin stays in range.
  if ( abs(adcmax-radcmax1) > 0.5*wadc ) adcmax = radcmax1;
  double adc1 = adcmax - 0.5*wadc;
  adc1 = 10*int(adc1/10);
  if ( adcmean > adcmax + 10) adc1 += 10;
  double adc2 = adc1 + wadc;
  TH1* phf = new TH1F(hname.c_str(), htitl.c_str(), wadc, adc1, adc2);
  phf->SetDirectory(nullptr);
  for ( Index isam=0; isam<nsam; ++isam ) {
    phf->Fill(acd.raw[isam]);
  }
  phf->SetStats(0);
  phf->SetLineWidth(2);
  // Fetch the peak bin and suppress it for the fit if more than 20% (but not
  // all) the data is in it.
  int binmax = phf->GetMaximumBin();
  double valmax = phf->GetBinContent(binmax);
  double rangeIntegral = phf->Integral(1, phf->GetNbinsX());
  double peakBinFraction = valmax/rangeIntegral;
  bool allBin = peakBinFraction > 0.99;
  bool dropBin = valmax > 0.2*phf->Integral() && !allBin;
  int nbinsRemoved = 0;
  if ( dropBin ) {
    double tmpval = 0.5*(phf->GetBinContent(binmax-1)+phf->GetBinContent(binmax+1));
    phf->SetBinContent(binmax, tmpval);
    rangeIntegral = phf->Integral(1, phf->GetNbinsX());
    ++nbinsRemoved;
  }
  double amean = phf->GetMean() + 0.5;
  double alim1 = amean - 25.0;
  double alim2 = amean + 25.0;
  TF1 fitter("pedgaus", "gaus", adc1, adc2, TF1::EAddToList::kNo);
  fitter.SetParameters(0.1*rangeIntegral, amean, 5.0);
  fitter.SetParLimits(0, 0.01*rangeIntegral, rangeIntegral);
  fitter.SetParLimits(1, alim1, alim2);
  if ( allBin ) fitter.FixParameter(1, amean);  // Fix posn.
  fitter.SetParLimits(2, 3.0, 10.0);
  TF1* pfinit = dynamic_cast<TF1*>(fitter.Clone("pedgaus0"));
  pfinit->SetLineColor(3);
  pfinit->SetLineStyle(2);
  string fopt = "0";
  if ( m_LogLevel < 2 ) fopt += "WWQB";
  phf->Fit(&fitter, fopt.c_str());
  phf->GetListOfFunctions()->AddLast(pfinit, "0");
  phf->GetListOfFunctions()->Last()->SetBit(TF1::kNotDraw, true);
  if ( dropBin ) phf->SetBinContent(binmax, valmax);
  res.setHist("pedestal", phf);
  res.setFloat("fitPedestal", fitter.GetParameter(1) - 0.5);
  res.setFloat("fitPedestalRms", fitter.GetParameter(2));
  res.setFloat("fitChiSquare", fitter.GetChisquare());
  res.setFloat("fitPeakBinFraction", peakBinFraction);
  res.setInt("fitChannel", acd.channel);
  res.setInt("fitNBinsRemoved", nbinsRemoved);
  if ( m_LogLevel >= 3 ) cout << myname << "Exiting..." << endl;
  return res;
}

//**********************************************************************
