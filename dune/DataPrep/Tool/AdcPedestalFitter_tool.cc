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
  m_HistManager(ps.get<string>("HistManager")) {
  const string myname = "AdcPedestalFitter::ctor: ";
  if ( m_HistManager.size() ) {
    DuneToolManager* ptm = DuneToolManager::instance();
    m_phm = ptm->getShared<HistogramManager>(m_HistManager);
    if ( m_phm == nullptr ) {
      cout << myname << "WARNING: Histogram manager not found: " << m_HistManager << endl;
    }
  }
  if ( m_LogLevel > 0 ) {
    cout << myname << "      LogLevel: [" << m_LogLevel << endl;
    cout << myname << "      HistName: " << m_HistName << endl;
    cout << myname << "     HistTitle: " << m_HistTitle << endl;
    cout << myname << "   HistManager: " << m_HistManager << endl;
  }
}

//**********************************************************************

int AdcPedestalFitter::view(const AdcChannelData& acd) const {
  const string myname = "AdcPedestalFitter::view: ";
  string hnameBase = m_HistName;
  string htitlBase = m_HistTitle;
  Index nsam = acd.raw.size();
  if ( nsam == 0 ) {
    cout << myname << "WARNING: Raw data is empty." << endl;
    return 1;
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
  int binmax = phr->GetMaximumBin();
  double adcmax = phr->GetBinCenter(binmax);
  double adcmean = phr->GetMean();
  //double adcrms = phr->GetRMS();
  delete phr;
  double wadc = 100.0;
  double adc1 = adcmax - 0.5*wadc;
  adc1 = 10*int(adc1/10);
  if ( adcmean > adcmax + 10) adc1 += 10;
  double adc2 = adc1 + wadc;
  TH1* phf = new TH1F(hname.c_str(), htitl.c_str(), wadc, adc1, adc2);
  for ( Index isam=0; isam<nsam; ++isam ) {
    phf->Fill(acd.raw[isam]);
  }
  phf->SetStats(0);
  phf->SetLineWidth(2);
  phf->Fit("gaus");
  int rstat = m_phm != nullptr;
  if ( rstat != 0 ) {
    if ( m_phm->manage(phf) ) {
      cout << myname << "WARNING: Attempt to manage histogram " << phf->GetName()
           << " returned error " << rstat << endl;
      phf->SetDirectory(nullptr);
      delete phf;
    }
  }
  return 0;
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
