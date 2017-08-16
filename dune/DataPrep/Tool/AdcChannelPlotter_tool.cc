// AdcChannelPlotter_tool.cc

#include "AdcChannelPlotter.h"
#include "dune/DuneCommon/StringManipulator.h"
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

AdcChannelPlotter::AdcChannelPlotter(fhicl::ParameterSet const& ps)
: m_LogLevel(ps.get<int>("LogLevel")),
  m_HistTypes(ps.get<NameVector>("HistTypes")),
  m_HistName(ps.get<string>("HistName")),
  m_HistTitle(ps.get<string>("HistTitle")),
  m_RootFileName(ps.get<string>("RootFileName")) {
  const string myname = "AdcChannelPlotter::ctor: ";
  if ( m_HistTypes.size() == 0 ) {
    cout << myname << "WARNING: No histogram types are specified." << endl;
    return;
  }
  if ( m_LogLevel > 0 ) {
    cout << myname << "      LogLevel: [" << m_LogLevel << endl;
    cout << myname << "     HistTypes: [";
    bool first = true;
    for ( string name : m_HistTypes ) {
      if ( ! first ) cout << ", ";
      first = false;
      cout << " " << name;
    }
    cout << "]" << endl;
    cout << myname << "      HistName: " << m_HistName << endl;
    cout << myname << "     HistTitle: " << m_HistTitle << endl;
    cout << myname << "  RootFileName: " << m_RootFileName << endl;
  }
}

//**********************************************************************

int AdcChannelPlotter::view(const AdcChannelData& acd) const {
  const string myname = "AdcChannelPlotter::view: ";
  if ( m_HistTypes.size() == 0 ) {
    cout << myname << "WARNING: No histogram types are specified." << endl;
    return 1;
  }
  string hnameBase = m_HistName;
  if ( hnameBase == "" ) hnameBase = "%TYPE%";
  string htitlBase = m_HistTitle;
  if ( htitlBase == "" ) htitlBase = "%TYPE%";
  TDirectory* polddir = gDirectory;
  TFile* pfile = nullptr;
  if ( m_RootFileName.size() ) {
    string fname = nameReplace(m_RootFileName, acd, m_HistTypes[0]);
    pfile = TFile::Open(fname.c_str(), "UPDATE");
  }
  vector<TH1*> hists;
  for ( string type : m_HistTypes ) {
    if ( type == "raw" ) {
      Index nsam = acd.raw.size();
      if ( nsam == 0 ) {
        cout << myname << "WARNING: Raw data is empty." << endl;
        continue;
      }
      string hname = nameReplace(hnameBase, acd, type);
      string htitl = nameReplace(htitlBase, acd, type);
      htitl += "; Tick; ADC count";
      TH1* ph = new TH1F(hname.c_str(), htitl.c_str(), nsam, 0, nsam);
      hists.push_back(ph);
      for ( Index isam=0; isam<nsam; ++isam ) {
        ph->SetBinContent(isam+1, acd.raw[isam]);
      }
    } else if ( type == "rawdist" ) {
      Index nsam = acd.raw.size();
      if ( nsam == 0 ) {
        cout << myname << "WARNING: Raw data is empty." << endl;
        continue;
      }
      string hname = nameReplace(hnameBase, acd, type);
      string htitl = nameReplace(htitlBase, acd, type);
      htitl += "; ADC count; # samples";
      unsigned int nadc = 4096;
      TH1* ph = new TH1F(hname.c_str(), htitl.c_str(), nadc, 0, nadc);
      hists.push_back(ph);
      for ( Index isam=0; isam<nsam; ++isam ) {
        ph->Fill(acd.raw[isam]);
      }
    } else {
      cout << myname << "WARNING: Unknown type: " << type << endl;
    }
  }
  for ( TH1* ph : hists ) {
    ph->SetStats(0);
    ph->SetLineWidth(2);
  }
  if ( pfile != nullptr ) {
    pfile->Write();
    pfile->Close();
    delete pfile;
    gDirectory = polddir;
  }
  return 0;
}

//**********************************************************************

string AdcChannelPlotter::
nameReplace(string name, const AdcChannelData& acd, string type) const {
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
  if ( type.size() ) sman.replace("%TYPE%", type);
  return nameout;
}

//**********************************************************************
