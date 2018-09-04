// AdcEventViewer_tool.cc

#include "AdcEventViewer.h"
#include "dune/DuneInterface/Tool/AdcChannelStringTool.h"
#include "dune/ArtSupport/DuneToolManager.h"
#include "dune/DuneInterface/Tool/RunDataTool.h"
#include "dune/DuneCommon/TPadManipulator.h"
#include "TH1.h"
#include <iostream>
#include <sstream>
#include <iomanip>

using std::string;
using std::to_string;
using std::cout;
using std::endl;
using std::istringstream;
using fhicl::ParameterSet;
using std::setw;
using std::fixed;
using std::setprecision;

using Index = AdcEventViewer::Index;

//**********************************************************************
// Class methods.
//**********************************************************************

AdcEventViewer::AdcEventViewer(fhicl::ParameterSet const& ps)
: m_LogLevel(ps.get<int>("LogLevel")),
  m_EventHists(ps.get<NameVector>("EventHists")),
  m_state(new AdcEventViewer::State) {
  const string myname = "AdcEventViewer::ctor: ";
  // Display the configuration.
  if ( m_LogLevel>= 1 ) {
    cout << myname << "         LogLevel: " << m_LogLevel << endl;
    cout << myname << "       EventHists: [" << endl;
    for ( Index ievh=0; ievh<m_EventHists.size(); ++ievh ) {
      cout << myname << "                     " << m_EventHists[ievh] << endl;
    }
    cout << myname << "                   " << "]" << endl;
  }
  state().run = 0;
  state().event = 0;
  const string::size_type& npos = string::npos;
  for ( string hspec : m_EventHists ) {
    string name;
    int nbin = 0;
    float xmin = 0.0;
    float xmax = 0.0;
    bool ok = false;
    string::size_type ipos = 0;
    string::size_type jpos = 0;
    istringstream ssin;
    jpos = hspec.find(":");
    if ( jpos != npos ) {
      name = hspec.substr(ipos, jpos-ipos);
      ipos = jpos + 1;
      jpos = hspec.find(":", ipos);
      if ( jpos != npos && jpos > ipos ) {
        istringstream ssnbin(hspec.substr(ipos, jpos-ipos));
        ssnbin >> nbin;
        ipos = jpos + 1;
        jpos = hspec.find(":", ipos);
        if ( jpos != npos && jpos > ipos ) {
          istringstream ssxmin(hspec.substr(ipos, jpos-ipos));
          ssxmin >> xmin;
          ipos = jpos + 1;
          jpos = hspec.size();
          istringstream ssxmax(hspec.substr(ipos, jpos-ipos));
          ssxmax >> xmax;
          ok = true;
        }
      }
    }
    if ( ! ok ) {
      cout << "WARNING: Invalid histogram configuration string: " << hspec << endl;
      continue;
    }
    string sttl;
    if ( name.find("nfemb") != string::npos ) {
      sttl = "FEMB counts;# FEMB;# event";
    } else {
      cout << myname << "ERROR: No variable for histogram name " << name << endl;
      continue;
    }
    if ( m_LogLevel >= 2 ) {
      cout << myname << "Creating in histogram " << name << ", nbin=" << nbin
           << ", range=(" << xmin << ", " << xmax << ")" << endl;
    }
    TH1* ph = new TH1F(name.c_str(), sttl.c_str(), nbin, xmin, xmax);
    ph->SetDirectory(nullptr);
    ph->SetStats(0);
    ph->SetLineWidth(2);
    state().hists.push_back(ph);
    if ( m_LogLevel>= 1 ) cout << myname << "Created histogram " << name << endl;
  }
}

//**********************************************************************

AdcEventViewer::~AdcEventViewer() {
  const string myname = "AdcEventViewer::dtor: ";
  endEvent();
  displayHists();
  cout << myname << "Exiting." << endl;
}

//**********************************************************************

DataMap AdcEventViewer::view(const AdcChannelData& acd) const {
  DataMap res;
  Index ievt = acd.event;
  Index irun = acd.run;
  if ( ievt != state().event || irun != state().run ) {
    startEvent(acd);
  }
  state().fembIDSet.insert(acd.fembID);
  ++state().nchan;
  return res;
}

//**********************************************************************

DataMap AdcEventViewer::viewMap(const AdcChannelDataMap& acds) const {
  const string myname = "AdcEventViewer::viewMap: ";
  DataMap ret;
  if ( acds.size() == 0 ) {
    if ( m_LogLevel >=2 ) cout << myname << "Skipping group with no data" << endl;
    return ret;
  }
  const AdcChannelData& acd = acds.begin()->second;
  if ( acd.event != state().event ) startEvent(acd);
  ++state().ngroup;
  for ( const AdcChannelDataMap::value_type& iacd : acds ) view(iacd.second);
  return ret;
}

//**********************************************************************

void AdcEventViewer::startEvent(const AdcChannelData& acd) const {
  const string myname = "AdcEventViewer::startEvent: ";
  endEvent();
  if ( acd.run != state().run ) {
    if ( state().run != 0 ) {
      if ( m_LogLevel >= 2 ) {
        cout << myname << "Encountered new run " << acd.run
             << ". Ending old run " << state().run << endl;
      }
      displayHists();
      state().events.clear();
      state().eventSet.clear();
    }
    state().run = acd.run;
  }
  Index ievt = acd.event;
  state().event = ievt;
  state().events.push_back(ievt);
  state().eventSet.insert(ievt);
  state().ngroup = 0;
  state().fembIDSet.clear();
  state().nchan = 0;
}

//**********************************************************************

void AdcEventViewer::endEvent() const {
  const string myname = "AdcEventViewer::endEvent: ";
  printReport();
}

//**********************************************************************

void AdcEventViewer::printReport() const {
  const string myname = "AdcEventViewer::printReport: ";
  if ( state().event == 0 ) return;
  Index nevt = state().events.size();
  Index ndup = nevt - state().eventSet.size();
  Index nfmb = state().fembIDSet.size();
  Index nchn = state().nchan;
  double chanPerFemb = nfmb > 0 ? double(nchn)/nfmb : 0.0;
  if ( m_LogLevel >= 1 ) {
    const int w = 5;
    cout << myname << "               event: " << setw(w) << state().event << endl;
    cout << myname << "            # events: " << setw(w) << state().events.size() << endl;
    cout << myname << "  # duplicate events: " << setw(w) << ndup << endl;
    cout << myname << "            # groups: " << setw(w) << state().ngroup << endl;
    cout << myname << "             # FEMBs: " << setw(w) << nfmb << endl;
    cout << myname << "          # channels: " << setw(w) << nchn << endl;
    cout << myname << "     # channels/FEMB: " << setw(w+5) << fixed << setprecision(4) << chanPerFemb << endl;
  }
  for ( TH1* ph : state().hists ) {
    string name = ph->GetName();
    if ( name.find("nfemb") != string::npos ) ph->Fill(nfmb);
    else {
      cout << myname << "ERROR: No variable for histogram name " << name << endl;
    }
  }
}

//**********************************************************************

void AdcEventViewer::displayHists() const {
  const string myname = "AdcEventViewer::displayHists: ";
  string sttlSuf = " for run " + to_string(state().run);
  Index nplt = state().hists.size();
  if ( m_LogLevel >= 1 ) cout << myname << "Creating " << nplt << " plot"
                              << (nplt == 1 ? "" : "s") << sttlSuf
                              << (nplt > 0 ? ":" : "") << endl;
  Index nevt = state().eventSet.size();
  if ( nevt == 0 ) sttlSuf += " with no events.";
  else if ( nevt == 1 ) sttlSuf += " event " + to_string(*state().eventSet.begin());
  else sttlSuf += " events " + to_string(*state().eventSet.begin()) + "-" + to_string(*state().eventSet.rbegin());
  for ( TH1* ph : state().hists ) {
    TPadManipulator man;
    man.add(ph);
    string sttl = ph->GetTitle() + sttlSuf;
    man.setTitle(sttl.c_str());
    string fname = string("eview_") + ph->GetName() + "_run" + to_string(state().run) + ".png";
    man.showUnderflow();
    man.showOverflow();
    man.addAxis();
    man.print(fname);
    if ( m_LogLevel >= 1 ) cout << myname << "  " << fname << endl;
  }
}

//**********************************************************************
