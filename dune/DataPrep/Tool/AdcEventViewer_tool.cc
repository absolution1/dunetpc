// AdcEventViewer_tool.cc

#include "AdcEventViewer.h"
#include "dune/DuneInterface/Tool/AdcChannelStringTool.h"
#include "dune/ArtSupport/DuneToolManager.h"
#include "dune/DuneInterface/Tool/RunDataTool.h"
#include <iostream>
#include <iomanip>

using std::string;
using std::cout;
using std::endl;
using fhicl::ParameterSet;
using std::setw;

using Index = AdcEventViewer::Index;

//**********************************************************************
// Class methods.
//**********************************************************************

AdcEventViewer::AdcEventViewer(fhicl::ParameterSet const& ps)
: m_LogLevel(ps.get<int>("LogLevel")),
  m_state(new AdcEventViewer::State) {
  const string myname = "AdcEventViewer::ctor: ";
  // Display the configuration.
  if ( m_LogLevel>= 1 ) {
    cout << myname << "         LogLevel: " << m_LogLevel << endl;
  }
  state().event = 0;
}

//**********************************************************************

AdcEventViewer::~AdcEventViewer() {
  const string myname = "AdcEventViewer::dtor: ";
  printReport();
  cout << myname << "Exiting." << endl;
}

//**********************************************************************

DataMap AdcEventViewer::view(const AdcChannelData& acd) const {
  DataMap res;
  Index ievt = acd.event;
  if ( ievt != state().event ) {
    initializeState(ievt);
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
  Index ievt = acds.begin()->second.event;
  if ( ievt != state().event ) initializeState(ievt);
  ++state().ngroup;
  for ( const AdcChannelDataMap::value_type& iacd : acds ) view(iacd.second);
  return ret;
}

//**********************************************************************

void AdcEventViewer::initializeState(Index ievt) const {
  printReport();
  state().event = ievt;
  state().events.push_back(ievt);
  state().eventSet.insert(ievt);
  state().ngroup = 0;
  state().fembIDSet.clear();
  state().nchan = 0;
}

//**********************************************************************

void AdcEventViewer::printReport() const {
  const string myname = "AdcEventViewer::printReport: ";
  if ( state().event == 0 ) return;
  if ( m_LogLevel >= 1 ) {
    const int w = 6;
    Index nevt = state().events.size();
    Index ndup = nevt - state().eventSet.size();
    cout << myname << "               event: " << setw(w) << state().event << endl;
    cout << myname << "            # events: " << setw(w) << state().events.size() << endl;
    cout << myname << "  # duplicate events: " << setw(w) << ndup << endl;
    cout << myname << "            # groups: " << setw(w) << state().ngroup << endl;
    cout << myname << "             # FEMBs: " << setw(w) << state().fembIDSet.size() << endl;
    cout << myname << "          # channels: " << setw(w) << state().nchan << endl;
  }
}

//**********************************************************************
