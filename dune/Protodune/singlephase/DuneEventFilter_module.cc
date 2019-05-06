// DuneEventFilter_module.cc
//
// Module to select or reject events based on event number.
//
// David Adams
// October 2018
//
// Configuration:
//   LogLevel: 0=quiet, 1=ctor message, 2=status of each event
//   SelectEvents: If this has entries, the event must be among them
//   RejectEvents: Event is rejected if it is in this list.
//   EventBegin, EventEnd: if EventEnd > EventBegin, event must be in (EventBegin, EventEnd].
//   EventModFreq, EventModVal: keep events for which ievt%EventModFreq == EventModVal

#include <iostream>
#include <vector>
#include <set>

#include "art/Framework/Core/EDFilter.h" 
#include "art/Framework/Core/ModuleMacros.h" 
#include "art/Framework/Principal/Event.h" 

class DuneEventFilter : public art::EDFilter {

public:

  using Index = unsigned int;
  using IndexVector = std::vector<Index>;
  using IndexSet = std::set<Index>;

  explicit DuneEventFilter(fhicl::ParameterSet const & pset);
  virtual ~DuneEventFilter();
  virtual bool filter(art::Event& e);

private:

  // Configuration parameters.
  Index m_LogLevel;
  IndexVector m_SelectEventVector;
  IndexVector m_RejectEventVector;
  Index m_EventBegin;
  Index m_EventEnd;
  Index m_EventModFreq;
  Index m_EventModVal;

  // Derived from configuration.
  IndexSet m_SelectEvents;
  IndexSet m_RejectEvents;

  // Counters.
  Index m_nproc;
  Index m_nsel;

};

//*****************************************************************************

DuneEventFilter::DuneEventFilter(fhicl::ParameterSet const & pset)
: EDFilter(pset),
  m_LogLevel(pset.get<Index>("LogLevel")),
  m_SelectEventVector(pset.get<IndexVector>("SelectEvents")),
  m_RejectEventVector(pset.get<IndexVector>("RejectEvents")),
  m_EventBegin(pset.get<Index>("EventBegin")),
  m_EventEnd(pset.get<Index>("EventEnd")),
  m_EventModFreq(pset.get<Index>("EventModFreq")),
  m_EventModVal(pset.get<Index>("EventModVal")),
  m_nproc(0), m_nsel(0) {
  using std::cout;
  using std::endl;
  using std::string;
  const string myname = "DuneEventFilter::ctor: ";
  for ( Index ievt : m_SelectEventVector ) m_SelectEvents.insert(ievt);
  for ( Index ievt : m_RejectEventVector ) m_RejectEvents.insert(ievt);
  if ( m_LogLevel >= 1 ) {
    cout << myname << "      LogLevel: " << m_LogLevel << endl;
    cout << myname << "  SelectEvents: [";
    bool first = true;
    for ( unsigned int ievt : m_SelectEvents ) {
      if ( first ) first = false;
      else cout << ", ";
      cout << ievt;
    }
    cout << "]" << endl;
    cout << myname << "  RejectEvents: [";
    first = true;
    for ( unsigned int ievt : m_RejectEvents ) {
      if ( first ) first = false;
      else cout << ", ";
      cout << ievt;
    }
    cout << "]" << endl;
    cout << myname << "    EventBegin: " << m_EventBegin << endl;
    cout << myname << "      EventEnd: " << m_EventEnd << endl;
    cout << myname << "   EventModVal: " << m_EventModVal << endl;
    cout << myname << "  EventModFreq: " << m_EventModFreq << endl;
  }
}

//*****************************************************************************

bool DuneEventFilter::filter(art::Event & evt) {
  using std::cout;
  using std::endl;
  using std::string;
  const string myname = "DuneEventFilter::filter: ";
  Index ievt = evt.event();
  ++m_nproc;
  bool keep = true;
  if ( keep && m_SelectEvents.size() ) keep = m_SelectEvents.count(ievt);
  if ( keep ) keep = m_RejectEvents.count(ievt) == 0;
  if ( keep && m_EventEnd > m_EventBegin ) keep = ievt >= m_EventBegin && ievt < m_EventEnd;
  if ( keep && m_EventModFreq ) keep = (ievt % m_EventModFreq) == m_EventModVal;
  if ( m_LogLevel >= 2) {
    cout << myname << (keep ? "Sel" : "Rej") << "ecting event " << ievt << endl;
  }
  if ( keep ) ++m_nsel;
  return keep;
}

//*****************************************************************************

DuneEventFilter::~DuneEventFilter() {
  using std::cout;
  using std::endl;
  using std::string;
  const string myname = "DuneEventFilter::dtor: ";
  if ( m_LogLevel >= 1 ) {
    float fsel = double(m_nsel)/double(m_nproc);
    cout << myname << " Events processed: " << m_nproc << endl;
    cout << myname << "  Events selected: " << m_nsel << " (" << fsel << ")" << endl;
  }
}

//*****************************************************************************

DEFINE_ART_MODULE(DuneEventFilter)
