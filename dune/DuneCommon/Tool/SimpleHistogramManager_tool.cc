// SimpleHistogramManager.cxx

#include "SimpleHistogramManager.h"
#include "TH1.h"
#include <iostream>

using std::string;
using std::cout;
using std::endl;

using NameVector = HistogramManager::NameVector;

//**********************************************************************

SimpleHistogramManager::SimpleHistogramManager(fhicl::ParameterSet const& ps)
:  m_LogLevel(ps.get<int>("LogLevel")) {
  const string myname = "SimpleHistogramManager::ctor: ";
  if ( m_LogLevel > 0 ) {
    cout << myname << "LogLevel: " << m_LogLevel << endl;
  }
}

//**********************************************************************

SimpleHistogramManager::~SimpleHistogramManager() {
  release("*");
}

//**********************************************************************

int SimpleHistogramManager::manage(TH1* ph) {
  const string myname = "SimpleHistogramManager::manage: ";
  if ( ph == nullptr ) {
    if ( m_LogLevel > 1 ) {
      cout << myname << "Histogram pointer is null." << endl;
    }
    return 1;
  }
  string hname = ph->GetName();
  //if ( ph->GetDirectory() == nullptr ) {
  //  if ( m_LogLevel > 1 ) {
  //    cout << myname << "Histogram " << hname << " is not in a directory." << endl;
  //  }
  //  return 2;
  //}
  if ( m_hists.find(hname) != m_hists.end() ) {
    if ( m_LogLevel > 1 ) {
      cout << myname << "Histogram  with name " << hname << " is already managed here." << endl;
    }
    return 3;
  }
  ph->SetDirectory(nullptr);
  m_hists[hname] = ph;
  m_names.push_back(hname);
  return 0;
}
  
//**********************************************************************

NameVector SimpleHistogramManager::names() const {
  return m_names;
}

//**********************************************************************

TH1* SimpleHistogramManager::get(Name hname) const {
  std::map<Name, TH1*>::const_iterator ihst = m_hists.find(hname);
  if ( ihst == m_hists.end() ) return nullptr;
  return ihst->second;
}
  
//**********************************************************************

TH1* SimpleHistogramManager::last() const {
  if ( m_hists.size() == 0 ) return nullptr;
  return get(names().back());
}

//**********************************************************************

int SimpleHistogramManager::release(Name hname) {
  if ( hname == "*" ) {
    for ( auto ent : m_hists ) delete ent.second;
    m_names.clear();
    m_hists.clear();
  }
  std::map<Name, TH1*>::const_iterator ihst = m_hists.find(hname);
  if ( ihst == m_hists.end() ) return 1;
  NameVector::iterator inam = find(m_names.begin(), m_names.end(), hname );
  if ( inam == m_names.end() ) return 2;
  m_names.erase(inam);
  m_hists.erase(ihst);
  return 0;
}

//**********************************************************************
