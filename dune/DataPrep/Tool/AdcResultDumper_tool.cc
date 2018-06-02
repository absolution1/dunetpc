// AdcResultDumper_tool.cc

#include "AdcResultDumper.h"
#include "dune/ArtSupport/DuneToolManager.h"
#include <iostream>

using std::string;
using std::cout;
using std::endl;

using AdcChannelToolPtr = std::unique_ptr<AdcChannelTool>;

//**********************************************************************
// Class methods.
//**********************************************************************

AdcResultDumper::AdcResultDumper(fhicl::ParameterSet const& ps)
: m_LogLevel(ps.get<int>("LogLevel")), 
  m_Tool(ps.get<string>("Tool"))
{
  const string myname = "AdcResultDumper::ctor: ";
  if ( m_LogLevel ) {
    cout << myname << "Configuration: " << endl;
    cout << myname << "   LogLevel: " << m_LogLevel << endl;
    cout << myname << "       Tool: " << m_Tool << endl;
  }
  if ( m_LogLevel ) cout << myname << "Retrieving tool." << endl;
  DuneToolManager* ptm = DuneToolManager::instance();
  m_ptool = ptm->getPrivate<AdcChannelTool>(m_Tool);
  if ( m_ptool == nullptr ) {
    cout << myname << "Tool retrieval failed." << endl;
  } else if ( m_LogLevel ) {
    cout << myname << "Tool retrieval succeeded." << endl;
  }
}

//**********************************************************************

DataMap AdcResultDumper::view(const AdcChannelData& acd) const {
  const string myname = "AdcResultDumper::view: ";
  if ( m_ptool == nullptr ) return DataMap(101);
  if ( m_LogLevel >= 2 ) cout << myname << "Calling tool " << m_Tool << endl;
  DataMap ret = m_ptool->view(acd);
  ret.print();
  return ret;
}

//**********************************************************************

DataMap AdcResultDumper::update(AdcChannelData& acd) const {
  const string myname = "AdcResultDumper::updateMap: ";
  if ( m_ptool == nullptr ) return DataMap(101);
  if ( m_LogLevel >= 2 ) cout << myname << "Calling tool " << m_Tool << endl;
  DataMap ret = m_ptool->update(acd);
  ret.print();
  return ret;
}

//**********************************************************************

DataMap AdcResultDumper::viewMap(const AdcChannelDataMap& acds) const {
  const string myname = "AdcResultDumper::viewMap: ";
  if ( m_ptool == nullptr ) return DataMap(101);
  if ( m_LogLevel >= 2 ) cout << myname << "Calling tool " << m_Tool << endl;
  DataMap ret = m_ptool->viewMap(acds);
  ret.print();
  return ret;
}

//**********************************************************************

DataMap AdcResultDumper::updateMap(AdcChannelDataMap& acds) const {
  const string myname = "AdcResultDumper::updateMap: ";
  if ( m_ptool == nullptr ) return DataMap(101);
  if ( m_LogLevel >= 2 ) cout << myname << "Calling tool " << m_Tool << endl;
  DataMap ret = m_ptool->updateMap(acds);
  ret.print();
  return ret;
}

//**********************************************************************
