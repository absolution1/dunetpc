// ToolBasedRawDigitPrepService_service.cc

#include "ToolBasedRawDigitPrepService.h"
#include "dune/ArtSupport/DuneToolManager.h"
#include "dune/DuneInterface/Tool/AdcChannelTool.h"
#include "dune/DuneInterface/AdcWireBuildingService.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include <iostream>

using std::string;
using std::cout;
using std::endl;
using std::vector;
using std::ostringstream;
using raw::RawDigit;

//**********************************************************************

ToolBasedRawDigitPrepService::
ToolBasedRawDigitPrepService(fhicl::ParameterSet const& pset, art::ActivityRegistry&)
: m_LogLevel(pset.get<int>("LogLevel")),
  m_DoWires(pset.get<bool>("DoWires")),
  m_AdcChannelToolNames(pset.get<vector<string>>("AdcChannelToolNames")),
  m_pWireBuildingService(nullptr) {
  const string myname = "ToolBasedRawDigitPrepService::ctor: ";
  pset.get_if_present<int>("LogLevel", m_LogLevel);
  // Fetch the tools.
  if ( m_AdcChannelToolNames.size() ) {
    DuneToolManager* ptm = DuneToolManager::instance("");
    if ( ptm == nullptr ) {
      cout << myname << "ERROR: Unable to retrieve tool manaager." << endl;
    } else {
      for ( string tname : m_AdcChannelToolNames ) {
        if ( m_LogLevel ) cout << myname << "     Fetching " << tname << endl;
        AdcChannelToolPtr ptool = ptm->getPrivate<AdcChannelTool>(tname);
        NamedTool nt(tname, ptool.get());
        if ( nt.tool ) {
          if ( m_LogLevel ) cout << myname << "    Found tool " << tname << " @ " << nt.tool << endl;
          m_AdcChannelTools.push_back(std::move(ptool));
          m_AdcChannelNamedTools.push_back(nt);
        } else {
          cout << myname << "ERROR: Unable to retrieve display tool " << tname << endl;
        }
      }
    }
  }
  if ( m_DoWires ) {
    if ( m_LogLevel ) cout << myname << "Fetching wire building service." << endl;
    m_pWireBuildingService = &*art::ServiceHandle<AdcWireBuildingService>();
    if ( m_LogLevel ) cout << myname << "  Wire building service: @" <<  m_pWireBuildingService << endl;
  }
  if ( m_LogLevel >=1 ) print(cout, myname);
}

//**********************************************************************

int ToolBasedRawDigitPrepService::
prepare(AdcChannelDataMap& datamap,
        std::vector<recob::Wire>* pwires, WiredAdcChannelDataMap* pintStates) const {
  const string myname = "ToolBasedRawDigitPrepService:prepare: ";
  // Loop over tools.
  if ( m_LogLevel >= 2 ) cout << myname << "Processing " << datamap.size() << " channels with "
                              << m_AdcChannelNamedTools.size() << " tools." << endl;
  if ( m_AdcChannelNamedTools.size() ) {
    for ( NamedTool nt : m_AdcChannelNamedTools ) {
      if ( m_LogLevel >= 3 ) cout << myname << "  Running tool " << nt.name << endl;
      DataMap ret = nt.tool->updateMap(datamap);
      if ( ret ) cout << myname << "WARNING: Tool " << nt.name << " failed with error " << ret.status() << endl;
      if ( m_LogLevel >= 4 ) {
        cout << myname << "----------------------------" << endl;
        ret.print();
        cout << myname << "----------------------------" << endl;
      }
    }
  }
  if ( m_DoWires ) {
    for ( auto& chdata : datamap ) {
      AdcChannelData& acd = chdata.second;
      m_pWireBuildingService->build(acd, pwires);
    }
  }
  return 0;
}

//**********************************************************************

std::ostream& ToolBasedRawDigitPrepService::
print(std::ostream& out, std::string prefix) const {
  out << prefix << "ToolBasedRawDigitPrepService:"                      << endl;
  out << prefix << "                    LogLevel: " << m_LogLevel             << endl;
  out << prefix << "                     DoWires: " << m_DoWires              << endl;
  if ( m_AdcChannelNamedTools.size() ) {
    cout << prefix << "     ADC channel tools:";
    for ( const NamedTool& nm : m_AdcChannelNamedTools ) {
       out << "\n" << prefix << "           " << nm.name;
    }
    cout << endl;
  } else {
    out << prefix << "    No ADC channel tools." << endl;
  }
  return out;
}

//**********************************************************************

DEFINE_ART_SERVICE_INTERFACE_IMPL(ToolBasedRawDigitPrepService, RawDigitPrepService)

//**********************************************************************

