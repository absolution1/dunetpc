// ToolBasedRawDigitPrepService_service.cc

#include "ToolBasedRawDigitPrepService.h"
#include "dune/ArtSupport/DuneToolManager.h"
#include "dune/DuneInterface/Tool/AdcChannelTool.h"
#include "dune/DuneInterface/Service/AdcWireBuildingService.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Registry/ServiceDefinitionMacros.h"
#include <iostream>
#include <iomanip>

#include "valgrind/callgrind.h"

using std::string;
using std::cout;
using std::endl;
using std::vector;
using std::ostringstream;
using std::setw;
using std::setprecision;

using raw::RawDigit;

using Index = unsigned int;

//**********************************************************************

ToolBasedRawDigitPrepService::
ToolBasedRawDigitPrepService(fhicl::ParameterSet const& pset, art::ActivityRegistry&)
: m_LogLevel(pset.get<int>("LogLevel")),
  m_DoWires(pset.get<bool>("DoWires")),
  m_ToolNames(pset.get<vector<string>>("ToolNames")),
  m_CallgrindToolNames(pset.get<vector<string>>("CallgrindToolNames")),
  m_pWireBuildingService(nullptr),
  m_cgset(m_CallgrindToolNames.begin(), m_CallgrindToolNames.end()),
  m_pstate(new State) {
  const string myname = "ToolBasedRawDigitPrepService::ctor: ";
  pset.get_if_present<int>("LogLevel", m_LogLevel);
  // Fetch the tools.
  if ( m_ToolNames.size() ) {
    DuneToolManager* ptm = DuneToolManager::instance("");
    if ( ptm == nullptr ) {
      cout << myname << "ERROR: Unable to retrieve tool manaager." << endl;
    } else {
      for ( string tname : m_ToolNames ) {
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
  state().toolTimes.resize(m_ToolNames.size());
  if ( m_DoWires ) {
    if ( m_LogLevel ) cout << myname << "Fetching wire building service." << endl;
    m_pWireBuildingService = &*art::ServiceHandle<AdcWireBuildingService>();
    if ( m_LogLevel ) cout << myname << "  Wire building service: @" <<  m_pWireBuildingService << endl;
  }
  if ( m_LogLevel >=1 ) print(cout, myname);
}

//**********************************************************************

ToolBasedRawDigitPrepService::~ToolBasedRawDigitPrepService() {
  const string myname = "ToolBasedRawDigitPrepService:dtor: ";
  if ( state().nevtBegin != state().nevtEnd ) {
    cout << myname << "WARNING: Event counts are inconsistent: " << state().nevtBegin
         << " != " << state().nevtEnd << endl;
  }
  Index ntoo = m_ToolNames.size();
  Index nevt = state().nevtEnd;
  if ( m_LogLevel >= 1 ) {
    cout << myname << "Event count: " << nevt << endl;
    cout << myname << " Call count: " << state().ncall << endl;
    cout << myname << "Time report for " << ntoo << " tools." << endl;
    string sunit = "sec/event";
    float xnevt = float(nevt);
    if ( nevt == 0 ) {
      xnevt = 1.0;
      sunit = "sec";
    }
    for ( Index itoo=0; itoo<ntoo; ++itoo ) {
      string name = m_ToolNames[itoo];
      double time = state().toolTimes[itoo].count();
      cout << myname << setw(30) << name << ":"
           << setw(7) << std::fixed << setprecision(2)
           << time/xnevt << " " << sunit << endl;
    }
  }
}

//**********************************************************************

int ToolBasedRawDigitPrepService::beginEvent(const DuneEventInfo& devt) const {
  const string myname = "ToolBasedRawDigitPrepService:beginEvent: ";
  if ( m_LogLevel >= 2 ) {
    cout << myname << "Begin processing run " << devt.runString();
    cout << " event " << devt.event;
    cout << " with " << m_AdcChannelNamedTools.size() << " tools." << endl;
  }
  if ( state().nevtBegin != state().nevtEnd ) {
    cout << myname << "WARNING: Event counts are inconsistent: " << state().nevtBegin
         << " != " << state().nevtEnd << endl;
  }
  ++state().nevtBegin;
  Index nfail = 0;
  if ( m_AdcChannelNamedTools.size() ) {
    for ( NamedTool nt : m_AdcChannelNamedTools ) {
      DataMap ret = nt.tool->beginEvent(devt);
      if ( ret.status() ) {
        ++nfail;
        cout << myname << "WARNING: Iniitalization for tool " << nt.name
             << " failed for event " << devt.event << " with status code " << ret.status() << endl;
      }
    }
  }
  return nfail;
}

//**********************************************************************

int ToolBasedRawDigitPrepService::endEvent(const DuneEventInfo& devt) const {
  const string myname = "ToolBasedRawDigitPrepService:endEvent: ";
  if ( m_LogLevel >= 2 ) {
    cout << myname << "End processing run " << devt.runString();
    cout << " event " << devt.event;
    cout << " with " << m_AdcChannelNamedTools.size() << " tools." << endl;
  }
  Index nfail = 0;
  ++state().nevtEnd;
  if ( state().nevtBegin != state().nevtEnd ) {
    cout << myname << "WARNING: Event counts are inconsistent: " << state().nevtBegin
         << " != " << state().nevtEnd << endl;
  }
  if ( m_AdcChannelNamedTools.size() ) {
    for ( NamedTool nt : m_AdcChannelNamedTools ) {
      DataMap ret = nt.tool->endEvent(devt);
      if ( ret.status() ) {
        ++nfail;
        cout << myname << "WARNING: Event finalization for tool " << nt.name
             << " failed for event " << devt.event << " with status code " << ret.status() << endl;
      }
    }
  }
  return nfail;
}

//**********************************************************************

int ToolBasedRawDigitPrepService::
prepare(detinfo::DetectorClocksData const& clockData,
        AdcChannelDataMap& datamap,
        std::vector<recob::Wire>* pwires, WiredAdcChannelDataMap* pintStates) const {
  const string myname = "ToolBasedRawDigitPrepService:prepare: ";
  // Loop over tools.
  ++state().ncall;
  if ( m_LogLevel >= 2 ) cout << myname << "Processing " << datamap.size() << " channels with "
                              << m_AdcChannelNamedTools.size() << " tools." << endl;
  if ( m_AdcChannelNamedTools.size() ) {
    Index itoo = 0;
    for ( NamedTool nt : m_AdcChannelNamedTools ) {
      if ( m_LogLevel >= 3 ) cout << myname << "  Running tool " << nt.name << endl;
      bool useCallgrind = m_cgset.count(nt.name);
      if ( useCallgrind ) {
        CALLGRIND_START_INSTRUMENTATION;
        CALLGRIND_TOGGLE_COLLECT;
      }
      auto start = Clock::now();
      DataMap ret = nt.tool->updateMap(datamap);
      auto stop = Clock::now();
      if ( useCallgrind ) {
        CALLGRIND_TOGGLE_COLLECT;
        CALLGRIND_STOP_INSTRUMENTATION;
      }
      Duration dtim =stop - start;
      state().toolTimes[itoo] += dtim;
      if ( ret ) {
        cout << myname << "WARNING: Tool " << nt.name << " failed";
        if ( ret.haveIntVector("failedChannels") ) {
          Index ncha = ret.getIntVector("failedChannels").size();
          cout << " for " << ncha << " channel" << (ncha == 1 ? "" : "s");
        }
        cout << " with error " << ret.status();
        if ( ret.haveIntVector("failedCodes") ) {
          Index ncod = ret.getIntVector("failedCodes").size();
          if ( ncod > 1 ) cout << " and "  << ncod - 1 << " other errors";
        }
        cout << "." << endl;
      }
      if ( m_LogLevel >= 4 ) {
        cout << myname << "----------------------------" << endl;
        ret.print();
        cout << myname << "----------------------------" << endl;
      }
      ++itoo;
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
       if ( m_cgset.count(nm.name) ) cout << " (callgrind enabled)";
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
