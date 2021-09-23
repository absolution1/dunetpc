// TpcToolBasedRawDigitPrepService_service.cc

#include "TpcToolBasedRawDigitPrepService.h"
#include "dune/ArtSupport/DuneToolManager.h"
#include "dune/DuneInterface/Tool/TpcDataTool.h"
#include "dune/DuneInterface/Service/AdcWireBuildingService.h"
#include "art/Framework/Services/Registry/ServiceDefinitionMacros.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
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

TpcToolBasedRawDigitPrepService::
TpcToolBasedRawDigitPrepService(fhicl::ParameterSet const& pset, art::ActivityRegistry&)
: m_LogLevel(pset.get<int>("LogLevel")),
  m_DoWires(pset.get<bool>("DoWires")),
  m_ToolNames(pset.get<vector<string>>("ToolNames")),
  m_CallgrindToolNames(pset.get<vector<string>>("CallgrindToolNames")),
  m_pWireBuildingService(nullptr),
  m_cgset(m_CallgrindToolNames.begin(), m_CallgrindToolNames.end()),
  m_pstate(new State) {
  const string myname = "TpcToolBasedRawDigitPrepService::ctor: ";
  pset.get_if_present<int>("LogLevel", m_LogLevel);
  // Fetch the tools.
  if ( m_ToolNames.size() ) {
    DuneToolManager* ptm = DuneToolManager::instance("");
    if ( ptm == nullptr ) {
      cout << myname << "ERROR: Unable to retrieve tool manaager." << endl;
    } else {
      for ( string tname : m_ToolNames ) {
        if ( m_LogLevel ) cout << myname << "     Fetching " << tname << endl;
        TpcDataToolPtr ptool = ptm->getPrivate<TpcDataTool>(tname);
        NamedTool nt(tname, ptool.get());
        if ( nt.tool ) {
          if ( m_LogLevel ) cout << myname << "    Found tool " << tname << " @ " << nt.tool << endl;
          m_TpcDataTools.push_back(std::move(ptool));
          m_TpcDataNamedTools.push_back(nt);
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

TpcToolBasedRawDigitPrepService::~TpcToolBasedRawDigitPrepService() {
  const string myname = "TpcToolBasedRawDigitPrepService:dtor: ";
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

int TpcToolBasedRawDigitPrepService::beginEvent(const DuneEventInfo& devt) const {
  const string myname = "TpcToolBasedRawDigitPrepService:beginEvent: ";
  if ( m_LogLevel >= 2 ) {
    cout << myname << "Begin processing run " << devt.runString();
    cout << " event " << devt.event;
    cout << " with " << m_TpcDataNamedTools.size() << " tools." << endl;
  }
  if ( state().nevtBegin != state().nevtEnd ) {
    cout << myname << "WARNING: Event counts are inconsistent: " << state().nevtBegin
         << " != " << state().nevtEnd << endl;
  }
  ++state().nevtBegin;
  Index nfail = 0;
  if ( m_TpcDataNamedTools.size() ) {
    for ( NamedTool nt : m_TpcDataNamedTools ) {
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

int TpcToolBasedRawDigitPrepService::endEvent(const DuneEventInfo& devt) const {
  const string myname = "TpcToolBasedRawDigitPrepService:endEvent: ";
  if ( m_LogLevel >= 2 ) {
    cout << myname << "End processing run " << devt.runString();
    cout << " event " << devt.event;
    cout << " with " << m_TpcDataNamedTools.size() << " tools." << endl;
  }
  Index nfail = 0;
  ++state().nevtEnd;
  if ( state().nevtBegin != state().nevtEnd ) {
    cout << myname << "WARNING: Event counts are inconsistent: " << state().nevtBegin
         << " != " << state().nevtEnd << endl;
  }
  if ( m_TpcDataNamedTools.size() ) {
    for ( NamedTool nt : m_TpcDataNamedTools ) {
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

int TpcToolBasedRawDigitPrepService::
prepare(detinfo::DetectorClocksData const& clockData,
        AdcChannelDataMap& datamap,
        std::vector<recob::Wire>* pwires, WiredAdcChannelDataMap* pintStates) const {
  const string myname = "TpcToolBasedRawDigitPrepService:prepare: ";
  // Loop over tools.
  ++state().ncall;
  if ( m_LogLevel >= 2 ) cout << myname << "Processing " << datamap.size() << " channels with "
                              << m_TpcDataNamedTools.size() << " tools." << endl;
  if ( m_TpcDataNamedTools.size() ) {
    Index itoo = 0;
    for ( NamedTool nt : m_TpcDataNamedTools ) {
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

std::ostream& TpcToolBasedRawDigitPrepService::
print(std::ostream& out, std::string prefix) const {
  out << prefix << "TpcToolBasedRawDigitPrepService:"                      << endl;
  out << prefix << "                    LogLevel: " << m_LogLevel             << endl;
  out << prefix << "                     DoWires: " << m_DoWires              << endl;
  if ( m_TpcDataNamedTools.size() ) {
    cout << prefix << "     ADC channel tools:";
    for ( const NamedTool& nm : m_TpcDataNamedTools ) {
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

DEFINE_ART_SERVICE_INTERFACE_IMPL(TpcToolBasedRawDigitPrepService, RawDigitPrepService)

//**********************************************************************
