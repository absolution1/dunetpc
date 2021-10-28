// AdcEventViewer_tool.cc

#include "AdcEventViewer.h"
#include "dune/DuneInterface/Tool/AdcChannelStringTool.h"
#include "dune/ArtSupport/DuneToolManager.h"
#include "dune/DuneInterface/Tool/RunDataTool.h"
#include "dune/DuneInterface/Tool/IndexRangeTool.h"
#include "dune/DuneCommon/Utility/TPadManipulator.h"
#include "dune/DuneCommon/Utility/StringManipulator.h"
#include "dune/ArtSupport/DuneToolManager.h"
#include "TH1.h"
#include <iostream>
#include <sstream>
#include <iomanip>

using std::string;
using std::to_string;
using std::cout;
using std::endl;
using std::istringstream;
using std::ostringstream;
using fhicl::ParameterSet;
using std::setw;
using std::fixed;
using std::setprecision;

using Index = AdcEventViewer::Index;

//**********************************************************************
// Local methods.
//**********************************************************************

namespace {

struct VarInfo {
  string name;
  string vname;
  string label;
  string unit;
  VarInfo(string aname, const IndexRange& cr, string clockUnit);
  bool isValid() const { return vname.size(); }
};

// Extract variable info from a name.
VarInfo::VarInfo(string aname, const IndexRange& cr, string clockUnit) : name(aname) {
  if ( name.find("nfemb") != string::npos ) {
    vname = "nfemb";
    label = "FEMB count";
  } else if ( name.find("event") != string::npos ) {
    vname = "event";
    label = "Event";
  } else if ( name.find("clock") != string::npos ) {
    vname = "clock";
    label = "Timing clock";
    unit = clockUnit;
  } else if ( name.find("time") != string::npos ) {
    vname = "time";
    label = "Time";
    unit = "sec";
  } else if ( name.find("rmPedPower") != string::npos ) {
    vname = "rmPedPower";
    label = "Pedestal noise RMS";
    unit = "ADC counts";
  } else if ( name.find("meanPed") != string::npos ) {
    vname = "meanPed";
    label = "Pedestal mean";
    unit = "ADC counts";
  }
}

string formatValue(float val, string styp) {
  ostringstream ssval;
  ssval << val;
  string sval = ssval.str();
  if ( styp == "event" ) {
    while ( sval.size() < 8 ) sval = "0" + sval;
  }
  return sval;
}

}  // end unnamed namespace

//**********************************************************************
// Subclass methods.
//**********************************************************************

std::string AdcEventViewer::State::runString() const {
  const string myname = "AdcEventViewer::State::runString: ";
  if ( eventInfo.isValid() ) return eventInfo.runString();
  if ( runSet.size() == 0 ) return "<no-runs>";
  if ( runSet.size() > 1 ) {
    cout << myname << "WARNING: Multiple runs have been processed:";
    bool first = true;
    for ( Index irun : runSet ) {
      if ( first ) first = false;
      else cout << ", ";
      cout << irun;
    }
    cout << endl;
    return "<multiple-runs>";
  }
  DuneEventInfo evi;
  evi.run = *runSet.begin();
  return evi.runString();
}

//**********************************************************************
// Class methods.
//**********************************************************************

AdcEventViewer::AdcEventViewer(fhicl::ParameterSet const& ps)
: m_LogLevel(ps.get<int>("LogLevel")),
  m_EventHists(ps.get<NameVector>("EventHists")),
  m_EventGraphs(ps.get<NameVector>("EventGraphs")),
  m_ChannelRanges(ps.get<NameVector>("ChannelRanges")),
  m_ChannelRangeLabel(ps.get<Name>("ChannelRangeLabel")),
  m_ClockUnit(ps.get<Name>("ClockUnit")),
  m_ClockRate(ps.get<double>("ClockRate")),
  m_state(new AdcEventViewer::State) {
  const string myname = "AdcEventViewer::ctor: ";
  // Fetch the channel ranges.
  bool toolNotFound = false;
  const IndexRangeTool* pcrt = nullptr;
  DuneToolManager* ptm = DuneToolManager::instance();
  for ( Name crn : m_ChannelRanges.size() ? m_ChannelRanges : NameVector(1, "") ) {
    if ( crn.size() == 0 || crn == "all" ) {
      m_crs.emplace_back("all", 0, 0, "All");
    } else {
      if ( pcrt == nullptr && !toolNotFound ) {
        pcrt = ptm->getShared<IndexRangeTool>("channelRanges");
        if ( pcrt == nullptr ) {
          cout << myname << "ERROR: IndexRangeTool not found: channelRanges" << endl;
        }
      }
      if ( pcrt != nullptr ) {
        IndexRange ran = pcrt->get(crn);
        if ( ran.isValid() ) {
          m_crs.push_back(ran);
        } else {
          cout << myname << "WARNING: Channel range not found: " << crn << endl;
        }
      }
    }
  }
  // Display the configuration.
  if ( m_LogLevel>= 1 ) {
    cout << myname << "         LogLevel: " << m_LogLevel << endl;
    cout << myname << "       EventHists: [" << endl;
    for ( Index ievh=0; ievh<m_EventHists.size(); ++ievh ) {
      cout << myname << "                     " << m_EventHists[ievh] << endl;
    }
    cout << myname << "                   " << "]" << endl;
    cout << myname << "      EventGraphs: [" << endl;
    for ( Index ievg=0; ievg<m_EventGraphs.size(); ++ievg ) {
      cout << myname << "                     " << m_EventGraphs[ievg] << endl;
    }
    cout << myname << "                   " << "]" << endl;
    cout << myname << "       ChannelRanges: [";
    bool first = true;
    for ( const IndexRange& ran : m_crs ) {
      if ( ! first ) cout << ", ";
      else first = false;
      cout << ran.name;
    }
    cout << "]" << endl;
    cout << myname << "   ChannelRangeLabel: " << m_ChannelRangeLabel << endl;
    cout << myname << "           ClockUnit: " << m_ClockUnit << endl;
    cout << myname << "           ClockRate: " << m_ClockRate << " tick/sec" << endl;
  }
  state().firstClock = 0;
  state().minClock = 0;
  const string::size_type& npos = string::npos;
  // Make sure clock unit is valid.
  if ( m_ClockUnit != "tick" &&
       m_ClockUnit != "ktick" &&
       m_ClockUnit != "Mtick" &&
       m_ClockUnit != "Gtick" ) {
    cout << myname << "WARNING: invalid clock unit changed to tick" << endl;
    m_ClockUnit = "tick";
  }
  // Loop over channel ranges
  for ( const IndexRange& cr : m_crs ) {
    // Create the substate;
    Name crn = cr.name;
    if ( state().crstates.find(crn) != state().crstates.end() ) {
      cout << myname << "WARNING: Skipping duplicate CR name: " << crn << endl;
      continue;
    }
    ChannelRangeState& crstate = state().crstates[crn];
    // Create the histograms.
    for ( string hspec : m_EventHists ) {
      string basename;
      int nbin = 0;
      float xmin = 0.0;
      float xmax = 0.0;
      bool ok = false;
      string::size_type ipos = 0;
      string::size_type jpos = 0;
      istringstream ssin;
      jpos = hspec.find(":");
      if ( jpos != npos ) {
        basename = hspec.substr(ipos, jpos-ipos);
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
      string hname = basename;
      if ( crn != "all " ) hname += "_" + crn;
      string vname;
      if ( basename.find("nfemb") != string::npos ) {
        vname = "nfemb";
      } else if ( basename.find("rmPedPower") != string::npos ) {
        vname = "rmPedPower";
      } else if ( basename.find("meanPed") != string::npos ) {
        vname = "meanPed";
      } else {
        cout << myname << "ERROR: No variable for histogram name " << hname << endl;
        continue;
      }
      VarInfo vinfo(vname, cr, m_ClockUnit);
      string sttl = vinfo.label;
      sttl += ";" + vinfo.label;
      Name unit = vinfo.unit;
      if ( unit.size() ) sttl += " [" + unit + "]";
      sttl += ";# event";
      if ( m_LogLevel >= 2 ) {
        cout << myname << "Creating in histogram " << hname << ", nbin=" << nbin
             << ", range=(" << xmin << ", " << xmax << ")" << endl;
      }
      TH1* ph = new TH1F(hname.c_str(), sttl.c_str(), nbin, xmin, xmax);
      ph->SetDirectory(nullptr);
      ph->SetStats(0);
      ph->SetLineWidth(2);
      crstate.hists.push_back(ph);
      if ( m_LogLevel>= 1 ) cout << myname << "Created histogram " << hname << endl;
    }
    for ( string gspec : m_EventGraphs ) {
      if ( m_LogLevel >= 2 ) cout << myname << "Building graph info" << endl;
      string xname;
      string yname;
      float xmin = 0.0;
      float xmax = 0.0;
      float xoff = 0.0;
      float ymin = 0.0;
      float ymax = 0.0;
      float yoff = 0.0;
      string::size_type ipos = 0;
      string::size_type jpos = 0;
      istringstream ssin;
      jpos = gspec.find(":");
      // x name
      if ( jpos != npos ) {
        xname = gspec.substr(ipos, jpos-ipos);
        ipos = jpos + 1;
        jpos = gspec.find(":", ipos);
        if ( m_LogLevel >= 2 ) cout << myname << "  xname: " << xname << endl;
      }
      // x range
      if ( jpos != npos && jpos > ipos ) {
        // Range from job range tool:   name:ranName:xoff
        if ( isalpha(gspec[ipos]) ) {
          if ( m_LogLevel >= 2 ) cout << myname << "  Taking x range from tool." << endl;
          istringstream ssnran(gspec.substr(ipos, jpos-ipos));
          Name ranName = gspec.substr(ipos, jpos-ipos);
          const IndexRangeTool* pcrt = ptm->getShared<IndexRangeTool>("jobRanges");
          if ( pcrt == nullptr ) {
            cout << myname << "ERROR: Unable to find index range tool jobRanges" << endl;
          } else {
            IndexRange ran = pcrt->get(ranName);
            if ( ran.isValid() ) {
              xmin = ran.begin;
              xmax = ran.end;
            } else {
              cout << myname << "ERROR: Unable to find job range " << ranName << endl;
            }
          }
          ipos = jpos + 1;
          jpos = gspec.find(":", ipos);
          if ( jpos != npos && jpos > ipos ) {
            istringstream ssxoff(gspec.substr(ipos, jpos-ipos));
            ssxoff >> xoff;
          }
        // Explicit range
        } else {
          if ( m_LogLevel >= 2 ) cout << myname << "  Using explicit x range." << endl;
          istringstream ssxmin(gspec.substr(ipos, jpos-ipos));
          ssxmin >> xmin;
          ipos = jpos + 1;
          jpos = gspec.find(":", ipos);
          if ( jpos != npos && jpos > ipos ) {
            istringstream ssxmax(gspec.substr(ipos, jpos-ipos));
            ssxmax >> xmax;
          }
        }
      }
      ipos = jpos + 1;
      jpos = gspec.find(":", ipos);
      // y name
      if ( jpos != npos && jpos > ipos ) {
        yname = gspec.substr(ipos, jpos-ipos);
        ipos = jpos + 1;
        jpos = gspec.find(":", ipos);
        if ( m_LogLevel >= 2 ) cout << myname << "  yname: " << yname << endl;
      }
      // y range
      if ( jpos != npos && jpos > ipos ) {
        // Range from job range tool:   name:ranName:yoff
        if ( isalpha(gspec[ipos]) ) {
          if ( m_LogLevel >= 2 ) cout << myname << "  Taking y range from tool." << endl;
          istringstream ssnran(gspec.substr(ipos, jpos-ipos));
          Name ranName = gspec.substr(ipos, jpos-ipos);
          cout << myname << "Taking y-range from job range " << ranName << endl;
          const IndexRangeTool* pcrt = ptm->getShared<IndexRangeTool>("jobRanges");
          if ( pcrt == nullptr ) {
            cout << myname << "ERROR: Unable to find index range tool jobRanges" << endl;
          } else {
            IndexRange ran = pcrt->get(ranName);
            if ( ran.isValid() ) {
              ymin = ran.begin;
              ymax = ran.end;
            } else {
              cout << myname << "ERROR: Unable to fund job range " << ranName << endl;
            }
          }
          ipos = jpos + 1;
          jpos = gspec.find(":", ipos);
          if ( jpos != npos && jpos > ipos ) {
            istringstream ssyoff(gspec.substr(ipos, jpos-ipos));
            ssyoff >> yoff;
          }
        // Explicit range
        } else {
          if ( m_LogLevel >= 2 ) cout << myname << "  Using explicit y range." << endl;
          istringstream ssymin(gspec.substr(ipos, jpos-ipos));
          ssymin >> ymin;
          ipos = jpos + 1;
          jpos = gspec.find(":", ipos);
          if ( jpos > ipos ) {
            istringstream ssymax(gspec.substr(ipos, jpos-ipos));
            ssymax >> ymax;
          }
        }
      }
      VarInfo xvin(xname, cr, m_ClockUnit);
      VarInfo yvin(yname, cr, m_ClockUnit);
      bool ok = xvin.isValid() && yvin.isValid();
      if ( ! ok ) {
        cout << myname << "WARNING: Invalid graph configuration string: " << gspec << endl;
        continue;
      }
      string sttl = yvin.label + " vs. " + xvin.label;
      if ( m_LogLevel >= 2 ) {
        cout << myname << "Creating graph of " << yname;
        if ( ymax > ymin ) cout << " range=(" << ymin - yoff << ", " << ymax + yoff << ")";
        else if ( ymax >= 0.0 ) cout << " range from data with border " << yoff << endl;
        cout << " vs. " << xname;
        if ( xmax > xmin ) cout << " range=(" << xmin << ", " << xmax << ")";
        else if ( xmax >= 0.0 ) cout << " range from data with border " << xoff << endl;
        cout << endl;
      } 
      crstate.graphs.emplace_back(xname, xvin.label, xvin.unit, xmin, xmax, xoff,
                                  yname, yvin.label, yvin.unit, ymin, ymax, yoff );
      if ( m_LogLevel>= 1 ) cout << myname << "Created graph info for " << yname << " vs. " << xname << endl;
    }
  }  // end loop over channel ranges
}

//**********************************************************************

AdcEventViewer::~AdcEventViewer() {
  const string myname = "AdcEventViewer::dtor: ";
  if ( state().beginEventCount == 0 ) endEventState(state().eventInfo);
  cout << myname << " Begin event count: " <<  state().beginEventCount << endl;
  cout << myname << "Unique event count: " <<  state().events.size() << endl;
  displayHists();
  displayGraphs();
  cout << myname << "Exiting." << endl;
}

//**********************************************************************

DataMap AdcEventViewer::beginEvent(const DuneEventInfo& devt) const {
  const string myname = "AdcEventViewer::beginEvent: ";
  DataMap ret;
  if ( m_LogLevel >= 2 ) {
    cout << myname << "Begin processing event " << devt.event << " for run " << devt.run;
    if ( devt.subRun ) cout << "-" << devt.subRun;
    cout << endl;
  }
  if ( state().haveEvent() ) {
    cout << myname << "WARNING: Begin called wth event " << devt.event
         << " during processing of event " << state().eventInfo.event
         << ". Ending current event." << endl;
    endEvent(state().eventInfo);
  }
  ++state().beginEventCount;
  beginEventState(devt);
  return ret;
}

//**********************************************************************

DataMap AdcEventViewer::endEvent(const DuneEventInfo& devt) const {
  const string myname = "AdcEventViewer::endEvent: ";
  if ( m_LogLevel >= 2 ) {
    cout << myname << "End processing event " << devt.event << " for run " << devt.run;
    if ( devt.subRun ) cout << "-" << devt.subRun;
    cout << endl;
  }
  DataMap ret;
  endEventState(devt);
  return ret;
}

//**********************************************************************

DataMap AdcEventViewer::view(const AdcChannelData& acd) const {
  const string myname = "AdcEventViewer::view: ";
  DataMap res;
  if ( m_LogLevel >= 4 ) cout << myname << "Processing channel " << acd.channel()
                              << " (FEMB " << acd.fembID() << ")"
                              << " in run " << acd.run() << " event " << acd.event() << endl;
  if ( acd.event() != state().event() || acd.run() != state().run() ) {
    if ( state().beginEventCount ) {
      cout << myname << "ERROR: Run:event from data " << acd.run() << ":" << acd.event()
           << " does not match that from beginEvent "
           << state().runString() << ":" << state().event() << endl;
      return res.setStatus(1);
    }
    beginEventState(acd.getEventInfo());
  }
  for ( const IndexRange& cr : m_crs ) {
    if ( cr.name != "all" && !cr.contains(acd.channel()) ) continue;
    ChannelRangeState& crstate = state().crstates[cr.name];
    crstate.fembIDSet.insert(acd.fembID());
    ++crstate.nchan;
    crstate.pedSum += acd.pedestal;
    float pedNoise = acd.pedestalRms;
    crstate.pedPower += pedNoise*pedNoise;
  }
  return res;
}

//**********************************************************************

DataMap AdcEventViewer::viewMap(const AdcChannelDataMap& acds) const {
  const string myname = "AdcEventViewer::viewMap: ";
  if ( m_LogLevel >= 3 ) {
    Index ncha = acds.size();
    cout << myname << "Processing " << ncha << " channel" << (ncha == 1 ? "" : "s");
    if ( ncha > 0 ) {
      const AdcChannelData& acd = acds.begin()->second;
      cout << " in run " << acd.run() << " event " << acd.event();
      cout << ". Count is " << state().beginEventCount << "." << endl;
    }
  }
  DataMap ret;
  if ( state().beginEventCount == 0 ) {
    if ( acds.size() == 0 ) {
      if ( m_LogLevel >=2 ) cout << myname << "Skipping group with no data" << endl;
      return ret;
    }
    const AdcChannelData& acd = acds.begin()->second;
    if ( acd.event() != state().event() ) beginEventState(acd.getEventInfo());
  }
  ++state().ngroup;
  for ( const AdcChannelDataMap::value_type& iacd : acds ) view(iacd.second);
  return ret;
}

//**********************************************************************

AdcEventViewer::ChannelRangeState& AdcEventViewer::crstate(Name crn) const {
  const string myname = "AdcEventViewer::crstate: ";
  ChannelRangeStates::iterator icrs = state().crstates.find(crn);
  if ( icrs == state().crstates.end() ) {
    cout << myname << "ERROR: There is no substate for channel range " << crn << endl;
    static ChannelRangeState crsbad;
    return crsbad;
  }
  return icrs->second;
}

//**********************************************************************

void AdcEventViewer::beginEventState(const DuneEventInfo& devt) const {
  const string myname = "AdcEventViewer::beginEventState: ";
  if ( state().haveEvent() ) {
    cout << "WARNING: Begin called during event processing. Ending current event." << endl;
    endEventState(devt);
  }
  if ( m_LogLevel >= 4 ) cout << myname << "Starting run" << devt.runString()
                              << " event " << devt.event << endl;
  if ( state().haveEvent() && state().run() != devt.run ) {
    cout << myname << "ERROR: change in run number is not (yet) supported." << endl;
    return;
  }
  Index ievt = devt.event;
  Index irun = devt.run;
  state().eventInfo = devt;
  state().events.push_back(ievt);
  state().eventSet.insert(ievt);
  state().runSet.insert(irun);
  LongIndex clk = devt.triggerClock;
  if ( state().firstClock == 0 ) state().firstClock = clk;
  if ( clk < state().minClock ) state().minClock = clk;
  state().ngroup = 0;
  for ( ChannelRangeStates::iterator icr=state().crstates.begin();
        icr!=state().crstates.end(); ++icr ) {
    ChannelRangeState& crstate = icr->second;
    crstate.fembIDSet.clear();
    crstate.nchan = 0;
    crstate.pedSum = 0.0;
    crstate.pedPower = 0.0;
  }
}

//**********************************************************************

void AdcEventViewer::endEventState(const DuneEventInfo& devt) const {
  const string myname = "AdcEventViewer::endEventState: ";
  if ( ! devt.isValid() ) {
    cout << myname << "ERROR: End called without an active event." << endl;
    return;
  }
  if ( devt != state().eventInfo ) {
    cout << myname << "ERROR: End called without incorrect event info." << endl;
    return;
  }
  if ( state().event() == 0 ) return;
  Index nevt = state().events.size();
  Index ndup = nevt - state().eventSet.size();
  const int w = 5;
  long iclk = state().clock();
  iclk -= state().firstClock;
  // Plotted clock has zero 1 sec after the clock for the first event.
  float xclk = iclk + m_ClockRate;
  float time = m_ClockRate == 0.0 ? 0.0 : xclk/m_ClockRate;
  Name clkunit = m_ClockUnit;
  if      ( clkunit == "ktick" ) xclk *= 0.001;
  else if ( clkunit == "Mtick" ) xclk *= 0.000001;
  else if ( clkunit == "Gtick" ) xclk *= 0.000000001;
  if ( m_LogLevel >= 2 ) {
    cout << myname << "                 run: " << setw(w) << state().runString() << endl;
    cout << myname << "               event: " << setw(w) << state().event() << endl;
    cout << myname << "            # events: " << setw(w) << state().events.size() << endl;
    cout << myname << "  # duplicate events: " << setw(w) << ndup << endl;
    cout << myname << "           Min clock: " << setw(w) << state().minClock << endl;
    cout << myname << "               clock: " << setw(w+3) << fixed << setprecision(2) << xclk << " " << clkunit << endl;
    cout << myname << "                time: " << setw(w+3) << fixed << setprecision(2) << time << " sec" << endl;
    cout << myname << "            # groups: " << setw(w) << state().ngroup << endl;
  }
  for ( const IndexRange& cr : m_crs ) {
    ChannelRangeState& crstate = state().crstates[cr.name];
    Index nfmb = crstate.fembIDSet.size();
    Index nchn = crstate.nchan;
    float meanPed = nchn > 0 ? crstate.pedSum/nchn : 0.0;
    float meanPedPower = nchn > 0 ? crstate.pedPower/nchn : 0.0;
    float rmPedPower = sqrt(meanPedPower);
    double chanPerFemb = nfmb > 0 ? double(nchn)/nfmb : 0.0;
    if ( m_LogLevel >= 2 ) {
      cout << myname << "----- Channel range " << cr.label() << ":" << endl;
      cout << myname << "             # FEMBs: " << setw(w) << nfmb << endl;
      cout << myname << "          # channels: " << setw(w) << nchn << endl;
      cout << myname << "     # channels/FEMB: " << setw(w+5) << fixed << setprecision(4) << chanPerFemb << endl;
      cout << myname << "       Mean pedestal: " << setw(w+3) << fixed << setprecision(2) << meanPed << endl;
      cout << myname << "   RM pedestal power: " << setw(w+3) << fixed << setprecision(2) << rmPedPower << endl;
    }
    for ( TH1* ph : crstate.hists ) {
      string name = ph->GetName();
      if ( name.find("nfemb") != string::npos ) ph->Fill(nfmb);
      else if ( name.find("rmPedPower") != string::npos ) ph->Fill(rmPedPower);
      else if ( name.find("meanPed") != string::npos ) ph->Fill(meanPed);
      else {
        cout << myname << "ERROR: No variable for histogram name " << name << endl;
      }
    }
    for ( GraphInfo& gin : crstate.graphs ) {
      gin.add("event", state().event());
      gin.add("clock", xclk);
      gin.add("time", time);
      gin.add("nfemb", nfmb);
      gin.add("rmPedPower", rmPedPower);
      gin.add("meanPed", meanPed);
    }
  }
  state().eventInfo.clear();
}

//**********************************************************************

void AdcEventViewer::displayHists() const {
  const string myname = "AdcEventViewer::displayHists: ";
  string sttlSuf = " for run " + state().runString();
  Index nevt = state().eventSet.size();
  if ( nevt == 0 ) sttlSuf += " with no events.";
  else if ( nevt == 1 ) sttlSuf += " event " + to_string(*state().eventSet.begin());
  else sttlSuf += " events " + to_string(*state().eventSet.begin()) + "-" + to_string(*state().eventSet.rbegin());
  string evstring;
  if ( state().eventSet.size() ) {
    evstring += "_event" + formatValue(*state().eventSet.cbegin(), "event");
    if ( state().eventSet.size() > 1 ) {
      evstring += "-" + formatValue(*state().eventSet.crbegin(), "event");
    }
  }
  for ( const IndexRange& cr : m_crs ) {
    ChannelRangeState& crstate = state().crstates[cr.name];
    Index nplt = crstate.hists.size();
    if ( m_LogLevel >= 1 ) cout << myname << "Creating " << nplt << " plot"
                                << (nplt == 1 ? "" : "s") << sttlSuf
                                << (nplt > 0 ? ":" : "") << endl;
    for ( TH1* ph : crstate.hists ) {
      TPadManipulator man;
      man.add(ph);
      string sttl = ph->GetTitle() + sttlSuf;
      string crlab = m_ChannelRangeLabel;
      StringManipulator smlab(crlab, false);
      smlab.replace("%CRNAME%", cr.name);
      smlab.replace("%CRLABEL%", cr.label());
      smlab.replace("%CRLABEL1%", cr.label(1));
      smlab.replace("%CRLABEL2%", cr.label(2));
      if ( crlab.size() ) sttl += " " + crlab;
      man.setTitle(sttl.c_str());
      
      string fname = string("eviewh_") + ph->GetName() + "_run" + state().runString()
                     + evstring + "_" + cr.name + ".{png,tpad}";
      man.showUnderflow();
      man.showOverflow();
      man.addAxis();
      man.print(fname);
      if ( m_LogLevel >= 1 ) cout << myname << "  " << fname << endl;
    }
  }
}

//**********************************************************************

void AdcEventViewer::displayGraphs() const {
  const string myname = "AdcEventViewer::displayGraphs: ";
  string sttlSufBase = " for run " + state().runString();
  for ( const IndexRange& cr : m_crs ) {
    ChannelRangeState& crstate = state().crstates[cr.name];
    Index nplt = crstate.graphs.size();
    Index nevt = state().eventSet.size();
    string sttlSuf = sttlSufBase;
    if ( nevt == 0 ) sttlSuf += " with no events.";
    else if ( nevt == 1 ) sttlSuf += " event " + to_string(*state().eventSet.begin());
    else sttlSuf += " events " + to_string(*state().eventSet.begin()) + "-" + to_string(*state().eventSet.rbegin());
    if ( m_LogLevel >= 1 ) cout << myname << "Creating " << nplt << " graph"
                                << (nplt == 1 ? "" : "s") << sttlSuf
                                << (nplt > 0 ? ":" : "") << endl;
    string crlab = m_ChannelRangeLabel;
    StringManipulator smlab(crlab, false);
    smlab.replace("%CRNAME%", cr.name);
    smlab.replace("%CRLABEL%", cr.label());
    smlab.replace("%CRLABEL1%", cr.label(1));
    smlab.replace("%CRLABEL2%", cr.label(2));
    if ( crlab.size() ) sttlSuf += " " + crlab;
    for ( GraphInfo& gin : crstate.graphs ) {
      if ( m_LogLevel >= 1 ) cout << myname << "Creating graph of " << gin.vary << " vs. " << gin.varx << endl;
      Index npt = gin.xvals.size();
      if ( npt == 0 ) {
        cout << myname << "Skipping graph with no entries." << endl;
        continue;
      }
      if ( gin.yvals.size() != npt ) {
        cout << myname << "Skipping graph with inconsistent entries." << endl;
        continue;
      }
      TGraph* pg = new TGraph(npt, &gin.xvals[0], &gin.yvals[0]);
      pg->GetXaxis()->SetTitle(gin.xAxisLabel().c_str());
      pg->GetYaxis()->SetTitle(gin.yAxisLabel().c_str());
      pg->SetMarkerStyle(2);
      double xmin = gin.xmin;
      double xmax = gin.xmax;
      double xoff = gin.xoff;
      double ymin = gin.ymin;
      double ymax = gin.ymax;
      double yoff = gin.yoff;
      bool setXLimitsFromData = pg->GetN() && xmax == xmin && xmax >= 0.0;
      bool setYLimitsFromData = pg->GetN() && ymax == ymin && ymax >= 0.0;
      if ( setXLimitsFromData || setYLimitsFromData ) {
        double xdmin = 0.0;
        double xdmax = 0.0;
        double ydmin = 0.0;
        double ydmax = 0.0;
        double xpt, ypt;
        for ( int ipt=0; ipt<pg->GetN(); ++ipt ) {
          pg->GetPoint(ipt, xpt, ypt);
          if ( ipt == 0 ) {
            xdmin = xpt;
            xdmax = xpt;
            ydmin = ypt;
            ydmax = ypt;
          } else {
            if ( xpt < xdmin ) xdmin = xpt;
            if ( xpt > xdmax ) xdmax = xpt;
            if ( ypt < ydmin ) ydmin = ypt;
            if ( ypt > ydmax ) ydmax = ypt;
          }
        }
        if ( setXLimitsFromData ) {
          xoff = xmax;
          xmin = xdmin;
          xmax = xdmax;
          if ( m_LogLevel >= 2 ) {
            cout << myname << "  Setting graph X range from data (" << xdmin << ", " << xdmax << ") to ("
                 << xmin-xoff << ", " << xmax+xoff << ")" << endl;
          }
        }
        if ( setYLimitsFromData ) {
          yoff = ymax;
          ymin = ydmin;
          ymax = ydmax;
          if ( m_LogLevel >= 2 ) {
            cout << myname << "  Setting graph Y range from data (" << ydmin << ", " << ydmax << ") to ("
                 << ymin-yoff << ", " << ymax+yoff << ")" << endl;
          }
        }
      }
      if ( xmax > xmin ) pg->GetXaxis()->SetRangeUser(xmin-xoff, xmax+xoff);
      if ( ymax > ymin ) pg->GetYaxis()->SetRangeUser(ymin-yoff, ymax+yoff);
      TPadManipulator man(1400, 500);
      man.add(pg, "P");
      string sttl = gin.ylab + " vs. " + gin.xlab + sttlSuf;
      man.setTitle(sttl.c_str());
      man.showUnderflow();
      man.showOverflow();
      man.addAxis();
      man.setGridX();
      man.setGridY();
      if ( xmax > xmin ) man.setRangeX(xmin-xoff, xmax+xoff);
      if ( ymax > ymin ) man.setRangeY(ymin-yoff, ymax+yoff);
      ostringstream ssfname;
      ssfname << "eviewg_" << gin.varx;
      if ( gin.xmax > gin.xmin ) ssfname << "-" << formatValue(gin.xmin, gin.varx)
                                         << "-" << formatValue(gin.xmax, gin.varx);
      ssfname << "_" << gin.vary;
      if ( gin.ymax > gin.ymin ) ssfname << "-" << formatValue(gin.ymin, gin.vary)
                                         << "-" << formatValue(gin.ymax, gin.vary);
      ssfname << "_run" << state().runString();
      ssfname << "_" << cr.name;
      ssfname << ".{png,tpad}";
      Name fname = ssfname.str();
      man.print(fname);
      if ( m_LogLevel >= 1 ) cout << myname << "  " << fname << endl;
    }
  }
}

//**********************************************************************

DEFINE_ART_CLASS_TOOL(AdcEventViewer)
