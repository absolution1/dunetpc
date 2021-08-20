// AdcDetectorPlotter_tool.cc

#include "AdcDetectorPlotter.h"
#include <iostream>
#include <sstream>
#include "dune/DuneCommon/Utility/TPadManipulator.h"
#include "dune/DuneCommon/Utility/RootPalette.h"
#include "dune/DuneCommon/Utility/LineColors.h"
#include "dune/ArtSupport/DuneToolManager.h"
#include "dune/DuneInterface/Tool/AdcChannelStringTool.h"
#include "larcore/Geometry/Geometry.h"
#include "larevt/CalibrationDBI/Interface/ChannelStatusService.h"
#include "larevt/CalibrationDBI/Interface/ChannelStatusProvider.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "TColor.h"
#include "TStyle.h"
#include "TDirectory.h"
#include "TFile.h"
#include "TGraph.h"

using std::string;
using std::cout;
using std::cin;
using std::endl;
using std::ostringstream;

using Tick = AdcSignalVector::size_type;
using geo::GeometryCore;
using std::vector;

//**********************************************************************
// Local methods.
//**********************************************************************

namespace {

void initializeState(AdcDetectorPlotter::State& state, const AdcChannelData& acd) {
  state.reportCount = 0;
  state.channelCount = 0;
  state.run    = acd.run();
  state.subrun = acd.subRun();
  state.event  = acd.event();
  state.ppad.reset(nullptr);
}

}  // end unnamed namespace

//**********************************************************************
// Class methods.
//**********************************************************************

AdcDetectorPlotter::AdcDetectorPlotter(fhicl::ParameterSet const& ps)
: m_LogLevel(ps.get<int>("LogLevel")), 
  m_WireAngle(ps.get<float>("WireAngle")),
  m_DataType(ps.get<int>("DataType")),
  m_Tick0(ps.get<float>("Tick0")),
  m_DriftSpeed(ps.get<float>("DriftSpeed")),
  m_XMin(ps.get<float>("XMin")),
  m_XMax(ps.get<float>("XMax")),
  m_ZMin(ps.get<float>("ZMin")),
  m_ZMax(ps.get<float>("ZMax")),
  m_SignalThreshold(ps.get<float>("SignalThreshold")),
  m_SkipBadChannels(ps.get<bool>("SkipBadChannels")),
  m_ShowAllTicks(ps.get<bool>("ShowAllTicks")),
  m_FirstTick(ps.get<unsigned long>("FirstTick")),
  m_LastTick(ps.get<unsigned long>("LastTick")),
  m_ShowWires(ps.get<bool>("ShowWires")),
  m_ShowCathode(ps.get<bool>("ShowCathode")),
  m_ShowTpcSets(ps.get<IndexVector>("ShowTpcSets")),
  m_ShowGrid(ps.get<bool>("ShowGrid")),
  m_Title(ps.get<string>("Title")),
  m_PlotTitle(ps.get<string>("PlotTitle")),
  m_FileName(ps.get<string>("FileName")),
  m_pChannelStatusProvider(nullptr),
  m_state(new State) {
  const string myname = "AdcDetectorPlotter::ctor: ";
  DuneToolManager* ptm = DuneToolManager::instance();
  string stringBuilder = "adcStringBuilder";
  m_adcStringBuilder = ptm->getShared<AdcChannelStringTool>(stringBuilder);
  if ( m_adcStringBuilder == nullptr ) {
    cout << myname << "WARNING: AdcChannelStringTool not found: " << stringBuilder << endl;
  }
  if ( m_SkipBadChannels ) {
    if ( m_LogLevel >= 1 ) cout << myname << "Fetching channel status service." << endl;
    m_pChannelStatusProvider = &art::ServiceHandle<lariov::ChannelStatusService>()->GetProvider();
    if ( m_pChannelStatusProvider == nullptr ) {
      cout << myname << "WARNING: Channel status provider not found." << endl;
    }
  }
  if ( m_LogLevel ) {
    cout << myname << "Configuration: " << endl;
    cout << myname << "         LogLevel: " << m_LogLevel << endl;
    cout << myname << "        WireAngle: " << m_WireAngle << endl;
    cout << myname << "         DataType: " << m_DataType << endl;
    cout << myname << "            Tick0: " << m_Tick0 << endl;
    cout << myname << "       DriftSpeed: " << m_DriftSpeed << " cm/tick" << endl;
    cout << myname << "             XMin: " << m_XMin << " cm" << endl;
    cout << myname << "             XMax: " << m_XMax << " cm" << endl;
    cout << myname << "             ZMin: " << m_ZMin << " cm" << endl;
    cout << myname << "             ZMax: " << m_ZMax << " cm" << endl;
    cout << myname << "  SignalThreshold: " << m_SignalThreshold << endl;
    cout << myname << "  SkipBadChannels: " << (m_SkipBadChannels ? "true" : "false") << endl;
    cout << myname << "     ShowAllTicks: " << m_ShowAllTicks << endl;
    cout << myname << "        FirstTick: " << m_FirstTick << endl;
    cout << myname << "         LastTick: " << m_LastTick << endl;
    cout << myname << "        ShowWires: " << m_ShowWires << endl;
    cout << myname << "      ShowCathode: " << m_ShowCathode << endl;
    cout << myname << "      ShowTpcSets: [";
    bool first = true;
    for ( Index itps : m_ShowTpcSets ) {
      if ( first ) first = false;
      else cout << ", ";
      cout << itps;
    }
    cout << "]" << endl;
    cout << myname << "         ShowGrid: " << m_ShowGrid << endl;
    cout << myname << "            Title: " << m_Title << endl;
    cout << myname << "        PlotTitle: " << m_PlotTitle << endl;
    cout << myname << "         FileName: " << m_FileName << endl;
  }
  WireSelector& sel = getState()->sel;
  sel.selectWireAngle(m_WireAngle);
  sel.selectTpcSets(m_ShowTpcSets);
  //for ( Index itps : m_ShowTpcSets ) sel.selectTpcSet(itps);
  const WireSelector::WireInfoVector& wdat = sel.fillData();
  const WireSelector::WireInfoMap& wmap = sel.fillDataMap();
  const WireSelector::WireSummary& wsum = sel.fillWireSummary();
  const GeometryCore* pgeo = sel.geometry();
  string gname = pgeo == nullptr ? "NONE" : pgeo->DetectorName();
  if ( m_LogLevel ) {
    cout << myname << "Derived: " << endl;
    cout << myname << "           Geometry name: " << gname << endl;
    cout << myname << "       # selected planes: " << sel.planeIDs().size() << endl;
    cout << myname << "        # selected wires: " << wdat.size() << endl;
    cout << myname << "  # selected wire points: " << wsum.size() << endl;
    cout << myname << "     # selected channels: " << wmap.size() << endl;
  }
}

//**********************************************************************

void AdcDetectorPlotter::initialize() {
}

//**********************************************************************

DataMap AdcDetectorPlotter::viewMap(const AdcChannelDataMap& acds) const {
  const string myname = "AdcDetectorPlotter::view: ";
  DataMap ret;
  State& state = *getState();
  if ( m_LogLevel >= 2 ) cout << myname << "Begin call " << state.reportCount
                              << "/" << state.jobCount << "." << endl;
  if ( acds.size() == 0 ) {
    cout << myname << "WARNING: Channel map is empty. No data extracted." << endl;
    return ret.setStatus(1);
  }
  const AdcChannelData& acdFirst = acds.begin()->second;
  const AdcChannelData& acdLast = acds.rbegin()->second;
  string hname = "hdet";
  int npadx = 1200;
  int npady = 1000;
  double xmin = m_XMin;
  double xmax = m_XMax;
  string sttlx = "Drift coordinate [cm]";
  double xsign = 1.0;
  if ( xmax < xmin ) {
    xmin = m_XMax;
    xmax = m_XMin;
    sttlx = "-" + sttlx;
    xsign = -1.0;
  }
  string sttly = "Wire coordinate [cm]";
  if ( state.reportCount ) {
    if ( acdFirst.run() != state.run ||
         acdFirst.subRun() != state.subrun ||
         acdFirst.event() != state.event ) {
      cout << myname << "ERROR: Received unexpected event ID. Clearing data." << endl;
      cout << myname << "State: " << state.event << "-" << state.subrun << "-" << state.event << endl;
      cout << myname << " Data: " << acdFirst.event() << "-" << acdFirst.subRun() << "-" << acdFirst.event() << endl;
      state.reportCount = 0;
    }
  }
  if ( state.reportCount == 0 ) {
    if ( m_LogLevel >= 2 ) cout << myname << "  Starting new event." << endl;
    initializeState(state, acdFirst);
    string sttl = AdcChannelStringTool::build(m_adcStringBuilder, acdFirst, m_Title);
    string spttl = AdcChannelStringTool::build(m_adcStringBuilder, acdFirst, m_PlotTitle);
    state.ofname = AdcChannelStringTool::build(m_adcStringBuilder, acdFirst, m_FileName);
    // Create graph.
    state.ppad.reset(new TPadManipulator(npadx, npady));
    TGraph* pg = new TGraph;
    pg->SetTitle(sttl.c_str());
    pg->GetXaxis()->SetTitle(sttlx.c_str());
    pg->GetYaxis()->SetTitle(sttly.c_str());
    state.ppad->add(pg, "P");
    state.ppad->graph()->Expand(512);  // Allocate 1024 points.
    state.ppad->setRangeX(xmin, xmax);
    state.ppad->setRangeY(m_ZMin, m_ZMax);
    if ( m_ShowGrid ) state.ppad->setGrid();
    LineColors cols;
    const WireSelector::WireInfoVector& wins = state.sel.fillData();
    if ( m_ShowWires ) {
      TGraph* pgw = new TGraph;
      pgw->SetMarkerColor(cols.red());
      pgw->Expand(wins.size());
      for ( Index iwin=0; iwin<wins.size(); ++iwin ) {
        const WireSelector::WireInfo& win = wins[iwin];
        pgw->SetPoint(iwin, xsign*win.x, win.z);
      }
      state.ppad->add(pgw, "P");
    }
    if ( m_ShowCathode ) {
      TGraph* pgc = new TGraph;
      pgc->SetMarkerColor(cols.green());
      pgc->Expand(wins.size());
      for ( Index iwin=0; iwin<wins.size(); ++iwin ) {
        const WireSelector::WireInfo& win = wins[iwin];
        pgc->SetPoint(iwin, xsign*(win.x + win.driftMax), win.z);
      }
      state.ppad->add(pgc, "P");
    }
    // Add lower left label.
    if ( spttl.size() ) {
      state.pttl.reset(new TLatex(0.01, 0.015, spttl.c_str()));
      state.pttl->SetNDC();
      state.pttl->SetTextFont(42);
      state.pttl->SetTextSize(0.030);
      state.ppad->add(state.pttl.get());
    }
    ++state.jobCount;
  } else {
    TGraph* pgr = state.ppad->graph();
    if ( pgr == nullptr ) {
      cout << "ERROR: Graph not found." << endl;
      return ret;
    }
    if ( m_LogLevel >= 2 ) {
      cout << myname << "  Adding to existing event. Graph point count is " << pgr->GetN() << endl;
    }
  }
  ++state.reportCount;
  Tick maxtick = 0;
  for ( const AdcChannelDataMap::value_type& iacd : acds ) {
    if ( iacd.first == AdcChannelData::badChannel() ) {
      cout << myname << "WARNING: Channel map has invalid channels. No plot is created." << endl;
    }
    Tick ntick = iacd.second.samples.size();
    if ( ntick > maxtick ) maxtick = ntick;
  }
  AdcIndex chanFirst = acdFirst.channel();
  AdcIndex chanLast = acdLast.channel();
  AdcIndex nchan = chanLast + 1 - chanFirst;
  if ( m_LogLevel >= 2 ) cout << myname << "  Input channel count is " << nchan << endl;
  // Fill graph.
  for ( const AdcChannelDataMap::value_type& iacd : acds ) {
    if ( m_LogLevel >= 3 ) cout << myname << "    Filling with channel " << iacd.first << endl;
    const AdcChannelData& acd = iacd.second;
    if ( m_SkipBadChannels && m_pChannelStatusProvider != nullptr &&
         m_pChannelStatusProvider->IsBad(acd.channel()) ) {
      if ( m_LogLevel >= 3 ) cout << myname << "      Skipping bad channel " << acd.channel() << endl;
    } else {
      if ( m_LogLevel >= 3 ) cout << myname << "      Adding channel " << acd.channel() << endl;
      addChannel(acd, xsign);
    }
  }
  if ( state.ppad->graph()->GetN() == 0 ) {
    cout << myname << "Graph has no points. Adding one to avoid root exception." << endl;
    state.ppad->graph()->SetPoint(0, m_XMin, m_ZMin);
  }
  if ( m_LogLevel >= 2 ) cout << myname << "  Graph point count: " << state.ppad->graph()->GetN() << endl;
  return ret;
}

//**********************************************************************

DataMap AdcDetectorPlotter::endEvent(const DuneEventInfo& evi) const {
  const string myname = "AdcDetectorPlotter::endEvent: ";
  DataMap ret;
  State& state = *getState();
  if ( state.reportCount == 0 ) {
    cout << myname << "WARNING: No data recorded." << endl;
    return ret;
  }
  if ( evi.run != state.run ) {
    cout << myname << "ERROR: Received request for unexpected run: " << evi.run << endl;
    return ret.setStatus(1);
  }
  if ( evi.event != state.event ) {
    cout << myname << "ERROR: Received request for unexpected event: " << evi.event << endl;
    return ret.setStatus(2);
  }
  state.ppad->print(state.ofname);
  state.reportCount = 0;
  return ret;
}

//**********************************************************************

int AdcDetectorPlotter::addChannel(const AdcChannelData& acd, double xsign) const {
  const string myname = "AdcDetectorPlotter::addChannel: ";
  bool isRaw = m_DataType == 1;
  bool isPrep = m_DataType == 0;
  if ( !isRaw && !isPrep ) {
    cout << myname << "ERROR: Invalid data type: " << m_DataType << endl;
    return 2;
  }
  State& state = *getState();
  TGraph* pg = state.ppad->graph();
  if ( pg == nullptr ) {
    cout << myname << "ERROR: Graph is missing." << endl;
    return 1;
  }
  AdcChannel icha = acd.channel();
  auto rng = state.sel.dataMap().equal_range(icha);
  Index nsam = isRaw ? acd.raw.size() : acd.samples.size();
  if ( m_LogLevel >= 4 ) {
    cout << myname << "  Adding " << nsam << (isRaw ? " raw" : " processed")
         << " sample" << (nsam==1 ? "" : "s") << endl;
  } 
  for ( auto ient=rng.first; ient!=rng.second; ++ient) {
    const WireSelector::WireInfo& win = *(ient->second);
    float z = win.z;
    float x1 = win.x1();
    float x2 = win.x2();
    float driftVelocity = win.driftSign()*m_DriftSpeed;
    Index isam1 = 0;
    Index isam2 = nsam;
    if ( m_LastTick > m_FirstTick ) {
      if ( m_FirstTick > isam1 ) isam1 = m_FirstTick;
      if ( m_LastTick < isam2 ) isam2 = m_LastTick;
    }
    for ( Index isam=isam1; isam<isam2; ++isam ) {
      float sig = isRaw ? acd.raw[isam] - acd.pedestal : acd.samples[isam];
      if ( sig > m_SignalThreshold ) {
        float x = win.x + driftVelocity*(isam - m_Tick0);
        if ( ! m_ShowAllTicks ) {
          if ( x < x1 ) continue;
          if ( x > x2 ) continue;
        }
        Index ipt = pg->GetN();
        pg->SetPoint(ipt, xsign*x, z);
        if ( m_LogLevel >= 4 ) {
          ostringstream sout;
          sout.precision(2);
          sout << myname << "Added point " << ipt << " (" << x << ", " << z << ")";
          cout << sout.str() << endl;
        }
      } else if ( m_LogLevel >= 5 ) {
        cout << myname << "Skipped sample " << isam << " with signal " << sig << endl;
      }
    }
  }
/*
  AdcChannel chan = acd.channel();
  const AdcSignalVector& sams = acd.samples;
  const AdcCountVector& raw = acd.raw;
  AdcSignal ped = 0.0;
  bool isRawPed = false;
  if ( isRaw ) {
    ped = acd.pedestal;
    isRawPed = ped != AdcChannelData::badSignal;
  }
    unsigned int ibin = ph->GetBin(1, chan-chanFirst+1);
    for ( Tick isam=0; isam<sams.size(); ++isam, ++ibin ) {
      unsigned int isig = isam + m_FirstTick;
      if ( isPrep ) {
        if ( isig < sams.size() ) ph->SetBinContent(ibin, sams[isig]);
      } else if ( isRawPed ) {
        if ( isig < raw.size() ) ph->SetBinContent(ibin, raw[isig] - ped);
      } else {
        cout << myname << "Fill failed for bin " << ibin << endl;
      }
  }
*/
  return 0;
}

//**********************************************************************

DEFINE_ART_CLASS_TOOL(AdcDetectorPlotter)
