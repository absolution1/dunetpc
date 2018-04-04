// AdcDetectorPlotter_tool.cc

#include "AdcDetectorPlotter.h"
#include <iostream>
#include <sstream>
#include "dune/DuneCommon/TPadManipulator.h"
#include "dune/DuneCommon/StringManipulator.h"
#include "dune/DuneCommon/RootPalette.h"
#include "larcore/Geometry/Geometry.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "TColor.h"
#include "TStyle.h"
#include "TDirectory.h"
#include "TFile.h"

using std::string;
using std::cout;
using std::cin;
using std::endl;

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
  state.run    = acd.run;
  state.subrun = acd.subRun;
  state.event  = acd.event;
  state.ppad.reset(new TPadManipulator);
/*
  // Create title and file names.
  string htitl = m_Title;
  string ofname = m_FileName;
  vector<string*> strs = {&hname, &htitl, &ofname, &ofrname};
  for ( string* pstr : strs ) {
    string& str = *pstr;
    StringManipulator sman(str);
    sman.replace("%PAT%", fpat);
    if ( acdFirst.run != AdcChannelData::badIndex ) sman.replace("%RUN%", acdFirst.run);
    else sman.replace("%RUN%", "RunNotFound");
    if ( acdFirst.subRun != AdcChannelData::badIndex ) sman.replace("%SUBRUN%", acdFirst.subRun);
    else sman.replace("%SUBRUN%", "SubRunNotFound");
    if ( acdFirst.event != AdcChannelData::badIndex ) sman.replace("%EVENT%", acdFirst.event);
    else sman.replace("%EVENT%", "EventNotFound");
    sman.replace("%CHAN1%", chanFirst);
    sman.replace("%CHAN2%", chanLast);
  }
  string szunits =  isRaw ? "(ADC counts)" : acdFirst.sampleUnit;
  htitl += "; Tick; Channel; Signal [" + szunits + "/Tick/Channel]";
  // Create histogram.
  TH2* ph = new TH2F(hname.c_str(), htitl.c_str(), ntick, tick1, tick2, nchan, chanFirst, chanLast+1);
  ph->SetDirectory(nullptr);
  ph->SetStats(0);
  double zmax = m_MaxSignal;
  if ( zmax <= 0.0 ) zmax = 100.0;
  ph->GetZaxis()->SetRangeUser(-zmax, zmax);
  ph->SetContour(40);
  man.add(ph, "colz");
  man.addAxis();
  man.print(ofname);
  if ( 0 ) {
    string line;
    cout << myname;
    cout.flush();
    std::getline(cin, line);
  }
  if ( m_LogLevel > 1 ) {
    cout << myname << "Created plot ";
    if ( label.size() ) cout << label << " ";
    cout << "for channels " << acds.begin()->first << " - "
                            << acds.rbegin()->first
         << ": " << ofname << endl;
  }
*/
}

}  // end unnamed namespace

//**********************************************************************
// Class methods.
//**********************************************************************

AdcDetectorPlotter::AdcDetectorPlotter(fhicl::ParameterSet const& ps)
: m_LogLevel(ps.get<int>("LogLevel")), 
  m_WireAngle(ps.get<float>("WireAngle")),
  m_DataType(ps.get<int>("DataType")),
  m_XMin(ps.get<float>("XMin")),
  m_XMax(ps.get<float>("XMax")),
  m_ZMin(ps.get<float>("ZMin")),
  m_ZMax(ps.get<float>("ZMax")),
  m_SignalThreshold(ps.get<float>("SignalThreshold")),
  m_Title(ps.get<string>("Title")),
  m_FileName(ps.get<string>("FileName")),
  m_state(new State) {
  const string myname = "AdcDetectorPlotter::ctor: ";
  if ( m_LogLevel ) {
    cout << myname << "Configuration: " << endl;
    cout << myname << "         LogLevel: " << m_LogLevel << endl;
    cout << myname << "        WireAngle: " << m_WireAngle << endl;
    cout << myname << "         DataType: " << m_DataType << endl;
    cout << myname << "             XMin: " << m_XMin << endl;
    cout << myname << "             XMax: " << m_XMax << endl;
    cout << myname << "             ZMin: " << m_ZMin << endl;
    cout << myname << "             ZMax: " << m_ZMax << endl;
    cout << myname << "  SignalThreshold: " << m_SignalThreshold << endl;
    cout << myname << "            Title: " << m_Title << endl;
    cout << myname << "         FileName: " << m_FileName << endl;
  }
  WireSelector& sel = getState()->sel;
  sel.selectWireAngle(m_WireAngle);
  const WireSelector::WireInfoVector& wdat = sel.fillData();
  const WireSelector::WireInfoMap& wmap = sel.fillDataMap();
  const WireSelector::WireSummary& wsum = sel.fillWireSummary();
  const GeometryCore* pgeo = sel.geometry();
  string gname = pgeo == nullptr ? "NONE" : pgeo->DetectorName();
  if ( m_LogLevel ) {
    cout << myname << "Derived: " << endl;
    cout << myname << "  Geometry name :" << gname << endl;
    cout << myname << "       # selected planes :" << sel.planeIDs().size() << endl;
    cout << myname << "        # selected wires :" << wdat.size() << endl;
    cout << myname << "  # selected wire points :" << wsum.size() << endl;
    cout << myname << "     # selected channels :" << wmap.size() << endl;
  }
}

//**********************************************************************

void AdcDetectorPlotter::initialize() {
}

//**********************************************************************

int AdcDetectorPlotter::view(const AdcChannelDataMap& acds, string label, string fpat) const {
  const string myname = "AdcDetectorPlotter::view: ";
  State& state = *getState();
  if ( m_LogLevel >= 2 ) cout << myname << "Begin call " << state.reportCount
                              << "/" << state.jobCount << "." << endl;
  if ( acds.size() == 0 ) {
    cout << myname << "WARNING: Channel map is empty. No data extracted." << endl;
    return 1;
  }
  const AdcChannelData& acdFirst = acds.begin()->second;
  const AdcChannelData& acdLast = acds.rbegin()->second;
  bool isRaw = m_DataType == 1;
  bool isPrep = m_DataType == 0;
  if ( !isRaw && !isPrep ) {
    cout << myname << "ERROR: Invalid data type: " << m_DataType << endl;
    return 2;
  }
  string hname = "hdet";
  if ( state.jobCount == 0 ||
       acdFirst.run != state.run || acdFirst.subRun != state.subrun || acdFirst.event != state.event ) {
    initializeState(state, acdFirst);
    string htitl = m_Title;
    state.ofname = m_FileName;
    vector<string*> strs = {&htitl, &state.ofname};
    for ( string* pstr : strs ) {
      string& str = *pstr;
      StringManipulator sman(str);
      sman.replace("%PAT%", fpat);
      if ( acdFirst.run != AdcChannelData::badIndex ) sman.replace("%RUN%", acdFirst.run);
      else sman.replace("%RUN%", "RunNotFound");
      if ( acdFirst.subRun != AdcChannelData::badIndex ) sman.replace("%SUBRUN%", acdFirst.subRun);
      else sman.replace("%SUBRUN%", "SubRunNotFound");
      if ( acdFirst.event != AdcChannelData::badIndex ) sman.replace("%EVENT%", acdFirst.event);
      else sman.replace("%EVENT%", "EventNotFound");
    }
    //string szunits =  isRaw ? "(ADC counts)" : acdFirst.sampleUnit;
    //htitl += "; Tick; Channel; Signal [" + szunits + "/Tick/Channel]";
    htitl += "; Tick; Channel";
    // Create histogram.
    TH2* ph = new TH2F(hname.c_str(), htitl.c_str(), 10, m_ZMin, m_ZMax, 10, m_XMin, m_XMax);
    ph->SetDirectory(nullptr);
    ph->SetStats(0);
    state.ppad->add(ph);
    delete ph;
  }
  ++state.jobCount;
  ++state.reportCount;
  Tick maxtick = 0;
  for ( const AdcChannelDataMap::value_type& iacd : acds ) {
    if ( iacd.first == AdcChannelData::badChannel ) {
      cout << myname << "WARNING: Channel map has invalid channels. No plot is created." << endl;
    }
    Tick ntick = iacd.second.samples.size();
    if ( ntick > maxtick ) maxtick = ntick;
  }
  AdcIndex chanFirst = acdFirst.channel;
  AdcIndex chanLast = acdLast.channel;
  AdcIndex nchan = chanLast + 1 - chanFirst;
  if ( m_LogLevel >= 2 ) cout << myname << "Input channel count is " << nchan << endl;
/*
  // Fill histogram.
  for ( const AdcChannelDataMap::value_type& iacd : acds ) {
    AdcChannel chan = iacd.first;
    const AdcChannelData& acd = iacd.second;
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
  }
*/
  state.ppad->print(state.ofname);
  return 0;
}

