// AdcDataPlotter_tool.cc

#include "AdcDataPlotter.h"
#include <iostream>
#include <sstream>
#include "dune/DuneCommon/TPadManipulator.h"
#include "dune/DuneCommon/RootPalette.h"
#include "dune/ArtSupport/DuneToolManager.h"
#include "dune/DuneInterface/Tool/AdcChannelStringTool.h"
#include "dune/DuneInterface/Tool/IndexMapTool.h"
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
using std::vector;

//**********************************************************************
// Class methods.
//**********************************************************************

AdcDataPlotter::AdcDataPlotter(fhicl::ParameterSet const& ps)
: m_LogLevel(ps.get<int>("LogLevel")), 
  m_DataType(ps.get<int>("DataType")),
  m_FirstTick(ps.get<unsigned long>("FirstTick")),
  m_LastTick(ps.get<unsigned long>("LastTick")),
  m_FirstChannel(ps.get<unsigned int>("FirstChannel")),
  m_LastChannel(ps.get<unsigned int>("LastChannel")),
  m_FembTickOffsets(ps.get<IntVector>("FembTickOffsets")),
  m_MaxSignal(ps.get<double>("MaxSignal")),
  m_ChannelLineModulus(ps.get<Index>("ChannelLineModulus")),
  m_ChannelLinePattern(ps.get<IndexVector>("ChannelLinePattern")),
  m_Palette(ps.get<int>("Palette")),
  m_HistName(ps.get<string>("HistName")),
  m_HistTitle(ps.get<string>("HistTitle")),
  m_PlotSizeX(ps.get<Index>("PlotSizeX")),
  m_PlotSizeY(ps.get<Index>("PlotSizeY")),
  m_PlotFileName(ps.get<string>("PlotFileName")),
  m_RootFileName(ps.get<string>("RootFileName")),
  m_pOnlineChannelMapTool(nullptr) {
  const string myname = "AdcDataPlotter::ctor: ";
  DuneToolManager* ptm = DuneToolManager::instance();
  string stringBuilder = "adcStringBuilder";
  m_adcStringBuilder = ptm->getShared<AdcChannelStringTool>(stringBuilder);
  if ( m_adcStringBuilder == nullptr ) {
    cout << myname << "WARNING: AdcChannelStringTool not found: " << stringBuilder << endl;
  }
  if ( m_FembTickOffsets.size() ) {
    m_OnlineChannelMapTool = ps.get<string>("OnlineChannelMapTool");
    m_pOnlineChannelMapTool = ptm->getShared<const IndexMapTool>(m_OnlineChannelMapTool);
  }
  if ( m_LogLevel ) {
    cout << myname << "Configuration: " << endl;
    cout << myname << "              LogLevel: " << m_LogLevel << endl;
    cout << myname << "              DataType: " << m_DataType << endl;
    cout << myname << "             FirstTick: " << m_FirstTick << endl;
    cout << myname << "              LastTick: " << m_LastTick << endl;
    cout << myname << "          FirstChannel: " << m_FirstChannel << endl;
    cout << myname << "           LastChannel: " << m_LastChannel << endl;
    cout << myname << "       FembTickOffsets: [";
    bool first = true;
    for ( int ioff : m_FembTickOffsets ) {
      if ( first ) first = false;
      else cout << ", ";
      cout << ioff;
    }
    cout << "]" << endl;
    if ( m_FembTickOffsets.size() ) {
      cout << myname << "  OnlineChannelMapTool: " << m_OnlineChannelMapTool << " @ "
           << m_pOnlineChannelMapTool << endl;
    }
    cout << myname << "             MaxSignal: " << m_MaxSignal << endl;
    cout << myname << "    ChannelLineModulus: " << m_ChannelLineModulus << endl;
    cout << myname << "    ChannelLinePattern: {";
    first = true;
    for ( Index icha : m_ChannelLinePattern ) {
      if ( ! first ) cout << ", ";
      first = false;
      cout << icha;
    }
    cout << "}" << endl;
    cout << myname << "             Palette: " << m_Palette << endl;
    cout << myname << "            HistName: " << m_HistName << endl;
    cout << myname << "           HistTitle: " << m_HistTitle << endl;
    cout << myname << "           PlotSizeX: " << m_PlotSizeX << endl;
    cout << myname << "           PlotSizeY: " << m_PlotSizeY << endl;
    cout << myname << "        PlotFileName: " << m_PlotFileName << endl;
    cout << myname << "        RootFileName: " << m_RootFileName << endl;
  }
}

//**********************************************************************

DataMap AdcDataPlotter::viewMap(const AdcChannelDataMap& acds) const {
  const string myname = "AdcDataPlotter::view: ";
  DataMap ret;
  if ( acds.size() == 0 ) {
    cout << myname << "WARNING: Channel map is empty. No plot is created." << endl;
    return ret.setStatus(1);
  }
  const AdcChannelData& acdFirst = acds.begin()->second;
  const AdcChannelData& acdLast = acds.rbegin()->second;
  bool isPrep = m_DataType == 0;
  bool isRaw = m_DataType == 1;
  bool isSig = m_DataType == 2;
  Tick maxtick = 0;
  for ( const AdcChannelDataMap::value_type& iacd : acds ) {
    if ( iacd.first == AdcChannelData::badChannel ) {
      cout << myname << "WARNING: Channel map has invalid channels. No plot is created." << endl;
    }
    Tick ntick = isRaw ? iacd.second.raw.size() : iacd.second.samples.size();
    if ( ntick > maxtick ) maxtick = ntick;
  }
  AdcIndex acdChanFirst = acdFirst.channel;
  AdcIndex acdChanLast = acdLast.channel;
  AdcIndex chanFirst = acdChanFirst;
  AdcIndex chanLast = acdChanLast;
  // If the parameters specify a channel range, we use it.
  // No action if the map does not have channels in this range.
  if ( m_LastChannel > m_FirstChannel ) {
    chanFirst = m_FirstChannel;
    chanLast = m_LastChannel - 1;
    if ( acdChanFirst > chanLast || acdChanLast < chanFirst ) return ret;
  }
  if ( chanLast < chanFirst ) {
    if ( m_LogLevel >= 3 ) cout << myname << "No channels in view range for data range ("
                                << acdChanFirst << ", " << acdChanLast << ")" << endl;
    if ( m_LogLevel >= 2 ) cout << myname << "Skipping plot for " << acds.size() << " channels." << endl;
    return ret;
  }
  if ( m_LogLevel >= 2 ) cout << myname << "Creating plot for " << acds.size() << " channels." << endl;
  unsigned long tick1 = m_FirstTick;
  unsigned long tick2 = m_LastTick;
  if ( tick2 <= tick1 ) {
    tick1 = 0;
    tick2 = maxtick;
  }
  if ( tick2 <= tick1 ) {
    cout << myname << "WARNING: Invalid tick range: (" << tick1 << ", " << tick2 << ")." << endl;
    cout << myname << "           Configured range: (" << m_FirstTick << ", " << m_LastTick << ")." << endl;
    cout << myname << "Data size: " << maxtick << " ticks" << endl;
    return ret.setStatus(2);
  }
  Tick ntick = tick2 - tick1;
  AdcIndex nchan = chanLast + 1 - chanFirst;
  // Create title and file names.
  DataMap dm;
  dm.setInt("chan1", acdChanFirst);
  dm.setInt("chan2", acdChanLast);
  string hname =   AdcChannelStringTool::build(m_adcStringBuilder, acdFirst, dm, m_HistName);
  string htitl =   AdcChannelStringTool::build(m_adcStringBuilder, acdFirst, dm, m_HistTitle);
  string ofname =  AdcChannelStringTool::build(m_adcStringBuilder, acdFirst, dm, m_PlotFileName);
  string ofrname = AdcChannelStringTool::build(m_adcStringBuilder, acdFirst, dm, m_RootFileName);
  string szunits = "(ADC counts)";
  if ( ! isRaw ) {
    szunits = acdFirst.sampleUnit;
    if ( szunits.find(" ") != string::npos ) szunits = "(" + szunits + ")";
  }
  htitl += "; Tick; Channel; Signal [" + szunits + "/Tick/Channel]";
  // Create histogram.
  ret.setInt("ntick", ntick);
  ret.setInt("nchan", nchan);
  TH2* ph = new TH2F(hname.c_str(), htitl.c_str(), ntick, tick1, tick2, nchan, chanFirst, chanLast+1);
  ph->SetDirectory(nullptr);
  ph->SetStats(0);
  double zmax = m_MaxSignal;
  if ( zmax <= 0.0 ) zmax = 100.0;
  ph->GetZaxis()->SetRangeUser(-zmax, zmax);
  ph->SetContour(40);
  // Fill histogram.
  const bool doZero = false;
  for ( AdcChannel chan=chanFirst; chan<chanLast; ++chan ) {
    unsigned int ibin = ph->GetBin(1, chan-chanFirst+1);
    AdcChannelDataMap::const_iterator iacd = acds.find(chan);
    if ( iacd == acds.end() ) {
      if ( doZero ) {
        if ( m_LogLevel >= 3 ) cout << myname << "Filling channel-tick histogram with zero for channel " << chan << endl;
        unsigned int ibin = ph->GetBin(1, chan-chanFirst+1);
        for ( Tick isam=tick1; isam<tick2; ++isam, ++ibin ) ph->SetBinContent(ibin, 0.0);
      } else {
        if ( m_LogLevel >= 3 ) cout << myname << "Not filling channel-tick histogram for channel " << chan << endl;
      }
      continue;
    }
    const AdcChannelData& acd = iacd->second;
    const AdcSignalVector& sams = acd.samples;
    const AdcFilterVector& keep = acd.signal;
    const AdcCountVector& raw = acd.raw;
    AdcSignal ped = 0.0;
    bool isRawPed = false;
    if ( isRaw ) {
      ped = acd.pedestal;
      isRawPed = ped != AdcChannelData::badSignal;
    }
    Tick nsam = isRaw ? raw.size() : sams.size();
    if ( m_LogLevel >= 3 ) {
      cout << myname << "Filling channel-tick histogram with " << nsam << " samples for channel " << chan << endl;
    }
    int tickOffset = 0;
    if ( m_FembTickOffsets.size() ) {
      if ( m_pOnlineChannelMapTool == nullptr ) {
        cout << myname << "  FEMB tick offsets provided without online channel mapping tool." << endl;
        break;
      }
      Index chanOn = m_pOnlineChannelMapTool->get(chan);
      Index ifmb = chanOn/128;
      if ( ifmb < m_FembTickOffsets.size() ) tickOffset = m_FembTickOffsets[ifmb];
    }
    Index dsam = tickOffset < 0 ? -tickOffset : tickOffset;
    bool addOffset = tickOffset > 0;
    bool subtractOffset = tickOffset < 0;
    for ( Tick itck=tick1; itck<tick2; ++itck, ++ibin ) {
      bool haveSam = true;
      if ( subtractOffset ) haveSam = itck >= dsam;
      float sig = 0.0;
      if ( haveSam ) {
        Index isam = itck;
        if ( addOffset ) isam += dsam;
        if ( subtractOffset ) isam -= dsam;
        if ( isSig && isam >= acd.signal.size() ) {
          if ( m_LogLevel >= 3 ) {
            cout << myname << "  Signal array not filled for sample " << isam << " and above--stopping fill." << endl;
          }
          break;
        }
        if ( isPrep ) {
          if ( isam < sams.size() ) sig = sams[isam];
        } else if ( isRawPed ) {
          if ( isam < raw.size() ) sig = raw[isam] - ped;
        } else if ( isSig ) {
          if ( isam < sams.size() && isam < keep.size() && keep[isam] ) sig = sams[isam];
        } else {
          cout << myname << "Fill failed for bin " << ibin << endl;
        }
      }
      ph->SetBinContent(ibin, sig);
    }
  }
  // Save the original color map.
  RootPalette oldPalette;
  RootPalette::set(m_Palette);
  const RootPalette* ppal = RootPalette::find(m_Palette);
  if ( ppal == nullptr ) {
    cout << myname << "ERROR: Unable to find palette " << m_Palette << endl;
    return ret.setStatus(3);
  }
  TPadManipulator man;
  if ( m_PlotSizeX && m_PlotSizeY ) man.setCanvasSize(m_PlotSizeX, m_PlotSizeY);
  man.add(ph, "colz");
  man.addAxis();
  man.setFrameFillColor(ppal->colorVector()[0]);
  if ( m_ChannelLineModulus ) {
    for ( Index icha : m_ChannelLinePattern ) {
      man.addHorizontalModLines(m_ChannelLineModulus, icha);
    }
  } else {
    for ( Index icha : m_ChannelLinePattern ) {
      if ( icha > chanFirst && icha < chanLast ) {
        man.addHorizontalLine(icha, 1.0, 3);
      }
    }
  }
  man.print(ofname);
  if ( 0 ) {
    string line;
    cout << myname;
    cout.flush();
    std::getline(cin, line);
  }
  if ( m_LogLevel > 1 ) {
    cout << myname << "Created plot ";
    cout << "for channels " << acds.begin()->first << " - "
                            << acds.rbegin()->first
         << ": " << ofname << endl;
  }
  if ( ofrname.size() ) {
    TFile* pfile = TFile::Open(ofrname.c_str(), "UPDATE");
    ph->Write();
    if ( m_LogLevel > 1 ) cout << myname << "Wrote " << ph->GetName() << " to " << ofrname << endl;
    delete pfile;
  }
  oldPalette.setRootPalette();
  return ret;
}

//**********************************************************************
