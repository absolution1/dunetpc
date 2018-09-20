// AdcDataPlotter_tool.cc

#include "AdcDataPlotter.h"
#include <iostream>
#include <sstream>
#include "dune/DuneCommon/TPadManipulator.h"
#include "dune/DuneCommon/StringManipulator.h"
#include "dune/DuneCommon/RootPalette.h"
#include "dune/ArtSupport/DuneToolManager.h"
#include "dune/DuneInterface/Tool/AdcChannelStringTool.h"
#include "dune/DuneInterface/Tool/IndexMapTool.h"
#include "dune/DuneInterface/Tool/IndexRangeTool.h"
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
  m_TickRange(ps.get<string>("TickRange")),
  m_ChannelRanges(ps.get<NameVector>("ChannelRanges")),
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
  // Fetch tick range.
  string descTickRange;
  if ( m_TickRange.size() ) {
    const string tnam = "tickRanges";
    const IndexRangeTool* ptool = ptm->getShared<IndexRangeTool>(tnam);
    if ( ptool == nullptr ) { 
      cout << myname << "WARNING: Tick range tool not found: " << tnam << endl;
    } else {
      m_tickRange = ptool->get(m_TickRange);
    }
    if ( ! m_tickRange.isValid() ) {
      cout << myname << "WARNING: Tick range not found: " << m_TickRange << endl;
    }
    descTickRange = m_TickRange + " " + m_tickRange.rangeString();
  }
  // Fetch channel ranges.
  const IndexRangeTool* pcrt = nullptr;
  for ( Name crn : m_ChannelRanges.size() ? m_ChannelRanges : NameVector(1, "") ) {
    if ( crn.size() == 0 || crn == "data" ) {
      m_crs.emplace_back("data", 0, 0, "All data");
    } else {
      if ( pcrt == nullptr ) {
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
  // Display configuration.
  if ( m_LogLevel ) {
    cout << myname << "Configuration: " << endl;
    cout << myname << "              LogLevel: " << m_LogLevel << endl;
    cout << myname << "              DataType: " << m_DataType << endl;
    cout << myname << "             TickRange: " << descTickRange << endl;
    cout << myname << "       ChannelRanges: [";
    bool first = true;
    for ( const IndexRange& ran : m_crs ) {
      if ( ! first ) cout << ", ";
      else first = false;
      cout << ran.name;
    }
    cout << "]" << endl;
    cout << myname << "       FembTickOffsets: [";
    first = true;
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
  // Find the tick range.
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
  unsigned long tick1 = 0;
  unsigned long tick2 = maxtick;
  if ( m_tickRange.isValid() ) {
    tick1 = m_tickRange.begin;
    tick2 = m_tickRange.end;
  }
  Tick ntick = tick2 - tick1;
  ret.setInt("ntick", ntick);
  if ( ntick <= 0 ) {
    cout << myname << "ERROR: Invalid tick range." << endl;
    return ret.setStatus(2);
  }
  // Loop over channel ranges.
  Index nhist = 0;
  for ( IndexRange ran : m_crs ) {
    if ( ran.name == "data" ) {
      ran.begin = acdChanFirst;
      ran.end = acdChanLast + 1;
    }
    if ( ! ran.isValid() ) {
      cout << myname << "ERROR: Skipping invalid channel range " << ran.name << endl;
      continue;
    }
    Index chanBegin = ran.first();
    Index chanEnd = ran.end;
    // Find the range of data channels in the plot range: [ichanDataBegin, ichanDataEnd)
    AdcChannelDataMap::const_iterator ichanDataBegin = acds.lower_bound(chanBegin);
    AdcChannelDataMap::const_iterator   ichanDataEnd = acds.upper_bound(ran.last());
    Index chanDataEnd   =   ichanDataEnd == acds.end() ?     chanEnd :   ichanDataEnd->first;
    Index chanDataBegin = ichanDataBegin == acds.end() ? chanDataEnd : ichanDataBegin->first;
    // Skip plot if no data channels are in range.
    if ( chanDataEnd <= chanDataBegin ) {
      if ( m_LogLevel >= 2 ) cout << myname << "Data has no entries in channel range " << ran.name << endl;
      if ( m_LogLevel >= 3 ) {
        cout << myname << "      chanBegin: " << chanBegin << endl;
        cout << myname << "        chanEnd: " << chanEnd << endl;
        cout << myname << "  chanDataBegin: " << chanDataBegin << endl;
        cout << myname << "    chanDataEnd: " << chanDataEnd << endl;
      }
      continue;
    }
    if ( m_LogLevel >= 2 ) cout << myname << "Creating histogram for channel range " << ran.name << endl;
    AdcIndex nchan = chanEnd - chanBegin;
    // Create title and file names.
    DataMap dm;
    dm.setInt("chan1", acdChanFirst);
    dm.setInt("chan2", acdChanLast-1);
    string   hname = nameReplace(    m_HistName, acdFirst, ran);
    string   htitl = nameReplace(   m_HistTitle, acdFirst, ran);
    string  ofname = nameReplace(m_PlotFileName, acdFirst, ran);
    string ofrname = nameReplace(m_RootFileName, acdFirst, ran);
    string szunits = "(ADC counts)";
    if ( ! isRaw ) {
      szunits = acdFirst.sampleUnit;
      if ( szunits.find(" ") != string::npos ) szunits = "(" + szunits + ")";
    }
    htitl += "; Tick; Channel; Signal [" + szunits + "/Tick/Channel]";
    // Create histogram.
    TH2* ph = new TH2F(hname.c_str(), htitl.c_str(), ntick, tick1, tick2, nchan, chanBegin, chanEnd);
    ph->SetDirectory(nullptr);
    ph->SetStats(0);
    if ( m_LogLevel >= 2 ) cout << myname << "Created histogram " << hname << endl;
    double zmax = m_MaxSignal;
    if ( zmax <= 0.0 ) zmax = 100.0;
    ph->GetZaxis()->SetRangeUser(-zmax, zmax);
    ph->SetContour(40);
    // Fill histogram.
    const bool doZero = false;
    for ( AdcChannel chan=chanDataBegin; chan<chanDataEnd; ++chan ) {
      unsigned int ibin = ph->GetBin(1, chan-chanBegin+1);
      AdcChannelDataMap::const_iterator iacd = acds.find(chan);
      if ( iacd == acds.end() ) {
        if ( doZero ) {
          if ( m_LogLevel >= 3 ) cout << myname << "Filling channel-tick histogram with zero for channel " << chan << endl;
          unsigned int ibin = ph->GetBin(1, chan-chanBegin+1);
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
        if ( icha > chanBegin && icha < chanEnd ) {
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
      cout << "for channels [" << chanBegin << " - "
                              << chanEnd << "): " << ofname << endl;
    }
    if ( ofrname.size() ) {
      TFile* pfile = TFile::Open(ofrname.c_str(), "UPDATE");
      ph->Write();
      if ( m_LogLevel > 1 ) cout << myname << "Wrote " << ph->GetName() << " to " << ofrname << endl;
      delete pfile;
    }
    oldPalette.setRootPalette();
    ++nhist;
  }
  ret.setInt("nhist", nhist);
  return ret;
}

//**********************************************************************

string AdcDataPlotter::
nameReplace(string name, const AdcChannelData& acd, const IndexRange& ran) const {
  StringManipulator sman(name);
  sman.replace("%CRNAME%", ran.name);
  sman.replace("%CRLABEL%", ran.label(0));
  sman.replace("%CRLABEL1%", ran.label(1));
  sman.replace("%CRLABEL2%", ran.label(2));
  const AdcChannelStringTool* pnbl = m_adcStringBuilder;
  if ( pnbl == nullptr ) return name;
  DataMap dm;
  dm.setInt("chan1", ran.first());
  dm.setInt("chan2", ran.last());
  return pnbl->build(acd, dm, name);
}

//**********************************************************************
