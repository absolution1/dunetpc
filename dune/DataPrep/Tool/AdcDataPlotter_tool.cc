// AdcDataPlotter_tool.cc

#include "AdcDataPlotter.h"
#include <iostream>
#include <sstream>
#include "dune/DuneInterface/Data/RunData.h"
#include "dune/DuneInterface/Tool/RunDataTool.h"
#include "dune/DuneCommon/Utility/TPadManipulator.h"
#include "dune/DuneCommon/Utility/StringManipulator.h"
#include "dune/DuneCommon/Utility/RootPalette.h"
#include "dune/DuneCommon/Utility/RootParFormula.h"
#include "dune/ArtSupport/DuneToolManager.h"
#include "dune/DuneInterface/Tool/AdcChannelStringTool.h"
#include "dune/DuneInterface/Tool/IndexMapTool.h"
#include "dune/DuneInterface/Tool/IndexRangeTool.h"
#include "larevt/CalibrationDBI/Interface/ChannelStatusService.h"
#include "larevt/CalibrationDBI/Interface/ChannelStatusProvider.h"
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
  m_DataView(ps.get<string>("DataView")),
  m_TickRange(ps.get<string>("TickRange")),
  m_TickRebin(ps.get<Index>("TickRebin")),
  m_ChannelRanges(ps.get<NameVector>("ChannelRanges")),
  m_ClockFactor(ps.get<float>("ClockFactor")),
  m_ClockOffset(ps.get<float>("ClockOffset")),
  m_FembTickOffsets(ps.get<IntVector>("FembTickOffsets")),
  m_MinSignal(nullptr),
  m_MaxSignal(new RootParFormula("MaxSignal", ps.get<Name>("MaxSignal"))),
  m_SkipBadChannels(ps.get<bool>("SkipBadChannels")),
  m_EmptyColor(ps.get<double>("EmptyColor")),
  m_ChannelLineModulus(ps.get<Index>("ChannelLineModulus")),
  m_ChannelLinePattern(ps.get<IndexVector>("ChannelLinePattern")),
  m_Palette(ps.get<int>("Palette")),
  m_HistName(ps.get<string>("HistName")),
  m_HistTitle(ps.get<string>("HistTitle")),
  m_PlotTitle(ps.get<string>("PlotTitle")),
  m_PlotSizeX(ps.get<Index>("PlotSizeX")),
  m_PlotSizeY(ps.get<Index>("PlotSizeY")),
  m_PlotFileName(ps.get<string>("PlotFileName")),
  m_RootFileName(ps.get<string>("RootFileName")),
  m_needRunData(false),
  m_pOnlineChannelMapTool(nullptr),
  m_pChannelStatusProvider(nullptr),
  m_prdtool(nullptr) {
  const string myname = "AdcDataPlotter::ctor: ";
  string stmp;
  if ( m_MaxSignal->npar() ) m_needRunData = true;
  ps.get_if_present("MinSignal", stmp);
  if ( stmp.size() ) {
    m_MinSignal = new RootParFormula("MinSignal", stmp);
    if ( m_MinSignal->npar() ) m_needRunData = true;
  }
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
      // Issue warning if the range does not have  label.
      if ( m_tickRange.label(0) == "" ) {
        cout << myname << "WARNING: Tick range not found: " << m_TickRange << endl;
      }
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
  // Fetch the channel status service.
  if ( m_SkipBadChannels ) {
    if ( m_LogLevel >= 1 ) cout << myname << "Fetching channel status service." << endl;
    m_pChannelStatusProvider = &art::ServiceHandle<lariov::ChannelStatusService>()->GetProvider();
    if ( m_pChannelStatusProvider == nullptr ) {
      cout << myname << "WARNING: Channel status provider not found." << endl;
    }
  }
  // Fetch the run data tool.
  if ( m_needRunData ) {
    string stnam = "runDataTool";
    m_prdtool = ptm->getShared<RunDataTool>(stnam);
    if ( m_prdtool == nullptr ) {
      cout << myname << "ERROR: RunDataTool " << stnam
           << " not found. Metric limits will not be evaluated." << endl;
    } else {
      cout << myname << "RunDataTool retrieved." << endl;
    }
  }
  // Display configuration.
  if ( m_LogLevel ) {
    cout << myname << "Configuration: " << endl;
    cout << myname << "              LogLevel: " << m_LogLevel << endl;
    cout << myname << "              DataType: " << m_DataType << endl;
    cout << myname << "              DataView: " << m_DataView << endl;
    cout << myname << "             TickRange: " << descTickRange << endl;
    cout << myname << "             TickRebin: " << m_TickRebin << endl;
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
    cout << myname << "           ClockFactor: " << m_ClockFactor << endl;
    cout << myname << "           ClockOffset: " << m_ClockOffset << endl;
    if ( m_FembTickOffsets.size() ) {
      cout << myname << "  OnlineChannelMapTool: " << m_OnlineChannelMapTool << " @ "
           << m_pOnlineChannelMapTool << endl;
    }
    if ( m_MinSignal != nullptr )
      cout << myname << "             MinSignal: " << m_MinSignal->formulaString() << endl;
    cout << myname << "             MaxSignal: " << m_MaxSignal->formulaString() << endl;
    cout << myname << "            EmptyColor: " << m_EmptyColor << endl;
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
    cout << myname << "           PlotTitle: " << m_PlotTitle << endl;
    cout << myname << "           PlotSizeX: " << m_PlotSizeX << endl;
    cout << myname << "           PlotSizeY: " << m_PlotSizeY << endl;
    cout << myname << "        PlotFileName: " << m_PlotFileName << endl;
    cout << myname << "        RootFileName: " << m_RootFileName << endl;
  }
}

//**********************************************************************

DataMap AdcDataPlotter::viewMap(const AdcChannelDataMap& acds) const {
  const string myname = "AdcDataPlotter::viewMap: ";
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
  bool isNsg = m_DataType == 3;
  // Find the tick range.
  Tick maxtick = 0;
  for ( const AdcChannelDataMap::value_type& iacd : acds ) {
    if ( iacd.first == AdcChannelData::badChannel() ) {
      cout << myname << "WARNING: Channel map has invalid channels. No plot is created." << endl;
    }
    Tick ntick = isRaw ? iacd.second.raw.size() : iacd.second.samples.size();
    if ( ntick > maxtick ) maxtick = ntick;
  }
  AdcIndex acdChanFirst = acdFirst.channel();
  AdcIndex acdChanLast = acdLast.channel();
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
  // Evaluate formulas.
  if ( m_prdtool != nullptr ) {
    Index run = acdFirst.run();
    RunData rdat = m_prdtool->runData(run);
    if ( ! rdat.isValid() ) {
      cout << myname << "WARNING: RunData not found for run " << run << "." << endl;
    } else {
      rdat.setFormulaPars(*m_MaxSignal);
      if ( m_MinSignal != nullptr ) rdat.setFormulaPars(*m_MinSignal);
    }
  }
  double zmax = m_MaxSignal->eval();
  double zmin = m_MinSignal == nullptr ? -zmax : m_MinSignal->eval();
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
    string   ptitl = nameReplace(   m_PlotTitle, acdFirst, ran);
    string  ofname = nameReplace(m_PlotFileName, acdFirst, ran);
    string ofrname = nameReplace(m_RootFileName, acdFirst, ran);
    string szunits = "(ADC counts)";
    if ( ! isRaw ) {
      szunits = acdFirst.sampleUnit;
      if ( szunits.find(" ") != string::npos ) szunits = "(" + szunits + ")";
    }
    string zunit = szunits;
    if ( zunit.find(" ") != string::npos ) zunit = "(" + zunit + ")";
    bool zhasTick = zunit.find("Tick") != string::npos || zunit.find("tick") != string::npos;
    if ( ! zhasTick ) zunit += "/tick";
    zunit += "/channel";
    htitl += "; Tick; Channel; Signal [" + zunit + "]";
    // Set flag indicating we want to show empty bins with the color m_EmptyColor.
    // We initialize all bins below zmin and fill with zmin where the value would be lower.
    // We do not attempt this this where rebinning is done.
    //bool colorEmptyBins = m_EmptyColor >= 0 && m_TickRebin <= 1;
    bool colorEmptyBins = m_TickRebin <= 1;
    // Create histogram.
    TH2* ph = new TH2F(hname.c_str(), htitl.c_str(), ntick, tick1, tick2, nchan, chanBegin, chanEnd);
    ph->SetDirectory(nullptr);
    ph->SetStats(0);
    if ( m_LogLevel >= 2 ) cout << myname << "Created histogram " << hname << endl;
    ph->GetZaxis()->SetRangeUser(zmin, zmax);
    ph->SetContour(200);
    double zempty = colorEmptyBins ? zmin - 1000.0 : 0.0;
    // Initialize bins to zmin.
    for ( Index icha=1; icha<=nchan; ++icha ) {
      Index ibin0 = (ntick+2)*icha + 1;
      for ( Index ibin=ibin0; ibin<ibin0+ntick; ++ibin ) ph->SetBinContent(ibin, zempty);
    }
    // Fill histogram.
    for ( AdcChannel chan=chanDataBegin; chan<chanDataEnd; ++chan ) {
      AdcChannelDataMap::const_iterator iacd = acds.find(chan);
      if ( iacd == acds.end() ) continue;
      const AdcChannelData& acdtop = iacd->second;
      if ( m_SkipBadChannels && m_pChannelStatusProvider != nullptr &&
           m_pChannelStatusProvider->IsBad(acdtop.channel()) ) {
        if ( m_LogLevel >= 3 ) cout << myname << "Skipping bad channel " << acdtop.channel() << endl;
        continue;
      }
      Index ibiny = chan-chanBegin + 1;
      Index nent = acdtop.viewSize(m_DataView);
      if ( m_LogLevel >= 2 && nent == 0 ) {
        cout << myname << "WARNING: Unable to find data view name " << m_DataView << endl;
      }
      for ( Index ient=0; ient<nent; ++ient ) {
        if ( m_LogLevel >= 3 ) {
          cout << myname << "Processing data view " << ient << "/" << nent << endl;
        }
        const AdcChannelData* pacd = acdtop.viewEntry(m_DataView, ient);
        if ( pacd == nullptr ) {
          cout << myname << "ERROR: Skipping null view entry " << m_DataView << "[" << ient
               << "]." << endl;
          continue;
        }
        if ( pacd->channel() != chan ) {
          cout << myname << "ERROR: Skipping view entry " << m_DataView << "[" << ient
               << "] with the wrong the wrong channel: "
               << pacd->channel() << " != " << chan <<"." << endl;
          continue;
        }
        const AdcSignalVector& sams = pacd->samples;
        const AdcFilterVector& keep = pacd->signal;
        const AdcCountVector& raw = pacd->raw;
        AdcSignal ped = 0.0;
        bool isRawPed = false;
        if ( isRaw ) {
          ped = pacd->pedestal;
          isRawPed = ped != AdcChannelData::badSignal();
        }
        Tick nsam = isRaw ? raw.size() : sams.size();
        AdcInt dsam = pacd->tick0;
        if ( m_ClockFactor > 0.0 ) {
          double dclk = 0.0;
          if ( pacd->channelClock >= pacd->triggerClock() ) dclk = pacd->channelClock - pacd->triggerClock();
          else dclk = -double(pacd->triggerClock() - pacd->channelClock);
          dsam += m_ClockFactor*(dclk + m_ClockOffset);
        }
        if ( m_FembTickOffsets.size() ) {
          if ( m_pOnlineChannelMapTool == nullptr ) {
            cout << myname << "  FEMB tick offsets provided without online channel mapping tool." << endl;
            break;
          }
          Index chanOn = m_pOnlineChannelMapTool->get(chan);
          Index ifmb = chanOn/128;
          if ( ifmb < m_FembTickOffsets.size() ) dsam += m_FembTickOffsets[ifmb];
          
        }
        if ( m_LogLevel >= 3 ) {
          cout << myname << "Filling channel-tick histogram with " << nsam << " samples and offset "
               << dsam << " for channel " << chan << endl;
        }
        Index ibin = ph->GetBin(1, ibiny);
        for ( Tick itck=tick1; itck<tick2; ++itck, ++ibin ) {
          AdcInt iisam = itck;
          iisam -= dsam;
          if ( iisam >= 0 ) {
            AdcIndex isam = iisam;
            float sig = 0.0;
            if ( (isSig || isNsg) && isam >= pacd->signal.size() ) {
              if ( m_LogLevel >= 3 ) {
                cout << myname << "  Signal array not filled for sample " << isam
                     << " and above--stopping fill." << endl;
              }
              break;
            }
            if ( isPrep ) {
              if ( isam < sams.size() ) sig = sams[isam];
              else break;
            } else if ( isRawPed ) {
              if ( isam < raw.size() ) sig = raw[isam] - ped;
              else break;
            } else if ( isRaw && isam >= raw.size() ) {
              // Ped missing is ok if we have no raw
              break;
            } else if ( isSig ) {
              if ( isam < sams.size() && isam < keep.size() ) {
                if ( keep[isam] ) sig = sams[isam];
              } else {
                break;
              }
            } else if ( isNsg ) {
              if ( isam < sams.size() ) {
                if ( isam >= keep.size() || !keep[isam] ) sig = sams[isam];
              } else {
                break;
              }
            } else {
              cout << myname << "Fill failed for bin " << ibin << endl;
            }
            if ( colorEmptyBins && sig < zmin ) sig = zmin;
            ph->SetBinContent(ibin, sig);
          }
        }
      }
    }
    // Rebin.
    if ( m_TickRebin > 1 ) {
      TH2* ph0 = ph;
      ph = ph0->RebinX(m_TickRebin, ph->GetName());
      ph->SetDirectory(nullptr);
      delete ph0;
      ph->Scale(1.0/m_TickRebin);
      ph->GetZaxis()->SetRangeUser(-zmax, zmax);
    }
/*
    // Save the original color map.
    RootPalette oldPalette;
    RootPalette::set(m_Palette);
*/
    const RootPalette* ppal = RootPalette::find(m_Palette);
    if ( ppal == nullptr ) {
      cout << myname << "ERROR: Unable to find palette " << m_Palette << endl;
      return ret.setStatus(3);
    }
    TPadManipulator man;
    man.setPalette(m_Palette);
    if ( m_PlotSizeX && m_PlotSizeY ) man.setCanvasSize(m_PlotSizeX, m_PlotSizeY);
    man.add(ph, "colz");
    man.setRangeZ(zmin, zmax);
    man.addAxis();
    // Root uses the frame color for underflows.
    // If we are coloring empty bins, we use that color.
    // If not, we use the color for the first (lowest) bin.
    int frameColor = colorEmptyBins ? m_EmptyColor : ppal->colorVector()[0];
    man.setFrameFillColor(frameColor);   // Otherwise Root uses white for underflows
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
    TLatex* pptl = nullptr;
    if ( ptitl.size() ) {
      pptl = new TLatex(0.01, 0.015, ptitl.c_str());
      pptl->SetNDC();
      pptl->SetTextFont(42);
      pptl->SetTextSize(0.030);
      man.add(pptl);
    }
    man.print(ofname);
    delete pptl;
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
/*
    oldPalette.setRootPalette();
*/
    ++nhist;
  }
  ret.setInt("nhist", nhist);
  return ret;
}

//**********************************************************************

string AdcDataPlotter::
nameReplace(string name, const AdcChannelData& acd, const IndexRange& ran) const {
  StringManipulator sman(name, false);
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

DEFINE_ART_CLASS_TOOL(AdcDataPlotter)
