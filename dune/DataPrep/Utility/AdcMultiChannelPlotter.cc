// AdcMultiChannelPlotter_t.cc

#include "AdcMultiChannelPlotter.h"
#include "dune/DuneInterface/Tool/AdcChannelStringTool.h"
#include "dune/DuneInterface/Tool/IndexRangeTool.h"
#include "dune/DuneInterface/Tool/IndexRangeGroupTool.h"
#include "dune/DuneCommon/Utility/TPadManipulator.h"
#include "dune/DuneCommon/Utility/StringManipulator.h"
#include "dune/ArtSupport/DuneToolManager.h"
#include <iostream>
#include <sstream>

using std::cout;
using std::endl;
using std::ostringstream;

//**********************************************************************

AdcMultiChannelPlotter::AdcMultiChannelPlotter(fhicl::ParameterSet const& ps, Name prefix)
: m_LogLevel(ps.get<Index>("LogLevel")),
  m_PlotChannelRanges(ps.get<NameVector>(prefix + "ChannelRanges")),
  m_PlotChannelGroups(ps.get<NameVector>(prefix + "ChannelGroups")),
  m_PlotOverlayGroups(ps.get<Index>(prefix + "OverlayGroups")),
  m_PlotDataView(ps.get<Name>(prefix + "DataView")),
  m_PlotName(ps.get<Name>(prefix + "Name")),
  m_PlotSummaryNames(ps.get<NameVector>(prefix + "SummaryNames")),
  m_PlotSizeX(ps.get<Index>(prefix + "SizeX")),
  m_PlotSizeY(ps.get<Index>(prefix + "SizeY")),
  m_PlotSplitX(ps.get<Index>(prefix + "SplitX")),
  m_PlotSplitY(ps.get<Index>(prefix + "SplitY")) {
  const Name myname = "AdcMultiChannelPlotter::ctor: ";
  if ( m_LogLevel >= 1 ) {
    cout << myname << "           LogLevel: " << m_LogLevel << endl;
    bool first = true;
    cout << myname << "   " << prefix + "ChannelRanges: [";
    for ( Name crn : m_PlotChannelRanges ) {
      if ( first ) first = false;
      else cout << ", ";
      cout << crn;
    }
    cout << "]" << endl;
    cout << myname << "   " << prefix + "ChannelGroups: [";
    first = true;
    for ( Name cgn : m_PlotChannelGroups ) {
      if ( first ) first = false;
      else cout << ", ";
      cout << cgn;
    }
    cout << "]" << endl;
    cout << myname << "  " << prefix << "OverlayGroups: " << m_PlotOverlayGroups << endl;
    cout << myname << "       " << prefix << "DataView: " << m_PlotDataView << endl;
    cout << myname << "           " << prefix << "Name: " << m_PlotName << endl;
    cout << myname << "   " << prefix << "SummaryNames: [";
    first = true;
    for ( Name name : m_PlotSummaryNames ) {
      if ( first ) first = false;
      else cout << ", ";
      cout << name;
    }
    cout << "]" << endl;
    cout << myname << "          " << prefix << "SizeX: " << m_PlotSizeX << endl;
    cout << myname << "          " << prefix << "SizeY: " << m_PlotSizeY << endl;
    cout << myname << "         " << prefix << "SplitX: " << m_PlotSplitX << endl;
    cout << myname << "         " << prefix << "SplitY: " << m_PlotSplitY << endl;
  }
  Name stringBuilder = "adcStringBuilder";
  DuneToolManager* ptm = DuneToolManager::instance();
  if ( ptm == nullptr ) {
    cout << myname << "ERROR: Tool manager not found." << endl;
  } else {
    // Fetch the string builder.
    m_adcStringBuilder = ptm->getShared<AdcChannelStringTool>(stringBuilder);
    if ( m_adcStringBuilder == nullptr ) {
      cout << myname << "WARNING: AdcChannelStringTool not found: " << stringBuilder << endl;
    }
    // Build pad descriptions from groups.
    if ( haveChannelGroups() && ptm != nullptr ) {
      const IndexRangeGroupTool* pcgt = ptm->getShared<IndexRangeGroupTool>("channelGroups");
      if ( pcgt == nullptr ) {
        cout << myname << "ERROR: IndexRangeGroupTool not found: channelGroups" << endl;
      } else {
        for ( Name cgn : m_PlotChannelGroups ) {
          IndexRangeGroup grp = pcgt->get(cgn);
          if ( ! grp.isValid() ) {
            cout << myname << "WARNING: Skipping invalid group " << cgn << endl;
            continue;
          }
          m_cgmap[cgn] = grp;
          Pad* ppad = nullptr;
          for ( const IndexRange& ran : grp.ranges ) {
            Name crn = ran.name;
            if ( ! ran.isValid() ) {
              cout << myname << "WARNING: Skipping invalid range " << crn << " in group " << cgn << endl;
              continue;
            }
            if ( ppad == nullptr ) {
              m_pads.emplace_back(cgn);
              ppad = &m_pads.back();
            }
            if ( ppad->crmap.find(crn) != ppad->crmap.end() ) {
              cout << myname << "WARNING: Ignoring duplicate channel range group " << cgn << endl;
              continue;
            }
            ppad->crnames.push_back(crn);
            ppad->crmap[crn] = ran;
            if ( ! overlayGroups() ) ppad = nullptr;
          }  // end loop over ranges in group
        }  // end loop over groups
      }
    }
    // Build pad descriptions from ranges.
    if ( haveChannelRanges() && ptm != nullptr ) {
      const IndexRangeTool* pcrt = ptm->getShared<IndexRangeTool>("channelRanges");
      if ( pcrt == nullptr ) {
        cout << myname << "ERROR: IndexRangeTool not found: channelRanges" << endl;
      } else {
        for ( Name crn : m_PlotChannelRanges ) {
          IndexRange ran = pcrt->get(crn);
          if ( ! ran.isValid() ) {
            cout << myname << "WARNING: Skipping invalid range " << crn << endl;
            continue;
          }
          m_pads.emplace_back(crn);    // Use range name as group name
          Pad& pad = m_pads.back();
          pad.crnames.push_back(crn);
          pad.crmap[crn] = ran;
        }  // end loop over ranges
      }
    }
  }
  if ( m_LogLevel >= 1 ) cout << myname << "Drawing pad count is " << m_pads.size() << endl;
}

//**********************************************************************

AdcMultiChannelPlotter::~AdcMultiChannelPlotter() {
  // viewSummary();   // Can't call this here b/c it uses virtual calls.
}

//**********************************************************************

DataMap AdcMultiChannelPlotter::viewMap(const AdcChannelDataMap& acds) const {
  const Name myname = "AdcMultiChannelPlotter::viewMap: ";
  DataMap ret;
  if ( acds.size() == 0 ) return ret;
  if ( ! getBaseState().hasRun() ) {
    getBaseState().setRun(acds.begin()->second.run());
  }
  getBaseState().event = acds.begin()->second.event();
  Index npadx = 0;
  Index npady = 0;
  Index npadOnPage = 0;
  if ( getPlotName().size() && m_PlotSplitX > 0 ) {
    npadx = m_PlotSplitX;
    npady = m_PlotSplitY ? m_PlotSplitY : m_PlotSplitX;
    npadOnPage = npadx*npady;
  }
  if ( getLogLevel() >= 2 ) {
    cout << myname << "Per page pad count is " << npadOnPage << " (" << npady << " x " << npadx << ")" << endl;
  }
  Name plotName;
  Name cgrn;
  TPadManipulator* pmantop = nullptr;
  bool doRanges = haveChannelRanges() || haveChannelGroups();
  PadVector localPads;
  // Build pad vector for single channels.
  if ( ! doRanges ) {
    for ( AdcChannelDataMap::value_type iacd : acds ) {
      const AdcChannelData& acd = iacd.second;
      Index icha = acd.channel();
      getBaseState().channels.insert(icha);
      ostringstream sscha;
      sscha << icha;
      Name scha = sscha.str();
      while ( scha.size() < 6 ) scha = "0" + scha;
      scha = "chan" + scha;
      localPads.emplace_back(scha);
      Pad& pad = localPads.back();
      pad.crnames.push_back(scha);
      pad.crmap[scha] = IndexRange(icha, icha+1);
    }
  }
  const PadVector& pads = doRanges ? m_pads : localPads;
  Index npad = pads.size();
  if ( getLogLevel() >= 4 ) cout << myname << "Pad count is " << npad << endl;
  Index ncha = 0;
  Index nplt = 0;
  Index ipadOnPage = 0;
  for ( Index ipad=0; ipad<npad; ++ipad ) {
    const Pad& pad = pads[ipad];
    Name cgn = pad.cgname;
    if ( getLogLevel() >= 4 ) {
      cout << myname << "  Processing pad " << ipad << " for channel group " << cgn << endl;
    }
    Index ncr = pad.crnames.size();
    for ( Index icr=0; icr<ncr; ++icr ) {
      Name crn = pad.crnames[icr];
      const IndexRange& ran = pad.crmap.at(crn);
      if ( ! ran.isValid() ) {
        cout << "ERROR: Range " << crn << " in pad group " << cgn << " is invalid. Aborting." << endl;
        abort();
      }
      if ( getLogLevel() >= 5 ) {
        cout << myname << "    Processing range " << ran << endl;
      }
      AdcChannelDataMap::const_iterator iacd1 = acds.lower_bound(ran.first());
      if ( iacd1 == acds.end() ) {
        if ( getLogLevel() >= 5 ) cout << myname << "   Skipping range above data." << endl;
        continue;
      }
      AdcChannelDataMap::const_iterator iacd2 = acds.upper_bound(ran.last());
      if ( iacd2 == iacd1 ) {
        if ( getLogLevel() >= 5 ) cout << myname << "   Skipping range below data." << endl;
        continue;
      }
      AcdVector acdvec;
      Index ncha = 0;
      for ( AdcChannelDataMap::const_iterator iacd=iacd1; iacd!=iacd2; ++iacd ) {
        const AdcChannelData& acd = iacd->second;
        Index icha = iacd->first;
        if ( ! acd.hasView(getDataView()) ) {
          cout << myname << "WARNING: View " << getDataView() << " is missing for channel "
               << icha << "." << endl;
        }
        Index ndve = acd.viewSize(getDataView());
        for ( Index idve=0; idve<ndve; ++idve ) {
           acdvec.push_back(acd.viewEntry(getDataView(), idve));
        }
        ++ncha;
      }
      if ( getLogLevel() >= 4 ) {
        Index nacd = acdvec.size();
        cout << myname << "    Processing channel range " << crn << " with " << nacd
             << " entr" << (nacd == 1 ? "y" : "ies" ) << " in " << ncha << " channel"
             << (ncha == 1 ? "" : "s") << "."  << endl;
      }
      // If needed, create a new canvas.
      AdcChannelData acdTmp;   // Temporary channel data to hold event and (someday?) channel info.
      acdTmp.setEventInfo(acds.begin()->second.getEventInfoPtr());
      if ( pmantop == nullptr ) {
        pmantop = new TPadManipulator;
        if ( m_PlotSizeX && m_PlotSizeY ) pmantop->setCanvasSize(m_PlotSizeX, m_PlotSizeY);
        if ( npadOnPage > 1 ) pmantop->split(npadx, npady);
      }
      // Assign group/range name using the first plot on the canvas.
      if ( cgrn.size() == 0 ) {
        cgrn = overlayGroups() ? cgn : crn;
      }
      // Assign plot name using the first plot with data on the canvas.
      // Thepreviously-assigned group/range name is used.
      if ( plotName.size() == 0 && acdvec.size() ) {
        const AdcChannelData& acdFirst = *acdvec[0];
        plotName = AdcChannelStringTool::build(m_adcStringBuilder, acdFirst, getPlotName());
        StringManipulator sman(plotName, false);
        sman.replace("%CGNAME%", cgn);
        sman.replace("%CRNAME%", crn);
        sman.replace("%CGRNAME%", cgrn);
        sman.replace("%VIEW%", getDataView());
        plotName = sman.str();
        if ( getLogLevel() >= 4 ) {
          if ( plotName.size() ) cout << myname << "    Created pad  with plot name " << plotName << endl;
          else cout << myname << "    Created pad with no plot name." << endl;
        }
        if ( getLogLevel() >= 3 ) {
          cout << myname << "    Run " << acdFirst.run() << ", event " << acdFirst.event()
               << " group/range " << cgrn << " ==> plot name " << plotName << endl;
        }
      }
      // View this channel range.
      TPadManipulator* ppadPage = pmantop->man(ipadOnPage);
      if ( ppadPage == nullptr ) {
        cout << myname << "ERROR: View pad " << ipadOnPage << " not found for pad " << ipad
             << ", channel range " << crn << "." << endl;
      } else {
        viewMapChannels(crn, acdvec, *ppadPage, ncr, icr);
        ncha += acds.size();
      }
    }
    // Handle the end of a plot file.
    ++ipadOnPage;
    bool lastpad = (npadOnPage == 0) || (ipadOnPage == npadOnPage) || (ipad+1 == npad);
    if ( lastpad && pmantop != nullptr ) {
      if ( plotName.size() ) {
        if ( getLogLevel() >= 3 ) cout << myname << "  Printing canvas to " << plotName << endl;
        pmantop->print(plotName);
      } else {
        if ( getLogLevel() >= 3 ) cout << myname << "  Not printing canvas for channel group/range "
                                       << cgrn << "." << endl;
      }
      delete pmantop;
      pmantop = nullptr;
      plotName = "";
      cgrn = "";
      ipadOnPage = 0;
      ++nplt;
    }
  }
  delete pmantop;
  ret.setInt("multiChannelNChannel", ncha);
  ret.setInt("multiChannelNPlot", nplt);
  return ret;
}

//**********************************************************************

void AdcMultiChannelPlotter::viewSummary(Index ilev) const {
  const Name myname = "AdcMultiChannelPlotter::viewSummary: ";
  if ( getPlotSummaryName(ilev).size() == 0 ) {
    if ( getLogLevel() >= 3 ) cout << myname << "Summmary plots not named for level " << ilev << endl;
    return;
  }
  if ( getLogLevel() >= 3 ) cout << myname << "Creating summmary plots for level " << ilev << endl;
  Index npadx = 0;
  Index npady = 0;
  Index npadOnPage = 1;
  if ( m_PlotSplitX > 0 ) {
    npadx = m_PlotSplitX;
    npady = m_PlotSplitY ? m_PlotSplitY : m_PlotSplitX;
    npadOnPage = npadx*npady;
  }
  if ( getLogLevel() >= 2 ) {
    cout << myname << "Pad count/page is " << npadOnPage << " (" << npady << " x " << npadx << ")" << endl;
  }
  TPadManipulator* pmantop = nullptr;
  PadVector localPads;
  // Build pad vector for single channels.
  bool doRanges = haveChannelRanges() || haveChannelGroups();
  if ( ! doRanges ) {
    for ( Index icha : getBaseState().channels ) {
      ostringstream sscha;
      sscha << icha;
      Name scha = sscha.str();
      while ( scha.size() < 6 ) scha = "0" + scha;
      scha = "chan" + scha;
      localPads.emplace_back(scha);
      Pad& pad = localPads.back();
      pad.crnames.push_back(scha);
      pad.crmap[scha] = IndexRange(icha, icha+1);
    }
  }
  const PadVector& pads = doRanges ? m_pads : localPads;
  Index npad = pads.size();
  Name plotName;
  Name cgrn;
  AdcChannelData acdPrint;
  using EventInfo = AdcChannelData::EventInfo;
  acdPrint.setEventInfo(new EventInfo(getBaseState().run(), getBaseState().event));
  if ( getLogLevel() >= 3 ) cout << myname << "Pad count is " << npad << endl;
  for ( Index ipad =0; ipad<npad; ++ipad ) {
    const Pad& pad = pads[ipad];
    Name cgn = pad.cgname;
    Index ipadOnPage = npad == 0 ? 0 : ipad % npadOnPage;
    if ( getLogLevel() >= 4 ) cout << myname << "  Channel group " << cgn
                                   << " channel range count is " << pad.crnames.size() << endl;
    Index ncrn = pad.crnames.size();
    // # CRNs included in plot.
    for ( Index icrn=0; icrn<ncrn; ++icrn ) {
      Name crn = pad.crnames[icrn];
      if ( getLogLevel() >= 4 ) {
        cout << myname << "    Processing group/range " << cgn << "/" << crn << endl;
      }
      // Assign group/range name for the canvas.
      if ( cgrn.size() == 0 ) {
        cgrn = overlayGroups() ? cgn : crn;
      }
      // If needed, create a new canvas and a name.
      acdPrint.setChannelInfo(pad.crmap.at(crn).begin);
      if ( pmantop == nullptr ) {
        if ( getLogLevel() >= 3 ) cout << myname << "    Creating canvas for run " << acdPrint.run()
                                       << ", event " << acdPrint.event() << "." << endl;
        pmantop = new TPadManipulator;
        if ( m_PlotSizeX && m_PlotSizeY ) pmantop->setCanvasSize(m_PlotSizeX, m_PlotSizeY);
        if ( npadOnPage > 1 ) pmantop->split(npadx, npady);
        DataMap dmPrint;
        //dmPrint.setInt("CHAN1", ran.first());
        //dmPrint.setInt("CHAN2", ran.last());
        plotName = AdcChannelStringTool::build(m_adcStringBuilder, acdPrint, dmPrint, getPlotSummaryName(ilev));
        StringManipulator sman(plotName, false);
        sman.replace("%CRNAME%", crn);
        sman.replace("%CGNAME%", cgn);
        sman.replace("%CGRNAME%", cgrn);
        sman.replace("%VIEW%", getDataView());
      }
      // View this channel range.
      int rstat = viewMapSummary(ilev, cgn, crn, *pmantop->man(ipadOnPage), ncrn, icrn);
      if ( rstat >= 1 ) cout << myname << "WARNING: viewMapSummary returned error code " << rstat << endl;
      if ( getLogLevel() >= 5 ) cout << myname << "    Pad " << ipadOnPage << " extra object count: "
                                     << pmantop->man(ipadOnPage)->objects().size() << endl;
    }
    // Handle the end of a plot file.
    bool lastpad = (npadOnPage == 0) || (ipadOnPage+1 == npadOnPage) || (ipad+1 == npad);
    if ( lastpad && pmantop != nullptr ) {
      if ( getLogLevel() >= 3 ) cout << myname << "  Printing summary canvas to " << plotName << endl;
      pmantop->print(plotName);
      delete pmantop;
      pmantop = nullptr;
      plotName = "";
      cgrn = "";
    }
  }
}

//**********************************************************************

const IndexRangeGroup& AdcMultiChannelPlotter::getChannelGroup(Name cgn) const {
  const Name myname = "AdcMultiChannelPlotter::getChannelGroup: ";
  ChannelGroupMap::const_iterator icgr = m_cgmap.find(cgn);
  if ( icgr == m_cgmap.end() ) {
    // Not an error -- we may index by range or channel name instead of group name.
    //cout << myname << "ERROR: Group not found: " << cgn << endl;
    static const IndexRangeGroup badcgr;
    return badcgr;
  }
  return icgr->second;
}

//**********************************************************************
