// AdcMultiChannelPlotter_t.cc

#include "AdcMultiChannelPlotter.h"
#include "dune/DuneInterface/Tool/AdcChannelStringTool.h"
#include "dune/DuneInterface/Tool/IndexRangeTool.h"
#include "dune/DuneCommon/TPadManipulator.h"
#include "dune/DuneCommon/StringManipulator.h"
#include "dune/ArtSupport/DuneToolManager.h"
#include <iostream>
#include <sstream>

using std::cout;
using std::endl;
using std::ostringstream;

//**********************************************************************

AdcMultiChannelPlotter::AdcMultiChannelPlotter(fhicl::ParameterSet const& ps, Name prefix)
: m_PlotChannelRanges(ps.get<NameVector>(prefix + "ChannelRanges")),
  m_PlotName(ps.get<Name>(prefix + "Name")),
  m_PlotSizeX(ps.get<Index>(prefix + "SizeX")),
  m_PlotSizeY(ps.get<Index>(prefix + "SizeY")),
  m_PlotSplitX(ps.get<Index>(prefix + "SplitX")),
  m_PlotSplitY(ps.get<Index>(prefix + "SplitY")) {
  const Name myname = "AdcMultiChannelPlotter::ctor: ";
  Name stringBuilder = "adcStringBuilder";
  DuneToolManager* ptm = DuneToolManager::instance();
  m_adcStringBuilder = ptm->getShared<AdcChannelStringTool>(stringBuilder);
  if ( m_adcStringBuilder == nullptr ) {
    cout << myname << "WARNING: AdcChannelStringTool not found: " << stringBuilder << endl;
  }
  if ( m_PlotChannelRanges.size() ) {
    const IndexRangeTool* pcrt = ptm->getShared<IndexRangeTool>("channelRanges");
    if ( pcrt == nullptr ) {
      cout << myname << "ERROR: IndexRangeTool not found: channelRanges" << endl;
    } else {
      for ( Name crn : m_PlotChannelRanges ) {
        IndexRange ran = pcrt->get(crn);
        if ( ran.isValid() ) {
          if ( m_crmap.find(crn) == m_crmap.end() ) {
            m_crmap[crn] = ran;
          } else {
            cout << myname << "WARNING: Ignoring duplicate channel range " << crn << endl;
          }
        } else {
          cout << myname << "WARNING: Ignoring invalid include channel range " << crn << endl;
        }
      }
    }
  }
}

//**********************************************************************

DataMap AdcMultiChannelPlotter::viewMap(const AdcChannelDataMap& acds) const {
  const Name myname = "AdcMultiChannelPlotter::viewMap: ";
  DataMap ret;
  if ( acds.size() == 0 ) return ret;
  Index npadx = 0;
  Index npady = 0;
  Index npad = 0;
  if ( getPlotName().size() && m_PlotSplitX > 0 ) {
    npadx = m_PlotSplitX;
    npady = m_PlotSplitY ? m_PlotSplitY : m_PlotSplitX;
    npad = npadx*npady;
  }
  if ( getLogLevel() >= 2 ) {
    cout << myname << "Pad count is " << npad << " (" << npady << " x " << npadx << ")" << endl;
  }
  Name plotName;
  TPadManipulator* pmantop = nullptr;
  bool doRanges = haveChannelRanges();
  // Build the vector of channel data objects for each range.
  // Also build ordered list of retained CR names.
  using AcdMap = std::map<Name, AcdVector>;
  AcdMap acdmap;
  NameVector crNames;
  if ( doRanges ) {
    // Make a plot for each range that overlaps the data.
    if ( getLogLevel() >= 4 ) cout << myname << "Creating plots for channel ranges." << endl;
    for ( ChannelRangeMap::value_type namedRange : m_crmap ) {
      Name crn = namedRange.first;
      const IndexRange& ran = namedRange.second;
      if ( getLogLevel() >= 5 ) cout << myname << "  Range " << crn << ": " << ran << endl;
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
      for ( AdcChannelDataMap::const_iterator iacd=iacd1; iacd!=iacd2; ++iacd ) {
        const AdcChannelData* pacd = &iacd->second;
        if ( acdmap[crn].size() == 0 ) crNames.push_back(crn);
        acdmap[crn].push_back(pacd);
      }
    }
  } else {
    // Make a plot for each channel.
    if ( getLogLevel() >= 4 ) cout << myname << "Creating plots for channels." << endl;
    for ( const AdcChannelDataMap::value_type& chanAcd : acds ) {
      Index icha = chanAcd.first;
      const AdcChannelData* pacd = &chanAcd.second;
      ostringstream sscha;
      sscha << icha;
      Name scha = sscha.str();
      while ( scha.size() < 6 ) scha = "0" + scha;
      scha = "chan" + scha;
      acdmap[scha].push_back(pacd);
      crNames.push_back(scha);
    }
  }
  // Loop over channels or channel ranges.
  Index ncan = crNames.size();
  if ( getLogLevel() >= 4 ) cout << myname << "Channel range count is " << ncan << endl;
  Index ican = 0;
  Index ncha = 0;
  Index nplt = 0;
  for ( Name crn : crNames ) {
    const AcdVector& crAcds = acdmap[crn];
    if ( getLogLevel() >= 4 ) {
      Index nacd = crAcds.size();
      cout << myname << "  Processing channel range " << crn << " with " << nacd
           << " channel" << (nacd == 1 ? "" : "s") << "."  << endl;
    }
    const AdcChannelData& acdFirst = *crAcds.front();
    // If needed, create a new canvas and a name.
    if ( pmantop == nullptr ) {
      if ( getLogLevel() >= 3 ) cout << myname << "  Creating canvas." << endl;
      pmantop = new TPadManipulator;
      if ( m_PlotSizeX && m_PlotSizeY ) pmantop->setCanvasSize(m_PlotSizeX, m_PlotSizeY);
      if ( npad > 1 ) pmantop->split(npadx, npady);
      plotName = AdcChannelStringTool::build(m_adcStringBuilder, acdFirst, getPlotName());
      StringManipulator sman(plotName);
      sman.replace("%CRNAME%", crn);
      plotName = sman.string();
    }
    // View this channel range.
    Index ipad = npad == 0 ? 0 : ican % npad;
    viewMapChannels(crn, crAcds, ret, *pmantop->man(ipad));
    // Handle the end of a plot file.
    ++ican;
    bool lastpad = (npad == 0) || (++ipad == npad) || (ican == ncan);
    if ( lastpad && pmantop != nullptr ) {
      if ( getLogLevel() >= 3 ) cout << myname << "  Printing canvas to " << plotName << endl;
      pmantop->print(plotName);
      delete pmantop;
      pmantop = nullptr;
      ++nplt;
    }
    ++ncha;
  }
  ret.setInt("multiChannelNChannel", ncha);
  ret.setInt("multiChannelNPlot", nplt);
  return ret;
}

//**********************************************************************
