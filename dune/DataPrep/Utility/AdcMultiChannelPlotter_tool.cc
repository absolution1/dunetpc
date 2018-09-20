// AdcMultiChannelPlotter_t.cc

#include "AdcMultiChannelPlotter.h"
#include "dune/DuneInterface/Tool/AdcChannelStringTool.h"
#include "dune/DuneCommon/TPadManipulator.h"
#include "dune/DuneCommon/StringManipulator.h"
#include "dune/ArtSupport/DuneToolManager.h"
#include <iostream>

using std::cout;
using std::endl;

//**********************************************************************

AdcMultiChannelPlotter::AdcMultiChannelPlotter(fhicl::ParameterSet const& ps, Name prefix)
: m_PlotName(ps.get<Name>(prefix + "Name")),
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
}

//**********************************************************************

DataMap AdcMultiChannelPlotter::viewMap(const AdcChannelDataMap& acds) const {
  const Name myname = "AdcMultiChannelPlotter::viewMap: ";
  DataMap ret;
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
  // Loop over channels.
  Index nacd = acds.size();
  Index iacd = 0;
  Index ncha = 0;
  Index nplt = 0;
  for ( auto& acdPair : acds ) {
    const AdcChannelData& acd = acdPair.second;
    // If needed, create a new canvas and a name.
    if ( pmantop == nullptr ) {
      if ( getLogLevel() >= 3 ) cout << myname << "  Creating canvas." << endl;
      pmantop = new TPadManipulator;
      if ( m_PlotSizeX && m_PlotSizeY ) pmantop->setCanvasSize(m_PlotSizeX, m_PlotSizeY);
      if ( npad > 1 ) pmantop->split(npady, npady);
      plotName = AdcChannelStringTool::build(m_adcStringBuilder, acd, getPlotName());
    }
    // View this channel.
    Index ipad = npad == 0 ? 0 : iacd % npad;
    viewMapChannel(acd, ret, *pmantop->man(ipad));
    // Handle the end of a plot file.
    ++iacd;
    bool lastpad = (npad == 0) || (++ipad == npad) || (iacd == nacd);
    if ( lastpad && pmantop != nullptr ) {
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
