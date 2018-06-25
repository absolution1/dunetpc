// AdcTickModViewer_tool.cc

#include "AdcTickModViewer.h"
#include "dune/DataPrep/Utility/StickyCodeMetrics.h"
#include "dune/DuneCommon/TPadManipulator.h"
#include "dune/DuneCommon/StringManipulator.h"
#include "dune/DuneInterface/Tool/AdcChannelStringTool.h"
#include "dune/DuneInterface/Tool/TimeOffsetTool.h"
#include "dune/ArtSupport/DuneToolManager.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include "TDirectory.h"
#include "TFile.h"
#include "TH1F.h"
#include "TF1.h"
#include "TROOT.h"
#include "TError.h"

using std::string;
using std::cout;
using std::endl;
using std::ofstream;
using std::ostream;
using std::ostringstream;
using std::setw;
using std::fixed;
using std::setprecision;
using std::vector;

using Index = AdcTickModViewer::Index;
using Name = AdcTickModViewer::Name;

//**********************************************************************
// Class methods.
//**********************************************************************

AdcTickModViewer::AdcTickModViewer(fhicl::ParameterSet const& ps)
: m_LogLevel(ps.get<int>("LogLevel")),
  m_TickModPeriod(ps.get<Index>("TickModPeriod")),
  m_FitRmsMin(ps.get<float>("FitRmsMin")),
  m_FitRmsMax(ps.get<float>("FitRmsMax")),
  m_HistName(ps.get<string>("HistName")),
  m_HistTitle(ps.get<string>("HistTitle")),
  m_HistChannelCount(ps.get<Index>("HistChannelCount")),
  m_PlotFileName(ps.get<string>("PlotFileName")),
  m_RootFileName(ps.get<string>("RootFileName")),
  m_PlotChannels(ps.get<IndexVector>("PlotChannels")),
  m_PlotSizeX(ps.get<Index>("PlotSizeX")),
  m_PlotSizeY(ps.get<Index>("PlotSizeY")),
  m_PlotShowFit(ps.get<Index>("PlotShowFit")),
  m_PlotSplitX(ps.get<Index>("PlotSplitX")),
  m_PlotSplitY(ps.get<Index>("PlotSplitY")),
  m_state(new State) 
{
  const string myname = "AdcTickModViewer::ctor: ";
  DuneToolManager* ptm = DuneToolManager::instance();
  string stringBuilder = "adcStringBuilder";
  m_adcStringBuilder = ptm->getShared<AdcChannelStringTool>(stringBuilder);
  if ( m_adcStringBuilder == nullptr ) {
    cout << myname << "WARNING: AdcChannelStringTool not found: " << stringBuilder << endl;
  }
  string tname = "tickOffsetFinder";
  m_tickOffsetTool = ptm->getShared<TimeOffsetTool>(tname);
  if ( m_tickOffsetTool == nullptr ) {
    cout << myname << "WARNING: TimeOffsetTool not found: " << tname << endl;
  }
  if ( m_LogLevel >= 1 ) {
    cout << myname << "Configuration parameters:" << endl;
    cout << myname << "          LogLevel: " << m_LogLevel << endl;
    cout << myname << "     TickModPeriod: " << m_TickModPeriod << endl;
    cout << myname << "          HistName: " << m_HistName << endl;
    cout << myname << "         HistTitle: " << m_HistTitle << endl;
    cout << myname << "  HistChannelCount: " << m_HistChannelCount << endl;
    cout << myname << "      PlotFileName: " << m_PlotFileName << endl;
    cout << myname << "      RootFileName: " << m_RootFileName << endl;
    cout << myname << "      PlotChannels: [";
    bool first = true;
    for ( Index icha : m_PlotChannels ) {
      if ( !first ) cout << ", ";
      else first = false;
      cout << icha;
    }
    cout << "]" << endl;
    cout << myname << "         PlotSizeX: " << m_PlotSizeX << endl;
    cout << myname << "         PlotSizeY: " << m_PlotSizeY << endl;
    cout << myname << "       PlotShowFit: " << m_PlotShowFit << endl;
    cout << myname << "        PlotSplitX: " << m_PlotSplitX << endl;
    cout << myname << "        PlotSplitY: " << m_PlotSplitY << endl;
  }
  if ( m_LogLevel >=4 ) {
    cout << myname << "INFO: Checking state." << endl;
    cout << myname << "INFO:   Hist map size: " << state().ChannelTickModHists.size() << endl;
  }
}

//**********************************************************************

DataMap AdcTickModViewer::view(const AdcChannelData& acd) const {
  const string myname = "AdcTickModViewer::view: ";
  DataMap res;
  Index icha = acd.channel;
  if ( m_LogLevel >= 3 ) cout << myname << "Processing channel " << icha << endl;
  HistVector& tmhs = state().ChannelTickModHists[icha];
  Index ntkm = m_TickModPeriod;
  if ( tmhs.size() == 0 ) tmhs.resize(ntkm, nullptr);
  Index itkm0 = 0;
  for ( Index itkm=0; itkm<ntkm; ++itkm ) {
    fillChannelTickMod(acd, itkm0, itkm);
  }
  // Process histograms with the sticky code utility.
  HistVector tkmProcs;
  Index chmod = 10;
  for ( Index itkm=0; itkm<ntkm; ++itkm ) {
    TH1* phi = tmhs[itkm].get();
    StickyCodeMetrics scm(phi->GetName(), phi->GetTitle(), m_HistChannelCount, chmod);
    if ( scm.evaluate(phi) ) {
      cout << myname << "Sticky code evaluation failed for channel " << icha
           << " tickmod " << itkm << endl;
      tkmProcs.emplace_back(phi);
    } else {
      tkmProcs.push_back(scm.getSharedHist());
    }
  }
  // Draw the histograms.
  Index nplot = 0;
  if ( m_PlotFileName.size() &&
       ( m_PlotChannels.size() == 0 ||
         find(m_PlotChannels.begin(), m_PlotChannels.end(), icha) != m_PlotChannels.end() )
     ) {
    Index npad = 0;
    Index npadx = 0;
    Index npady = 0;
    if ( m_PlotFileName.size() && m_PlotSplitX > 0 ) {
      npadx = m_PlotSplitX;
      npady = m_PlotSplitY ? m_PlotSplitY : m_PlotSplitX;
      npad = npadx*npady;
    }
    TPadManipulator* pmantop = nullptr;
    Name plotFileName;
    Index itkm = 0;
    Index ipad = 0;
    for ( HistPtr ph : tkmProcs ) {
      if ( pmantop == nullptr ) {
        if ( m_LogLevel >= 3 ) cout << myname << "  Creating canvas." << endl;
        pmantop = new TPadManipulator;
        if ( m_PlotSizeX && m_PlotSizeY ) pmantop->setCanvasSize(m_PlotSizeX, m_PlotSizeY);
        if ( npad > 1 ) pmantop->split(npady, npady);
        plotFileName = nameReplace(m_PlotFileName, acd, itkm);
      }
      TPadManipulator* pman = pmantop->man(ipad);
      pman->add(ph.get(), "hist", false);
      if ( m_PlotShowFit > 1 ) pman->addHistFun(1);
      if ( m_PlotShowFit ) pman->addHistFun(0);
      pman->addVerticalModLines(64);
      pman->showUnderflow();
      pman->showOverflow();
      ++itkm;
      if ( ++ipad == npad || itkm == ntkm ) {
        pmantop->print(plotFileName);
        ++nplot;
        ipad = 0;
        delete pmantop;
        pmantop = nullptr;
      }
    }
  }
  res.setHistVector("tmHists", tkmProcs);
  res.setHistVector("tmWideHists", tmhs);   // Passing out hist sthat will be updated!
  res.setInt("tmCount", ntkm);
  res.setInt("tmPlotCount", nplot);
  return res;
}

//**********************************************************************

Name AdcTickModViewer::nameReplace(Name nameIn, const AdcChannelData& acd, Index itkm) const {
  const string myname = "AdcTickModViewer::nameReplace: ";
  DataMap dm;
  string nameOut = m_adcStringBuilder->build(acd, dm, nameIn);
  string stkm = std::to_string(itkm);
  string::size_type stkm0Len = std::to_string(m_TickModPeriod-1).size();
  string stkm0 = stkm;
  while ( stkm0.size() < stkm0Len ) stkm0 = "0" + stkm0;
  StringManipulator sman(nameOut);
  sman.replace("%TICKMOD%", stkm);
  sman.replace("%0TICKMOD%", stkm0);
  return nameOut;
}

//**********************************************************************

int
AdcTickModViewer::fillChannelTickMod(const AdcChannelData& acd, Index itkm0, Index itkm) const {
  const string myname = "AdcTickModViewer::fillChannelTickMod: ";
  DataMap res;
  Index nsam = acd.raw.size();
  if ( nsam == 0 ) {
    if ( m_LogLevel >= 2 ) cout << myname << "WARNING: Raw data is empty." << endl;
    return 1;
  }
  HistPtr& ph = state().ChannelTickModHists[acd.channel][itkm];
  if ( ph == nullptr ) {
    string hname = nameReplace(m_HistName, acd, itkm);
    string htitl = nameReplace(m_HistTitle, acd, itkm);
    htitl += "; ADC count; # samples";
    unsigned int nadc = 4096;
    if ( m_LogLevel >= 2 ) cout << myname << "Creating histogram " << hname
                                << " for channel " << acd.channel
                                << " tickmod " << itkm << endl;
    ph.reset(new TH1F(hname.c_str(), htitl.c_str(), nadc, 0, nadc));
    ph->SetDirectory(0);
  }
  Index period = m_TickModPeriod;
  if ( m_LogLevel >= 2 ) cout << myname << "Filling hist " << ph->GetName()
                              << " with channel " << acd.channel
                              << " tickmod " << itkm << endl;
  for ( Index isam = (itkm0 + itkm) % period; isam<nsam; isam += period ) ph->Fill(acd.raw[isam]);
  return 0;
}

//**********************************************************************
