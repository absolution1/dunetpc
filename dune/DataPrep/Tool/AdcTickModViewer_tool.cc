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
#include "TTree.h"
#include "TROOT.h"
#include "TError.h"
#include "TGraph.h"

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
// Private definitions.
//**********************************************************************

namespace {

void copyAcd(const AdcChannelData& acdin, AdcChannelData& acdout) {
  acdout.run         = acdin.run;
  acdout.subRun      = acdin.subRun;
  acdout.event       = acdin.event;
  acdout.channel     = acdin.channel;
  acdout.fembID      = acdin.fembID;
  acdout.fembChannel = acdin.fembChannel;
  acdout.pedestal    = acdin.pedestal;
}

// Class to hold the adc counts for a tickmod.
class TickModData {
public:
  AdcCountVector data;
  AdcCount adcMin = 0;
  AdcCount adcMax = 0;
  long adcSum = 0;
  bool dbg = false;
  void add(AdcCount adc) {
    if ( data.size() ) {
      if ( adc < adcMin ) adcMin = adc;
      if ( adc > adcMax ) adcMax = adc;
    } else {
      adcMin = adc;
      adcMax = adc;
    }
    data.push_back(adc);
    adcSum += adc;
    if ( dbg ) cout << "... " << adc << " (" << adcMin << ", " << adcMax << ")" << endl;
  }
  AdcCount min() const { return adcMin; }
  AdcCount max() const { return adcMax; }
  double mean() const {
    return double(adcSum)/double(data.size());
  }
};

}  // end unnamed namespace

//**********************************************************************
// Class methods.
//**********************************************************************

AdcTickModViewer::AdcTickModViewer(fhicl::ParameterSet const& ps)
: m_LogLevel(ps.get<int>("LogLevel")),
  m_TickModPeriod(ps.get<Index>("TickModPeriod")),
  m_TimeOffsetTool(ps.get<Name>("TimeOffsetTool")),
  m_FitSigmaMin(ps.get<float>("FitSigmaMin")),
  m_FitSigmaMax(ps.get<float>("FitSigmaMax")),
  m_HistName(ps.get<string>("HistName")),
  m_HistTitle(ps.get<string>("HistTitle")),
  m_HistChannelCount(ps.get<Index>("HistChannelCount")),
  m_AllPlotFileName(ps.get<string>("AllPlotFileName")),
  m_MinPlotFileName(ps.get<string>("MinPlotFileName")),
  m_MaxPlotFileName(ps.get<string>("MaxPlotFileName")),
  m_PhasePlotFileName(ps.get<string>("PhasePlotFileName")),
  m_PhaseVariable(ps.get<string>("PhaseVariable")),
  m_RootFileName(ps.get<string>("RootFileName")),
  m_TreeFileName(ps.get<string>("TreeFileName")),
  m_PlotChannels(ps.get<IndexVector>("PlotChannels")),
  m_PlotSizeX(ps.get<Index>("PlotSizeX")),
  m_PlotSizeY(ps.get<Index>("PlotSizeY")),
  m_PlotShowFit(ps.get<Index>("PlotShowFit")),
  m_PlotSplitX(ps.get<Index>("PlotSplitX")),
  m_PlotSplitY(ps.get<Index>("PlotSplitY")),
  m_PlotFrequency(ps.get<Index>("PlotFrequency")),
  m_PhaseGrouping(ps.get<Name>("PhaseGrouping")),
  m_PhasePlotSizeX(ps.get<Index>("PhasePlotSizeX")),
  m_PhasePlotSizeY(ps.get<Index>("PhasePlotSizeY")),
  m_PhasePlotSplitX(ps.get<Index>("PhasePlotSplitX")),
  m_PhasePlotSplitY(ps.get<Index>("PhasePlotSplitY")),
  m_tickOffsetTool(nullptr),
  m_varPhase(m_PhaseVariable == "phase"),
  m_varEvent(m_PhaseVariable == "event"),
  m_varTick0(m_PhaseVariable == "tick0"),
  m_varNchan(m_PhaseVariable == "nchan"),
  m_plotAll(m_AllPlotFileName.size()),
  m_plotMin(m_MinPlotFileName.size()),
  m_plotMax(m_MaxPlotFileName.size()),
  m_plotAny(m_plotAll || m_plotMin || m_plotMax),
  m_plotPhase(m_PhasePlotFileName.size()),
  m_makeTree(m_TreeFileName.size()),
  m_groupByChannel(m_PhaseGrouping == "channel"),
  m_groupByFemb(m_PhaseGrouping == "femb"),
  m_state(new State)
{
  const string myname = "AdcTickModViewer::ctor: ";
  DuneToolManager* ptm = DuneToolManager::instance();
  string stringBuilder = "adcStringBuilder";
  m_adcStringBuilder = ptm->getShared<AdcChannelStringTool>(stringBuilder);
  if ( m_adcStringBuilder == nullptr ) {
    cout << myname << "WARNING: AdcChannelStringTool not found: " << stringBuilder << endl;
  }
  string tname = m_TimeOffsetTool;
  if ( tname.size() ) {
    m_tickOffsetTool = ptm->getShared<TimeOffsetTool>(tname);
    if ( m_tickOffsetTool == nullptr ) {
      cout << myname << "WARNING: Requested TimeOffsetTool not found: " << tname << endl;
    }
  }
  if ( m_plotPhase ) {
    if ( ! m_groupByChannel && ! m_groupByFemb ) {
      cout << myname << "ERROR: Invalid phase grouping: " << m_PhaseGrouping << endl;
      m_plotPhase = false;
    }
    if ( !m_varPhase && !m_varEvent && !m_varTick0 && !m_varNchan ) {
      cout << myname << "ERROR: Invalid phase variable: " << m_PhaseVariable << endl;
      m_plotPhase = false;
    }
  }
  if ( m_LogLevel >= 1 ) {
    cout << myname << "Configuration parameters:" << endl;
    cout << myname << "          LogLevel: " << m_LogLevel << endl;
    cout << myname << "     TickModPeriod: " << m_TickModPeriod << endl;
    cout << myname << "    TimeOffsetTool: " << m_TimeOffsetTool << endl;
    cout << myname << "       FitSigmaMin: " << m_FitSigmaMin << endl;
    cout << myname << "       FitSigmaMax: " << m_FitSigmaMax << endl;
    cout << myname << "          HistName: " << m_HistName << endl;
    cout << myname << "         HistTitle: " << m_HistTitle << endl;
    cout << myname << "  HistChannelCount: " << m_HistChannelCount << endl;
    cout << myname << "   AllPlotFileName: " << m_AllPlotFileName << endl;
    cout << myname << "   MinPlotFileName: " << m_MinPlotFileName << endl;
    cout << myname << "   MaxPlotFileName: " << m_MaxPlotFileName << endl;
    cout << myname << " PhasePlotFileName: " << m_PhasePlotFileName << endl;
    cout << myname << "     PhaseVariable: " << m_PhaseVariable << endl;
    cout << myname << "     PhaseGrouping: " << m_PhaseGrouping << endl;
    cout << myname << "      RootFileName: " << m_RootFileName << endl;
    cout << myname << "      TreeFileName: " << m_TreeFileName << endl;
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
    cout << myname << "     PlotFrequency: " << m_PlotFrequency << endl;
    cout << myname << "    PhasePlotSizeX: " << m_PhasePlotSizeX << endl;
    cout << myname << "    PhasePlotSizeY: " << m_PhasePlotSizeY << endl;
    cout << myname << "   PhasePlotSplitX: " << m_PhasePlotSplitX << endl;
    cout << myname << "   PhasePlotSplitY: " << m_PhasePlotSplitY << endl;
    cout << myname << "      NTimingPhase: " << m_NTimingPhase << endl;
  }
  if ( m_LogLevel >=4 ) {
    cout << myname << "INFO: Checking state." << endl;
    cout << myname << "INFO: Full hist map size: " << state().ChannelTickModFullHists.size() << endl;
    cout << myname << "INFO: Proc hist map size: " << state().ChannelTickModProcHists.size() << endl;
  }
}

//**********************************************************************

AdcTickModViewer::~AdcTickModViewer() {
  const string myname = "AdcTickModViewer::dtor: ";
  if ( m_LogLevel >= 1 ) cout << myname << "Closing." << endl;
  Index nplot = 0;
  processAccumulation(nplot);
  if ( m_LogLevel >= 1 ) cout << myname << "Plot count: " << nplot << endl;
  if ( m_LogLevel >= 1 ) cout << myname << "TM initial hist count: " << state().tickModHistogramInitialCount << endl;
  if ( m_LogLevel >= 1 ) cout << myname << "TM rebuild hist count: " << state().tickModHistogramRebuildCount << endl;
  if ( m_LogLevel >= 1 ) cout << myname << "Plot count: " << nplot << endl;
}
  
//**********************************************************************

DataMap AdcTickModViewer::viewMap(const AdcChannelDataMap& acds) const {
  DataMap ret;
  if ( acds.size() == 0 ) return ret.setStatus(0);
  // Save channel count.
  //Index eventID = acds.begin()->second.event;
  //state().eventIDs.push_back(eventID);
  //state().nchans.push_back(acds.size());
  return AdcChannelTool::viewMap(acds);
}
 
//**********************************************************************

DataMap AdcTickModViewer::view(const AdcChannelData& acd) const {
  const string myname = "AdcTickModViewer::view: ";
  DataMap res;
  Index icha = acd.channel;
  if ( m_LogLevel >= 3 ) cout << myname << "Processing channel " << icha << endl;
  setChannelData(acd);
  HistVector& tmhsFull = state().ChannelTickModFullHists[icha];
  Index ntkm = m_TickModPeriod;
  if ( tmhsFull.size() == 0 ) tmhsFull.resize(ntkm, nullptr);
  Index eventID = acd.event;
  Index ndigi = 0;
  if ( m_varNchan ) {
    if ( acd.hasMetadata("ndigi") ) {
      ndigi = acd.getMetadata("ndigi");
    } else {
      cout << myname << "WARNING: Metadata not found in ADC data: " << "digi" << endl;
    }
  }
  double tick0 = 0.0;   // Tick number in this run/job
  Index itkm0 = 0;
  Index timingPhase = 0;
  if ( m_tickOffsetTool != nullptr ) {
    TimeOffsetTool::Data dat;
    dat.run = acd.run;
    dat.subrun = acd.subRun;
    dat.event = acd.event;
    dat.channel = acd.channel;
    dat.triggerClock = acd.triggerClock;
    TimeOffsetTool::Offset off = m_tickOffsetTool->offset(dat);
    if ( ! off.isValid() ) {
      cout << myname << "Error finding tick offset: " << off.status << endl;
      return res.setStatus(1);
    }
    if ( off.unit != "tick" ) {
      cout << myname << "Time offset has wrong unit: " << off.unit << endl;
      return res.setStatus(2);
    }
    long toff = off.value;
    if ( ! state().haveTick0Job ) {
      state().tick0Job = toff;
      state().haveTick0Job = true;
    }
    tick0 = toff - state().tick0Job;
    while ( toff < 0.0 ) toff += m_TickModPeriod;
    itkm0 = toff % m_TickModPeriod;
    if ( m_LogLevel >= 3 ) cout << myname << "     Tick0: " << tick0 << endl;
    if ( m_LogLevel >= 3 ) cout << myname << "  Tickmod0: " << itkm0 << endl;
    // Get timing phase.
    if ( m_NTimingPhase ) {
      double delta = 0.01/m_NTimingPhase;
      timingPhase = m_NTimingPhase*(off.rem + delta);
      if ( m_LogLevel >= 3 ) cout << myname << "  Timing phase: " << timingPhase << endl;
    }
  }
  // Set the event variables, the first time an event ID is encountered.
  if ( state().eventIDs.size() == 0 || state().eventIDs.back() != eventID ) {
    state().eventIDs.push_back(eventID);
    state().nchans.push_back(ndigi);
    state().tick0s.push_back(tick0);
    if ( m_LogLevel >= 3 ) cout << myname << "Creating event entry: EventID=" << eventID
                                << ", ndigi=" << ndigi << ", tick0=" << std::fixed << tick0 << endl;
  }
  // Process each tickmod for this channel.
  FloatVector sigs(ntkm, 0.0);
  Index itkmMax = 0;
  float sigMax = -1.e20;
  FloatVector sigMeans(ntkm, 0.0);
  for ( Index itkm=0; itkm<ntkm; ++itkm ) {
    float& sig = sigMeans[itkm];
    processChannelTickMod(acd, itkm0, itkm, sig);
    if ( sig > sigMax || itkm == 0 ) {
      itkmMax = itkm;
      sigMax = sig;
    }
  }
  // Find and record the tickmod position of the signal peak.
  // The peak position is found with three-point interpolation.
  if ( true ) {
    FloatVVector& mtmsForAllPhases = state().MaxTickMods[icha];
    Index ivar = timingPhase;
    if ( m_varPhase ) {
      if ( mtmsForAllPhases.size() < m_NTimingPhase ) mtmsForAllPhases.resize(m_NTimingPhase);
    } else {
      ivar = state().eventIDs.size() - 1;  // For now, assume no overlap in event processing.
      if ( mtmsForAllPhases.size() < ivar+1 ) mtmsForAllPhases.resize(ivar+1);
    }
    FloatVector& mtms = mtmsForAllPhases[ivar];
    float ym = itkmMax > 0 ? sigMeans[itkmMax-1] : sigMeans[ntkm-1];
    float y0 = sigMeans[itkmMax];
    float yp = itkmMax+1 < ntkm ? sigMeans[itkmMax+1] : sigMeans[0];
    float den = 2.0*y0 - yp - ym;
    if ( den <= 0.0 ) {
      cout << myname << "WARNING: " << "Peak not found for channel " << icha << endl;
      cout << myname << "itkmMax: " << itkmMax << endl;
      cout << myname << "sigs: " << ym << ", " << y0 << ", " << yp << endl;
    } else {
      float xrem = 0.5*(yp - ym)/den;
      float itkmPeak = itkmMax + xrem;
//cout << "XXX: " << itkmMax << ": (" << ym << ", " << y0 << ", " << yp << "): " << xrem << ", " << itkmPeak << endl;
      mtms.push_back(itkmPeak);
    }
  }
  // If requested, make the plots of accumulated data for this channel.
  Index nplot = 0;
  if ( m_PlotFrequency ) {
    processAccumulation(nplot);
    HistVectorMap& ctmprocs = state().ChannelTickModProcHists;
    if ( ctmprocs.find(icha) == ctmprocs.end() ) {
      static HistVector empty;
      res.setHistVector("tmHists", empty);
    } else {
      res.setHistVector("tmHists", ctmprocs[icha]);
    }
  }
  res.setHistVector("tmHists", tmhsFull);   // Passing out hist that will be updated!
  res.setInt("tmCount", ntkm);
  return res;
}

//**********************************************************************

void AdcTickModViewer::setChannelData(const AdcChannelData& acd) const {
  const string myname = "AdcTickModViewer::setChannelData: ";
  Index icha = acd.channel;
  if ( icha == AdcChannelData::badIndex ) {
    cout << myname << "WARNING: Invalid channel index." << endl;
    icha = 0;
  }
  // Find the index of the variable that provides the grouping
  // for phase plots.
  Index igrp = AdcChannelData::badIndex;
  if ( m_groupByChannel ) igrp = icha;
  if ( m_groupByFemb ) igrp = acd.fembID;
  if ( igrp == AdcChannelData::badIndex ) {
    cout << myname << "WARNING: Invalid phase channel grouping: " << m_PhaseGrouping << endl;
  }
  state().channel = icha;
  // First call with this channel, we copy the channel data and record
  // the channel-grouping mappings.
  if ( state().acdMap.find(icha) == state().acdMap.end() ) {
    if ( m_LogLevel >= 3 ) {
      cout << myname << "Adding ADC channel state for channel " << icha << endl;
    }
    AdcChannelData& newacd = state().acdMap[icha];
    copyAcd(acd, newacd);
    state().phaseIndexMap[icha] = igrp;
    state().phaseChannels[igrp].push_back(icha);
  }
}
  
//**********************************************************************

Index AdcTickModViewer::setChannel(Index icha) const {
  const string myname = "AdcTickModViewer::setChannel: ";
  if ( icha == AdcChannelData::badIndex ) {
    cout << myname << "ERROR: Invalid channel index." << endl;
    icha = 0;
  }
  if ( state().acdMap.find(icha) == state().acdMap.end() ) {
    cout << myname << "ERROR: There is no description for channel " << icha << endl;
    state().currentAcd.clear();
  } else {
    copyAcd(state().acdMap[icha], state().currentAcd);
  }
  state().channel = icha;
  return icha;
}

//**********************************************************************

void AdcTickModViewer::clearChannelIndex() const {
  const string myname = "AdcTickModViewer::getChannelIndex: ";
  state().channel = AdcChannelData::badIndex;
  state().currentAcd.clear();
}

//**********************************************************************

Index AdcTickModViewer::getChannelIndex() const {
  const string myname = "AdcTickModViewer::getChannelIndex: ";
  Index icha = state().channel;
  if ( icha == AdcChannelData::badIndex ) {
    cout << myname << "ERROR: Channel index requested before being set." << endl;
    return setChannel(0);
  }
  return icha;
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
AdcTickModViewer::processChannelTickMod(const AdcChannelData& acd, Index itkm0, Index itkm, float& sigMean) const {
  const string myname = "AdcTickModViewer::processChannelTickMod: ";
  const AdcCount border = 20;
  const AdcCount adcLim = 4096;
  DataMap res;
  Index nsam = acd.raw.size();
  if ( nsam == 0 ) {
    if ( m_LogLevel >= 2 ) cout << myname << "WARNING: Raw data is empty." << endl;
    return 1;
  }
  // Extract the data for this tickmod.
  Index period = m_TickModPeriod;
  TickModData adcData;
  Index isam0 = (itkm + period - itkm0) % period;
  for ( Index isam=isam0; isam<nsam; isam+=period ) adcData.add(acd.raw[isam]);
  // Fetch the histogram pointer.
  Index icha = getChannelIndex();
  HistPtr& ph = state().ChannelTickModFullHists[icha][itkm];
  bool haveHist(ph);
  // Check if we need new histogram.  I.e. if not present or insufficient range.
  bool needHist = ! haveHist;
  if ( ! needHist ) {
    AdcCount hstMin = ph->GetXaxis()->GetXmin() + 0.001;
    AdcCount hstMax = ph->GetXaxis()->GetXmax() + -0.999;
    needHist |= adcData.min() < hstMin;
    needHist |= adcData.max() > hstMax;
  }
  // If needed, create histogram.
  if ( needHist ) {
    string hname = nameReplace(m_HistName, acd, itkm);
    string htitl = nameReplace(m_HistTitle, acd, itkm);
    htitl += "; ADC count; # samples";
    AdcCount adcMin = adcData.min() < 2*border ? 0.0 : adcData.min() - border;
    AdcCount adcMax = (adcData.max() + 2*border > adcLim) ? adcLim : (adcData.max() + border);
    if ( haveHist ) {
      AdcCount hstMin = ph->GetXaxis()->GetXmin() + 0.001;
      AdcCount hstMax = ph->GetXaxis()->GetXmax() + -0.999;
      if ( hstMin < adcMin ) adcMin = hstMin;
      if ( hstMax > adcMax ) adcMax = hstMax;
    }
    float xmin = adcMin;
    float xmax = adcMax + 1.0;
    int nbin = xmax - xmin + 0.001;
    HistPtr phold = ph;
    if ( m_LogLevel >= 4 ) cout << myname << (haveHist ? "Re-c" : "C")
                                << "reating histogram " << hname
                                << " for channel " << acd.channel
                                << " tickmod " << itkm
                                << ": " << nbin << ": [" << xmin << ", " << xmax << ")" << endl;
    if ( nbin <= 0 ) {
      cout << myname << "  adcData.min(): " << adcData.min() << endl;
      cout << myname << "  adcData.max(): " << adcData.max() << endl;
      cout << myname << "         adcMin: " << adcMin << endl;
      cout << myname << "         adcMax: " << adcMax << endl;
      abort();
    }
    ph.reset(new TH1F(hname.c_str(), htitl.c_str(), nbin, xmin, xmax));
    ph->SetDirectory(0);
    if ( haveHist ) {
      for ( int ibin=1; ibin<=phold->GetNbinsX(); ++ibin ) {
        int adc = phold->GetXaxis()->GetXmin() - 0.999 + ibin;
        for ( int icnt=0; icnt<phold->GetBinContent(ibin); ++icnt ) {
          ph->Fill(adc);
        }
      }
      if ( m_LogLevel > 4 ) cout << myname << "Check new hist: " << phold->GetMean() << " ?= " << ph->GetMean() << endl;
      ++state().tickModHistogramRebuildCount;
    } else {
      ++state().tickModHistogramInitialCount;
    }
  }
  // Add the new data to the histogram.
  for ( AdcCount adc : adcData.data ) ph->Fill(adc);
  sigMean = adcData.mean();
  return 0;
}

//**********************************************************************

// Process the accumulated data for the current channel/femb.
// A StickyCodeMetric is created from the full ADC-range frequency histogram
// for each tickmod and from each of these:
// 1. The limited-range histogram is cached and optionally plotted and/or written
//    to the root file.
// 2. An entry with metrics is added to the metric tree.
int AdcTickModViewer::processAccumulatedChannel(Index& nplot) const {
  const string myname = "AdcTickModViewer::processAccumulatedChannel: ";
  nplot = 0;
  Index icha = getChannelIndex();
  if ( state().ChannelTickModFullHists.find(icha) == state().ChannelTickModFullHists.end() ) return 1;
  const HistVector& tmhsFull = state().ChannelTickModFullHists[icha];
  HistVector& tmhsProc = state().ChannelTickModProcHists[icha];
  Index ntkm = tmhsFull.size();
  if ( m_LogLevel >= 3 ) {
    cout << myname << "Tickmod hist count: " << ntkm << endl;
  }
  if ( tmhsProc.size() == 0 ) tmhsProc.resize(ntkm, nullptr);
  if ( m_plotAny || m_makeTree ) {
    // Fetch the tree.
    TTree* ptree = state().tickmodTree;
    TickModTreeData& data = state().treedata;
    if ( ptree != nullptr ) {
      data.run = state().currentAcd.run;
      data.chan = state().currentAcd.channel;
      data.femb = state().currentAcd.fembID;
      data.fembChan = state().currentAcd.fembChannel;
      data.pedestal = state().currentAcd.pedestal;
    }
    // Loop over tickmods and, for each, create metrics and limited-range ADC
    // frequency histo and fill metric tree.
    for ( Index itkm=0; itkm<ntkm; ++itkm ) {
      const HistPtr& ph = tmhsFull[itkm];
      // Process histograms with the sticky code utility.
      Index chmod = 10;
      StickyCodeMetrics scm(ph->GetName(), ph->GetTitle(), m_HistChannelCount, chmod, m_FitSigmaMin, m_FitSigmaMax);
      if ( scm.evaluate(ph.get()) ) {
        cout << myname << "Sticky code evaluation failed for channel " << icha
             << " tickmod " << itkm << endl;
        //tmhsProc[itkm].reset(ph);
      } else {
        tmhsProc[itkm] = scm.getSharedHist();
        if ( ptree != nullptr ) {
          data.itkm = itkm;
          data.fill(scm);
          ptree->Fill();
        }
      }
    }
  }
  if ( m_plotAny ) {
    // Draw the ADC frequency histogram for each tickmod.
    makeTickModPlots(nplot);
    if ( m_LogLevel >= 3 ) {
      cout << myname << "  Plot file count for channel " << icha
           << ": " << nplot << endl;
    }
  }
  return 0;
}

//**********************************************************************

// Process the accumulated data for all channels.
// The metric tree is created, processAccumulatedChannel is called for
// each channel and phase plots are created.

int AdcTickModViewer::processAccumulation(Index& nplot) const {
  const string myname = "AdcTickModViewer::processAccumulation: ";
  Index ncha = state().ChannelTickModFullHists.size();
  if ( ncha == 0 ) return 0;
  // Set currentAcd for tree file name.
  setChannel(state().ChannelTickModFullHists.begin()->first);
  // Create tree to hold results.
  TTree*& ptree = state().tickmodTree;
  TFile*& pfile = state().pfile;
  if ( m_TreeFileName.size() ) {
    if ( m_LogLevel >= 2 ) cout << myname << "Creating tickmod tree." << endl;
    TDirectory* psavdir = gDirectory;
    string tfname = nameReplace(m_TreeFileName, state().currentAcd, 0);
    pfile = TFile::Open(tfname.c_str(), "RECREATE");
    if ( pfile->IsOpen() ) {
      ptree = new TTree("tickmod", "TickMod tree");
      ptree->Branch("data", &(state().treedata), 64000, 1);
    } else {
      cout << myname << "Unable to open file " << m_TreeFileName << endl;
    }
    psavdir->cd();
  }
  Index nplotTot = 0;
  int rstat = 0;
  if ( m_LogLevel >= 2 ) cout << myname << "Processing " << ncha << " channels" << endl;
  for ( HistVectorMap::value_type icvm : state().ChannelTickModFullHists ) {
    Index icha = setChannel(icvm.first);
    if ( m_LogLevel >= 3 ) {
      cout << myname << "Processing channel " << icha << endl;
    }
    Index nplot = 0;
    rstat += processAccumulatedChannel(nplot);
    nplotTot += nplot;
  }
  clearChannelIndex();
  if ( m_LogLevel >= 2 ) {
    cout << myname << "  Total plot file count: " << nplotTot << endl;
  }
  if ( pfile != nullptr ) {
    if ( m_LogLevel >= 2 ) {
      cout << myname << "Tree size: " << ptree->GetEntries() << endl;
    }
    pfile->Write();
    pfile->Close();
    delete pfile;
    pfile = nullptr;
    ptree = nullptr;
  }
  // Make phase graphs and plots.
  makePhaseGraphs();
  plotPhaseGraphs();
  return rstat;
}

//**********************************************************************

int AdcTickModViewer::makeTickModPlots(Index& nplot) const {
  const string myname = "AdcTickModViewer::makeTickModPlots: ";
  nplot = 0;
  // Exit if no plots are requested.
  if ( !m_plotAll && !m_plotMin && !m_plotMax ) return 0;
  if ( !m_plotAll && !m_plotMin && !m_plotMax && !m_plotPhase ) return 0;
  Index icha = getChannelIndex();
  // Exit if this channel should not be plotted.
  if ( m_PlotChannels.size() ) {
    if ( find(m_PlotChannels.begin(), m_PlotChannels.end(), icha) == m_PlotChannels.end() ) {
      if ( m_LogLevel >= 3 )  cout << myname << "Skipping channel not in PlotChannels: " << icha << endl;
      return 0;
    }
  }
  // Find the pad counts. Plot has npady x npadx pads with one tickmod plot per pad.
  Index npad = 0;
  Index npadx = 0;
  Index npady = 0;
  if ( m_PlotSplitX > 0 ) {
    npadx = m_PlotSplitX;
    npady = m_PlotSplitY ? m_PlotSplitY : m_PlotSplitX;
    npad = npadx*npady;
  }
  // Find the min and max ticks.
  Index ntkm = m_TickModPeriod;
  const HistVector& tmhs = state().ChannelTickModProcHists[icha];
  Index itkmMin = 9999;
  Index itkmMax = 9999;
  if ( true ) {
    double meanMin = tmhs[0]->GetMean();
    double meanMax = meanMin;
    for ( Index itkm=0; itkm<ntkm; ++itkm ) {
      double mean = tmhs[itkm]->GetMean();
      if ( mean > meanMax ) {
        meanMax = mean;
        itkmMax = itkm;
      }
      if ( mean < meanMin ) {
        meanMin = mean;
        itkmMin = itkm;
      }
    }
  }
  if ( m_LogLevel >= 3 ) cout << myname << "Processing channel " << icha << endl;
  // Build the vectors describing the plots.
  //   showTickModVectors - Vector of tickmods included in each tickmod plot.
  //   showPlotNames - Name for each tickmod plot file.
  vector<IndexVector> showTickModVectors;  // Vectors of tickmods to plot
  vector<Name> showPlotNames;
  // Add plot descriptions for all tickmods.
  if ( m_plotAll ) {
    if ( m_LogLevel >= 3 ) cout << myname << "Add full vector of "
                                << ntkm << " tickmods." << endl;
    showTickModVectors.emplace_back(ntkm);
    for ( Index itkm=0; itkm<ntkm; ++itkm ) showTickModVectors.back()[itkm] = itkm;
    showPlotNames.push_back(m_AllPlotFileName);
  }
  // Add plot descriptions for min and max plots.
  if ( m_plotMin || m_plotMax ) {
    IndexVector tkmsMin;
    IndexVector tkmsMax;
    int idel1 = npad/2;
    idel1 *= -1;
    int idel2 = idel1 + npad;
    if ( npad > ntkm ) {
      idel1 = 0;
      idel2 = ntkm;
    }
    if ( m_LogLevel >= 3 ) cout << myname << "Tick delta range: [" << idel1
                                << ", " << idel2 << ")" << endl;
    for ( int idel=idel1; idel<idel2; ++idel ) {
      Index itkm = (itkmMin + ntkm + idel) % ntkm;
      tkmsMin.push_back(itkm);
      itkm = (itkmMax + ntkm + idel) % ntkm;
      tkmsMax.push_back(itkm);
    }
    if ( m_plotMin ) {
      if ( m_LogLevel >= 3 ) cout << myname << "Adding min vector of " << tkmsMin.size()
                                  << " tickmods." << endl;
      showTickModVectors.push_back(tkmsMin);
      showPlotNames.push_back(m_MinPlotFileName);
    }
    if ( m_plotMax ) {
      if ( m_LogLevel >= 3 ) cout << myname << "Adding max vector of " << tkmsMax.size()
                                  << " tickmods." << endl;
      showTickModVectors.push_back(tkmsMax);
      showPlotNames.push_back(m_MaxPlotFileName);
    }
  }
  // Loop over descriptions and build plots.
  TPadManipulator* pmantop = nullptr;
  for ( Index ihv=0; ihv<showTickModVectors.size(); ++ihv ) {
    Name plotFileName;
    const IndexVector tkms = showTickModVectors[ihv];
    Name pfname = showPlotNames[ihv];
    if ( m_LogLevel >= 3 ) cout << myname << "Plotting " << tkms.size() << " tickmods with name "
                                << pfname << endl;
    Index ipad = 0;
    Index icount = 0;
    vector<TObject*> managedObjects;
    for ( Index itkm : tkms ) {
      HistPtr ph = tmhs[itkm];
      if ( pmantop == nullptr ) {
        pmantop = new TPadManipulator;
        if ( m_PlotSizeX && m_PlotSizeY ) pmantop->setCanvasSize(m_PlotSizeX, m_PlotSizeY);
        if ( npad > 1 ) pmantop->split(npadx, npady);
        plotFileName = nameReplace(pfname, state().currentAcd, itkm);
      }
      TPadManipulator* pman = pmantop->man(ipad);
      pman->add(ph.get(), "hist", false);
      if ( m_PlotShowFit > 1 ) pman->addHistFun(1);
      if ( m_PlotShowFit ) {
        pman->addHistFun(0);
        TF1* pfun = dynamic_cast<TF1*>(ph->GetListOfFunctions()->At(0));
        if ( pfun != nullptr ) {
          ostringstream ssout;
          ssout.precision(1);
          ssout << "Mean: " << std::fixed << pfun->GetParameter(1);
          TLatex* ptxt = new TLatex(0.72, 0.86, ssout.str().c_str());
          ptxt->SetNDC();
          ptxt->SetTextFont(42);
          pman->add(ptxt);
          managedObjects.push_back(ptxt);
          ssout.str("");
          ssout << "Sigma: " << std::fixed << pfun->GetParameter(2);
          ptxt = new TLatex(0.72, 0.80, ssout.str().c_str());
          ptxt->SetNDC();
          ptxt->SetTextFont(42);
          pman->add(ptxt);
          managedObjects.push_back(ptxt);
        }
      }
      pman->addVerticalModLines(64);
      pman->showUnderflow();
      pman->showOverflow();
      ++icount;
      if ( ++ipad == npad || icount == tkms.size() ) {
        if ( m_LogLevel >= 2 ) cout << myname << "Saving plot " << plotFileName << endl;
        pmantop->print(plotFileName);
        ++nplot;
        ipad = 0;
        delete pmantop;
        pmantop = nullptr;
        for ( TObject* pobj : managedObjects ) delete pobj;
        managedObjects.clear();
      }
    }
  }
  return 0;
}

//**********************************************************************

// Make the phase graphs.

int AdcTickModViewer::makePhaseGraphs() const {
  const string myname = "AdcTickModViewer::makePhaseGraph: ";
  if ( ! m_plotPhase ) return 0;
  // Plot phase vs. peak tick.
  Index nvar = m_NTimingPhase;
  if ( ! m_varPhase ) nvar = state().eventIDs.size();
  // Smearing of data points.
  Index ntkm = m_TickModPeriod;
  if ( ntkm == 0 ) return 0;
  // Loop over group indices (group = channel. femb, ...).
  Index ncan = state().phaseChannels.size();
  if ( m_LogLevel >= 2 ) {
    cout << myname << "Building phase-peak graphs for " << ncan
         << " " << m_PhaseGrouping << " indices" << endl;
  }
  for ( IndexVectorMap::value_type iich : state().phaseChannels ) {
    Index igrp = iich.first;
    const IndexVector ichas = iich.second;
    setChannel(ichas.front());
    // Collect the data for this channel from all contributing channels.
    FloatVVector xpksByPhase(nvar);
    for ( Index icha : ichas ) {
      const FloatVVector& chanXpksByPhase = state().MaxTickMods[icha];
      for ( Index ivar=0; ivar<nvar; ++ivar ) {
        FloatVector& out = xpksByPhase[ivar];
        const FloatVector& in = chanXpksByPhase[ivar];
        out.insert(out.end(), in.begin(), in.end());
      }
    }
    // Build the graph.
    string hnam = "hphtm";
    string vname = m_varPhase ? "Phase" :
                   m_varEvent ? "Event" :
                   m_varTick0 ? "Tick0" :
                   m_varNchan ? "Nchan" : "UNKNOWN";
    string httl = vname + " vs. tickmod peak for run %RUN% " + m_PhaseGrouping
                  + " " + std::to_string(igrp);
    httl = nameReplace(httl, state().currentAcd, 0);
    TGraph* pg = new TGraph;
    pg->SetName(hnam.c_str());
    pg->SetTitle(httl.c_str());
    pg->SetMarkerStyle(2);
    pg->SetMarkerColor(602);
    pg->GetXaxis()->SetTitle("Peak position [tick]");
    pg->GetXaxis()->SetTitle("Peak position [tick]");
    pg->GetXaxis()->SetTitle("Peak position [tick]");
    if ( m_varPhase ) pg->GetYaxis()->SetTitle("Timing phase");
    if ( m_varEvent ) pg->GetYaxis()->SetTitle("Event ID");
    if ( m_varTick0 ) pg->GetYaxis()->SetTitle("Event time [Tick]");
    if ( m_varNchan ) pg->GetYaxis()->SetTitle("Channel count");
    if ( m_LogLevel >=3 ) cout << myname << "Plotting " << m_PhaseVariable << " vs. peak tick for "
                               << m_PhaseGrouping << " " << igrp << endl;
    // Find range of tickmods with no change and with smaller values shifted
    // up by one period. The distribution with smaller range is used.
    double xmin = 1.e20;
    double xmax = -1.e20;
    double xminShift = 1.e20;
    double xmaxShift = -1.e20;
    float xhalf = 0.5*ntkm;
    for ( Index ivar=0; ivar<nvar; ++ivar ) {
      for ( float x : xpksByPhase[ivar] ) {
        if ( x < xmin ) xmin = x;
        if ( x > xmax ) xmax = x;
        float xs = x < xhalf ? x + ntkm : x;
        if ( xs < xminShift ) xminShift = xs;
        if ( xs > xmaxShift ) xmaxShift = xs;
      }
    }
    bool doShift = (xmax - xmin) > xhalf;
    doShift &= xmaxShift - xminShift < xmax - xmin;  // Should we shift points?
    if ( doShift ) {
      xmin = xminShift;
      xmax = xmaxShift;
    }
    // Fill the graph.
    for ( Index ivar=0; ivar<nvar; ++ivar ) {
      float y0 = ivar;
      if ( m_varEvent ) y0 = state().eventIDs[ivar];
      if ( m_varTick0 ) y0 = state().tick0s[ivar];
      if ( m_varNchan ) y0 = state().nchans[ivar];
      FloatVector& xpks = xpksByPhase[ivar];
      if ( m_LogLevel >= 4 ) cout << myname << "Phase " << ivar << " peak count: " << xpks.size() << endl;
      Index nxpk = xpks.size();
      for ( Index ixpk=0; ixpk<nxpk; ++ixpk ) {
        float x = xpks[ixpk];
        if ( doShift && x < xhalf ) x += ntkm;
        float y = y0;
        if ( m_varPhase || m_varEvent ) y+= -0.25 + 0.5*(ixpk+0.5)/nxpk;  // Smear integers
        pg->SetPoint(pg->GetN(), x, y);
      }
    }
    pg->GetXaxis()->SetRangeUser(xmin, xmax);
    // Save the graph in the state for this tool.
    state().phaseGraphs[igrp].reset(pg);
  }
  return 0;
}

//**********************************************************************

// Plot the phase graphs.

int AdcTickModViewer::plotPhaseGraphs() const {
  const string myname = "AdcTickModViewer::plotPhaseGraphs: ";
  if ( ! m_plotPhase ) return 0;
  if ( m_LogLevel >=3 ) cout << myname << "Making phase graph plots. Count: "
                             << state().phaseGraphs.size() << endl;
  Index npha = m_NTimingPhase;
  // Find the pad counts. Plot has npady x npadx pads with one tickmod plot per pad.
  Index npad = 0;
  Index npadx = 0;
  Index npady = 0;
  if ( m_PhasePlotSplitX > 0 ) {
    npadx = m_PhasePlotSplitX;
    npady = m_PhasePlotSplitY ? m_PhasePlotSplitY : m_PhasePlotSplitX;
    npad = npadx*npady;
  }
  // Find the min and max ticks.
  TPadManipulator* pmantop = nullptr;
  Index igrpLast = state().phaseGraphs.rbegin()->first;
  Index ipad = 0;
  Name plotFileName;
  for ( GraphMap::value_type icgr : state().phaseGraphs ) {
    Index igrp = icgr.first;
    setChannel(state().phaseChannels[igrp].front());
    TGraph* pg = icgr.second.get();
    // If needed, start a new plot.
    if ( pmantop == nullptr ) {
      pmantop = new TPadManipulator;
      if ( m_PhasePlotSizeX && m_PhasePlotSizeY )
        pmantop->setCanvasSize(m_PhasePlotSizeX, m_PhasePlotSizeY);
      if ( npad > 1 ) pmantop->split(npadx, npady);
      plotFileName = nameReplace(m_PhasePlotFileName, state().currentAcd, 0);
      ipad = 0;
    }
    TPadManipulator* pman = pmantop->man(ipad);
    // Evaluate plot range.
    double xmin = int(pg->GetXaxis()->GetXmin());
    double xmax = int(pg->GetXaxis()->GetXmax() + 1.0);
    while ( xmax - xmin < 10.0 ) {
      xmin -= 0.5;
      xmax += 0.5;
    }
    // Create plot.
    pman->add(pg->Clone(), "P");
    pman->addAxis();
    pman->setRangeX(xmin, xmax);
    if ( m_varPhase ) pman->setRangeY(0, npha);
    pman->setGridX();
    pman->setGridY();
    // Print the plot.
    if ( ++ipad == npad || igrp == igrpLast ) {
      if ( m_LogLevel >= 2 ) cout << myname << "Saving plot " << plotFileName << endl;
      pmantop->print(plotFileName);
      pmantop = nullptr;
    }
  }
  clearChannelIndex();
  return 0;
}

//**********************************************************************
