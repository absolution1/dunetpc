// AdcChannelMetric_tool.cc

#include "AdcChannelMetric.h"
#include <iostream>
#include <sstream>
#include <iomanip>
#include <set>
#include "dune/DuneInterface/Data/RunData.h"
#include "dune/DuneInterface/Tool/RunDataTool.h"
#include "dune/DuneInterface/Tool/AdcChannelStringTool.h"
#include "dune/DuneInterface/Tool/IndexRangeTool.h"
#include "dune/DuneInterface/Tool/FloatArrayTool.h"
#include "dune/DuneCommon/Utility/TPadManipulator.h"
#include "dune/DuneCommon/Utility/StringManipulator.h"
#include "dune/DuneCommon/Utility/LineColors.h"
#include "dune/DuneCommon/Utility/RootParFormula.h"
#include "dune/ArtSupport/DuneToolManager.h"
#include "larevt/CalibrationDBI/Interface/ChannelStatusService.h"
#include "larevt/CalibrationDBI/Interface/ChannelStatusProvider.h"
#include "TH1F.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TCanvas.h"
#include "TColor.h"
#include "TStyle.h"
#include "TDirectory.h"
#include "TFile.h"

using std::string;
using std::cout;
using std::setw;
using std::cin;
using std::endl;
using std::vector;
using std::istringstream;
using Index = AdcChannelMetric::Index;
using TGraphVector = std::vector<TGraph*>;
using TGraphErrorsVector = std::vector<TGraphErrors*>;

//**********************************************************************
// local definitiions.
//**********************************************************************

//**********************************************************************
// Sublass methods.
//**********************************************************************

bool AdcChannelMetric::AdcChannelMetric::State::update(Index run, Index event) {
  bool runHasChanged = false;
  if ( callCount == 0 ) {
    firstRun = run;
    firstEvent = event;
    runCount = 1;
    eventCount = 1;
    runHasChanged = true;
  } else {
    if ( run != lastRun ) {
      ++runCount;
      runHasChanged = true;
    }
    if ( event != lastEvent ) ++eventCount;
  }
  ++callCount;
  lastEvent = event;
  lastRun = run;
  return runHasChanged;
}

//**********************************************************************
// Class methods.
//**********************************************************************

AdcChannelMetric::AdcChannelMetric(fhicl::ParameterSet const& ps)
: m_LogLevel(ps.get<int>("LogLevel")), 
  m_Metric(ps.get<Name>("Metric")),
  m_DataView(ps.get<Name>("DataView")),
  m_PedestalReference(ps.get<Name>("PedestalReference")),
  m_MetricSummaryView(ps.get<Name>("MetricSummaryView")),
  m_ChannelRanges(ps.get<NameVector>("ChannelRanges")),
  m_MetricMin(new RootParFormula("MetricMin", ps.get<Name>("MetricMin"))),
  m_MetricMax(new RootParFormula("MetricMax", ps.get<Name>("MetricMax"))),
  m_MetricBins(ps.get<Index>("MetricBins")),
  m_ChannelLineModulus(ps.get<Index>("ChannelLineModulus")),
  m_ChannelLinePattern(ps.get<IndexVector>("ChannelLinePattern")),
  m_HistName(ps.get<Name>("HistName")),
  m_HistTitle(ps.get<Name>("HistTitle")),
  m_MetricLabel(ps.get<Name>("MetricLabel")),
  m_PlotSizeX(ps.get<Index>("PlotSizeX")),
  m_PlotSizeY(ps.get<Index>("PlotSizeY")),
  m_PlotFileName(ps.get<Name>("PlotFileName")),
  m_PlotUsesStatus(ps.get<int>("PlotUsesStatus")),
  m_RootFileName(ps.get<Name>("RootFileName")),
  m_MetadataFlags(ps.get<NameVector>("MetadataFlags")),
  m_mdRead(false), m_mdWrite(false), m_mdWarnAbsent(false), m_mdWarnPresent(false),
  m_prdtool(nullptr),
  m_doSummary(false),
  m_doSummaryError(false),
  m_pPedestalReference(nullptr),
  m_pChannelStatusProvider(nullptr),
  m_state(new State) {
  const string myname = "AdcChannelMetric::ctor: ";
  string stringBuilder = "adcStringBuilder";
  DuneToolManager* ptm = DuneToolManager::instance();
  // Process metadata flags.
  for ( string mflag : m_MetadataFlags ) {
    if ( mflag == "read" ) m_mdRead = true;
    else if ( mflag == "write" ) m_mdWrite = true;
    else if ( mflag == "warnabsent" ) m_mdWarnAbsent = true;
    else if ( mflag == "warnpresent" ) m_mdWarnPresent = true;
    else {
      cout << myname << "WARNING: Ignoring invalid metadata flag: " << mflag << endl;
    }
  }
  // Fetch the channel ranges.
  bool toolNotFound = false;
  const IndexRangeTool* pcrt = nullptr;
  for ( Name crn : m_ChannelRanges.size() ? m_ChannelRanges : NameVector(1, "") ) {
    if ( crn.size() == 0 ) {
      m_crs.emplace_back("all", 0, 0, "All");
    } else {
      if ( pcrt == nullptr && !toolNotFound ) {
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
  // Initialize the state.
  for ( const IndexRange& cr : m_crs ) getState().crsums[cr].resize(cr.size());
  // Fetch the pedestalreference tool.
  if ( m_PedestalReference.size() ) {
    m_pPedestalReference = ptm->getShared<FloatArrayTool>(m_PedestalReference);
    if ( m_pPedestalReference == nullptr ) {
      cout << "WARNING: Pedestal reference tool not found: " << m_PedestalReference << endl;
    } else {
      Index nref = m_pPedestalReference->size();
      Index off = m_pPedestalReference->offset();
      cout << myname << "Pedestal reference array has " << nref << " value"
           << ( nref == 1 ? "" : "s" );
      if ( nref ) cout << " starting at channel " << off;
      cout << "." << endl;
    }
  }
  // Fetch the naming tool.
  m_adcStringBuilder = ptm->getShared<AdcChannelStringTool>(stringBuilder);
  if ( m_adcStringBuilder == nullptr ) {
    cout << myname << "WARNING: AdcChannelStringTool not found: " << stringBuilder << endl;
  }
  // Fetch the run data tool.
  if ( m_MetricMin->npar() || m_MetricMax->npar() ) {
    string stnam = "runDataTool";
    m_prdtool = ptm->getShared<RunDataTool>(stnam);
    if ( m_prdtool == nullptr ) {
      cout << myname << "ERROR: RunDataTool " << stnam
           << " not found. Metric limits will not be evaluated." << endl;
    } else {
      cout << myname << "RunDataTool retrieved." << endl;
    }
  }
  // Set summary fields.
  m_doSummary = m_HistName.find("EVENT%") == string::npos;
  if ( m_doSummary ) {
    const std::set<Name> sumVals = {"count", "mean", "rms", "drms"};
    if ( m_MetricSummaryView.size() == 0 ) {
      m_MetricSummaryView = "mean:dmean";
      cout << myname << "WARNING: Missing metric summary view set to \"" << m_MetricSummaryView
           << "\"." << endl;
    }
    Name vnam = m_MetricSummaryView;
    Name enam;
    Name::size_type ipos = vnam.find(":");
    if ( ipos != Name::npos ) {
      enam = vnam.substr(ipos+1);
      vnam = vnam.substr(0, ipos);
    }
    if ( ! MetricSummary::isValueName(vnam) ) {
      cout << myname << "WARNING: Invalid value for metric summary view reset from " << vnam
           << " to mean." << endl;
      vnam = "mean";
    }
    if ( enam.size() ) {
      if ( ! MetricSummary::isValueName(enam) ) {
        cout << myname << "WARNING: Ignoring invalid error for metric summary view: " << enam << endl;
      } else {
        m_doSummaryError = true;
      }
    }
    m_summaryValue = vnam;
    m_summaryError = enam;
  }
  // Fetch the channel status service.
  if ( m_PlotUsesStatus ) {
    if ( m_LogLevel >= 1 ) cout << myname << "Fetching channel status service." << endl;
    m_pChannelStatusProvider = &art::ServiceHandle<lariov::ChannelStatusService>()->GetProvider();
    if ( m_pChannelStatusProvider == nullptr ) {
      cout << myname << "WARNING: Channel status provider not found." << endl;
      m_PlotUsesStatus = false;
    }
  }
  // Display the configuration.
  if ( m_LogLevel ) {
    cout << myname << "Configuration: " << endl;
    cout << myname << "            LogLevel: " << m_LogLevel << endl;
    cout << myname << "              Metric: " << m_Metric << endl;
    cout << myname << "            DataView: " << m_DataView << endl;
    cout << myname << "   PedestalReference: " << m_PedestalReference << endl;
    cout << myname << "   MetricSummaryView: " << m_MetricSummaryView;
    if ( m_summaryValue.size() ) {
      cout << " (" << m_summaryValue;
      if ( m_summaryError.size() ) cout << " +/- " << m_summaryError;
      cout << ")";
    }
    cout << endl;
    cout << myname << "       ChannelRanges: [";
    bool first = true;
    for ( const IndexRange& ran : m_crs ) {
      if ( ! first ) cout << ", ";
      else first = false;
      cout << ran.name;
    }
    cout << "]" << endl;
    cout << myname << "           MetricMin: " << m_MetricMin->formulaString() << endl;
    cout << myname << "           MetricMax: " << m_MetricMax->formulaString() << endl;
    cout << myname << "          MetricBins: " << m_MetricBins << endl;
    cout << myname << "  ChannelLineModulus: " << m_ChannelLineModulus << endl;
    cout << myname << "  ChannelLinePattern: {";
    first = true;
    for ( Index icha : m_ChannelLinePattern ) {
      if ( ! first ) cout << ", ";
      first = false;
      cout << icha;
    }
    cout << "}" << endl;
    cout << myname << "  ChannelCounts: {";
    first = true;
    for ( Index icha : m_ChannelCounts ) {
      if ( ! first ) cout << ", ";
      first = false;
      cout << icha;
    }
    cout << "}" << endl;
    cout << myname << "           PlotSizeX: " << m_PlotSizeX << endl;
    cout << myname << "           PlotSizeY: " << m_PlotSizeY << endl;
    cout << myname << "            HistName: " << m_HistName << endl;
    cout << myname << "           HistTitle: " << m_HistTitle << endl;
    cout << myname << "         MetricLabel: " << m_MetricLabel << endl;
    cout << myname << "        PlotFileName: " << m_PlotFileName << endl;
    cout << myname << "      PlotUsesStatus: " << m_PlotUsesStatus << endl;
    cout << myname << "        RootFileName: " << m_RootFileName << endl;
    cout << myname << "       MetadataFlags: [";
    first = true;
    for ( string mflag : m_MetadataFlags ) {
      if ( ! first ) cout << ", ";
      first = false;
      cout << mflag;
    }
    cout << "]" << endl;
  }
  // We might have to move this to getMetric.
  initialize();
}

//**********************************************************************

AdcChannelMetric::~AdcChannelMetric() {
  const string myname = "AdcChannelMetric::dtor: ";
  Index ncha = 0;
  Index nchaData = 0;
  Index nchaDataMax = 0;
  Index countMax = 0;
  Index weightSumMax = 0.0;
  for ( const MetricSummaryMap::value_type& imsm : getState().crsums ) {
    IndexRange cr = imsm.first;
    const MetricSummaryVector& msums = imsm.second;
    if ( m_LogLevel >= 3 ) {
      cout << myname << "Channel range " << cr.name << endl;
    }
    for ( Index kcha=0; kcha<cr.size(); ++kcha ) {
      const MetricSummary& ms = msums[kcha];
      ++ncha;
      if ( ms.eventCount ) {
        ++nchaData;
        if ( ms.eventCount > countMax ) {
          countMax = ms.eventCount;
          nchaDataMax = 1;
        } else if ( ms.eventCount == countMax ) {
          ++nchaDataMax;
        }
        if ( ms.weightSum > weightSumMax ) {
          weightSumMax = ms.weightSum;
        }
      }
    }
    if ( m_doSummary ) {
      MetricMap mets;
      const MetricSummaryVector& msums = getState().crsums[cr];
      for ( Index kcha=0; kcha<cr.size(); ++kcha ) {
        Index icha = cr.first() + kcha;
        const MetricSummary& msum = msums[kcha];
        if ( msum.eventCount ) {
          Metric& met = mets[icha];
          if ( m_summaryValue.size() ) {
            met.setValue(msum.getValue(m_summaryValue));
            if ( m_LogLevel >=3 ) cout << myname << "Channel: summary value: " << icha << ": " << met.value;
            if ( m_summaryError.size() ) {
              met.setError(msum.getValue(m_summaryError));
              if ( m_LogLevel >=3 ) cout << " +/- " << met.error;
            }
            if ( m_LogLevel >=3 ) cout << endl;
          }
        }
      }
      AdcChannelData acd;
      if ( getState().runCount == 1 ) {
        acd.setEventInfo(
          getState().firstRun,
          getState().eventCount==1 ? getState().firstEvent : AdcChannelData::badIndex()
        );
      }
      Name ofpname = nameReplace(m_PlotFileName, acd, cr);
      Name ofrname = nameReplace(m_RootFileName, acd, cr);
      TH1* ph = createHisto(acd, cr);
      processMetricsForOneRange(cr, mets, ph, ofpname, ofrname, true);
    }
  }
  if ( m_LogLevel >= 1 ) {
    Index w = 1;
    Index valmax = std::max(ncha, getState().callCount);
    if ( ncha ) w = log10(valmax) + 1.01;
    cout << myname << "Summary for metric " << m_Metric << endl;
    cout << myname << "                           # inits: " << setw(w) << getState().initCount << endl;
    cout << myname << "                           # calls: " << setw(w) << getState().callCount << endl;
    cout << myname << "                          # events: " << setw(w) << getState().eventCount << endl;
    cout << myname << "                            # runs: " << setw(w) << getState().runCount << endl;
    cout << myname << "   Maximum # entries for a channel: " << setw(w) << countMax << endl;
    cout << myname << "  Maximum weight sum for a channel: " << setw(w) << weightSumMax << endl;
    cout << myname << "        Total # channels in ranges: " << setw(w) << ncha << endl;
    cout << myname << "           # channels without data: " << setw(w) << ncha - nchaData << endl;
    cout << myname << "              # channels with data: " << setw(w) << nchaData << endl;
    cout << myname << "     # channels with max # entries: " << setw(w) << nchaDataMax << endl;
    if ( m_LogLevel >= 3 ) {
      Index icnt = 0;
      for ( const MetricSummaryMap::value_type& imsm : getState().crsums ) {
        IndexRange cr = imsm.first;
        cout << "============= Channel range " << cr.name << endl;
        const MetricSummaryVector& msums = imsm.second;
        for ( const MetricSummary& msum : msums ) {
          if ( msum.eventCount ) {
            cout << "------------- " << m_Metric << "[" << icnt++ << "]" << endl;
            cout << myname << "        weight flag: " << setw(w) << msum.weightFlag << endl;
            cout << myname << "           # events: " << setw(w) << msum.eventCount << endl;
            cout << myname << "  # weighted events: " << setw(w) << msum.weightedEventCount << endl;
            cout << myname << "               Neff: " << setw(w) << msum.neff() << endl;
            cout << myname << "               mean: " << setw(w) << msum.mean() << " +/- "
                                                      << setw(w) << msum.dmean() << endl;
            cout << myname << "                RMS: " << setw(w) << msum.rms() << " +/- "
                                                      << setw(w) << msum.drms() << endl;
            cout << myname << "             center: " << setw(w) << msum.center() << " +/- "
                                                      << setw(w) << 0.5*msum.range() << endl;
          }
        }
        cout << "===================================" << endl;
      }
    }
  }
}

//**********************************************************************

void AdcChannelMetric::initialize(bool force) {
  const string myname = "AdcChannelMetric::initialize: ";
  if ( !force && getState().initCount ) return;
  if ( m_LogLevel >= 1 ) cout << myname << "Initializing " << m_crs.size()
                              << " channel ranges." << endl;
  // Loop over channels and fetch status for each.
  Index ncha = 0;
  for ( const IndexRange& ran : m_crs ) {
    for ( Index icha =ran.begin; icha<ran.end; ++icha ) { 
      channelStatus(icha);
      ++ncha;
    }
  }
  if ( m_LogLevel >= 1 ) cout << myname << "Initialized " << ncha << " channels." << endl;
  ++getState().initCount;
}

//**********************************************************************

DataMap AdcChannelMetric::view(const AdcChannelData& acd) const {
  const string myname = "AdcChannelMetric::view: ";
  DataMap ret;
  double val = 0.0;
  double wt = 0.0;
  Name sunits;
  int rstat = getMetric(acd, m_Metric, val, sunits, wt);
  if ( rstat ) return ret.setStatus(rstat);
  ret.setString("metricName", m_Metric);
  ret.setFloat("metricValue", val);
  ret.setFloat("metricWeight", wt);
  ret.setString("metricUnits", sunits);
  return ret;
}

//**********************************************************************

DataMap AdcChannelMetric::viewMap(const AdcChannelDataMap& acds) const {
  const string myname = "AdcChannelMetric::viewMap: ";
  MetricMap mets;
  return viewMapLocal(acds, mets);
}

//**********************************************************************

DataMap AdcChannelMetric::update(AdcChannelData& acd) const {
  const string myname = "AdcChannelMetric::update: ";
  DataMap ret = view(acd);
  string mnam = m_Metric;
  if ( m_mdWrite ) {
    if ( m_mdWarnPresent && acd.hasMetadata(mnam) ) {
      cout << "WARNING: Run-event-channel " << acd.run() << "-" << acd.event() << "-"
           << acd.channel() << " already has metadata field " << mnam << endl;
    }
    float mval = ret.getFloat("metricValue");
    acd.setMetadata(mnam, mval);
  }
  return ret;
}

//**********************************************************************

DataMap AdcChannelMetric::updateMap(AdcChannelDataMap& acds) const {
  const string myname = "AdcChannelMetric::updateMap: ";
  MetricMap mets;
  DataMap ret = viewMapLocal(acds, mets);
  if ( m_mdWrite ) {
    string mnam = m_Metric;
    for ( const auto& imet : mets ) {
      Index icha = imet.first;
      AdcChannelDataMap::iterator iacd = acds.find(icha);
      if ( iacd == acds.end() ) {
        cout << myname << "ERROR: Channel map is missing channel " << icha << endl;
      } else {
        AdcChannelData& acd = iacd->second;
        if ( m_mdWarnPresent && acd.hasMetadata(mnam) ) {
          cout << "WARNING: Run-event-channel " << acd.run() << "-" << acd.event() << "-"
               << acd.channel() << " already has metadata field " << mnam << endl;
        }
        float mval = imet.second.value;
        acd.setMetadata(mnam, mval);
      }
    }
  }
  return ret;
}

//**********************************************************************

DataMap AdcChannelMetric::viewMapLocal(const AdcChannelDataMap& acds, MetricMap& mets) const {
  const string myname = "AdcChannelMetric::viewMapLocal: ";
  DataMap ret;
  if ( acds.size() == 0 ) {
    cout << myname << "Input channel map is empty." << endl;
    return ret.setStatus(1);
  }
  Index chanFirst = acds.begin()->first;
  Index chanLast = acds.rbegin()->first;
  if ( m_LogLevel >= 2 ) cout << myname << "Processing " << acds.size() << " channels: ["
                              << chanFirst << ", " << chanLast << "]" << endl;
  for ( IndexRange ran : m_crs ) {
    if ( ran.name == "all" ) {
      ran.begin = chanFirst;
      ran.end = chanLast + 1;
    }
    Index chanLo = std::max(chanFirst, ran.first());
    Index chanHi = std::min(chanLast, ran.last());
    if ( chanHi >= chanLo ) ret += viewMapForOneRange(acds, ran, mets);
  }
  return ret;
}

//**********************************************************************

DataMap AdcChannelMetric::
viewMapForOneRange(const AdcChannelDataMap& acds, const IndexRange& ran, MetricMap& mets) const {
  const string myname = "AdcChannelMetric::viewMapForOneRange: ";
  DataMap ret;
  vector<int> ichas;
  vector<float> vals;
  vector<float> wgts;
  // Extract the metrics for the subset of this range included in the data.
  Index icha0 = ran.begin;
  AdcChannelDataMap::const_iterator iacd1=acds.lower_bound(icha0);
  AdcChannelDataMap::const_iterator iacd2=acds.upper_bound(ran.last());
  if ( getState().update(iacd1->second.run(), iacd1->second.event()) ) {
    evaluateFormulas();
  }
  MetricSummaryVector& metricSums = getState().crsums[ran];
  if ( metricSums.size() < ran.size() ) metricSums.resize(ran.size());
  for ( AdcChannelDataMap::const_iterator iacd=iacd1; iacd!=iacd2; ++iacd ) {
    const AdcChannelData& acd = iacd->second;
    double met;
    double wt;
    Name sunits;
    int rstat = getMetric(acd, m_Metric, met, sunits, wt);
    if ( rstat ) {
      cout << myname << "WARNING: Metric evaluation failed for metric " << m_Metric
           << " channel " << acd.channel() << endl;
      continue;
    }
    Index icha = iacd->first;
    mets[icha].setValue(met);
    MetricSummary& metricSum = metricSums[icha-icha0];
    metricSum.add(met, wt);
    if ( m_LogLevel >=4 ) {
      cout << myname << "    Chan: met, wt (wsum, range): " << icha <<  ": " << met << ", " << wt
           << " (" << metricSum.weightSum << ", " << metricSum.range() << ")" << endl;
    }
    ichas.push_back(icha);
    vals.push_back(met);
    wgts.push_back(wt);
  }
  // Create the histogram for this data and this range.
  const AdcChannelData& acdFirst = acds.begin()->second;
  string ofpname;
  string ofrname;
  if ( ! m_doSummary ) {
    ofpname = nameReplace(m_PlotFileName, acdFirst, ran);
    ofrname = nameReplace(m_RootFileName, acdFirst, ran);
  }
  TH1* ph = createHisto(acdFirst, ran);
  // Fill the histogram and create the plots for this data and this range.
  processMetricsForOneRange(ran, mets, ph, ofpname, ofrname, false);
  ret.setHist(ph, true);
  ret.setString("metricName", m_Metric);
  ret.setIntVector("metricChannels", ichas);
  ret.setFloatVector("metricValues", vals);
  ret.setFloatVector("metricWeights", wgts);
  return ret;
}

//**********************************************************************

int AdcChannelMetric::getMetric(const AdcChannelData& acdtop, Name met, double& val,
                                Name& sunits, double& weight) const {
  const string myname = "AdcChannelMetric::getMetric: ";
  val = 0.0;
  weight = 0.0;
  sunits = "";
  AdcIndex nent = acdtop.viewSize(m_DataView);
  const AdcChannelData* pacd0 = acdtop.viewEntry(m_DataView, 0);
  if ( nent && pacd0 == nullptr ) {
    cout << myname << "ERROR: Channel " << acdtop.channel() << " view " << m_DataView
         << " is missing its first entry." << endl;
    return 4;
  }
  if ( nent ) weight = 1.0;
  if ( met == "pedestal" ) {
    val = nent ? pacd0->pedestal : 0;
    sunits = "ADC count";
  } else if ( met == "pedestalDiff" ) {
    val = nent ? pacd0->pedestal : 0.0;
    if ( m_pPedestalReference != nullptr ) {
      double pedRef = m_pPedestalReference->value(acdtop.channel(), 0.0);
      val -= pedRef;
    } else if ( m_PedestalReference == "first" && nent ) {
      Index icha = acdtop.channel();
      MetricMap& pedRefs = getState().pedRefs;
      MetricMap::const_iterator ipdr = pedRefs.find(icha);
      if ( ipdr == pedRefs.end() ) {
        pedRefs[icha].value = val;
        val = 0.0;
      } else {
        val -= ipdr->second.value;
      }
    } else {
      cout << myname << "WARNING: Pedestal difference requested without reference." << endl;
    }
    sunits = "ADC count";
  } else if ( met == "pedestalRms" ) {
    val = nent ? pacd0->pedestalRms : 0.0;
    sunits = "ADC count";
  } else if ( met == "time" ) {
    val = acdtop.time();
    sunits = "sec";
  } else if ( met == "fembID" ) {
    val = acdtop.fembID();
  } else if ( met == "fembChannel" ) {
    val = acdtop.fembChannel();
  } else if ( met == "apaFembID" ) {
    val = acdtop.fembID()%20;
  } else if ( met == "nraw" ) {
    val = 0.0;
    weight = 1.0;
    for ( Index ient=0; ient<nent; ++ient ) val += acdtop.viewEntry(m_DataView, ient)->raw.size();
  } else if ( met == "nsam" ) {
    val = 0.0;
    weight = 1.0;
    for ( Index ient=0; ient<nent; ++ient ) val += acdtop.viewEntry(m_DataView, ient)->samples.size();
  } else if ( met == "rawRms" ) {
    double sum = 0.0;
    AdcIndex nsam = 0;
    for ( Index ient=0; ient<nent; ++ient ) {
      const AdcChannelData* pacd = acdtop.viewEntry(m_DataView, ient);
      double ped = pacd->pedestal;
      for ( AdcSignal sig : pacd->raw ) {
        double dif = double(sig) - ped;
        sum += dif*dif;
        ++nsam;
      }
    }
    val = nsam == 0 ? 0.0 : sqrt(sum/nsam);
    weight = nsam;
  } else if ( met == "samRms" ) {
    double sum = 0.0;
    Index nsam = 0;
    for ( Index ient=0; ient<nent; ++ient ) {
      const AdcChannelData* pacd = acdtop.viewEntry(m_DataView, ient);
      for ( double sig : pacd->samples ) {
        sum += sig*sig;
        ++nsam;
      }
    }
    val = nsam == 0 ? 0.0 : sqrt(sum/nsam);
    weight = nsam;
  } else if ( met.substr(0, 6) == "samRms" ) {   // samRmsNN, NN is integer
    val = 0.0;
    weight = 0.0;
    istringstream sscnt(met.substr(6));
    Index ncnt = 0;
    sscnt >> ncnt;
    std::vector<float> samSums;  // Sum over samples for each group of cnt samples.
    if ( ncnt == 0 ) {
      cout << myname << "WARNING: Invalid metric: " << met << endl;
    } else {
      float samSum = 0.0;
      for ( Index ient=0; ient<nent; ++ient ) {
        const AdcChannelData* pacd = acdtop.viewEntry(m_DataView, ient);
        Index nsam = pacd->samples.size();
        Index samCount = 0;
        for ( Index isam=0; isam<nsam; ++isam ) {
          float sam = pacd->samples[isam];
          samSum += sam;
          ++samCount;
          if ( samCount == ncnt ) {
            samSums.push_back(samSum);
            samSum = 0.0;
            samCount = 0;
          }
        }
      }
    }
    if ( samSums.size() ) {
      float sum = 0.0;
      for ( float samSum : samSums ) sum += samSum*samSum;
      val = sqrt(sum/samSums.size());
      weight = samSums.size();
    }
  } else if ( met == "sigRms" || met == "nsgRms" ) {
    double sum = 0.0;
    Index nsum = 0;
    Index nsamtot = 0;
    bool doSignal = met == "sigRms";
    for ( Index ient=0; ient<nent; ++ient ) {
      const AdcChannelData* pacd = acdtop.viewEntry(m_DataView, ient);
      Index nsam = pacd->samples.size();
      if ( pacd->signal.size() != nsam ) {
        cout << myname << "WARNING: signal and sample sizes differ";
        if ( m_DataView.size() ) cout << " for view entry " << m_DataView << "[" << ient << "]";
        cout << ": " << pacd->signal.size() << " != " << nsam << "." << endl;
      } else {
        for ( Index isam=0; isam<nsam; ++isam ) {
          if ( pacd->signal[isam] == doSignal ) {
            double sig = pacd->samples[isam];
            sum += sig*sig;
            ++nsum;
          }
          ++nsamtot;
        }
      }
    }
    val = nsum == 0 ? 0.0 : sqrt(sum/nsum);
    weight = nsum;
    if ( m_LogLevel >= 4 ) {
      cout << myname << "Sample count for " << met << " for channel " << acdtop.channel();
      if ( m_DataView.size() ) cout << " view " << m_DataView;
      cout << ": " << nsum << "/" << nsamtot << "."
           << " " << met << " = " << val << endl;
    }
  // nsgRmsNN - Integrated noise with NN samples.
  } else if ( met.substr(0, 6) == "nsgRms" ) {
    val = 0.0;
    weight = 0.0;
    istringstream sscnt(met.substr(6));
    Index ncnt = 0;
    sscnt >> ncnt;
    std::vector<float> samSums;  // Sum over samples for each group of ncnt samples.
    if ( ncnt == 0 ) {
      cout << myname << "WARNING: Invalid metric: " << met << endl;
    } else {
      for ( Index ient=0; ient<nent; ++ient ) {
        const AdcChannelData* pacd = acdtop.viewEntry(m_DataView, ient);
        Index nsam = pacd->samples.size();
        if ( pacd->signal.size() != nsam ) {
          cout << myname << "WARNING: signal and sample sizes differ";
          if ( m_DataView.size() ) cout << " for view entry " << m_DataView << "[" << ient << "]";
          cout << " metric " << met << ": " << pacd->signal.size() << " != " << nsam << "." << endl;
          continue;
        }
        for ( Index isam0=0; isam0+ncnt<=nsam; ++isam0 ) {
          bool foundSignal = false;
          float samSum = 0.0;
          for ( Index icnt=0; icnt<ncnt; ++icnt ) {
            Index isam = isam0 + icnt;
            if ( pacd->signal[isam] ) {
              foundSignal = true;
              break;
            }
            float sam = pacd->samples[isam];
            samSum += sam;
          }
          if ( foundSignal ) break;
          samSums.push_back(samSum);
        }
      }
    }
    if ( samSums.size() ) {
      float sum = 0.0;
      for ( float samSum : samSums ) sum += samSum*samSum;
      val = sqrt(sum/samSums.size());
      weight = samSums.size();
    }
  } else if ( met == "sigFrac" ) {
    Index nsamtot = 0;
    Index nsig = 0;
    for ( Index ient=0; ient<nent; ++ient ) {
      const AdcChannelData* pacd = acdtop.viewEntry(m_DataView, ient);
      Index nsam = pacd->samples.size();
      if ( pacd->signal.size() != nsam ) {
        cout << myname << "WARNING: Signal and sample sizes differ for sigFrac";
        if ( m_DataView.size() ) cout << " view entry " << m_DataView << "[" << ient << "]";
        cout << ": " << pacd->signal.size() << " != " << nsam << "." << endl;
      } else {
        for ( Index isam=0; isam<nsam; ++isam ) {
          ++nsamtot;
          if ( pacd->signal[isam] ) ++nsig;
        }
      }
    }
    val = nsamtot > 0 ? float(nsig)/float(nsamtot) : 0.0;
  } else if ( met == "rawTailFraction" ) {
    Index ntail = 0;
    Index nsam = 0;
    for ( Index ient=0; ient<nent; ++ient ) {
      const AdcChannelData* pacd = acdtop.viewEntry(m_DataView, ient);
      double lim = 3.0*pacd->pedestalRms;
      double ped = pacd->pedestal;
      for ( AdcSignal sig : pacd->raw ) {
        double dif = double(sig) - ped;
        if ( fabs(dif) > lim ) ++ntail;
        ++nsam;
      }
    }
    val = nsam > 0 ? double(ntail)/nsam : 0.0;
  } else if ( nent == 1 && pacd0->hasMetadata(met) ) {
    if ( nent != 1 ) {
      cout << myname << "WARNING: Channel " << acdtop.channel()
           << " filling metadata only for the first of " << nent
           << " entries." << endl;
    }
    val = pacd0->metadata.find(met)->second;
  // Compound metric: met1+met2
  // TODO: Move this to ctor.
  } else if ( met.find("+") != string::npos ) {
    vector<string> nams;
    string metsrem = met;
    string::size_type ipos = 0;
    val = 0.0;
    while ( ipos != string::npos ) {
      ipos = metsrem.find("+");
      string newmet = metsrem.substr(0, ipos);
      double newval = 0.0;
      double newwt = 0.0;
      int sstat = getMetric(acdtop, newmet, newval, sunits, newwt);
      if ( sstat ) {
        cout << myname << "ERROR: Invalid sub-metric name: " << newmet << endl;
        return 2;
      }
      if ( newwt != 1.0 ) {
        cout << myname << "Evaluation of sub-metric " << newmet << " for metric " << met
             << " returned non-unity weight " << newwt << endl;
        return 3;
      }
      val += newval;
      if ( ipos == string::npos ) break;
      metsrem = metsrem.substr(ipos + 1);
    }
  } else {
    cout << myname << "ERROR: Invalid metric name: " << met << endl;
    return 1;
  }
  if ( m_LogLevel >= 4 ) {
    cout << myname << setw(20) << met << ": " << val << endl;
  }
  return 0;
}

//**********************************************************************

string AdcChannelMetric::
nameReplace(string name, const AdcChannelData& acd, const IndexRange& ran) const {
  StringManipulator sman(name, false);
  sman.replace("%CRNAME%", ran.name);
  sman.replace("%CRLABEL%", ran.label());
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

void AdcChannelMetric::
processMetricsForOneRange(const IndexRange& ran, const MetricMap& mets, TH1* ph,
                          Name ofpname, Name ofrname, bool useErrors) const {
  const string myname = "AdcChannelMetric::processMetricsForOneRange: ";
  unsigned int ngraph = m_PlotUsesStatus ? 4 : 1;
  NameVector statNames = {"All", "Good", "Bad", "Noisy"};
  LineColors lc;
  std::vector<int> statCols = {lc.blue(), lc.green(), lc.red(), lc.brown()};
  Name hname = ph->GetName();
  Name htitl = ph->GetTitle();
  // # channels vs. metric
  if ( m_MetricBins > 0 ) {
    if ( m_LogLevel >= 2 ) cout << myname << "Plotting # channels vs. metric. Count is " << mets.size() << endl;
    for ( MetricMap::value_type imet : mets ) {
      double met = imet.second.value;
      ph->Fill(met);
    }
    if ( ofpname.size() ) {
      TPadManipulator man;
      if ( m_PlotSizeX && m_PlotSizeY ) man.setCanvasSize(m_PlotSizeX, m_PlotSizeY);
      man.add(ph);
      man.addAxis();
      man.showUnderflow();
      man.showOverflow();
      man.print(ofpname);
    }
  } else {
    // Metric vs. channel.
    Name slaby = ph->GetYaxis()->GetTitle();
    TGraphVector graphs(ngraph, nullptr);
    TGraphErrorsVector egraphs(ngraph, nullptr);
    if ( m_LogLevel >= 2 ) cout << myname << "Plotting metric vs. channel. Count is " << mets.size() << endl;
    for ( Index igra=0; igra<ngraph; ++igra ) {
      string gname = hname;
      string gtitl = htitl;
      StringManipulator smanName(gname, false);
      smanName.replace("%STATUS%", statNames[igra]);
      StringManipulator smanTitl(gtitl, false);
      smanTitl.replace("%STATUS%", statNames[igra]);
      if ( useErrors ) {
        egraphs[igra] = new TGraphErrors;
        graphs[igra] = egraphs[igra];
      } else {
        graphs[igra] = new TGraph;
      }
      TGraph* pg = graphs[igra];
      pg->SetName(gname.c_str());
      pg->SetTitle(gtitl.c_str());
      if ( ! useErrors ) pg->SetMarkerStyle(2);
      pg->SetMarkerColor(statCols[igra]);
      pg->SetLineColor(statCols[igra]);
      pg->GetXaxis()->SetTitle("Channel");
      pg->GetYaxis()->SetTitle(slaby.c_str());
    }
    TGraph* pgAll = graphs[0];
    TGraphErrors* pgeAll = egraphs[0];
    double ex = 0.25;
    Index nfill = 0;
    Index icha0 = ran.begin;
    for ( MetricMap::value_type imet : mets ) {
      Index icha = imet.first;
      double met = imet.second.value;
      double err = imet.second.error;
      if ( m_LogLevel >= 3 ) {
        cout << myname << "  " << met << ": " << met << " +/- " << err << endl;
      }
      Index bin = (icha + 1) - icha0;
      ph->SetBinContent(bin, met);
      if ( err ) ph->SetBinError(bin, err);
      double gval = met;
      Index iptAll = pgAll->GetN();
      float xcha = icha + 0.5;
      pgAll->SetPoint(iptAll, xcha, gval);
      if ( pgeAll != nullptr ) pgeAll->SetPointError(iptAll, ex, err);
      if ( m_PlotUsesStatus ) {
        Index stat = channelStatus(icha);
        if ( stat > 0 ) {
          TGraph* pgStat = graphs[stat];
          TGraphErrors* pgeStat = egraphs[stat];
          Index iptStat = pgStat->GetN();
          pgStat->SetPoint(iptStat, xcha, gval);
          if ( pgeStat != nullptr ) pgeStat->SetPointError(iptStat, ex, err);
        }
      }
      ++nfill;
    }
    if ( m_LogLevel >= 3 ) cout << myname << "Filled " << nfill << " channels." << endl;
    if ( ofpname.size() && nfill > 0 ) {
      TPadManipulator man;
      if ( m_PlotSizeX && m_PlotSizeY ) man.setCanvasSize(m_PlotSizeX, m_PlotSizeY);
      //man.add(ph, "hist");
      //man.add(ph, "axis");
      man.add(pgAll, "P");
      if ( m_PlotUsesStatus ) {
        for ( int igra : {1, 3, 2} ) {
          TGraph* pgra = graphs[igra];
          if ( pgra->GetN() ) man.add(pgra, "P");
        }
      }
      man.addAxis();
      if ( m_ChannelLineModulus ) {
        for ( Index icha : m_ChannelLinePattern ) {
          man.addVerticalModLines(m_ChannelLineModulus, icha);
        }
      } else {
        for ( Index icha : m_ChannelLinePattern ) {
          if ( icha > icha0 && icha < ran.last() ) {
            man.addVerticalLine(icha, 1.0, 3);
          }
        }
      }
      man.setRangeX(ran.begin, ran.end);
      float mmin = getState().metricMin;
      float mmax = getState().metricMax;
      if ( mmax > mmin ) man.setRangeY(mmin, mmax);
      man.showGraphOverflow("BLTR", 2, statCols[0]);
      man.setGridY();
      man.print(ofpname);
    }
    if ( m_LogLevel > 1 ) {
      cout << myname << "Created plot ";
      cout << "for " << nfill << " channels in range " << ran.name << endl;
      cout << myname << " Output file: " << ofpname << endl;
    }
  }
  if ( ofrname.size() ) {
    TFile* pfile = TFile::Open(ofrname.c_str(), "UPDATE");
    ph->Write();
    if ( m_LogLevel > 1 ) cout << myname << "Wrote " << ph->GetName() << " to " << ofrname << endl;
    delete pfile;
  }
}

//**********************************************************************

TH1* AdcChannelMetric::createHisto(const AdcChannelData& acd, const IndexRange& ran) const {
  string   hname = nameReplace(    m_HistName, acd, ran);
  string   htitl = nameReplace(   m_HistTitle, acd, ran);
  string   slabm = nameReplace( m_MetricLabel, acd, ran);
  TH1* ph = nullptr;
  Index nbins = m_MetricBins;
  if ( nbins == 0 ) {
    ph = new TH1F(hname.c_str(), htitl.c_str(), ran.size(), ran.begin, ran.end);
    ph->GetXaxis()->SetTitle("Channel");
    ph->GetYaxis()->SetTitle(slabm.c_str());
  } else {
    ph = new TH1F(hname.c_str(), htitl.c_str(), nbins, getState().metricMin, getState().metricMax);
    ph->GetXaxis()->SetTitle(slabm.c_str());
    ph->GetYaxis()->SetTitle("# channels");
  }
  ph->SetDirectory(nullptr);
  ph->SetLineWidth(2);
  ph->SetStats(0);
  return ph;
}

//**********************************************************************

Index AdcChannelMetric::channelStatus(Index icha) const {
  IndexVector& stats = getState().channelStatuses;
  if ( icha >= stats.size() ) stats.resize(icha + 1, 0);
  Index& stat = stats[icha];
  if ( stat == 0 && m_pChannelStatusProvider != nullptr ) {
    if      ( m_pChannelStatusProvider->IsBad(icha) )   stat = 2;
    else if ( m_pChannelStatusProvider->IsNoisy(icha) ) stat = 3;
    else                                                stat = 1;
  }
  return stat;
}

//**********************************************************************

void AdcChannelMetric::evaluateFormulas() const {
  const string myname = "AdcChannelMetric::evaluateFormulas: ";
  if ( m_prdtool != nullptr ) {
    Index run = getState().lastRun;
    RunData rdat = m_prdtool->runData(run);
    if ( ! rdat.isValid() ) {
      cout << myname << "WARNING: RunData not found for run " << run << "." << endl;
    } else {
      rdat.setFormulaPars(*m_MetricMin);
      rdat.setFormulaPars(*m_MetricMax);
    }
  }
  float mmin = 0.0;
  float mmax = 100.0;
  if ( m_MetricMin->ready() ) {
    mmin = m_MetricMin->eval();
  } else {
    cout << myname << "WARNING: Using default metric min " << mmin << endl;
  } 
  getState().metricMin = mmin;
  if ( m_MetricMax->ready() ) {
    mmax = m_MetricMax->eval();
  } else {
    cout << myname << "WARNING: Using default metric max " << mmax << endl;
  } 
  getState().metricMax = mmax;
}

//**********************************************************************

DEFINE_ART_CLASS_TOOL(AdcChannelMetric)
