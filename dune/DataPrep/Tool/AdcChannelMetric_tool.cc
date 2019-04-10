// AdcChannelMetric_tool.cc

#include "AdcChannelMetric.h"
#include <iostream>
#include <sstream>
#include <iomanip>
#include <set>
#include "dune/DuneInterface/Tool/AdcChannelStringTool.h"
#include "dune/DuneInterface/Tool/IndexRangeTool.h"
#include "dune/DuneCommon/TPadManipulator.h"
#include "dune/DuneCommon/StringManipulator.h"
#include "dune/DuneCommon/LineColors.h"
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
using Index = AdcChannelMetric::Index;
using TGraphVector = std::vector<TGraph*>;
using TGraphErrorsVector = std::vector<TGraphErrors*>;

//**********************************************************************
// local definitiions.
//**********************************************************************

//**********************************************************************
// Sublass methods.
//**********************************************************************

void AdcChannelMetric::AdcChannelMetric::State::update(Index run, Index event) {
  if ( callCount == 0 ) {
    firstRun = run;
    firstEvent = event;
    runCount = 1;
    eventCount = 1;
  } else {
    if ( run != lastRun ) ++runCount;
    if ( event != lastEvent ) ++eventCount;
  }
  ++callCount;
  lastEvent = event;
  lastRun = run;
}

//**********************************************************************
// Class methods.
//**********************************************************************

AdcChannelMetric::AdcChannelMetric(fhicl::ParameterSet const& ps)
: m_LogLevel(ps.get<int>("LogLevel")), 
  m_Metric(ps.get<Name>("Metric")),
  m_MetricSummaryView(ps.get<Name>("MetricSummaryView")),
  m_ChannelRanges(ps.get<NameVector>("ChannelRanges")),
  m_MetricMin(ps.get<float>("MetricMin")),
  m_MetricMax(ps.get<float>("MetricMax")),
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
  m_doSummary(false),
  m_doSummaryError(false),
  m_pChannelStatusProvider(nullptr),
  m_state(new State) {
  const string myname = "AdcChannelMetric::ctor: ";
  string stringBuilder = "adcStringBuilder";
  DuneToolManager* ptm = DuneToolManager::instance();
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
  // Fetch the naming tool.
  m_adcStringBuilder = ptm->getShared<AdcChannelStringTool>(stringBuilder);
  if ( m_adcStringBuilder == nullptr ) {
    cout << myname << "WARNING: AdcChannelStringTool not found: " << stringBuilder << endl;
  }
  // Set summary fields.
  m_doSummary = m_HistName.find("EVENT%") == string::npos;
  if ( m_doSummary ) {
    const std::set<Name> sumVals = {"count", "mean", "rms", "drms"};
    if ( m_MetricSummaryView.size() == 0 ) {
      m_MetricSummaryView = "mean:rms";
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
    cout << myname << "           MetricMin: " << m_MetricMin << endl;
    cout << myname << "           MetricMax: " << m_MetricMax << endl;
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
  for ( const MetricSummaryMap::value_type& imsm : getState().crsums ) {
    IndexRange cr = imsm.first;
    const MetricSummaryVector& msums = imsm.second;
    if ( m_LogLevel >= 3 ) {
      cout << myname << "Channel range " << cr.name << endl;
    }
    for ( Index kcha=0; kcha<cr.size(); ++ kcha ) {
      Index icha = cr.first() + kcha;
      const MetricSummary& ms = msums[kcha];
      ++ncha;
      if ( ms.count ) {
        ++nchaData;
        if ( ms.count > countMax ) {
          countMax = ms.count;
          nchaDataMax = 1;
        } else if ( ms.count == countMax ) {
          ++nchaDataMax;
        }
      }
      if ( m_LogLevel >= 3 ) {
        cout << myname << setw(8) << icha << ":" << ms.mean() << " +/- " << ms.dmean() << endl;
      }
    }
    if ( m_doSummary ) {
      MetricMap mets;
      const MetricSummaryVector& msums = getState().crsums[cr];
      for ( Index kcha=0; kcha<cr.size(); ++kcha ) {
        Index icha = cr.first() + kcha;
        const MetricSummary& msum = msums[kcha];
        if ( msum.count ) {
          Metric& met = mets[icha];
          if ( m_summaryValue.size() ) {
            met.setValue(msum.getValue(m_summaryValue));
            if ( m_summaryError.size() ) met.setError(msum.getValue(m_summaryError));
          }
        }
      }
      AdcChannelData acd;
      acd.run = getState().firstRun;
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
    cout << myname << "                          # inits: " << setw(w) << getState().initCount << endl;
    cout << myname << "                          # calls: " << setw(w) << getState().callCount << endl;
    cout << myname << "                         # events: " << setw(w) << getState().eventCount << endl;
    cout << myname << "                           # runs: " << setw(w) << getState().runCount << endl;
    cout << myname << "  Maximum # entries for a channel: " << setw(w) << countMax << endl;
    cout << myname << "       Total # channels in ranges: " << setw(w) << ncha << endl;
    cout << myname << "          # channels without data: " << setw(w) << ncha - nchaData << endl;
    cout << myname << "             # channels with data: " << setw(w) << nchaData << endl;
    cout << myname << "    # channels with max # entries: " << setw(w) << nchaDataMax << endl;
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
  float val = 0.0;
  Name sunits;
  int rstat = getMetric(acd, val, sunits);
  if ( rstat ) return ret.setStatus(rstat);
  ret.setString("metricName", m_Metric);
  ret.setFloat("metricValue", val);
  ret.setString("metricUnits", sunits);
  return ret;
}

//**********************************************************************

DataMap AdcChannelMetric::viewMap(const AdcChannelDataMap& acds) const {
  const string myname = "AdcChannelMetric::viewMap: ";
  DataMap ret;
  if ( acds.size() == 0 ) {
    cout << myname << "Input channel map is empty." << endl;
    return ret.setStatus(1);
  }
  Index chanFirst = acds.begin()->first;
  Index chanLast = acds.rbegin()->first;
  if ( m_LogLevel >= 2 ) cout << "Processing " << acds.size() << " channels: ["
                              << chanFirst << ", " << chanLast << "]" << endl;
  for ( IndexRange ran : m_crs ) {
    if ( ran.name == "all" ) {
      ran.begin = chanFirst;
      ran.end = chanLast + 1;
    }
    Index chanLo = std::max(chanFirst, ran.first());
    Index chanHi = std::min(chanLast, ran.last());
    if ( chanHi >= chanLo ) ret += viewMapForOneRange(acds, ran);
  }
  return ret;
}

//**********************************************************************

DataMap AdcChannelMetric::viewMapForOneRange(const AdcChannelDataMap& acds, const IndexRange& ran) const {
  const string myname = "AdcChannelMetric::viewMapForOneRange: ";
  DataMap ret;
  // Extract the metrics for the subset of this range included in the data.
  MetricMap mets;
  Index icha0 = ran.begin;
  AdcChannelDataMap::const_iterator iacd1=acds.lower_bound(icha0);
  AdcChannelDataMap::const_iterator iacd2=acds.upper_bound(ran.last());
  getState().update(iacd1->second.run, iacd1->second.event);
  MetricSummaryVector& metricSums = getState().crsums[ran];
  if ( metricSums.size() < ran.size() ) metricSums.resize(ran.size());
  for ( AdcChannelDataMap::const_iterator iacd=iacd1; iacd!=iacd2; ++iacd ) {
    const AdcChannelData& acd = iacd->second;
    float met;
    Name sunits;
    int rstat = getMetric(acd, met, sunits);
    if ( rstat ) {
      cout << myname << "WARNING: Metric evaluation failed for channel " << acd.channel << endl;
      continue;
    }
    Index icha = iacd->first;
    mets[icha].setValue(met);
    MetricSummary& metricSum = metricSums[icha-icha0];
    metricSum.add(met);
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
  return ret;
}

//**********************************************************************

int AdcChannelMetric::getMetric(const AdcChannelData& acd, float& val, Name& sunits) const {
  const string myname = "AdcChannelMetric::getMetric: ";
  val = 0.0;
  sunits = "";
  if ( m_Metric == "pedestal" ) {
    val = acd.pedestal;
    sunits = "ADC count";
  } else if ( m_Metric == "pedestalDiff" ) {
    float ped = acd.pedestalRms;
    val = 0.0;
    Index icha = acd.channel;
    MetricMap& pedRefs = getState().pedRefs;
    MetricMap::const_iterator ipdr = pedRefs.find(icha);
    if ( ipdr == pedRefs.end() ) {
      pedRefs[icha].value = ped;
    } else {
      val = ped - ipdr->second.value;
    }
    sunits = "ADC count";
  } else if ( m_Metric == "pedestalRms" ) {
    val = acd.pedestalRms;
    sunits = "ADC count";
  } else if ( m_Metric == "fembID" ) {
    val = acd.fembID;
  } else if ( m_Metric == "apaFembID" ) {
    val = acd.fembID%20;
  } else if ( m_Metric == "nraw" ) {
    val = acd.raw.size();
  } else if ( m_Metric == "nsam" ) {
    val = acd.samples.size();
  } else if ( m_Metric == "fembChannel" ) {
    val = acd.fembChannel;
  } else if ( m_Metric == "rawRms" ) {
    double sum = 0.0;
    double ped = acd.pedestal;
    double nsam = acd.raw.size();
    for ( AdcSignal sig : acd.raw ) {
      double dif = double(sig) - ped;
      sum += dif*dif;
    }
    val = acd.raw.size() == 0 ? 0.0 : sqrt(sum/nsam);
  } else if ( m_Metric == "rawTailFraction" ) {
    Index ntail = 0;
    double lim = 3.0*acd.pedestalRms;
    double ped = acd.pedestal;
    double nsam = acd.raw.size();
    for ( AdcSignal sig : acd.raw ) {
      double dif = double(sig) - ped;
      if ( fabs(dif) > lim ) ++ntail;
    }
    val = acd.raw.size() == 0 ? 0.0 : double(ntail)/nsam;
  } else if ( acd.hasMetadata(m_Metric) ) {
    val = acd.metadata.find(m_Metric)->second;
  // Compound metric: met1+met2
  // TODO: Move this to ctor.
  } else if ( m_Metric.find("+") != string::npos ) {
    vector<string> nams;
    string metsrem = m_Metric;
    string::size_type ipos = 0;
    val = 0.0;
    while ( ipos != string::npos ) {
      ipos = metsrem.find("+");
      string met = metsrem.substr(0, ipos);
      if ( ! acd.hasMetadata(met) ) {
        cout << myname << "ERROR: Invalid sub-metric name: " << met << endl;
        return 2;
      }
      val += acd.metadata.find(met)->second;
      if ( ipos == string::npos ) break;
      metsrem = metsrem.substr(ipos + 1);
    }
  } else {
    cout << myname << "ERROR: Invalid metric name: " << m_Metric << endl;
    return 1;
  }
  return 0;
}

//**********************************************************************

string AdcChannelMetric::
nameReplace(string name, const AdcChannelData& acd, const IndexRange& ran) const {
  StringManipulator sman(name);
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
    if ( m_LogLevel >= 2 ) cout << "Plotting # channels vs. metric. Count is " << mets.size() << endl;
    for ( MetricMap::value_type imet : mets ) {
      float met = imet.second.value;
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
    if ( m_LogLevel >= 2 ) cout << "Plotting metric vs. channel. Count is " << mets.size() << endl;
    for ( Index igra=0; igra<ngraph; ++igra ) {
      string gname = hname;
      string gtitl = htitl;
      StringManipulator smanName(gname);
      smanName.replace("%STATUS%", statNames[igra]);
      StringManipulator smanTitl(gtitl);
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
      float met = imet.second.value;
      float err = imet.second.error;
      Index bin = (icha + 1) - icha0;
      ph->SetBinContent(bin, met);
      if ( err ) ph->SetBinError(bin, err);
      float gval = met;
      if ( m_MetricMax > m_MetricMin ) {
        if ( met < m_MetricMin ) gval = m_MetricMin;
        if ( met > m_MetricMax ) gval = m_MetricMax;
      }
      Index iptAll = pgAll->GetN();
      pgAll->SetPoint(iptAll, icha, gval);
      if ( pgeAll != nullptr ) pgeAll->SetPointError(iptAll, ex, err);
      if ( m_PlotUsesStatus ) {
        Index stat = channelStatus(icha);
        if ( stat > 0 ) {
          TGraph* pgStat = graphs[stat];
          TGraphErrors* pgeStat = egraphs[stat];
          Index iptStat = pgStat->GetN();
          pgStat->SetPoint(iptStat, icha, gval);
          if ( pgeStat != nullptr ) pgeStat->SetPointError(iptStat, ex, err);
        }
      }
      ++nfill;
    }
    if ( m_LogLevel >= 3 ) cout << myname << "Filled " << nfill << " channels." << endl;
    if ( ofpname.size() ) {
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
      if ( m_MetricMax > m_MetricMin ) man.setRangeY(m_MetricMin, m_MetricMax);
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
    ph = new TH1F(hname.c_str(), htitl.c_str(), nbins, m_MetricMin, m_MetricMax);
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
