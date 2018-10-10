// AdcChannelMetric_tool.cc

#include "AdcChannelMetric.h"
#include <iostream>
#include <sstream>
#include "dune/DuneInterface/Tool/AdcChannelStringTool.h"
#include "dune/DuneInterface/Tool/IndexRangeTool.h"
#include "dune/DuneCommon/TPadManipulator.h"
#include "dune/DuneCommon/StringManipulator.h"
#include "dune/ArtSupport/DuneToolManager.h"
#include "TH1F.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TColor.h"
#include "TStyle.h"
#include "TDirectory.h"
#include "TFile.h"

using std::string;
using std::cout;
using std::cin;
using std::endl;
using std::vector;

//**********************************************************************
// local definitiions.
//**********************************************************************

//**********************************************************************
// Class methods.
//**********************************************************************

AdcChannelMetric::AdcChannelMetric(fhicl::ParameterSet const& ps)
: m_LogLevel(ps.get<int>("LogLevel")), 
  m_Metric(ps.get<Name>("Metric")),
  m_ChannelRanges(ps.get<NameVector>("ChannelRanges")),
  m_MetricMin(ps.get<float>("MetricMin")),
  m_MetricMax(ps.get<float>("MetricMax")),
  m_ChannelLineModulus(ps.get<Index>("ChannelLineModulus")),
  m_ChannelLinePattern(ps.get<IndexVector>("ChannelLinePattern")),
  m_HistName(ps.get<Name>("HistName")),
  m_HistTitle(ps.get<Name>("HistTitle")),
  m_MetricLabel(ps.get<Name>("MetricLabel")),
  m_PlotSizeX(ps.get<Index>("PlotSizeX")),
  m_PlotSizeY(ps.get<Index>("PlotSizeY")),
  m_PlotFileName(ps.get<Name>("PlotFileName")),
  m_RootFileName(ps.get<Name>("RootFileName")) {
  const string myname = "AdcChannelMetric::ctor: ";
  string stringBuilder = "adcStringBuilder";
  DuneToolManager* ptm = DuneToolManager::instance();
  // Fetch the channel ranges.
  bool toolNotFound = false;
  const IndexRangeTool* pcrt = nullptr;
  for ( Name crn : m_ChannelRanges.size() ? m_ChannelRanges : NameVector(1, "") ) {
    if ( crn.size() == 0 || crn == "all" ) {
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
  // Fetch the naming tool.
  m_adcStringBuilder = ptm->getShared<AdcChannelStringTool>(stringBuilder);
  if ( m_adcStringBuilder == nullptr ) {
    cout << myname << "WARNING: AdcChannelStringTool not found: " << stringBuilder << endl;
  }
  // Display the configuration.
  if ( m_LogLevel ) {
    cout << myname << "Configuration: " << endl;
    cout << myname << "            LogLevel: " << m_LogLevel << endl;
    cout << myname << "              Metric: " << m_Metric << endl;
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
    cout << myname << "        RootFileName: " << m_RootFileName << endl;
  }
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
  // At this point, there is only one range to plot.
  const AdcChannelData& acdFirst = acds.begin()->second;
  string   hname = nameReplace(    m_HistName, acdFirst, ran);
  string   htitl = nameReplace(   m_HistTitle, acdFirst, ran);
  string   slaby = nameReplace( m_MetricLabel, acdFirst, ran);
  string  ofname = nameReplace(m_PlotFileName, acdFirst, ran);
  string ofrname = nameReplace(m_RootFileName, acdFirst, ran);
  Index nchan = ran.size();
  Index ich0 = ran.begin;
  TH1* ph = new TH1F(hname.c_str(), htitl.c_str(), nchan, ich0, ran.end);
  ph->SetDirectory(nullptr);
  ph->SetLineWidth(2);
  ph->SetStats(0);
  ph->GetXaxis()->SetTitle("Channel");
  ph->GetYaxis()->SetTitle(slaby.c_str());
  TGraph* pg = new TGraph;
  pg->SetName(hname.c_str());
  pg->SetTitle(htitl.c_str());
  pg->SetMarkerStyle(2);
  pg->SetMarkerColor(602);
  pg->GetXaxis()->SetTitle("Channel");
  pg->GetYaxis()->SetTitle(slaby.c_str());
  Index nfill = 0;
  float val = 0.0;
  Name sunits;
  AdcChannelDataMap::const_iterator iacd1=acds.lower_bound(ich0);
  AdcChannelDataMap::const_iterator iacd2=acds.upper_bound(ran.last());
  for ( AdcChannelDataMap::const_iterator iacd=iacd1; iacd!=iacd2; ++iacd ) {
    const AdcChannelData& acd = iacd->second;
    int rstat = getMetric(acd, val, sunits);
    if ( rstat ) {
      cout << myname << "WARNING: Metric evaluation failed for channel " << acd.channel << endl;
      continue;
    }
    Index icha = iacd->first;
    Index bin = (icha + 1) - ich0;
    ph->SetBinContent(bin, val);
    float gval = val;
    if ( m_MetricMax > m_MetricMin ) {
      if ( val < m_MetricMin ) gval = m_MetricMin;
      if ( val > m_MetricMax ) gval = m_MetricMax;
    }
    pg->SetPoint(pg->GetN(), icha, gval);
    ++nfill;
  }
  if ( m_LogLevel >= 3 ) cout << myname << "Filled " << nfill << " channels." << endl;
  if ( ofname.size() ) {
    TPadManipulator man;
    if ( m_PlotSizeX && m_PlotSizeY ) man.setCanvasSize(m_PlotSizeX, m_PlotSizeY);
    //man.add(ph, "hist");
    //man.add(ph, "axis");
    man.add(pg, "P");
    man.addAxis();
    if ( m_ChannelLineModulus ) {
      for ( Index icha : m_ChannelLinePattern ) {
        man.addVerticalModLines(m_ChannelLineModulus, icha);
      }
    } else {
      for ( Index icha : m_ChannelLinePattern ) {
        if ( icha > ich0 && icha < ran.last() ) {
          man.addVerticalLine(icha, 1.0, 3);
        }
      }
    }
    man.setRangeX(ran.begin, ran.end);
    if ( m_MetricMax > m_MetricMin ) man.setRangeY(m_MetricMin, m_MetricMax);
    man.setGridY();
    man.print(ofname);
  }
  if ( m_LogLevel > 1 ) {
    cout << myname << "Created plot ";
    cout << "for " << nfill << " channels in range " << ran.name << endl;
    cout << myname << " Output file: " << ofname << endl;
  }
  if ( ofrname.size() ) {
    TFile* pfile = TFile::Open(ofrname.c_str(), "UPDATE");
    ph->Write();
    if ( m_LogLevel > 1 ) cout << myname << "Wrote " << ph->GetName() << " to " << ofrname << endl;
    delete pfile;
  }
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
    sunits = "ADC counts";
  } else if ( m_Metric == "pedestalRms" ) {
    val = acd.pedestalRms;
    sunits = "ADC counts";
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
