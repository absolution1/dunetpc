// AdcChannelMetric_tool.cc

#include "AdcChannelMetric.h"
#include <iostream>
#include <sstream>
#include "dune/DuneCommon/TPadManipulator.h"
#include "dune/DuneCommon/StringManipulator.h"
#include "TH1F.h"
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
  m_FirstChannel(ps.get<unsigned int>("FirstChannel")),
  m_LastChannel(ps.get<unsigned int>("LastChannel")),
  m_ChannelCounts(ps.get<IndexVector>("ChannelCounts")),
  m_MetricMin(ps.get<float>("MetricMin")),
  m_MetricMax(ps.get<float>("MetricMax")),
  m_ChannelLineModulus(ps.get<Index>("ChannelLineModulus")),
  m_ChannelLinePattern(ps.get<IndexVector>("ChannelLinePattern")),
  m_HistName(ps.get<Name>("HistName")),
  m_HistTitle(ps.get<Name>("HistTitle")),
  m_PlotSizeX(ps.get<Index>("PlotSizeX")),
  m_PlotSizeY(ps.get<Index>("PlotSizeY")),
  m_PlotFileName(ps.get<Name>("PlotFileName")),
  m_RootFileName(ps.get<Name>("RootFileName"))
  {

  const string myname = "AdcChannelMetric::ctor: ";
  if ( m_LogLevel ) {
    cout << myname << "Configuration: " << endl;
    cout << myname << "            LogLevel: " << m_LogLevel << endl;
    cout << myname << "              Metric: " << m_Metric << endl;
    cout << myname << "        FirstChannel: " << m_FirstChannel << endl;
    cout << myname << "         LastChannel: " << m_LastChannel << endl;
    cout << myname << "           MetricMin: " << m_MetricMin << endl;
    cout << myname << "           MetricMax: " << m_MetricMax << endl;
    cout << myname << "  ChannelLineModulus: " << m_ChannelLineModulus << endl;
    cout << myname << "  ChannelLinePattern: {";
    bool first = true;
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
  const AdcChannelData& acdFirst = acds.begin()->second;
  //const AdcChannelData& acdLast = acds.rbegin()->second;
  Index acdChannelFirst = acds.begin()->first;
  Index acdChannelLast = acds.rbegin()->first;
  AdcIndex chanFirst = m_FirstChannel;
  AdcIndex chanLast = m_LastChannel;
  if ( chanLast <= chanFirst ) {
    chanFirst = acdChannelFirst;
    chanLast = acdChannelLast + 1;
  }
  if ( m_LogLevel >= 2 ) cout << "Processing " << acds.size() << " channels in the range ("
                              << chanFirst << ", " << chanLast << ")" << endl;
  // If channel limits are defined, then we process each channel range separately.
  Index ncl = m_ChannelCounts.size();
  if ( ncl ) {
    Index icl = 0;
    Index ich1 = 0;
    while ( ich1 <= acdChannelLast ) {
      Index ich2 = ich1 + m_ChannelCounts[icl];
      icl = (icl + 1) % ncl;
      AdcChannelMetric acd(*this);
      acd.m_ChannelCounts.clear();
      acd.m_FirstChannel = ich1;
      acd.m_LastChannel = ich2;
      ret = acd.viewMap(acds);
      ich1 = ich2;
    }
    return ret;
  }
  // Create title and file names.
  string hname = m_HistName;
  string htitl = m_HistTitle;
  string ofname = m_PlotFileName;
  string ofrname = m_RootFileName;
  vector<string*> strs = {&hname, &htitl, &ofname, &ofrname};
  for ( string* pstr : strs ) {
    string& str = *pstr;
    StringManipulator sman(str);
    if ( acdFirst.run != AdcChannelData::badIndex ) sman.replace("%RUN%", acdFirst.run);
    else sman.replace("%RUN%", "RunNotFound");
    if ( acdFirst.subRun != AdcChannelData::badIndex ) sman.replace("%SUBRUN%", acdFirst.subRun);
    else sman.replace("%SUBRUN%", "SubRunNotFound");
    if ( acdFirst.event != AdcChannelData::badIndex ) sman.replace("%EVENT%", acdFirst.event);
    else sman.replace("%EVENT%", "EventNotFound");
    sman.replace("%CHAN1%", chanFirst);
    sman.replace("%CHAN2%", chanLast);
  }
  Index nchan = chanLast - chanFirst;
  string sttl;
  string slaby;
  if ( m_Metric == "pedestal" ) {
    sttl = "ADC pedestal";
    slaby = "Pedestal";
  } else if ( m_Metric == "pedestalRms" ) {
    sttl = "ADC pedestal RMS";
    slaby = "Pedestal RMS";
  } else {
    cout << myname << "ERROR: Invalid metric name: " << m_Metric << endl;
    return ret.setStatus(1);
  }
  if ( htitl == "AUTO" ) htitl = sttl;
  TH1* ph = new TH1F(hname.c_str(), htitl.c_str(), nchan, chanFirst, chanLast);
  ph->SetDirectory(nullptr);
  ph->SetLineWidth(2);
  ph->SetStats(0);
  Index nfill = 0;
  AdcChannelDataMap::const_iterator iacd1=acds.lower_bound(chanFirst);
  AdcChannelDataMap::const_iterator iacd2=acds.upper_bound(chanLast);
  float val = 0.0;
  Name sunits;
  for ( AdcChannelDataMap::const_iterator iacd=iacd1; iacd!=iacd2; ++iacd ) {
    const AdcChannelData& acd = iacd->second;
    int rstat = getMetric(acd, val, sunits);
    if ( rstat ) {
      cout << myname << "WARNING: Metric evaluation failed for channel " << acd.channel << endl;
      continue;
    }
    Index icha = iacd->first;
    Index bin = icha - chanFirst + 1;
    ph->SetBinContent(bin, val);
    if ( nfill == 0 ) {
      if ( sunits.size() ) {
        slaby = slaby + " [" + sunits + "]";
      }
      ph->GetXaxis()->SetTitle("Channel");
      ph->GetYaxis()->SetTitle(slaby.c_str());
    }
    ++nfill;
  }
  if ( m_LogLevel >= 3 ) cout << myname << "Filled " << nfill << " channels." << endl;
  if ( ofname.size() ) {
    TPadManipulator man;
    if ( m_PlotSizeX && m_PlotSizeY ) man.setCanvasSize(m_PlotSizeX, m_PlotSizeY);
    man.add(ph, "hist");
    man.addAxis();
    if ( m_ChannelLineModulus ) {
      for ( Index icha : m_ChannelLinePattern ) {
        man.addVerticalModLines(m_ChannelLineModulus, icha);
      }
    } else {
      for ( Index icha : m_ChannelLinePattern ) {
        if ( icha > chanFirst && icha < chanLast ) {
          man.addVerticalLine(icha, 1.0, 3);
        }
      }
    }
    if ( m_MetricMax > m_MetricMin ) man.setRangeY(m_MetricMin, m_MetricMax);
    man.print(ofname);
  }
  if ( m_LogLevel > 1 ) {
    cout << myname << "Created plot ";
    cout << "for " << nfill << " channels in range " << chanFirst << " - " << chanLast << endl;
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
  } else {
    cout << myname << "ERROR: Invalid metric name: " << m_Metric << endl;
    return 1;
  }
  return 0;
}

//**********************************************************************
