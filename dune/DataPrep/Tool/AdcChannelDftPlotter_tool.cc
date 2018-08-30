// AdcChannelDftPlotter_tool.cc

#include "AdcChannelDftPlotter.h"
#include "dune/DuneInterface/Tool/AdcChannelStringTool.h"
#include "dune/DuneCommon/TPadManipulator.h"
#include "dune/DuneCommon/StringManipulator.h"
#include "dune/ArtSupport/DuneToolManager.h"
#include "canvas/Utilities/Exception.h"
#include <iostream>
#include <sstream>
#include <vector>
#include <iomanip>
#include "TH1F.h"
#include "TGraph.h"

using std::string;
using std::cout;
using std::cin;
using std::endl;
using std::vector;
using std::ostringstream;
using std::setw;
using std::fixed;

//**********************************************************************
// Class methods.
//**********************************************************************

AdcChannelDftPlotter::AdcChannelDftPlotter(fhicl::ParameterSet const& ps)
: AdcMultiChannelPlotter(ps, "Plot"),
  m_LogLevel(ps.get<int>("LogLevel")), 
  m_Variable(ps.get<Name>("Variable")),
  m_SampleFreq(ps.get<float>("SampleFreq")),
  m_YMax(0.0),
  m_YMinLog(ps.get<float>("YMinLog")),
  m_NBinX(0),
  m_HistName(ps.get<Name>("HistName")),
  m_HistTitle(ps.get<Name>("HistTitle")),
  m_PlotName(ps.get<Name>("PlotName")) {
  const string myname = "AdcChannelDftPlotter::ctor: ";
  bool doMag = m_Variable == "magnitude";
  bool doPha = m_Variable == "phase";
  bool doPwr = m_Variable == "power";
  bool doPwt = m_Variable == "power/tick";
  // Check variable and get optional fields.
  if ( doPwr || doPwt ) {
    m_YMax = ps.get<float>("YMax");
    m_NBinX = ps.get<Index>("NBinX");
  } else if ( doMag ) {
    m_YMax = ps.get<float>("YMax");
  } else if ( doPha ) {
  } else {
    cout << myname << "Invalid Variable: " << m_Variable << endl;
    throw art::Exception(art::errors::Configuration, "InvalidFclValue");
  }
  // Fetch renaming tools.
  string snameBuilder = "adcStringBuilder";
  DuneToolManager* ptm = DuneToolManager::instance();
  m_adcStringBuilder = ptm->getShared<AdcChannelStringTool>(snameBuilder);
  if ( m_adcStringBuilder == nullptr ) {
    cout << myname << "WARNING: AdcChannelStringTool not found: " << snameBuilder << endl;
  }
  // Display the configuration.
  if ( m_LogLevel ) {
    cout << myname << "Configuration: " << endl;
    cout << myname << "         LogLevel: " << m_LogLevel << endl;
    cout << myname << "         Variable: " << m_Variable << endl;
    cout << myname << "       SampleFreq: " << m_SampleFreq << endl;
    if ( doMag || doPwr || doPwt ) cout << myname << "             YMax: " << m_YMax << endl;
    if ( doPwr || doPwt )          cout << myname << "            NBinX: " << m_NBinX << endl;
    cout << myname << "          YMinLog: " << m_YMinLog << endl;
    cout << myname << "         HistName: " << m_HistName << endl;
    cout << myname << "        HistTitle: " << m_HistTitle << endl;
    cout << myname << "         PlotName: " << m_PlotName << endl;
  }
}


//**********************************************************************

int AdcChannelDftPlotter::
viewMapChannel(const AdcChannelData& acd, DataMap&, TPadManipulator& man) const {
  const string myname = "AdcChannelDftPlotter::viewMapChannel: ";
  DataMap chret = viewLocal(acd);
  fillChannelPad(chret, man);
  return 0;
}

//**********************************************************************

DataMap AdcChannelDftPlotter::view(const AdcChannelData& acd) const {
  const string myname = "AdcChannelDftPlotter::view: ";
  DataMap chret = viewLocal(acd);
  if ( getPlotName().size() ) {
    string pname = AdcChannelStringTool::build(m_adcStringBuilder, acd, getPlotName());
    TPadManipulator man;
    fillChannelPad(chret, man);
    man.print(pname);
  }
  return chret;
}

//**********************************************************************

DataMap AdcChannelDftPlotter::viewLocal(const AdcChannelData& acd) const {
  const string myname = "AdcChannelDftPlotter::viewLocal: ";
  DataMap ret;
  bool doMag = m_Variable == "magnitude";
  bool doPha = m_Variable == "phase";
  bool doPwr = m_Variable == "power";
  bool doPwt = m_Variable == "power/tick";
  bool haveFreq = m_SampleFreq > 0.0;
  if ( ! doMag && !doPha && !doPwr && !doPwt ) {
    cout << myname << "Invalid plot variable: " << m_Variable << endl;
    return ret.setStatus(1);
  }
  Index nmag = acd.dftmags.size();
  Index npha = acd.dftphases.size();
  Index nsam = nmag + npha - 1;
  if ( nmag == 0 ) {
    cout << myname << "DFT is not present." << endl;
    return ret.setStatus(2);
  }
  if ( npha > nmag || nmag - npha > 1 ) {
    cout << myname << "DFT is not valid." << endl;
    return ret.setStatus(3);
  }
  string hname = AdcChannelStringTool::build(m_adcStringBuilder, acd, m_HistName);
  string htitl = AdcChannelStringTool::build(m_adcStringBuilder, acd, m_HistTitle);
  float pi = acos(-1.0);
  double xFac = haveFreq ? m_SampleFreq/nsam : 1.0;
  double xmin = 0.0;
  double xmax = (nmag-1)*xFac;
  float yValMax = 0.0;
  string dopt;
  string xtitl = haveFreq ? "Frequency [kHz]" : "Frequency index";
  if ( doMag || doPha ) {  
    string ytitl = "Phase";
    if ( doMag ) {
      ytitl = AdcChannelStringTool::build(m_adcStringBuilder, acd, "Amplitude% [SUNIT]%");
    }
    TGraph* pg = new TGraph;
    pg->SetName(hname.c_str());
    pg->SetTitle(htitl.c_str());
    pg->SetMarkerStyle(2);
    pg->SetMarkerColor(602);
    //ph = new TH1F(hname.c_str(), htitl.c_str(), nbin, 0, nbin);
    //ph->SetDirectory(nullptr);
    //ph->SetLineWidth(2);
    string xtitl = haveFreq ? "Frequency [kHz]" : "Frequency index";
    pg->GetXaxis()->SetTitle(xtitl.c_str());
    pg->GetYaxis()->SetTitle(ytitl.c_str());
    for ( Index ipha=0; ipha<nmag; ++ipha ) {
      float mag = acd.dftmags[ipha];
      float pha = ipha < npha ? acd.dftphases[ipha] : 0.0;
      if ( mag < 0 ) {
        pha += ( pha > 0 ? -pi : pi);
        mag = -mag;
      }
      float x = ipha*xFac;
      float y = doPha ? pha : mag;
      if ( y > yValMax ) yValMax = y;
      pg->SetPoint(ipha, x, y);
    }
    ret.setGraph("dftGraph", pg);
    ret.setString("dftDopt", "P");
  } else if ( doPwr || doPwt ) {
    string ytitl = doPwr ? "Power" : "Power/tick";
    if ( acd.sampleUnit.size() ) {
      ytitl = AdcChannelStringTool::build(m_adcStringBuilder, acd, ytitl + " [(%SUNIT%)^{2}]");
    }
    if ( m_NBinX == 0 ) {
      cout << myname << "Invalid bin count: " << m_NBinX << endl;
      return ret.setStatus(2);
    }
    // Shift bins sightly so f=0 is an underflow and last frequency is not an overflow.
    double delx = 1.e-5*xmax;
    double x1 = xmin + delx;
    double x2 = xmax + delx;
    xmax = xmin;
    TH1* ph = new TH1F(hname.c_str(), htitl.c_str(), m_NBinX, x1, x2);
    ph->SetDirectory(nullptr);
    ph->SetLineWidth(2);
    ph->GetXaxis()->SetTitle(xtitl.c_str());
    ph->GetYaxis()->SetTitle(ytitl.c_str());
    float pwrFac = doPwr ? 1.0 : 1.0/nsam;
    for ( Index ipha=0; ipha<nmag; ++ipha ) {
      float mag = acd.dftmags[ipha];
      float x = ipha*xFac;
      float y = pwrFac*mag*mag;
      ph->Fill(x, y);
    }
    if ( ph->GetBinContent(m_NBinX+1) ) {
      cout << myname << "WARNING: Power histogram has overflow." << endl;
    }
    for ( Index ibin=0; ibin<m_NBinX+2; ++ibin ) {
      double y = ph->GetBinContent(ibin);
      if ( y > yValMax ) yValMax = y;
    }
    dopt = "hist";
    ret.setHist("dftHist", ph, true);
    ret.setString("dftDopt", "hist");
  }
  ret.setFloat("dftYValMax", yValMax);
  return ret;
}

//**********************************************************************

int AdcChannelDftPlotter::fillChannelPad(DataMap& dm, TPadManipulator& man) const {
  const string myname = "AdcChannelDftPlotter::fillChannelPad: ";
  TGraph* pg = dm.getGraph("dftGraph");
  TH1* ph = dm.getHist("dftHist");
  float yValMax = dm.getFloat("dftYValMax");
  bool doPha = m_Variable == "phase";
  bool doPwr = m_Variable == "power";
  bool doPwt = m_Variable == "power/tick";
  string dopt = dm.getString("dftDopt");
  bool logy = false;
  // Assign y limits.
  double ymin = 0.0;
  double ymax = 0.0;
  if ( doPha ) {
    ymax =  3.2;
    ymin = -ymax;
  } else {
    ymin = 0.0;
    if ( m_YMax > 0 ) ymax = m_YMax;
    else if ( m_YMax < 0 && -m_YMax > yValMax ) ymax = -m_YMax;
    else ymax = yValMax*1.02;
  }
  if ( m_YMinLog ) {
    ymin = m_YMinLog;
    logy = true;
  }
  double xmin = 0.0;
  double xmax = 0.0;
  if ( pg != nullptr ) {
    xmax = pg->GetXaxis()->GetXmax();
    xmin = -0.02*xmax;
    xmax *= 1.02;
    man.add(pg, dopt);
  } else if ( ph != nullptr ) {
    man.add(ph, dopt);
  } else {
    cout << myname << "ERROR: Neither hist or graph is defined." << endl;
    return 1;
  }
  man.addAxis();
  if ( xmax > xmin ) man.setRangeX(xmin, xmax);
  if ( ymax > ymin ) man.setRangeY(ymin, ymax);
  if ( doPwr || doPwt ) man.showUnderflow();
  if ( logy ) man.setLogY();
  if ( logy ) man.setGridY();
  if ( doPwt ) {
    ostringstream ssout;
    ssout.precision(2);
    double sum = ph->Integral(0, ph->GetNbinsX()+1);
    ssout << "#sqrt{#Sigma} = " << fixed << setw(2) << sqrt(sum);
    TLatex* ptxt = new TLatex(0.75, 0.80, ssout.str().c_str());
    ptxt->SetNDC();
    ptxt->SetTextFont(42);
    man.add(ptxt);
  }
  return 0;
}

//**********************************************************************
