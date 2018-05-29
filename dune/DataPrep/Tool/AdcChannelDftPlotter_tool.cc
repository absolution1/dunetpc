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
#include "TH1F.h"
#include "TGraph.h"

using std::string;
using std::cout;
using std::cin;
using std::endl;
using std::vector;

//**********************************************************************
// Class methods.
//**********************************************************************

AdcChannelDftPlotter::AdcChannelDftPlotter(fhicl::ParameterSet const& ps)
: m_LogLevel(ps.get<int>("LogLevel")), 
  m_Variable(ps.get<Name>("Variable")),
  m_HistName(ps.get<Name>("HistName")),
  m_HistTitle(ps.get<Name>("HistTitle")),
  m_PlotName(ps.get<Name>("PlotName")) {
  const string myname = "AdcChannelDftPlotter::ctor: ";
  // Fetch optional vlaues and check Variable.
  if ( m_Variable == "power" ) {
    m_SampleFreq = ps.get<float>("SampleFreq");
    m_Rebin = ps.get<Index>("Rebin");
  } else if ( m_Variable == "magnitude" ) {
  } else if ( m_Variable == "phase" ) {
  } else {
    cout << myname << "Invalid Variable: " << m_Variable << endl;
    throw art::Exception(art::errors::Configuration, "InvalidFclValue");
  }
  // Fetch renaming tools.
  string snameBuilder = "adcNameBuilder";
  DuneToolManager* ptm = DuneToolManager::instance();
  m_adcNameBuilder = ptm->getShared<AdcChannelStringTool>(snameBuilder);
  if ( m_adcNameBuilder == nullptr ) {
    cout << myname << "WARNING: AdcChannelStringTool not found: " << snameBuilder << endl;
  }
  string stitlBuilder = "adcTitleBuilder";
  m_adcTitleBuilder = ptm->getShared<AdcChannelStringTool>(stitlBuilder);
  if ( m_adcTitleBuilder == nullptr ) {
    cout << myname << "WARNING: AdcChannelStringTool not found: " << stitlBuilder << endl;
  }
  if ( m_LogLevel ) {
    cout << myname << "Configuration: " << endl;
    cout << myname << "         LogLevel: " << m_LogLevel << endl;
    cout << myname << "         Variable: " << m_Variable << endl;
    cout << myname << "       SampleFreq: " << m_SampleFreq << endl;
    cout << myname << "            Rebin: " << m_Rebin << endl;
    cout << myname << "         HistName: " << m_HistName << endl;
    cout << myname << "        HistTitle: " << m_HistTitle << endl;
    cout << myname << "         PlotName: " << m_PlotName << endl;
  }
}

//**********************************************************************

DataMap AdcChannelDftPlotter::view(const AdcChannelData& acd) const {
  const string myname = "AdcChannelDftPlotter::view: ";
  DataMap ret;
  bool doMag = m_Variable == "magnitude";
  bool doPha = m_Variable == "phase";
  bool doPwr = m_Variable == "power";
  if ( ! doMag && !doPha && !doPwr ) {
    cout << myname << "Invalid plot variable: " << m_Variable << endl;
    return ret.setStatus(1);
  }
  Index nmag = acd.dftmags.size();
  Index npha = acd.dftphases.size();
  if ( nmag == 0 ) {
    cout << myname << "DFT is not present." << endl;
    return ret.setStatus(2);
  }
  if ( npha > nmag || nmag - npha > 1 ) {
    cout << myname << "DFT is not valid." << endl;
    return ret.setStatus(3);
  }
  TObject* pobj = nullptr;
  //string hname = nameReplace(m_HistName,  acd, false);
  string hname = AdcChannelStringTool::build(m_adcNameBuilder, acd, m_HistName);
  string htitl = nameReplace(m_HistTitle, acd, true);
  string pname = nameReplace(m_PlotName,  acd, false);
  float pi = acos(-1.0);
  double xmin = 0.0;
  double xmax = 0.0;
  double ymin = 0.0;
  double ymax = 0.0;
  string dopt;
  if ( doMag || doPha ) {  
    string ytitl = "Phase";
    if ( doMag ) {
      ytitl = AdcChannelStringTool::build(m_adcTitleBuilder, acd, "Amplitude% [SUNIT]%");
    }
    TGraph* pg = new TGraph;
    pg->SetName(hname.c_str());
    pg->SetTitle(htitl.c_str());
    pg->SetMarkerStyle(2);
    pg->SetMarkerColor(602);
    //ph = new TH1F(hname.c_str(), htitl.c_str(), nbin, 0, nbin);
    //ph->SetDirectory(nullptr);
    //ph->SetLineWidth(2);
    pg->GetXaxis()->SetTitle("Frequency index");
    pg->GetYaxis()->SetTitle(ytitl.c_str());
    float magMax = 0.0;
    for ( Index ipha=0; ipha<nmag; ++ipha ) {
      float mag = acd.dftmags[ipha];
      float pha = ipha < npha ? acd.dftphases[ipha] : 0.0;
      if ( mag < 0 ) {
        pha += ( pha > 0 ? -pi : pi);
        mag = -mag;
      }
      if ( mag > magMax ) magMax = mag;
      if ( doPha ) pg->SetPoint(ipha, ipha, pha);
      else pg->SetPoint(ipha, ipha, mag);
    }
    float delx = 0.02*nmag;
    xmin = -delx;
    xmax = nmag -1 + delx;
    ymax =  3.2;
    ymin = -ymax;
    if ( doMag ) {
      ymin = 0.0;
      ymax = magMax*1.03;
    }
    pobj = pg;
    dopt = "P";
  }
  if ( pname.size() ) {
    TPadManipulator man;
    //if ( m_PlotSizeX && m_PlotSizeY ) man.setCanvasSize(m_PlotSizeX, m_PlotSizeY);
    man.add(pobj, dopt);
    man.addAxis();
    if ( xmax > xmin ) man.setRangeX(xmin, xmax);
    if ( ymax > ymin ) man.setRangeY(ymin, ymax);
    man.print(pname);
  }
  return ret;
}

//**********************************************************************

string AdcChannelDftPlotter::
nameReplace(string name, const AdcChannelData& acd, bool isTitle) const {
  const AdcChannelStringTool* pnbl = nullptr;
  if ( isTitle ) pnbl = m_adcTitleBuilder;
  else {
    pnbl = m_adcNameBuilder == nullptr ? m_adcTitleBuilder : m_adcNameBuilder;
  }
  if ( pnbl == nullptr ) return name;
  DataMap dm;
  return pnbl->build(acd, dm, name);
}

//**********************************************************************
