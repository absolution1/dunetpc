// AdcDataPlotter_tool.cc

#include "AdcDataPlotter.h"
#include <iostream>
#include <sstream>
#include "dune/DuneCommon/TH1Manipulator.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "TColor.h"

using std::string;
using std::cout;
using std::cin;
using std::endl;
using std::ostringstream;

using Tick = AdcSignalVector::size_type;

//**********************************************************************
// Helper.
//**********************************************************************

// Substitute xsubin for string subout in a string.
template<typename T>
void replace(string& str, string ssubout, const T& xsubin) {
  string ssubin;
  bool havesub = false;
  string::size_type ipos = str.find(ssubout);
  while ( ipos != string::npos ) {
    if ( ! havesub ) {
      ostringstream sssubin;
      sssubin << xsubin;
      ssubin = sssubin.str();
    }
    string::size_type lout = ssubout.size();
    str.replace(ipos, lout, ssubin);
    ipos = str.find(ssubout, ipos+lout);
  }
}

//**********************************************************************
// Class methods.
//**********************************************************************

AdcDataPlotter::AdcDataPlotter(fhicl::ParameterSet const& ps)
: m_LogLevel(ps.get<int>("LogLevel")), 
  m_FileName(ps.get<string>("FileName")),
  m_FirstTick(ps.get<unsigned long>("FirstTick")),
  m_LastTick(ps.get<unsigned long>("LastTick")),
  m_MaxSignal(ps.get<double>("MaxSignal")) {
  const string myname = "AdcDataPlotter::ctor: ";
  if ( m_FileName == "" ) m_FileName = "adcplot_%PAT%evt%EVENT%.png";
  if ( m_LogLevel ) {
    cout << myname << "Configuration: " << endl;
    cout << myname << "   LogLevel: " << m_LogLevel << endl;
    cout << myname << "   FileName: " << m_FileName << endl;
    cout << myname << "  FirstTick: " << m_FirstTick << endl;
    cout << myname << "   LastTick: " << m_LastTick << endl;
    cout << myname << "  MaxSignal: " << m_MaxSignal << endl;
  }
}

//**********************************************************************

int AdcDataPlotter::view(const AdcChannelDataMap& acds, string label, string fpat) const {
  const string myname = "AdcDataPlotter::view: ";
  if ( acds.size() == 0 ) {
    cout << myname << "WARNING: Channel map is empty. No plot is created." << endl;
    return 1;
  }
  const AdcChannelData& acdFirst = acds.begin()->second;
  const AdcChannelData& acdLast = acds.rbegin()->second;
  Tick maxtick = 0;
  for ( const AdcChannelDataMap::value_type& iacd : acds ) {
    if ( iacd.first == AdcChannelData::badChannel ) {
      cout << myname << "WARNING: Channel map has invalid channels. No plot is created." << endl;
    }
    Tick ntick = iacd.second.samples.size();
    if ( ntick > maxtick ) maxtick = ntick;
  }
  AdcIndex chanFirst = acdFirst.channel;
  AdcIndex chanLast = acdLast.channel;
  unsigned long tick1 = m_FirstTick;
  unsigned long tick2 = m_LastTick;
  if ( tick2 <= tick1 ) {
    tick1 = 0;
    tick2 = maxtick;
  }
  Tick ntick = tick2 - tick1;
  AdcIndex nchan = chanLast + 1 - chanFirst;
  string htitl = label;
  htitl += "; Tick; Channel";
  TH2* ph = new TH2F("hadp", htitl.c_str(), ntick, tick1, tick2, nchan, chanFirst, chanLast+1);
  ph->SetDirectory(nullptr);
  ph->SetStats(0);
  double zmax = m_MaxSignal;
  if ( zmax <= 0.0 ) zmax = 100.0;
  ph->GetZaxis()->SetRangeUser(-zmax, zmax);
  ph->SetContour(40);
  for ( const AdcChannelDataMap::value_type& iacd : acds ) {
    AdcChannel chan = iacd.first;
    const AdcSignalVector& sams = iacd.second.samples;
    unsigned int ibin = ph->GetBin(1, chan-chanFirst+1);
    for ( Tick isam=0; isam<sams.size(); ++isam, ++ibin ) {
      ph->SetBinContent(ibin, sams[isam + m_FirstTick]);
    }
  }
  if ( 1 ) {
    double alpha = 1.0;
    const int nRGBs = 8;
    Double_t stops[nRGBs] = { 0.00, 0.48, 0.50, 0.53, 0.56, 0.62, 0.80, 1.00};
    Double_t red[nRGBs]   = { 0.09, 0.75, 1.00, 1.00, 1.00, 1.00, 0.70, 0.00};
    Double_t green[nRGBs] = { 0.60, 0.80, 1.00, 1.00, 0.75, 0.55, 0.20, 0.00};
    Double_t blue[nRGBs]  = { 0.48, 0.93, 1.00, 1.00, 0.00, 0.00, 0.10, 0.00};
    //colout = 40;
    TColor::CreateGradientColorTable(nRGBs, stops, red, green, blue, 255, alpha);
  }
  TCanvas* pcan = new TCanvas;;
  pcan->SetRightMargin(0.12);
  ph->Draw("colz");
  TH1Manipulator::fixFrameFillColor();
  TH1Manipulator::addaxis(ph);
  string ofname = m_FileName;
  replace(ofname, "%PAT%", fpat);
  if ( acdFirst.event != AdcChannelData::badIndex ) replace(ofname, "%EVENT%", acdFirst.event);
  else replace(ofname, "%EVENT%", "EventNotFound");
  replace(ofname, "%CHAN1%", chanFirst);
  replace(ofname, "%CHAN2%", chanLast);
  pcan->Print(ofname.c_str());
  if ( 0 ) {
    string line;
    cout << myname;
    cout.flush();
    std::getline(cin, line);
  }
  if ( m_LogLevel > 1 ) {
    cout << myname << "Created plot ";
    if ( label.size() ) cout << label << " ";
    cout << "for channels " << acds.begin()->first << " - "
                            << acds.rbegin()->first
         << ": " << ofname << endl;
  }
  delete ph;
  delete pcan;
  return 0;
}

//**********************************************************************
