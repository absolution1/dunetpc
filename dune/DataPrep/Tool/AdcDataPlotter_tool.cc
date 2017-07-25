// AdcDataPlotter_tool.cc

#include "AdcDataPlotter.h"
#include <iostream>
#include "dune/DuneCommon/TH1Manipulator.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "TColor.h"

using std::string;
using std::cout;
using std::cin;
using std::endl;

using Tick = AdcSignalVector::size_type;

//**********************************************************************
// Class methods.
//**********************************************************************

AdcDataPlotter::AdcDataPlotter(fhicl::ParameterSet const& ps)
: m_FileName(ps.get<string>("FileName")),
  m_FirstTick(ps.get<unsigned long>("FirstTick")),
  m_LastTick(ps.get<unsigned long>("LastTick")),
  m_MaxSignal(ps.get<unsigned long>("MaxSignal")) { }

//**********************************************************************

int AdcDataPlotter::view(const AdcChannelDataMap& acds, string label, string fpat) const {
  const string myname = "AdcDataPlotter::view: ";
  if ( acds.size() == 0 ) {
    cout << myname << "WARNING: Channel map is empty. No plot is created." << endl;
    return 1;
  }
  Tick maxtick = 0;
  for ( const AdcChannelDataMap::value_type& iacd : acds ) {
    if ( iacd.first == AdcChannelData::badChannel ) {
      cout << myname << "WARNING: Channel map has invalid channels. No plot is created." << endl;
    }
    Tick ntick = iacd.second.samples.size();
    if ( ntick > maxtick ) maxtick = ntick;
  }
  AdcIndex chan1 = acds.cbegin()->first;
  AdcIndex chan2 = acds.crbegin()->first + 1;
  unsigned long tick1 = m_FirstTick;
  unsigned long tick2 = m_LastTick;
  if ( tick2 <= tick1 ) {
    tick1 = 0;
    tick2 = maxtick;
  }
  Tick ntick = tick2 - tick1;
  AdcIndex nchan = chan2 - chan1;
  string htitl = "Prepared ADC data";
  if ( label.size() ) htitl += " for " + label;
  htitl += "; Tick; Channel";
  TH2* ph = new TH2F("hadp", htitl.c_str(), ntick, tick1, tick2, nchan, chan1, chan2);
  ph->SetDirectory(nullptr);
  ph->SetStats(0);
  double zmax = m_MaxSignal;
  if ( zmax <= 0.0 ) zmax = 100.0;
  ph->GetZaxis()->SetRangeUser(-zmax, zmax);
  ph->SetContour(40);
  for ( const AdcChannelDataMap::value_type& iacd : acds ) {
    AdcChannel chan = iacd.first;
    const AdcSignalVector& sams = iacd.second.samples;
    unsigned int ibin = ph->GetBin(1, chan-chan1+1);
    for ( Tick isam=0; isam<sams.size(); ++isam, ++ibin ) {
      ph->SetBinContent(ibin, sams[isam]);
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
    cout << TColor::GetColorPalette(0) << endl;
  }
  TCanvas* pcan = new TCanvas;;
  pcan->SetRightMargin(0.12);
  ph->Draw("colz");
  TH1Manipulator::fixFrameFillColor();
  TH1Manipulator::addaxis(ph);
  string ofname = "adc-prepared";
  if ( fpat.size() ) ofname += "_" + fpat;
  ofname += ".png";
  pcan->Print(ofname.c_str());
  if ( 0 ) {
    string line;
    cout << myname;
    cout.flush();
    std::getline(cin, line);
  }
  delete ph;
  delete pcan;
  return 0;
}

//**********************************************************************
