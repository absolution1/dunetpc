// AdcDataPlotter_tool.cc

#include "AdcDataPlotter.h"
#include <iostream>
#include <sstream>
#include "dune/DuneCommon/TH1Manipulator.h"
#include "dune/DuneCommon/StringManipulator.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "TColor.h"
#include "TStyle.h"
#include "TDirectory.h"
#include "TFile.h"

using std::string;
using std::cout;
using std::cin;
using std::endl;

using Tick = AdcSignalVector::size_type;
using std::vector;

//**********************************************************************
// Class methods.
//**********************************************************************

AdcDataPlotter::AdcDataPlotter(fhicl::ParameterSet const& ps)
: m_LogLevel(ps.get<int>("LogLevel")), 
  m_DataType(ps.get<int>("DataType")),
  m_FirstTick(ps.get<unsigned long>("FirstTick")),
  m_LastTick(ps.get<unsigned long>("LastTick")),
  m_MaxSignal(ps.get<double>("MaxSignal")),
  m_HistName(ps.get<string>("HistName")),
  m_PlotFileName(ps.get<string>("PlotFileName")),
  m_RootFileName(ps.get<string>("RootFileName")) {
  const string myname = "AdcDataPlotter::ctor: ";
  if ( m_LogLevel ) {
    cout << myname << "Configuration: " << endl;
    cout << myname << "      LogLevel: " << m_LogLevel << endl;
    cout << myname << "      DataType: " << m_DataType << endl;
    cout << myname << "     FirstTick: " << m_FirstTick << endl;
    cout << myname << "      LastTick: " << m_LastTick << endl;
    cout << myname << "     MaxSignal: " << m_MaxSignal << endl;
    cout << myname << "      HistName: " << m_HistName << endl;
    cout << myname << "     HistTitle: " << m_HistTitle << endl;
    cout << myname << "  PlotFileName: " << m_PlotFileName << endl;
    cout << myname << "  RootFileName: " << m_RootFileName << endl;
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
  // Create title and file names.
  string hname = m_HistName;
  string htitl = m_HistTitle;
  string ofname = m_PlotFileName;
  string ofrname = m_RootFileName;
  vector<string*> strs = {&hname, &htitl, &ofname, &ofrname};
  for ( string* pstr : strs ) {
    string& str = *pstr;
    StringManipulator sman(str);
    sman.replace("%PAT%", fpat);
    if ( acdFirst.event != AdcChannelData::badIndex ) sman.replace("%EVENT%", acdFirst.event);
    else sman.replace("%EVENT%", "EventNotFound");
    sman.replace("%CHAN1%", chanFirst);
    sman.replace("%CHAN2%", chanLast);
  }
  htitl += "; Tick; Channel";
  // Create histogram.
  TH2* ph = new TH2F(hname.c_str(), htitl.c_str(), ntick, tick1, tick2, nchan, chanFirst, chanLast+1);
  ph->SetDirectory(nullptr);
  ph->SetStats(0);
  double zmax = m_MaxSignal;
  if ( zmax <= 0.0 ) zmax = 100.0;
  ph->GetZaxis()->SetRangeUser(-zmax, zmax);
  ph->SetContour(40);
  // Fill histogram.
  bool isPrep = m_DataType == 0;
  bool isRaw = m_DataType == 1;
  for ( const AdcChannelDataMap::value_type& iacd : acds ) {
    AdcChannel chan = iacd.first;
    const AdcChannelData& acd = iacd.second;
    const AdcSignalVector& sams = acd.samples;
    const AdcCountVector& raw = acd.raw;
    AdcSignal ped = 0.0;
    bool isRawPed = false;
    if ( isRaw ) {
      ped = acd.pedestal;
      isRawPed = ped != AdcChannelData::badSignal;
    }
    unsigned int ibin = ph->GetBin(1, chan-chanFirst+1);
    for ( Tick isam=0; isam<sams.size(); ++isam, ++ibin ) {
      unsigned int isig = isam + m_FirstTick;
      if ( isPrep ) {
        if ( isig < sams.size() ) ph->SetBinContent(ibin, sams[isig]);
      } else if ( isRawPed ) {
        if ( isig < raw.size() ) ph->SetBinContent(ibin, raw[isig] - ped);
      } else {
        cout << myname << "Fill failed for bin " << ibin << endl;
      }
    }
  }
  // Save the original color map.
  vector<int> saveColors;
  if ( saveColors.size() == 0 ) {
    for ( int icol=0; icol<gStyle->GetNumberOfColors(); ++icol ) {
      saveColors.push_back(gStyle->GetColorPalette(icol));
    }
  }
  TStyle* poldstyle = gStyle;
  string styleName = "adcstyle";
  TStyle* pnewstyle = dynamic_cast<TStyle*>(gStyle->Clone(styleName.c_str()));
  pnewstyle->cd();
  if ( 1 ) {
    const double alpha = 1.0;
    const int nRGBs = 8;
    const int ncol = 255;
    Double_t stops[nRGBs] = { 0.00, 0.48, 0.50, 0.53, 0.56, 0.62, 0.80, 1.00};
    Double_t red[nRGBs]   = { 0.09, 0.75, 1.00, 1.00, 1.00, 1.00, 0.70, 0.00};
    Double_t green[nRGBs] = { 0.60, 0.80, 1.00, 1.00, 0.75, 0.55, 0.20, 0.00};
    Double_t blue[nRGBs]  = { 0.48, 0.93, 1.00, 1.00, 0.00, 0.00, 0.10, 0.00};
    static vector<int> colors;
    colors.reserve(ncol);
    // First pass we define the color table (which also sets the palette).
    // On subsequent passes, we set the palette with that table.
    if ( colors.size() == 0 ) {
      int coloff = TColor::CreateGradientColorTable(nRGBs, stops, red, green, blue, 255, alpha);
      for ( unsigned int icol=0; icol<ncol; ++icol ) colors.push_back(coloff + icol);
    } else {
      gStyle->SetPalette(ncol, &colors[0]);
    }
  }
  TCanvas* pcan = new TCanvas;;
  pcan->SetRightMargin(0.12);
  ph->Draw("colz");
  TH1Manipulator::fixFrameFillColor();
  TH1Manipulator::addaxis(ph);
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
  if ( ofrname.size() ) {
    TFile* pfile = TFile::Open(ofrname.c_str(), "UPDATE");
    ph->Write();
    if ( m_LogLevel > 1 ) cout << myname << "Wrote " << ph->GetName() << " to " << ofrname << endl;
    delete pfile;
  }
  poldstyle->cd();
  // Restore the original color map.
  gStyle->SetPalette(saveColors.size(), &saveColors[0]);
  delete pnewstyle;
  delete ph;
  delete pcan;
  return 0;
}

//**********************************************************************
