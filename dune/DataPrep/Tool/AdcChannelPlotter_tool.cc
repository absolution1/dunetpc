// AdcChannelPlotter_tool.cc

#include "AdcChannelPlotter.h"
#include "dune/DuneCommon/StringManipulator.h"
#include "dune/DuneInterface/Tool/AdcChannelStringTool.h"
#include "dune/DuneCommon/TPadManipulator.h"
#include "dune/DuneInterface/Tool/HistogramManager.h"
#include "dune/ArtSupport/DuneToolManager.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include "TDirectory.h"
#include "TFile.h"
#include "TH1F.h"

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

using Index = unsigned int;

//**********************************************************************
// Class methods.
//**********************************************************************

AdcChannelPlotter::AdcChannelPlotter(fhicl::ParameterSet const& ps)
: m_LogLevel(ps.get<int>("LogLevel")),
  m_HistTypes(ps.get<NameVector>("HistTypes")),
  m_HistName(ps.get<string>("HistName")),
  m_HistTitle(ps.get<string>("HistTitle")),
  m_RootFileName(ps.get<string>("RootFileName")),
  m_PlotFileName(ps.get<string>("PlotFileName")),
  m_HistManager(ps.get<string>("HistManager")),
  m_phm(nullptr) {
  const string myname = "AdcChannelPlotter::ctor: ";
  if ( m_HistTypes.size() == 0 ) {
    cout << myname << "WARNING: No histogram types are specified." << endl;
    return;
  }
  DuneToolManager* ptm = DuneToolManager::instance();
  if ( m_HistManager.size() ) {
    m_phm = ptm->getShared<HistogramManager>(m_HistManager);
    if ( m_phm == nullptr ) {
      cout << myname << "WARNING: Histoggram manager not found: " << m_HistManager << endl;
    }
  }
  string stringBuilder = "adcStringBuilder";
  m_adcStringBuilder = ptm->getShared<AdcChannelStringTool>(stringBuilder);
  if ( m_adcStringBuilder == nullptr ) {
    cout << myname << "WARNING: AdcChannelStringTool not found: " << stringBuilder << endl;
  }
  if ( m_LogLevel > 0 ) {
    cout << myname << "      LogLevel: " << m_LogLevel << endl;
    cout << myname << "     HistTypes: [";
    bool first = true;
    for ( string name : m_HistTypes ) {
      if ( ! first ) cout << ", ";
      first = false;
      cout << name;
    }
    cout << "]" << endl;
    cout << myname << "      HistName: " << m_HistName << endl;
    cout << myname << "     HistTitle: " << m_HistTitle << endl;
    cout << myname << "  RootFileName: " << m_RootFileName << endl;
    cout << myname << "  PlotFileName: " << m_PlotFileName << endl;
    cout << myname << "   HistManager: " << m_HistManager << endl;
  }
}

//**********************************************************************

DataMap AdcChannelPlotter::view(const AdcChannelData& acd) const {
  const string myname = "AdcChannelPlotter::view: ";
  DataMap res;
  if ( m_HistTypes.size() == 0 ) {
    cout << myname << "WARNING: No histogram types are specified." << endl;
    return res.setStatus(1);
  }
  string hnameBase = m_HistName;
  if ( hnameBase == "" ) hnameBase = "%TYPE%";
  string htitlBase = m_HistTitle;
  if ( htitlBase == "" ) htitlBase = "%TYPE%";
  TDirectory* polddir = gDirectory;
  TFile* pfile = nullptr;
  if ( m_RootFileName.size() ) {
    string fname = nameReplace(m_RootFileName, acd, m_HistTypes[0]);
    pfile = TFile::Open(fname.c_str(), "UPDATE");
  }
  vector<TH1*> hists;
  for ( string type : m_HistTypes ) {
    TH1* ph = nullptr;
    string hname = nameReplace(hnameBase, acd, type);
    string htitl = nameReplace(htitlBase, acd, type);
    if ( type == "raw" ) {
      Index nsam = acd.raw.size();
      if ( nsam == 0 ) {
        cout << myname << "WARNING: Raw data is empty." << endl;
        continue;
      }
      htitl += "; Tick; ADC count";
      ph = new TH1F(hname.c_str(), htitl.c_str(), nsam, 0, nsam);
      hists.push_back(ph);
      float sigMin = acd.raw[0];
      float sigMax = sigMin;
      for ( Index isam=0; isam<nsam; ++isam ) {
        float sig = acd.raw[isam];
        ph->SetBinContent(isam+1, sig);
        if ( isam >= m_plotSamMin && isam < m_plotSamMax ) {
          if ( sig < sigMin ) sigMin = sig;
          if ( sig > sigMax ) sigMax = sig;
        }
      }
      res.setFloat("plotSigMin_" + type, sigMin);
      res.setFloat("plotSigMax_" + type, sigMax);
    } else if ( type == "rawdist" ) {
      Index nsam = acd.raw.size();
      if ( nsam == 0 ) {
        cout << myname << "WARNING: Raw data is empty." << endl;
        continue;
      }
      htitl += "; ADC count; # samples";
      unsigned int nadc = 4096;
      ph = new TH1F(hname.c_str(), htitl.c_str(), nadc, 0, nadc);
      hists.push_back(ph);
      float sigMin = acd.raw[0];
      float sigMax = sigMin;
      for ( Index isam=0; isam<nsam; ++isam ) {
        float sig = acd.raw[isam];
        ph->Fill(sig);
        if ( isam >= m_plotSamMin && isam < m_plotSamMax ) {
          if ( sig < sigMin ) sigMin = sig;
          if ( sig > sigMax ) sigMax = sig;
        }
      }
      res.setFloat("plotSigMin_" + type, sigMin);
      res.setFloat("plotSigMax_" + type, sigMax);
    } else if ( type == "prepared" ) {
      Index nsam = acd.samples.size();
      if ( nsam == 0 ) {
        cout << myname << "WARNING: Prepared data is empty." << endl;
        continue;
      }
      htitl += "; Tick; Signal";
      if ( acd.sampleUnit.size() ) {
        htitl += " [" + acd.sampleUnit + "]";
      }
      ph = new TH1F(hname.c_str(), htitl.c_str(), nsam, 0, nsam);
      hists.push_back(ph);
      float sigMin = acd.raw[0];
      float sigMax = sigMin;
      for ( Index isam=0; isam<nsam; ++isam ) {
        float sig = acd.samples[isam];
        ph->SetBinContent(isam+1, sig);
        if ( isam >= m_plotSamMin && isam < m_plotSamMax ) {
          if ( sig < sigMin ) sigMin = sig;
          if ( sig > sigMax ) sigMax = sig;
        }
      }
      res.setFloat("plotSigMin_" + type, sigMin);
      res.setFloat("plotSigMax_" + type, sigMax);
    } else {
      cout << myname << "WARNING: Unknown type: " << type << endl;
    }
    if ( ph == nullptr ) continue;
    ph->SetStats(0);
    ph->SetLineWidth(2);
    bool resManage = m_phm == nullptr;
    if ( ! resManage ) {
      int rstat = m_phm->manage(ph);
      if ( rstat ) {
        cout << myname << "WARNING: Attempt to manage histogram " << ph->GetName()
             << " with manager " << m_HistManager << " returned error " << rstat << endl;
        resManage = true;
      }
    }
    res.setHist(type, ph, resManage);
  }
  if ( pfile != nullptr ) {
    if ( m_LogLevel >= 2 ) cout << myname << "Writing to  " << pfile->GetName() << endl;
    for ( TH1* ph : hists ) ph->Write();
    if ( m_LogLevel >= 3 ) {
      cout << myname << "File listing: " << endl;
      pfile->ls();
    }
    pfile->Close();
    delete pfile;
    gDirectory = polddir;
  }
  return res;
}

//**********************************************************************

DataMap AdcChannelPlotter::viewMap(const AdcChannelDataMap& acds) const {
  const string myname = "AdcChannelPlotter::view: ";
  DataMap resall;
  bool doPlots = m_PlotFileName.size();
  Index nx = 1;
  Index ny = 0;
  Index iplt = 0;
  using ManMap = std::map<string, TPadManipulator>;
  using NameMap = std::map<string, string>;
  ManMap mans;
  if ( doPlots ) ny = 8;  // For now 8 histos/plot
  Index nplt = nx*ny;
  NameMap pfnames;
  std::vector<TLatex*> labs;
  for ( const AdcChannelDataMap::value_type& iacd : acds ) {
    Index icha = iacd.first;
    string schan = std::to_string(icha);
    TLatex* ptxt = new TLatex(0.98, 0.14, schan.c_str());
    ptxt->SetNDC();
    ptxt->SetTextFont(42);
    ptxt->SetTextSize(0.16);
    ptxt->SetTextAlign(31);
    labs.push_back(ptxt);
    const AdcChannelData& acd = iacd.second;
    DataMap res = view(acd);
    if ( doPlots ) {
      if ( mans.size() == 0 ) {
        for ( string type : m_HistTypes ) {
          TPadManipulator& man = mans[type];
          man.setCanvasSize(1400, 1000);
          man.split(nx,ny);
          for ( Index ipad=0; ipad<nplt; ++ipad ) {
            man.man(ipad)->addHorizontalModLines(64);
            man.man(ipad)->setRangeX(m_plotSamMin, m_plotSamMax);
          }
          pfnames[type] = nameReplace(m_PlotFileName, acd, type);
        }
      }
      for ( string type : m_HistTypes ) {
        TH1* ph = res.getHist(type);
        TPadManipulator& man = *mans[type].man(iplt);
        man.add(ph, "hist", false);
        if ( type == "raw" ) {
          Index dSigMin = 100;
          Index gSigMin = res.getFloat("plotSigMin_" + type);
          Index gSigMax = res.getFloat("plotSigMax_" + type) + 0.999;
          if ( gSigMax - gSigMin < dSigMin ) {
            while ( gSigMax - gSigMin < dSigMin ) {
              if ( gSigMin > 0 ) --gSigMin;
              if ( gSigMax - gSigMin < dSigMin ) ++gSigMax;
            }
            man.setRangeY(gSigMin, gSigMax);
          }
        }
        man.add(ptxt);
      }
      if ( ++iplt >= nplt ) {
        for ( string type : m_HistTypes ) mans[type].print(pfnames[type]);
        mans.clear();
        pfnames.clear();
        iplt = 0;
      }
    }
  }
  if ( mans.size() ) {
    for ( string type : m_HistTypes ) mans[type].print(pfnames[type]);
  }
  for ( TLatex* ptxt : labs ) delete ptxt;
  return resall;
}

//**********************************************************************

string AdcChannelPlotter::
nameReplace(string name, const AdcChannelData& acd, string type) const {
  const AdcChannelStringTool* pnbl = m_adcStringBuilder;
  string nameout = name;
  StringManipulator sman(nameout);
  if ( type.size() ) sman.replace("%TYPE%", type);
  if ( pnbl == nullptr ) return nameout;
  DataMap dm;
  return pnbl->build(acd, dm, nameout);
}

//**********************************************************************
