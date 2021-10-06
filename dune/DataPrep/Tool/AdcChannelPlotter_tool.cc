// AdcChannelPlotter_tool.cc

#include "AdcChannelPlotter.h"
#include "dune/DuneCommon/Utility/StringManipulator.h"
#include "dune/DuneInterface/Tool/AdcChannelStringTool.h"
#include "dune/DuneCommon/Utility/TPadManipulator.h"
#include "dune/ArtSupport/DuneToolManager.h"
#include "larevt/CalibrationDBI/Interface/ChannelStatusService.h"
#include "larevt/CalibrationDBI/Interface/ChannelStatusProvider.h"
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
  m_PlotSamMin(ps.get<Index>("PlotSamMin")),
  m_PlotSamMax(ps.get<Index>("PlotSamMax")),
  m_PlotSigOpt(ps.get<string>("PlotSigOpt")),
  m_PlotSigMin(ps.get<float>("PlotSigMin")),
  m_PlotSigMax(ps.get<float>("PlotSigMax")),
  m_PlotDistMin(ps.get<float>("PlotDistMin")),
  m_PlotDistMax(ps.get<float>("PlotDistMax")),
  m_ColorBad(ps.get<Index>("ColorBad")),
  m_ColorNoisy(ps.get<Index>("ColorNoisy")),
  m_LabelSize(ps.get<float>("LabelSize")),
  m_SkipFlags(ps.get<IndexVector>("SkipFlags")) {
  const string myname = "AdcChannelPlotter::ctor: ";
  if ( m_HistTypes.size() == 0 ) {
    cout << myname << "WARNING: No histogram types are specified." << endl;
    return;
  }
  DuneToolManager* ptm = DuneToolManager::instance();
  string stringBuilder = "adcStringBuilder";
  m_adcStringBuilder = ptm->getShared<AdcChannelStringTool>(stringBuilder);
  if ( m_adcStringBuilder == nullptr ) {
    cout << myname << "WARNING: AdcChannelStringTool not found: " << stringBuilder << endl;
  }
  if ( m_ColorBad || m_ColorNoisy ) {
    if ( m_LogLevel >= 1 ) cout << myname << "Fetching channel status service." << endl;
    m_pChannelStatusProvider = &art::ServiceHandle<lariov::ChannelStatusService>()->GetProvider();
    if ( m_pChannelStatusProvider == nullptr ) {
      cout << myname << "WARNING: Channel status provider not found." << endl;
      m_ColorBad = 0;
      m_ColorNoisy = 0;
    }
  }
  for ( Index flg : m_SkipFlags ) m_skipFlags.insert(flg);
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
    cout << myname << "    PlotSamMin: " << m_PlotSamMin << endl;
    cout << myname << "    PlotSamMax: " << m_PlotSamMax << endl;
    cout << myname << "    PlotSigOpt: " << m_PlotSigOpt << endl;
    cout << myname << "    PlotSigMin: " << m_PlotSigMin << endl;
    cout << myname << "    PlotSigMax: " << m_PlotSigMax << endl;
    cout << myname << "   PlotDistMin: " << m_PlotDistMin << endl;
    cout << myname << "   PlotDistMax: " << m_PlotDistMax << endl;
    cout << myname << "      ColorBad: " << m_ColorBad << endl;
    cout << myname << "    ColorNoisy: " << m_ColorNoisy << endl;
    cout << myname << "     LabelSize: " << m_LabelSize << endl;
    cout << myname << "      SkipFlags: [";
    first = true;
    for ( Index flg : m_SkipFlags ) {
       if ( first ) first = false;
       else cout << ", ";
       cout << flg;
    }
    cout << "]" << endl;
  }
}

//**********************************************************************

AdcChannelPlotter::~AdcChannelPlotter() {
}

//**********************************************************************

DataMap AdcChannelPlotter::view(const AdcChannelData& acd) const {
  const string myname = "AdcChannelPlotter::view: ";
  if ( m_LogLevel >= 3 ) cout << myname << "Processing channel " << acd.channel() << endl;
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
  bool resManage = true;   // Does returned data map manage the histogram
  for ( string type : m_HistTypes ) {
    TH1* ph = nullptr;
    string hname = nameReplace(hnameBase, acd, type);
    string htitl = nameReplace(htitlBase, acd, type);
    bool isRawDist = type == "rawdist" || type == "rawdistlog";
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
        if ( isam >= m_PlotSamMin && isam < m_PlotSamMax ) {
          if ( sig < sigMin ) sigMin = sig;
          if ( sig > sigMax ) sigMax = sig;
        }
      }
      res.setFloat("plotSigMin_" + type, sigMin);
      res.setFloat("plotSigMax_" + type, sigMax);
    } else if ( isRawDist ) {
      Index nsam = acd.raw.size();
      if ( nsam == 0 ) {
        cout << myname << "WARNING: Raw data is empty." << endl;
        continue;
      }
      // Flag samples to keep in pedestal fit.
      vector<bool> keep(nsam, true);
      Index nkeep = 0;
      for ( Index isam=0; isam<nsam; ++isam ) {
        if ( isam >= acd.flags.size() ) {
          if ( m_LogLevel >= 2 ) cout << myname << "WARNING: flags are missing." << endl;
          break;
        }
        Index flg = acd.flags[isam];
        if ( m_skipFlags.count(flg) ) keep[isam] = false;
        else ++nkeep;
      }
      if ( nkeep == 0 ) {
        if ( m_LogLevel >= 2 ) cout << myname << "WARNING: No raw data is selected." << endl;
        return res.setStatus(2);
      }
      // Create a new histogram if %EVENT% appears in the name.
      bool useExistingHist = m_HistName.find("%EVENT%") == string::npos;
      if ( useExistingHist ) {
        HistMap::iterator ihst = getState().hists.find(hname);
        if ( ihst != getState().hists.end() ) ph = ihst->second;
        resManage = false;
      }
      if ( ph == nullptr ) {
        htitl += "; ADC count; # samples";
        unsigned int nadc = 4096;
        ph = new TH1F(hname.c_str(), htitl.c_str(), nadc, 0, nadc);
        if ( ! useExistingHist ) hists.push_back(ph);
        if ( useExistingHist ) getState().hists[hname] = ph;
      }
      float sigMin = acd.raw[0];
      float sigMax = sigMin;
      for ( Index isam=0; isam<nsam; ++isam ) {
        float sig = acd.raw[isam];
        if ( keep[isam] ) ph->Fill(sig);
        if ( isam >= m_PlotSamMin && isam < m_PlotSamMax ) {
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
      float sigMin = acd.samples[0];
      float sigMax = sigMin;
      for ( Index isam=0; isam<nsam; ++isam ) {
        float sig = acd.samples[isam];
        ph->SetBinContent(isam+1, sig);
        if ( isam >= m_PlotSamMin && isam < m_PlotSamMax ) {
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
    if ( m_ColorBad && m_pChannelStatusProvider->IsBad(acd.channel()) ) {
      ph->SetLineColor(m_ColorBad);
    }
    if ( m_ColorNoisy && m_pChannelStatusProvider->IsNoisy(acd.channel()) ) {
      ph->SetLineColor(m_ColorNoisy);
    }
    res.setHist(type, ph, resManage);
  }
  if ( pfile != nullptr ) {
    if ( m_LogLevel >= 2 ) cout << myname << "Writing to  " << pfile->GetName() << endl;
    for ( TH1* ph : hists ) ph->Write();
    if ( m_LogLevel >= 4 ) {
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
  const string myname = "AdcChannelPlotter::viewMap: ";
  DataMap resall;
  bool doPlots = m_PlotFileName.size();
  Index nx = 1;
  Index ny = 8;
  Index nplt = nx*ny;
  Index ndx = 4;
  Index ndy = 4;
  Index ndplt = ndx*ndy;
  using ManMap = std::map<string, TPadManipulator>;
  using NameMap = std::map<string, string>;
  using IndexMap = std::map<string, Index>;
  ManMap mans;
  IndexMap iplts;
  IndexMap nplts;
  NameMap pfnames;
  std::vector<TLatex*> labs;
  bool useViewPort = true;
  for ( const AdcChannelDataMap::value_type& iacd : acds ) {
    Index icha = iacd.first;
    string schan = std::to_string(icha);
    TLatex* ptxt = new TLatex(0.98, 0.025, schan.c_str());
    ptxt->SetNDC();
    ptxt->SetTextFont(42);
    ptxt->SetTextAlign(31);
    labs.push_back(ptxt);
    const AdcChannelData& acd = iacd.second;
    DataMap res = view(acd);
    if ( doPlots ) {
      for ( string type : m_HistTypes ) {
        bool isRaw = type == "raw";
        bool isRawDist = type == "rawdist" || type == "rawdistlog";
        bool isLogY = type == "rawdistlog";
        float marginTop    = 0.0;
        float marginBottom = isRawDist ? 0.12 : 0.09;
        float marginLeft   = isRawDist ? 0.12 : 0.05;
        float marginRight  = isRawDist ? 0.02 : 0.01;
        float xlab = isRawDist ? 0.95 : 0.98;
        float ylab = 0.05 + marginBottom;
        float hlab = isRawDist ? 0.08 : 0.16;
        ptxt->SetX(xlab);
        ptxt->SetY(ylab);
        ptxt->SetTextSize(hlab);
        if ( mans.find(type) == mans.end() ) {
          if ( m_LogLevel >= 3 ) cout << "Creating new top-level plot of type " << type << "." << endl;
          TPadManipulator& topman = mans[type];
          topman.setCanvasSize(1400, 1000);
          if ( useViewPort ) {
            float xview1 = 0.0;
            float yview1 = isRawDist ? 0.0 : 0.00;
            float xview2 = 1.0;
            float yview2 = isRawDist ? 0.96 : 0.96;
            if ( topman.addPad(xview1, yview1, xview2, yview2) ) {
              cout << myname << "ERROR: Unable to add subpad." << endl;
              abort();
            }
          }
          TPadManipulator& man = *topman.man(0);
          if ( isRaw || type == "prepared" ) {
            man.split(nx,ny);
            for ( Index ipad=0; ipad<nplt; ++ipad ) {
              if ( isRaw && (m_PlotSigMax - m_PlotSigMin) < 1001 ) man.man(ipad)->addHorizontalModLines(64);
              man.man(ipad)->setRangeX(m_PlotSamMin, m_PlotSamMax);
              man.addAxis();
            }
            nplts[type] = nx*ny;
          } else if ( isRawDist ) {
            man.split(ndx, ndy);
            for ( Index ipad=0; ipad<ndplt; ++ipad ) {
              man.man(ipad)->addVerticalModLines(64);
              man.man(ipad)->setRangeX(m_PlotSigMin, m_PlotSigMax);
              man.man(ipad)->showUnderflow();
              man.man(ipad)->showOverflow();
              if ( type == "rawdistlog" ) man.man(ipad)->setLogY();
              man.addAxis();
            }
            nplts[type] = ndx*ndy;
          }
          pfnames[type] = nameReplace(m_PlotFileName, acd, type);
          iplts[type] = 0;
          if ( m_LabelSize ) {
            man.setLabelSizeX(m_LabelSize);
            man.setLabelSizeY(m_LabelSize);
            //man.setTitleSize(m_LabelSize);
          }
          for ( Index ipsm=0; ipsm<man.npad(); ++ipsm ) {
            TPadManipulator* psman = man.man(ipsm);
             psman->setMarginTop(marginTop);
            psman->setMarginBottom(marginBottom);
            psman->setMarginLeft(marginLeft);
            psman->setMarginRight(marginRight);
          }
        }
        Index& iplt = iplts[type];
        Index nplt = nplts[type];
        TH1* ph = res.getHist(type);
        if ( m_LogLevel >= 3 ) cout << myname << "  Adding subplot " << iplt << " for type " << type << "." << endl;
        TPadManipulator& topman = mans[type];
        string sttl = ph->GetTitle();
        Index ipos = sttl.find(" channel");
        sttl = sttl.substr(0, ipos);
        topman.setTitle(sttl, 0.025);
        TPadManipulator& man = useViewPort ? *topman.man(0)->man(iplt) : *topman.man(iplt);
        man.add(ph, "hist", false);
        man.setTitle("");
        if ( type == "raw" || type == "prepared" ) {
          float ymin = m_PlotSigMin;
          float ymax = m_PlotSigMax;
          if ( m_PlotSigOpt == "pedestal" ) {
            if ( type == "raw" ) {
              ymin += acd.pedestal;
              ymax += acd.pedestal;
            } else {
              cout << myname << "Invalid raw PlotSigOpt = " << m_PlotSigOpt << ". Using fixed." << endl;
            }
          } else if ( m_PlotSigOpt == "full" ) {
            int dSigMin = m_PlotSigMin;
            int gSigMin = res.getFloat("plotSigMin_" + type);
            int gSigMax = res.getFloat("plotSigMax_" + type) + 0.999;
            if ( gSigMax - gSigMin < dSigMin ) {
              while ( gSigMax - gSigMin < dSigMin ) {
                if ( gSigMin > 0 ) --gSigMin;
                if ( gSigMax - gSigMin < dSigMin ) ++gSigMax;
              }
            }
            ymin = gSigMin;
            ymax = gSigMax;
          } else if ( m_PlotSigOpt != "fixed" ) {
            cout << myname << "Invalid " << type << " PlotSigOpt = " << m_PlotSigOpt << ". Using fixed." << endl;
          }
          man.setRangeY(ymin, ymax);
          man.add(ptxt);
          // For raw data, add line showing the pedestal.
          if ( type == "raw" && acd.pedestal > ymin && acd.pedestal < ymax ) { 
            man.addHorizontalLine(acd.pedestal);
          }
          // For prepared data, add line showing zero.
          if ( type == "prepared" ) { 
            man.addHorizontalLine(0.0);
          }
          
        } else if ( isRawDist ) {
          if ( m_PlotSigOpt == "fixed" ) {
            float xmin = m_PlotSigMin;
            float xmax = m_PlotSigMax;
            man.setRangeX(xmin, xmax);
          } else if ( m_PlotSigOpt == "pedestal" ) {
            float xmin = acd.pedestal + m_PlotSigMin;
            float xmax = acd.pedestal + m_PlotSigMax;
            man.setRangeX(xmin, xmax);
          } else if ( m_PlotSigOpt != "full" ) {
            cout << myname << "Invalid rawdist PlotSigOpt = " << m_PlotSigOpt << ". Using full." << endl;
          }
          if ( m_PlotDistMax > m_PlotDistMin ) {
            if ( isLogY ) man.setLogRangeY(m_PlotDistMin, m_PlotDistMax);
            else          man.setRangeY(m_PlotDistMin, m_PlotDistMax);
          }
          man.add(ptxt);
        }
        man.addAxis();
        if ( ++iplt >= nplt ) {
          if ( m_LogLevel >= 3 ) cout << "Printing plot " << pfnames[type] << endl;
          for ( string type : m_HistTypes ) mans[type].print(pfnames[type]);
          mans.erase(type);
          pfnames.erase(type);
        }
      }
    }
  }
  // Print any left over (i.e. partial) plots.
  for ( ManMap::value_type& iman : mans ) iman.second.print(pfnames[iman.first]);
  for ( TLatex* ptxt : labs ) delete ptxt;
  return resall;
}

//**********************************************************************

string AdcChannelPlotter::
nameReplace(string name, const AdcChannelData& acd, string type) const {
  const AdcChannelStringTool* pnbl = m_adcStringBuilder;
  string nameout = name;
  StringManipulator sman(nameout, false);
  if ( type.size() ) sman.replace("%TYPE%", type);
  if ( pnbl == nullptr ) return nameout;
  DataMap dm;
  return pnbl->build(acd, dm, nameout);
}

//**********************************************************************

DEFINE_ART_CLASS_TOOL(AdcChannelPlotter)
