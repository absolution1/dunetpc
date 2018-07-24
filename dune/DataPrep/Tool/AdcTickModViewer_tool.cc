// AdcTickModViewer_tool.cc

#include "AdcTickModViewer.h"
#include "dune/DataPrep/Utility/StickyCodeMetrics.h"
#include "dune/DuneCommon/TPadManipulator.h"
#include "dune/DuneCommon/StringManipulator.h"
#include "dune/DuneInterface/Tool/AdcChannelStringTool.h"
#include "dune/DuneInterface/Tool/TimeOffsetTool.h"
#include "dune/ArtSupport/DuneToolManager.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include "TDirectory.h"
#include "TFile.h"
#include "TH1F.h"
#include "TF1.h"
#include "TTree.h"
#include "TROOT.h"
#include "TError.h"

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

using Index = AdcTickModViewer::Index;
using Name = AdcTickModViewer::Name;

//**********************************************************************
// Class methods.
//**********************************************************************

AdcTickModViewer::AdcTickModViewer(fhicl::ParameterSet const& ps)
: m_LogLevel(ps.get<int>("LogLevel")),
  m_TickModPeriod(ps.get<Index>("TickModPeriod")),
  m_TimeOffsetTool(ps.get<Name>("TimeOffsetTool")),
  m_FitRmsMin(ps.get<float>("FitRmsMin")),
  m_FitRmsMax(ps.get<float>("FitRmsMax")),
  m_HistName(ps.get<string>("HistName")),
  m_HistTitle(ps.get<string>("HistTitle")),
  m_HistChannelCount(ps.get<Index>("HistChannelCount")),
  m_PlotFileName(ps.get<string>("PlotFileName")),
  m_RootFileName(ps.get<string>("RootFileName")),
  m_TreeFileName(ps.get<string>("TreeFileName")),
  m_PlotChannels(ps.get<IndexVector>("PlotChannels")),
  m_PlotSizeX(ps.get<Index>("PlotSizeX")),
  m_PlotSizeY(ps.get<Index>("PlotSizeY")),
  m_PlotShowFit(ps.get<Index>("PlotShowFit")),
  m_PlotSplitX(ps.get<Index>("PlotSplitX")),
  m_PlotSplitY(ps.get<Index>("PlotSplitY")),
  m_PlotWhich(ps.get<Index>("PlotWhich")),
  m_PlotFrequency(ps.get<Index>("PlotFrequency")),
  m_tickOffsetTool(nullptr),
  m_state(new State) 
{
  const string myname = "AdcTickModViewer::ctor: ";
  DuneToolManager* ptm = DuneToolManager::instance();
  string stringBuilder = "adcStringBuilder";
  m_adcStringBuilder = ptm->getShared<AdcChannelStringTool>(stringBuilder);
  if ( m_adcStringBuilder == nullptr ) {
    cout << myname << "WARNING: AdcChannelStringTool not found: " << stringBuilder << endl;
  }
  string tname = m_TimeOffsetTool;
  if ( tname.size() ) {
    m_tickOffsetTool = ptm->getShared<TimeOffsetTool>(tname);
    if ( m_tickOffsetTool == nullptr ) {
      cout << myname << "WARNING: Requested TimeOffsetTool not found: " << tname << endl;
    }
  }
  if ( m_LogLevel >= 1 ) {
    cout << myname << "Configuration parameters:" << endl;
    cout << myname << "          LogLevel: " << m_LogLevel << endl;
    cout << myname << "     TickModPeriod: " << m_TickModPeriod << endl;
    cout << myname << "    TimeOffsetTool: " << m_TimeOffsetTool << endl;
    cout << myname << "          HistName: " << m_HistName << endl;
    cout << myname << "         HistTitle: " << m_HistTitle << endl;
    cout << myname << "  HistChannelCount: " << m_HistChannelCount << endl;
    cout << myname << "      PlotFileName: " << m_PlotFileName << endl;
    cout << myname << "      RootFileName: " << m_RootFileName << endl;
    cout << myname << "      TreeFileName: " << m_TreeFileName << endl;
    cout << myname << "      PlotChannels: [";
    bool first = true;
    for ( Index icha : m_PlotChannels ) {
      if ( !first ) cout << ", ";
      else first = false;
      cout << icha;
    }
    cout << "]" << endl;
    cout << myname << "         PlotSizeX: " << m_PlotSizeX << endl;
    cout << myname << "         PlotSizeY: " << m_PlotSizeY << endl;
    cout << myname << "       PlotShowFit: " << m_PlotShowFit << endl;
    cout << myname << "        PlotSplitX: " << m_PlotSplitX << endl;
    cout << myname << "         PlotWhich: " << m_PlotWhich << endl;
    cout << myname << "     PlotFrequency: " << m_PlotFrequency << endl;
  }
  if ( m_LogLevel >=4 ) {
    cout << myname << "INFO: Checking state." << endl;
    cout << myname << "INFO: Full hist map size: " << state().ChannelTickModFullHists.size() << endl;
    cout << myname << "INFO: Proc hist map size: " << state().ChannelTickModProcHists.size() << endl;
  }
}

//**********************************************************************

AdcTickModViewer::~AdcTickModViewer() {
  const string myname = "AdcTickModViewer::dtor: ";
  if ( m_LogLevel >= 1 ) cout << myname << "Closing." << endl;
  Index nplot = 0;
  processAccumulation(nplot);
  if ( m_LogLevel >= 1 ) cout << myname << "Plot count: " << nplot << endl;
}
  
//**********************************************************************

DataMap AdcTickModViewer::view(const AdcChannelData& acd) const {
  const string myname = "AdcTickModViewer::view: ";
  DataMap res;
  Index icha = acd.channel;
  if ( m_LogLevel >= 3 ) cout << myname << "Processing channel " << icha << endl;
  HistVector& tmhsFull = state().ChannelTickModFullHists[icha];
  Index ntkm = m_TickModPeriod;
  if ( tmhsFull.size() == 0 ) tmhsFull.resize(ntkm, nullptr);
  Index itkm0 = 0;
  if ( m_tickOffsetTool != nullptr ) {
    TimeOffsetTool::Data dat;
    dat.run = acd.run;
    dat.subrun = acd.subRun;
    dat.event = acd.event;
    dat.channel = acd.channel;
    TimeOffsetTool::Offset off = m_tickOffsetTool->offset(dat);
    if ( ! off.isValid() ) {
      cout << myname << "Error finding tick offset: " << off.status << endl;
      return res.setStatus(1);
    }
    if ( off.unit != "tick" ) {
      cout << myname << "Time offset has wrong unit: " << off.unit << endl;
      return res.setStatus(2);
    }
    long toff = off.value;
    while ( toff < 0.0 ) toff += m_TickModPeriod;
    itkm0 = toff % m_TickModPeriod;
    if ( m_LogLevel >= 3 ) cout << myname << "Using tick offset " << itkm0 << endl;
  }
  if ( state().run < 0 && acd.run != acd.badIndex ) state().run = acd.run;
  for ( Index itkm=0; itkm<ntkm; ++itkm ) {
    fillChannelTickMod(acd, itkm0, itkm);
  }
  Index nplot = 0;
  if ( m_PlotFrequency ) {
    processAccumulation(nplot);
    HistVectorMap& ctmprocs = state().ChannelTickModProcHists;
    if ( ctmprocs.find(icha) == ctmprocs.end() ) {
      static HistVector empty;
      res.setHistVector("tmHists", empty);
    } else {
      res.setHistVector("tmHists", ctmprocs[icha]);
    }
  }
  res.setHistVector("tmWideHists", tmhsFull);   // Passing out hist sthat will be updated!
  res.setInt("tmCount", ntkm);
  res.setInt("tmPlotCount", nplot);
  return res;
}

//**********************************************************************

Name AdcTickModViewer::nameReplace(Name nameIn, const AdcChannelData& acd, Index itkm) const {
  const string myname = "AdcTickModViewer::nameReplace: ";
  DataMap dm;
  string nameOut = m_adcStringBuilder->build(acd, dm, nameIn);
  string stkm = std::to_string(itkm);
  string::size_type stkm0Len = std::to_string(m_TickModPeriod-1).size();
  string stkm0 = stkm;
  while ( stkm0.size() < stkm0Len ) stkm0 = "0" + stkm0;
  StringManipulator sman(nameOut);
  sman.replace("%TICKMOD%", stkm);
  sman.replace("%0TICKMOD%", stkm0);
  return nameOut;
}

//**********************************************************************

int
AdcTickModViewer::fillChannelTickMod(const AdcChannelData& acd, Index itkm0, Index itkm) const {
  const string myname = "AdcTickModViewer::fillChannelTickMod: ";
  DataMap res;
  Index nsam = acd.raw.size();
  if ( nsam == 0 ) {
    if ( m_LogLevel >= 2 ) cout << myname << "WARNING: Raw data is empty." << endl;
    return 1;
  }
  HistPtr& ph = state().ChannelTickModFullHists[acd.channel][itkm];
  if ( ph == nullptr ) {
    string hname = nameReplace(m_HistName, acd, itkm);
    string htitl = nameReplace(m_HistTitle, acd, itkm);
    htitl += "; ADC count; # samples";
    unsigned int nadc = 4096;
    if ( m_LogLevel >= 2 ) cout << myname << "Creating histogram " << hname
                                << " for channel " << acd.channel
                                << " tickmod " << itkm << endl;
    ph.reset(new TH1F(hname.c_str(), htitl.c_str(), nadc, 0, nadc));
    ph->SetDirectory(0);
  }
  Index period = m_TickModPeriod;
  if ( m_LogLevel >= 4 ) cout << myname << "Filling hist " << ph->GetName()
                              << " with channel " << acd.channel
                              << " tickmod " << itkm << endl;
  Index isam0 = (itkm + period - itkm0) % period;
  for ( Index isam=isam0; isam<nsam; isam+=period ) ph->Fill(acd.raw[isam]);
  return 0;
}

//**********************************************************************

int AdcTickModViewer::processAccumulatedChannel(Index icha, Index& nplot) const {
  const string myname = "AdcTickModViewer::processAccumulatedChannel: ";
  nplot = 0;
  if ( state().ChannelTickModFullHists.find(icha) == state().ChannelTickModFullHists.end() ) return 1;
  const HistVector& tmhsFull = state().ChannelTickModFullHists[icha];
  HistVector& tmhsProc = state().ChannelTickModProcHists[icha];
  Index ntkm = tmhsFull.size();
  if ( m_LogLevel >= 3 ) {
    cout << myname << "Tickmod hist count: " << ntkm << endl;
  }
  if ( tmhsProc.size() == 0 ) tmhsProc.resize(ntkm, nullptr);
  // Fetch the tree.
  TTree* ptree = state().tickmodTree;
  TickModTreeData& data = state().treedata;
  if ( ptree != nullptr ) {
    data.run = state().run;
    data.chan = icha;
  }
  // Loop over tickmods.
  for ( Index itkm=0; itkm<ntkm; ++itkm ) {
    const HistPtr& ph = tmhsFull[itkm];
    // Process histograms with the sticky code utility.
    Index chmod = 10;
    StickyCodeMetrics scm(ph->GetName(), ph->GetTitle(), m_HistChannelCount, chmod);
    if ( scm.evaluate(ph.get()) ) {
      cout << myname << "Sticky code evaluation failed for channel " << icha
           << " tickmod " << itkm << endl;
      //tmhsProc[itkm].reset(ph);
    } else {
      tmhsProc[itkm] = scm.getSharedHist();
      if ( ptree != nullptr ) {
        data.itkm = itkm;
        data.fill(scm);
        ptree->Fill();
      }
    }
  }
  // Draw the histograms.
  makeTickModPlots(icha, nplot);
  if ( m_LogLevel >= 3 ) {
    cout << myname << "  Plot file count for channel " << icha << ": " << nplot << endl;
  }
  return 0;
}

//**********************************************************************

int AdcTickModViewer::processAccumulation(Index& nplot) const {
  const string myname = "AdcTickModViewer::processAccumulation: ";
  Index nhst = state().ChannelTickModFullHists.size();
  if ( m_LogLevel >= 2 ) {
    cout << myname << "Channel count: " << nhst << endl;
  }
  if ( nhst == 0 ) return 0;
  // Create tree to hold results.
  TTree*& ptree = state().tickmodTree;
  TFile*& pfile = state().pfile;
  if ( m_TreeFileName.size() ) {
    TDirectory* psavdir = gDirectory;
    pfile = TFile::Open(m_TreeFileName.c_str(), "CREATE");
    if ( pfile->IsOpen() ) {
      ptree = new TTree("tickmod", "TickMod tree");
      //ptree->Branch("data", &mydat, "Tickmod data");
      //ptree->Branch("run", &(state().run));
      //ptree->Branch("which", &m_PlotWhich, "i");
/*
      static int myint = 123;
      ptree->Branch("myint", &myint, "myint/I");
      static TickModTreeData* pdumdat = new TickModTreeData;
      ptree->Branch("dumdat", pdumdat, "TickModTreeData");
      ptree->Branch("pdumdat", &pdumdat, "TickModTreeData");
      ptree->Branch("pdumdat2", "TickModTreeData", &pdumdat);
      TickModTreeData& mydat = state().treedata;
      ptree->Branch("mydata", &mydat);
*/
      ptree->Branch("data", &(state().treedata), 64000, 1);
    } else {
      cout << myname << "Unable to open file " << m_TreeFileName << endl;
    }
    psavdir->cd();
  }
  Index nplotTot = 0;
  int rstat = 0;
  for ( HistVectorMap::value_type icvm : state().ChannelTickModFullHists ) {
    Index icha = icvm.first;
    Index nplot = 0;
    rstat += processAccumulatedChannel(icha, nplot);
    nplotTot += nplot;
  }
  if ( m_LogLevel >= 2 ) {
    cout << myname << "  Total plot file count: " << nplotTot << endl;
  }
  if ( pfile != nullptr ) {
    if ( m_LogLevel >= 2 ) {
      cout << myname << "Tree size: " << ptree->GetEntries() << endl;
    }
    pfile->Write();
    pfile->Close();
    delete pfile;
    pfile = nullptr;
    ptree = nullptr;
  }
  return rstat;
}

//**********************************************************************

int AdcTickModViewer::makeTickModPlots(Index icha, Index& nplot) const {
  const string myname = "AdcTickModViewer::makeTickModPlots: ";
  nplot = 0;
  if ( m_PlotFileName.size() == 0 ) return 0;
  if ( m_PlotChannels.size() ) {
    if ( find(m_PlotChannels.begin(), m_PlotChannels.end(), icha) == m_PlotChannels.end() ) {
      if ( m_LogLevel >= 3 )  cout << myname << "Skipping channel not in PlotChannels: " << icha << endl;
      return 0;
    }
  }
  Index npad = 0;
  Index npadx = 0;
  Index npady = 0;
  if ( m_PlotFileName.size() && m_PlotSplitX > 0 ) {
    npadx = m_PlotSplitX;
    npady = m_PlotSplitY ? m_PlotSplitY : m_PlotSplitX;
    npad = npadx*npady;
  }
  Name plotFileName;
  Index ntkm = m_TickModPeriod;
  const HistVector& tmhs = state().ChannelTickModProcHists[icha];
  vector<IndexVector> showTickModVectors;  // Vectors of tickmods to plot
  vector<Name> showPlotNames;
  if ( m_PlotWhich & 1 ) {
    if ( m_LogLevel >= 3 ) cout << myname << "Add full vector of "
                                << ntkm << " tickmods." << endl;
    showTickModVectors.emplace_back(ntkm);
    for ( Index itkm=0; itkm<ntkm; ++itkm ) showTickModVectors.back()[itkm] = itkm;
    showPlotNames.push_back(m_PlotFileName);
  }
  IndexVector tkmsMin;
  IndexVector tkmsMax;
  // If needed, find min and max tickmods and build histogram vectors(s).
  if ( m_PlotWhich & 6 ) {
    Index itkmMin = 0;
    Index itkmMax = 0;
    double meanMin = tmhs[0]->GetMean();
    double meanMax = meanMin;
    for ( Index itkm=0; itkm<ntkm; ++itkm ) {
      double mean = tmhs[itkm]->GetMean();
      if ( mean > meanMax ) {
        meanMax = mean;
        itkmMax = itkm;
      }
      if ( mean < meanMin ) {
        meanMin = mean;
        itkmMin = itkm;
      }
    }
    int idel1 = npad/2;
    idel1 *= -1;
    int idel2 = idel1 + npad;
    if ( npad > ntkm ) {
      idel1 = 0;
      idel2 = ntkm;
    }
    if ( m_LogLevel >= 3 ) cout << myname << "Tick delta range: [" << idel1
                                << ", " << idel2 << ")" << endl;
    for ( int idel=idel1; idel<idel2; ++idel ) {
      Index itkm = (itkmMin + ntkm + idel) % ntkm;
      tkmsMin.push_back(itkm);
      itkm = (itkmMax + ntkm + idel) % ntkm;
      tkmsMax.push_back(itkm);
    }
    if ( m_PlotWhich & 2 ) {
      if ( m_LogLevel >= 3 ) cout << myname << "Add min vector of " << tkmsMin.size()
                                  << " tickmods." << endl;
      showTickModVectors.push_back(tkmsMin);
      Name pnam = m_PlotFileName;
      StringManipulator sman(pnam);
      sman.replace("%TICKMOD%", "Min");
      sman.replace("%0TICKMOD%", "Min");
      showPlotNames.push_back(pnam);
    }
    if ( m_PlotWhich & 4 ) {
      if ( m_LogLevel >= 3 ) cout << myname << "Add max vector of " << tkmsMax.size()
                                  << " tickmods." << endl;
      showTickModVectors.push_back(tkmsMax);
      Name pnam = m_PlotFileName;
      StringManipulator sman(pnam);
      sman.replace("%TICKMOD%", "Max");
      sman.replace("%0TICKMOD%", "Max");
      showPlotNames.push_back(pnam);
    }
  }
  AdcChannelData acd;  // For building plot file name
  acd.channel = icha;
  if ( state().run >= 0 ) acd.run = state().run;
  TPadManipulator* pmantop = nullptr;
  for ( Index ihv=0; ihv<showTickModVectors.size(); ++ihv ) {
    const IndexVector tkms = showTickModVectors[ihv];
    Name pfname = showPlotNames[ihv];
    if ( m_LogLevel >= 3 ) cout << "Plotting " << tkms.size() << " tickmods with name "
                                << pfname << endl;
    Index ipad = 0;
    Index icount = 0;
    for ( Index itkm : tkms ) {
      HistPtr ph = tmhs[itkm];
      if ( pmantop == nullptr ) {
        pmantop = new TPadManipulator;
        if ( m_PlotSizeX && m_PlotSizeY ) pmantop->setCanvasSize(m_PlotSizeX, m_PlotSizeY);
        if ( npad > 1 ) pmantop->split(npady, npady);
        plotFileName = nameReplace(pfname, acd, itkm);
      }
      TPadManipulator* pman = pmantop->man(ipad);
      pman->add(ph.get(), "hist", false);
      if ( m_PlotShowFit > 1 ) pman->addHistFun(1);
      if ( m_PlotShowFit ) pman->addHistFun(0);
      pman->addVerticalModLines(64);
      pman->showUnderflow();
      pman->showOverflow();
      ++icount;
      if ( ++ipad == npad || icount == tkms.size() ) {
        if ( m_LogLevel >= 3 ) cout << myname << "  Creating " << plotFileName << endl;
        pmantop->print(plotFileName);
        ++nplot;
        ipad = 0;
        delete pmantop;
        pmantop = nullptr;
      }
    }
  }
  return 0;
}

//**********************************************************************
