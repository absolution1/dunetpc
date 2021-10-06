// AdcPedestalFitter_tool.cc

#include "AdcPedestalFitter.h"
#include "dune/DuneCommon/Utility/TPadManipulator.h"
#include "dune/DuneInterface/Tool/AdcChannelStringTool.h"
#include "dune/DuneInterface/Tool/RunDataTool.h"
#include "dune/ArtSupport/DuneToolManager.h"
#include "dune/DuneCommon/Utility/RootParFormula.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include "TDirectory.h"
#include "TFile.h"
#include "TH1F.h"
#include "TF1.h"
#include "TROOT.h"
#include "TError.h"
#include "TFitResult.h"
#include "TFitter.h"
#include "TLatex.h"

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
using NameSet = std::set<string>;

//**********************************************************************
// Local definitions.
//**********************************************************************

namespace {

void copyMetadata(const DataMap& res, AdcChannelData& acd) {
  acd.metadata["fitPedFractionLow"] = res.getFloat("fitFractionLow");
  acd.metadata["fitPedFractionHigh"] = res.getFloat("fitFractionHigh");
  acd.metadata["fitPedestal"] = res.getFloat("fitPedestal");
  acd.metadata["fitPedRms"] = res.getFloat("fitPedestalRms");
  acd.metadata["fitPedChiSquare"] = res.getFloat("fitChiSquare");
  acd.metadata["fitPedReducedChiSquare"] = res.getFloat("fitReducedChiSquare");
  acd.metadata["fitPedPeakBinFraction"] = res.getFloat("fitPeakBinFraction");
  acd.metadata["fitPedPeakBinExcess"] = res.getFloat("fitPeakBinExcess");
  acd.metadata["fitPedNBinsRemoved"] = res.getFloat("fitNBinsRemoved");
}

//void handleRootError(int level, bool doAbort, const char* clocation, const char* cmsg) {
//  cout << "AdcPedestalFitter::handleRootError: "
//       << "Received Root error: level=" << level
//       << ", doAbort=" << (doAbort ? "true" : "false")
//       << ", location=" << clocation
//       << ": " << cmsg << endl;
//}

}  // end unnamed namespace

//**********************************************************************
// Class methods.
//**********************************************************************

AdcPedestalFitter::AdcPedestalFitter(fhicl::ParameterSet const& ps)
: m_LogLevel(ps.get<int>("LogLevel")),
  m_AdcRange(ps.get<Name>("AdcRange")),
  m_FitOpt(ps.get<Index>("FitOpt")),
  m_FitPrecision(ps.get<float>("FitPrecision")),
  m_SkipFlags(ps.get<IndexVector>("SkipFlags")),
  m_AdcFitRange(ps.get<Name>("AdcFitRange")),
  m_FitRmsMin(ps.get<Name>("FitRmsMin")),
  m_FitRmsMax(ps.get<Name>("FitRmsMax")),
  m_RemoveStickyCode(ps.get<bool>("RemoveStickyCode")),
  m_HistName(ps.get<string>("HistName")),
  m_HistTitle(ps.get<string>("HistTitle")),
  m_PlotFileName(ps.get<string>("PlotFileName")),
  m_RootFileName(ps.get<string>("RootFileName")),
  m_PlotSizeX(ps.get<Index>("PlotSizeX")),
  m_PlotSizeY(ps.get<Index>("PlotSizeY")),
  m_PlotShowFit(ps.get<Index>("PlotShowFit")),
  m_PlotSplitX(ps.get<Index>("PlotSplitX")),
  m_PlotSplitY(ps.get<Index>("PlotSplitY")),
  m_prdtool(nullptr),
  m_pstate(new State) {
  const string myname = "AdcPedestalFitter::ctor: ";
  DuneToolManager* ptm = DuneToolManager::instance();
  string stringBuilder = "adcStringBuilder";
  m_adcStringBuilder = ptm->getShared<AdcChannelStringTool>(stringBuilder);
  if ( m_adcStringBuilder == nullptr ) {
    cout << myname << "WARNING: AdcChannelStringTool not found: " << stringBuilder << endl;
  }
  for ( Index flg : m_SkipFlags ) m_skipFlags.insert(flg);
  // Set the vector of fit options.
  // These are tried in turn until one succeeds.
  {
    Name foptLsq = "WWB";
    // Nov2018: Likelihood fit better handles pdsp run 5803 event 86 channel 7309.
    Name foptLik = "LWB";
    if ( m_FitOpt == 0 ) {
      if ( m_LogLevel > 0 ) cout << myname << "Mean used in place of pedestal fit." << endl;
    } else if ( m_FitOpt == 1 ) {
      m_fitOpts.push_back(foptLsq);
      if ( m_LogLevel > 0 ) cout << myname << "Chi-square fit used for pedestal." << endl;
    } else if ( m_FitOpt == 2 ) {
      m_fitOpts.push_back(foptLik);
      if ( m_LogLevel > 0 ) cout << myname << "Likelihood fit used for pedestal." << endl;
    } else if ( m_FitOpt == 3 ) {
      if ( m_LogLevel > 0 ) cout << myname << "Chi-square/likelihood fit used for pedestal." << endl;
      m_fitOpts.push_back(foptLsq);
      m_fitOpts.push_back(foptLik);
    } else if ( m_FitOpt != 0.0 ) {
      cout << "WARNING: Invalid FitOpt: " << m_FitOpt << ". Not fit will be performed." << endl;
    }
  }
  // Fetch the formula parameters.
  m_tfs["AdcRange"]    = new RootParFormula("adcRange", m_AdcRange);
  m_tfs["AdcFitRange"] = new RootParFormula("adcFitRange", m_AdcFitRange);
  m_tfs["FitRmsMin"]   = new RootParFormula("fitRmsMin", m_FitRmsMin);
  m_tfs["FitRmsMax"]   = new RootParFormula("fitRmsMax", m_FitRmsMax);
  m_haveFormulaParams = false;
  for ( const auto& itf : m_tfs ) {
    if ( itf.second->npar() > 0 ) m_haveFormulaParams = true;
  }
  if ( m_haveFormulaParams ) {
    string stnam = "runDataTool";
    m_prdtool = ptm->getShared<RunDataTool>(stnam);
    if ( m_prdtool == nullptr ) {
      cout << myname << "ERROR: RunDataTool " << stnam
           << " not found. Formulas will not be evaluated." << endl;
    } else {
      cout << myname << "RunDataTool retrieved." << endl;
    }
  } else {
    cout << myname << "No formula parameters. RunDataTool is not retrieved." << endl;
  }
  if ( m_LogLevel >= 1 ) {
    cout << myname << "Configuration parameters:" << endl;
    cout << myname << "       LogLevel: " << m_LogLevel << endl;
    cout << myname << "       AdcRange: " << m_AdcRange << endl;
    cout << myname << "         FitOpt: " << m_FitOpt << endl;
    cout << myname << "   FitPrecision: " << m_FitPrecision << endl;
    cout << myname << "      SkipFlags: [";
    bool first = true;
    for ( Index flg : m_SkipFlags ) {
       if ( first ) first = false;
       else cout << ", ";
       cout << flg;
    }
    cout << "]" << endl;
    cout << myname << "       AdcFitRange: " << m_AdcFitRange << endl;
    cout << myname << "         FitRmsMin: " << m_FitRmsMin << endl;
    cout << myname << "         FitRmsMax: " << m_FitRmsMax << endl;
    cout << myname << "  RemoveStickyCode: " << m_RemoveStickyCode << endl;
    cout << myname << "          HistName: " << m_HistName << endl;
    cout << myname << "         HistTitle: " << m_HistTitle << endl;
    cout << myname << "      PlotFileName: " << m_PlotFileName << endl;
    cout << myname << "      RootFileName: " << m_RootFileName << endl;
    cout << myname << "         PlotSizeX: " << m_PlotSizeX << endl;
    cout << myname << "         PlotSizeY: " << m_PlotSizeY << endl;
    cout << myname << "       PlotShowFit: " << m_PlotShowFit << endl;
    cout << myname << "        PlotSplitX: " << m_PlotSplitX << endl;
    cout << myname << "        PlotSplitY: " << m_PlotSplitY << endl;
  }
}

//**********************************************************************

DataMap AdcPedestalFitter::view(const AdcChannelData& acd) const {
  const string myname = "AdcPedestalFitter::view: ";
  if ( m_LogLevel >= 3 ) cout << myname << "Calling " << endl;
  TPadManipulator* pman = nullptr;
  if ( m_PlotFileName.size() ) {
    pman = new TPadManipulator;
    if ( m_PlotSizeX && m_PlotSizeY ) pman->setCanvasSize(m_PlotSizeX, m_PlotSizeY);
  }
  DataMap res = getPedestal(acd);
  fillChannelPad(res, acd, pman);
  if ( pman != nullptr ) {
    string pfname = nameReplace(m_PlotFileName, acd, false);
    if ( m_LogLevel >= 3 ) cout << myname << "Creating plot " << pfname << endl;
    pman->print(pfname);
    delete pman;
  }
  if ( m_LogLevel >= 3 ) cout << myname << "Called " << endl;
  //TH1* phped = res.getHist("pedestal");
  //if ( res.haveHist("pedestal") ) res.setHist("pedestal", phped, true);
  return res;
}

//**********************************************************************

DataMap AdcPedestalFitter::update(AdcChannelData& acd) const {
  const string myname = "AdcPedestalFitter::update: ";
  DataMap res = view(acd);
  if ( res.status() != 0 ) return res;
  if ( m_LogLevel >= 3 ) cout << myname << "Old pedestal: " << acd.pedestal << endl;
  acd.pedestal = res.getFloat("fitPedestal");
  acd.pedestalRms = res.getFloat("fitPedestalRms");
  copyMetadata(res, acd);
  if ( m_LogLevel >= 3 ) cout << myname << "New pedestal: " << acd.pedestal << endl;
  return res;
}

//**********************************************************************

DataMap AdcPedestalFitter::updateMap(AdcChannelDataMap& acds) const {
  const string myname = "AdcPedestalFitter::updateMap: ";
  DataMap ret;
  Index nPedFitGood = 0;
  Index nPedFitFail = 0;
  if ( m_LogLevel >= 2 ) cout << myname << "Fitting " << acds.size() << " channels." << endl;
  Index npad = 0;
  Index npadx = 0;
  Index npady = 0;
  bool doPlot = m_PlotFileName.size() && m_PlotSplitX > 0;
  if ( doPlot ) {
    npadx = m_PlotSplitX;
    npady = m_PlotSplitY ? m_PlotSplitY : m_PlotSplitX;
    npad = npadx*npady;
  }
  if ( m_LogLevel >= 2 ) {
    cout << myname << "Pad count is " << npad << " (" << npady << " x " << npadx << ")" << endl;
  }
  TPadManipulator* pmantop = nullptr;
  // Loop over channels.
  Index nacd = acds.size();
  Index iacd = 0;
  vector<int> fitStats(nacd, 999);
  vector<float> fitPedestals(nacd, 0.0);
  vector<float> fitPedestalRmss(nacd, 0.0);
  string plotFileName = "";
  for ( auto& acdPair : acds ) {
    const AdcChannelData& acd = acdPair.second;
    AdcChannelData& acdMutable = acdPair.second;
    if ( m_LogLevel >= 3 ) cout << myname << "  " << iacd << ": Processing channel " << acd.channel() << endl;
    // If needed, create a new canvas.
    if ( npad > 0 && pmantop == nullptr ) {
      if ( m_LogLevel >= 3 ) cout << myname << "  Creating canvas." << endl;
      pmantop = new TPadManipulator;
      if ( m_PlotSizeX && m_PlotSizeY ) pmantop->setCanvasSize(m_PlotSizeX, m_PlotSizeY);
      if ( npad > 1 ) pmantop->split(npady, npady);
      plotFileName = nameReplace(m_PlotFileName, acd, false);
    }
    Index ipad = npad == 0 ? 0 : iacd % npad;
    TPadManipulator* pman = pmantop == nullptr ? nullptr : pmantop->man(ipad);
    DataMap tmpres = getPedestal(acd);
    fillChannelPad(tmpres, acd, pman);
    fitStats[iacd] = tmpres.status();
    float fitPedestal = 0.0;
    float fitPedestalRms = 0.0;
    if ( tmpres.status() == 0 ) {
      ++nPedFitGood;
      fitPedestal = tmpres.getFloat("fitPedestal");
      fitPedestalRms = tmpres.getFloat("fitPedestalRms");
      copyMetadata(tmpres, acdMutable);
    } else {
      ++nPedFitFail;
    }
    fitPedestals[iacd] = fitPedestal;
    fitPedestalRmss[iacd] = fitPedestalRms;
    // If needed, print and delete the canvas.
    ++iacd;
    bool lastpad = (++ipad == npad) || (iacd == nacd);
    if ( lastpad && pmantop != nullptr ) {
      if ( m_LogLevel >= 3 ) cout << myname << "  Creating plot " << plotFileName << endl;
      pmantop->print(plotFileName);
      delete pmantop;
      pmantop = nullptr;
    }
  }
  iacd = 0;
  for ( auto& acdPair : acds ) {
    AdcChannelData& acd = acdPair.second;
    acd.pedestal = fitPedestals[iacd];
    acd.pedestalRms = fitPedestalRmss[iacd];
    ++iacd;
  }
  if ( m_LogLevel >= 2 ) {
    cout << myname <<   " # good pedestal fits: " << nPedFitGood << endl;
    cout << myname <<   " # fail pedestal fits: " << nPedFitFail << endl;
  }
  ret.setInt("nPedFitGood", nPedFitGood);
  ret.setInt("nPedFitFail", nPedFitFail);
  ret.setIntVector("fitStats", fitStats);
  ret.setIntVector("fitStats", fitStats);
  ret.setFloatVector("fitPedestals", fitPedestals);
  ret.setFloatVector("fitPedestalRmss", fitPedestalRmss);
  return ret;
}

//**********************************************************************

DataMap AdcPedestalFitter::beginEvent(const DuneEventInfo& evi) const {
  const string myname = "AdcPedestalFitter::beginEvent: ";
  DataMap res;
  ++state().nevt;
  if ( evi.run == state().run ) return res;
  if ( m_LogLevel >= 2 ) {
    cout << myname << "Setting run " << evi.run << endl;
  }
  state().run = evi.run;
  string msg;
  if ( m_haveFormulaParams ) {
    RunData rdat;
    if ( m_prdtool == nullptr ) {
      msg = "WARNING: RunData tool not found. Using default parameters.";
    } else {
      rdat = m_prdtool->runData(evi.run);
    }
    if ( msg.size() == 0 && ! rdat.isValid() ) {
      msg = "WARNING: RunData not found. Using default parameters.";
    }
    if ( msg.size() ) {
      cout << myname << msg << endl;
      rdat.setGain(14.0);
      rdat.setShaping(2.0);
    }
    for ( auto& itr : m_tfs ) {
      ParFormula& form = *itr.second;
      form.unsetParValues();
      rdat.setFormulaPars(form);
      if ( form.unsetPars().size() ) {
        cout << "ERROR: Formula " << form.name() << " has unset parameters: [";
        bool first = true;
        for ( string spar : form.unsetPars() ) {
          if ( first ) first = false;
          else cout << ", ";
          cout << spar;
        }
        cout << "]" << endl;
      }
      if ( form.resetPars().size() ) {
        cout << "WARNING: Formula " << form.name() << " has reset parameters" << endl;
      }
    }
  }
  return res;
}

//**********************************************************************

string AdcPedestalFitter::
nameReplace(string name, const AdcChannelData& acd, bool isTitle) const {
  const AdcChannelStringTool* pnbl = m_adcStringBuilder;
  if ( pnbl == nullptr ) return name;
  DataMap dm;
  return pnbl->build(acd, dm, name);
}

//**********************************************************************

DataMap
AdcPedestalFitter::getPedestal(const AdcChannelData& acd) const {
  const string myname = "AdcPedestalFitter::getPedestal: ";
  DataMap res;
  if ( m_LogLevel >= 2 ) cout << myname << "Fitting pedestal for channel " << acd.channel() << endl;
  string hnameBase = m_HistName;
  string htitlBase = m_HistTitle;
  Index nsam = acd.raw.size();
  if ( nsam == 0 ) {
    if ( m_LogLevel >= 2 ) cout << myname << "WARNING: Raw data is empty." << endl;
    return res.setStatus(1);
  }
  // If we need RunData parameters and Make sure that beginEvent has been called to
  // evaluate formulas.
  if ( m_haveFormulaParams && state().nevt == 0 ) {
    cout << myname << "WARNING: Calling beginEvent to evaluate formulas." << endl;
    beginEvent(acd.getEventInfo());
  }
  // Flag samples to keep in pedestal fit.
  //return res;
  vector<bool> keep(nsam, true);
  Index nkeep = 0;
  Index nskip = 0;
  for ( Index isam=0; isam<nsam; ++isam ) {
    if ( isam >= acd.flags.size() ) {
      if ( m_LogLevel >= 2 ) cout << myname << "WARNING: flags are missing." << endl;
      nkeep += nsam - isam;    // 2021-04-30: Keep unflagged samples
      break;
    }
    Index flg = acd.flags[isam];
    if ( m_skipFlags.count(flg) ) {
      keep[isam] = false;
      ++nskip;
    } else {
      ++nkeep;
    }
  }
  if ( nkeep == 0 ) {
    if ( m_LogLevel >= 2 ) cout << myname << "WARNING: No raw data is selected." << endl;
    return res.setStatus(2);
  }
  string hname = nameReplace(hnameBase, acd, false);
  string htitl = nameReplace(htitlBase, acd, true);
  htitl += "; ADC count; # samples";
  Index adcRange = m_tfs.at("AdcRange")->eval();
  if ( state().pfitter == nullptr || adcRange != state().adcRange ) {
    cout << myname << "New ADC range is " << adcRange << endl;
    state().pfitter = new TF1("pedgaus", "gaus", 0, adcRange, TF1::EAddToList::kNo);
  }
  state().adcRange = adcRange;
  unsigned int nadc = adcRange;
  unsigned int rebin = 10;
  unsigned int nbin = (nadc + rebin - 0.01)/rebin;
  double xmax = rebin*nbin;
  string hnamr = hname + "_rebin";
  TH1* phr = new TH1F(hnamr.c_str(), htitl.c_str(), nbin, 0, xmax);
  phr->SetDirectory(0);
  IndexVector rcounts(nbin, 0);
  for ( Index isam=0; isam<nsam; ++isam ) {
    if ( keep[isam] ) {
      Index iadc = acd.raw[isam];
      Index ibin = iadc/10;
      if ( ibin >= nbin ) cout << myname << "ERROR: Too many ADC counts for channel " << acd.channel()
                               << " sample " << isam << ": " << iadc << endl;
      else ++rcounts[ibin];
    }
  }
  for ( Index ibin=0; ibin<nbin; ++ibin ) {
    phr->SetBinContent(ibin+1, rcounts[ibin]);
  }
  double wadc = m_tfs.at("AdcFitRange")->eval();
  if ( wadc < 10 ) {
    cout << myname << "INFO: Invalid fit range: " << wadc << " < 10 ADC counts" << endl;
    return res.setStatus(3);
  }
  if ( m_LogLevel >= 4 ) cout << myname << "INFO: Width = " << wadc << " ADC counts" << endl;
  int rbinmax1 = phr->GetMaximumBin();
  double adcmax = phr->GetBinCenter(rbinmax1);
  double adc1 = adcmax - 0.5*wadc;
  double adc2 = adc1 + wadc;
  if ( m_RemoveStickyCode ) {
    double radcmax1 = phr->GetBinCenter(rbinmax1);
    double radcmean = phr->GetMean();
    double radcsum1 = phr->Integral();
    // Max may be due to a sticky code. Reduce it and find the next maximum.
    double tmpval = 0.5*(phr->GetBinContent(rbinmax1-1)+phr->GetBinContent(rbinmax1+1));
    phr->SetBinContent(rbinmax1, tmpval);
    int rbinmax2 = phr->GetMaximumBin();
    double radcsum2 = phr->Integral();
    // Evaluate the histogram mean and peak position.
    // If the peak removal has not removed too much data, these values are
    // re-evaluated using the peak-removed histogram.
    double adcmean = radcmean;        // Mean position.
    int rbinmax = rbinmax1;           // Peak position.
    if ( radcsum2 > 0.01*radcsum1 ) {
      // Define the max to be the first value if the two maxima are close or the
      // average if they are far part.
      if ( abs(rbinmax2-rbinmax1) > 1 ) {
        rbinmax = (rbinmax1 + rbinmax2)/2;
        adcmean = phr->GetMean();
      }
    }
    adcmax = phr->GetBinCenter(rbinmax);
    // Make sure the peak bin stays in range.
    if ( abs(adcmax-radcmax1) > 0.45*wadc ) adcmax = radcmax1;
    adc1 = adcmax - 0.5*wadc;
    adc1 = 10*int(adc1/10);
    if ( adcmean > adcmax + 10) adc1 += 10;
    adc2 = adc1 + wadc;
    if ( radcmax1 < adc1 || radcmax1+1.0 > adc2 ) {
      cout << myname << "WARNING: Histogram range (" << adc1 << ", " << adc2
           << ") does not include peak at " << radcmax1 << "." << endl;
    }
  }
  delete phr;
  TH1* phf = new TH1F(hname.c_str(), htitl.c_str(), wadc, adc1, adc2);
  phf->SetDirectory(0);
  //phf->Sumw2();
  Index countLo = 0;
  Index countHi = 0;
  Index count = 0;
  IndexVector fcounts(wadc, 0);
  for ( Index isam=0; isam<nsam; ++isam ) {
    AdcIndex val = acd.raw[isam];
    ++count;
    if ( val < adc1 ) {
      ++countLo;
    } else if ( val >= adc2 ) {
      ++countHi;
    } else {
      if ( keep[isam] ) ++fcounts[val-adc1];
    }
  }
  for ( Index iadc=0; iadc<wadc; ++iadc ) {
    phf->SetBinContent(iadc+1, fcounts[iadc]);
    //float err = sqrt(fcounts[iadc]);
    //if ( err <= 0.0 ) err = 1.0;
    //phf->SetBinError(iadc+1, err);
  }
  float fracLo = count > 0 ? float(countLo)/count : 0.0;
  float fracHi = count > 0 ? float(countHi)/count : 0.0;
  phf->SetStats(0);
  phf->SetLineWidth(2);
  // Fetch the peak bin and suppress it for the fit if more than 20% (but not
  // all) the data is in it.
  // May 2021: We expect the peak bin to hold 20% of the samples if the RMS is near
  // or below 2.0. Don't drop bin for low RMS and log a warning if bin is dropped
  // and fit sigm < 2.5.
  int binmax = phf->GetMaximumBin();
  double valmax = phf->GetBinContent(binmax);
  float peakFrac = valmax/phf->Integral();
  double xcomax = phf->GetBinLowEdge(binmax);
  double rangeIntegral = phf->Integral(1, phf->GetNbinsX());
  double peakBinFraction = (rangeIntegral > 0) ? valmax/rangeIntegral : 1.0;
  bool allBin = peakBinFraction > 0.99;
  float rawRms = phf->GetRMS();
  bool isNarrow = rawRms < 2.5;
  bool dropBin = m_RemoveStickyCode && peakFrac > 0.2 && !allBin && !isNarrow;
  int nbinsRemoved = 0;
  if ( dropBin ) {
    double tmpval = 0.5*(phf->GetBinContent(binmax-1)+phf->GetBinContent(binmax+1));
    phf->SetBinContent(binmax, tmpval);
    rangeIntegral = phf->Integral(1, phf->GetNbinsX());
    ++nbinsRemoved;
  }
  double amean = phf->GetMean() + 0.5;
  double arms = phf->GetRMS();
  double fitRmsMin = m_tfs.at("FitRmsMin")->eval();
  double fitRmsMax = m_tfs.at("FitRmsMax")->eval();
  double ameanWin = fitRmsMax > fitRmsMin ? fitRmsMax : 0.0;
  bool doFit = peakBinFraction < 0.99 && m_fitOpts.size() > 0;
  float pedestal = phf->GetMean();
  //float pedestalRms = 1.0/sqrt(12.0);
  float pedestalRms = arms;
  if ( fitRmsMax > fitRmsMin ) {
    if ( arms < fitRmsMin ) arms = fitRmsMin;
    if ( arms > fitRmsMax ) arms = fitRmsMax;
  }
  float fitChiSquare = 0.0;
  float fitReducedChiSquare = 0.0;
  float peakBinExcess = 0.0;
  if ( doFit ) {
    // Block Root info message for new Canvas produced in fit.
    int levelSave = gErrorIgnoreLevel;
    gErrorIgnoreLevel = 1001;
    gErrorIgnoreLevel = 2001;   // Block warning in Fit
    // Block non-default (e.g. art) from handling the Root "error".
    // We switch to the Root default handler while making the call to Print.
    ErrorHandlerFunc_t pehSave = nullptr;
    ErrorHandlerFunc_t pehDefault = DefaultErrorHandler;
    //pehDefault = handleRootError;
    if ( GetErrorHandler() != pehDefault ) {
      pehSave = SetErrorHandler(pehDefault);
    }
    // Set fit precision. Least-squares only?
    float fprec = TFitter::GetPrecision();
    if ( m_FitPrecision > 0.0 ) TFitter::SetPrecision(m_FitPrecision);
    int fitStat = 99;
    TF1* pfinit = nullptr;
    //TF1 fitter("pedgaus", "gaus", adc1, adc2, TF1::EAddToList::kNo);
    TF1& fitter = *state().pfitter;
    for ( Name fopt : m_fitOpts ) {
      fitter.SetParameters(0.1*rangeIntegral, amean, pedestalRms);
      fitter.SetParLimits(0, 0.01*rangeIntegral, rangeIntegral);
      if ( ameanWin > 0.0 ) fitter.SetParLimits(1, amean - ameanWin, amean + ameanWin);
      if ( allBin ) fitter.FixParameter(1, amean);  // Fix posn.
      if ( fitRmsMin < fitRmsMax ) {
        fitter.SetParLimits(2, fitRmsMin, fitRmsMax);
      }
      fitter.SetParError(0, 0.0);
      if ( m_PlotShowFit >= 2 ) {
        pfinit = dynamic_cast<TF1*>(fitter.Clone("pedgaus0"));
        pfinit->SetLineColor(3);
        pfinit->SetLineStyle(2);
      }
      if ( m_LogLevel < 3 ) fopt += "Q";
      // Root calls error handler but returns 0 if the histo has no data in range
      // so we check that before calling fitter.
      fitStat = (phf->Integral() == 0.0) ? 999 : int(phf->Fit(&fitter, fopt.c_str(), "", adc1, adc2));
      //fopt += " S";  // So we can retrieve fit status
      //TFitResultPtr pfres = phf->Fit(&fitter, fopt.c_str());
      //bool haveFitStatus = pfres.Get();
      //int fitStat = haveFitStatus ? pfres->Status() : 999;
      if ( fitStat == 0 ) break;
    }
    if ( m_FitPrecision > 0.0 ) TFitter::SetPrecision(fprec);
    if ( pehSave != nullptr ) SetErrorHandler(pehSave);
    gErrorIgnoreLevel = levelSave;
    if ( pfinit != nullptr ) {
      phf->GetListOfFunctions()->AddLast(pfinit, "0");
      phf->GetListOfFunctions()->Last()->SetBit(TF1::kNotDraw, true);
    }
    if ( fitStat ) {
      Index istat = acd.channelStatus();
      string sstat = istat == 0 ? "good" : istat == 1 ? "bad" : istat == 2 ? "noisy" : "unknown";
      cout << myname << "WARNING: Fit status is " << fitStat << " for " << sstat << " channel "
           << acd.channel() << endl;
      cout << myname << "  Errors[0]: " << fitter.GetParErrors()[0] << endl;
      //cout << myname << "  radcmax1 = " << radcmax1 << endl;
      //cout << myname << "  radcmean = " << radcmean << endl;
      //cout << myname << "   adcmean = " << adcmean << endl;
      cout << myname << "    adcmax = " << adcmax << endl;
      //cout << myname << "  rbinmax1,2: " << rbinmax1 << ", " << rbinmax2 << endl;
      cout << myname << "  peakBinFraction = " << valmax << "/" << rangeIntegral
           << " = " << peakBinFraction << endl;
      cout << myname << "  allBin = " << allBin << endl;
      cout << myname << "  dropBin = " << dropBin << endl;
      cout << myname << "  amean = " << amean << " +/- " << ameanWin << endl;
    } else {
      double valEval = fitter.Eval(xcomax);
      peakBinExcess = (valmax - valEval)/rangeIntegral;
      pedestal = fitter.GetParameter(1) - 0.5;
      pedestalRms = fitter.GetParameter(2);
      // Chi-square from the fitter.
      fitChiSquare = fitter.GetChisquare();
      // Estimate of chi-square/DOF
      float ibin1 = int(pedestal - 3.0*pedestalRms) - adc1;
      if ( ibin1 < 1 ) ibin1 = 1;
      float ibin2 = int(pedestal + 3.0*pedestalRms + 2.0) - adc1;
      if ( ibin2 > wadc ) ibin2 = wadc;
      int nbin = 0;
      double sumsq = 0.0;
      for ( int ibin=ibin1; ibin<=ibin2; ++ibin ) {
        float xfun = adc1 + ibin - 0.5;
        float yfun = fitter.Eval(xfun);
        if ( yfun < 1.0 ) continue;
        float yhst = phf->GetBinContent(ibin);
        float dely = yhst - yfun;
        float vary = yfun > 1.0 ? yfun : 1.0;
        sumsq += dely*dely/vary;
        ++nbin;
      }
      float estndof = nbin > 4 ? nbin - 3 : 1;
      fitReducedChiSquare = sumsq/estndof;
      // Warn if bin was dropped for a too-narrow distribution.
      if ( dropBin && pedestalRms < 2.5 ) {
        int precSave = cout.precision();
        cout.precision(2);
        cout << myname << "WARNING: Bin was dropped for channel " << acd.channel()
             << " with fit sigma=" << fixed << pedestalRms
             << ", raw RMS=" << fixed << rawRms
             << " and peak fraction=" << fixed << peakFrac << endl;
        cout.precision(precSave);
      }
    }
  }
  if ( dropBin ) phf->SetBinContent(binmax, valmax);
  res.setHist("pedestal", phf, true);
  res.setFloat("fitFractionLow", fracLo);
  res.setFloat("fitFractionHigh", fracHi);
  res.setFloat("fitPedestal", pedestal);
  res.setFloat("fitPedestalRms", pedestalRms);
  res.setFloat("fitChiSquare", fitChiSquare);
  res.setFloat("fitPeakBinFraction", peakBinFraction);
  res.setFloat("fitPeakBinExcess", peakBinExcess);
  res.setFloat("fitChiSquare", fitChiSquare);
  res.setFloat("fitReducedChiSquare", fitReducedChiSquare);
  res.setInt("fitChannel", acd.channel());
  res.setInt("fitNSkip", nskip);
  res.setInt("fitNBinsRemoved", nbinsRemoved);
  string rfname = nameReplace(m_RootFileName, acd, false);
  if ( rfname.size() ) {
    if ( m_LogLevel >=2 ) cout << myname << "Write histogram " << phf->GetName() << " to " << rfname << endl;
    TFile* pf = TFile::Open(rfname.c_str(), "UPDATE");
    TH1* phfCopy = (TH1*) phf->Clone();
    phfCopy->GetListOfFunctions()->Clear();
    phfCopy->Write();
    delete pf;
  }
  if ( m_LogLevel >= 3 ) cout << myname << "Exiting..." << endl;
  return res;
}

//**********************************************************************

int AdcPedestalFitter::fillChannelPad(DataMap& dm, const AdcChannelData& acd, TPadManipulator* pman) const {
  if ( pman == nullptr ) return 1;
  TH1* phf = dm.getHist("pedestal");
  pman->add(phf, "hist");
  if ( m_PlotShowFit > 1 ) pman->addHistFun(1);
  if ( m_PlotShowFit ) pman->addHistFun(0);
  pman->addVerticalModLines(64);
  pman->showUnderflow();
  pman->showOverflow();
  pman->addAxis();
  Index istat = acd.channelStatus();
  string sstat = istat == 0 ? "Good" : istat == 1 ? "Bad" : istat == 2 ? "Noisy" : "No-status";
  NameVector slabs(3);
  slabs[0] = sstat + " channel";
  slabs[1] = "# skipped bins: " + std::to_string(dm.getInt("fitNSkip"));
  slabs[2] = "# dropped bins: " + std::to_string(dm.getInt("fitNBinsRemoved"));
  float xlab = 0.13;
  float ylab = 0.86;
  float dylab = 0.06;
  for ( Name slab : slabs ) {
    TLatex* ptxt = new TLatex(xlab, ylab, slab.c_str());
    ptxt->SetNDC();
    ptxt->SetTextFont(42);
    pman->add(ptxt);
    ylab -= dylab;
  }
  slabs.clear();
  ostringstream sslab;
  sslab << "Pedestal: " << std::fixed << std::setprecision(1) << dm.getFloat("fitPedestal");
  slabs.push_back(sslab.str());
  sslab.str("");
  sslab << "Ped RMS: " << std::fixed << std::setprecision(1) << dm.getFloat("fitPedestalRms");
  slabs.push_back(sslab.str());
  sslab.str("");
  //sslab << "Fit #chi^{2}: " << std::fixed << std::setprecision(1) << dm.getFloat("fitChiSquare");
  //slabs.push_back(sslab.str());
  //sslab.str("");
  sslab << "#chi^{2}/DOF: " << std::fixed << std::setprecision(1) << dm.getFloat("fitReducedChiSquare");
  slabs.push_back(sslab.str());
  sslab.str("");
  float frac = dm.getFloat("fitFractionLow");
  int prec = frac > 0 ? int(-log10(frac)) + 3 : 0;
  sslab << "f_{lo}: " << std::fixed << std::setprecision(prec) << frac;
  slabs.push_back(sslab.str());
  sslab.str("");
  frac = dm.getFloat("fitFractionHigh");
  prec = frac > 0 ? int(-log10(frac)) + 3 : 0;
  sslab << "f_{hi}: " << std::fixed << std::setprecision(prec) << frac;
  slabs.push_back(sslab.str());
  xlab = 0.94;
  ylab = 0.86;
  for ( Name slab : slabs ) {
    TLatex* ptxt = new TLatex(xlab, ylab, slab.c_str());
    ptxt->SetNDC();
    ptxt->SetTextFont(42);
    ptxt->SetTextAlign(31);
    pman->add(ptxt);
    ylab -= dylab;
  }
  return 0;
}

//**********************************************************************

DEFINE_ART_CLASS_TOOL(AdcPedestalFitter)
