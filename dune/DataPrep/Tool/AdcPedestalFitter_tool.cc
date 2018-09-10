// AdcPedestalFitter_tool.cc

#include "AdcPedestalFitter.h"
#include "dune/DuneCommon/TPadManipulator.h"
#include "dune/DuneCommon/StringManipulator.h"
#include "dune/DuneInterface/Tool/AdcChannelStringTool.h"
#include "dune/DuneInterface/Tool/HistogramManager.h"
#include "dune/ArtSupport/DuneToolManager.h"
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
// Local definitions.
//**********************************************************************

namespace {

void copyMetadata(const DataMap& res, AdcChannelData& acd) {
  acd.metadata["fitPedFractionLow"] = res.getFloat("fitFractionLow");
  acd.metadata["fitPedFractionHigh"] = res.getFloat("fitFractionHigh");
  acd.metadata["fitPedestal"] = res.getFloat("fitPedestal");
  acd.metadata["fitPedRms"] = res.getFloat("fitPedestalRms");
  acd.metadata["fitPedChiSquare"] = res.getFloat("fitChiSquare");
  acd.metadata["fitPedPeakBinFraction"] = res.getFloat("fitPeakBinFraction");
  acd.metadata["fitPedPeakBinExcess"] = res.getFloat("fitPeakBinExcess");
  acd.metadata["fitPedNBinsRemoved"] = res.getFloat("fitNBinsRemoved");
}

}  // end unnamed namespace

//**********************************************************************
// Class methods.
//**********************************************************************

AdcPedestalFitter::AdcPedestalFitter(fhicl::ParameterSet const& ps)
: m_LogLevel(ps.get<int>("LogLevel")),
  m_FitRmsMin(ps.get<float>("FitRmsMin")),
  m_FitRmsMax(ps.get<float>("FitRmsMax")),
  m_HistName(ps.get<string>("HistName")),
  m_HistTitle(ps.get<string>("HistTitle")),
  m_HistManager(ps.get<string>("HistManager")),
  m_PlotFileName(ps.get<string>("PlotFileName")),
  m_RootFileName(ps.get<string>("RootFileName")),
  m_PlotSizeX(ps.get<Index>("PlotSizeX")),
  m_PlotSizeY(ps.get<Index>("PlotSizeY")),
  m_PlotShowFit(ps.get<Index>("PlotShowFit")),
  m_PlotSplitX(ps.get<Index>("PlotSplitX")),
  m_PlotSplitY(ps.get<Index>("PlotSplitY")),
  m_phm(nullptr) {
  const string myname = "AdcPedestalFitter::ctor: ";
  DuneToolManager* ptm = DuneToolManager::instance();
  if ( m_HistManager.size() ) {
    m_phm = ptm->getShared<HistogramManager>(m_HistManager);
    if ( m_phm == nullptr ) {
      cout << myname << "WARNING: Histogram manager not found: " << m_HistManager << endl;
    }
  }
  string stringBuilder = "adcStringBuilder";
  m_adcStringBuilder = ptm->getShared<AdcChannelStringTool>(stringBuilder);
  if ( m_adcStringBuilder == nullptr ) {
    cout << myname << "WARNING: AdcChannelStringTool not found: " << stringBuilder << endl;
  }
  if ( m_LogLevel >= 1 ) {
    cout << myname << "Configuration parameters:" << endl;
    cout << myname << "       LogLevel: " << m_LogLevel << endl;
    cout << myname << "       HistName: " << m_HistName << endl;
    cout << myname << "      HistTitle: " << m_HistTitle << endl;
    cout << myname << "   PlotFileName: " << m_PlotFileName << endl;
    cout << myname << "   RootFileName: " << m_RootFileName << endl;
    cout << myname << "    HistManager: " << m_HistManager << endl;
    cout << myname << "      PlotSizeX: " << m_PlotSizeX << endl;
    cout << myname << "      PlotSizeY: " << m_PlotSizeY << endl;
    cout << myname << "    PlotShowFit: " << m_PlotShowFit << endl;
    cout << myname << "     PlotSplitX: " << m_PlotSplitX << endl;
    cout << myname << "     PlotSplitY: " << m_PlotSplitY << endl;
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
  fillChannelPad(res, pman);
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
  if ( m_PlotFileName.size() && m_PlotSplitX > 0 ) {
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
    if ( m_LogLevel >= 3 ) cout << myname << "  " << iacd << ": Processing channel " << acd.channel << endl;
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
    fillChannelPad(tmpres, pman);
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
  if ( m_LogLevel >= 2 ) cout << myname << "Fitting pedestal for channel " << acd.channel << endl;
  string hnameBase = m_HistName;
  string htitlBase = m_HistTitle;
  Index nsam = acd.raw.size();
  if ( nsam == 0 ) {
    if ( m_LogLevel >= 2 ) cout << myname << "WARNING: Raw data is empty." << endl;
    return res.setStatus(1);
  }
  string hname = nameReplace(hnameBase, acd, false);
  string htitl = nameReplace(htitlBase, acd, true);
  string pfname = nameReplace(m_PlotFileName, acd, false);
  string rfname = nameReplace(m_RootFileName, acd, false);
  htitl += "; ADC count; # samples";
  unsigned int nadc = 4096;
  unsigned int rebin = 10;
  unsigned int nbin = (nadc + rebin - 0.01)/rebin;
  double xmax = rebin*nbin;
  string hnamr = hname + "_rebin";
  TH1* phr = new TH1F(hname.c_str(), htitl.c_str(), nbin, 0, xmax);
  phr->SetDirectory(0);
  for ( Index isam=0; isam<nsam; ++isam ) {
    phr->Fill(acd.raw[isam]);
  }
  int rbinmax1 = phr->GetMaximumBin();
  double radcmax1 = phr->GetBinCenter(rbinmax1);
  double adcmean = phr->GetMean();
  // Max may just be a sticky code. Reduce it and find the next maximum.
  double tmpval = 0.5*(phr->GetBinContent(rbinmax1-1)+phr->GetBinContent(rbinmax1+1));
  phr->SetBinContent(rbinmax1, tmpval);
  int rbinmax2 = phr->GetMaximumBin();
  // Define the max to be the first value if the two maxima are close or the
  // average if they are far part.
  int rbinmax = rbinmax1;
  if ( abs(rbinmax2-rbinmax1) > 1 ) {
    rbinmax = (rbinmax1 + rbinmax2)/2;
    adcmean = phr->GetMean();
  }
  double adcmax = phr->GetBinCenter(rbinmax);
  delete phr;
  double wadc = 100.0;
  // Make sure the peak bin stays in range.
  if ( abs(adcmax-radcmax1) > 0.5*wadc ) adcmax = radcmax1;
  double adc1 = adcmax - 0.5*wadc;
  adc1 = 10*int(adc1/10);
  if ( adcmean > adcmax + 10) adc1 += 10;
  double adc2 = adc1 + wadc;
  TH1* phf = new TH1F(hname.c_str(), htitl.c_str(), wadc, adc1, adc2);
  phf->SetDirectory(nullptr);
  Index countLo = 0;
  Index countHi = 0;
  Index count = 0;
  for ( Index isam=0; isam<nsam; ++isam ) {
    AdcIndex val = acd.raw[isam];
    if ( val < adc1 ) ++countLo;
    if ( val >= adc2 ) ++countHi;
    ++count;
    phf->Fill(acd.raw[isam]);
  }
  float fracLo = count > 0 ? float(countLo)/count : 0.0;
  float fracHi = count > 0 ? float(countHi)/count : 0.0;
  phf->SetStats(0);
  phf->SetLineWidth(2);
  // Fetch the peak bin and suppress it for the fit if more than 20% (but not
  // all) the data is in it.
  int binmax = phf->GetMaximumBin();
  double valmax = phf->GetBinContent(binmax);
  double xcomax = phf->GetBinLowEdge(binmax);
  double rangeIntegral = phf->Integral(1, phf->GetNbinsX());
  double peakBinFraction = valmax/rangeIntegral;
  bool allBin = peakBinFraction > 0.99;
  bool dropBin = valmax > 0.2*phf->Integral() && !allBin;
  int nbinsRemoved = 0;
  if ( dropBin ) {
    double tmpval = 0.5*(phf->GetBinContent(binmax-1)+phf->GetBinContent(binmax+1));
    phf->SetBinContent(binmax, tmpval);
    rangeIntegral = phf->Integral(1, phf->GetNbinsX());
    ++nbinsRemoved;
  }
  double amean = phf->GetMean() + 0.5;
  double alim1 = amean - 25.0;
  double alim2 = amean + 25.0;
  TF1 fitter("pedgaus", "gaus", adc1, adc2, TF1::EAddToList::kNo);
  fitter.SetParameters(0.1*rangeIntegral, amean, 5.0);
  fitter.SetParLimits(0, 0.01*rangeIntegral, rangeIntegral);
  fitter.SetParLimits(1, alim1, alim2);
  if ( allBin ) fitter.FixParameter(1, amean);  // Fix posn.
  if ( m_FitRmsMin < m_FitRmsMax ) {
    fitter.SetParLimits(2, m_FitRmsMin, m_FitRmsMax);
  }
  TF1* pfinit = dynamic_cast<TF1*>(fitter.Clone("pedgaus0"));
  pfinit->SetLineColor(3);
  pfinit->SetLineStyle(2);
  string fopt = "0";
  fopt = "WWB";
  if ( m_LogLevel < 3 ) fopt += "Q";
  // Block Root info message for new Canvas produced in fit.
  int levelSave = gErrorIgnoreLevel;
  gErrorIgnoreLevel = 1001;
  // Block non-default (e.g. art) from handling the Root "error".
  // We switch to the Root default handler while making the call to Print.
  ErrorHandlerFunc_t pehSave = nullptr;
  ErrorHandlerFunc_t pehDefault = DefaultErrorHandler;
  if ( GetErrorHandler() != pehDefault ) {
    pehSave = SetErrorHandler(pehDefault);
  }
  phf->Fit(&fitter, fopt.c_str());
  if ( pehSave != nullptr ) SetErrorHandler(pehSave);
  gErrorIgnoreLevel = levelSave;
  phf->GetListOfFunctions()->AddLast(pfinit, "0");
  phf->GetListOfFunctions()->Last()->SetBit(TF1::kNotDraw, true);
  double valEval = fitter.Eval(xcomax);
  double peakBinExcess = (valmax - valEval)/rangeIntegral;
  if ( dropBin ) phf->SetBinContent(binmax, valmax);
  res.setHist("pedestal", phf, true);
  res.setFloat("fitFractionLow", fracLo);
  res.setFloat("fitFractionHigh", fracHi);
  res.setFloat("fitPedestal", fitter.GetParameter(1) - 0.5);
  res.setFloat("fitPedestalRms", fitter.GetParameter(2));
  res.setFloat("fitChiSquare", fitter.GetChisquare());
  res.setFloat("fitPeakBinFraction", peakBinFraction);
  res.setFloat("fitPeakBinExcess", peakBinExcess);
  res.setInt("fitChannel", acd.channel);
  res.setInt("fitNBinsRemoved", nbinsRemoved);
/*
  if ( pman != nullptr ) {
    pman->add(phf, "hist", false);
    if ( m_PlotShowFit > 1 ) pman->addHistFun(1);
    if ( m_PlotShowFit ) pman->addHistFun(0);
    pman->addVerticalModLines(64);
    pman->showUnderflow();
    pman->showOverflow();
  }
*/
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

int AdcPedestalFitter::fillChannelPad(DataMap& dm, TPadManipulator* pman) const {
  if ( pman == nullptr ) return 1;
  TH1* phf = dm.getHist("pedestal");
  pman->add(phf, "hist", false);
  if ( m_PlotShowFit > 1 ) pman->addHistFun(1);
  if ( m_PlotShowFit ) pman->addHistFun(0);
  pman->addVerticalModLines(64);
  pman->showUnderflow();
  pman->showOverflow();
  return 0;
}

//**********************************************************************
