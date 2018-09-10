// StickyCodeMetrics.cxx

#include "StickyCodeMetrics.h"
#include <map>
#include <iostream>
#include <sstream>
#include <vector>

#include "TH1F.h"
#include "TF1.h"
//#include "TFitResult.h"

using std::map;
using std::string;
using std::cout;
using std::endl;
using std::ostringstream;
using std::vector;

namespace {
using Index = unsigned int;
}

//**********************************************************************

StickyCodeMetrics::
StickyCodeMetrics(Name hnam, Name httl, Index nbin, Index lowbin,
                  float sigmaMin, float sigmaMax)
: m_hnam(hnam), m_httl(httl), m_nbin(nbin), m_lowbin(lowbin),
  m_sigmaMin(sigmaMin), m_sigmaMax(sigmaMax) { }

//**********************************************************************

int StickyCodeMetrics::evaluate(const BinCounter& counts) {
  m_counts = counts;
  return evaluateMetrics();
}

//**********************************************************************

int StickyCodeMetrics::evaluate(const AdcCountVector& adcs) {
  m_counts.clear();
  for ( AdcCount adc : adcs ) {
    if ( m_counts.find(adc) == m_counts.end() ) m_counts[adc] = 0;
    ++m_counts[adc];
  }
  return evaluateMetrics();
}

//**********************************************************************

int StickyCodeMetrics::evaluate(const TH1* pha) {
  const string myname = "StickyCodeMetrics::evaluate: ";
  Index nbin = pha->GetNbinsX();
  Index iadc0 = pha->GetBinLowEdge(1);
  Index iadcLast = pha->GetBinLowEdge(nbin);
  Index nhadc = iadcLast - iadc0 + 1;
  m_counts.clear();
  if ( nhadc != nbin ) {
    cout << myname << "Histogram " << pha->GetName() << " has inconsistent binning." << endl;
  } else {
    for ( Index ibin=0; ibin<nbin; ++ibin ) {
      double count = pha->GetBinContent(ibin+1);
      if ( count != 0.0 ) m_counts[iadc0+ibin] = count;
    }
  }
  return evaluateMetrics();
}

//**********************************************************************

void StickyCodeMetrics::clear() {
  m_counts.clear();
  evaluateMetrics();
}

//**********************************************************************

int StickyCodeMetrics::evaluateMetrics() {
  const string myname = "StickyCodeMetrics::evaluateMetrics: ";
  const int dbg = 0;  // Nonzero give debug logging.
  if ( dbg ) cout << myname << "Debugging level: " << dbg << endl;
  const BinCounter& counts = m_counts;
  m_nsample = 0;
  Index badadc = 999999;
  Index maxadc = badadc;
  Index maxCount = 0;
  Index maxCount2 = 0;
  Index nmod0 = 0;
  Index nmod1 = 0;
  Index nmod63 = 0;
  double adcSum = 0.0;
  double countSum = 0.0;
  for ( BinCounter::value_type ibco : counts ) {
    Index iadc = ibco.first;
    double count = ibco.second;
    m_nsample += count;
    if ( count > maxCount ) {
      maxCount = count;
      maxadc = iadc;
    }
    AdcCount adcmod = iadc%64;
    if ( adcmod ==  0 ) nmod0 += count;
    if ( adcmod ==  1 ) nmod1 += count;
    if ( adcmod == 63 ) nmod63 += count;
    countSum += count;
    adcSum += count*iadc;
  }
  double countSum2 = 0;
  double adcSum2 = 0;
  Index maxadc2 = badadc;
  for ( BinCounter::value_type ibco : counts ) {
    Index iadc = ibco.first;
    double count = ibco.second;
    if ( iadc != maxadc ) {
      countSum2 += count;
      adcSum2 += count*iadc;
      if ( count > maxCount2 ) {
        maxCount2 = count;
        maxadc2 = iadc;
      }
    }
  }
  m_maxAdc = maxadc;
  m_maxAdc2 = maxadc2;
  m_meanAdc = countSum>0 ? adcSum/countSum : -1.0;
  m_meanAdc2 = countSum2>0 ? adcSum2/countSum2 : -1.0;
  m_maxFraction = maxCount/countSum;
  m_zeroFraction = nmod0/countSum;
  m_oneFraction = nmod1/countSum;
  m_highFraction = nmod63/countSum;
  // Create the histogram if configured.
  if ( counts.size() == 0 ) return 0;
  if ( m_hnam.size() == 0 ) return 0;
  if ( m_nbin == 0 ) return 1;
  if ( m_lowbin == 0 ) return 2;
  Index nbinHist = m_nbin;           // # bins in the histogram
  // Create vectors with counts and rebinned counts.
  // binCounts[ibin] holds the count for iadc = binOffset + ibin
  // rebinCounts[ireb] holds the count sum for nrebin bins starting at
  // iadc = binOffset + ireb*rebin
  // Not use of fabs in case someone passes negative values.
  Index idadc1 = counts.begin()->first;
  Index idadc2 = counts.rbegin()->first;
  Index nrebin = m_lowbin;
  vector<double> binCounts;
  // Pad the count vector with half the histo range so we can center the
  // data in the histogram.
  Index halfHist = nbinHist/2;
  Index iadc1 = idadc1 < halfHist ? 0 : idadc1 - halfHist;
  iadc1 = nrebin*(iadc1/nrebin);
  Index binOffset = iadc1;
  Index iadc2= idadc2 + halfHist;
  if ( iadc2 > 4096 && idadc2 < 4096 ) iadc2 = 4096;
  countSum = 0.0;
  for ( Index iadc=iadc1; iadc<=iadc2; ++iadc ) {
    BinCounter::const_iterator icnt = counts.find(iadc);
    double count = icnt != counts.end() ? fabs(icnt->second) : 0;
    countSum += count;
    binCounts.push_back(count);
  }
  Index nbin = binCounts.size();
  Index nreb = (nbin+nrebin-1)/nrebin;
  vector<double> rebinCounts(nreb, 0.0);
  vector<double> rebinAdcCounts(nreb, 0.0);
  for ( Index ibin=0; ibin<nbin; ++ibin ) {
    double count = binCounts[ibin];
    Index iadc = binOffset + ibin;
    if ( count > 0 ) {
      Index ireb = ibin/nrebin;
      rebinCounts[ireb] += count;
      rebinAdcCounts[ireb] += iadc*count;
    }
  }
  if ( dbg >= 2 ) {
    cout << myname << "   Data ADC range: [" << idadc1 << ", " << idadc2 << "]" << endl;
    cout << myname << "       Bin offset: " << binOffset << endl;
    cout << myname << "   Data bin count: " << nbin << endl;
    cout << myname << " Data rebin count: " << rebinCounts.size() << endl;
    cout << myname << "   Hist bin count: " << nbinHist << endl;
  }
  // Initialize the histogram under and overflows.
  double countUnder = 0.0;
  double countOver = 0.0;
  // Evaluate the hist ADC range [ihadc1, ihadc2)
  // First check if all data fit in histogram.
  Index ihadc1 = binOffset;
  while ( ihadc1 + nrebin <= idadc1 ) ihadc1 += nrebin;
  Index ihadc2 = ihadc1 + nbinHist;
  // Data fits in histogram range.
  if ( ihadc2 > idadc2 ) {
    Index nbinLeft = idadc1 - ihadc1;
    Index nbinRight = ihadc2 - idadc2 - 1;
    if ( dbg >= 3 ) cout << myname << "All data are in histogram range." << endl;
    if ( dbg >= 4 ) cout << myname << "nl, nr, hleft: " << nbinLeft << ", "
                         << nbinRight << ", " << ihadc1 << endl;
    while ( nbinRight > nbinLeft + nrebin ) {
      ihadc1 -= nrebin;
      ihadc2 -= nrebin;
      nbinRight -= nrebin;
      nbinLeft += nrebin;
      if ( dbg >= 4 ) cout << myname << "nl, nr, hleft: " << nbinLeft << ", "
                           << nbinRight << ", " << ihadc1 << endl;
      if ( ihadc1 < nrebin ) break;
    }
  // Data extends beyond histogram range.
  } else { 
    ihadc1 = binOffset;
    ihadc2 = ihadc1 + nbinHist;
    if ( dbg >= 3 ) cout << myname << "Some data are outside histogram range." << endl;
    // Preferred ADC count for the center of the histogram.
    double adcHistCenter = maxAdc();
    if ( fabs(maxAdc2() - maxAdc()) > 0.1*nbinHist ) {
      adcHistCenter = 0.5*(maxAdc2() + maxAdc());
    }
    if ( dbg >= 3 ) {
      cout << myname << "   max ADC: " << maxAdc() << endl;
      cout << myname << "  max ADC2: " << maxAdc() << endl;
      cout << myname << "  mean ADC: " << meanAdc() << endl;
      cout << myname << "  Preferred hist center: " << adcHistCenter << endl;
    }
    // Get the total counts in the starting histogram range.
    double countSumHist = 0.0;
    double adcSumHist = 0.0;
    for ( Index iadc=ihadc1; iadc<ihadc2; ++iadc ) {
      Index ibin = iadc - binOffset;
      countSumHist += binCounts[ibin];
      adcSumHist += binCounts[ibin]*iadc;
    }
    // Find the value of ihadc1 that gives the maximum area.
    // If areas are similar, center the distribution.
    Index ihadc1Selected = ihadc1;
    double countSumHistSelected = countSumHist;
    double hmean = countSumHist == 0.0 ? 0.0 : adcSumHist/countSumHist;
    double hmeanSelected = hmean;
    if ( dbg >= 4 ) cout << "  hleft, hmean, area: " << ihadc1 << ", " << hmean
                         << ", " << countSumHist << endl;
    while ( ihadc1 <= idadc2 ) {
      // Add the new area to the right.
      for ( Index iadc=ihadc2; iadc<ihadc2+nrebin; ++iadc ) {
        if ( iadc > idadc2 ) break;
        Index ibin = iadc - binOffset;
        countSumHist += binCounts[ibin];
        adcSumHist += iadc*binCounts[ibin];
      }
      // Subtract the area to the left.
      Index ireb = (ihadc1-binOffset)/nrebin;
      if ( dbg >= 4 ) cout << "  ireb/nreb: " << ireb << "/" << nrebin << endl;
      countSumHist -= rebinCounts[ireb];
      adcSumHist -= rebinAdcCounts[ireb];
      hmean = countSumHist == 0.0 ? 0.0 : adcSumHist/countSumHist;
      // Update the range.
      ihadc1 += nrebin;
      ihadc2 += nrebin;
      if ( dbg >= 4 ) cout << "  hleft, hmean, area: " << ihadc1 << ", " << hmean
                           << ", " << countSumHist << endl;
      // Update the histogram range if this area is bigger or
      // if it is about the same and is closer to the mean.
      double diffCount = countSumHist - countSumHistSelected;
      bool sameCount = diffCount == 0.0;
      if ( ! sameCount ) {
        double rat = diffCount/(countSumHist + countSumHistSelected);
        sameCount = fabs(rat) < 0.001;
      }
      if ( ihadc1Selected > ihadc1 ) {
        cout << myname << "Unexpected value for ihadc1Selected: " << ihadc1Selected << endl;
        cout << myname << "                        with ihadc1: " << ihadc1 << endl;
        abort();
      }
      bool haveOverlap = ihadc1Selected + nbinHist > ihadc1;
      bool updateSel = false;
      if ( sameCount && haveOverlap ) {
        if ( dbg >= 4 ) cout << "  Same count." << endl;
        if ( true ) {
          double diffOld = fabs(ihadc1Selected + nbinHist/2 - hmeanSelected);
          double diffNew = fabs(ihadc1         + nbinHist/2 - hmean);
          updateSel = diffNew < diffOld;
        } else {
          double diffOld = fabs(ihadc1Selected + nbinHist/2 - adcHistCenter);
          double diffNew = fabs(ihadc1         + nbinHist/2 - adcHistCenter);
          updateSel = diffNew < diffOld;
        }
      } else {
        updateSel = diffCount > 0.0;
      }
      if ( updateSel ) {
        countSumHistSelected = countSumHist;
        hmeanSelected = hmean;
        ihadc1Selected = ihadc1;
        if ( dbg >= 4 ) cout << "  Updated selected hleft, hmean to " << ihadc1Selected
                             << ", " << hmeanSelected << endl;
      }
    }
    ihadc1 = ihadc1Selected;
    ihadc2 = ihadc1 + nbinHist;
    // Get the under and overflows.
    for ( BinCounter::value_type icnt : counts ) {
      Index iadc = icnt.first;
      double count = icnt.second;
      if ( iadc < ihadc1 ) countUnder += count;
      if ( iadc >= ihadc2 ) countOver += count;
    }
  }
  // Create histogram.
  TH1* ph = new TH1F(m_hnam.c_str(), m_httl.c_str(), nbinHist, ihadc1, ihadc2);
  string::size_type iposx = m_httl.find(";");
  bool haveXlab = iposx != string::npos;
  bool haveYlab = haveXlab && m_httl.find(";", iposx+1) != string::npos;
  if ( ! haveXlab ) ph->GetXaxis()->SetTitle("ADC count");
  if ( ! haveYlab ) ph->GetYaxis()->SetTitle("# samples");
  ph->SetDirectory(0);
  ph->SetStats(0);
  ph->SetLineWidth(2);
  for ( Index iadc=binOffset; iadc<=idadc2; ++iadc ) {
    BinCounter::const_iterator icnt = counts.find(iadc);
    Index count = icnt != counts.end() ? fabs(icnt->second) : 0;
    ph->SetBinContent(iadc+1-ihadc1, count);
  }
  ph->SetBinContent(0, countUnder);
  ph->SetBinContent(nbinHist+1, countOver);
  m_ph.reset(ph);
  // Fit the histogram.
  TF1 fitter("pedgaus", "gaus", ihadc1, ihadc2, TF1::EAddToList::kNo);
  fitter.SetParameters(0.1*ph->Integral(), ph->GetMean(), 5.0);
  if ( m_sigmaMax > m_sigmaMin ) {
    fitter.SetParLimits(2, m_sigmaMin, m_sigmaMax);
  } else if ( m_sigmaMax == m_sigmaMin ) {
    fitter.FixParameter(2, m_sigmaMin);
  } else {
    fitter.SetParLimits(2, 0.1, 100.0);
  }
  TF1* pfinit = dynamic_cast<TF1*>(fitter.Clone("pedgaus0"));
  pfinit->SetLineColor(3);
  pfinit->SetLineStyle(2);
  string fopt = "0";
  fopt = "WWBQ";
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
  //TFitResultPtr pres = ph->Fit(&fitter, fopt.c_str());
  //m_fitStatus = pres->Status();
  m_fitStatus = ph->Fit(&fitter, fopt.c_str());
  if ( pehSave != nullptr ) SetErrorHandler(pehSave);
  gErrorIgnoreLevel = levelSave;
  ph->GetListOfFunctions()->AddLast(pfinit, "0");
  ph->GetListOfFunctions()->Last()->SetBit(TF1::kNotDraw, true);
  m_fitMean = fitter.GetParameter(1);
  m_fitSigma = fitter.GetParameter(2);
  m_fitExcess = (maxCount - fitter.Eval(maxAdc()))/countSum;
  return 0;
}

//**********************************************************************

DataMap StickyCodeMetrics::getMetrics(string prefix) const {
  DataMap res;
  res.setInt(  prefix + "MaxAdc", maxAdc());
  res.setInt(  prefix + "MaxAdc2", maxAdc2());
  res.setFloat(prefix + "MeanAdc", meanAdc());
  res.setFloat(prefix + "MeanAdc2", meanAdc2());
  res.setFloat(prefix + "MaxFraction", maxFraction());
  res.setFloat(prefix + "ZeroFraction", zeroFraction());
  res.setFloat(prefix + "OneFraction", oneFraction());
  res.setFloat(prefix + "HighFraction", highFraction());
  res.setFloat(prefix + "ClassicFraction", classicFraction());
  res.setInt(prefix + "FitStatus", fitStatus());
  res.setFloat(prefix + "FitMean", fitMean());
  res.setFloat(prefix + "FitSigma", fitSigma());
  res.setFloat(prefix + "FitExcess", fitExcess());
  return res;
}

//**********************************************************************

void StickyCodeMetrics::print(string prefix) const {
  ostringstream sout;
  sout.precision(2);
  sout << prefix << "             # samples: " << nsample();
  sout << "\n";
  sout << prefix << "  Most common ADC code: " << maxAdc();
  sout << "\n";
  sout << prefix << " Next most common code: " << maxAdc2();
  sout << "\n";
  sout << prefix << "              Mean ADC: " << std::fixed << meanAdc();
  sout << "\n";
  sout << prefix << "      Mean ADC w/o max: " << std::fixed << meanAdc2();
  sout << "\n";
  sout.precision(3);
  sout << prefix << "       Frac in max bin: " << maxFraction();
  sout << "\n";
  sout << prefix << "            Frac LSB=0: " << zeroFraction();
  sout << "\n";
  sout << prefix << "            Frac LSB=1: " << oneFraction();
  sout << "\n";
  sout << prefix << "           Frac LSB=64: " << highFraction();
  sout << "\n";
  sout << prefix << "         Frac LSB=0,64: " << classicFraction();
  sout << "\n";
  sout << prefix << "            Fit status: " << fitStatus();
  sout << "\n";
  sout << prefix << "              Fit mean: " << fitMean();
  sout << "\n";
  sout << prefix << "             Fit sigma: " << fitSigma();
  sout << "\n";
  sout << prefix << "            Fit excess: " << fitExcess();
  cout << sout.str() << endl;
}

//**********************************************************************
