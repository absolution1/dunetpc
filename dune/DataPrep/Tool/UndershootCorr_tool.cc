// UndershootCorr_tool.cc

#include "UndershootCorr.h"
#include <iostream>
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "lardata/Utilities/LArFFT.h"
#include "dune/Utilities/SignalShapingServiceDUNE.h"
#include "TMath.h"
#include "dune-raw-data/Services/ChannelMap/PdspChannelMapService.h"


using std::string;
using std::cout;
using std::endl;

//**********************************************************************
// Class methods.
//**********************************************************************

UndershootCorr::UndershootCorr(fhicl::ParameterSet const& ps)
: m_LogLevel(ps.get<int>("LogLevel"))
, m_CorrectFlag(ps.get<std::vector<bool> >("CorrectFlag"))
, m_TDecayConst(ps.get<std::vector<double> >("TDecayConst")) 
, m_FSubConst(ps.get<std::vector<double> >("FSubConst"))
, m_RestoreBaseline(ps.get<std::vector<bool> >("RestoreBaseline")) 
, m_SignalThreshold(ps.get<std::vector<double> >("SignalThreshold")) 
, m_SignalUnit(ps.get<std::string>("SignalUnit")) 
, m_LCA(ps.get<std::vector<double> >("LCA"))
, m_LCB(ps.get<std::vector<double> >("LCB"))
, m_LCC(ps.get<std::vector<double> >("LCC"))
, m_LCD(ps.get<std::vector<double> >("LCD")) {
  const string myname = "UndershootCorr::ctor: ";
  if ( m_LogLevel >= 1 ) {
    cout << myname << "Parameters:" << endl;
    cout << myname << "  LogLevel: " << m_LogLevel << endl;
    for ( unsigned int i=0; i<3; ++i ) {
      cout << myname << "  For view " << i << endl;
      cout << myname << "         CorrectFlag: " << (m_CorrectFlag[i] ? "true" : "false" ) << endl;
      cout << myname << "         TDecayConst: " << m_TDecayConst[i] << endl;
      cout << myname << "           FSubConst: " << m_FSubConst[i] << endl;
      cout << myname << "     RestoreBaseline: " << (m_RestoreBaseline[i] ? "true" : "false" ) << endl;
      cout << myname << "     SignalThreshold: " << m_SignalThreshold[i] << endl;
      cout << myname << "          SignalUnit: " << m_SignalUnit << endl;
      cout << myname << "                 LCA: " << m_LCA[i] << endl;
      cout << myname << "                 LCB: " << m_LCB[i] << endl;
      cout << myname << "                 LCC: " << m_LCC[i] << endl;
      cout << myname << "                 LCD: " << m_LCD[i] << endl;
    }
  }
}

//**********************************************************************

DataMap UndershootCorr::update(AdcChannelData& acd) const {
  const string myname = "UndershootCorr::view: ";
  if ( m_LogLevel >= 2 ) cout << "Processing run " << acd.run() << " event " << acd.event()
                              << " channel " << acd.channel() << endl;
  DataMap ret;

  AdcSignalVector& samples = acd.samples;

  // Check input data size.
  size_t nticks = samples.size();
  if ( nticks < 10 ) {
    cout << myname << "Data for channel " << acd.channel() << " has "
         << ( nticks==0 ? "no" : "too few" ) << " ticks." << endl;
    return ret;
  }

  // Check input data unit.
  if ( m_SignalUnit.size() ) {
    if ( acd.sampleUnit.size() == 0 ) {
      cout << myname << "WARNING: Input data does not have a sample unit." << endl;
    } else if ( acd.sampleUnit != m_SignalUnit ) {
      cout << myname << "WARNING: Unexpected input data unit: " << acd.sampleUnit
           << " != " << m_SignalUnit << endl;
    }
  }

  art::ServiceHandle<dune::PdspChannelMapService> channelMap;
  size_t offlineChannel = acd.channel();
  size_t plane = channelMap->PlaneFromOfflineChannel(offlineChannel);
  if (plane >= m_CorrectFlag.size()) return ret;  
  double pedfit = 0;
  double csifit = 0;
  if ( m_CorrectFlag[plane] ) {

    double median = TMath::Median(nticks,samples.data());

    std::vector<double> x(nticks);
    std::vector<double> yorig(nticks);
    std::vector<double> ycorr(nticks);

    bool allempty = true;
    for (size_t i=0;i<nticks;++i) {
      double bc = samples[i]-median;
      if (bc != 0) allempty=false;
      yorig[i] = bc;
      x[i] = i;
    }
    if (allempty) return ret;
    estimatepars(x,yorig,pedfit,csifit,plane);
    crc(yorig,ycorr,pedfit,csifit,plane);
    double offset = m_RestoreBaseline[plane] ? median : 0.0;
    for ( size_t i=0; i<nticks; ++i ) samples[i] = ycorr[i] + offset;
    if ( m_LogLevel >= 3 ) {
      cout << "  Median: " << median
           << " (" << (m_RestoreBaseline[plane] ? "" : "not ") << "restored)" << endl;
      cout << "  Fit pedestal: " << pedfit << endl;
      cout << "    Fit offset: " << csifit << endl;
    }
  }
  acd.metadata["uscPedestal"] = pedfit;
  acd.metadata["uscPedestalOffset"] = csifit;

  return ret;
}

//**********************************************************************

void UndershootCorr::estimatepars(std::vector<double> &x, std::vector<double> &y, double &pedi,
                                  double &csi, size_t plane) const {
  pedi = 0;
  csi = 0;
  int nbinsx = x.size();
  std::vector<double> corr1(nbinsx);
  std::vector<double> corr1e(nbinsx);
  crc(y,corr1,0,0,plane);  // initial try -- just take out undershoot, no pedestal or initial charge sum

  // use the same convention for error bars the event display added -- seems to result in a better fit
  for (int i=0;i<nbinsx;++i) {
    corr1e[i] = TMath::Max(1.0,TMath::Abs(corr1[i]));
  }
  double slope=0;
  double intercept=0;
  wlinfit(x,corr1,corr1e,slope,intercept);

  // de-weight bins that differ from linear prediction
  for (int i=0;i<nbinsx;++i) {
    double xcent = x[i];
    double yval = corr1[i];
    float pred = intercept + slope*xcent;
    if ( TMath::Abs(yval - pred) > m_SignalThreshold[plane] ) {
      corr1e[i] = 1000;
    } else {
      corr1e[i] = 1;
    }
  }
  // fit again with outliers de-weighted
  wlinfit(x,corr1,corr1e,slope,intercept);

  // parameters tuned for the collection plane Nov. 2018 for ProtoDUNE-SP, now read from fcl.

  //pedi = 1980*slope - 0.08277*intercept;
  //csi = -2043.65*slope + 1.03*intercept;

  pedi = m_LCA[plane]*slope + m_LCB[plane]*intercept;
  csi =  m_LCC[plane]*slope + m_LCD[plane]*intercept;
}

//**********************************************************************

void UndershootCorr::crc(std::vector<double> &orig, std::vector<double> &corr,
                         double pedi, double csi, size_t plane) const {
  int nbinsx = orig.size();
  double csum=csi;
  double trueped=pedi;
  for (int i=0;i<nbinsx;++i) {
    float bc = orig[i] - trueped;
    float cv = bc - csum;
    corr[i] = cv;

    // parameters tuned for the collection plane Nov. 2018 for ProtoDUNE-SP, now read from fcl.
    //csum -= cv*0.0005;
    //csum *= 0.99955;

    csum -= cv*m_FSubConst[plane];
    csum *= m_TDecayConst[plane];
  }  
}

//**********************************************************************

void UndershootCorr::wlinfit(std::vector<double> x, std::vector<double> &y,
                             std::vector<double> &e, double &slope, double &intercept) const {
  slope = 0;
  intercept = 0;

  double sumx=0;
  double sumy=0;
  double sumxy=0;
  double sum=0;
  double sumxx=0;
  size_t npts = x.size();

  for ( size_t i=0; i<npts; ++i ) {
    double ooe2 = 1.0/TMath::Sq(e[i]);
    sumx += x[i]*ooe2;
    sumy += y[i]*ooe2;
    sumxy += x[i]*y[i]*ooe2;
    sumxx += TMath::Sq(x[i])*ooe2;
    sum += ooe2;
  }

  double denom = TMath::Sq(sumx) - sumxx*sum;
  if ( denom != 0 ) {
    slope = (sumx*sumy - sumxy*sum)/denom;
    intercept = (sumx*sumxy - sumy*sumxx)/denom;
  }
}

//**********************************************************************

DEFINE_ART_CLASS_TOOL(UndershootCorr)
