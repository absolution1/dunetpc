// ExpTailRemover_tool.cc

#include "ExpTailRemover.h"
#include <iostream>
#include <iomanip>
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "dune/Utilities/SignalShapingServiceDUNE.h"
#include "TMath.h"
#include "dune-raw-data/Services/ChannelMap/PdspChannelMapService.h"

using std::string;
using std::cout;
using std::endl;
using std::setw;

using Index = unsigned int;

//**********************************************************************
// Class methods.
//**********************************************************************

ExpTailRemover::ExpTailRemover(fhicl::ParameterSet const& ps)
: m_LogLevel(ps.get<int>("LogLevel")),
  m_SignalUnit(ps.get<string>("SignalUnit")),
  m_DecayTime(ps.get<double>("DecayTime")) ,
  m_SignalThreshold(ps.get<double>("SignalThreshold")),
  m_CorrectFlag(ps.get<std::vector<bool> >("CorrectFlag")) {
  const string myname = "ExpTailRemover::ctor: ";
  if ( m_LogLevel >= 1 ) {
    cout << myname << "Parameters:" << endl;
    cout << myname << "         LogLevel: " << m_LogLevel << endl;
    cout << myname << "       SignalUnit: " << m_SignalUnit << endl;
    cout << myname << "        DecayTime: " << m_DecayTime << endl;
    cout << myname << "  SignalThreshold: " << m_SignalThreshold << endl;
    cout << myname << "      CorrectFlag: [";
    for ( Index iori=0; iori<m_CorrectFlag.size(); ++iori ) {
      if ( iori ) cout << ", ";
      cout  << (m_CorrectFlag[iori] ? "true" : "false" );
    }
    cout << endl;
  }
}

//**********************************************************************

DataMap ExpTailRemover::update(AdcChannelData& acd) const {
  const string myname = "ExpTailRemover::view: ";
  if ( m_LogLevel >= 2 ) cout << "Processing run " << acd.run << " event " << acd.event
                              << " channel " << acd.channel << endl;
  DataMap ret;

  AdcSignalVector& samples = acd.samples;

  // Check input data size.
  size_t nsam = samples.size();
  if ( nsam < 10 ) {
    cout << myname << "WARNING: Data for channel " << acd.channel << " has "
         << ( nsam==0 ? "no" : "too few" ) << " ticks." << endl;
    return ret.setStatus(1);;
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

  // Check plane.
  if ( m_CorrectFlag.size() ) {
    art::ServiceHandle<dune::PdspChannelMapService> channelMap;
    size_t offlineChannel = acd.channel;
    size_t plane = channelMap->PlaneFromOfflineChannel(offlineChannel);
    if ( plane >= m_CorrectFlag.size() ) {
      cout << myname << "WARNING: Unexpected plane index: " << plane << "." << endl;
      return ret.setStatus(2);;  
    }
    if ( ! m_CorrectFlag[plane] ) return ret.setStatus(3);
  }

  if ( m_LogLevel >= 2 ) cout << myname << "Correcting channel " << acd.channel << endl;

  double pedfit = 0;
  double csifit = 0;
  Index nsamFit = 0;
  std::vector<double> qdats(nsam);
  for ( size_t isam=0; isam<nsam; ++isam ) qdats[isam] = samples[isam];
  estimatepars(qdats, pedfit, csifit, nsamFit);
  std::vector<double> qsigs(nsam);
  getSignal(qdats, pedfit, csifit, qsigs);
  for ( size_t isam=0; isam<nsam; ++isam ) samples[isam] = qsigs[isam];
  if ( m_LogLevel >= 3 ) {
    cout << "  Fit residual pedestal: " << pedfit << endl;
    cout << "       Fit initial tail: " << csifit << endl;
  }
  if ( m_LogLevel >= 4 ) {
    cout << myname << "fit data:" << endl;
    cout << myname << "               dat         sig" << endl;
    for ( unsigned int isam=0; isam<nsam; ++isam ) {
      cout << myname << setw(6) << isam << setw(12) << qdats[isam]  << setw(12) << qsigs[isam] << endl;
    }
  }
  acd.metadata["uscPedestal"] = pedfit;
  acd.metadata["uscTail"] = csifit;
  ret.setFloat("uscPedestal", pedfit);
  ret.setFloat("uscTail", csifit);
  ret.setInt("uscNsamFit", nsamFit);

  return ret;
}

//**********************************************************************

// Fit here for the pedestal ped and initial charge sum csi. If these values
// change by Dped and and Dcsi, the charge in bin i will shift by
//
//  DQ_i = Dcsi alpha beta^i + Dped
//
// where alpha = 1/tdecay and  beta = exp(-alpha);
//
// So if we fit Q_i vs. (-alpha beta^i) for signal-free bins, the slope and
// intercept directly give the shifts in csi and ped.

void ExpTailRemover::
estimatepars(const Vector& qdats, double& ped, double& csi, Index& nsamFit) const {
  const string myname = "ExpTailRemover::estimatepars: ";
  ped = 0;
  csi = 0;
  nsamFit = 0;
  unsigned int nsam = qdats.size();
  double alpha = 1.0/m_DecayTime;
  double beta = exp(-alpha);
  alpha *= sqrt(beta);
  std::vector<double> xfit(nsam);
  std::vector<double> corr1(nsam);
  std::vector<double> corr1e(nsam);
  getSignal(qdats, 0, 0, corr1);  // initial try -- just take out undershoot, no pedestal or initial charge sum

  // Set x-values for fit: alpha beta^i
  double xval = alpha;
  xval = 1.0;
  for ( unsigned int isam=0; isam<nsam; ++isam ) {
    xfit[isam] = xval;
    xval *= beta;
  }
  // use the same convention for error bars the event display added -- seems to result in a better fit
  for ( Index isam=0; isam<nsam; ++isam ) {
    corr1e[isam] = TMath::Max(1.0, TMath::Abs(corr1[isam]));
  }
  double slope=0;
  double intercept=0;
  wlinfit(xfit, corr1, corr1e, slope, intercept);
  if ( m_LogLevel >= 5 ) {
    cout << myname << "Initial fit data:" << endl;
    cout << myname << "             xfit,           dat" << endl;
    for ( unsigned int isam=0; isam<nsam; ++isam ) {
      cout << myname << setw(6) << isam << setw(12) << xfit[isam]  << ", " << setw(12) << corr1[isam]
           << " +/- " << setw(12) << corr1e[isam] << endl;
    }
    cout << "Intercept/pedestal: " << intercept << endl;
    cout << "         Slope/csi: " << slope << endl;
  }

  // de-weight bins that differ from linear prediction
  for ( unsigned int isam=0; isam<nsam; ++isam ) {
    double xcent = xfit[isam];
    double yval = corr1[isam];
    float pred = intercept + slope*xcent;
    if ( TMath::Abs(yval - pred) > m_SignalThreshold ) {
      corr1e[isam] = 1000;
    } else {
      corr1e[isam] = 1;
      ++nsamFit;
    }
  }
  // fit again with outliers de-weighted
  wlinfit( xfit, corr1, corr1e, slope, intercept);
  if ( m_LogLevel >= 5 ) {
    cout << myname << "Second fit data:" << endl;
    cout << myname << "             xfit,           dat" << endl;
    for ( unsigned int isam=0; isam<nsam; ++isam ) {
      cout << myname << setw(6) << isam << setw(12) << xfit[isam] << ", " << setw(12) << corr1[isam]
           << " +/- " << setw(12) << corr1e[isam] << endl;
    }
    cout << "Intercept/pedestal: " << intercept << endl;
    cout << "         Slope/csi: " << slope << endl;
  }
  ped = intercept;
  csi =  slope;
}

//**********************************************************************

void ExpTailRemover::
getSignal(const Vector& qdats, double ped, double csi, Vector& qsigs) const {
  double alpha = m_DecayTime > 0.0 ? 1.0/m_DecayTime : 0.0;
  double beta = exp(-alpha);
  alpha *= sqrt(beta);  // Improves tail cancellation but is a tiny change
  Index nsam = qdats.size();
  double qtai = csi;  // This is the charge in the current tail sample
  alpha *= sqrt(beta);  // Improves tail cancellation but is a tiny change
  for ( Index isam=0; isam<nsam; ++isam ) {
    float qcor = qdats[isam] - ped;
    float qsig = qcor - qtai;
    qsigs[isam] = qsig;
    qtai -= alpha*qsig;    // Add signal contribution to the tail.
    qtai *= beta;          // Decay the tail by one tick
  }  
}

//**********************************************************************

void ExpTailRemover::
wlinfit(const Vector& x, const Vector& y, const Vector& e,
        double &slope, double &intercept) const {
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
