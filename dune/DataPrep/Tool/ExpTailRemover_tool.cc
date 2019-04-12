// ExpTailRemover_tool.cc

#include "ExpTailRemover.h"
#include <iostream>
#include <iomanip>
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "dune/DuneCommon/SampleTailer.h"
#include "dune/ArtSupport/DuneToolManager.h"
#include "TMath.h"
#include "dune-raw-data/Services/ChannelMap/PdspChannelMapService.h"

using std::string;
using std::cout;
using std::endl;
using std::setw;
using std::copy;

using Index = unsigned int;
using DoubleVector = std::vector<double>;

//**********************************************************************
// Class methods.
//**********************************************************************

ExpTailRemover::ExpTailRemover(fhicl::ParameterSet const& ps)
: m_LogLevel(ps.get<int>("LogLevel")),
  m_SignalFlag(ps.get<Index>("SignalFlag")),
  m_SignalIterationLimit(ps.get<Index>("SignalIterationLimit")),
  m_SignalTool(ps.get<string>("SignalTool")),
  m_DecayTime(ps.get<double>("DecayTime")) ,
  m_CorrectFlag(ps.get<std::vector<bool> >("CorrectFlag")),
  m_pSignalTool(nullptr) {
  const string myname = "ExpTailRemover::ctor: ";
  if ( m_SignalFlag > 3 ) {
    Index signalFlag = 2;
    cout << myname << "WARNING: Invalid SignalFlag value " << m_SignalFlag
         << " reset to " << signalFlag << "." << endl;
    m_SignalFlag = signalFlag;
  }
  if ( m_SignalFlag >= 2  ) {
    DuneToolManager* ptm = DuneToolManager::instance();
    if ( ptm == nullptr ) {
      cout << myname << "ERROR: Unable to retrieve tool manager." << endl;
    } else {
      m_pSignalTool = ptm->getShared<AdcChannelTool>(m_SignalTool);
      if ( m_pSignalTool == nullptr ) {
        cout << myname << "ERROR: Signal finding tool not found: " << m_SignalTool << endl;
      }
    }
  }
  if ( m_LogLevel >= 1 ) {
    cout << myname << "Parameters:" << endl;
    cout << myname << "              LogLevel: " << m_LogLevel << endl;
    cout << myname << "            SignalFlag: " << m_SignalFlag << endl;
    cout << myname << "  SignalIterationLimit: " << m_SignalIterationLimit << endl;
    cout << myname << "            SignalTool: " << m_SignalTool << endl;
    cout << myname << "             DecayTime: " << m_DecayTime << endl;
    cout << myname << "           CorrectFlag: [";
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
  DataMap ret;

  // Save the data before tail removal.
  AdcSignalVector samples = acd.samples;
  Index nsam = samples.size();

  // Check input data size.
  if ( nsam < 10 ) {
    cout << myname << "WARNING: Data for channel " << acd.channel << " has "
         << ( nsam==0 ? "no" : "too few" ) << " ticks." << endl;
    return ret.setStatus(1);;
  }

  // Check plane.
  if ( m_CorrectFlag.size() ) {
    art::ServiceHandle<dune::PdspChannelMapService> channelMap;
    size_t offlineChannel = acd.channel;
    size_t plane = channelMap->PlaneFromOfflineChannel(offlineChannel);
    if ( plane >= m_CorrectFlag.size() ) {
      cout << myname << "WARNING: Unexpected plane index: " << plane << "." << endl;
      return ret.setStatus(2);  
    }
    if ( ! m_CorrectFlag[plane] ) return ret.setStatus(3);
  }

  if ( m_LogLevel >= 2 ) cout << myname << "Correcting run " << acd.run << " event " << acd.event
                              << " channel " << acd.channel << endl;

  // Build the initial signal selection.
  bool checkSignal = true;   // Whether to use only non-signal in fit.
  bool findSignal = false;   // Whether to find signals each iteration.
  if ( m_SignalFlag == 0 ) {
    checkSignal = false;
  } else if ( m_SignalFlag == 1  ) {
    if ( acd.signal.size() < nsam ) {
      cout << myname << "WARNING: Data is missing signal flags--padding from " << acd.signal.size()
           << " to " << nsam << " samples." << endl;
    }
  } else if ( m_SignalFlag >= 2  ) {
    if ( m_pSignalTool == nullptr ) {
      cout << myname << "WARNING: Signal-finding tool is missing. Using all signals." << endl;
      checkSignal = false;
    } else {
      findSignal = true;
    }
  }

  Index niter = 0;  // # fit iterations
  float ped = 0.0;  // Fitted pedestal
  float tau = 0.0;  // Fitted tail0
  Index nsamKeep = 0;
  Index maxiter = findSignal ? m_SignalIterationLimit : 1;
  double noise = 0.0;
  while ( niter < maxiter ) {
    // Do signal finding.
    if ( findSignal ) {
      AdcFilterVector signalLast = acd.signal;
      DataMap fret = m_pSignalTool->update(acd);
      if ( fret ) {
        cout << myname << "WARNING: Signal-finding failed for event " << acd.event
             << " channel " << acd.channel << endl;
        break;
      }
      if ( acd.signal == signalLast ) {
        if ( m_LogLevel >=3 ) cout << myname << "Signal is unchanged. Exiting loop." << endl;
        break;
      }
    }
    //
    // Evaluate the signal coefficients.
    SampleTailer sta(m_DecayTime);
    sta.setTail0(1.0);
    sta.setDataZero(nsam);
    DoubleVector ctau(nsam);
    copy(sta.signal().begin(), sta.signal().end(), ctau.begin());
    sta.setTail0(0.0);
    sta.setPedestal(1.0);
    sta.setDataZero(nsam);
    DoubleVector cped(nsam);
    copy(sta.signal().begin(), sta.signal().end(), cped.begin());
    sta.setPedestal(0.0);
    sta.setData(samples);
    DoubleVector cdat(nsam);
    copy(sta.signal().begin(), sta.signal().end(), cdat.begin());
    // Evaluate the K-params.
    double kdd = 0.0;
    double kdt = 0.0;
    double kdp = 0.0;
    double ktt = 0.0;
    double ktp = 0.0;
    double kpp = 0.0;
    nsamKeep = 0;
    for ( Index  isam=0; isam<nsam; ++isam ) {
      if ( checkSignal &&  isam < acd.signal.size() && acd.signal[isam] ) continue;
      double cd = cdat[isam];
      double ct = ctau[isam];
      double cp = cped[isam];
      kdd += cd*cd;
      kdt += cd*ct;
      kdp += cd*cp;
      ktt += ct*ct;
      ktp += ct*cp;
      kpp += cp*cp;
      ++nsamKeep;
    }
    // Invert matrix and solve for (ped, tau).
    double den = ktt*kpp - ktp*ktp;
    if ( den == 0.0 ) {
      cout << myname << "WARNING: Unable to invert K-matrix with "
           << nsamKeep << " samples--stopping iteration." << endl;
      break;
    }
    double deninv = 1.0/den;
    tau = deninv*(kdp*ktp-kdt*kpp);
    ped = deninv*(kdt*ktp-kdp*ktt);
    double chsq = kdd + ktt*tau*tau + kpp*ped*ped + 2.0*(kdt*tau + kdp*ped + ktp*tau*ped);
    noise = 0.0;
    if ( nsamKeep > 2 && chsq > 0.0 ) {
      noise =  sqrt(chsq/(nsamKeep-2));
    }
    if ( m_LogLevel >= 3 ) cout << myname << "Iteration " << niter << ": ped, tau, noise: "
                                << ped << ", " << tau << ", " << noise
                                << " (" << nsamKeep << " samples)." << endl;
    // Use the parameters to remove tail from the original data and store
    // these tail-subtracted samples in the channel data.
    sta.setTail0(tau);
    sta.setPedestal(ped);
    sta.setData(samples);
    acd.samples = sta.signal();
    // Update iteration count.
    ++niter;
  }

  // Log result of interation.
  if ( m_LogLevel >= 2 ) {
    cout << myname << "Iteration count: " << niter << endl;
    cout << myname << "Final ped, tau: " << ped << ", " << tau << endl;
  }

  // Use the fitted
  acd.metadata["uscPedestal"] = ped;
  acd.metadata["uscTail"] = tau;
  acd.metadata["uscNoise"] = noise;
  ret.setFloat("uscPedestal", ped);
  ret.setFloat("uscTail", tau);
  ret.setFloat("uscNoise", noise);
  ret.setInt("uscNsamFit", nsamKeep);
  ret.setInt("uscNiteration", niter);

  return ret;
}

//**********************************************************************
