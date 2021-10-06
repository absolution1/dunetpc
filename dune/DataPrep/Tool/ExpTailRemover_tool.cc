// ExpTailRemover_tool.cc

#include "ExpTailRemover.h"
#include <iostream>
#include <iomanip>
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "dune/DuneCommon/Utility/SampleTailer.h"
#include "dune/ArtSupport/DuneToolManager.h"
#include "dune/DuneInterface/Tool/IndexRangeTool.h"
#include "dune-raw-data/Services/ChannelMap/PdspChannelMapService.h"
//#include "TMath.h"

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
  m_IncludeChannelRanges(ps.get<NameVector>("IncludeChannelRanges")),
  m_ExcludeChannelRanges(ps.get<NameVector>("ExcludeChannelRanges")),
  m_useChannelRanges(false),
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
  Index nchaCheck = 0;
  if ( m_IncludeChannelRanges.size() ) {
    DuneToolManager* ptm = DuneToolManager::instance();
    const IndexRangeTool* pcrt = ptm->getShared<IndexRangeTool>("channelRanges");
    if ( pcrt == nullptr ) {
      cout << myname << "ERROR: IndexRangeTool not found: channelRanges" << endl;
    } else {
      for ( Name crn : m_IncludeChannelRanges ) {
        IndexRange ran = pcrt->get(crn);
        if ( ran.isValid() ) {
          if ( ran.end > m_checkChannels.size() ) m_checkChannels.resize(ran.end, false);
          for ( Index icha=ran.begin; icha<ran.end; ++icha ) m_checkChannels[icha] = true;
        } else {
          cout << myname << "WARNING: Ignoring invalid include channel range " << crn << endl;
        }
      }
      if ( m_ExcludeChannelRanges.size() ) {
        for ( Name crn : m_ExcludeChannelRanges ) {
          if ( crn == "all" ) {
            m_checkChannels.clear();
            break;
          }
          IndexRange ran = pcrt->get(crn);
          if ( ran.isValid() ) {
            Index end = ran.end < m_checkChannels.size() ? ran.end : m_checkChannels.size();
            for ( Index icha=ran.begin; icha<end; ++icha ) m_checkChannels[icha] = false;
          } else {
            cout << myname << "WARNING: Ignoring invalid exclude channel range " << crn << endl;
          }
        }
      }
      m_useChannelRanges = true;
      for ( bool keep : m_checkChannels ) if ( keep ) ++nchaCheck;
    }
  }
  if ( m_LogLevel >= 1 ) {
    cout << myname << "Parameters:" << endl;
    cout << myname << "              LogLevel: " << m_LogLevel << endl;
    cout << myname << "            SignalFlag: " << m_SignalFlag << endl;
    cout << myname << "  SignalIterationLimit: " << m_SignalIterationLimit << endl;
    cout << myname << "            SignalTool: " << m_SignalTool << endl;
    cout << myname << "             DecayTime: " << m_DecayTime << endl;
    cout << myname << "  IncludeChannelRanges: [";
    bool first = true;
    for ( Name crn : m_IncludeChannelRanges ) {
      if ( first ) first = false;
      else cout << ", ";
      cout  << crn;
    }
    cout << "]" << endl;
    cout << myname << "  ExcludeChannelRanges: [";
    first = true;
    for ( Name crn : m_ExcludeChannelRanges ) {
      if ( first ) first = false;
      else cout << ", ";
      cout  << crn;
    }
    cout << "]" << endl;
    if ( m_useChannelRanges ) {
      cout << myname << "Channel checking enabled for " << nchaCheck << " channel"
           << ( nchaCheck == 1 ? "" : "s") << "." << endl;
    } else {
      cout << myname << "Channel checking disabled." << endl;
    }
  }
}

//**********************************************************************

DataMap ExpTailRemover::update(AdcChannelData& acd) const {
  const string myname = "ExpTailRemover::update: ";
  DataMap ret;

  // Save the data before tail removal.
  AdcSignalVector samples = acd.samples;
  Index nsam = samples.size();

  // Check the channel.
  if ( m_useChannelRanges ) {
    if ( acd.channel() >= m_checkChannels.size() || ! m_checkChannels[acd.channel()] ) {
      if ( m_LogLevel >= 2 ) cout << myname << "Skipping channel " << acd.channel() << endl;
      return ret;
    }
  }

  // Check input data size.
  if ( nsam < 10 ) {
    cout << myname << "WARNING: Data for channel " << acd.channel() << " has "
         << ( nsam==0 ? "no" : "too few" ) << " ticks." << endl;
    return ret.setStatus(1);
  }

  if ( m_LogLevel >= 2 ) cout << myname << "Correcting run " << acd.run() << " event " << acd.event()
                              << " channel " << acd.channel() << endl;

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
        cout << myname << "WARNING: Signal-finding failed for event " << acd.event()
             << " channel " << acd.channel() << endl;
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
      if ( acd.channelStatus() == 0 || m_LogLevel >= 2 ) {
        cout << myname << "WARNING: Unable to invert K-matrix with "
             << nsamKeep << " of " << nsam << " samples--stopping iteration for channel "
             << acd.channel() << " with status " << acd.channelStatus() << "." << endl;
      }
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

DEFINE_ART_CLASS_TOOL(ExpTailRemover)
