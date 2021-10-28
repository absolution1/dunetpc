// AdcNoiseSignalFinder_tool.cc

#include "AdcNoiseSignalFinder.h"
#include <iostream>

using std::string;
using std::cout;
using std::endl;

//**********************************************************************
// Local defintitions.
//**********************************************************************

namespace {

string boolToString(bool val) {
  return val ? "true" : "false";
}

}  // end unnamed namespace

//**********************************************************************
// Class methods.
//**********************************************************************

AdcNoiseSignalFinder::AdcNoiseSignalFinder(fhicl::ParameterSet const& ps)
: m_LogLevel(ps.get<int>("LogLevel")),
  m_SigFracMax(ps.get<float>("SigFracMax")),
  m_ThresholdMin(ps.get<float>("ThresholdMin")),
  m_ThresholdRatio(ps.get<float>("ThresholdRatio")),
  m_ThresholdRatioTol(ps.get<float>("ThresholdRatioTol")),
  m_MaxLoop(ps.get<unsigned int>("MaxLoop")),
  m_BinsBefore(ps.get<unsigned int>("BinsBefore")),
  m_BinsAfter(ps.get<unsigned int>("BinsAfter")),
  m_FlagPositive(ps.get<bool>("FlagPositive")),
  m_FlagNegative(ps.get<bool>("FlagNegative")) {
  const string myname = "AdcNoiseSignalFinder::ctor: ";
  if ( m_LogLevel >= 1 ) {
    cout << myname << "Configuration parameters:" << endl;
    cout << myname << "          LogLevel: " << m_LogLevel << endl;
    cout << myname << "        SigFracMax: " << m_SigFracMax << endl;
    cout << myname << "      ThresholdMin: " << m_ThresholdMin << endl;
    cout << myname << "    ThresholdRatio: " << m_ThresholdRatio << endl;
    cout << myname << " ThresholdRatioTol: " << m_ThresholdRatioTol << endl;
    cout << myname << "        BinsBefore: " << m_BinsBefore << endl;
    cout << myname << "           MaxLoop: " << m_MaxLoop << endl;
    cout << myname << "         BinsAfter: " << m_BinsAfter << endl;
    cout << myname << "      FlagPositive: " << boolToString(m_FlagPositive) << endl;
    cout << myname << "      FlagNegative: " << boolToString(m_FlagNegative) << endl;
  }
}

//**********************************************************************

DataMap AdcNoiseSignalFinder::update(AdcChannelData& acd) const {
  const string myname = "AdcNoiseSignalFinder::update: ";
  DataMap ret;
  AdcIndex nsam = acd.samples.size();
  if ( m_ThresholdRatio <= 0 ) {
    cout << myname << "Invalid noise ratio: " << m_ThresholdRatio << endl;
    return ret.setStatus(2);
  }
  if ( m_ThresholdRatioTol <= 0 ) {
    cout << myname << "Invalid noise ratio tolerance: " << m_ThresholdRatioTol << endl;
    return ret.setStatus(2);
  }
  if ( nsam == 0 ) {
    cout << myname << "ERROR: No samples found in channel " << acd.channel() << endl;
    acd.signal.clear();
    acd.rois.clear();
    return ret.setStatus(1);
  }
  if ( m_LogLevel >= 3 ) cout << myname << "Finding ROIs for channel " << acd.channel() << endl;
  float thr = m_ThresholdMin;
  float noise = 0.0;
  float sigfrac = 0.0;
  Index nloop = 0;
  float trtgt = m_ThresholdRatio;
  float trmin = trtgt - m_ThresholdRatioTol;
  float trmax = trtgt + m_ThresholdRatioTol;
  float thrTooLow = m_ThresholdMin;  // We know threshold is at or above this value
  float thrTooHigh = 1.e20;          // We know threshold is below this value
  while ( true ) {
    // Find the and flag the signal.
    AdcIndex nsamlo = m_BinsBefore;
    AdcIndex nsamhi = m_BinsAfter;
    AdcIndex nsamAbove = 0;
    AdcIndex isamUnknown = 0;  // First sample not known to be in or outside a ROI.
    acd.signal.clear();
    acd.signal.resize(nsam, false);
    for ( AdcIndex isam=0; isam<nsam; ++isam ) {
      bool keep = ( m_FlagPositive && acd.samples[isam] >  thr ) ||
                  ( m_FlagNegative && acd.samples[isam] < -thr );
      if ( keep ) {
        ++nsamAbove;
        AdcIndex jsam1 = isam > nsamlo ? isam - nsamlo : 0;
        if ( jsam1 < isamUnknown ) jsam1 = isamUnknown;
        AdcIndex jsam2 = isam + nsamhi + 1;
        if ( jsam2 > nsam ) jsam2 = nsam;
        if ( m_LogLevel >= 5 ) cout << myname << "Trigger: " << isam << ", range: ["
                                    << jsam1 << ", " << jsam2 << ")" << endl;
        for ( AdcIndex jsam=jsam1; jsam<jsam2; ++jsam ) acd.signal[jsam] = true;
        isamUnknown = jsam2;
      }
    }
    // Evaluate the noise and signal fraction.
    Index nsig = 0;
    Index nnsg = 0;
    float ssqsum = 0.0;
    for ( AdcIndex isam=0; isam<nsam; ++isam ) {
      if ( acd.signal[isam] ) {
        ++nsig;
      } else {
        ++nnsg;
        float val = acd.samples[isam];
        ssqsum += val*val;
      }
    }
    noise = nnsg ? sqrt(ssqsum/nnsg) : 0.0;
    sigfrac = float(nsig)/nsam;
    // Check is we need andother loop and rest threshold accordingly.
    ++nloop;
    if ( m_LogLevel >= 4 ) {
      cout << myname << "  Loop " << nloop << endl;
      cout << myname << "    Threshold: " << thr << endl;
      cout << myname << "    Noise: " << noise << endl;
      cout << myname << "    Thr/noise: " << (noise > 0.0 ? thr/noise : -999) << endl;
      cout << myname << "    # ticks above threshold: " << nsamAbove << endl;
      cout << myname << "    SigFrac: " << sigfrac << endl;
      acd.roisFromSignal();
      cout << myname << "    # ROI: " << acd.rois.size() << endl;
    }
    if ( nloop >= m_MaxLoop ) {
      cout << myname << "WARNING: Channel " << acd.channel() << " exiting after "
           << nloop << " loops." << endl;
      break;
    }
    // Use current noise estimate to evaluate the target threshold and range.
    float thrtgt = trtgt*noise;
    float thrmin = trmin*noise;
    float thrmax = trmax*noise;
    float thrOld = thr;
    // Noise level too high ==> raise threshold.
    if ( sigfrac > m_SigFracMax || thrmin > thr ) {
      thrTooLow = thr;
      // Double the threshold unless that goes past the known upper limit.
      // In that case, go halfway.
      float thrEst = thrtgt > thr ? 2.0*thrtgt : 2.0*thr;
      if ( thrEst < thrTooHigh ) thr = thrEst;
      else thr = 0.5*(thr + thrTooHigh);
    } else if ( thr <= m_ThresholdMin ) {
      break;
    } else if ( thrmax >= thr ) {
      break;
    } else {
      // Noise level too low ==> lower threshold.
      thrTooHigh = thr;
      // Use the target as the next estimate unless it below the thw known lower
      // limit. In the latter case, go halfway to that limit.
      if ( thrtgt > thrTooLow ) thr = thrtgt;
      else thr = 0.5*(thr + thrTooLow);
    }
    float dthrmax = 0.001*thrOld;
    if ( m_LogLevel >= 1 ) {
      if ( fabs(thr - thrOld) < dthrmax ) {
        Index chstat = acd.channelStatus();
        if ( m_LogLevel >= 2 || (m_LogLevel == 1 && chstat == 0) ) {
          string schstat = chstat==0 ? "Good" :
                           chstat==1 ? "Bad" :
                           chstat==2 ? "Noisy" : "UnkownStatus";
          cout << myname << "WARNING: " << schstat << " channel " << acd.channel()
             << " exiting prematurely due to threshold convergence." << endl;
        }
        break;
      }
    }
  }
  acd.roisFromSignal();
  ret.setFloat(  "nsfSigFrac", sigfrac);
  ret.setFloat(    "nsfNoise", noise);
  ret.setFloat("nsfThreshold", thr);
  ret.setInt(  "nsfLoopCount", nloop);
  ret.setInt(   "nsfRoiCount", acd.rois.size());
  if ( m_LogLevel >= 4 ) ret.print(myname);
  acd.metadata[  "nsfSigFrac"] = sigfrac;
  acd.metadata[    "nsfNoise"] = noise;
  acd.metadata["nsfThreshold"] = thr;
  acd.metadata["nsfLoopCount"] = nloop;
  acd.metadata[ "nsfRoiCount"] = acd.rois.size();

  return ret;
}
//**********************************************************************

DataMap AdcNoiseSignalFinder::view(const AdcChannelData& acd) const {
  AdcChannelData acdtmp;
  acdtmp.samples = acd.samples;
  return update(acdtmp);
}

//**********************************************************************

DEFINE_ART_CLASS_TOOL(AdcNoiseSignalFinder)
