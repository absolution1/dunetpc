// AdcChannelFFT_tool.cc

#include "AdcChannelFFT.h"
#include "dune/DuneCommon/DuneFFT.h"
#include <iostream>
#include <sstream>
#include <vector>
#include <iomanip>
#include "TComplex.h"

using std::string;
using std::cout;
using std::cin;
using std::endl;
using std::vector;
using std::setw;
using std::fixed;
using DFT = DuneFFT::DFT;

//**********************************************************************
// Class methods.
//**********************************************************************

AdcChannelFFT::AdcChannelFFT(fhicl::ParameterSet const& ps)
: m_LogLevel(ps.get<int>("LogLevel")), 
  m_FirstTick(ps.get<Index>("FirstTick")),
  m_NTick(ps.get<Index>("NTick")),
  m_Action(ps.get<Index>("Action")),
  m_ReturnOpt(ps.get<Index>("ReturnOpt")) {
  const string myname = "AdcChannelFFT::ctor: ";
  if ( m_LogLevel ) {
    cout << myname << "Configuration: " << endl;
    cout << myname << "            LogLevel: " << m_LogLevel << endl;
    cout << myname << "           FirstTick: " << m_FirstTick << endl;
    cout << myname << "               NTick: " << m_NTick << endl;
    cout << myname << "              Action: " << m_Action << endl;
    cout << myname << "           ReturnOpt: " << m_ReturnOpt << endl;
  }
}

//**********************************************************************

DataMap AdcChannelFFT::view(const AdcChannelData& acd) const {
  const string myname = "AdcChannelFFT::view: ";
  DataMap ret;
  DataMap::FloatVector sams;
  DataMap::FloatVector mags;
  DataMap::FloatVector phas;
  internalView(acd, sams, mags, phas, ret);
  return ret;
}

//**********************************************************************

DataMap AdcChannelFFT::update(AdcChannelData& acd) const {
  const string myname = "AdcChannelFFT::update: ";
  DataMap ret;
  DataMap::FloatVector sams;
  DataMap::FloatVector mags;
  DataMap::FloatVector phas;
  internalView(acd, sams, mags, phas, ret);
  if ( ret ) return ret;
  if ( m_Action == 3 || m_Action == 4 ) {
    if ( mags.size() ) {
      if ( m_LogLevel >= 2 ) cout << myname << "Saving DFT." << endl;
      acd.dftmags = mags;
      acd.dftphases = phas;
    }
  }
  if ( m_Action == 13 || m_Action == 14 ) {
    if ( sams.size() ) {
      if ( m_LogLevel >= 2 ) cout << myname << "Saving samples." << endl;
      acd.samples = sams;
    }
  }
  return ret;
}

//**********************************************************************

void AdcChannelFFT::
internalView(const AdcChannelData& acd, FloatVector& sams, FloatVector& xams, FloatVector& xphs, DataMap& ret) const {
  const string myname = "AdcChannelFFT::internalView: ";
  bool doForward = false;
  bool doInverse = false;
  if ( m_Action <= 4 ) {
    if ( m_Action > 0 ) {
      if ( m_Action == 2 || m_Action == 4 ) {
        doForward = acd.dftmags.size() == 0;
      } else {
        doForward = true;
      }
    }
  } else if ( m_Action >= 10  && m_Action <= 14 ) {
    if ( m_Action != 10 ) {
      if ( m_Action == 12 || m_Action == 14 ) {
        doInverse = acd.samples.size() == 0;
      } else {
        doInverse = true;
      }
    }
  } else {
    cout << myname << "ERROR: Invalid action: " << m_Action << endl;
    ret.setStatus(1);
    return;
  }
  Index isam0 = 0;
  Index nsam = 0;
  DFT dft(DFT::Norm(AdcChannelData::dftNormalization()));
  if ( doForward ) {
    isam0 = m_FirstTick;
    if ( isam0 >= acd.samples.size() ) {
      cout << myname << "WARNING: No data in range." << endl;
      ret.setStatus(11);
      return;
    }
    nsam = acd.samples.size() - isam0;
    if ( m_LogLevel >= 3 ) cout << myname << "Forward FFT with " << nsam << " samples." << endl;
    if ( m_NTick > 0 && m_NTick < nsam ) nsam = m_NTick;
    int rstat = DuneFFT::fftForward(nsam, &acd.samples[isam0], dft, m_LogLevel);
    if ( rstat ) {
      ret.setStatus(10+rstat);
      cout << myname << "WARNING: Forward FFT failed." << endl;
      return;
    }
  } else if ( doInverse ) {
    dft.copyIn(acd.dftmags, acd.dftphases);
    if ( ! dft.isValid() ) {
      ret.setStatus(20);
      cout << "ERROR: Unable to find DFT in AdcChannelData." << endl;
      return;
    }
    int rstat = DuneFFT::fftInverse(dft, sams, m_LogLevel);
    xams = acd.dftmags;
    xphs = acd.dftphases;
    if ( m_LogLevel >= 3 ) cout << myname << "Inverse FFT for " << dft.size() << " samples." << endl;
    if ( rstat ) {
      ret.setStatus(20+rstat);
      cout << myname << "WARNING: Inverse FFT failed." << endl;
      return;
    }
  }
  // Fetch return data stored in the DFT.
  Index dftRet = m_ReturnOpt % 10;
  if ( dftRet >= 3 ) {
    FloatVector fftres(nsam);
    FloatVector fftims(nsam);
    for ( Index ifrq=0; ifrq<nsam; ++ifrq ) {
      fftres[ifrq] = dft.real(ifrq);
      fftims[ifrq] = dft.imag(ifrq);
    }
    ret.setFloatVector("fftReals", fftres);
    ret.setFloatVector("fftImags", fftims);
  }
  // Move data out of the DFT (leaving it invalid).
  if ( doForward ) dft.moveOut(xams, xphs);
  // Fetch return data not stored in the DFT.
  if ( m_ReturnOpt >= 1 ) {
    ret.setInt("fftTick0", isam0);
    ret.setInt("fftNTick", nsam);
  }
  if ( dftRet >= 1 ) {
    ret.setInt("fftNMag",   xams.size());
    ret.setInt("fftNPhase", xphs.size());
  }
  if ( dftRet >= 2 ) {
    ret.setFloatVector("fftMags",   xams);
    ret.setFloatVector("fftPhases", xphs);
  }
  if ( m_ReturnOpt >= 10 ) {
    if ( sams.size() ) ret.setFloatVector("fftSamples", sams);
    else ret.setFloatVector("fftSamples", acd.samples);
  }
  return;
}

//**********************************************************************
