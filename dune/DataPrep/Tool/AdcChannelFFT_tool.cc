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

//**********************************************************************
// Class methods.
//**********************************************************************

AdcChannelFFT::AdcChannelFFT(fhicl::ParameterSet const& ps)
: m_LogLevel(ps.get<int>("LogLevel")), 
  m_FirstTick(ps.get<Index>("FirstTick")),
  m_NTick(ps.get<Index>("NTick")),
  m_NormOpt(ps.get<Index>("NormOpt")),
  m_Action(ps.get<Index>("Action")),
  m_ReturnOpt(ps.get<Index>("ReturnOpt")) {
  const string myname = "AdcChannelFFT::ctor: ";
  if ( m_LogLevel ) {
    cout << myname << "Configuration: " << endl;
    cout << myname << "            LogLevel: " << m_LogLevel << endl;
    cout << myname << "           FirstTick: " << m_FirstTick << endl;
    cout << myname << "               NTick: " << m_NTick << endl;
    cout << myname << "             NormOpt: " << m_NormOpt << endl;
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
internalView(const AdcChannelData& acd, FloatVector& sams, FloatVector& xmgs, FloatVector& xphs, DataMap& ret) const {
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
  DataMap::FloatVector xres;
  DataMap::FloatVector xims;
  Index isam0 = 0;
  Index nsam = 0;
  if ( doForward ) {
    isam0 = m_FirstTick;
    if ( isam0 >= acd.samples.size() ) {
      cout << myname << "No data in range." << endl;
      ret.setStatus(11);
      return;
    }
    nsam = acd.samples.size() - isam0;
    if ( m_NTick > 0 && m_NTick < nsam ) nsam = m_NTick;
    int rstat = DuneFFT::fftForward(m_NormOpt, nsam, &acd.samples[isam0], xres, xims, xmgs, xphs, m_LogLevel);
    if ( rstat ) {
      ret.setStatus(10+rstat);
      return;
    }
  } else if ( doInverse ) {
    int rstat = DuneFFT::fftInverse(m_NormOpt, acd.dftmags, acd.dftphases, xres, xims, sams, m_LogLevel);
    xmgs = acd.dftmags;
    xphs = acd.dftphases;
    if ( rstat ) {
      ret.setStatus(20+rstat);
      return;
    }
  }
  if ( m_ReturnOpt >= 1 ) {
    ret.setInt("fftTick0", isam0);
    ret.setInt("fftNTick", nsam);
  }
  Index dftRet = m_ReturnOpt % 10;
  if ( dftRet >= 1 ) {
    ret.setInt("fftNMag",   xmgs.size());
    ret.setInt("fftNPhase", xphs.size());
  }
  if ( dftRet >= 2 ) {
    ret.setFloatVector("fftMags",   xmgs);
    ret.setFloatVector("fftPhases", xphs);
  }
  if ( dftRet >= 3 ) {
    ret.setFloatVector("fftReals",  xres);
    ret.setFloatVector("fftImags",  xims);
  }
  if ( m_ReturnOpt >= 10 ) {
    if ( sams.size() ) ret.setFloatVector("fftSamples", sams);
    else ret.setFloatVector("fftSamples", acd.samples);
  }
  return;
}

//**********************************************************************
