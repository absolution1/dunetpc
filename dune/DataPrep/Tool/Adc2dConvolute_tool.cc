// Adc2dConvolute_tool.cc

#include "Adc2dConvolute.h"
#include <iostream>
#include <iomanip>

using std::string;
using std::cout;
using std::endl;
using std::setw;
using std::fixed;
using std::setprecision;

using Index = unsigned int;
using FloatVector = AdcSignalVector;
using DoubleVector = std::vector<double>;
using DoubleVectorMap = std::map<Index, DoubleVector>;
using Name = std::string;

//**********************************************************************
// Class methods.
//**********************************************************************

Adc2dConvolute::Adc2dConvolute(fhicl::ParameterSet const& ps)
: m_LogLevel(ps.get<int>("LogLevel")),
  m_BinsPerTick(ps.get<Index>("BinsPerTick")),
  m_BinsPerWire(ps.get<Index>("BinsPerWire")),
  m_ResponseVectors(ps.get<ResponseVectorVector>("ResponseVectors")),
  m_ResponseCenter(ps.get<Index>("ResponseCenter")),
  m_MinChannel(ps.get<Index>("MinChannel")),
  m_MaxChannel(ps.get<Index>("MaxChannel")) {
  const string myname = "Adc2dConvolute::ctor: ";
  if ( m_LogLevel >= 1 ) {
    std::ios cout_state(nullptr);
    cout_state.copyfmt(std::cout);
    cout << myname << "Parameters:" << endl;
    cout << myname << "         LogLevel: " << m_LogLevel << endl;
    cout << myname << "      BinsPerTick: " << m_BinsPerTick << endl;
    cout << myname << "      BinsPerWire: " << m_BinsPerWire << endl;
    cout << myname << "  ResponseVectors: [";
    Index indent = 22;
    bool first = true;
    for ( const ResponseVector& vec : m_ResponseVectors ) {
      if ( first ) first = false;
      else cout << ", ";
      for ( Index isam=0; isam<vec.size(); ++isam ) {
        float val = vec[isam];
        if ( isam != 0 ) cout << ", ";
        if ( isam == 0 ) cout << "\n" << myname << setw(indent) << "[";
        else if ( 10*(isam/10) == isam ) cout << "\n" << myname << setw(indent) << " ";
        cout << setw(11) << fixed << setprecision(7) << val;
      }
      cout << "]";
    }
    cout << endl;
    cout << myname << setw(indent-2) << "]" << endl;
    cout << myname << "   ResponseCenter: " << m_ResponseCenter << endl;
    cout << myname << "       MinChannel: " << m_MinChannel << endl;
    cout << myname << "       MaxChannel: " << m_MaxChannel << endl;
  }
}

//**********************************************************************

DataMap Adc2dConvolute::updateMap(AdcChannelDataMap& acds) const {
  const string myname = "Adc2dConvolute::update: ";
  DataMap ret;
  if ( acds.size() == 0 ) return ret;
  const AdcChannelData& acd0 = acds.begin()->second;
  // Find the range of channels.
  Index icha1 = acds.begin()->first;
  if ( m_MinChannel > icha1 ) icha1 = m_MinChannel;
  Index icha2 = acds.rbegin()->first;
  if ( m_MaxChannel < icha1 ) icha2 = m_MaxChannel;
  ret.setInt("channelMin", icha1);
  ret.setInt("channelMax", icha2);
  if ( icha2 < icha1 ) return ret.setStatus(1);
  if ( m_LogLevel >= 2 ) cout << myname << "Processing run() " << acd0.run() << " event " << acd0.event()
                              << " channel " << icha1 << " to " << icha2 << endl;

  bool readBins = m_BinsPerWire;       // Read data from channel rather than local cache.
  Index nbinw = readBins ? m_BinsPerWire : 1;  // # bins per channel
  Index nbint = m_BinsPerTick;                 // # bins per output tick
  bool readSamples = nbint == 0;               // If true, input is samples instead of binSamples
  if ( readSamples ) nbint = 1;
  Index nrv = m_ResponseVectors.size();        // # of response vector (each acts on a channel bin)
  Index jofft = m_ResponseCenter;
  ret.setInt("responseVectorCount", nrv);

  // If we are reading and writing to samples, make a copy of the former.
  std::map<AdcChannel, AdcSignalVector> inSamples;
  for ( Index icha=icha1; icha<=icha2; ++icha ) {
    if ( acds.count(icha) ) inSamples[icha] = acds[icha].samples;
  }

  // Loop over output channels.
  Index nmult = 0;
  for ( Index icha=icha1; icha<=icha2; ++icha ) {
    if ( m_LogLevel >= 3 ) {
      cout << myname << "Creating samples for channel " << icha << endl;
    }
    if ( acds.count(icha) == 0 ) {
      cout << myname << "Skipping missing output channel " << icha << endl;
      continue;
    }
    AdcChannelData& acdi = acds[icha];
    AdcSignalVector& samsi = acdi.samples;
    samsi.clear();
    // Loop over response vectors.
    for ( Index irv=0; irv<nrv; ++irv ) {
      if ( m_LogLevel >= 4 ) {
        cout << myname << "  Using Response vector " << irv << endl;
      }
      const ResponseVector& resvec = m_ResponseVectors[irv];
      Index nbinrv = resvec.size();
      Index dcha = (irv + nbinw/2)/nbinw;
      Index jbinRight = (irv + nbinw/2)%nbinw;
      Index jbinLeft = nbinw - jbinRight - 1;
      using ChannelBin = std::pair<Index, Index>;
      // Use set so we don't double count when the bin count is odd.
      std::set<ChannelBin> jchbs;
      if ( icha + dcha <= icha2 ) {  // Right side.
        jchbs.emplace(icha + dcha, jbinRight);
      }
      if ( icha >= icha1 + dcha ) {  // Left side.
        jchbs.emplace(icha - dcha, jbinLeft);
      }
      // Loop over left and right channel bins.
      for ( ChannelBin jchb : jchbs ) {
        Index jcha = jchb.first;
        Index jbin = jchb.second;
        if ( m_LogLevel >= 5 ) {
          cout << myname << "    Filling from channel/bin " << jcha << "/" << jbin << endl;
        }
        if ( acds.count(jcha) == 0 ) {
          cout << myname << "      Skipping missing input channel " << jcha << endl;
          continue;
        }
        const AdcChannelData& acdj = acds[jcha];
        if ( ! readSamples ) {
          if ( acdj.binSamples.size() != m_BinsPerWire ) {
            cout << myname << "ERROR: Input channel " << jcha
                 << " has an unexpected number of sample bins: "
                 << acdj.binSamples.size() << " != " << m_BinsPerWire << endl;
            continue;
          }
        }
        const AdcSignalVector& samsj = readSamples ? inSamples[jcha] : acdj.binSamples[jbin];
        Index nsamj = samsj.size();
        if ( nsamj == 0 ) {
          cout << myname << "      Skipping empty channel/bin " << jcha << "/" << jbin << endl;
          continue;
        }
        Index nsami = (nsamj-1)/nbint + 1;
        if ( nsami > samsi.size() ) {
          if ( m_LogLevel >= 5 ) {
            cout << myname << "      Increasing output sample size to " << nsami << endl;
          }
          samsi.resize(nsami, 0.0);
        }
        // Loop over output tick bins.
        for ( Index isami=0; isami<nsami; ++isami ) {
          if ( m_LogLevel >= 6 ) {
            cout << myname << "        Filling output channel " << icha << " sample " << isami << endl;
          }
          Index ibint = nbint*isami;
          float& sami = samsi[isami];
          // Loop over input tick bins.
          // Note jbint = ibint + jofft - ijbinrv
          Index ibintOff = ibint + jofft;
          Index jbint1 = ibintOff >= nbinrv ? ibintOff - nbinrv + 1: 0;
          Index jbint2 = ibintOff < nsamj ? ibintOff + 1 : nsamj;
          if ( m_LogLevel >= 6 ) {
            cout << myname << "          Input channel tick range: (" << jbint1 << ", "
                 << jbint2 - 1 << ")" << endl;
            cout << myname << "          Response vector " << irv << " tick range: ("
                 << ibintOff + 1 - jbint2 << ", "
                 << ibintOff - jbint1 << ")" << endl;
          }
          for ( Index jbint=jbint1; jbint<jbint2; ++jbint ) {
            Index ijbinrv = ibintOff - jbint;
            sami += resvec[ijbinrv]*samsj[jbint];
            ++nmult;
          }  // End loop over input ticks
          if ( m_LogLevel >= 6 ) {
            cout << myname << "          Channel " << icha << " sample[" << isami
                 << "] = " << sami << endl;
          }
        }  // End loop over output ticks
      }  // End loop over left and right 
    }  // End loop over response vectors
    if ( m_LogLevel >= 3 ) {
      cout << myname << "Sample count for channel " << icha << " is " << samsi.size() << endl;
    }
  }  // End loop over output channels

  ret.setInt("nmult", nmult);

  return ret;
}

//**********************************************************************

DEFINE_ART_CLASS_TOOL(Adc2dConvolute)
