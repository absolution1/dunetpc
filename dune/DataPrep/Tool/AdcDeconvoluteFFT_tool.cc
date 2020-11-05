// AdcDeconvoluteFFT_tool.cc

#include "AdcDeconvoluteFFT.h"
#include "dune/DuneCommon/DuneFFT.h"
#include "dune/ArtSupport/DuneToolManager.h"
#include "TMath.h"
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
using Name = std::string;
using DFT = DuneFFT::DFT;

//**********************************************************************
// Class methods.
//**********************************************************************

AdcDeconvoluteFFT::AdcDeconvoluteFFT(fhicl::ParameterSet const& ps)
: m_LogLevel(ps.get<int>("LogLevel")),
  m_Action(ps.get<Index>("Action")),
  m_ResponseVectors(ps.get<ResponseVectorVector>("ResponseVectors")),
  m_GausFilterSigmas(ps.get<FloatVector>("GausFilterSigmas")),
  m_IndexMapTool(ps.get<Name>("IndexMapTool")),
  m_useResponse(false), m_useFilter(false),
  m_doFFTConvolute(false), m_doDirectConvolute(false), m_doDeconvolute(false) {
  const Name myname = "AdcDeconvoluteFFT:ctor: ";
  // Set the action flags.
  Name amsg;
  if ( m_Action == 0 ) {
    amsg = "No action.";
  } else if ( m_Action == 1 ) {
    m_doDeconvolute = true;
    m_useResponse = true;
    m_useFilter = true;
    amsg = "Deconvolution.";
  } else if ( m_Action == 2 ) {
    m_doFFTConvolute = true;
    m_useResponse = true;
    amsg = "FFT response convolution.";
  } else if ( m_Action == 3 ) {
    m_doFFTConvolute = true;
    m_useFilter = true;
    amsg = "FFT filter convolution.";
  } else if ( m_Action == 4 ) {
    m_doDirectConvolute = true;
    m_useResponse = true;
    amsg = "Direct response convolution.";
  } else if ( m_Action == 5 ) {
    m_doDirectConvolute = true;
    m_useFilter = true;
    amsg = "Direct filter convolution.";
  } else {
    abort();
    cout << myname << "ERROR: Invalid action flag: " << m_Action << endl;
  }
  // Fetch index mapper.
  if ( m_IndexMapTool.size() ) {
    DuneToolManager* ptm = DuneToolManager::instance();
    m_channelToIndex = ptm->getPrivate<IndexMapTool>(m_IndexMapTool);
    if ( m_channelToIndex ) {
      cout << myname << "Using channel mapping tool " << m_IndexMapTool << endl;
    } else {
      cout << myname << "WARNING: Channel mapping tool not found: " << m_IndexMapTool << endl;
    }
  }
  if ( m_LogLevel >= 1 ) {
    std::ios cout_state(nullptr);
    cout_state.copyfmt(std::cout);
    cout << myname << "Parameters:" << endl;
    cout << myname << "         LogLevel: " << m_LogLevel << endl;
    cout << myname << "           Action: " << m_Action << " (" << amsg << ")" << endl;
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
    cout << myname << "  GausFilterSigmas: [";
    first = true;
    for ( float sgm : m_GausFilterSigmas ) {
      if ( first ) first = false;
      else cout << ", ";
      cout << sgm;
    }
    cout << "]" << endl;
    cout << myname << "    ChannelMapTool: " << m_IndexMapTool << endl;;
    std::cout.copyfmt(cout_state);
  }
}

//**********************************************************************

DataMap AdcDeconvoluteFFT::update(AdcChannelData& acd) const {
  const string myname = "AdcDeconvoluteFFT::update: ";
  DataMap ret;
  Index icha = acd.channel;
  if ( m_LogLevel >= 2 ) cout << myname << "Processing run " << acd.run << " event " << acd.event
                              << " channel " << icha << endl;
  // Retrive the response vector and sigma.
  Index ivec = 0;
  if ( m_channelToIndex ) {
    ivec = m_channelToIndex->get(icha);
  }
  if ( ivec >= m_ResponseVectors.size() ) {
    cout << "WARNING: There is no response vector for index " << ivec
         << " (channel " << icha << ")." << endl;
    return ret;
  }
  const ResponseVector& responseVector = m_ResponseVectors[ivec];
  if ( ivec >= m_GausFilterSigmas.size() ) {
    cout << "WARNING: There is no Gaus filter sigma for index "
         << ivec << " (channel " << icha << ")." << endl;
    return ret;
  }
  float gausFilterSigma = m_GausFilterSigmas[ivec];
  // Do action.
  DFT::Norm fnormConv(DFT::convolutionNormalization());
  DFT::Norm fnormData(AdcChannelData::dftNormalization());
  Index rstat = 0;
  Index fftLogLevel = m_LogLevel > 2 ? m_LogLevel - 2 : 0.0;

  if ( !m_doDeconvolute && !m_doFFTConvolute && !m_doDirectConvolute ) return ret;

  bool useFFT = m_doDeconvolute || m_doFFTConvolute;

  // Input data.
  FloatVector& xdasInp = acd.samples;
  Index nsam = xdasInp.size();
  if ( nsam == 0 ) return ret;

  // Build the response sequence extended to the data length.
  FloatVector xdasRes(nsam, 0.0);
  Index nres = responseVector.size();
  if ( m_useResponse ) {
    if ( m_LogLevel >= 3 ) cout << myname << "Building time-domain response vector." << endl;
    if ( nres == 0 ) {
      cout << myname << "ERROR: Response function is undefined." << endl;
      return ret.setStatus(1);
    }
    if ( nsam < nres ) {
      cout << myname << "WARNING: Data is shorter than the response function. No action taken." << endl;
      return ret.setStatus(2);
    }
    for ( Index isam=0; isam<nres; ++isam ) {
      xdasRes[isam] = responseVector[isam];
    }
  }

  // Build response deconvolution.
  // Transform the response.
  // We may want to cache this result mapped by nsam.
  DFT dftRes(fnormConv);
  if ( m_useResponse && useFFT ) {
    if ( m_LogLevel >= 3 ) cout << myname << "Using FFT to build frequency-domain response vectors." << endl;
    rstat = DuneFFT::fftForward(xdasRes, dftRes, fftLogLevel);
    if ( rstat ) {
      cout << myname << "WARNING: Unable to xform response function. No action taken." << endl;
      return ret.setStatus(3);
    }
    if ( m_LogLevel >= 4 ) {
      cout << myname << "Frequency components:" << endl;
      for ( Index ifrq=0; ifrq<dftRes.nCompact(); ++ifrq ) {
        cout << myname << setw(4) << ifrq << ": " << setw(10) << fixed << dftRes.amplitude(ifrq);
        if ( ifrq < dftRes.nPhase() ) cout << " @ " << setw(10) << fixed << dftRes.phase(ifrq);
        cout << endl;
      }
    }
  }
  
  // Transform the data to frequency domain.
  DFT dftInp(fnormData);
  FloatVector xamsInp, xphsInp;
  if ( useFFT ) {
    if ( m_LogLevel >= 3 ) cout << myname << "Using FFT to build frequency-domain input vectors." << endl;
    FloatVector xresInp, ximsInp;
    rstat = DuneFFT::fftForward(xdasInp, dftInp, fftLogLevel);
    if ( rstat ) {
      cout << myname << "WARNING: Unable to xform input data. No action taken." << endl;
    return ret.setStatus(4);
    }
    if ( m_LogLevel >= 4 ) {
      cout << myname << "Frequency components:" << endl;
      for ( Index ifrq=0; ifrq<dftRes.nCompact(); ++ifrq ) {
        cout << myname << setw(4) << ifrq << ": " << setw(10) << fixed << dftInp.amplitude(ifrq);
        if ( ifrq < dftRes.nPhase() ) cout << " @ " << setw(10) << fixed << dftInp.phase(ifrq);
        cout << endl;
      }
    }
  }
  Index namp = dftInp.nAmplitude();
  Index npha = dftInp.nPhase();
  
  // Build the filter.
  // In frequncy space, this is just a Gaussian with sima = N/(2*pi*sigma_t)
  const float xmin = 1.e-20;
  FloatVector xdasFil(nsam, 1.0);
  DFT dftFil(fnormConv);
  if ( m_useFilter ) {
    FloatVector xamsFil(namp);
    FloatVector xphsFil(npha, 0.0);
    if ( m_LogLevel >= 3 ) cout << myname << "Building frequency-domain filter vectors." << endl;
    if ( gausFilterSigma == 0.0 ) {
      xamsFil[0] = 1.0;
    } else {
      static float novertwopi = 0.5*nsam/acos(-1);
      float sigmaFreq = novertwopi/gausFilterSigma;
      for ( Index ifrq=0; ifrq<xamsFil.size(); ++ifrq ) {
        float val = TMath::Gaus(ifrq, 0.0, sigmaFreq, false);
        xamsFil[ifrq] = val;
      }
    }
    dftFil.moveIn(xamsFil, xphsFil);
    cout << myname << "Frequency components:" << endl;
    for ( Index ifrq=0; ifrq<dftFil.nCompact(); ++ifrq ) {
      cout << myname << setw(4) << ifrq << ": " << setw(10) << fixed << dftFil.amplitude(ifrq);
      if ( ifrq < dftFil.nPhase() ) cout << " @ " << setw(10) << fixed << dftFil.phase(ifrq);
      cout << endl;
    }
  }

  // Create the output vectors.
  FloatVector xdasOut;
  FloatVector xamsOut, xphsOut;

  // Deconvolute/convolute (divide/multiply) in frequency space.
  if ( useFFT ) {
    if ( m_LogLevel >= 3 ) cout << myname << "Building frequency-domain output vectors." << endl;
    DFT dftOut(fnormData, nsam);
    for ( Index ifrq=0; ifrq<namp; ++ifrq ) {
      float xamOut = dftInp.amplitude(ifrq);
      float xphOut = dftInp.phase(ifrq);
      if ( m_doFFTConvolute ) {
        if ( m_useFilter ) {
          xamOut *= dftFil.amplitude(ifrq);
          xphOut += dftFil.phase(ifrq);
        } else if ( m_useResponse ) {
          xamOut *= dftRes.amplitude(ifrq);
          xphOut += dftRes.phase(ifrq);
        }
      } else if ( m_doDeconvolute ) {
        if ( m_useFilter ) {
          xamOut *= dftFil.amplitude(ifrq);
          xphOut += dftFil.phase(ifrq);
        }
        if ( xamOut > xmin ) {
          float den = dftRes.amplitude(ifrq);
          if ( fabs(den) < xmin ) {
            cout << myname << "WARNING: Ignoring near-zero division in deconvolution. " << endl;
          } else {
            xamOut /= den;
          }
          xphOut -= dftRes.phase(ifrq);
        }
      } else {
        cout << myname << "ERROR: Unxpected error in FFT (de)convolution." << endl;
      }
      if ( m_LogLevel >= 4 ) {
        cout << myname << setw(4) << ifrq << ": " << setw(10) << fixed << xamOut;
        if ( ifrq < npha ) cout << " @ " << setw(10) << fixed << xphOut;
        cout << endl;
      }
      dftOut.setAmplitude(ifrq, xamOut);
      if ( ifrq < npha ) dftOut.setPhase(ifrq, xphOut);
    }
  
    // Transform back.
    if ( m_LogLevel >= 3 ) cout << myname << "Using FFT to build time-domain output vector." << endl;
    rstat = DuneFFT::fftInverse(dftOut, xdasOut, fftLogLevel);
    if ( rstat ) {
      cout << myname << "WARNING: Inverse xform failed. No action taken." << endl;
      return ret.setStatus(4);
    }
  }

  // Direct convolution.
  if ( m_doDirectConvolute ) {
    if ( m_LogLevel >= 3 ) cout << myname << "Directly building time-domain output vector." << endl;
    xdasOut.resize(nsam);
    FloatVector empty;
    const FloatVector& xdasFld = m_useResponse ? xdasRes : m_useFilter ? xdasFil : empty;
    for ( Index isam=0; isam<nsam; ++isam ) {
      double xsam = 0.0;
      for ( Index jsam=0; jsam<nsam; ++jsam ) {
        Index imj = (isam + nsam - jsam) % nsam;
        xsam += xdasInp[jsam]*xdasFld[imj];
      }
      xdasOut[isam] = xsam;
    }
  }

  // Record results.
  acd.samples = xdasOut;
  //acd.dftmags = xamsOut;
  //acd.dftphases = xphsOut;

  return ret;
}

//**********************************************************************
