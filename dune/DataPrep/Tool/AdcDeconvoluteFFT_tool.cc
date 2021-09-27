// AdcDeconvoluteFFT_tool.cc

#include "AdcDeconvoluteFFT.h"
#include "dune/DuneCommon/Utility/DuneFFT.h"
#include "dune/ArtSupport/DuneToolManager.h"
#include "TMath.h"
#include "TDecompLU.h"
#include "TDecompSVD.h"
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
using Name = std::string;
using DFT = DuneFFT::DFT;

//**********************************************************************
// Local functions.
//**********************************************************************

namespace {

// Fill a matrix with a vector.
// The matrix is near diagonal with the near-diagonal columns populated with
// the entries in the vector vrcol.
// The diagonal element is vres[ishift].
//
// vec - Input smearing vector.
// ishift - vector offset in matrix
// mres - Output matrix.
// smsgpre - If not blank, a message prefixed with this string is logged for each row.
template<class RootMatrix>
int responseVectorToMatrix(const DoubleVector& vrcol, Index ishiftin, RootMatrix& mres,
                           int logLevel, string smsgpre ="") {
  string myname = smsgpre + "reponseVectorToMatrix: ";
  // Flip vector so we can fill rows instead of columns.
  Index nvec = vrcol.size();
  Index nrow = mres.GetNrows();
  Index ncol = mres.GetNcols();
  if ( nvec == 0 ) {
    if ( logLevel > 0 ) cout << myname << "WARNING: Response vector is empty." << endl;
    return 1;
  }
  if ( ishiftin >= nvec ) {
    if ( logLevel > 0 ) cout << myname << "WARNING: Shift " << ishiftin
                             << " is outside vector range (0, " << nvec << "]" << endl;
    return 2;
  }
  if ( logLevel >= 4 ) {
    cout << myname << "Matrix is " << nrow << " x " << ncol << endl;
    cout << myname << "Vector has length " << nvec << " and shift " << ishiftin << endl;
  }
  Index ishift = nvec - 1 - ishiftin;
  DoubleVector vrrow(nvec, 0.0);
  for ( Index ircol=0, irrow=nvec-1; ircol<nvec; ++ircol, --irrow ) {
    vrrow[irrow] = vrcol[ircol];
  }
  // Populate matrix with vector.
  if ( nvec == 0 ) return 0;
  for ( Index irow=0; irow<ncol; ++irow ) {
    Index icol = irow > ishift ? irow - ishift : 0;
    Index ivec = icol + ishift - irow;
    Index nvecInsert = icol + nvec < nrow ? nvec - ivec : nrow - icol;
    if ( logLevel >= 4 ) {
      cout << smsgpre << "Column " << icol << ": inserting " << nvecInsert
           << " entries at row " << irow << endl;
    }
    mres.InsertRow(irow, icol, &vrrow[ivec], nvecInsert);
  }
  return 0;
}

template<class RootMatrix>
int responseVectorToMatrix(const FloatVector& vresf, Index ishift, RootMatrix& mres,
                           int logLevel, string smsgpre ="") {
  DoubleVector vresd(vresf.begin(), vresf.end());
  return responseVectorToMatrix(vresd, ishift, mres, logLevel, smsgpre);
}

}  // End unnamed namepsace.

//**********************************************************************
// Class methods.
//**********************************************************************

AdcDeconvoluteFFT::AdcDeconvoluteFFT(fhicl::ParameterSet const& ps)
: m_LogLevel(ps.get<int>("LogLevel")),
  m_Action(ps.get<Index>("Action")),
  m_ResponseVectors(ps.get<ResponseVectorVector>("ResponseVectors")),
  m_ResponseCenters(ps.get<IndexVector>("ResponseCenters")),
  m_SmoothVectors(ps.get<ResponseVectorVector>("SmoothVectors")),
  m_SmoothScales(ps.get<FloatVector>("SmoothScales")),
  m_GausFilterSigmas(ps.get<FloatVector>("GausFilterSigmas")),
  m_LowFilterWidths(ps.get<FloatVector>("LowFilterWidths")),
  m_LowFilterPowers(ps.get<FloatVector>("LowFilterPowers")),
  m_IndexMapTool(ps.get<Name>("IndexMapTool")),
  m_useResponse(false), m_useFilter(false),
  m_doFFTConvolute(false), m_doDirectConvolute(false),
  m_doDeconvolute(false), m_directDeconvolute(0) {
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
  } else if ( m_Action == 6 ) {
    m_directDeconvolute = 1;
    amsg = "Direct matrix deconvolution.";
  } else if ( m_Action == 7 ) {
    m_directDeconvolute = 2;
    amsg = "Filtered matrix deconvolution.";
  } else if ( m_Action == 8 ) {
    m_directDeconvolute = 3;
    amsg = "Chi-square matrix deconvolution.";
  } else if ( m_Action == 9 ) {
    m_directDeconvolute = 4;
    amsg = "Direct deconvolution.";
  } else {
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
    cout << myname << "   ResponseShifts: [";
    first = true;
    for ( Index val : m_ResponseCenters ) {
      if ( first ) first = false;
      else cout << ", ";
      cout << val;
    }
    cout << "]" << endl;
    cout << myname << "    SmoothVectors: [";
    indent = 22;
    first = true;
    for ( const ResponseVector& vec : m_SmoothVectors ) {
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
    cout << myname << "     SmoothScales: [";
    first = true;
    for ( float val : m_SmoothScales ) {
      if ( first ) first = false;
      else cout << ", ";
      cout << setw(11) << fixed << setprecision(7) << val;
    }
    cout << "]" << endl;
    cout << myname << " GausFilterSigmas: [";
    first = true;
    for ( float sgm : m_GausFilterSigmas ) {
      if ( first ) first = false;
      else cout << ", ";
      cout << sgm;
    }
    cout << "]" << endl;
    cout << myname << "  LowFilterWidths: [";
    first = true;
    for ( float val : m_LowFilterWidths ) {
      if ( first ) first = false;
      else cout << ", ";
      cout << val;
    }
    cout << "]" << endl;
    cout << myname << "  LowFilterPowers: [";
    first = true;
    for ( float val : m_LowFilterPowers ) {
      if ( first ) first = false;
      else cout << ", ";
      cout << val;
    }
    cout << "]" << endl;
    cout << myname << "   ChannelMapTool: " << m_IndexMapTool << endl;;
    std::cout.copyfmt(cout_state);
  }
}

//**********************************************************************

DataMap AdcDeconvoluteFFT::update(AdcChannelData& acd) const {
  const string myname = "AdcDeconvoluteFFT::update: ";
  DataMap ret;
  Index icha = acd.channel();
  if ( m_LogLevel >= 2 ) cout << myname << "Processing run " << acd.run() << " event " << acd.event()
                              << " channel " << icha << endl;
  // Retrieve the response vector and sigma.
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
  Index ishift = ivec < m_ResponseCenters.size() ? m_ResponseCenters[ivec] : 0;
  if ( ivec >= m_GausFilterSigmas.size() ) {
    cout << "WARNING: There is no Gaus filter sigma for index "
         << ivec << " (channel " << icha << ")." << endl;
    return ret;
  }
  const ResponseVector empty;
  const ResponseVector& smoothVector = ivec < m_SmoothVectors.size() ? m_SmoothVectors[ivec] : empty;
  float smoothScale = ivec < m_SmoothScales.size() ? m_SmoothScales[ivec] : 1.0;
  float gausFilterSigma = m_GausFilterSigmas[ivec];
  float lfwidth = -1.0;
  float lfpower = 0;
  if ( ivec >= m_LowFilterWidths.size() ) {
    cout << "WARNING: There is no low-filter width for index "
         << ivec << " (channel " << icha << ")." << endl;
  } else {
    lfwidth = m_LowFilterWidths[ivec];
    if ( lfwidth > 0.0 && ivec >= m_LowFilterPowers.size() ) {
      cout << "WARNING: There is no low-filter width for index "
           << ivec << " (channel " << icha << ")." << endl;
      return ret;
    }
    lfpower = m_LowFilterPowers[ivec];
  }
  // Do action.
  DFT::Norm fnormConv(RealDftNormalization::convolutionNormalization());
  DFT::Norm fnormData(AdcChannelData::dftNormalization());
  Index rstat = 0;
  Index fftLogLevel = m_LogLevel > 2 ? m_LogLevel - 2 : 0.0;

  if ( !m_doDeconvolute && !m_doFFTConvolute && !m_doDirectConvolute && !m_directDeconvolute ) return ret;

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
      Index jsam = isam >= ishift ? isam - ishift : nsam + isam - ishift;
      xdasRes[jsam] = responseVector[isam];
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
  // In frequncy space, this is just a Gaussian with sigma = N/(2*pi*sigma_t)
  const float xmin = 1.e-20;
  FloatVector xdasFil(nsam, 1.0);
  DFT dftFil(fnormConv);
  if ( m_useFilter ) {
    FloatVector xamsFil(namp, 1.0);
    FloatVector xphsFil(npha, 0.0);
    if ( m_LogLevel >= 3 ) cout << myname << "Building frequency-domain filter vectors." << endl;
    if ( gausFilterSigma > 0.0 ) {
      static float novertwopi = 0.5*nsam/acos(-1);
      float sigmaFreq = novertwopi/gausFilterSigma;
      for ( Index ifrq=0; ifrq<xamsFil.size(); ++ifrq ) {
        float val = TMath::Gaus(ifrq, 0.0, sigmaFreq, false);
        xamsFil[ifrq] = val;
      }
    }
    if ( lfwidth == 0.0 ) xamsFil[0] = 0.0;
    if ( lfwidth > 0.0 ) {
      double fac = lfwidth/nsam;
      bool square = lfpower == 2.0;
      for ( Index ifrq=0; ifrq<xamsFil.size(); ++ifrq ) {
        double term = fac*ifrq;
        double termp = square ? term*term : std::pow(term, lfpower);
        double fac = termp/(1 + termp);
        xamsFil[ifrq] *= fac;
      }
    }
    dftFil.moveIn(xamsFil, xphsFil);
    if ( m_LogLevel >= 4 ) {
      cout << myname << "Filter frequency components:" << endl;
      for ( Index ifrq=0; ifrq<dftFil.nCompact(); ++ifrq ) {
        cout << myname << setw(4) << ifrq << ": " << setw(10) << fixed << dftFil.amplitude(ifrq);
        if ( ifrq < dftFil.nPhase() ) cout << " @ " << setw(10) << fixed << dftFil.phase(ifrq);
        cout << endl;
      }
    }
  }

  // Create the output vectors.
  FloatVector xdasOut;
  FloatVector xamsOut, xphsOut;

  // Deconvolute/convolute (divide/multiply) in frequency space.
  if ( useFFT ) {
    xamsOut.resize(namp, 0.0);
    xphsOut.resize(npha, 0.0);
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
      xamsOut[ifrq] = xamOut;
      if ( ifrq < npha ) {
        dftOut.setPhase(ifrq, xphOut);
        xphsOut[ifrq] = xphOut;
      }
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

  // Deconvolution with matrix solve.
  if ( m_directDeconvolute ) {
    if ( m_LogLevel >= 1 ) cout << myname << "Deconvoluting with Root: option "
                                << m_directDeconvolute << "." << endl;
    // Construct the response matrix from the response vector.
    TMatrixDSparse mres(nsam, nsam);
    if ( responseVectorToMatrix(responseVector, ishift, mres, m_LogLevel, myname) ) return ret.setStatus(6);
    // Build data vector.
    TVectorD vdat(nsam);
    for ( Index isam=0; isam<nsam; ++isam ) {
      vdat[isam] = xdasInp[isam];
    }
    // The following are alternatives to direct solution of vdat = mres q.
    // The direct solution doesn't work because mres is ill conditioned.
    // Either mres or mres and vdat are updated.
    if ( m_directDeconvolute == 1 ) {
      cout << "Performing direct matrix deconvolution." << endl;
    } else if ( m_directDeconvolute == 2 ) {
      cout << "Performing filtered matrix deconvolution." << endl;
      // Filtered deconvolution.
      // vdat = (mres mfil_inv) (mfil q)
      // Doesn't work because mfil is ill conditioned.
      // Build smearing vector vfil and matrix mfil.
      DoubleVector vfil(smoothVector.begin(), smoothVector.end());
      Index nfil = vfil.size();
      TMatrixDSparse mfil(nsam, nsam);
      if ( responseVectorToMatrix(vfil, nfil/2, mfil, m_LogLevel, myname) ) return ret.setStatus(7);
      // Solve R_t = S_t (R Sinv)_t to get RSinv.
      TMatrixDSparse mfil_t(nsam,nsam);
      mfil_t.Transpose(mfil);
      TMatrixD rsinv_t(nsam,nsam);
      rsinv_t.Transpose(mres);
      TDecompSVD ssolver(mfil_t);
      if ( ssolver.Decompose() ) {
        cout << myname << "Filtering condition is " << ssolver.Condition() << endl;
        cout << myname << "S_t is " << (ssolver.GetMatrix().IsValid() ? "" : "not ") << "valid."  << endl;
        cout << myname << "R_t is " << (rsinv_t.IsValid() ? "" : "not ") << "valid."  << endl;
        ssolver.MultiSolve(rsinv_t);
        cout << myname << "(R Sinv)_t is " << (rsinv_t.IsValid() ? "" : "not ") << "valid."  << endl;
        // Replace R with RSinv.
        mres.Transpose(rsinv_t);
      } else {
        cout << myname << "Unable to decompose smearing matrix!" << endl;
        return ret.setStatus(8);
      }
    // y = (R_t R + S_t S_t S) x
    } else if ( m_directDeconvolute == 3 ) {
      cout << "Performing chi-square matrix deconvolution." << endl;
      TMatrixDSparse mres_t(nsam, nsam);
      mres_t.Transpose(mres);
      TMatrixDSparse mresres(mres_t);
      mresres *= mres;
      vdat *= mres_t;
      vdat.Print();
      // Fetch smearing vector vsmr and build matrix msmr.
      DoubleVector vsmr(smoothVector.begin(), smoothVector.end());
      for ( double& val : vsmr ) val*= smoothScale;
      Index nsmr = vsmr.size();
      TMatrixDSparse msmr(nsam, nsam);
      if ( responseVectorToMatrix(vsmr, nsmr/2, msmr, m_LogLevel, myname) ) return ret.setStatus(9);
      TMatrixDSparse msmr_t(nsam, nsam);
      msmr_t.Transpose(msmr);
      TMatrixDSparse msmrsmr(msmr_t);
      msmrsmr *= msmr;
      // Add smearing to reponse.
      mres.SetSparseIndexAB(mresres, msmrsmr);
      mres = mresres + msmrsmr;
    } else if ( m_directDeconvolute == 4 ) {
      DoubleVector vsmr = {0.05, 0.20, 0.5, 0.20, 0.05};
      Index nsmr = vsmr.size();
      TMatrixDSparse msmr(nsam, nsam);
      if ( responseVectorToMatrix(vsmr, nsmr/2, msmr, m_LogLevel, myname) ) return ret.setStatus(9);
      //vdat *= msmr;
      msmr *= mres;
      mres.SetSparseIndex(msmr);
      mres = msmr;
    } else if ( m_directDeconvolute != 1 ) {
      return ret.setStatus(10);
    }
    // Solve vdat = mres q
    //TDecompLU solver(mres);
    TDecompSVD solver(mres);
    if ( solver.Decompose() ) {
      cout << myname << "Final matrix condition is " << solver.Condition() << endl;
      solver.Solve(vdat);
      for ( Index isam=0; isam<nsam; ++isam ) {
        acd.samples[isam] = vdat[isam];
      }
    } else {
      cout << myname << "Unable to decompose response matrix!" << endl;
      return ret.setStatus(11);
    }
  } else {
    // Record results from any other operation.
    acd.samples = xdasOut;
    if ( xamsOut.size() ) {
      acd.dftmags = xamsOut;
      acd.dftphases = xphsOut;
    }
  }

  return ret;
}

//**********************************************************************

DEFINE_ART_CLASS_TOOL(AdcDeconvoluteFFT)
