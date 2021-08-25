// Tpc2dDeconvolute_tool.cc

#include "Tpc2dDeconvolute.h"
#include "TMath.h"
#include <iostream>
#include <iomanip>
#include <sstream>

using std::string;
using std::cout;
using std::endl;
using std::setw;
using std::fixed;
using std::setprecision;
using std::ostringstream;

using Index = unsigned int;
using IndexArray = Fw2dFFT::IndexArray;
using FloatVector = AdcSignalVector;
using DoubleVector = std::vector<double>;
using DoubleVectorMap = std::map<Index, DoubleVector>;
using Name = std::string;
using RealData = Tpc2dRoi::DataArray;
using DftData = Tpc2dRoi::Dft;

namespace {

// Complex to string:  mag@thtdeg
string comstr(Tpc2dRoi::Dft::Complex val) {
  float mag = std::abs(val);
  float tht = 180/acos(-1.0)*std::arg(val);
  ostringstream ssout;
  if ( mag > 0.01 ) {
    ssout.setf(std::ios::fixed);
    ssout.precision(3);
  }
  ssout << "(" << mag << "@" << setw(4) << int(tht) << ")";
  return ssout.str();
}

}  // end unnamed namespace

//**********************************************************************
// Class methods.
//**********************************************************************

Tpc2dDeconvolute::Tpc2dDeconvolute(fhicl::ParameterSet const& ps)
: m_LogLevel(ps.get<int>("LogLevel")),
  m_ResponseVectors(ps.get<ResponseVectorVector>("ResponseVectors")),
  m_ResponseCenter(ps.get<int>("ResponseCenter")),
  m_FftSize(ps.get<Index>("FftSize")),
  m_InPath(ps.get<Name>("InPath")),
  m_OutPath(ps.get<Name>("OutPath")),
  m_SampleSigma(ps.get<float>("SampleSigma")),
  m_ChannelSigma(ps.get<float>("ChannelSigma")),
  m_pfft(new Fw2dFFT(m_FftSize, 1)) {
  const string myname = "Tpc2dDeconvolute::ctor: ";
  if ( m_LogLevel >= 1 ) {
    std::ios cout_state(nullptr);
    cout_state.copyfmt(std::cout);
    cout << myname << "Parameters:" << endl;
    cout << myname << "         LogLevel: " << m_LogLevel << endl;
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
    cout << myname << "          FftSize: " << m_FftSize << endl;
    cout << myname << "           InPath: " << m_InPath << endl;
    cout << myname << "          OutPath: " << m_OutPath << endl;
  }
}

//**********************************************************************

DataMap Tpc2dDeconvolute::updateTpcData(TpcData& tpd) const {
  const string myname = "Tpc2dDeconvolute::updateTpcData: ";
  DataMap ret;
  // Find input data.
  TpcData* ptpdIn = tpd.getTpcData(m_InPath);
  if ( ptpdIn == nullptr ) {
    cout << myname << "WARNING: Input TPC data not found: m_InPath." << endl;
    return ret.setStatus(1);
  }
  TpcData::Tpc2dRoiVector& inrois = ptpdIn->get2dRois();
  // Create output data object.
  TpcData* ptpdOut = tpd.addTpcData(m_OutPath);
  if ( ptpdOut == nullptr ) {
    cout << myname << "ERROR: Unable to create output TPC data at m_InPath." << endl;
    return ret.setStatus(2);
  }
  TpcData::Tpc2dRoiVector& outrois = ptpdOut->get2dRois();
  // Define DFT normalization.
  RealDftNormalization dftNorm(12);
  // Loop over 2d ROIs.
  Index xflog = 0;
  Index nroi = 0;
  for ( Tpc2dRoi& roi : inrois ) {
    // Check fft size.
    Index ncha = roi.channelSize();
    Index nsam = roi.sampleSize();
    IndexArray ndxs = {ncha, nsam};
    if ( fft().checkDataSize(ndxs) ) {
      cout << myname << "ERROR: Skipping too-large 2D ROI: " << ncha << " x " << nsam << endl;
      continue;
    }
    // Get ROI DFT.
    if ( roi.dft() == nullptr ) {
      const RealData& rdat = roi.data();
      DftData* pdft = new DftData(dftNorm, ndxs);
      fft().fftForward(rdat, *pdft, xflog);
      roi.resetDft(pdft);
    }
    if ( roi.dft() == nullptr ) {
      cout << myname << "ERROR: Forward FFT of ROI failed." << endl;
      continue;
    }
    DftData& roiDft = *roi.dft();
    // Build ncha x nsam response matrix.
    RealData resDat(ndxs);
    ResponseVector empty;
    // ... Response is periodic. Choose tick offset in first period.
    // ... Response is symmetric in channel.
    int nsamInt = nsam;
    int samoffInt = m_ResponseCenter;
    while ( samoffInt < 0 ) samoffInt += nsamInt;
    while ( samoffInt >= nsamInt ) samoffInt -= nsamInt;
    Index samoff = samoffInt;
    for ( Index dcha1=0; dcha1<(ncha+1)/2; ++dcha1 ) {
      const ResponseVector& res1d =
        dcha1 < m_ResponseVectors.size() ? m_ResponseVectors[dcha1] : empty;
      IndexVector dchas(1, dcha1);
      if ( dcha1 ) {
        Index dcha2 = ncha - dcha1;
        if ( dcha2 != dcha1 ) dchas.push_back(dcha2);
      }
      for ( Index dsam=0; dsam<nsam; ++dsam ) {
        Index isam = dsam + samoff;
        if ( isam >= nsam ) isam -= nsam;
        float rval = isam < res1d.size() ? res1d[isam] : 0.0;
        for ( Index dcha : dchas ) {
          IndexArray idxs = {dcha, dsam};
          Index idat = resDat.setValue(idxs, rval);
          if ( m_LogLevel >= 3 ) {
            cout << myname << "  res[" << dcha << ", " << dsam << "]: " << rval << endl;
          }
          // This shouldn't happen...
          if ( idat >= resDat.size() ) {
            cout << myname << "ERROR: Unexpected invalid response indices." << endl;
          }
        }
      }
    }
    // Get response DFT.
    RealDftNormalization respNorm(RealDftNormalization::convolutionNormalization());
    DftData resDft(respNorm, ndxs);
    fft().fftForward(resDat, resDft, xflog);
    // Build sample filter.
    // User supplies Gaussian sigma in sample (tick). Xform is a Gaussian
    // with sigma_freq = Nsam/(2 pi sigma_time).
    std::vector<double> samFilt(nsam, 1.0);
    if ( m_SampleSigma > 0.0 ) {
      float freqSigma = nsam/(2.0*acos(-1.0)*m_SampleSigma);
      for ( Index ifrq=0; ifrq<nsam; ++ifrq ) {
        float xf = ifrq/freqSigma;
        float fac = exp(-xf*xf/2.0);
        samFilt[ifrq] = fac;
        if ( m_LogLevel >= 3 ) cout << myname << "  samfil[" << ifrq << "]: " << fac << endl;
      }
    }
    // Build channel filter.
    // User supplies Gaussian sigma in channel. Xform is a Gaussian
    // with sigma_freq = Ncha/(2 pi sigma_time).
    std::vector<double> chaFilt(ncha, 1.0);
    if ( m_ChannelSigma > 0.0 ) {
      float freqSigma = nsam/(2.0*acos(-1.0)*m_SampleSigma);
      for ( Index ifrq=0; ifrq<ncha; ++ifrq ) {
        float xf = ifrq/freqSigma;
        float fac = exp(-xf*xf/2.0);
        chaFilt[ifrq] = fac;
        if ( m_LogLevel >= 3 ) cout << myname << "  chafil[" << ifrq << "]: " << fac << endl;
      }
    }
    // Evaluate deconvoluted DFT.
    // Copy the input ROI DFT, multiply by filter and divide by response.
    DftData* poutDft = new DftData(roiDft);
    bool logEntry = m_LogLevel>=4;
    for ( Index idat=0; idat<poutDft->size(); ++idat ) {
      auto idxs = poutDft->indexArrays(idat);
      Index dcha = idxs[0][0];
      Index dsam = idxs[0][1];
      Tpc2dRoi::Dft::Complex& val = poutDft->data()[idat];
      // For the filters, use the smaller of the wave number and its conjugate.
      Index dchaFilt = std::min(dcha, ncha-dcha);
      Index dsamFilt = std::min(dsam, nsam-dsam);
      if ( logEntry ) {
        cout << myname << "Decon[" << dcha << "," << dsam << "]";
        cout << "[" << dchaFilt << "," << dsamFilt << "]";
        cout << ": " << comstr(val);
      }
      val *= chaFilt[dchaFilt];
      val *= samFilt[dsamFilt];
      Tpc2dRoi::Dft::Complex resp = resDft.data()[idat];
      if ( resp != 0.0 ) val /= resp;
      if ( logEntry ) {
        cout << " * " << chaFilt[dchaFilt];
        cout << " * " << samFilt[dsamFilt];
        cout << " / " << comstr(resp) << " = " << comstr(val) << endl;
      }
    }
    if ( m_LogLevel >=3 ) {
      cout << myname << "   Input power: " << roiDft.power() << endl;
      cout << myname << "Response power: " << resDft.power() << endl;
      cout << myname << "  Output power: " << poutDft->power() << endl;
    }
    // Create the output ROI and attach the DFt data.
    outrois.emplace_back(ncha, nsam, roi.channelOffset(), roi.sampleOffset());
    Tpc2dRoi& outroi = outrois.back();
    outroi.resetDft(poutDft);
    // Evaluate deconvoluted reponse.
    int dstat = fft().fftBackward(*outroi.dft(), outroi.data());
    if ( dstat ) {
      cout << myname << "ERROR: Reverse xform failed with error " << dstat << endl;
      continue;
    }
    ++nroi;
  }
  ret.setInt("dcoNroiIn", inrois.size());
  ret.setInt("dcoNroiOut", nroi);
  return ret;
}

//**********************************************************************

DataMap Tpc2dDeconvolute::viewTpcData(const TpcData&) const {
  const string myname = "Tpc2dDeconvolute::updateTpcData: ";
  cout << myname << "ERROR: View of TPC data is not supported here." << endl;
  DataMap ret;
  return ret.setStatus(1);
}

//**********************************************************************

Fw2dFFT& Tpc2dDeconvolute::fft() const {
  return *m_pfft;
}

//**********************************************************************
