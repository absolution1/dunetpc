// DuneDPhase3x1x1NoiseRemovalService_service.cc

#include "DuneDPhase3x1x1NoiseRemovalService.h"
#include <iostream>
#include <sstream>
#include <iomanip>
#include <algorithm>
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "larcore/Geometry/Geometry.h"
#include "lardata/Utilities/LArFFT.h"

#include "larevt/CalibrationDBI/Interface/ChannelStatusService.h"
#include "larevt/CalibrationDBI/Interface/ChannelStatusProvider.h"

//**********************************************************************

DuneDPhase3x1x1NoiseRemovalService::
DuneDPhase3x1x1NoiseRemovalService(fhicl::ParameterSet const& pset, art::ActivityRegistry&) :
    fDoTwoPassFilter( pset.get<bool>("DoTwoPassFilter") ),
    fCoherent32( pset.get<bool>("Coherent32") ),
    fCoherent16( pset.get<bool>("Coherent16") ),
    fLowPassFlt( pset.get<bool>("LowPassFlt") ),
    fFlatten( pset.get<bool>("Flatten") ),
    fCoherent32Groups( pset.get<std::vector<size_t>>("Coherent32Groups") ),
    fCoherent16Groups( pset.get<std::vector<size_t>>("Coherent16Groups") ),
    fRoiStartThreshold( pset.get<float>("RoiStartThreshold") ),
    fRoiEndThreshold( pset.get<float>("RoiEndThreshold") ),
    fRoiPadLow( pset.get<int>("RoiPadLow") ),
    fRoiPadHigh( pset.get<int>("RoiPadHigh") ),
    fBinsToSkip( pset.get<int>("BinsToSkip") ),
    fGeometry( &*art::ServiceHandle<geo::Geometry>() ),
    fFFT( &*art::ServiceHandle<util::LArFFT>() )
{
    if (pset.get<std::string>("CorrMode") == "mean") { fMode = 1; }
    else if (pset.get<std::string>("CorrMode") == "median") { fMode = 2; }
    else
    {
        std::cout << "DuneDPhase3x1x1NoiseRemovalService WARNING: correction set to mean value." << std::endl;
        fMode = 1;
    }

    fLowPassCoeffs.resize(fFFT->FFTSize() / 2 + 1);

    const float fcut = 0.5; // [MHz]
    for (size_t i = 0; i < fLowPassCoeffs.size(); ++i)
    {
        float f = 0.0015 * i; // [MHz]
        fLowPassCoeffs[i] = 1.0 / sqrt(1.0 + pow(f/fcut, 8));
    }
}
//**********************************************************************

int DuneDPhase3x1x1NoiseRemovalService::update(AdcChannelDataMap& datamap) const
{
  const std::string myname = "DuneDPhase3x1x1NoiseRemovalService:update: ";
  if ( datamap.size() == 0 ) {
    std::cout << myname << "WARNING: No channels found." << std::endl;
    return 1;
  }

  unsigned int nsig = 0;
  bool first = true;
  for (const AdcChannelDataMap::value_type& ent : datamap)
  {
    const AdcChannelData& data = ent.second;
    if (first) { nsig = data.samples.size(); first = false; }
    else if (data.samples.size() != nsig)
    {
      std::cout << myname << "WARNING: Channels have inconsistent sample counts." << std::endl;
      return 2;
    }
  }
  
  if (nsig == 0)
  {
    std::cout << myname << "WARNING: No ADC samples found." << std::endl;
    return 3;
  }

  std::cout << myname << "Processing noise removal..." << std::endl;

  if(fDoTwoPassFilter)
  {
    if (fFlatten)
    {
      removeSlopePolynomial(datamap);
    }

    if (fCoherent16)
    {
      auto ch_groups = makeDaqGroups(16, fCoherent16Groups);
      removeCoherent(ch_groups, datamap);
    }

    if (fCoherent32)
    {
      auto ch_groups = makeDaqGroups(32, fCoherent32Groups);
      removeCoherent(ch_groups, datamap);
    }

    if (fLowPassFlt)
    {
      removeHighFreq(datamap);
    }
  }
  else
  {
    if (fCoherent32)
    {
      auto ch_groups = makeDaqGroups(32, fCoherent32Groups);
      removeCoherent(ch_groups, datamap);
    }

    if (fCoherent16)
    {
      auto ch_groups = makeDaqGroups(16, fCoherent16Groups);
      removeCoherent(ch_groups, datamap);
    }

    if (fLowPassFlt)
    {
      removeHighFreq(datamap);
    }

    if (fFlatten)
    {
      removeSlope(datamap);
    }
  }
  std::cout << myname << "...done." << std::endl;

  return 0;
}
//**********************************************************************

std::vector<float> DuneDPhase3x1x1NoiseRemovalService::getMeanCorrection(
    const std::vector<unsigned int> & channels,
    const AdcChannelDataMap & datamap) const
{
  size_t n_samples = datamap.begin()->second.samples.size();
  std::vector<size_t> ch_averaged(n_samples, 0);
  std::vector<float> correction(n_samples, 0);

  for (unsigned int ch : channels)
  {
      auto iacd = datamap.find(ch);
      if (iacd == datamap.end()) continue;

      const AdcChannelData & adc = iacd->second;
      auto mask = roiMask(adc);

      for (size_t s = 0; s < n_samples; ++s)
      {
          if (!mask[s]) { continue; }

          AdcFlag flag = adc.flags.size() ? adc.flags[s] : AdcGood;
          if (flag != AdcGood) { continue; }

          correction[s] += adc.samples[s];
          ch_averaged[s]++;
      }
  }
  for (size_t s = 0; s < n_samples; ++s)
  {
      if (ch_averaged[s] > 0) { correction[s] /= ch_averaged[s]; }
  }
  return correction;
}

std::vector<float> DuneDPhase3x1x1NoiseRemovalService::getMedianCorrection(
    const std::vector<unsigned int> & channels,
    const AdcChannelDataMap & datamap) const
{
  size_t n_samples = datamap.begin()->second.samples.size();
  std::vector< std::vector<float> > samples(n_samples);

  for (unsigned int ch : channels)
  {
      auto iacd = datamap.find(ch);
      if (iacd == datamap.end()) continue;

      const AdcChannelData & adc = iacd->second;
      auto mask = roiMask(adc);

      for (size_t s = 0; s < n_samples; ++s)
      {
          if (!mask[s]) { continue; }

          AdcFlag flag = adc.flags.size() ? adc.flags[s] : AdcGood;
          if (flag != AdcGood) { continue; }

          samples[s].push_back(adc.samples[s]);
      }
  }

  std::vector<float> correction(n_samples);
  for (size_t s = 0; s < n_samples; ++s)
  {
      size_t n = samples[s].size();
      if (n < 2) { correction[s] = 0; continue; }

      std::sort(samples[s].begin(), samples[s].end());

      if ((n % 2) == 0) { correction[s] = 0.5 * (samples[s][n/2] + samples[s][(n/2)-1]); }
      else              { correction[s] = samples[s][(n-1)/2]; }
  }
  return correction;
}

void DuneDPhase3x1x1NoiseRemovalService::removeCoherent(const GroupChannelMap & ch_groups, AdcChannelDataMap& datamap) const
{
  if (datamap.empty()) return;

  size_t n_samples = datamap.begin()->second.samples.size();
  std::vector<double> correction(n_samples);

  for (const auto & entry : ch_groups)
  {
    const auto & channels = entry.second;
    std::vector<float> correction;

    if (fMode == 1) // mean
    {
        correction = getMeanCorrection(channels, datamap);
    }
    else if (fMode == 2) // median
    {
        correction = getMedianCorrection(channels, datamap);
    }

    for (unsigned int ch : channels)
    {
        auto iacd = datamap.find(ch);
        if (iacd == datamap.end()) continue;

        AdcChannelData & adc = iacd->second;
        for (size_t s = 0; s < n_samples; ++s)
        {
            adc.samples[s] -= correction[s];
        }

        if (adc.samples[2] - adc.samples[1] > 20) // ugly fix of the two ticks in plane0
        {
            adc.samples[0] = adc.samples[2];
            adc.samples[1] = adc.samples[2];
        }
    }
  }
}
//**********************************************************************

void DuneDPhase3x1x1NoiseRemovalService::removeHighFreq(AdcChannelDataMap& datamap) const
{
    auto const & chStatus = art::ServiceHandle< lariov::ChannelStatusService >()->GetProvider();

    for (auto & entry : datamap)
    {
        if (chStatus.IsPresent(entry.first) && !chStatus.IsNoisy(entry.first)) { fftFltInPlace(entry.second.samples, fLowPassCoeffs); }
    }
}
//**********************************************************************

void DuneDPhase3x1x1NoiseRemovalService::removeSlope(AdcChannelDataMap& datamap) const
{
  size_t n_samples = datamap.begin()->second.samples.size(); 
  std::vector< float > slope(n_samples);

  auto const & chStatus = art::ServiceHandle< lariov::ChannelStatusService >()->GetProvider();

  for (auto & entry : datamap)
  {
    if (!chStatus.IsPresent(entry.first) || chStatus.IsNoisy(entry.first)) { continue; }

    auto & adc = entry.second.samples;
    auto mask = roiMask(entry.second);

    bool start = true;
    float i = 0, s0 = 0, s1 = 0;
    while (i < n_samples)
    {
        if (mask[i])
        {
            if (start) { s0 = adc[i]; start = false; }
            else { s0 = 0.8 * s0 + 0.2 * adc[i]; } // average over some past adc
            slope[i] = adc[i]; ++i;
        }
        else
        {
            size_t j = i;
            while ((j < n_samples) && !mask[j]) { ++j; }
            if (j < n_samples)
            {
                size_t k = 1;
                float f = 1, g = f;
                s1 = adc[j];
                while ((k < 25) && (j+k < n_samples) && mask[j+k])
                {
                    f *= 0.8; g += f; s1 += f * adc[j+k]; ++k;
                }
                s1 /= g; // average ofer following adc
            }
            else { s1 = s0; }

            float ds = (s1 - s0) / (j - i + 1);

            j = i;
            while ((j < n_samples) && !mask[j])
            {
                slope[j] = s0 + (j - i + 1) * ds;
                ++j;
            }
            start = true;
            i = j;
        }
    }

    double y, sx = 0, sy = 0, sxy = 0, sx2 = 0, sy2 = 0;
    for (size_t s = 0; s < n_samples; ++s)
    {
        y = slope[s]; sx += s; sy += y; sxy += s*y; sx2 += s*s; sy2 += y*y;
    }

    double ssx = sx2 - ((sx * sx) / n_samples);
    double c = sxy - ((sx * sy) / n_samples);
    double mx = sx / n_samples;
    double my = sy / n_samples;
    double b = my - ((c / ssx) * mx);
    double a = c / ssx;
    for (size_t s = 0; s < n_samples; ++s) { adc[s] -= (a*s + b); }
  }
}
//**********************************************************************

void DuneDPhase3x1x1NoiseRemovalService::removeSlopePolynomial(AdcChannelDataMap& datamap) const
{
  size_t n_samples = datamap.begin()->second.samples.size();  //1667
  std::vector< float > slope(n_samples);

  auto const & chStatus = art::ServiceHandle< lariov::ChannelStatusService >()->GetProvider();

  for (auto & entry : datamap)
  {
    if (!chStatus.IsPresent(entry.first) || chStatus.IsNoisy(entry.first)) { continue; }

    auto & adc = entry.second.samples;
    auto & signal = entry.second.signal;

    //preparation for plynomial fit
    int fOrder = 3; //order of the polynomial to be fitted to baseline
    std::vector<long double> line(fOrder+2,0.);
    std::vector< std::vector< long double> > matrix(fOrder+1,line);

    for(int i = 0; i < fOrder+1; i++)
    {
      for(int j = 0; j < fOrder+2; j++)
      {
        matrix[i][j] = 0.;
      }
    }

    double x, y;
    AdcIndex np = 0;
    std::vector<double> sol;
    sol.resize(fOrder+1, 0.);


    for(AdcIndex i = fBinsToSkip; i < n_samples; i++)
    {
      if(signal[i])
      { 
	continue;
      }

      x = i;
      y = adc[i];

      int j = 0;
      for( int l = 1; l <= fOrder+1; l++)
      {
        long double wt = pow(x, l-1);
        int k = 0;

        for(int m=1; m <= l; m++)
        {
          matrix[j][k] += wt*pow(x,m-1);
          k++;
        }

        matrix[j][fOrder+1] += y*wt;
        j++;
      }
      np++;
    }

    for(int i = 1; i < fOrder+1; i++){
      for(int j = 0; j<i; j++){
        matrix[j][i] = matrix[i][j];
      }
    }   

    //get parameters for polynomial with Gauss-Jordan
    if(np > n_samples/10) //need minimum number of bins for fit
    {
      sol = GaussJordanSolv(matrix);

       //subtract fit from waveform
      for(AdcIndex i = 0; i < n_samples; i++)
      {
        double corr = 0.;
        for(size_t a = 0; a < sol.size(); a++)
	{
          corr += sol[a]*pow(i,a);
	}
        adc[i] -= corr;
      }
    }
  }
}
//**********************************************************************


void DuneDPhase3x1x1NoiseRemovalService::fftFltInPlace(std::vector< float > & adc, const std::vector< float > & coeffs) const
{
  std::vector< TComplex > ch_spectrum(fFFT->FFTSize() / 2 + 1);
  std::vector< float > ch_waveform(fFFT->FFTSize(), 0);

  size_t n_samples = adc.size();

  std::copy(adc.begin(), adc.end(), ch_waveform.begin());
  for (size_t s = n_samples; s < ch_waveform.size(); ++s)
  {
      ch_waveform[s] = ch_waveform[2*n_samples - s - 1];
  }
  fFFT->DoFFT(ch_waveform, ch_spectrum);
  for (size_t c = 0; c < coeffs.size(); ++c)
  {
      ch_spectrum[c] *= coeffs[c];
  }
  fFFT->DoInvFFT(ch_spectrum, ch_waveform);

  std::copy(ch_waveform.begin(), ch_waveform.begin()+n_samples, adc.begin());
}
//**********************************************************************

std::vector< float > DuneDPhase3x1x1NoiseRemovalService::fftFlt(const std::vector< float > & adc, const std::vector< float > & coeffs) const
{
  std::vector< TComplex > ch_spectrum(fFFT->FFTSize() / 2 + 1);
  std::vector< float > ch_waveform(fFFT->FFTSize(), 0);

  size_t n_samples = adc.size();

  std::copy(adc.begin(), adc.end(), ch_waveform.begin());
  for (size_t s = n_samples; s < ch_waveform.size(); ++s)
  {
      ch_waveform[s] = ch_waveform[2*n_samples - s - 1];
  }
  fFFT->DoFFT(ch_waveform, ch_spectrum);
  for (size_t c = 0; c < coeffs.size(); ++c)
  {
      ch_spectrum[c] *= coeffs[c];
  }
  fFFT->DoInvFFT(ch_spectrum, ch_waveform);

  std::vector< float > flt_adc(n_samples);
  std::copy(ch_waveform.begin(), ch_waveform.begin()+n_samples, flt_adc.begin());
  return flt_adc;
}
//**********************************************************************

std::vector<bool> DuneDPhase3x1x1NoiseRemovalService::roiMask(const AdcChannelData & adc) const
{
  std::vector<bool> mask(adc.samples.size(), true);

  auto adc_flt = fftFlt(adc.samples, fLowPassCoeffs);

  bool inroi = false;
  for (int i = 0; i < (int)adc_flt.size(); ++i)
  {
    auto sig = adc_flt[i];
    if (inroi)
    {
      if (sig > fRoiEndThreshold) { mask[i] = false; }
      else
      {
        for (int p = 0; p <= fRoiPadHigh; ++p) { if ((i + p) < (int)mask.size()) { mask[i + p] = false; } }
        inroi = false;
      }
    }
    else
    {
      if (sig > fRoiStartThreshold )
      {
        for (int p = fRoiPadLow; p >= 0; --p) { if (i - p >= 0) { mask[i - p] = false; } }
        inroi = true;
      }
    }
  }
  return mask;
}
//**********************************************************************

GroupChannelMap DuneDPhase3x1x1NoiseRemovalService::makeDaqGroups(size_t gsize, const std::vector< size_t > & gidx) const
{
  GroupChannelMap groups;

  auto const & chStatus = art::ServiceHandle< lariov::ChannelStatusService >()->GetProvider();

  const unsigned int nchan = fGeometry->Nchannels();
  for (unsigned int ch = 0; ch < nchan; ++ch)
  {
    size_t g = get311Chan(ch) / gsize;
    //std::cout << "p:" << fGeometry->View(raw::ChannelID_t(ch)) << " g:" << g << " ch:" << ch << std::endl;
    if (gidx.empty() || has(gidx, g))
    {
        if (chStatus.IsPresent(ch) && !chStatus.IsNoisy(ch)) { groups[g].push_back(ch); }
    }
  }

  return groups;
}
//**********************************************************************

GroupChannelMap DuneDPhase3x1x1NoiseRemovalService::makeGroups(size_t gsize, const std::vector< size_t > & gidx) const
{
  GroupChannelMap groups;

  auto const & chStatus = art::ServiceHandle< lariov::ChannelStatusService >()->GetProvider();

  const unsigned int nchan = fGeometry->Nchannels();
  for (unsigned int ch = 0; ch < nchan; ++ch)
  {
    size_t g = ch / gsize;
    //std::cout << "p:" << fGeometry->View(raw::ChannelID_t(ch)) << " g:" << g << " ch:" << ch << std::endl;
    if (gidx.empty() || has(gidx, g))
    {
        if (chStatus.IsPresent(ch) && !chStatus.IsNoisy(ch)) { groups[g].push_back(ch); }
    }
  }

  return groups;
}
//**********************************************************************

size_t DuneDPhase3x1x1NoiseRemovalService::get311Chan(size_t LAr_chan)
{
  size_t crate = LAr_chan / 320;
  size_t Chan311;

  LAr_chan = 8*(LAr_chan/8+1)-LAr_chan%8 -1;

  if(crate == 0)
  {
    LAr_chan = 32*(LAr_chan/32+1)-LAr_chan%32 -1;
    size_t card = 4 - ((LAr_chan / 32) % 5);
    if(LAr_chan > 159)
    {
        size_t shift = 31 - (LAr_chan % 32);
        Chan311 = (2*card)*32 + shift;
    }
    else
    {
       size_t shift = 31 - (LAr_chan % 32);
       Chan311 = (2*card + 1)*32 + shift;
    }
  }
  else
  {
     size_t new_LAr_chan = LAr_chan - crate*320;
     size_t card = ((new_LAr_chan / 32) % 5);
     if(new_LAr_chan > 159)
     {
        size_t shift = new_LAr_chan % 32;
        Chan311 = (2*card)*32 + shift;
     }
     else
     {
       size_t shift = new_LAr_chan % 32;
       Chan311 = (2*card + 1)*32 + shift;
     }
     Chan311 = Chan311 + crate*320;
  } // end of if/else statementi

  return Chan311;
}
//**********************************************************************

std::vector<double> DuneDPhase3x1x1NoiseRemovalService::GaussJordanSolv(std::vector< std::vector<long double> > matrix) const
{

  int n = matrix.size();
  
  for (int i=0; i<n; i++) {
    // Search for maximum in this column
    double maxEl = std::abs(matrix[i][i]);
    int maxRow = i;
    for (int k=i+1; k<n; k++) {
      if (std::abs(matrix[k][i]) > maxEl) {
        maxEl = std::abs(matrix[k][i]);
        maxRow = k;
      }
    }

       // Swap maximum row with current row (column by column)
    for (int k=i; k<n+1;k++) {
      double tmp = matrix[maxRow][k];
      matrix[maxRow][k] = matrix[i][k];
      matrix[i][k] = tmp;
    }

    // Make all rows below this one 0 in current column
    for (int k=i+1; k<n; k++) {
      double c = -matrix[k][i]/matrix[i][i];
      for (int j=i; j<n+1; j++) {
        if (i==j) {
          matrix[k][j] = 0;
        } else {
          matrix[k][j] += c * matrix[i][j];
        }
      }
    }
  }

  // Solve equation Ax=b for an upper triangular matrix A
  std::vector<double> x(n);
  for (int i=n-1; i>=0; i--) {
    x[i] = matrix[i][n]/matrix[i][i];
    for (int k=i-1;k>=0; k--) {
      matrix[k][n] -= matrix[k][i] * x[i];
    }
  }
  return x;  
}
//**********************************************************************

std::ostream& DuneDPhase3x1x1NoiseRemovalService::print(std::ostream& out, std::string prefix) const
{
  out << prefix << "DuneDPhase3x1x1NoiseRemovalService:  ...info" << std::endl;
  return out;
}
//**********************************************************************

DEFINE_ART_SERVICE_INTERFACE_IMPL(DuneDPhase3x1x1NoiseRemovalService, AdcNoiseRemovalService)

//**********************************************************************
