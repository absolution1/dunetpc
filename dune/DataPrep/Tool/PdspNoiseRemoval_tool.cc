// PdspNoiseRemoval_tool.cc
// Jingbo Wang (UC Davis)
// Email: jiowang@ucdavis.edu
// January 2019

#include "PdspNoiseRemoval.h"
#include <iostream>
#include "lardataobj/RawData/raw.h" 
#include "dune/ArtSupport/DuneToolManager.h"
#include "larcore/Geometry/Geometry.h"
#include "lardata/Utilities/LArFFT.h"
#include "larevt/CalibrationDBI/Interface/ChannelStatusService.h"
#include "larevt/CalibrationDBI/Interface/ChannelStatusProvider.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "dune-raw-data/Services/ChannelMap/PdspChannelMapService.h"
using std::string;
using std::cout;
using std::endl;
	
//**********************************************************************
// Class methods.
//**********************************************************************

PdspNoiseRemoval::PdspNoiseRemoval(fhicl::ParameterSet const& ps)
: m_LogLevel(ps.get<int>("LogLevel")) {
  const string myname = "PdspNoiseRemoval::ctor: ";
  
  fRemoveHighFrequency     = ps.get<bool>("RemoveHighFrequency");
  fRemoveCoherent          = ps.get<bool>("RemoveCoherent");
  fCutoffFrequency         = ps.get<float>("CutoffFrequency");
  fCoherentOffline16       = ps.get<bool>("CoherentOffline16");
  fCoherentDaq8            = ps.get<bool>("CoherentDaq8");
  fCoherentDaq16           = ps.get<bool>("CoherentDaq16");
  fCoherentFEMB128         = ps.get<bool>("CoherentFEMB128");
  fCoherentOffline16Groups = ps.get<std::vector<size_t>>("CoherentOffline16Groups");
  fCoherentDaq8Groups      = ps.get<std::vector<size_t>>("CoherentDaq8Groups");
  fCoherentDaq16Groups     = ps.get<std::vector<size_t>>("CoherentDaq16Groups");
  fCoherentFEMB128Groups   = ps.get<std::vector<size_t>>("CoherentFEMB128Groups");
  fUseBasicROIForCNR       = ps.get<bool>("UseBasicROIForCNR");
  fRoiStartThreshold       = ps.get<float>("RoiStartThreshold");
  fRoiEndThreshold         = ps.get<float>("RoiEndThreshold");
  fRoiPadLow               = ps.get<int>("RoiPadLow");
  fRoiPadHigh              = ps.get<int>("RoiPadHigh");
  
  fGeometry = &*art::ServiceHandle<geo::Geometry>();
  fFFT = &*art::ServiceHandle<util::LArFFT>(); 
  if(ps.get<std::string>("CorrMode") == "mean") { fMode = 1; }
  else if(ps.get<std::string>("CorrMode") == "median") { fMode = 2; }
  else {
    std::cout << "PdspNoiseRemoval WARNING: correction set to mean value." << std::endl;
    fMode = 1;
  }
  auto const *fDetProp = lar::providerFrom<detinfo::DetectorPropertiesService>();
  double binWidth = 1.0/(fFFT->FFTSize()*fDetProp->SamplingRate()*1.0e-6);
  fLowPassCoeffs.resize(fFFT->FFTSize() / 2 + 1);
  for(size_t i = 0; i < fLowPassCoeffs.size(); ++i) {
    float f = binWidth * i;
    fLowPassCoeffs[i] = 1.0 / sqrt(1.0 + pow(f/fCutoffFrequency, 8));
  }
  
  if( m_LogLevel >= 1 ) {
    cout << myname << "               LogLevel: " << m_LogLevel << endl;
    cout << myname << "    RemoveHighFrequency: " << fRemoveHighFrequency << endl;
    cout << myname << "         RemoveCoherent: " << fRemoveCoherent << endl;
    cout << myname << "        CutoffFrequency: " << fCutoffFrequency << endl;
    cout << myname << "      CoherentOffline16: " << fCoherentOffline16 << endl;
    cout << myname << "           CoherentDaq8: " << fCoherentDaq8  << endl;
    cout << myname << "          CoherentDaq16: " << fCoherentDaq16  << endl;
    cout << myname << "        CoherentFEMB128: " << fCoherentFEMB128  << endl;
    cout << myname << "CoherentOffline16Groups:  [";
    for(size_t i=0; i<fCoherentOffline16Groups.size(); i++) {
       cout<< fCoherentOffline16Groups.at(i) <<", ";	
    }
    cout << " ]" << endl;
    cout << myname << "CoherentDaq8Groups:  [";
    for(size_t i=0; i<fCoherentDaq8Groups.size(); i++) {
       cout<< fCoherentDaq8Groups.at(i) <<", ";	
    }
    cout << " ]" << endl;
    cout << myname << "CoherentDaq16Groups:  [";
    for(size_t i=0; i<fCoherentDaq16Groups.size(); i++) {
       cout<< fCoherentDaq16Groups.at(i) <<", ";	
    }
    cout << " ]" << endl;
    cout << myname << "CoherentFEMB128Groups:  [";
    for(size_t i=0; i<fCoherentFEMB128Groups.size(); i++) {
       cout<< fCoherentFEMB128Groups.at(i) <<", ";	
    }
    cout << " ]" << endl;
    cout << myname << "      UseBasicROIForCNR: " << fUseBasicROIForCNR  << endl;
    cout << myname << "      RoiStartThreshold: " << fRoiStartThreshold  << endl;
    cout << myname << "        RoiEndThreshold: " << fRoiEndThreshold  << endl;
    cout << myname << "              RoiPadLow: " << fRoiPadLow  << endl;
    cout << myname << "             RoiPadHigh: " << fRoiPadHigh  << endl;
  }  
}

//**********************************************************************
DataMap PdspNoiseRemoval::updateMap(AdcChannelDataMap& acds) const {
	const string myname = "PdspNoiseRemoval::updateMap: ";
	if( updateWithView() ) return viewMap(acds);
	DataMap ret(0);
  if (acds.size() == 0) {
    std::cout << myname << "WARNING: No channels found." << std::endl;
    ret.setStatus(1);
    return ret;
  }

	if(fRemoveHighFrequency) {
    removeHighFreq(acds);
  }
  if(fRemoveCoherent) {
    if(fCoherentOffline16) {
    	auto ch_groups = makeGroupsByOfflineChannels(16, fCoherentOffline16Groups);
      removeCoherent(ch_groups, acds);
    }                 
	  if(fCoherentFEMB128) {
    	auto ch_groups = makeGroupsByFEMBPlaneType(128, fCoherentFEMB128Groups);
      removeCoherent(ch_groups, acds);
    }
	  if(fCoherentDaq16) {
      auto ch_groups = makeGroupsByDAQChannels(16, fCoherentDaq16Groups);
      removeCoherent(ch_groups, acds);
    }
    if(fCoherentDaq8) {
      auto ch_groups = makeGroupsByDAQChannels(8, fCoherentDaq8Groups);
      removeCoherent(ch_groups, acds);
    }	
  }
  ret.setStatus(0);
  return ret;
}

//**********************************************************************
void PdspNoiseRemoval::removeHighFreq(AdcChannelDataMap& datamap) const
{
    auto const & chStatus = art::ServiceHandle< lariov::ChannelStatusService >()->GetProvider();
    for(auto & entry : datamap) {
      if(chStatus.IsPresent(entry.first) && !chStatus.IsNoisy(entry.first)) {
      	fftFltInPlace(entry.second.samples, fLowPassCoeffs);
      }
    }
}

//**********************************************************************
void PdspNoiseRemoval::fftFltInPlace(std::vector< float > & adc, const std::vector< float > & coeffs) const {
  std::vector< TComplex > ch_spectrum(fFFT->FFTSize() / 2 + 1);
  std::vector< float > ch_waveform(fFFT->FFTSize(), 0);
  size_t n_samples = adc.size();
  std::copy(adc.begin(), adc.end(), ch_waveform.begin());
  for(size_t s = n_samples; s < ch_waveform.size(); ++s) {
    ch_waveform[s] = ch_waveform[2*n_samples - s - 1];
  }
  fFFT->DoFFT(ch_waveform, ch_spectrum);
  for(size_t c = 0; c < coeffs.size(); ++c) {
      ch_spectrum[c] *= coeffs[c]; 
  }
  fFFT->DoInvFFT(ch_spectrum, ch_waveform);
  std::copy(ch_waveform.begin(), ch_waveform.begin()+n_samples, adc.begin());
}

//**********************************************************************
std::vector< float > PdspNoiseRemoval::fftFlt(const std::vector< float > & adc, const std::vector< float > & coeffs) const {
  std::vector< TComplex > ch_spectrum(fFFT->FFTSize() / 2 + 1);
  std::vector< float > ch_waveform(fFFT->FFTSize(), 0);
  size_t n_samples = adc.size();
  std::copy(adc.begin(), adc.end(), ch_waveform.begin());
  for(size_t s = n_samples; s < ch_waveform.size(); ++s) {
    ch_waveform[s] = ch_waveform[2*n_samples - s - 1];
  }
  fFFT->DoFFT(ch_waveform, ch_spectrum);
  for(size_t c = 0; c < coeffs.size(); ++c) {
    ch_spectrum[c] *= coeffs[c];
  }
  fFFT->DoInvFFT(ch_spectrum, ch_waveform);
  std::vector< float > flt_adc(n_samples);
  std::copy(ch_waveform.begin(), ch_waveform.begin()+n_samples, flt_adc.begin());
  return flt_adc;
}

//**********************************************************************
GroupChannelMap PdspNoiseRemoval::makeGroupsByOfflineChannels(size_t gsize, const std::vector< size_t > & gidx) const {
  GroupChannelMap groups;
  auto const & chStatus = art::ServiceHandle< lariov::ChannelStatusService >()->GetProvider();
  const unsigned int nchan = fGeometry->Nchannels();
  for(unsigned int ch = 0; ch < nchan; ++ch) {
    size_t g = ch / gsize;
    if(gidx.empty() || has(gidx, g)) {
      if(chStatus.IsPresent(ch) && !chStatus.IsBad(ch) && !chStatus.IsNoisy(ch)) { groups[g].push_back(ch); }
    }
  }
  return groups;
}

//**********************************************************************
GroupChannelMap PdspNoiseRemoval::makeGroupsByDAQChannels(size_t gsize, const std::vector< size_t > & gidx) const {
  GroupChannelMap groups;
  auto const & chStatus = art::ServiceHandle< lariov::ChannelStatusService >()->GetProvider();
  const unsigned int nchan = fGeometry->Nchannels();
  for(unsigned int ch = 0; ch < nchan; ++ch) {
    size_t g = getDAQChan(ch) / gsize;
    if(gidx.empty() || has(gidx, g)) {
      if(chStatus.IsPresent(ch) && !chStatus.IsBad(ch) && !chStatus.IsNoisy(ch)) {
      	groups[g].push_back(ch); 
      }
    }
  }
  return groups;
}

//**********************************************************************
GroupChannelMap PdspNoiseRemoval::makeGroupsByFEMBPlaneType(size_t gsize, const std::vector< size_t > & gidx) const {
  // Get channel map
  art::ServiceHandle<dune::PdspChannelMapService> channelMap;
  GroupChannelMap groups;
  auto const & chStatus = art::ServiceHandle< lariov::ChannelStatusService >()->GetProvider();
  const unsigned int nchan = fGeometry->Nchannels();
  for(unsigned int ch = 0; ch < nchan; ++ch) {
    //size_t g = ch / gsize;
    size_t g = getDAQChan(ch) / gsize;
    size_t plane = channelMap->PlaneFromOfflineChannel(ch);
    if(gidx.empty() || has(gidx, g) ) {
      if(chStatus.IsPresent(ch) && !chStatus.IsBad(ch) && !chStatus.IsNoisy(ch)) { 
      	switch (plane) {
          case 0: groups[3*g].push_back(ch); 
            break;
          case 1: groups[3*g+1].push_back(ch); 
            break;
          case 2: groups[3*g+2].push_back(ch); 
          	break;
        }
      }
    }
  }
  return groups;
}

//**********************************************************************
void PdspNoiseRemoval::removeCoherent(const GroupChannelMap & ch_groups, AdcChannelDataMap& datamap) const {
  if(datamap.empty()) return;
  //PdspNoiseRemoval *th = const_cast <PdspNoiseRemoval *> (this); 
  size_t n_samples = datamap.begin()->second.samples.size();
  std::vector<float> correction(n_samples);   
  std::vector<float> correctionFFT(n_samples/2);
  //std::vector<float> correctionFFT(n_samples / 2);
  for(const auto & entry : ch_groups) {   
    const auto & channels = entry.second;
    std::vector<float> correction;
    if(fMode == 1)  { // mean
      correction = getMeanCorrection(channels, datamap);
    }
    else if(fMode == 2)  { // median
      correction = getMedianCorrection(channels, datamap);
    }
    for(unsigned int ch : channels) {
      auto iacd = datamap.find(ch);
      if(iacd == datamap.end()) continue;
      AdcChannelData & acd = iacd->second;
      if(acd.samples.size() == 0) continue;
      for(size_t s = 0; s < acd.samples.size(); ++s) {
        acd.samples[s] -= correction[s];  
      }
    }
  }
}

//**********************************************************************
std::vector<float> PdspNoiseRemoval::getMeanCorrection( const std::vector<unsigned int> & channels, const AdcChannelDataMap & datamap) const {
  size_t n_samples = datamap.begin()->second.samples.size();
  std::vector<size_t> ch_averaged(n_samples, 0);
  std::vector<float> correction(n_samples, 0);
  for(unsigned int ch : channels) {    
    auto iacd = datamap.find(ch);
    if(iacd == datamap.end()) continue;
    const AdcChannelData & acd = iacd->second; 
    
    if(fUseBasicROIForCNR) {
      auto mask = roiMask(acd);  
      for(size_t s = 0; s < n_samples; ++s) {
        if (!mask[s]) { continue; }
        AdcFlag flag = acd.flags.size() ? acd.flags[s] : AdcGood;
        if(flag != AdcGood) { continue; }
        correction[s] += acd.samples[s];
        ch_averaged[s]++;
      }
    }
    else {
      for(size_t s = 0; s < n_samples; ++s) {
	      if(acd.signal[s]) { continue; }
        AdcFlag flag = acd.flags.size() ? acd.flags[s] : AdcGood;
        if(flag != AdcGood) { continue; }
        correction[s] += acd.samples[s];
        ch_averaged[s]++;
      }
    }
  }
  for(size_t s = 0; s < n_samples; ++s) {
    if(ch_averaged[s] > 0) { 
    	correction[s] /= ch_averaged[s]; 
    }
  }
  return correction;
}
std::vector<float> PdspNoiseRemoval::getMedianCorrection( const std::vector<unsigned int> & channels, const AdcChannelDataMap & datamap) const {
  size_t n_samples = datamap.begin()->second.samples.size();
  std::vector< std::vector<float> > samples(n_samples);
  for(unsigned int ch : channels) {
    auto iacd = datamap.find(ch);
    if(iacd == datamap.end()) continue;
    const AdcChannelData & acd = iacd->second;
    if(fUseBasicROIForCNR) {
      auto mask = roiMask(acd); 
      for(size_t s = 0; s < n_samples; ++s) {
        if (!mask[s]) { continue; }
        AdcFlag flag = acd.flags.size() ? acd.flags[s] : AdcGood;
        if(flag != AdcGood) { continue; }
        samples[s].push_back(acd.samples[s]);
      }
    }
    else {
      for(size_t s = 0; s < n_samples; ++s) {
	      if(acd.signal[s]) { continue; }
	      AdcFlag flag = acd.flags.size() ? acd.flags[s] : AdcGood;
        if(flag != AdcGood) { continue; }
        samples[s].push_back(acd.samples[s]);
      }
    }
  }
  std::vector<float> correction(n_samples);
  for(size_t s = 0; s < n_samples; ++s) {
    size_t n = samples[s].size();
    if(n < 2) { correction[s] = 0; continue; }
    std::sort(samples[s].begin(), samples[s].end());
    if((n % 2) == 0) { correction[s] = 0.5 * (samples[s][n/2] + samples[s][(n/2)-1]); }
    else { correction[s] = samples[s][(n-1)/2]; }
  }
  return correction;
}

//**********************************************************************
std::vector<bool> PdspNoiseRemoval::roiMask(const AdcChannelData & acd) const {
  std::vector<bool> mask(acd.samples.size(), true);
  auto acd_flt = fftFlt(acd.samples, fLowPassCoeffs);
  bool inroi = false;
  for (int i = 0; i < (int)acd_flt.size(); ++i) {
    auto sig = acd_flt[i];
    if (inroi){
      if (sig > fRoiEndThreshold || sig < -1.0*fRoiEndThreshold) { mask[i] = false; }
      else {
        for (int p = 0; p <= fRoiPadHigh; ++p) { if ((i + p) < (int)mask.size()) { mask[i + p] = false; } }
        inroi = false;
      }
    }
    else {
      if (sig > fRoiStartThreshold || sig < -1.0*fRoiStartThreshold) {
        for (int p = fRoiPadLow; p >= 0; --p) { if (i - p >= 0) { mask[i - p] = false; } }
        inroi = true;
      }
    }
  }
  return mask;
}

//**********************************************************************
size_t PdspNoiseRemoval::getDAQChan(size_t LAr_chan) {
	// Get channel map
  art::ServiceHandle<dune::PdspChannelMapService> channelMap;
  size_t apa = channelMap->APAFromOfflineChannel(LAr_chan);
  size_t wib = channelMap->WIBFromOfflineChannel(LAr_chan);
  size_t femb = channelMap->FEMBFromOfflineChannel(LAr_chan);
  size_t fembChan = channelMap->FEMBChannelFromOfflineChannel(LAr_chan);
  size_t DaqChan = 2560*apa + 512* wib + 128* femb + fembChan; //does not depend on RCE or FELIX
  return DaqChan;
}
//**********************************************************************
