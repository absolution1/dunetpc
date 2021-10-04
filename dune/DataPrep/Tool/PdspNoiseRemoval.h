// PdspNoiseRemoval.h
// Jingbo Wang (UC Davis)
// Email: jiowang@ucdavis.edu
// January 2019
//
// Tool to remove noise
//
// At present, high frequency noise removal and coherent noise removal are supported. More noise removal functions will be added. 
// The high frequency noise removal is done first, followed by the coherent noise removal. Both are configurable. 
// The high frequency noise removal is done through a Butterworth low-pass filter. 
// The coherent noise is caculated as the mean or median waveform of a group of channels. 
// ProtoDUNE channel map service is used to make the online channel groups. 
// Different grouping methods are available. By default, every 128 online channels associated to the same FEMB is grouped together. 
// A simple threshold ROI finder is used to find the ROIs, and signals within the ROI are protected from the coherent noise calculation. 
//
// Configuration:
// Defined in dune/DataPrep/fcl/protodune_dataprep_tools.fcl
//  RemoveHighFrequency: false        # specify if the high frequency noise removal is enabled 
//  RemoveCoherent: true              # specify if the coherent noise removal is enabled. 
//                                      If RemoveHighFrequency is true, coherent noise removal is done afterward. 
//  CutoffFrequency: 300              # cutoff frequency in kHz for Butterworth low-pass filter
//  CorrMode: "median"                # mean or median waveform for coherent noise determination
//  CoherentOffline16: false          # remove coherent noise for groups consisting of 16 consecutive offline channels
//  CoherentOffline16Groups: []       # remove coherent noise for all groups if list is empty
//  CoherentDaq8: false               # remove coherent noise for groups consisting of 8 online DAQ channels on the same ASIC
//  CoherentDaq8Groups: []            # remove coherent noise for all groups if list is empty
//  CoherentDaq16: false              # remove coherent noise for groups consisting of 16 online DAQ channels on the same chip
//  CoherentDaq16Groups: []           # remove coherent noise for all groups if list is empty
//  CoherentFEMB128: true             # remove coherent noise for groups consisting of 128 online DAQ channels on the same FEMB
//  CoherentFEMB128Groups: []         # remove coherent noise for all groups if list is empty
//  UseBasicROIForCNR: true           # use simple threshold ROI finder
//  RoiStartThreshold: 20             # threshold on the leading edge
//  RoiEndThreshold: 20               # threshold on the trailing edge
//  RoiPadLow: 8                      # low bin extension  
//  RoiPadHigh: 20                    # high bin extension

// This tool was tested with ProtoDUNE-SP data Run4696. 


#ifndef PdspNoiseRemoval_H
#define PdspNoiseRemoval_H
#include "art/Utilities/ToolMacros.h"
#include "fhiclcpp/ParameterSet.h"
#include "dune/DuneInterface/Tool/TpcDataTool.h"
#include "TFFTRealComplex.h"
#include "TFFTComplexReal.h"
#include "dune-raw-data/Services/ChannelMap/PdspChannelMapService.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
namespace geo { class Geometry; }
namespace util { class LArFFT; }
using GroupChannelMap = std::unordered_map<unsigned int, std::vector<unsigned int> >;
	
class PdspNoiseRemoval : TpcDataTool {
	
public:
  PdspNoiseRemoval(fhicl::ParameterSet const& ps);
  ~PdspNoiseRemoval() override =default;
  DataMap updateMap(AdcChannelDataMap& acds) const override;
  
private:
  	
  // Configuration data.
  int  m_LogLevel;
  bool fRemoveHighFrequency, fRemoveCoherent;
  bool fCoherentOffline16, fCoherentDaq8, fCoherentDaq16, fCoherentFEMB128;
  std::vector< size_t > fCoherentOffline16Groups;
  std::vector< size_t > fCoherentDaq8Groups;
  std::vector< size_t > fCoherentDaq16Groups;
  std::vector< size_t > fCoherentFEMB128Groups;
  std::vector< float > fLowPassCoeffs;
  int fMode;
  bool fUseBasicROIForCNR;
  float fRoiStartThreshold;
  float fRoiEndThreshold;
  int fRoiPadLow;
  int fRoiPadHigh;
  float fCutoffFrequency;
  
  // Services.
  const geo::Geometry* fGeometry;
  mutable util::LArFFT* fFFT;
  	
  // Functions
  // Calculate mean wavefrom 
  std::vector<float> getMeanCorrection( const std::vector<unsigned int> & channels, const AdcChannelDataMap & datamap) const;
  // Calculate median wavefrom 
  std::vector<float> getMedianCorrection( const std::vector<unsigned int> & channels, const AdcChannelDataMap & datamap) const;
  // Remove coherent noise
  void removeCoherent(const GroupChannelMap & ch_groups, AdcChannelDataMap& datamap) const;
  // Remove high frequency noise
  void removeHighFreq(AdcChannelDataMap& datamap) const;
  // Do FFT and remove high frequency components for the input wavefrom
  void fftFltInPlace(std::vector< float > & adc, const std::vector< float > & coeffs) const;
  // Do FFT and remove high frequency components. Return the filtered wavefrom.
  std::vector< float > fftFlt(const std::vector< float > & adc, const std::vector< float > & coeffs) const;
  // Define ROI mask using a simple threshold ROI finder
  std::vector<bool> roiMask(const AdcChannelData & adc) const;
  // Make groups of channels using LArSoft numbering. Channels tagged as noisy are excluded at this stage.
  GroupChannelMap makeGroupsByOfflineChannels(size_t gsize, const std::vector< size_t > & gidx) const;
  // Make groups of channels by DAQ channel number.
  GroupChannelMap makeGroupsByDAQChannels(size_t gsize, const std::vector< size_t > & gidx) const;
  // Make groups of channels by FEMB and plane type.
  GroupChannelMap makeGroupsByFEMBPlaneType(size_t gsize, const std::vector< size_t > & gidx) const;
  bool has(const std::vector<size_t> & v, size_t idx) const {
    for (auto c : v) if (c == idx) return true;
    return false;
  }
  // Get DAQ channel number from the LArSoft's channel index.
  static size_t getDAQChan(size_t LAr_chan);
};
#endif
