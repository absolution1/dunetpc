// DuneDPhase3x1x1NoiseRemovalService.h
//
// Robert Sulej, Aug 2017
//
// Remove coherent noise from 3x1x1 data.
//

#ifndef DuneDPhase3x1x1NoiseRemovalService_H
#define DuneDPhase3x1x1NoiseRemovalService_H

#include "TFFTRealComplex.h"
#include "TFFTComplexReal.h"

#include "dune/DuneInterface/AdcNoiseRemovalService.h"
#include "dune/DuneInterface/AdcTypes.h"

namespace geo { class Geometry; }
namespace util { class LArFFT; }

using GroupChannelMap = std::unordered_map<unsigned int, std::vector<unsigned int> >;

class DuneDPhase3x1x1NoiseRemovalService : public AdcNoiseRemovalService {

public:

  DuneDPhase3x1x1NoiseRemovalService(fhicl::ParameterSet const& pset, art::ActivityRegistry&);
  ~DuneDPhase3x1x1NoiseRemovalService();

  int update(AdcChannelDataMap& datamap) const;

  std::ostream& print(std::ostream& out = std::cout, std::string prefix = "") const;

private:

  std::vector<float> getMeanCorrection(
    const std::vector<unsigned int> & channels,
    const AdcChannelDataMap & datamap) const;

  std::vector<float> getMedianCorrection(
    const std::vector<unsigned int> & channels,
    const AdcChannelDataMap & datamap) const;

  void fftFltInPlace(std::vector< float > & adc, const std::vector< float > & coeffs) const;
  std::vector< float > fftFlt(const std::vector< float > & adc, const std::vector< float > & coeffs) const;

  size_t fSize_2, fFreqSize_2;
  mutable TFFTRealComplex* fFFT_2;
  mutable TFFTComplexReal* fInvFFT_2;
  void doFFT_2(std::vector< float > & input, std::vector< TComplex > & output) const;
  void doInvFFT_2(std::vector< TComplex > & input, std::vector< float > & output) const;

  void removeCoherent(const GroupChannelMap & ch_groups, AdcChannelDataMap& datamap) const;
  void removeHighFreq(AdcChannelDataMap& datamap) const;
  void removeLowFreq(AdcChannelDataMap& datamap) const;

  std::vector<bool> roiMask(const AdcChannelData & adc) const;

  /// Make groups of channels using 3x1x1 DAQ numbering. Channels tagged as noisy are excluded at this stage.
  GroupChannelMap makeDaqGroups(size_t gsize, const std::vector< size_t > & gidx) const;
  /// Make groups of channels using LArSoft numbering. Channels tagged as noisy are excluded at this stage.
  GroupChannelMap makeGroups(size_t gsize, const std::vector< size_t > & gidx) const;

  bool has(const std::vector<size_t> & v, size_t idx) const
  {
    for (auto c : v) if (c == idx) return true;
    return false;
  }

  /// Get 3x1x1 DAQ channel number from the LArSoft's channel index.
  static size_t get311Chan(size_t LAr_chan);

  // Configuration parameters.
  bool fCoherent32, fCoherent16, fLowPassFlt, fFlattenLowFreq;
  std::vector< size_t > fCoherent32Groups;
  std::vector< size_t > fCoherent16Groups;
  std::vector< float > fFlattenCoeffs, fLowPassCoeffs;
  float fRoiStartThreshold;
  float fRoiEndThreshold;
  int fRoiPadLow;
  int fRoiPadHigh;
  int fMode;

  // Services.
  const geo::Geometry* fGeometry;
  mutable util::LArFFT* fFFT;
};

DECLARE_ART_SERVICE_INTERFACE_IMPL(DuneDPhase3x1x1NoiseRemovalService, AdcNoiseRemovalService, LEGACY)

#endif
