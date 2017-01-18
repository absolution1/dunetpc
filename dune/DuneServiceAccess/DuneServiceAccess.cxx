// DuneServiceAccess.cxx

#include "DuneServiceAccess.h"

//**********************************************************************

template<>
AdcSuppressService* ArtServicePointer<AdcSuppressService>() {
  return GenericArtServicePointer<AdcSuppressService>();
}

template<>
AdcDistortionService* ArtServicePointer<AdcDistortionService>() {
  return GenericArtServicePointer<AdcDistortionService>();
}

template<>
AdcCompressService* ArtServicePointer<AdcCompressService>() {
  return GenericArtServicePointer<AdcCompressService>();
}

//**********************************************************************

template<>
AdcMitigationService* ArtServicePointer<AdcMitigationService>() {
  return GenericArtServicePointer<AdcMitigationService>();
}

template<>
ChannelGroupService* ArtServicePointer<ChannelGroupService>() {
  return GenericArtServicePointer<ChannelGroupService>();
}

template<>
PedestalEvaluationService* ArtServicePointer<PedestalEvaluationService>() {
  return GenericArtServicePointer<PedestalEvaluationService>();
}

template<>
AdcSignalFindingService* ArtServicePointer<AdcSignalFindingService>() {
  return GenericArtServicePointer<AdcSignalFindingService>();
}

template<>
AdcNoiseRemovalService* ArtServicePointer<AdcNoiseRemovalService>() {
  return GenericArtServicePointer<AdcNoiseRemovalService>();
}

template<>
AdcChannelNoiseRemovalService* ArtServicePointer<AdcChannelNoiseRemovalService>() {
  return GenericArtServicePointer<AdcChannelNoiseRemovalService>();
}

template<>
AdcDeconvolutionService* ArtServicePointer<AdcDeconvolutionService>() {
  return GenericArtServicePointer<AdcDeconvolutionService>();
}

template<>
AdcChannelDataCopyService* ArtServicePointer<AdcChannelDataCopyService>() {
  return GenericArtServicePointer<AdcChannelDataCopyService>();
}

template<>
AdcRoiBuildingService* ArtServicePointer<AdcRoiBuildingService>() {
  return GenericArtServicePointer<AdcRoiBuildingService>();
}

//**********************************************************************

