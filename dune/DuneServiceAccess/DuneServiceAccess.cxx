// DuneServiceAccess.cxx

#include "DuneServiceAccess.h"

//**********************************************************************

template<>
ChannelGroupService* ArtServicePointer<ChannelGroupService>() {
  return GenericArtServicePointer<ChannelGroupService>();
}

template<>
ChannelMappingService* ArtServicePointer<ChannelMappingService>() {
  return GenericArtServicePointer<ChannelMappingService>();
}

//**********************************************************************

template<>
SimChannelExtractService* ArtServicePointer<SimChannelExtractService>() {
  return GenericArtServicePointer<SimChannelExtractService>();
}

template<>
PedestalAdditionService* ArtServicePointer<PedestalAdditionService>() {
  return GenericArtServicePointer<PedestalAdditionService>();
}

template<>
AdcDistortionService* ArtServicePointer<AdcDistortionService>() {
  return GenericArtServicePointer<AdcDistortionService>();
}

template<>
ChannelNoiseService* ArtServicePointer<ChannelNoiseService>() {
  return GenericArtServicePointer<ChannelNoiseService>();
}

template<>
AdcCompressService* ArtServicePointer<AdcCompressService>() {
  return GenericArtServicePointer<AdcCompressService>();
}

template<>
AdcSuppressService* ArtServicePointer<AdcSuppressService>() {
  return GenericArtServicePointer<AdcSuppressService>();
}

template<>
MultiAdcSuppressService* ArtServicePointer<MultiAdcSuppressService>() {
  return GenericArtServicePointer<MultiAdcSuppressService>();
}

template<>
BaseTriggerService* ArtServicePointer<BaseTriggerService>() {
  return GenericArtServicePointer<BaseTriggerService>();
}

//**********************************************************************

template<>
RawDigitPrepService* ArtServicePointer<RawDigitPrepService>() {
  return GenericArtServicePointer<RawDigitPrepService>();
}

template<>
RawDigitExtractService* ArtServicePointer<RawDigitExtractService>() {
  return GenericArtServicePointer<RawDigitExtractService>();
}

template<>
AdcMitigationService* ArtServicePointer<AdcMitigationService>() {
  return GenericArtServicePointer<AdcMitigationService>();
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

template<>
AdcWireBuildingService* ArtServicePointer<AdcWireBuildingService>() {
  return GenericArtServicePointer<AdcWireBuildingService>();
}

//**********************************************************************

#include "dune/Utilities/SignalShapingServiceDUNE.h"

template<>
SignalShapingService* ArtServicePointer<SignalShapingService>() {
  return GenericArtServicePointer<util::SignalShapingServiceDUNE>();
}

//**********************************************************************

