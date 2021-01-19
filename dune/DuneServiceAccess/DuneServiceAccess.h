// DuneServiceAccess.h

// David Adams
// January 2017
//
// Add ArtServicePointer specializations for DUNE services so they
// can be accessed from the Root command line.
//
// Usage example:
//   root> gSystem->SetBuildDir(".aclic");
//   root> gSystem->AddLinkedLibs(gSystem->ExpandPathName("$DUNETPC_LIB/libdune_ArtSupport.so"));
//   root> gSystem->AddLinkedLibs(gSystem->ExpandPathName("$DUNETPC_LIB/libdune_DuneServiceAccess.so"));
//   root> .L $DUNETPC_INC/dune/ArtSupport/ArtServiceHelper.h+
//   root> .L $DUNETPC_INC/dune/DuneServiceAccess/DuneServiceAccess.h+
//   root> ArtServiceHelper& ash = ArtServiceHelper::load("standard_reco_dune35tdata.fcl");
//   root> ChannelGroupService* pcgs = ArtServicePointer<ChannelGroupService>();
//   root> pcgs->print();

#ifndef DuneServiceAccess_H
#define DuneServiceAccess_H

#include "dune/ArtSupport/ArtServicePointer.h"

#include "dune/DuneInterface/Service/ChannelGroupService.h"
#include "dune/DuneInterface/Service/ChannelMappingService.h"

#include "dune/DuneInterface/Service/SimChannelExtractService.h"
#include "dune/DuneInterface/Service/PedestalAdditionService.h"
#include "dune/DuneInterface/Service/AdcDistortionService.h"
#include "dune/DuneInterface/Service/ChannelNoiseService.h"
#include "dune/DuneInterface/Service/AdcSuppressService.h"
#include "dune/DuneInterface/Service/MultiAdcSuppressService.h"
#include "dune/DuneInterface/Service/AdcCompressService.h"
#include "dune/DuneInterface/Service/BaseTriggerService.h"

#include "dune/DuneInterface/Service/RawDigitPrepService.h"
#include "dune/DuneInterface/Service/RawDigitExtractService.h"
#include "dune/DuneInterface/Service/AdcMitigationService.h"
#include "dune/DuneInterface/Service/PedestalEvaluationService.h"
#include "dune/DuneInterface/Service/AdcSignalFindingService.h"
#include "dune/DuneInterface/Service/AdcNoiseRemovalService.h"
#include "dune/DuneInterface/Service/AdcChannelNoiseRemovalService.h"
#include "dune/DuneInterface/Service/AdcDeconvolutionService.h"
#include "dune/DuneInterface/Service/AdcChannelDataCopyService.h"
#include "dune/DuneInterface/Service/AdcRoiBuildingService.h"
#include "dune/DuneInterface/Service/AdcWireBuildingService.h"

#include "dune/DuneInterface/Service/SignalShapingService.h"

template<class T> T* ArtServicePointer();

template<>           ChannelGroupService* ArtServicePointer<ChannelGroupService>();
template<>         ChannelMappingService* ArtServicePointer<ChannelMappingService>();

template<>      SimChannelExtractService* ArtServicePointer<SimChannelExtractService>();
template<>       PedestalAdditionService* ArtServicePointer<PedestalAdditionService>();
template<>          AdcDistortionService* ArtServicePointer<AdcDistortionService>();
template<>           ChannelNoiseService* ArtServicePointer<ChannelNoiseService>();
template<>            AdcSuppressService* ArtServicePointer<AdcSuppressService>();
template<>            AdcCompressService* ArtServicePointer<AdcCompressService>();
template<>            BaseTriggerService* ArtServicePointer<BaseTriggerService>();

template<>           RawDigitPrepService* ArtServicePointer<RawDigitPrepService>();
template<>        RawDigitExtractService* ArtServicePointer<RawDigitExtractService>();
template<>          AdcMitigationService* ArtServicePointer<AdcMitigationService>();
template<>     PedestalEvaluationService* ArtServicePointer<PedestalEvaluationService>();
template<>       AdcSignalFindingService* ArtServicePointer<AdcSignalFindingService>();
template<>        AdcNoiseRemovalService* ArtServicePointer<AdcNoiseRemovalService>();
template<> AdcChannelNoiseRemovalService* ArtServicePointer<AdcChannelNoiseRemovalService>();
template<>       AdcDeconvolutionService* ArtServicePointer<AdcDeconvolutionService>();
template<>     AdcChannelDataCopyService* ArtServicePointer<AdcChannelDataCopyService>();
template<>         AdcRoiBuildingService* ArtServicePointer<AdcRoiBuildingService>();
template<>        AdcWireBuildingService* ArtServicePointer<AdcWireBuildingService>();

template<>          SignalShapingService* ArtServicePointer<SignalShapingService>();

#endif
