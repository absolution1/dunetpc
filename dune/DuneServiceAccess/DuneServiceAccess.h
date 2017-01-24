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

#include "dune/DuneInterface/ChannelGroupService.h"
#include "dune/DuneInterface/ChannelMappingService.h"

#include "dune/DuneInterface/SimChannelExtractService.h"
#include "dune/DuneInterface/PedestalAdditionService.h"
#include "dune/DuneInterface/AdcDistortionService.h"
#include "dune/DuneInterface/ChannelNoiseService.h"
#include "dune/DuneInterface/AdcSuppressService.h"
#include "dune/DuneInterface/MultiAdcSuppressService.h"
#include "dune/DuneInterface/AdcCompressService.h"
#include "dune/DuneInterface/BaseTriggerService.h"

#include "dune/DuneInterface/RawDigitPrepService.h"
#include "dune/DuneInterface/RawDigitExtractService.h"
#include "dune/DuneInterface/AdcMitigationService.h"
#include "dune/DuneInterface/PedestalEvaluationService.h"
#include "dune/DuneInterface/AdcSignalFindingService.h"
#include "dune/DuneInterface/AdcNoiseRemovalService.h"
#include "dune/DuneInterface/AdcChannelNoiseRemovalService.h"
#include "dune/DuneInterface/AdcDeconvolutionService.h"
#include "dune/DuneInterface/AdcChannelDataCopyService.h"
#include "dune/DuneInterface/AdcRoiBuildingService.h"
#include "dune/DuneInterface/AdcWireBuildingService.h"

#include "dune/DuneInterface/SignalShapingService.h"

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
