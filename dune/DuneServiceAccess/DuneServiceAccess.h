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
//   root> ArtServiceHelper& ash = ArtServiceHelper::instance("standard_reco_dune35tdata.fcl");
//   root> ChannelGroupService* pcgs = ArtServicePointer<ChannelGroupService>();
//   root> pcgs->print();

#ifndef DuneServiceAccess_H
#define DuneServiceAccess_H

#include "dune/ArtSupport/ArtServicePointer.h"

#include "dune/DuneInterface/AdcSuppressService.h"
#include "dune/DuneInterface/AdcDistortionService.h"
#include "dune/DuneInterface/AdcCompressService.h"

#include "dune/DuneInterface/AdcMitigationService.h"
#include "dune/DuneInterface/ChannelGroupService.h"
#include "dune/DuneInterface/PedestalEvaluationService.h"
#include "dune/DuneInterface/AdcSignalFindingService.h"
#include "dune/DuneInterface/AdcNoiseRemovalService.h"
#include "dune/DuneInterface/AdcChannelNoiseRemovalService.h"
#include "dune/DuneInterface/AdcDeconvolutionService.h"
#include "dune/DuneInterface/AdcChannelDataCopyService.h"
#include "dune/DuneInterface/AdcRoiBuildingService.h"

template<class T> T* ArtServicePointer();

template<>            AdcCompressService* ArtServicePointer<AdcCompressService>();
template<>          AdcDistortionService* ArtServicePointer<AdcDistortionService>();
template<>           ChannelGroupService* ArtServicePointer<ChannelGroupService>();

template<>          AdcMitigationService* ArtServicePointer<AdcMitigationService>();
template<>            AdcSuppressService* ArtServicePointer<AdcSuppressService>();
template<>     PedestalEvaluationService* ArtServicePointer<PedestalEvaluationService>();
template<>       AdcSignalFindingService* ArtServicePointer<AdcSignalFindingService>();
template<>        AdcNoiseRemovalService* ArtServicePointer<AdcNoiseRemovalService>();
template<> AdcChannelNoiseRemovalService* ArtServicePointer<AdcChannelNoiseRemovalService>();
template<>       AdcDeconvolutionService* ArtServicePointer<AdcDeconvolutionService>();
template<>     AdcChannelDataCopyService* ArtServicePointer<AdcChannelDataCopyService>();
template<>         AdcRoiBuildingService* ArtServicePointer<AdcRoiBuildingService>();

#endif
