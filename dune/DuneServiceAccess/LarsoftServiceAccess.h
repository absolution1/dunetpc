// LarsoftServiceAccess.h

// David Adams
// January 2017
//
// Add ArtServicePointer specializations for LArSoft services so they
// can be accessed from the Root command line.
//
// The providers listed here may be accessed via provider type instead
// of service type.
//
// Usage example:
//   root> gSystem->SetBuildDir(".aclic");
//   root> gSystem->AddLinkedLibs(gSystem->ExpandPathName("$DUNETPC_LIB/libdune_ArtSupport.so"));
//   root> gSystem->AddLinkedLibs(gSystem->ExpandPathName("$DUNETPC_LIB/libdune_LarsoftServiceAccess.so"));
//   root> .L $DUNETPC_INC/dune/ArtSupport/ArtServiceHelper.h+
//   root> .L $DUNETPC_INC/dune/LarsoftServiceAccess/LarsoftServiceAccess.h+
//   root> ArtServiceHelper& ash = ArtServiceHelper::load("standard_reco_dune35tdata.fcl");
//   root> ChannelGroupService* pcgs = ArtServicePointer<ChannelGroupService>();
//   root> pcgs->print();

#ifndef LarsoftServiceAccess_H
#define LarsoftServiceAccess_H

#include "larcorealg/Geometry/GeometryCore.h"
#include "larevt/CalibrationDBI/Interface/ChannelStatusProvider.h"

using geo::GeometryCore;
using lariov::ChannelStatusProvider;

template<class P> const P* ArtServiceProvider();

template<>   const GeometryCore*          ArtServiceProvider<GeometryCore>();
template<>   const ChannelStatusProvider* ArtServiceProvider<ChannelStatusProvider>();

#endif
