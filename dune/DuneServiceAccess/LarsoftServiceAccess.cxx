// LarsoftServiceAccess.cxx

#include "LarsoftServiceAccess.h"
#include "dune/ArtSupport/ArtServiceProvider.h"

#include "larevt/CalibrationDBI/Interface/ChannelStatusService.h"

//**********************************************************************

template<>
const ChannelStatusProvider* ArtServiceProvider<lariov::ChannelStatusProvider>() {
  return DefaultArtServiceProvider<lariov::ChannelStatusService>();
}

//**********************************************************************
