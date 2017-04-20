// aclic_ArtServiceHelper.h

// David Adams
// May 2016
//
// Special version of ArtServiceHelper header to be used in Root, e.g.
//
// gSystem->AddLinkedLibs($DUNETPC_LIB/libdune_ArtSupport.so");
// gROOT->ProcessLine(".L $DUNETPC_INC/dune/ArtSupport/aclic_ArtServiceHelper.h+");

#ifndef aclic_ArtServiceHelper_H
#define aclic_ArtServiceHelper_H

namespace art {
  class ServiceToken { };
  class ServiceRegistry {
  public:
    class Operate {
    public:
      Operate(ServiceToken) { }
    };
  };
}

#define HideArtServices

#include "dune/ArtSupport/ArtServiceHelper.h"

#pragma link C++ class ArtServiceHelper;

#endif
