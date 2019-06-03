/*
   Raw data decoder from uTCA crate AMCs

 */
// art includes
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art_root_io/TFileService.h"

// artdaq and dune-raw-data includes
#include "dune-raw-data/Overlays/RceFragment.hh"
#include "artdaq-core/Data/Fragment.hh"
#include "artdaq-core/Data/ContainerFragment.hh"
#include "dune-raw-data/Overlays/FragmentType.hh"
#include "dune-raw-data/Services/ChannelMap/PdspChannelMapService.h"

// larsoft includes
#include "lardataobj/RawData/RawDigit.h"

// ROOT includes
#include "TH1.h"
#include "TH2F.h"
#include "TProfile.h"

// C++ Includes
#include <memory>
#include <iostream>

namespace dune {
  class AmcRawDecoder;
}
