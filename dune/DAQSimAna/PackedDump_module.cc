////////////////////////////////////////////////////////////////////////
// Class:       PackedDump
// Plugin Type: analyzer (art v2_10_03)
// File:        PackedDump_module.cc
//
// Generated at Thu Jun 14 08:21:34 2018 by Philip Rodrigues using cetskelgen
// from cetlib version v3_02_00.
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "dune-raw-data/Services/ChannelMap/PdspChannelMapService.h"

class PackedDump;


class PackedDump : public art::EDAnalyzer {
public:
  explicit PackedDump(fhicl::ParameterSet const & p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  PackedDump(PackedDump const &) = delete;
  PackedDump(PackedDump &&) = delete;
  PackedDump & operator = (PackedDump const &) = delete;
  PackedDump & operator = (PackedDump &&) = delete;

  // Required functions.
  void analyze(art::Event const & e) override;

private:

  art::ServiceHandle<dune::PdspChannelMapService> m_channelMap;
};


PackedDump::PackedDump(fhicl::ParameterSet const & p)
  :
  EDAnalyzer(p)  // ,
{}

void PackedDump::analyze(art::Event const & e)
{
  for(unsigned int i=0; i<2; ++i){
    for(unsigned int j=0; j<24; ++j){
      // args are: unsigned int crate, unsigned int slot, unsigned int fiber, unsigned int fembchannel
      unsigned int offlineChan=m_channelMap->GetOfflineNumberFromDetectorElements(0, 0, i, j);
      std::cout << "Fiber " << i << " FEMB channel " << j << " maps to offline " << offlineChan << std::endl;
    }
  }
}

DEFINE_ART_MODULE(PackedDump)
