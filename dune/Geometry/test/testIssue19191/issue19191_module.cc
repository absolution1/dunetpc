////////////////////////////////////////////////////////////////////////
// Class:       issue19191
// Plugin Type: analyzer (art v2_11_03)
// File:        issue19191_module.cc
//
// Generated at Thu Dec  6 12:02:48 2018 by Saba Sehrish using cetskelgen
// from cetlib version v3_03_01.
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

#include "larcore/Geometry/Geometry.h"


#include "art/Framework/Services/Registry/ServiceHandle.h"

#include <cassert>
#include <iostream>

class issue19191;


class issue19191 : public art::EDAnalyzer {
public:
  explicit issue19191(fhicl::ParameterSet const & p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  issue19191(issue19191 const &) = delete;
  issue19191(issue19191 &&) = delete;
  issue19191 & operator = (issue19191 const &) = delete;
  issue19191 & operator = (issue19191 &&) = delete;

  // Required functions.
  void analyze(art::Event const & e) override;

private:

  // Declare member data here.

};


issue19191::issue19191(fhicl::ParameterSet const & p)
  :
  EDAnalyzer(p)  // ,
 // More initializers here.
{}

void issue19191::analyze(art::Event const & )
{
  art::ServiceHandle<geo::Geometry> geom;
  auto const *gcore = geom.get();
  std::vector<raw::ChannelID_t> channels;
  channels.reserve(gcore->Nchannels());
  channels = gcore->ChannelsInTPCs();
  assert(std::is_sorted(channels.begin(), channels.end()));
  assert(std::is_unique(channels.begin(), channels.end()));
  assert(channels.size() == gcore->Nchannels());
}

DEFINE_ART_MODULE(issue19191)
