////////////////////////////////////////////////////////////////////////
// Class:       PDDPTPCRawDecoder
// Plugin Type: producer (art v3_02_05)
// File:        PDDPTPCRawDecoder_module.cc
//
// Generated at Tue Jun  4 11:04:58 2019 by Vyacheslav Galymov using cetskelgen
// from cetlib version v3_07_02.
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include <memory>

class PDDPTPCRawDecoder;


class PDDPTPCRawDecoder : public art::EDProducer {
public:
  explicit PDDPTPCRawDecoder(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  PDDPTPCRawDecoder(PDDPTPCRawDecoder const&) = delete;
  PDDPTPCRawDecoder(PDDPTPCRawDecoder&&) = delete;
  PDDPTPCRawDecoder& operator=(PDDPTPCRawDecoder const&) = delete;
  PDDPTPCRawDecoder& operator=(PDDPTPCRawDecoder&&) = delete;

  // Required functions.
  void produce(art::Event& e) override;

private:

  produce;
  // Declare member data here.

};


PDDPTPCRawDecoder::PDDPTPCRawDecoder(fhicl::ParameterSet const& p)
  : EDProducer{p}  // ,
  // More initializers here.
{
  // Call appropriate produces<>() functions here.
  // Call appropriate consumes<>() for any products to be retrieved by this module.
}

void PDDPTPCRawDecoder::produce(art::Event& e)
{
  // Implementation of required member function here.
}

DEFINE_ART_MODULE(PDDPTPCRawDecoder)
