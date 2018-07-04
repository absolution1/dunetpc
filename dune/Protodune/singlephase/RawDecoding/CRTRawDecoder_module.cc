////////////////////////////////////////////////////////////////////////
// Class:       CRTRawDecoder
// Plugin Type: producer (art v2_10_03)
// File:        CRTRawDecoder_module.cc
//
// Generated at Wed Jul  4 08:14:43 2018 by Andrew Olivier using cetskelgen
// from cetlib version v3_02_00.
////////////////////////////////////////////////////////////////////////

//Framework includes
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

//dune-raw-data includes
#include "dune-raw-data/Overlays/CRTFragment.hh"

//dunetpc includes
#include "dunetpc/dune/Protodune/singlephase/CRT/Trigger.h"

//c++ includes
#include <memory>

namespace CRT 
{
  class CRTRawDecoder;
}

//A CRTRawDecoder takess artdaq::Fragments made from Cosmic Ray Tagger input 
//and produces a CRT::Trigger for each time a CRT module triggered in this Event.  

class CRT::CRTRawDecoder : public art::EDProducer 
{
  public:
    explicit CRTRawDecoder(fhicl::ParameterSet const & p);
    // The compiler-generated destructor is fine for non-base
    // classes without bare pointers or other resource use.
  
    // Plugins should not be copied or assigned.
    CRTRawDecoder(CRTRawDecoder const &) = delete;
    CRTRawDecoder(CRTRawDecoder &&) = delete;
    CRTRawDecoder & operator = (CRTRawDecoder const &) = delete;
    CRTRawDecoder & operator = (CRTRawDecoder &&) = delete;
  
    // Required functions.
    void produce(art::Event & e) override;
  
    // Selected optional functions.  I might get services that I expect to change here.  So far, seems like only service I might need is timing.
    //void beginJob() override;
    //void beginRun(art::Run & r) override;
    //void beginSubRun(art::SubRun & sr) override;
  
  private:
  
    // Declare member data here.
    art::InputTag fFragLabel; //Label of the module that produced artdaq::Fragments 
                              //from CRT data that I will turn into CRT::Triggers.
  
};


CRT::CRTRawDecoder::CRTRawDecoder(fhicl::ParameterSet const & p): fFragLabel(p.get<art::InputTag>("CRTFragLabel"))
{
  // Call appropriate produces<>() functions here.
  produces<CRT::Trigger>();
}

//Read artdaq::Fragments produced by fFragLabel, and use CRT::Fragment to convert them to CRT::Triggers.  
void CRT::CRTRawDecoder::produce(art::Event & e)
{
  const auto& fragHandle = e.getValidHandle<std::vector<artdaq::Fragment>>(fFragLabel);
  auto triggers = std::make_unique<std::vector<CRT::Trigger>>();
  
  for(const auto& artFrag: *fragHandle)
  {
    CRT::Fragment frag(artFrag);
    
    std::vector<CRT::Hit> hits;
    for(size_t hitNum = 0; hitNum < frag.num_hits(); ++hitNum)
    {
      const auto hit = *(frag.hit(hitNum));
      hits.emplace_back(hit->channel, hit->adc);
    }

    triggers->emplace_back(frag.module_num(), frag.unixtime() << 32 | frag.fifty_mhz_time(), std::move(hits));
  } 

  e.put(std::move(triggers));
}

/*void CRT::CRTRawDecoder::beginJob()
{
  // Implementation of optional member function here.
}

void CRT::CRTRawDecoder::beginRun(art::Run & r)
{
  // Implementation of optional member function here.
}

void CRT::CRTRawDecoder::beginSubRun(art::SubRun & sr)
{
  // Implementation of optional member function here.
}*/

DEFINE_ART_MODULE(CRT::CRTRawDecoder)
