////////////////////////////////////////////////////////////////////////
// Class:       HelloAuxDet
// Plugin Type: analyzer (art v2_10_03)
// File:        HelloAuxDet_module.cc
// Brief:       Demonstration of how to access AuxDetGeos and AuxDetDigits 
//              in LArSoft.  Checks for presence of AuxDetGeos in geometry 
//              service for ProtoDUNE-SP.  Absence may be due to needing 
//              a dedicated ChannelMapAlg to register AuxDets with 
//              the Geometry service.  
//
// Generated at Wed May 16 09:02:34 2018 by Andrew Olivier (aolivier@ur.rochester.edu) using cetskelgen
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

namespace ex {
  class HelloAuxDet;
}


class ex::HelloAuxDet : public art::EDAnalyzer {
public:
  explicit HelloAuxDet(fhicl::ParameterSet const & p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  HelloAuxDet(HelloAuxDet const &) = delete;
  HelloAuxDet(HelloAuxDet &&) = delete;
  HelloAuxDet & operator = (HelloAuxDet const &) = delete;
  HelloAuxDet & operator = (HelloAuxDet &&) = delete;

  // Required functions.
  void analyze(art::Event const & e) override;

  // Selected optional functions.
  void beginJob() override;

private:

  // Declare member data here.
  std::vector<art::InputTag> fDigitLabels; //Search for AuxDetDigits produced by modules with these labels
  std::vector<art::InputTag> fSimLabels; //Search for AuxDetSimChannels produced by modules with these labels
};


ex::HelloAuxDet::HelloAuxDet(fhicl::ParameterSet const & p)
  :
  EDAnalyzer(p)  // ,
 // More initializers here.
{
  fDigitLabels = p.get<std::vector<art::InputTag>>("DigitLabels", {});
  fSimLabels = p.get<std::vector<art::InputTag>>("SimLabels", {});
}

void ex::HelloAuxDet::analyze(art::Event const & e)
{
  // Print out information about AuxDetSimChannels and AuxDetDigits.  
  art::Handle<std::vector<raw::AuxDetDigit>> digits; 
  if(fDigitLabels.empty())
  {
    e.getManyByType(digits);
  }
  else
  {
    e.getByLabel(fDigitLabels, digits); 
  }

  //TODO: Print AuxDetDigits
  for(const auto& digit: digits)
  {
  }

  art::Handle<std::vector<simb::AuxDetSimChannel>> channels;
  if(fSimLabels.empty())
  {
    e.getManyByType(channels);
  }
  else
  {
    e.getByLabel(fSimLabels, channels); 
  }

  for(const auto& channel: channels)
  { 
    //TODO: Print AuxDetSimChannels
  }
}

void ex::HelloAuxDet::beginJob()
{
  // Print all of the AuxDetGeos that the Geometry service knows about.  
  const auto geom = lar::providerFrom<geo::Geometry>();

  //TODO: Get AuxDetGeos and print their names.  

  //Tell the framework that I intend to consume AuxDetDigits and AuxDetSimChannels
  for(const auto& label: fDigitLabels) consumes<std::vector<raw::AuxDetDigit>>(label);
  for(const auto& label: fSimLabels) consumes<std::vector<simb::AuxDetSimChannel>>(label);
}

DEFINE_ART_MODULE(ex::HelloAuxDet)
