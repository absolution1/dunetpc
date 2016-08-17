////////////////////////////////////////////////////////////////////////
// Class:       TestHardwareMapper
// Module Type: analyzer
// File:        TestHardwareMapper_module.cc
//
// Jonathan Davies j.p.davies@sussex.ac.uk
// August 2016
//
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// DUNETPC specific includes
#include "dune/DAQTriggerSim/Service/HardwareMapperService.h"


class TestHardwareMapper;

class TestHardwareMapper : public art::EDAnalyzer {

public:

  explicit TestHardwareMapper(fhicl::ParameterSet const & pset);

  // Plugins should not be copied or assigned.
  TestHardwareMapper(TestHardwareMapper const &) = delete;
  TestHardwareMapper(TestHardwareMapper &&) = delete;
  TestHardwareMapper & operator = (TestHardwareMapper const &) = delete;
  TestHardwareMapper & operator = (TestHardwareMapper &&) = delete;

  // The main guts...
  void analyze(art::Event const& evt) override;

  void beginJob() override;

  void reconfigure(fhicl::ParameterSet const & pset);

private:
  art::ServiceHandle<HardwareMapperService> mapperService;

};

//......................................................
TestHardwareMapper::TestHardwareMapper(fhicl::ParameterSet const & pset)
  : EDAnalyzer(pset)
{
  INFO_FUNCTION << std::endl;
  this->reconfigure(pset);
}

//......................................................
void TestHardwareMapper::reconfigure(fhicl::ParameterSet const & pset){
  INFO_FUNCTION << std::endl;
}

//......................................................
void TestHardwareMapper::analyze(art::Event const& evt){
  INFO_FUNCTION << std::endl;
}

//......................................................
void TestHardwareMapper::beginJob(){
  INFO_FUNCTION << std::endl;
  mapperService->printGeometryInfo();
}

//......................................................

//......................................................

//......................................................

//......................................................

DEFINE_ART_MODULE(TestHardwareMapper)
