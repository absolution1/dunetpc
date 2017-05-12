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

//We need a hardware mapper service, so we should include the header
#include "dune/DAQTriggerSim/Service/HardwareMapperService.h"
#include "dune/DAQTriggerSim/TriggerDataProducts/HardwareElements.h"

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
  unsigned int fLogLevel;
  art::ServiceHandle<HardwareMapperService> mapperService;
  std::string fModuleName;
};

//......................................................
TestHardwareMapper::TestHardwareMapper(fhicl::ParameterSet const & pset)
  : EDAnalyzer(pset), fLogLevel(1), fModuleName("TestHardwareMapper")
{
  pset.get_if_present<unsigned int>("LogLevel", fLogLevel);
  std::string func_name("TestHardwareMapper");
  if(fLogLevel>1) mf::LogInfo(fModuleName) << " in Function: " << func_name;
  this->reconfigure(pset);
}

//......................................................
void TestHardwareMapper::reconfigure(fhicl::ParameterSet const & pset){
  std::string func_name("TestHardwareMapper");
  if(fLogLevel>1) mf::LogInfo(fModuleName) << " in Function: " << func_name;
}

//......................................................
void TestHardwareMapper::analyze(art::Event const& evt){
  std::string func_name("TestHardwareMapper");
  if(fLogLevel>1) mf::LogInfo(fModuleName) << " in Function: " << func_name;

  //jpd -- this is how to loop through TPCs, asking the hardware mapper for the channels associated with each one
  unsigned int num_tpcs = mapperService->getNTPCs();
  {
    mf::LogInfo loginfo(fModuleName);
    for(unsigned int tpc_id = 0; tpc_id < num_tpcs; tpc_id++){
      //jpd -- now get the vector of channels
      std::vector<raw::ChannelID_t> tpc_channel_vector = mapperService->getTPCChannels(tpc_id);
      loginfo << "Got TPC number: " << tpc_id 
              << " - it has " << tpc_channel_vector.size() << " channels" 
              << " - first channel is " << tpc_channel_vector.at(0)
              << "\n";
    }//loop over TPCs
  }//annonymous to wrap loginfo
  
  //jpd -- this is how to loop through APAs, asking the hardware mapper for the channels associated with each one
  unsigned int num_apas = mapperService->getNAPAs();
  {
    mf::LogInfo loginfo(fModuleName);
    for(unsigned int apa_id = 0; apa_id < num_apas; apa_id++){
      //jpd -- now get the vector of channels
      std::vector<raw::ChannelID_t> apa_channel_vector = mapperService->getAPAChannels(apa_id);
      loginfo << "Got APA number: " << apa_id 
              << " - it has " << apa_channel_vector.size() << " channels" 
              << " - first channel is " << apa_channel_vector.at(0)
              << "\n";
    }//loop over APAs
  }//annonymous to wrap loginfo
}

//......................................................
void TestHardwareMapper::beginJob(){
  std::string func_name("TestHardwareMapper");
  if(fLogLevel>1) mf::LogInfo(fModuleName) << " in Function: " << func_name;
  mapperService->printGeometryInfo();
}

//......................................................

//......................................................

//......................................................

//......................................................

DEFINE_ART_MODULE(TestHardwareMapper)
