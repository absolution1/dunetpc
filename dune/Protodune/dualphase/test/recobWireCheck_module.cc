////////////////////////////////////////////////////////////////////////
// Class:       recobWireCheck
// Plugin Type: analyzer (art v3_02_06)
// File:        recobWireCheck_module.cc
//
// Generated at Thu Oct  3 14:07:32 2019 by Vyacheslav Galymov using cetskelgen
// from cetlib version v3_07_02.
//
// Get some statistics on recob::Wire objects built by DataPrep
//
////////////////////////////////////////////////////////////////////////
#include <iostream>

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "art_root_io/TFileService.h"

#include "lardataobj/RawData/raw.h"
//#include "lardataobj/RawData/RawDigit.h"
#include "lardataobj/RecoBase/Wire.h"

#include <TTree.h>


#include <algorithm>

namespace pddp {
  class recobWireCheck;
}


class pddp::recobWireCheck : public art::EDAnalyzer {
public:
  explicit recobWireCheck(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  recobWireCheck(recobWireCheck const&) = delete;
  recobWireCheck(recobWireCheck&&) = delete;
  recobWireCheck& operator=(recobWireCheck const&) = delete;
  recobWireCheck& operator=(recobWireCheck&&) = delete;

  // Required functions.
  void analyze(art::Event const& e) override;

  // Selected optional functions.
  void beginJob() override;
  void endJob() override;

private:

  // Declare member data here.
  int fLogLevel;
  std::string fInputLabel;
  
  //
  TTree *fTree;
  
  float fAdcMax;
  float fAdcSum;
  unsigned fChanId;
  unsigned fTicks;
  int      fStartTick;
  int      fEndTick;
  //int      fMaxTick;  
};


pddp::recobWireCheck::recobWireCheck(fhicl::ParameterSet const& p)
  : EDAnalyzer{p},
  fLogLevel( p.get< int >("LogLevel") ),
  fInputLabel( p.get< std::string  >("InputLabel") )
{
  // Call appropriate consumes<>() for any products to be retrieved by this module.
  //fWireToken = consumes< std::vector<recob::Wire> >(fInputLabel);
}
  
//
void pddp::recobWireCheck::analyze(art::Event const& e)
{
  //auto const& wireHandle = e.getValidHandle(fWireToken);
  auto wireHandle = e.getHandle< std::vector<recob::Wire> >(fInputLabel);
  if( !wireHandle )
    {
      std::cerr<<"Product "<<fInputLabel<<" was not found\n";
    }

  for(size_t wireIter = 0; wireIter < wireHandle->size(); wireIter++)
    {
      art::Ptr<recob::Wire> wire(wireHandle, wireIter);
      fChanId = wire->Channel();
      const recob::Wire::RegionsOfInterest_t& signals = wire->SignalROI();

      if( fLogLevel >= 2 ) {
	std::cout<<"Channel "<<fChanId<<" has number of ROIs "<<signals.n_ranges()<<std::endl;
      }
      
      // loop over regions of interest      
      for (const auto& range : signals.get_ranges())
	{
	  fStartTick = range.begin_index();
	  fEndTick   = range.end_index();
	  fTicks     = range.size();
	  
	  fAdcMax    = -9999; //*(std::max_element( range.begin(), range.end() ));
	  fAdcSum    = 0; 
	  for( float adc: range.data() )
	    {
	      fAdcSum += adc;
	      if( adc > fAdcMax ) fAdcMax = adc;
	    }
	  
	  //
	  fTree->Fill();
	}
    }
}

void pddp::recobWireCheck::beginJob()
{
  art::ServiceHandle<art::TFileService> tfs;
  fTree = tfs->make<TTree>("recobWireTree","recobWire summary tree");
  fTree->Branch("fChanId", &fChanId, "fChanId/i");
  fTree->Branch("fTicks", &fTicks, "fTicks/i");
  fTree->Branch("fStartTick", &fStartTick, "fStartTick/I");
  fTree->Branch("fEndTick", &fEndTick, "fEndTick/I");
  fTree->Branch("fAdcMax", &fAdcMax, "fAdcMax/F");
  fTree->Branch("fAdcSum", &fAdcSum, "fAdcSum/F");
}

void pddp::recobWireCheck::endJob()
{
  // Implementation of optional member function here.
}

DEFINE_ART_MODULE(pddp::recobWireCheck)
