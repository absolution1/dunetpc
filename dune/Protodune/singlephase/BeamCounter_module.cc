////////////////////////////////////////////////////////////////////////
// Class:       BeamCounter
// Plugin Type: analyzer (art v2_08_03)
// File:        BeamCounter_module.cc
//
// Module to count the number of Triggers 
// received from the beamline
// and which cerenkov detectors were activated
//
// Jacob Calcutt Oct 9, 2018
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
#include "IFBeam_service.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "lardataobj/RecoBase/TrackingTypes.h"
#include "lardataobj/RecoBase/TrackTrajectory.h"
#include "lardataobj/RecoBase/Track.h"
#include "dune/Protodune/singlephase/CTB/data/pdspctb.h"
#include "lardataobj/RawData/RDTimeStamp.h"
#include "dune/DuneObj/ProtoDUNETimeStamp.h"
#include <bitset>
#include <iomanip>
#include <utility>
#include <algorithm>
#include <limits>

#include "TTree.h"
#include "TH2F.h"
#include "TVectorD.h"
#include "TPolyMarker.h"

namespace proto {
  class BeamCounter;
}



class proto::BeamCounter : public art::EDAnalyzer {
public:
  explicit BeamCounter(fhicl::ParameterSet const & p);
  //  virtual ~BeamCounter();

  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  BeamCounter(BeamCounter const &) = delete;
  BeamCounter(BeamCounter &&) = delete;
  BeamCounter & operator = (BeamCounter const &) = delete;
  BeamCounter & operator = (BeamCounter &&) = delete;

  // Required functions.
  void analyze(art::Event const & e) override;

  // Selected optional functions.
  void beginJob() override;
  
  void GetCTBInfo(art::Event const &);
  
private:
  TTree * fOutTree;

  int HLTWord;
  double HLTTS;
  int BITrigger;
  int C1;
  int C2;

  int eventNum;
  double eventTime;
  int runNum;
  int subRunNum;

};

proto::BeamCounter::BeamCounter(fhicl::ParameterSet const & p) : EDAnalyzer(p)
{}

//Gets the Timing and CTB raw-decoded info.
//Finds the triggers, and looks for a valid trigger word
//(i.e. coming from beam)
//
//Returns the timestamp of the high level trigger.
void proto::BeamCounter::GetCTBInfo(art::Event const & e){

  std::cout << std::endl;
  std::cout << "Getting Raw Decoder Info" << std::endl;

  auto const CTBHandle = e.getValidHandle< std::vector< raw::ctb::pdspctb > >("ctbrawdecoder:daq");
  std::cout << "CTB valid? " << CTBHandle.isValid() << std::endl;
  if(CTBHandle.isValid()){
    auto const & CTB = (*CTBHandle)[0];

    bool noHLT = true;

    std::cout << "NTriggers: " << CTB.GetNTriggers() << std::endl;
    for (size_t nTrig = 0; nTrig < CTB.GetNTriggers(); ++nTrig){

      raw::ctb::Trigger ctbTrig = CTB.GetTrigger(nTrig);      
      uint32_t  theType  = ctbTrig.word_type;
      ULong64_t theWord  = ctbTrig.trigger_word;
      ULong64_t theTS    = ctbTrig.timestamp;

      if (theType == 2){
        std::cout << "Found the High Level Trigger" << std::endl;
        
        HLTWord = theWord;
        HLTTS = 2.e-8*theTS;
        
        //The High Level Trigger consists of 8 bits
        //HLT7 -> HLT0
        std::bitset<8> theHLT(theWord);
        std::cout << "High Level Trigger: " << theHLT << std::endl;
        
        //HLT5 corresponds to excluding Low Level Triggers
        //from the Beamline
        //So return 0, we'll skip this event
        if (theHLT[5]) {         
          std::cout << "HLT 5 activated. Excluding beam events. Skipping this event" << std::endl;
          break;
        }
        else if (theHLT[0] && (theHLT.count() == 1)) {
          std::cout << "Only HLT 0 activated. This is just a random trigger. No beamline info was activated." << std::endl
                    << "Skipping this Event." << std::endl;
          break;
        }
        else{
          noHLT = false;
          std::cout << "Found valid beam event." << std::endl;
          break;
        }
      }       
    }

    if(noHLT){
      //This should never happen but...
      std::cout << "Error! No High Level Trigger Found!" << std::endl;
      return;
    }
    else{

      //Now check the channel statuses        
      std::cout << "ChStatuses: " << CTB.GetNChStatuses() << std::endl;
      raw::ctb::ChStatus theStatus = CTB.GetChStatuse(0);
      
      uint32_t the_beam_hi    = theStatus.beam_hi; 
      uint32_t the_beam_lo    = theStatus.beam_lo; 
      ULong64_t the_timestamp = theStatus.timestamp; 
      
      std::cout << "Timestamp : " << the_timestamp << std::endl;
      if( 2.e-8*the_timestamp > HLTTS ){
        std::cout << "Found Channel status > HLT timestamp. Skipping" << std::endl;
      }

      std::cout << "beam_hi   : " << std::bitset<5>(the_beam_hi) << std::endl; 
      std::cout << "beam_lo   : " << std::bitset<4>(the_beam_lo) << std::endl; 

      std::bitset<4> beam_lo(the_beam_lo);
      std::bitset<5> beam_hi(the_beam_hi);

      BITrigger  =  beam_lo[1];
      C1         =  beam_lo[3];
      C2         =  beam_hi[0];


      std::cout << "%%%%Decoding the beam channels%%%" << std::endl;
      std::cout << "BI Trigger: " << BITrigger << std::endl
                << "C1:         " << C1        << std::endl
                << "C2:         " << C2        << std::endl
                << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << std::endl << std::endl;        

      return;
    }
    
  }

  std::cout << "Error! Invalid CTB Handle!" << std::endl;
  std::cout << std::endl;
  return;

}



////////////////////////
// Analyzer Method (reads in the event and derives values)
void proto::BeamCounter::analyze(art::Event const & e)
{
  //Reset 
  HLTWord = 0;
  HLTTS      = -1;
  BITrigger  = -1;
  C1         = -1;
  C2         = -1;

  eventNum = e.event();
  eventTime = e.time().timeHigh() + 1.e-9*e.time().timeLow(); 
  runNum = e.run();
  subRunNum = e.subRun();
   
  GetCTBInfo(e);

  // Write out the to tree
  fOutTree->Fill();
}


void proto::BeamCounter::beginJob()
{
  art::ServiceHandle<art::TFileService> tfs;
  fOutTree = tfs->make<TTree>("tree", ""); 

  fOutTree->Branch("Time",   &eventTime);
  fOutTree->Branch("Event",  &eventNum);
  fOutTree->Branch("Run",    &runNum);
  fOutTree->Branch("Subrun", &subRunNum);

  fOutTree->Branch("C1",        &C1);
  fOutTree->Branch("C2",        &C2);
  fOutTree->Branch("BITrigger", &BITrigger);
}

DEFINE_ART_MODULE(proto::BeamCounter)
