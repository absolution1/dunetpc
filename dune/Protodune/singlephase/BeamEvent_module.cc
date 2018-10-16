////////////////////////////////////////////////////////////////////////
// Class:       BeamEvent
// Plugin Type: producer (art v2_08_03)
// File:        BeamEvent_module.cc
//
// Generated at Thu Nov  2 22:57:41 2017 by Jonathan Paley using cetskelgen
// from cetlib version v3_01_01.
//
// Edited by Jake Calcutt
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDProducer.h"
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
#include "dune/DuneObj/ProtoDUNEBeamEvent.h"
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
  class BeamEvent;
}

enum tofChan{
  k1A2A,
  k1B2A,
  k1A2B,
  k1B2B
};

typedef std::numeric_limits< double > dbl;

class proto::BeamEvent : public art::EDProducer {
public:
  explicit BeamEvent(fhicl::ParameterSet const & p);
  //  virtual ~BeamEvent();

  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  BeamEvent(BeamEvent const &) = delete;
  BeamEvent(BeamEvent &&) = delete;
  BeamEvent & operator = (BeamEvent const &) = delete;
  BeamEvent & operator = (BeamEvent &&) = delete;

  // Required functions.
  void reconfigure(fhicl::ParameterSet const & p);
  void produce(art::Event & e) override;

  // Selected optional functions.
  void beginJob() override;
  void beginRun(art::Run & r) override;
  void beginSubRun(art::SubRun & sr) override;
  void endJob() override;
  void endRun(art::Run & r) override;
  void endSubRun(art::SubRun & sr) override;


  //std::bitset<sizeof(double)*CHAR_BIT> toBinary(const long num);  
  uint64_t joinHighLow(double,double);

  TVector3 ConvertProfCoordinates(double x, double y, double z, double zOffset);
  void BeamMonitorBasisVectors(); 
  void RotateMonitorVector(TVector3 &vec); 

  uint64_t GetRawDecoderInfo(art::Event &);
  void TimeIn(art::Event &, uint64_t);
  void GetSpillInfo(art::Event &);
  bool MatchBeamToTPC(art::Event &, uint64_t);

  void MakeTrack(size_t);
  void MomentumSpec(size_t);
  double MomentumCosTheta(double, double, double);

  void GetPairedFBMInfo(beam::ProtoDUNEBeamEvent beamevt, double Time);
  void GetPairedStraightFBMInfo(beam::ProtoDUNEBeamEvent beamevt, double Time);
  void GetUnpairedFBMInfo(beam::ProtoDUNEBeamEvent beamevt, double Time);

  double GetPosition(std::string, int);
  TVector3 ProjectToTPC(TVector3, TVector3);
  double GetPairedPosition(std::string, size_t);
  TVector3 TranslateDeviceToDetector(TVector3);
 
  void  InitXBPFInfo(beam::ProtoDUNEBeamEvent *);
  void  parseGeneralXBPF(std::string, uint64_t, size_t);
  void  parseXBPF(uint64_t);
  void  parsePairedXBPF(uint64_t);
  void  parsePairedStraightXBPF(uint64_t);

  void  parseXTOF(uint64_t);
  void  parseXCET(uint64_t);

  template<class T> 
  T FetchWithRetries(uint64_t, std::string, int);
   
private:
  
  TTree * fOutTree;
  TTree * fTrackTree;

  std::vector<double> * trackX;
  std::vector<double> * trackY;
  std::vector<double> * trackZ;

  TH2F * fFirstBeamProf2D;
  TH2F * fSecondBeamProf2D;
  std::vector<TH2F *> fBeamProf2D;
  std::map<std::string, TH1F *> fBeamProf1D;
  std::map<std::string, std::vector<short> *> fActiveFibers;
  std::map<std::string, double> fProfTime;
  std::map<std::string, double> fProfTrigger1;
  std::map<std::string, double> fProfTrigger2;
  std::map<std::string, double> fProfTime1;
  std::map<std::string, double> fProfTime2;
  std::map<std::string, TTree*> fProfTree;
  TH1F  * fMomentumHist;
  TH1F  * fFullMomentum;
  TH1F  * fCutMomentum;
  TTree * fGenTrigTree;
  TTree * fXTOF1ATree;
  TTree * fXTOF1BTree;
  TTree * fXTOF2ATree;
  TTree * fXTOF2BTree;
  TTree * fXBTF022702Tree;
  TTree * fMatchedTriggers;

  TH1F  * fDeltaXHist;
  TH1F  * fDeltaYHist;
  TH1F  * fDeltaZHist;
  TH2F  * fDeltaXYHist;
  TH2F  * fDeltaYZHist;
  TH2F  * fDeltaZXHist;
  TTree * fDeltaTree;
  double  fDeltaX;
  double  fDeltaY;
  double  fDeltaZ;
  double  fBeamPosX;
  double  fBeamPosY;
  double  fBeamPosZ;
  double  fTrackPosX;
  double  fTrackPosY;
  double  fTrackPosZ;
  int trackNum;
  int beamNum;

  double matchedGen;
  double matchedTOF1;
  double matchedTOF2;
  int    matchedChan;
  int    matchedNom;
  std::map<std::string, double> matchedXBPF;
  

  double fGenTrigFrac;
  double fGenTrigSec;
  double fGenTrigCoarse; 
  double fXTOF1AFrac;
  double fXTOF1ACoarse; 
  double fXTOF1BFrac;
  double fXTOF1BCoarse; 
  double fXTOF2AFrac;
  double fXTOF2ACoarse; 
  double fXTOF2BFrac;
  double fXTOF2BCoarse; 

  double fXTOF1ASec;
  double fXTOF1BSec;
  double fXTOF2ASec;
  double fXTOF2BSec;

  double fXBTF022702Frac;
  double fXBTF022702Coarse; 
  TH1F * fTOFHist;
  TH1F * fCKovHist;
  recob::Track * theTrack;
  std::vector<recob::Track*> theTracks;
  long long int eventTime;
  double SpillStart;
  ULong_t SpillStart_alt;
  double PrevStart;
  double SpillEnd;
  double SpillOffset;
  double ActiveTriggerTime;
  double RDTSTime;
  std::vector< double > * GenTriggers;

  double acqTime;
  double acqStampMBPL;
  double cycleStampMBPL;
  int HLTWord;
  long long int HLTTS;
  int BeamOn;
  int BITrigger;
  int Upstream;
  int C1;
  int C2;
  int BP1;
  int BP2;
  int BP3;
  int BP4;

  uint64_t prev_fetch_time;
  long long int prev_event_time;

  int eventNum;
  int runNum;
  int subRunNum;
  double CKov1Pressure;
  double CKov2Pressure;
  double CKov1Efficiency;
  double CKov2Efficiency;
  
  TVector3 fBMBasisX;
  TVector3 fBMBasisY;
  TVector3 fBMBasisZ;

  // Declare member data here.
  //bool fLoadFromDB; // unused
  double  fTimeWindow;
  double fTolerance;
  std::string fCSVFileName;
  std::string fBundleName;
  std::string fOutputLabel;
  std::string fInputLabel;
  double fDummyEventTime;
  std::string fURLStr;
  double fValidWindow;
  uint64_t fFixedTime;
  std::vector< uint64_t > fMultipleTimes;

  std::string firstUpstreamName;
  std::string secondUpstreamName;
  std::string firstDownstreamName;
  std::string secondDownstreamName;

  std::string firstBPROF1;
  std::string secondBPROF1;
  std::string BPROF2;
  std::string BPROF3;
  double      fBeamBend;

  std::vector< std::string > fDevices;
  std::vector< std::pair<std::string, std::string> > fPairedDevices;
  std::vector< std::pair<std::string, std::string> > fPairedStraightDevices;
  std::map<std::string, std::string > fDeviceTypes;
//  std::vector< std::array<double, 3> > fCoordinates; 
  std::map< std::string, std::array<double,3> > fCoordinates;
  std::map< std::string, std::array<double, 3> > fRotations;
  TVector3 fGlobalDetCoords;
  std::array<double,3> fDetRotation;
  std::map< std::string, double > fFiberDimension;

  int fNRetries;

  // Names of the CERN Beam Devices
  // Names have the form of a prefix + device
  
  // Prefixes for different devices classes:
  // Beam Positions (BPF)
  // Time of flights (TOF)
  // and Cerenkov counters (CET)
  std::string fXBPFPrefix;
  std::string fXTOFPrefix;
  std::string fXCETPrefix;

  // Device Names for each device

  // Time of Flights
  std::string fTOF1;
  std::string fTOF1A, fTOF1B;
  std::string fTOF2;
  std::string fTOF2A, fTOF2B;

  // Cerenkovs
  std::string fCKov1;
  std::string fCKov2;

  double fRotateMonitorXZ;
  double fRotateMonitorYZ;

  double fFirstTrackingProfZ;
  double fSecondTrackingProfZ;
  double fNP04FrontZ;  
  double fBeamX, fBeamY, fBeamZ;

  bool   fForceNewFetch;
  bool   fMatchTime;
  bool   fForceRead;

  double fTimingCalibration;
  double fCalibrationTolerance;
  double fOffsetTAI;

  beam::ProtoDUNEBeamEvent * beamevt;
  beam::ProtoDUNEBeamEvent prev_beamevt;

  std::unique_ptr<ifbeam_ns::BeamFolder> bfp;
  art::ServiceHandle<ifbeam_ns::IFBeam> ifb;

  art::Handle< std::vector<raw::RDTimeStamp> > RDTimeStampHandle;
//  art::Handle< std::vector<raw::ctb::pdspctb> > CTBHandle;

  uint64_t validTimeStamp;

  double L1=1.980, L2=1.69472, L3=2.11666;
  double magnetLen, magnetField;
  std::vector<double> current;

  //Hardware Parameters for magnetic field stuff
  double mag_P1 = 5.82044830e-3;
  double mag_P2 = 0.;
  double mag_P3 = -4.68880000e-6;
  double mag_P4 = 324.573967;

};

// Constructor
proto::BeamEvent::BeamEvent(fhicl::ParameterSet const & p)
{
  // Declare products this module will provide
  produces<std::vector<beam::ProtoDUNEBeamEvent>>();  

  // Configure/Reconfigure
  this->reconfigure(p);
}
// END Constructor
////////////////////////

////////////////////////
// Fetch Method
template <class T> 
T proto::BeamEvent::FetchWithRetries(uint64_t time, std::string name, int nRetry){
  T theResult;
  
  std::cout << std::endl;
  uint64_t newTime;
  //Search at and above time given with nRetries
  //Will later search below, just in case the event time is actually greater than
  for(newTime = time; newTime < time + nRetry; ++newTime){
//    std::cout << "Trying to grab from folder: " << name << std::endl;
//    std::cout << "At Time: " << newTime << std::endl;    
      try{
        theResult = bfp->GetNamedVector(newTime, name);
        std::cout << "Successfully fetched" << std::endl;
        prev_fetch_time = newTime;
        return theResult;
      }
    catch(std::exception e){
//      std::cout << "Could not fetch with time " << newTime << std::endl;      
    }
  }
  //Now search below
  for(newTime = time - 1; newTime > time - nRetry - 1; --newTime){
//    std::cout << "Trying to grab from folder: " << name << std::endl;
//    std::cout << "At Time: " << newTime << std::endl;    
    try{
      theResult = bfp->GetNamedVector(newTime, name);
      std::cout << "Successfully fetched" << std::endl;
      prev_fetch_time = newTime;
      return theResult;
    }
    catch(std::exception e){
//      std::cout << "Could not fetch with time " << newTime << std::endl;      
    }
  }
  
  //Try a final time. Let it crash if it doesn't work
  std::cout << "Trying a final time to grab from folder: " << name << std::endl;
  std::cout << "At time: " << newTime << std::endl;
  theResult = bfp->GetNamedVector(newTime, name);
  std::cout << "Successfully fetched" << std::endl;
  std::cout << std::endl;
  return theResult; 
}
// END FetchWithRetries
////////////////////////


//Gets the Timing and CTB raw-decoded info.
//Finds the triggers, and looks for a valid trigger word
//(i.e. coming from beam)
//
//Returns the timestamp of the high level trigger.
uint64_t proto::BeamEvent::GetRawDecoderInfo(art::Event & e){
  std::cout << std::endl;
  std::cout << "Getting Raw Decoder Info" << std::endl;
  e.getByLabel("timingrawdecoder","daq",RDTimeStampHandle);
  std::cout << "RDTS valid? " << RDTimeStampHandle.isValid() << std::endl;
  for (auto const & RDTS : *RDTimeStampHandle){
    std::cout << "High: " << RDTS.GetTimeStamp_High() << std::endl;
    std::cout << "Low: " << RDTS.GetTimeStamp_Low() << std::endl; 

    std::bitset<64> high = RDTS.GetTimeStamp_High();
    std::bitset<64> low  = RDTS.GetTimeStamp_Low();
    high = high << 32; 
    std::bitset<64> joined = (high ^ low);

    std::cout << "Join: " << (joined).to_ullong() << std::endl;
    std::cout << "Join: " << 2.e-8 * (joined).to_ullong() << std::endl;

    RDTSTime = 2.e-8 * joined.to_ullong(); 

  }

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
        HLTTS = theTS;
        
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
      return 0;
    }
    else{

      //Now check the channel statuses        
      std::cout << "ChStatuses: " << CTB.GetNChStatuses() << std::endl;
      raw::ctb::ChStatus theStatus = CTB.GetChStatuse(0);
      
      uint32_t the_beam_hi    = theStatus.beam_hi; 
      uint32_t the_beam_lo    = theStatus.beam_lo; 
      ULong64_t the_timestamp = theStatus.timestamp; 
      
      std::cout << "Timestamp : " << the_timestamp << std::endl;
      if( the_timestamp > ULong64_t(HLTTS) ){
        std::cout << "Found Channel status > HLT timestamp. Skipping" << std::endl;
//        continue;
      }

      std::cout << "beam_hi   : " << std::bitset<5>(the_beam_hi) << std::endl; 
      std::cout << "beam_lo   : " << std::bitset<4>(the_beam_lo) << std::endl; 

      std::bitset<4> beam_lo(the_beam_lo);
      std::bitset<5> beam_hi(the_beam_hi);

      BeamOn     =  beam_lo[0];
      BITrigger  =  beam_lo[1];
      Upstream   =  beam_lo[2];
      C1         =  beam_lo[3];
      C2         =  beam_hi[0];
      BP1        =  beam_hi[1];
      BP2        =  beam_hi[2];
      BP3        =  beam_hi[3];
      BP4        =  beam_hi[4];


      std::cout << "%%%%Decoding the beam channels%%%" << std::endl;
      std::cout << "Beam On:    " << BeamOn    << std::endl
                << "BI Trigger: " << BITrigger << std::endl
                << "Upstream:   " << Upstream  << std::endl
                << "C1:         " << C1        << std::endl
                << "C2:         " << C2        << std::endl
                << "BP1:        " << BP1       << std::endl
                << "BP2:        " << BP2       << std::endl
                << "BP3:        " << BP3       << std::endl
                << "BP4:        " << BP4       << std::endl;
      std::cout << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << std::endl << std::endl;        

      fCKovHist->Fill(C1 + 2*C2);

      return HLTTS;
    }
    
  }
  std::cout << "Error! Invalid CTB Handle!" << std::endl;
  std::cout << std::endl;
  return 0;

}

void proto::BeamEvent::GetSpillInfo(art::Event & e){

    std::cout << "Attempting to Time In" << std::endl;   
    auto const PDTStampHandle = e.getValidHandle< std::vector< dune::ProtoDUNETimeStamp > > ("timingrawdecoder:daq");  
    std::vector< dune::ProtoDUNETimeStamp > PDTStampVec(*PDTStampHandle); 
    auto PDTStamp = PDTStampVec[0];            
          
    UInt_t    ver        = PDTStamp.getVersion();   
    double RunStart   = 2.e-8*PDTStamp.getLastRunStart();            
    SpillStart = 2.e-8*PDTStamp.getLastSpillStart();          
    SpillEnd   = 2.e-8*PDTStamp.getLastSpillEnd();            
    SpillStart_alt = PDTStamp.getLastSpillStart();
    std::cout << PDTStamp.getLastSpillStart() << " " << SpillStart_alt << std::endl;
          
    std::cout.precision(dbl::max_digits10);
    std::cout << "Version:     " << ver        << std::endl;      
    std::cout << "Run:         " << RunStart   << std::endl;      
          
    std::cout << "Spill Start: " << SpillStart <<  std::endl;      
    std::cout << "Spill End:   " << SpillEnd   <<  std::endl;      
          
    std::cout << std::endl;           
          
    if( SpillEnd > SpillStart ){      
      std::cout << "Outside of spill" << std::endl; 
      std::cout << "End - Start: "    << SpillEnd - SpillStart << std::endl;        
    }     
    else{ 
      std::cout << "Within a spill" << std::endl;   
    }     
          
    std::cout << std::endl;           
          
}

void proto::BeamEvent::TimeIn(art::Event & e, uint64_t time){

    std::cout << "Attempting to Time In" << std::endl;   

    /////Now look at the cycleStamp and acqStamp coming out of IFBeam
    std::vector<double> acqStamp = FetchWithRetries< std::vector<double> >(time, "dip/acc/NORTH/NP04/POW/MBPL022699:acqStamp[]",fNRetries); 
    std::vector<double> cycleStamp = FetchWithRetries< std::vector<double> >(time, "dip/acc/NORTH/NP04/POW/MBPL022699:cycleStamp[]",fNRetries); 

    acqStampMBPL   = 1.e-9 * joinHighLow(acqStamp[0],   acqStamp[1]); 
    cycleStampMBPL = 1.e-9 * joinHighLow(cycleStamp[0], cycleStamp[1]);

    std::cout << "%%%IFBeam%%%" << std::endl;
    std::cout << std::setw(15) << acqStampMBPL*1.e-9   << std::endl;
    std::cout << std::setw(15) << cycleStampMBPL*1.e-9 << std::endl;

    std::cout << std::endl;


    //Assign the calibration offset
    SpillOffset = SpillStart - acqStampMBPL;
}

bool proto::BeamEvent::MatchBeamToTPC(art::Event & e, uint64_t time){

  std::cout << "Matching in time between Beamline and TPC!!!" << std::endl; 
  for(size_t iT = 0; iT < beamevt->GetNT0(); ++iT){
    
    //GenTrig = sec + 1.e-9*ns portions
    double GenTrigTime = beamevt->GetT0(iT).first + 1.e-09*beamevt->GetT0(iT).second;

    //HLTTime in 50MHz ticks
    double HLTTime = 2.e-08*HLTTS;
    double diff = HLTTime - GenTrigTime - SpillOffset;
    std::cout.precision(dbl::max_digits10);
    std::cout << GenTrigTime << " " << HLTTime << " " << diff << " " << fTimingCalibration << std::endl << std::endl;

  
    if( ( fTimingCalibration - fCalibrationTolerance < diff ) && (fTimingCalibration + fCalibrationTolerance > diff) ){
      std::cout << "FOUND MATCHING TIME!!!" << std::endl;


      beamevt->SetActiveTrigger( iT ); 
      return true;
    }
  }

  std::cout << "Could not find matching time " << std::endl << std::endl;
  beamevt->SetUnmatched();
  return false;
}

////////////////////////
// Producer Method (reads in the event and derives values)
void proto::BeamEvent::produce(art::Event & e)
{
  //Reset 
  acqTime = 0;
  acqStampMBPL = 0;
  cycleStampMBPL = 0;
  HLTWord = 0;
  HLTTS = 0;
  BeamOn     = -1;
  BITrigger  = -1;
  Upstream   = -1;
  C1         = -1;
  C2         = -1;
  BP1        = -1;
  BP2        = -1;
  BP3        = -1;
  BP4        = -1; 
  SpillStart = -1;
  SpillEnd   = -1;
  SpillOffset = -1;
  ActiveTriggerTime = -1;
  RDTSTime   = -1;
  GenTriggers->clear();

  bool usedEventTime = false;


  eventNum = e.event();
  runNum = e.run();
  subRunNum = e.subRun();

  // Now for the current event retrieve the time (trigger) of the event
  // This will take the form a 64-bit timestamp with a high and low word
  // Break this up to have access to each word and the long long word
  std::cout <<"Event Time: " << uint64_t(e.time().timeLow()) << std::endl;
  std::cout << "Low: " << e.time().timeLow()  << std::endl;
  std::cout << "High: " << e.time().timeHigh()  << std::endl;
  std::cout << std::endl;
  //Need to use either for different versions
  if( e.time().timeHigh() == 0 ) eventTime = e.time().timeLow();
  else eventTime = e.time().timeHigh();


   
  // Create a new beam event (note the "new" here)  
  beamevt = new beam::ProtoDUNEBeamEvent();

  // Get the coordinate system conversions
  BeamMonitorBasisVectors();

  validTimeStamp = GetRawDecoderInfo(e);
  std::cout << std::endl;

  //Get Spill Information
  //This stores Spill Start, Spill End, 
  //And Prev Spill Start
  GetSpillInfo(e);

  if(~SpillStart_alt == 0ul){
    std::cout << "Invalid Spill Start time! Skipping Event" << std::endl << std::endl;
  }
   
  //Check if we have a valid beam trigger
  //If not, just place an empty beamevt
  //and move on
  //
  //Also check if we've gotten good spill info
  //
  //Or if we're forcing to read out the Beamline Info
  if( (validTimeStamp && ( 0ul != ~(SpillStart_alt) ) ) || fForceRead ){


    // Read in and cache the beam bundle folder for a specific time
    bfp = ifb->getBeamFolder(fBundleName,fURLStr,fTimeWindow);
    std::cerr << "%%%%%%%%%% Got beam folder %%%%%%%%%%" << std::endl << std::endl;


    //Use multiple times provided to fcl
    if( fMultipleTimes.size() ){
      std::cout << "Using multiple times: " << std::endl;
      for(size_t it = 0; it < fMultipleTimes.size(); ++it){
        std::cout << fMultipleTimes[it] << std::endl;
      }
    }
    //Use singular time provided to fcl
    else if(fFixedTime > 0.){
      std::cout << "Using singular time: " << fFixedTime << std::endl;
      fMultipleTimes.push_back(fFixedTime);
    }
    //Use time of event
    else{
      std::cout <<" Using Event Time: " << uint64_t(eventTime) << std::endl;
      fMultipleTimes.push_back(uint64_t(eventTime));
      usedEventTime = true;
    } 




    //Start getting beam event info
    // Loop over the different particle times 
    for(size_t it = 0; it < fMultipleTimes.size(); ++it){
      std::cout << "Time: " << fMultipleTimes[it] << std::endl;
      uint64_t fetch_time;
      if( abs( (long long)(fMultipleTimes[it]) - (long long)(prev_fetch_time) ) < 18 ){
        fetch_time = prev_fetch_time;
      }
      else{
        fetch_time = fMultipleTimes[it];
      }
     

      //Check if we are still using the same spill information.
      //
      //If it's a new spill: Get new info from the database 
      //Each 'parse' command fetches new data from the database
      //
      //If it's the same spill then just pass the old BeamEvent
      //Object. Its 'active trigger' info will be updated below
      //
      //Can be overridden with a flag from the fcl
      if(PrevStart != SpillStart || fForceNewFetch){
        std::cout << "New spill or forced new fetch. Getting new beamevt info" << std::endl << std::endl;


        // Parse the Time of Flight Counter data for the list
        // of times that we are using
        try{
          parseXTOF(fetch_time);
        }
        catch(std::exception e){
          std::cout << "COULD NOT GET INFO" << std::endl;
          std::cout << "SKIPPING EVENT" << std::endl;
          //break;
        }
        std::cout << std::endl;

        // Parse the Beam postion counter information for the list
        // of time that we are using
        InitXBPFInfo(beamevt);
        parseXBPF(fetch_time);
        parsePairedXBPF(fetch_time);
        parsePairedStraightXBPF(fetch_time);
        parseXCET(fetch_time);
   
        //Set PrevStart to SpillStart here
        //so that we don't skip in the case 
        //the first event in the spill did not
        //a good beamline trigger
        PrevStart = SpillStart;

      }
      else{
        std::cout << "Same spill. Reusing beamevt info" << std::endl << std::endl;
        std::cout << prev_beamevt.GetNT0() << std::endl;
        *beamevt = prev_beamevt;
        std::cout << beamevt->GetNT0() << std::endl;;
      }

      std::cout << "NGoodParticles: " << beamevt->GetNT0()            << std::endl;
      std::cout << "NTOF0: "          << beamevt->GetNTOF0Triggers()  << std::endl;
      std::cout << "NTOF1: "          << beamevt->GetNTOF1Triggers()  << std::endl;
      std::cout << "acqTime: "        << acqTime                      << std::endl;
      std::cout << "NXBPF: "          << beamevt->GetNFBMTriggers(fDevices[0]) << std::endl;



      if( fMatchTime ){
   
        //Now do the matching in Time:
        //Get the conversion from the ProtoDUNE Timing system
        //To the one in the SPS.
        //
        TimeIn(e, fetch_time);
        std::cout << "SpillOffset " << SpillOffset << std::endl;
        bool matched = MatchBeamToTPC(e, fetch_time);


        std::cout << matched << std::endl;
        if( beamevt->CheckIsMatched() ){
          std::pair<double,double> theTime = beamevt->GetT0(beamevt->GetActiveTrigger());
          ActiveTriggerTime = theTime.first + theTime.second*1.e-9;
          std::cout << "Trigger: " << beamevt->GetActiveTrigger() << " " << ActiveTriggerTime << std::endl << std::endl;       

          MakeTrack( beamevt->GetActiveTrigger() );
          for(size_t iTrack = 0; iTrack < theTracks.size(); ++iTrack){    
            beamevt->AddBeamTrack( *(theTracks[iTrack]) ); 
            
            auto thisTrack = theTracks[iTrack]; 
            const recob::TrackTrajectory & theTraj = thisTrack->Trajectory();
            trackX->clear(); 
            trackY->clear(); 
            trackZ->clear(); 
            std::cout << "npositison: " << theTraj.NPoints() << std::endl;
            for(auto const & pos : theTraj.Trajectory().Positions()){
              std::cout << pos.X() << " " << pos.Y() << " " << pos.Z() << std::endl;
              trackX->push_back(pos.X());
              trackY->push_back(pos.Y());
              trackZ->push_back(pos.Z());
            }

            fTrackTree->Fill();
          }
          std::cout << "Added " << beamevt->GetNBeamTracks() << " tracks to the beam event" << std::endl << std::endl;

          //Momentum
          //First, try getting the current from the magnet in IFBeam
          std::cout << "Trying to get the current" << std::endl;
          try{
            current = FetchWithRetries< std::vector<double> >(fetch_time, "dip/acc/NORTH/NP04/POW/MBPL022699:current",fNRetries);
            std::cout << "Current: " << current[0] << std::endl;
    
            MomentumSpec( beamevt->GetActiveTrigger() ); 
            std::cout << "Got NRecoBeamMomenta: " << beamevt->GetNRecoBeamMomenta() << std::endl << std::endl;
          }
          catch(std::exception e){
            std::cout << "Could not get the current from the magnet. Skipping spectrometry" << std::endl;
          }
          
        }
      }
    }

    //Pass beamevt to the next event;
    //Erase the Track and Reco Momentum info
    prev_beamevt = *beamevt;
    prev_beamevt.ClearBeamTracks();
    prev_beamevt.ClearRecoBeamMomenta();
    prev_beamevt.SetUnmatched();

  }
  //Start of a new spill, but the first event was 
  //Not a beam trigger. In this case, it would not
  //have been filled with info in the block above
  //So let's make it empty so we aren't putting 
  //old spill info in the new event
  if(!validTimeStamp && PrevStart != SpillStart){
    prev_beamevt = *beamevt;   
  }
  
  beamevt->SetBITrigger(BITrigger);
  beamevt->SetSpillStart(SpillStart);
  beamevt->SetSpillOffset(SpillOffset);
  beamevt->SetCTBTimestamp( (validTimeStamp == 0 ? -1. : 2.e-8*validTimeStamp ) );

  std::unique_ptr<std::vector<beam::ProtoDUNEBeamEvent> > beamData(new std::vector<beam::ProtoDUNEBeamEvent>);
  beamData->push_back(beam::ProtoDUNEBeamEvent(*beamevt));

  delete beamevt;
  e.put(std::move(beamData));
 
  // Write out the to tree
  fOutTree->Fill();
 
  prev_event_time = eventTime; 
  theTracks.clear();
  current.clear();
  if(usedEventTime) fMultipleTimes.clear();
}
// END BeamEvent::produce
////////////////////////

void proto::BeamEvent::InitXBPFInfo(beam::ProtoDUNEBeamEvent * beamevt){
  // Places a dummy trigger vector for each device

  // Make a vector the names of each of the devices being readout
  std::vector<std::string> monitors;
  size_t nDev = 0; 

  for(size_t id = 0; id < fDevices.size(); ++id){
    std::string name = fDevices[id];
    std::cout << fXBPFPrefix + fDevices[id] << std::endl;
    std::cout << "At: "      << fCoordinates[name][0] << " " << fCoordinates[name][1] << " " << fCoordinates[name][2] << std::endl;
    std::cout << "Rotated: " << fRotations[name][0]   << " " << fRotations[name][1]   << " " << fRotations[name][2]   << std::endl;

    // Put the current device name on the list
    monitors.push_back(name);
    nDev++;
  }

  for(size_t id = 0; id < fPairedStraightDevices.size(); ++id){

    std::string name = fPairedStraightDevices[id].first;
    std::cout << fXBPFPrefix + name << std::endl;
    std::cout << "At: " << fCoordinates[name][0] << " " << fCoordinates[name][1] << " " << fCoordinates[name][2] << std::endl;

    std::cout << nDev << std::endl;
    nDev++;
    monitors.push_back(name);

    name = fPairedStraightDevices[id].second;
    std::cout << fXBPFPrefix + name << std::endl;
    std::cout << "At: " << fCoordinates[name][0] << " " << fCoordinates[name][1] << " " << fCoordinates[name][2] << std::endl;
    monitors.push_back(name);
    std::cout << nDev << std::endl;
    nDev++;
  }

  for(size_t id = 0; id < fPairedDevices.size(); ++id){

    std::string name = fPairedDevices[id].first;
    std::cout << fXBPFPrefix + name << std::endl;
    std::cout << "At: " << fCoordinates[name][0] << " " << fCoordinates[name][1] << " " << fCoordinates[name][2] << std::endl;

    std::cout << nDev << std::endl;
    nDev++;
    monitors.push_back(name);

    name = fPairedDevices[id].second;
    std::cout << fXBPFPrefix + name << std::endl;
    std::cout << "At: " << fCoordinates[name][0] << " " << fCoordinates[name][1] << " " << fCoordinates[name][2] << std::endl;

    std::cout << nDev << std::endl;
    nDev++;
    monitors.push_back(name);
  }

//  beamevt->InitFBMs(nDev);
  std::cout << "Initializing monitors";
  for(size_t id = 0; id < nDev; ++id){
    std::cout << " " << monitors.at(id);
  }
  std::cout << std::endl;
  beamevt->InitFBMs(monitors);
}
// END BeamEvent::InitXBPFInfo
////////////////////////


void proto::BeamEvent::parseXTOF(uint64_t time){
  std::cout << "Getting General trigger info " << std::endl;
  std::vector<double> coarseGeneralTrigger = FetchWithRetries< std::vector<double> >(time, "dip/acc/NORTH/NP04/BI/XBTF/GeneralTrigger:coarse[]",fNRetries);
  std::vector<double> fracGeneralTrigger = FetchWithRetries< std::vector<double> >(time, "dip/acc/NORTH/NP04/BI/XBTF/GeneralTrigger:frac[]",fNRetries); 
  std::vector<double> acqStampGeneralTrigger = FetchWithRetries< std::vector<double> >(time, "dip/acc/NORTH/NP04/BI/XBTF/GeneralTrigger:acqStamp[]",fNRetries); 
  std::vector<double> secondsGeneralTrigger = FetchWithRetries< std::vector<double> >(time, "dip/acc/NORTH/NP04/BI/XBTF/GeneralTrigger:seconds[]",fNRetries); 
  std::vector<double> timestampCountGeneralTrigger = FetchWithRetries< std::vector<double> >(time, "dip/acc/NORTH/NP04/BI/XBTF/GeneralTrigger:timestampCount",fNRetries); 
  std::cout << "Size of coarse,frac: " << coarseGeneralTrigger.size() << " " << fracGeneralTrigger.size() << std::endl; 

  std::cout <<"Size of acqStamp: " << acqStampGeneralTrigger.size() << std::endl;
  std::cout <<"Size of seconds: " << secondsGeneralTrigger.size() << std::endl;
  std::cout << "Size of counts: " << timestampCountGeneralTrigger.size() << std::endl;
  std::cout << "timestampCounts: " << timestampCountGeneralTrigger[0] << std::endl;
  matchedNom = int(timestampCountGeneralTrigger[0]);
  
  double low = acqStampGeneralTrigger[1];
  uint32_t low32 = (uint32_t)low;
  std::bitset<64> lowbits = low32;

  double high = acqStampGeneralTrigger[0];
  uint32_t high32 = (uint32_t)high;
  std::bitset<64> highbits = high32;

  highbits = highbits << 32;
  std::bitset<64> joinedbits = highbits ^ lowbits;

  acqTime = joinedbits.to_ullong() / 1000000000.; 

  std::cout << "Getting TOF1A info: " << fTOF1 << std::endl;
  std::vector<double> coarseTOF1A = FetchWithRetries< std::vector<double> >(time, fXTOFPrefix + fTOF1A + ":coarse[]",fNRetries);
  std::vector<double> fracTOF1A =   FetchWithRetries< std::vector<double> >(time, fXTOFPrefix + fTOF1A + ":frac[]",fNRetries);
  std::vector<double> secondsTOF1A =   FetchWithRetries< std::vector<double> >(time, fXTOFPrefix + fTOF1A + ":seconds[]",fNRetries);
  std::cout << "Size of coarse,frac: " << coarseTOF1A.size() << " " << fracTOF1A.size() << std::endl; 

  std::cout << "Getting TOF1B info: " << fTOF1 << std::endl;
  std::vector<double> coarseTOF1B = FetchWithRetries< std::vector<double> >(time, fXTOFPrefix + fTOF1B + ":coarse[]",fNRetries);
  std::vector<double> fracTOF1B = FetchWithRetries< std::vector<double> >(time, fXTOFPrefix + fTOF1B + ":frac[]",fNRetries);
  std::vector<double> secondsTOF1B = FetchWithRetries< std::vector<double> >(time, fXTOFPrefix + fTOF1B + ":seconds[]",fNRetries);
  std::cout << "Size of coarse,frac: " << coarseTOF1B.size() << " " << fracTOF1B.size() << std::endl; 

  std::cout << "Getting TOF2A info: " << fTOF2 << std::endl;
  std::vector<double> coarseTOF2A = FetchWithRetries< std::vector<double> >(time, fXTOFPrefix + fTOF2A + ":coarse[]",fNRetries);
  std::vector<double> fracTOF2A = FetchWithRetries< std::vector<double> >(time, fXTOFPrefix + fTOF2A + ":frac[]",fNRetries);
  std::vector<double> secondsTOF2A = FetchWithRetries< std::vector<double> >(time, fXTOFPrefix + fTOF2A + ":seconds[]",fNRetries);
  std::cout << "Size of coarse,frac: " << coarseTOF2A.size() << " " << fracTOF2A.size() << std::endl; 

  std::cout << "Getting TOF2B info: " << fTOF2 << std::endl;
  std::vector<double> coarseTOF2B = FetchWithRetries< std::vector<double> >(time, fXTOFPrefix + fTOF2B + ":coarse[]",fNRetries);
  std::vector<double> fracTOF2B = FetchWithRetries< std::vector<double> >(time, fXTOFPrefix + fTOF2B + ":frac[]",fNRetries);
  std::vector<double> secondsTOF2B = FetchWithRetries< std::vector<double> >(time, fXTOFPrefix + fTOF2B + ":seconds[]",fNRetries);
  std::cout << "Size of coarse,frac: " << coarseTOF2B.size() << " " << fracTOF2B.size() << std::endl; 

  std::vector<std::pair<double,double>> unorderedGenTrigTime;
  std::vector<std::pair<double,double>> unorderedXBTF022702Time;
  std::vector<std::pair<double,double>> unorderedTOF1ATime;
  std::vector<std::pair<double,double>> unorderedTOF1BTime;
  std::vector<std::pair<double,double>> unorderedTOF2ATime;
  std::vector<std::pair<double,double>> unorderedTOF2BTime;


  for(size_t i = 0; i < timestampCountGeneralTrigger[0]; ++i){
//    std::cout << i << " " << secondsGeneralTrigger[2*i + 1] << " "  << 8.*coarseGeneralTrigger[i] + fracGeneralTrigger[i]/512. << std::endl;
//    std::cout << "\t" << std::setw(15) << secondsGeneralTrigger[2*i + 1] + 1.e-9*(8.*coarseGeneralTrigger[i] + fracGeneralTrigger[i]/512.) << std::endl;
    fGenTrigCoarse = coarseGeneralTrigger[i];
    fGenTrigFrac = fracGeneralTrigger[i];

    //2*i + 1 because the format is weird
    fGenTrigSec = secondsGeneralTrigger[2*i + 1];
    unorderedGenTrigTime.push_back( std::make_pair(fGenTrigSec, (fGenTrigCoarse*8. + fGenTrigFrac/512.)) );

    if (fGenTrigFrac == 0.0) break;
    fGenTrigTree->Fill();
  }

  for(size_t i = 0; i < timestampCountGeneralTrigger[0]; ++i){
 // for(size_t i = 0; i < coarseTOF1A.size(); ++i){
 //   std::cout << "TOF1A " << i << " " << secondsTOF1A[2*i+1] << " "  << 8.*coarseTOF1A[i] << " " <<  fracTOF1A[i]/512. << std::endl;
    fXTOF1ACoarse = coarseTOF1A[i];
    fXTOF1AFrac = fracTOF1A[i];
    fXTOF1ASec = secondsTOF1A[2*i + 1];

    if(fXTOF1ASec < secondsTOF1A[1]) break; 

    if (fXTOF1ACoarse == 0.0 && fXTOF1AFrac == 0.0 && fXTOF1ASec == 0.0) break;
    unorderedTOF1ATime.push_back(std::make_pair(fXTOF1ASec, (fXTOF1ACoarse*8. + fXTOF1AFrac/512.)) );
    fXTOF1ATree->Fill();
  }
    
  for(size_t i = 0; i < timestampCountGeneralTrigger[0]; ++i){
 // for(size_t i = 0; i < coarseTOF1B.size(); ++i){
 //   std::cout << "TOF1B " << i << " " << secondsTOF1B[2*i+1] << " "  << 8.*coarseTOF1B[i] << " " <<  fracTOF1B[i]/512. << std::endl;
    fXTOF1BCoarse = coarseTOF1B[i];
    fXTOF1BFrac = fracTOF1B[i];
    fXTOF1BSec = secondsTOF1B[2*i + 1];

    if(fXTOF1BSec < secondsTOF1B[1]) break; 

    if (fXTOF1BCoarse == 0.0 && fXTOF1BFrac == 0.0 && fXTOF1BSec == 0.0) break;
    unorderedTOF1BTime.push_back(std::make_pair(fXTOF1BSec, (fXTOF1BCoarse*8. + fXTOF1BFrac/512.)) );
    fXTOF1BTree->Fill();
  }

  for(size_t i = 0; i < timestampCountGeneralTrigger[0]; ++i){
 // for(size_t i = 0; i < coarseTOF2A.size(); ++i){
 //   std::cout << "TOF2A " << i << " " << secondsTOF2A[2*i+1] << " "  << 8.*coarseTOF2A[i] << " " <<  fracTOF2A[i]/512. << std::endl;
    fXTOF2ACoarse = coarseTOF2A[i];
    fXTOF2AFrac = fracTOF2A[i];
    fXTOF2ASec = secondsTOF2A[2*i + 1];

    if(fXTOF2ASec < secondsTOF2A[1]) break; 

    if (fXTOF2ACoarse == 0.0 && fXTOF2AFrac == 0.0 && fXTOF2ASec == 0.0) break;
    unorderedTOF2ATime.push_back(std::make_pair(fXTOF2ASec, (fXTOF2ACoarse*8. + fXTOF2AFrac/512.)) );
    fXTOF2ATree->Fill();
  }  

  for(size_t i = 0; i < timestampCountGeneralTrigger[0]; ++i){
 // for(size_t i = 0; i < coarseTOF2B.size(); ++i){
 //   std::cout << "TOF2B " << i << " " << secondsTOF2B[2*i+1] << " "  << 8.*coarseTOF2B[i] << " " <<  fracTOF2B[i]/512. << std::endl;
    fXTOF2BCoarse = coarseTOF2B[i];
    fXTOF2BFrac = fracTOF2B[i];
    fXTOF2BSec = secondsTOF2B[2*i + 1];

    if(fXTOF2BSec < secondsTOF2B[1]) break; 

    if (fXTOF2BCoarse == 0.0 && fXTOF2BFrac == 0.0 && fXTOF2BSec == 0.0) break;
    unorderedTOF2BTime.push_back(std::make_pair(fXTOF2BSec, (fXTOF2BCoarse*8. + fXTOF2BFrac/512.) ));
    fXTOF2BTree->Fill();
  }

  for(size_t iT = 0; iT < unorderedGenTrigTime.size(); ++iT){
    
//    bool found_TOF1 = false;
//    bool found_TOF2 = false;

    bool found_TOF = false;

    double the_gen_sec = unorderedGenTrigTime[iT].first;
    double the_gen_ns = unorderedGenTrigTime[iT].second;

    double the_TOF1_sec = -1.;
    double the_TOF2_sec = -1.;
    double the_TOF1_ns = -1.;
    double the_TOF2_ns = -1.;

    //bool TOF1A_passed = false;
    //bool TOF1B_passed = false;
    //bool TOF2A_passed = false;
    //bool TOF2B_passed = false;

    //1A2A = 0; 1B2A = 1, 1A2B = 2,  1B2B = 3
    //Add 1 for 1B, add 2 for 2B
    int channel = 0;

    for(size_t iT1A = 0; iT1A < unorderedTOF1ATime.size(); ++iT1A){
      double TOF1A_sec = unorderedTOF1ATime[iT1A].first;
      double TOF1A_ns = unorderedTOF1ATime[iT1A].second;

      double delta_1A = 1.e9*(the_gen_sec - TOF1A_sec) + the_gen_ns - TOF1A_ns;
      if(delta_1A < 0.){
       // std::cout << "Passed TOF1A" << std::endl;
        //TOF1A_passed = true;
        break;
      }

      if( delta_1A > 500. ) continue;     

      //If here, then 0 < delta_1A < 500ns
      //Check the TOF2 times. 
      
      for(size_t iT2A = 0; iT2A < unorderedTOF2ATime.size(); ++iT2A){
        double TOF2A_sec = unorderedTOF2ATime[iT2A].first;
        double TOF2A_ns = unorderedTOF2ATime[iT2A].second;

        double delta = 1.e9*(TOF2A_sec - TOF1A_sec) + TOF2A_ns - TOF1A_ns;
        
        //Overtaken 1A
        if(delta < 0.){
          continue; 
        }

        if( delta > 500. ){
          break;
        }
        else{

          //If here, then TOF1 is within 500 ns below TOF2

          std::cout << "Found matching TOF2A and TOF1A" << std::endl;
          found_TOF = true;

          the_TOF2_sec = TOF2A_sec;
          the_TOF2_ns  = TOF2A_ns;

          the_TOF1_sec = TOF1A_sec;
          the_TOF1_ns  = TOF1A_ns;

          channel = 0;

          break;
        }       
      }

      if(!found_TOF){
        for(size_t iT2B = 0; iT2B < unorderedTOF2BTime.size(); ++iT2B){
          double TOF2B_sec = unorderedTOF2BTime[iT2B].first;
          double TOF2B_ns = unorderedTOF2BTime[iT2B].second;

          double delta = 1.e9*(TOF2B_sec - TOF1A_sec) + TOF2B_ns - TOF1A_ns;
          
          //Overtaken 1A
          if(delta < 0.){
            continue;
          }

          if( delta > 500. ){
            break;
          }
          else{
            //If here, then TOF1 is within 500 ns below TOF2

            std::cout << "Found matching TOF2B and TOF1A" << std::endl;
            found_TOF = true;

            the_TOF2_sec = TOF2B_sec;
            the_TOF2_ns  = TOF2B_ns;
  
            the_TOF1_sec = TOF1A_sec;
            the_TOF1_ns  = TOF1A_ns;

            channel = 2;

            break;
          }       
        }
      }

      if(found_TOF) break;
    }

    //Now check 1B with 2A and 2B
    if(!found_TOF){
      
      for(size_t iT1B = 0; iT1B < unorderedTOF1BTime.size(); ++iT1B){
        double TOF1B_sec = unorderedTOF1BTime[iT1B].first;
        double TOF1B_ns = unorderedTOF1BTime[iT1B].second;

        double delta_1B = 1.e9*(the_gen_sec - TOF1B_sec) + the_gen_ns - TOF1B_ns;
        if(delta_1B < 0.){
         // std::cout << "Passed TOF1B" << std::endl;
          //TOF1B_passed = true;
          break;
        }

        if( delta_1B > 500. ) continue;     

        //If here, then 0 < delta_1B < 500ns
        //Check the TOF2 times. 
        
        for(size_t iT2A = 0; iT2A < unorderedTOF2ATime.size(); ++iT2A){
          double TOF2A_sec = unorderedTOF2ATime[iT2A].first;
          double TOF2A_ns = unorderedTOF2ATime[iT2A].second;

          double delta = 1.e9*(TOF2A_sec - TOF1B_sec) + TOF2A_ns - TOF1B_ns;
          
          //Overtaken 1B
          if(delta < 0.){
            continue;
          }

          if( delta > 500. ){
            break;
          }
          else{

            //If here, then TOF1 is within 500 ns below TOF2

            std::cout << "Found matching TOF2A and TOF1B" << std::endl;
            found_TOF = true;

            the_TOF2_sec = TOF2A_sec;
            the_TOF2_ns  = TOF2A_ns;

            the_TOF1_sec = TOF1B_sec;
            the_TOF1_ns  = TOF1B_ns;

            channel = 1;

            break;
          }       
        }

        if(!found_TOF){
          for(size_t iT2B = 0; iT2B < unorderedTOF2BTime.size(); ++iT2B){
            double TOF2B_sec = unorderedTOF2BTime[iT2B].first;
            double TOF2B_ns = unorderedTOF2BTime[iT2B].second;

            double delta = 1.e9*(TOF2B_sec - TOF1B_sec) + TOF2B_ns - TOF1B_ns;
            
            //Overtaken 1B
            if(delta < 0.){
              continue;
            }

            if( delta > 500. ){
              break;
            }
            else{
              //If here, then TOF1 is within 500 ns below TOF2

              std::cout << "Found matching TOF2B and TOF1B" << std::endl;
              found_TOF = true;

              the_TOF2_sec = TOF2B_sec;
              the_TOF2_ns  = TOF2B_ns;
  
              the_TOF1_sec = TOF1B_sec;
              the_TOF1_ns  = TOF1B_ns;

              channel = 3;

              break;
            }       
          }
        }
        
        if(found_TOF) break;
      }    
    }

    if(found_TOF){
      //Convert from TAI to UTC at this point

      std::cout << "Adding matched tof" << std::endl;

      beamevt->AddT0(std::make_pair(the_gen_sec - fOffsetTAI, the_gen_ns));
      beamevt->AddTOF0Trigger(std::make_pair(the_TOF1_sec - fOffsetTAI, the_TOF1_ns));
      beamevt->AddTOF1Trigger(std::make_pair(the_TOF2_sec - fOffsetTAI, the_TOF2_ns));
      beamevt->AddTOFChan(channel);        
    }
    else{
      //Add dummy

      std::cout << "Adding unmatched tof" << std::endl;

      beamevt->AddT0(std::make_pair(the_gen_sec - fOffsetTAI, the_gen_ns));
      beamevt->AddTOF0Trigger(std::make_pair(0., 0.));
      beamevt->AddTOF1Trigger(std::make_pair(0., 0.));
      beamevt->AddTOFChan(-1);        
    }

  }

}
// END BeamEvent::parseXTOF
////////////////////////

void proto::BeamEvent::parseXCET(uint64_t time){
  if(fCKov1 != ""){  
    std::cout << "Getting CKov1 info: " << fCKov1 << std::endl;

    std::vector< double > countsTrigCKov1 = FetchWithRetries< std::vector< double > >(time, fXCETPrefix + fCKov1 + ":countsTrig",fNRetries);
    std::vector< double > countsCKov1     = FetchWithRetries< std::vector< double > >(time, fXCETPrefix + fCKov1 + ":counts",fNRetries);
    std::vector< double > pressureCKov1   = FetchWithRetries< std::vector< double > >(time, fXCETPrefix + fCKov1 + ":pressure",fNRetries);
    std::cout << "countsTrig: " << countsTrigCKov1[0] << std::endl; 
    std::cout << "counts: "     << countsCKov1[0] << std::endl;
    std::cout << "pressure: "   << pressureCKov1[0] << std::endl;

    CKov1Pressure = pressureCKov1[0];
    if(countsTrigCKov1[0] > 0){
      CKov1Efficiency = countsCKov1[0]/countsTrigCKov1[0];
    }
  }
  
  if(fCKov2 != ""){
    std::cout << "Getting CKov2 info: " << fCKov2 << std::endl;

    std::vector< double > countsTrigCKov2 = FetchWithRetries< std::vector< double > >(time, fXCETPrefix + fCKov2 + ":countsTrig",fNRetries);
    std::vector< double > countsCKov2     = FetchWithRetries< std::vector< double > >(time, fXCETPrefix + fCKov2 + ":counts",fNRetries);
    std::vector< double > pressureCKov2   = FetchWithRetries< std::vector< double > >(time, fXCETPrefix + fCKov2 + ":pressure",fNRetries);
    std::cout << "countsTrig: " << countsTrigCKov2[0] << std::endl; 
    std::cout << "counts: "     << countsCKov2[0] << std::endl;
    std::cout << "pressure: "   << pressureCKov2[0] << std::endl;

    CKov2Pressure = pressureCKov2[0];
    if(countsTrigCKov2[0] > 0){
      CKov2Efficiency = countsCKov2[0]/countsTrigCKov2[0];
    }
  }

  beam::CKov CKov1Status, CKov2Status;

//  double triggerTime = beamevt->GetT0( beamevt->GetActiveTrigger() ).first;
//  triggerTime += 1.e-9*beamevt->GetT0( beamevt->GetActiveTrigger() ).second;

//  CKov1Status.timeStamp = triggerTime;
  CKov1Status.pressure  = CKov1Pressure;
  CKov1Status.trigger   = C1;
  beamevt->SetCKov0( CKov1Status );

//  CKov2Status.timeStamp = triggerTime;
  CKov2Status.pressure  = CKov2Pressure;
  CKov2Status.trigger   = C2;
  beamevt->SetCKov1( CKov2Status );

}
// END BeamEvent::parseXCET
////////////////////////



////////////////////////
// 
void proto::BeamEvent::parseGeneralXBPF(std::string name, uint64_t time, size_t ID){

  // Retrieve the number of counts in the BPF
  std::vector<double> counts;
  counts = FetchWithRetries< std::vector<double> >(time, fXBPFPrefix + name + ":countsRecords[]",fNRetries);
  std::cout << "Counts: " << counts.size() << std::endl;
  for(size_t i = 0; i < counts.size(); ++i){
    std::cout << counts[i] << std::endl;
  }

  std::vector<double> data;
  data = FetchWithRetries< std::vector<double> >(time, fXBPFPrefix + name + ":eventsData[]",fNRetries);
  std::cout << "Data: " << data.size() << std::endl;

  // If the number of counts is larger than the number of general triggers
  // make note
  if(counts[1] > beamevt->GetNT0()){
    std::cout << "WARNING MISMATCH " << counts[1] << " " << beamevt->GetNT0() << std::endl;
  }
  
  beam::FBM fbm;
  fbm.ID = ID;
    
  //Use this just in case any are out of sync?
  //Shouldn't be, but just to be safe...
  //Helps cut down on time
  std::vector<size_t> leftOvers;      
  for(size_t lo = 0; lo < beamevt->GetNT0(); ++lo){
    leftOvers.push_back(lo);
  }
 
  for(size_t i = 0; i < counts[1]; ++i){      
//      std::cout << "Count: " << i << std::endl;
    
    for(int j = 0; j < 10; ++j){
      double theData = data[20*i + (2*j + 1)];
//      std::cout << std::setw(15) << theData ;
      if(j < 4){
	fbm.timeData[j] = theData;           
      }else{
	fbm.fiberData[j - 4] = theData;
      } 
    } 

    // Check the time data for corruption
    if(fbm.timeData[1] < .0000001){
//      std::cout << "Skipping bad time" << std::endl;
      continue;
    } 

    // Timestamp is in units of sec + 8ns ticks
    fbm.timeStamp = fbm.timeData[3] + fbm.timeData[2]*8.e-9;  
    
    //Go through the valid Good Particles, and emplace the FBM 
//    std::cout << "Checking " << beamevt->GetNT0() << " triggers " << leftOvers.size() << std::endl;

    for(std::vector<size_t>::iterator ip = leftOvers.begin(); ip != leftOvers.end(); ++ip){
//       std::cout.precision(dbl::max_digits10);
//       std::cout << "\t" << fbm.timeStamp  - fOffsetTAI<< " " << beamevt->GetFullT0(*ip) << std::endl;
//       std::cout.precision(dbl::max_digits10);
//       std::cout << "\t" << beamevt->GetFullT0(*ip) - (fbm.timeStamp - fOffsetTAI)  << std::endl;

      // Compute the time delta between the timeStamp and the T0, see if it's less than 500ns away
//      std::cout << fbm.timeStamp << " " << beamevt->GetFullT0(*ip) << std::endl;
      if( fabs(beamevt->GetFullT0(*ip) - (fbm.timeStamp - fOffsetTAI) ) < 1000.e-9
       /*&& (beamevt->GetFullT0(*ip) - (fbm.timeStamp - fOffsetTAI) ) > 0.*/     ){

	if(beamevt->GetFBM(name, *ip).ID != -1){
	  std::cout << "Warning: Replacing non-dummy FBM at "
		    << name << " " << *ip << std::endl;
	} 
	
//	 std::cout << "Replacing at timestamp " << fbm.timeStamp << std::endl;
	beamevt->ReplaceFBMTrigger(name, fbm, *ip);
	leftOvers.erase(ip);
	break;
      } 
    } 

//    std::cout << std::endl;
  } 

  for(size_t i = 0; i < beamevt->GetNFBMTriggers(name); ++i){
    beamevt->DecodeFibers(name,i);
//    std::cout << name << " at time: "
//	      << beamevt->DecodeFiberTime(name, i) << " has active fibers: ";
//    for(size_t iFiber = 0; iFiber < beamevt->GetActiveFibers(name,i).size(); ++iFiber){
//      std::cout << beamevt->GetActiveFibers(name, i)[iFiber] << " ";
//    }
//    std::cout << std::endl;
 /*   
    *fActiveFibers[name] = beamevt->GetActiveFibers(name,i);
    fProfTime[name] = beamevt->DecodeFiberTime(name, i, fOffsetTAI);
    std::cout << beamevt->ReturnTriggerAndTime(name,i)[0] << " "
	      << beamevt->ReturnTriggerAndTime(name,i)[1] << " "
	      << beamevt->ReturnTriggerAndTime(name,i)[2] << " "
	      << beamevt->ReturnTriggerAndTime(name,i)[3] << std::endl;

    fProfTrigger1[name] = beamevt->ReturnTriggerAndTime(name,i)[0];
    fProfTrigger2[name] = beamevt->ReturnTriggerAndTime(name,i)[1];
    fProfTime1[name] = beamevt->ReturnTriggerAndTime(name,i)[2];
    fProfTime2[name] = beamevt->ReturnTriggerAndTime(name,i)[3];
    fProfTree[name]->Fill(); 
*/
    
  } 

}

////////////////////////
// 
void proto::BeamEvent::parseXBPF(uint64_t time){
  for(size_t d = 0; d < fDevices.size(); ++d){
    std::string name = fDevices[d];
    std::cout <<"Device: " << name << std::endl;
        parseGeneralXBPF(name, time, d);
  }  
}
// END BeamEvent::parseXBFP
////////////////////////

////////////////////////
// 
void proto::BeamEvent::parsePairedXBPF(uint64_t time){
  for(size_t d = 0; d < fPairedDevices.size(); ++d){
    std::string name = fPairedDevices[d].first;
    std::cout <<"Device: " << name << std::endl;
      parseGeneralXBPF(name, time, d);

    name = fPairedDevices[d].second;
    std::cout <<"Device: " << name << std::endl;
        parseGeneralXBPF(name, time, d);
  }
}
// END BeamEvent::parsePairedXBFP
////////////////////////

////////////////////////
// 
void proto::BeamEvent::parsePairedStraightXBPF(uint64_t time){
  for(size_t d = 0; d < fPairedStraightDevices.size(); ++d){
    std::string name = fPairedStraightDevices[d].first;
    std::cout <<"Device: " << name << std::endl;
        parseGeneralXBPF(name, time, d);

    name = fPairedStraightDevices[d].second;
    std::cout <<"Device: " << name << std::endl;
        parseGeneralXBPF(name, time, d);
  }
}
// END BeamEvent::parsePairedStrightXBFP
////////////////////////

void proto::BeamEvent::beginJob()
{
  art::ServiceHandle<art::TFileService> tfs;
  

  fTOFHist = tfs->make<TH1F>("TOF","",500,0,500);
  fCKovHist = tfs->make<TH1F>("CKov","",4,0,4);
  fMomentumHist = tfs->make<TH1F>("Momentum", "", 150, 0, 15.);

  fDeltaXHist  = tfs->make<TH1F>("DeltaX",  "",600,-300,300);
  fDeltaYHist  = tfs->make<TH1F>("DeltaY",  "",600,-300,300);
  fDeltaZHist  = tfs->make<TH1F>("DeltaZ",  "",600,-300,300);
  fDeltaXYHist = tfs->make<TH2F>("DeltaXY", "",600,-300,300,600,-300,300);
  fDeltaYZHist = tfs->make<TH2F>("DeltaYZ", "",600,-300,300,600,-300,300);
  fDeltaZXHist = tfs->make<TH2F>("DeltaZX", "",600,-300,300,600,-300,300);
  fDeltaTree   = tfs->make<TTree>("Delta", "");
  fDeltaX = 0.;
  fDeltaY = 0.;
  fDeltaZ = 0.;
  trackNum = 0;
  beamNum  = 0;
  fDeltaTree->Branch("DeltaX", &fDeltaX);
  fDeltaTree->Branch("DeltaY", &fDeltaY);
  fDeltaTree->Branch("DeltaZ", &fDeltaZ);
  fDeltaTree->Branch("BeamX",  &fBeamPosX);
  fDeltaTree->Branch("BeamY",  &fBeamPosY);
  fDeltaTree->Branch("BeamZ",  &fBeamPosZ);
  fDeltaTree->Branch("TrackX", &fTrackPosX);
  fDeltaTree->Branch("TrackY", &fTrackPosY);
  fDeltaTree->Branch("TrackZ", &fTrackPosZ);
  fDeltaTree->Branch("Event",  &eventNum);
  fDeltaTree->Branch("Track",  &trackNum);
  fDeltaTree->Branch("Beam",   &beamNum);


  fFullMomentum = tfs->make<TH1F>("FullMomentum", "", 150, 0, 15.);
  fCutMomentum = tfs->make<TH1F>("CutMomentum", "", 150, 0, 15.);

  fOutTree = tfs->make<TTree>("tree", "lines"); 
  theTrack = new recob::Track();
  eventTime = 0;
  eventNum = 0;
  runNum = 0;
  subRunNum = 0;
  CKov1Pressure = 0.;
  CKov2Pressure = 0.;
  CKov1Efficiency = 0.;
  CKov2Efficiency = 0.;
  fOutTree->Branch("Track", &theTrack);
  fOutTree->Branch("Time", &eventTime);
  fOutTree->Branch("RDTS", &RDTSTime);
  GenTriggers = new std::vector<double>;
  fOutTree->Branch("Gens", &GenTriggers);
  fOutTree->Branch("Event", &eventNum);
  fOutTree->Branch("SpillStart", &SpillStart);
  fOutTree->Branch("SpillEnd", &SpillEnd);
  fOutTree->Branch("SpillOffset", &SpillOffset);
  fOutTree->Branch("ActiveTriggerTime", &ActiveTriggerTime);
  fOutTree->Branch("Run",   &runNum);
  fOutTree->Branch("Subrun", &subRunNum);
  fOutTree->Branch("Pressure1", &CKov1Pressure);
  fOutTree->Branch("Eff1", &CKov1Efficiency);
  fOutTree->Branch("Pressure2", &CKov2Pressure);
  fOutTree->Branch("Eff2", &CKov2Efficiency);
  
  acqTime = 0;
  fOutTree->Branch("acqTime",        &acqTime);
  fOutTree->Branch("acqStampMBPL",   &acqStampMBPL);
  fOutTree->Branch("cycleStampMBPL", &cycleStampMBPL);
  HLTWord = 0;
  HLTTS = 0;
  fOutTree->Branch("HLTWord",     &HLTWord);
  fOutTree->Branch("HLTTS",       &HLTTS);
  fOutTree->Branch("BeamOn",      &BeamOn);
  fOutTree->Branch("BITrigger",   &BITrigger);
  fOutTree->Branch("Upstream",    &Upstream);
  fOutTree->Branch("C1",          &C1);
  fOutTree->Branch("C2",          &C2);
  fOutTree->Branch("BP1",         &BP1);
  fOutTree->Branch("BP2",         &BP2);
  fOutTree->Branch("BP3",         &BP3);
  fOutTree->Branch("BP4",         &BP4);

//  fOutTree->Branch("fetchTime", &fetchTime);

  fTrackTree = tfs->make<TTree>("tracks","");
  trackX = new std::vector<double>;
  trackY = new std::vector<double>;
  trackZ = new std::vector<double>;
  fTrackTree->Branch("X", &trackX);
  fTrackTree->Branch("Y", &trackY);
  fTrackTree->Branch("Z", &trackZ);

  fMatchedTriggers = tfs->make<TTree>("matched","");
  
  fMatchedTriggers->Branch("Gen", &matchedGen);
  fMatchedTriggers->Branch("TOF1", &matchedTOF1);
  fMatchedTriggers->Branch("TOF2", &matchedTOF2);
  fMatchedTriggers->Branch("Chan", &matchedChan);
  fMatchedTriggers->Branch("NominalTriggers", &matchedNom);

  for(size_t i = 0; i < fPairedDevices.size(); ++i){
    std::string name = "BeamProf2D_" + fPairedDevices[i].first + "_" + fPairedDevices[i].second;
    std::string title = fPairedDevices[i].first + ", " + fPairedDevices[i].second;
    fBeamProf2D.push_back( tfs->make<TH2F>(name.c_str(),title.c_str(),192,0,192,192,0,192) );

    name = "BeamProf1D_" + fPairedDevices[i].first;
    title = fPairedDevices[i].first;
    fBeamProf1D[fPairedDevices[i].first] = ( tfs->make<TH1F>(name.c_str(),title.c_str(),192,0,192) );

    name = "BeamProf1D_" + fPairedDevices[i].second;
    title = fPairedDevices[i].second;
    fBeamProf1D[fPairedDevices[i].second] = ( tfs->make<TH1F>(name.c_str(),title.c_str(),192,0,192) );
    
    name = "Fibers_" + fPairedDevices[i].first;
    fActiveFibers[fPairedDevices[i].first] = new std::vector<short>;

    fProfTime[fPairedDevices[i].first] = 0.;
    fProfTrigger1[fPairedDevices[i].first] = 0.;
    fProfTrigger2[fPairedDevices[i].first] = 0.;
    fProfTime1[fPairedDevices[i].first] = 0.;
    fProfTime2[fPairedDevices[i].first] = 0.;

    matchedXBPF[fPairedDevices[i].first] = 0.;
    fMatchedTriggers->Branch((fPairedDevices[i].first).c_str(), &matchedXBPF[fPairedDevices[i].first]);


    fProfTree[fPairedDevices[i].first] = ( tfs->make<TTree>(name.c_str(), "XBPF") );
    fProfTree[fPairedDevices[i].first]->Branch("time", &fProfTime[fPairedDevices[i].first]);
    fProfTree[fPairedDevices[i].first]->Branch("fibers", &fActiveFibers[fPairedDevices[i].first]);
    fProfTree[fPairedDevices[i].first]->Branch("trigger_1", &fProfTrigger1[fPairedDevices[i].first]);
    fProfTree[fPairedDevices[i].first]->Branch("trigger_2", &fProfTrigger2[fPairedDevices[i].first]);
    fProfTree[fPairedDevices[i].first]->Branch("time_1", &fProfTime1[fPairedDevices[i].first]);
    fProfTree[fPairedDevices[i].first]->Branch("time_2", &fProfTime2[fPairedDevices[i].first]);

    name = "Fibers_" + fPairedDevices[i].second;
    fActiveFibers[fPairedDevices[i].second] = new std::vector<short>;

    fProfTime[fPairedDevices[i].second] = 0.;
    fProfTrigger1[fPairedDevices[i].second] = 0.;
    fProfTrigger2[fPairedDevices[i].second] = 0.;
    fProfTime1[fPairedDevices[i].second] = 0.;
    fProfTime2[fPairedDevices[i].second] = 0.;

    matchedXBPF[fPairedDevices[i].second] = 0.;
    fMatchedTriggers->Branch((fPairedDevices[i].second).c_str(), &matchedXBPF[fPairedDevices[i].second]);

    fProfTree[fPairedDevices[i].second] = ( tfs->make<TTree>(name.c_str(), "XBPF") );
    fProfTree[fPairedDevices[i].second]->Branch("time", &fProfTime[fPairedDevices[i].second]);
    fProfTree[fPairedDevices[i].second]->Branch("fibers", &fActiveFibers[fPairedDevices[i].second]);
    fProfTree[fPairedDevices[i].second]->Branch("trigger_1", &fProfTrigger1[fPairedDevices[i].second]);
    fProfTree[fPairedDevices[i].second]->Branch("trigger_2", &fProfTrigger2[fPairedDevices[i].second]);
    fProfTree[fPairedDevices[i].second]->Branch("time_1", &fProfTime1[fPairedDevices[i].second]);
    fProfTree[fPairedDevices[i].second]->Branch("time_2", &fProfTime2[fPairedDevices[i].second]);
  }

  for(size_t i = 0; i < fPairedStraightDevices.size(); ++i){
    std::string name = "BeamProf2D_" + fPairedStraightDevices[i].first + "_" + fPairedStraightDevices[i].second;
    std::string title = fPairedStraightDevices[i].first + ", " + fPairedStraightDevices[i].second;
    fBeamProf2D.push_back( tfs->make<TH2F>(name.c_str(),title.c_str(),192,0,192,192,0,192) );

    name = "BeamProf1D_" + fPairedStraightDevices[i].first;
    title = fPairedStraightDevices[i].first;
    fBeamProf1D[fPairedStraightDevices[i].first] = ( tfs->make<TH1F>(name.c_str(),title.c_str(),192,0,192) );

    name = "BeamProf1D_" + fPairedStraightDevices[i].second;
    title = fPairedStraightDevices[i].second;
    fBeamProf1D[fPairedStraightDevices[i].second] = ( tfs->make<TH1F>(name.c_str(),title.c_str(),192,0,192) );
    
    name = "Fibers_" + fPairedStraightDevices[i].first;
    fActiveFibers[fPairedStraightDevices[i].first] = new std::vector<short>;

    fProfTime[fPairedStraightDevices[i].first] = 0.;
    fProfTrigger1[fPairedStraightDevices[i].first] = 0.;
    fProfTrigger2[fPairedStraightDevices[i].first] = 0.;
    fProfTime1[fPairedStraightDevices[i].first] = 0.;
    fProfTime2[fPairedStraightDevices[i].first] = 0.;

    matchedXBPF[fPairedStraightDevices[i].first] = 0.;
    fMatchedTriggers->Branch((fPairedStraightDevices[i].first).c_str(), &matchedXBPF[fPairedStraightDevices[i].first]);

    fProfTree[fPairedStraightDevices[i].first] = ( tfs->make<TTree>(name.c_str(), "XBPF") );
    fProfTree[fPairedStraightDevices[i].first]->Branch("time", &fProfTime[fPairedStraightDevices[i].first]);
    fProfTree[fPairedStraightDevices[i].first]->Branch("fibers", &fActiveFibers[fPairedStraightDevices[i].first]);
    fProfTree[fPairedStraightDevices[i].first]->Branch("trigger_1", &fProfTrigger1[fPairedStraightDevices[i].first]);
    fProfTree[fPairedStraightDevices[i].first]->Branch("trigger_2", &fProfTrigger2[fPairedStraightDevices[i].first]);
    fProfTree[fPairedStraightDevices[i].first]->Branch("time_1", &fProfTime1[fPairedStraightDevices[i].first]);
    fProfTree[fPairedStraightDevices[i].first]->Branch("time_2", &fProfTime2[fPairedStraightDevices[i].first]);

    name = "Fibers_" + fPairedStraightDevices[i].second;
    fActiveFibers[fPairedStraightDevices[i].second] = new std::vector<short>;

    fProfTime[fPairedStraightDevices[i].second] = 0.;
    fProfTrigger1[fPairedStraightDevices[i].second] = 0.;
    fProfTrigger2[fPairedStraightDevices[i].second] = 0.;
    fProfTime1[fPairedStraightDevices[i].second] = 0.;
    fProfTime2[fPairedStraightDevices[i].second] = 0.;

    matchedXBPF[fPairedStraightDevices[i].second] = 0.;
    fMatchedTriggers->Branch((fPairedStraightDevices[i].second).c_str(), &matchedXBPF[fPairedStraightDevices[i].second]);

    fProfTree[fPairedStraightDevices[i].second] = ( tfs->make<TTree>(name.c_str(), "XBPF") );
    fProfTree[fPairedStraightDevices[i].second]->Branch("time", &fProfTime[fPairedStraightDevices[i].second]);
    fProfTree[fPairedStraightDevices[i].second]->Branch("fibers", &fActiveFibers[fPairedStraightDevices[i].second]);
    fProfTree[fPairedStraightDevices[i].second]->Branch("trigger_1", &fProfTrigger1[fPairedStraightDevices[i].second]);
    fProfTree[fPairedStraightDevices[i].second]->Branch("trigger_2", &fProfTrigger2[fPairedStraightDevices[i].second]);
    fProfTree[fPairedStraightDevices[i].second]->Branch("time_1", &fProfTime1[fPairedStraightDevices[i].second]);
    fProfTree[fPairedStraightDevices[i].second]->Branch("time_2", &fProfTime2[fPairedStraightDevices[i].second]);
  }

  for(size_t i = 0; i < fDevices.size(); ++i){
    std::string name = "BeamProf1D_" + fDevices[i];
    std::string title = fDevices[i];
    fBeamProf1D[fDevices[i]] = ( tfs->make<TH1F>(name.c_str(),title.c_str(),192,0,192) );

    name = "Fibers_" + fDevices[i];
    fProfTree[fDevices[i]] = ( tfs->make<TTree>(name.c_str(), "XBPF") );

    fProfTrigger1[fDevices[i]] = 0.;
    fProfTrigger2[fDevices[i]] = 0.;
    fProfTime1[fDevices[i]] = 0.;
    fProfTime2[fDevices[i]] = 0.;
    fProfTime[fDevices[i]] = 0.;

    matchedXBPF[fDevices[i]] = 0.;
    fMatchedTriggers->Branch((fDevices[i]).c_str(), &matchedXBPF[fDevices[i]]);

    fActiveFibers[fDevices[i]] = new std::vector<short>;
    fProfTree[fDevices[i]]->Branch("time", &fProfTime[fDevices[i]]);
    fProfTree[fDevices[i]]->Branch("fibers", &fActiveFibers[fDevices[i]]);
    fProfTree[fDevices[i]]->Branch("trigger_1", &fProfTrigger1[fDevices[i]]);
    fProfTree[fDevices[i]]->Branch("trigger_2", &fProfTrigger2[fDevices[i]]);
    fProfTree[fDevices[i]]->Branch("time_1", &fProfTime1[fDevices[i]]);
    fProfTree[fDevices[i]]->Branch("time_2", &fProfTime2[fDevices[i]]);
  }

  fGenTrigTree = tfs->make<TTree>("GenTrig","");
  fGenTrigTree->Branch("coarse", &fGenTrigCoarse);
  fGenTrigTree->Branch("frac", &fGenTrigFrac);
  fGenTrigTree->Branch("sec", &fGenTrigSec);

  fXTOF1ATree = tfs->make<TTree>("TOF1A","");
  fXTOF1ATree->Branch("coarse", &fXTOF1ACoarse);
  fXTOF1ATree->Branch("frac", &fXTOF1AFrac);
  fXTOF1ATree->Branch("sec", &fXTOF1ASec);

  fXTOF1BTree = tfs->make<TTree>("TOF1B","");
  fXTOF1BTree->Branch("coarse", &fXTOF1BCoarse);
  fXTOF1BTree->Branch("frac", &fXTOF1BFrac);
  fXTOF1BTree->Branch("sec", &fXTOF1BSec);

  fXTOF2ATree = tfs->make<TTree>("TOF2A","");
  fXTOF2ATree->Branch("coarse", &fXTOF2ACoarse);
  fXTOF2ATree->Branch("frac", &fXTOF2AFrac);
  fXTOF2ATree->Branch("sec", &fXTOF2ASec);

  fXTOF2BTree = tfs->make<TTree>("TOF2B","");
  fXTOF2BTree->Branch("coarse", &fXTOF2BCoarse);
  fXTOF2BTree->Branch("frac", &fXTOF2BFrac);
  fXTOF2BTree->Branch("sec", &fXTOF2BSec);

  fXBTF022702Tree = tfs->make<TTree>("XBTF022702","");
  fXBTF022702Tree->Branch("coarse", &fXBTF022702Coarse);
  fXBTF022702Tree->Branch("frac", &fXBTF022702Frac);


}

void proto::BeamEvent::beginRun(art::Run & r)
{
  // Implementation of optional member function here.
}

void proto::BeamEvent::beginSubRun(art::SubRun & sr)
{
  // Implementation of optional member function here.
}

void proto::BeamEvent::endJob()
{
  // Implementation of optional member function here.
}

void proto::BeamEvent::endRun(art::Run & r)
{
  // Implementation of optional member function here.
}

void proto::BeamEvent::endSubRun(art::SubRun & sr)
{
  // Implementation of optional member function here.
}

void proto::BeamEvent::reconfigure(fhicl::ParameterSet const & p)
{
  // Implementation of optional member function here.
  fBundleName  = p.get<std::string>("BundleName");
  fOutputLabel = p.get<std::string>("OutputLabel");
  fInputLabel  = p.get<std::string>("InputLabel");
  fURLStr      = p.get<std::string>("URLStr");
  fNRetries    = p.get<int>("NRetries");
  fValidWindow = p.get<double>("ValidWindow");
  fTimeWindow  = p.get<double>("TimeWindow");
  fFixedTime   = p.get<uint64_t>("FixedTime");
  fMultipleTimes = p.get< std::vector<uint64_t> >("MultipleTimes");
  fTolerance   = p.get<double>("Tolerance");

  fDevices     = p.get< std::vector< std::string > >("Devices");    
  fPairedDevices = p.get< std::vector< std::pair<std::string, std::string> > >("PairedDevices");
  fPairedStraightDevices = p.get< std::vector< std::pair<std::string, std::string> > >("PairedStraightDevices");   

  //For Tracking/////
  firstUpstreamName    = p.get< std::string >("FirstUpstream");
  secondUpstreamName   = p.get< std::string >("SecondUpstream");
  firstDownstreamName  = p.get< std::string >("FirstDownstream");
  secondDownstreamName = p.get< std::string >("SecondDownstream");
  ///////////////////

  //For Momentum Spectrometry////
  firstBPROF1  = p.get< std::string >("FirstBPROF1");
  secondBPROF1 = p.get< std::string >("SecondBPROF1");
  BPROF2       = p.get< std::string >("BPROF2");
  BPROF3       = p.get< std::string >("BPROF3");
  fBeamBend    = p.get< double >("BeamBend");
/*  L1           = 2.004;//(m)
  L2           = 1.718*cos(fBeamBend);//(m)
  L3           = 2.728*cos(fBeamBend);//(m)
*/
  magnetLen    = 1.;//(m)
  magnetField  = 1000.;//()
  ///////////////////////////////

  std::vector< std::pair<std::string, std::string> >  tempTypes = p.get<std::vector< std::pair<std::string, std::string> >>("DeviceTypes");
  fDeviceTypes  = std::map<std::string, std::string>(tempTypes.begin(), tempTypes.end() );

  //Location of Device 
  std::vector< std::pair<std::string, std::array<double,3> > > tempCoords = p.get<std::vector< std::pair<std::string, std::array<double,3> > > >("Coordinates");
  fCoordinates = std::map<std::string, std::array<double,3> >(tempCoords.begin(), tempCoords.end());
//  fCoordinates = p.get< std::vector< std::array<double,3> > >("Coordinates");

  //Rotation of Device
  std::vector< std::pair<std::string, std::array<double,3> > > tempRots = p.get< std::vector< std::pair<std::string, std::array<double,3> > > >("Rotations"); 
  fRotations = std::map<std::string, std::array<double,3> >(tempRots.begin(), tempRots.end());
//  fRotations   = p.get< std::vector< std::array<double,3> > >("Rotations"); 

  //Deminsion of Fibers
  std::vector< std::pair<std::string, double> > tempFiberDims = p.get< std::vector<std::pair<std::string,double> > >("Dimension");
  fFiberDimension = std::map<std::string, double>(tempFiberDims.begin(), tempFiberDims.end());
//  fFiberDimension = p.get< std::vector<double> >("Dimension");

  std::array<double,3> detCoords = p.get<std::array<double,3>>("GlobalDetCoords");
  fGlobalDetCoords = TVector3(detCoords[0],detCoords[1],detCoords[2]);

  fDetRotation = p.get<std::array<double,3>>("DetRotation");
  
  //XTOF devices 
  fTOF1 = p.get< std::string >("TOF1");
  fTOF2 = p.get< std::string >("TOF2");

  fTOF1A = fTOF1 + "A";
  fTOF1B = fTOF1 + "B";
  fTOF2A = fTOF2 + "A";
  fTOF2B = fTOF2 + "B";

  fCKov1 = p.get< std::string >("CKov1");
  fCKov2 = p.get< std::string >("CKov2");


  fXBPFPrefix      = p.get<std::string>("XBPFPrefix");
  fXTOFPrefix      = p.get<std::string>("XTOFPrefix");
  fXCETPrefix      = p.get<std::string>("XCETPrefix");
  fDummyEventTime = p.get<double>("DummyEventTime");


  //New parameters to match Leigh's
  fRotateMonitorXZ = p.get<double>("RotateMonitorXZ");
  fRotateMonitorYZ = p.get<double>("RotateMonitorYZ");

  fFirstTrackingProfZ  = p.get<double>("FirstTrackingProfZ");
  fSecondTrackingProfZ = p.get<double>("SecondTrackingProfZ");
  fNP04FrontZ          = p.get<double>("NP04FrontZ");  

  fBeamX               = p.get<double>("BeamX");
  fBeamY               = p.get<double>("BeamY");
  fBeamZ               = p.get<double>("BeamZ");

  fForceNewFetch       = p.get<bool>("ForceNewFetch");
  fMatchTime           = p.get<bool>("MatchTime");
  fForceRead           = p.get<bool>("ForceRead");


  fTimingCalibration      = p.get<double>("TimingCalibration");
  fCalibrationTolerance   = p.get<double>("CalibrationTolerance");
  fOffsetTAI              = p.get<double>("OffsetTAI");

}

uint64_t proto::BeamEvent::joinHighLow(double high, double low){

  std::cout << "%%% Joining high and low %%%" << std::endl;

  uint32_t low32 = (uint32_t)low;
  std::bitset<64> lowbits = low32;

  uint32_t high32 = (uint32_t)high;
  std::bitset<64> highbits = high32;

  highbits = highbits << 32;
  std::bitset<64> joinedbits = highbits ^ lowbits;

  std::cout << low << " " << low32 << std::endl;
  std::cout << lowbits << std::endl;

  std::cout << high << " " << high32 << std::endl;
  std::cout << highbits << std::endl;

  std::cout << joinedbits.to_ullong() << std::endl;
  std::cout << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << std::endl << std::endl;

  return joinedbits.to_ullong(); 
}

/*std::bitset<sizeof(long)*CHAR_BIT> proto::BeamEvent::toBinary(long num){
   
  std::bitset<sizeof(double)*CHAR_BIT> mybits(num);  
  std::bitset<32> upper, lower;
  for(int i = 0; i < 32; ++i){
    lower[i] = mybits[i];
    upper[i] = mybits[i + 32];   
  }
  if(upper.any()) std::cout << "WARNING: NONZERO HALF" << std::endl;

  return mybits;
}
*/

TVector3 proto::BeamEvent::ConvertProfCoordinates(double x, double y, double z, double zOffset){
  double off = fNP04FrontZ - zOffset;

  TVector3 old(x,y,z);

  double newX = x*fBMBasisX.X() + y*fBMBasisY.X() + /*(z-zOffset)*fBMBasisZ.X()*/ + off*fabs(fBMBasisZ.X());
  double newY = x*fBMBasisX.Y() + y*fBMBasisY.Y() + /*(z-zOffset)*fBMBasisZ.Y()*/ + off*fabs(fBMBasisZ.Y());
  double newZ = x*fBMBasisX.Z() + y*fBMBasisY.Z() + /*(z-zOffset)              */ - off*fabs(fBMBasisZ.Z());

  newX += fBeamX*10.;
  newY += fBeamY*10.;
  newZ += fBeamZ*10.;

  TVector3 result(newX/10., newY/10., newZ/10.);
  return result;
}

void proto::BeamEvent::BeamMonitorBasisVectors(){
  std::cout << "Rotating" << std::endl;
  fBMBasisX = TVector3(1.,0.,0.);
  fBMBasisY = TVector3(0.,1.,0.);
  fBMBasisZ = TVector3(0.,0.,1.);
  RotateMonitorVector(fBMBasisX);
  std::cout << fBMBasisX.X() << " " << fBMBasisX.Y() << " " << fBMBasisX.Z() << std::endl;
  RotateMonitorVector(fBMBasisY);
  std::cout << fBMBasisY.X() << " " << fBMBasisY.Y() << " " << fBMBasisY.Z() << std::endl;
  RotateMonitorVector(fBMBasisZ);
  std::cout << fBMBasisZ.X() << " " << fBMBasisZ.Y() << " " << fBMBasisZ.Z() << std::endl;
}

void proto::BeamEvent::RotateMonitorVector(TVector3 &vec){
  vec.RotateY(fRotateMonitorXZ * TMath::Pi()/180.);
  vec.RotateX(fRotateMonitorYZ * TMath::Pi()/180.);
}

void proto::BeamEvent::MakeTrack(size_t theTrigger){
  
  std::cout << "Making Track for time: " << beamevt->GetFullT0(theTrigger) << std::endl;

  //Get the active fibers from the upstream tracking XBPF
  std::vector<short> firstUpstreamFibers  = beamevt->GetActiveFibers(firstUpstreamName, theTrigger);
  std::vector<short> secondUpstreamFibers = beamevt->GetActiveFibers(secondUpstreamName, theTrigger);

  std::cout << firstUpstreamName << " has " << firstUpstreamFibers.size() << " active fibers at time " << beamevt->GetFiberTime(firstUpstreamName,theTrigger) << std::endl;
  for(size_t i = 0; i < firstUpstreamFibers.size(); ++i){
    std::cout << firstUpstreamFibers[i] << " ";
  }
  std::cout << std::endl;

  std::cout << secondUpstreamName << " has " << secondUpstreamFibers.size() << " active fibers at time " << beamevt->GetFiberTime(secondUpstreamName,theTrigger) << std::endl;
  for(size_t i = 0; i < secondUpstreamFibers.size(); ++i){
    std::cout << secondUpstreamFibers[i] << " ";
  }
  std::cout << std::endl;
  //////////////////////////////////////////////

  //Get the active fibers from the downstream tracking XBPF
  std::vector<short> firstDownstreamFibers = beamevt->GetActiveFibers(firstDownstreamName, theTrigger);
  std::vector<short> secondDownstreamFibers = beamevt->GetActiveFibers(secondDownstreamName, theTrigger);

  std::cout << firstDownstreamName << " has " << firstDownstreamFibers.size() << " active fibers at time " << beamevt->GetFiberTime(firstDownstreamName,theTrigger) << std::endl;
  for(size_t i = 0; i < firstDownstreamFibers.size(); ++i){
    std::cout << firstDownstreamFibers[i] << " ";
  }
  std::cout << std::endl;

  std::cout << secondDownstreamName << " has " << secondDownstreamFibers.size() << " active fibers at time " << beamevt->GetFiberTime(secondDownstreamName,theTrigger) << std::endl;
  for(size_t i = 0; i < secondDownstreamFibers.size(); ++i){
    std::cout << secondDownstreamFibers[i] << " ";
  }
  std::cout << std::endl;
  //////////////////////////////////////////////

  if( (firstUpstreamFibers.size() < 1) || (secondUpstreamFibers.size() < 1) || (firstDownstreamFibers.size() < 1) || (secondDownstreamFibers.size() < 1) ){
    std::cout << "Warning, at least one empty Beam Profiler. Not making track" << std::endl;
    return;
  }
  else if( (firstUpstreamFibers.size() > 5) || (secondUpstreamFibers.size() > 5) || (firstDownstreamFibers.size() > 5) || (secondDownstreamFibers.size() > 5) ){
    std::cout << "Warning, too many (>5) active fibers in at least one Beam Profiler. Not making track" << std::endl;
    return;
  }

  //We have the active Fibers, now go through them.
  //Skip the second of any adjacents 
  std::vector< std::pair<size_t, size_t> > upstreamPairedFibers;
  std::vector< std::pair<size_t, size_t> > downstreamPairedFibers;
  std::string firstUpstreamType    = fDeviceTypes[firstUpstreamName]; 
  std::string secondUpstreamType   = fDeviceTypes[secondUpstreamName]; 
  std::string firstDownstreamType  = fDeviceTypes[firstDownstreamName]; 
  std::string secondDownstreamType = fDeviceTypes[secondDownstreamName]; 

  std::vector< TVector3 > upstreamPositions;
  std::vector< TVector3 > downstreamPositions; 
  
  //Pair the upstream fibers together
  std::cout << "Upstream" << std::endl;
  for(size_t iF1 = 0; iF1 < firstUpstreamFibers.size(); ++iF1){
    
    size_t firstFiber = firstUpstreamFibers[iF1];

    for(size_t iF2 = 0; iF2 < secondUpstreamFibers.size(); ++iF2){
      size_t secondFiber = secondUpstreamFibers[iF2];

      std::cout << "Paired: " << firstFiber << " " << secondFiber << std::endl; 
      upstreamPairedFibers.push_back(std::make_pair(firstFiber, secondFiber));

      if (iF2 < secondUpstreamFibers.size() - 1){
        if (secondUpstreamFibers[iF2] == (secondUpstreamFibers[iF2 + 1] - 1)) ++iF2;
      }
    }

    if (iF1 < firstUpstreamFibers.size() - 1){
      if (firstUpstreamFibers[iF1] == (firstUpstreamFibers[iF1 + 1] - 1)) ++iF1;
    }
  }
 
  for(size_t iF = 0; iF < upstreamPairedFibers.size(); ++iF){
    
    std::pair<size_t,size_t> thePair = upstreamPairedFibers.at(iF);
 
    if(firstUpstreamType == "horiz" && secondUpstreamType == "vert"){
      double xPos = GetPosition(firstUpstreamName, thePair.first);
      double yPos = GetPosition(secondUpstreamName, thePair.second);
      
      std::cout << "normal " << xPos << " " << yPos <<  std::endl;
      TVector3 posInDet = ConvertProfCoordinates(xPos,yPos,0.,fFirstTrackingProfZ);
      std::cout << posInDet.X() << " " << posInDet.Y() << " " << posInDet.Z() << std::endl;
      upstreamPositions.push_back( posInDet );
    }
    else if(firstUpstreamType == "vert" && secondUpstreamType == "horiz"){
      double yPos = GetPosition(firstUpstreamName, thePair.first);
      double xPos = GetPosition(secondUpstreamName, thePair.second);
      std::cout << "normal " << xPos << " " << yPos <<  std::endl;

      TVector3 posInDet = ConvertProfCoordinates(xPos,yPos,0.,fFirstTrackingProfZ);
      std::cout << posInDet.X() << " " << posInDet.Y() << " " << posInDet.Z() << std::endl;
      upstreamPositions.push_back( posInDet );
    }

  }
  
  std::cout << "Downstream" << std::endl;
  //Pair the downstream fibers together
  for(size_t iF1 = 0; iF1 < firstDownstreamFibers.size(); ++iF1){
    
    size_t firstFiber = firstDownstreamFibers[iF1];

    for(size_t iF2 = 0; iF2 < secondDownstreamFibers.size(); ++iF2){
      size_t secondFiber = secondDownstreamFibers[iF2];

      std::cout << "Paired: " << firstFiber << " " << secondFiber << std::endl; 
      downstreamPairedFibers.push_back(std::make_pair(firstFiber, secondFiber));

      if (iF2 < secondDownstreamFibers.size() - 1){
        if (secondDownstreamFibers[iF2] == (secondDownstreamFibers[iF2 + 1] - 1)) ++iF2;
      }
    }

    if (iF1 < firstDownstreamFibers.size() - 1){
      if (firstDownstreamFibers[iF1] == (firstDownstreamFibers[iF1 + 1] - 1)) ++iF1;
    }
  }

  for(size_t iF = 0; iF < downstreamPairedFibers.size(); ++iF){
    
    std::pair<size_t,size_t> thePair = downstreamPairedFibers.at(iF);
 
    if(firstDownstreamType == "horiz" && secondDownstreamType == "vert"){
      double xPos = GetPosition(firstDownstreamName, thePair.first);
      double yPos = GetPosition(secondDownstreamName, thePair.second);

      std::cout << "normal " << xPos << " " << yPos <<  std::endl;
      TVector3 posInDet = ConvertProfCoordinates(xPos,yPos,0.,fSecondTrackingProfZ);
      std::cout << posInDet.X() << " " << posInDet.Y() << " " << posInDet.Z() << std::endl;
      downstreamPositions.push_back( posInDet );
    }
    else if(firstDownstreamType == "vert" && secondDownstreamType == "horiz"){
      double yPos = GetPosition(firstDownstreamName, thePair.first);
      double xPos = GetPosition(secondDownstreamName, thePair.second);

      std::cout << "normal " << xPos << " " << yPos <<  std::endl;
      TVector3 posInDet = ConvertProfCoordinates(xPos,yPos,0.,fSecondTrackingProfZ);
      std::cout << posInDet.X() << " " << posInDet.Y() << " " << posInDet.Z() << std::endl;
      downstreamPositions.push_back( posInDet );
    }

  }
 
  //Just for creating tracks
  std::vector< std::vector<double> > dummy;
  std::vector<double> mom(3, util::kBogusD);
  /// 

  for(size_t iU = 0; iU < upstreamPositions.size(); ++iU){
    for(size_t iD = 0; iD < downstreamPositions.size(); ++iD){
      std::vector<TVector3> thePoints;
      thePoints.push_back(upstreamPositions.at(iU));
      thePoints.push_back(downstreamPositions.at(iD));

      //Now project the last point to the TPC face
      thePoints.push_back( ProjectToTPC(thePoints[0],thePoints[1]) );    

      
      std::vector<TVector3> theMomenta;
      //Just push back the unit vector for each point 
      //For now.
      //Eventually, use momentum from curvature?
      theMomenta.push_back( ( downstreamPositions.at(iD) - upstreamPositions.at(iU) ).Unit() );
      theMomenta.push_back( ( downstreamPositions.at(iD) - upstreamPositions.at(iU) ).Unit() );
      theMomenta.push_back( ( downstreamPositions.at(iD) - upstreamPositions.at(iU) ).Unit() );

      recob::Track * tempTrack = new recob::Track(thePoints, theMomenta, dummy, mom, 1);      
      theTracks.push_back(tempTrack);
    }
  }

}

void proto::BeamEvent::MomentumSpec(size_t theTrigger){
  
  std::cout << "Doing momentum spectrometry for trigger " << beamevt->GetFullT0(theTrigger) << std::endl;

  //Get the active fibers from the upstream tracking XBPF
  std::string firstBPROF1Type    = fDeviceTypes[firstBPROF1]; 
  std::string secondBPROF1Type   = fDeviceTypes[secondBPROF1]; 
  std::vector<short> BPROF1Fibers;
  std::string BPROF1Name;

  if (firstBPROF1Type == "horiz" && secondBPROF1Type == "vert"){
    BPROF1Fibers = beamevt->GetActiveFibers(firstBPROF1, theTrigger);
    BPROF1Name = firstBPROF1;

    std::cout << firstBPROF1 << " has " << BPROF1Fibers.size() << " active fibers at time " << beamevt->GetFiberTime(firstBPROF1,theTrigger) << std::endl;
    for(size_t i = 0; i < BPROF1Fibers.size(); ++i){
      std::cout << BPROF1Fibers[i] << " ";
    }
    std::cout << std::endl;

  }
  else if(secondBPROF1Type == "horiz" && firstBPROF1Type == "vert"){
    BPROF1Fibers = beamevt->GetActiveFibers(secondBPROF1, theTrigger);
    BPROF1Name = secondBPROF1;

    std::cout << secondBPROF1 << " has " << BPROF1Fibers.size() << " active fibers at time " << beamevt->GetFiberTime(secondBPROF1,theTrigger) << std::endl;
    for(size_t i = 0; i < BPROF1Fibers.size(); ++i){
      std::cout << BPROF1Fibers[i] << " ";
    }
    std::cout << std::endl;
  }
  else{
    std::cout << "Error: type is not correct" << std::endl;
    return;
  }


  //////////////////////////////////////////////

  if( (BPROF1Fibers.size() < 1) ){
    std::cout << "Warning, at least one empty Beam Profiler. Not checking momentum" << std::endl;
    return;
  }
 /* else if( (BPROF1Fibers.size() > 5) ){
    std::cout << "Warning, too many (>5) active fibers in at least one Beam Profiler. Not checking momentum" << std::endl;
    return;
  }
*/
  //We have the active Fibers, now go through them.
  //Skip the second of any adjacents 
  std::vector< short > strippedFibers; 
  std::cout << "BPROF1" << std::endl;
  for(size_t iF1 = 0; iF1 < BPROF1Fibers.size(); ++iF1){
    
    size_t Fiber = BPROF1Fibers[iF1];
    strippedFibers.push_back(Fiber);

    if (iF1 < BPROF1Fibers.size() - 1){
      if (BPROF1Fibers[iF1] == (BPROF1Fibers[iF1 + 1] - 1)) ++iF1;
    }
  }
  
  double X1, X2, X3;

  //X-direction convention is reversed for the spectrometer   
  X1 = -1.*GetPosition(BPROF1Name, strippedFibers[0])/1.E3;

  //BPROF2////
  //
  std::vector<short>  BPROF2Fibers = beamevt->GetActiveFibers(BPROF2, theTrigger);
  if( (BPROF2Fibers.size() < 1) ){
    std::cout << "Warning, at least one empty Beam Profiler. Not checking momentum" << std::endl;
    return;
  }
  std::cout << BPROF2 << " has " << BPROF2Fibers.size() << " active fibers at time " << beamevt->GetFiberTime(BPROF2,theTrigger) << std::endl;
  for(size_t i = 0; i < BPROF2Fibers.size(); ++i){
    std::cout << BPROF2Fibers[i] << " ";
  }
  std::cout << std::endl;

  strippedFibers.clear(); 
  std::cout << "BPROF2" << std::endl;
  for(size_t iF1 = 0; iF1 < BPROF2Fibers.size(); ++iF1){
    
    size_t Fiber = BPROF2Fibers[iF1];
    strippedFibers.push_back(Fiber);

    if (iF1 < BPROF2Fibers.size() - 1){
      if (BPROF2Fibers[iF1] == (BPROF2Fibers[iF1 + 1] - 1)) ++iF1;
    }
  }
 
  X2 = -1.*GetPosition(BPROF2, strippedFibers[0])/1.E3; 
  //X2 = X2*cos(TMath::Pi() - fBeamBend/1000.);
  ////////////

  //BPROF3////
  //
  std::vector<short>  BPROF3Fibers = beamevt->GetActiveFibers(BPROF3, theTrigger);
  if( (BPROF3Fibers.size() < 1) ){
    std::cout << "Warning, at least one empty Beam Profiler. Not checking momentum" << std::endl;
    return;
  }
  std::cout << BPROF3 << " has " << BPROF3Fibers.size() << " active fibers at time " << beamevt->GetFiberTime(BPROF3,theTrigger) << std::endl;
  for(size_t i = 0; i < BPROF3Fibers.size(); ++i){
    std::cout << BPROF3Fibers[i] << " ";
  }
  std::cout << std::endl;

  strippedFibers.clear(); 
  std::cout << "BPROF3" << std::endl;
  for(size_t iF1 = 0; iF1 < BPROF3Fibers.size(); ++iF1){
    
    size_t Fiber = BPROF3Fibers[iF1];
    strippedFibers.push_back(Fiber);

    if (iF1 < BPROF3Fibers.size() - 1){
      if (BPROF3Fibers[iF1] == (BPROF3Fibers[iF1 + 1] - 1)) ++iF1;
    }
  }
 
  X3 = -1.*GetPosition(BPROF3, strippedFibers[0])/1.E3; 
  //X3 = X3*cos(TMath::Pi() - fBeamBend/1000.);
  ////////////
  
  std::cout << "fBeamBend " << fBeamBend  << std::endl;
  double cosTheta = MomentumCosTheta(X1,X2,X3);
  std::cout <<"CosTheta " << cosTheta << std::endl;
  std::cout <<"Theta " << acos(cosTheta) << std::endl;
  double LB = mag_P1*fabs(current[0]);
  double deltaI = fabs(current[0]) - mag_P4;
  if(deltaI>0) LB+= mag_P3*deltaI*deltaI;

  double momentum = 299792458*LB/(1.E9 * acos(cosTheta));
  fMomentumHist->Fill(momentum);
  if( (BPROF1Fibers.size() == 1) && (BPROF2Fibers.size() == 1) && (BPROF3Fibers.size() == 1) ){
    std::cout << "Filling Cut Momentum Spectrum" << std::endl;
    fCutMomentum->Fill(momentum);

    beamevt->AddRecoBeamMomentum(momentum);
  }

  std::cout << "Momentum: " << 299792458*LB/(1.E9 * acos(cosTheta)) << std::endl; 


  std::cout << "Getting all trio-wise hits" << std::endl;
  std::cout << "N1,N2,N3 " << BPROF1Fibers.size()
            << " "         << BPROF2Fibers.size() 
            << " "         << BPROF3Fibers.size() << std::endl;
  for(size_t i1 = 0; i1 < BPROF1Fibers.size(); ++i1){
    double x1,x2,x3;

    x1 = -1.*GetPosition(BPROF1Name, BPROF1Fibers[i1])/1.E3;
    if (i1 < BPROF1Fibers.size() - 1){
      if (BPROF1Fibers[i1] == (BPROF1Fibers[i1 + 1] - 1)){
        //Add .5 mm
        x1 += .0005;
      }
    }

    for(size_t i2 = 0; i2 < BPROF2Fibers.size(); ++i2){
      x2 = -1.*GetPosition(BPROF2, BPROF2Fibers[i2])/1.E3;
      if (i2 < BPROF2Fibers.size() - 1){
        if (BPROF2Fibers[i2] == (BPROF2Fibers[i2 + 1] - 1)){
          //Add .5 mm
          x2 += .0005;
        }
      }

      for(size_t i3 = 0; i3 < BPROF3Fibers.size(); ++i3){
        std::cout << "\t" << i1 << " " << i2 << " " << i3 << std::endl;
        x3 = -1.*GetPosition(BPROF3, BPROF3Fibers[i3])/1.E3;
        if (i3 < BPROF3Fibers.size() - 1){
          if (BPROF3Fibers[i3] == (BPROF3Fibers[i3 + 1] - 1)){
            //Add .5 mm
            x3 += .0005;
          }
        }

        double cosTheta_full = MomentumCosTheta(x1,x2,x3);        
        double momentum_full = 299792458*LB/(1.E9 * acos(cosTheta_full));
        fFullMomentum->Fill(momentum_full);

        if (i3 < BPROF3Fibers.size() - 1){
          if (BPROF3Fibers[i3] == (BPROF3Fibers[i3 + 1] - 1)){
            //Skip the next
            ++i3;
          }
        }
        
      }

      if (i2 < BPROF2Fibers.size() - 1){
        if (BPROF2Fibers[i2] == (BPROF2Fibers[i2 + 1] - 1)){
          //Skip the next
          ++i2;
        }
      }
    }

    if (i1 < BPROF1Fibers.size() - 1){
      if (BPROF1Fibers[i1] == (BPROF1Fibers[i1 + 1] - 1)){
        //Skip the next
        ++i1;
      }
    }
  }
}

double proto::BeamEvent::MomentumCosTheta(double X1, double X2, double X3){
/*
  double a = cos(fBeamBend)*(X3*L2 - X2*L3) / (L3 - L2);

  double cosTheta = (a + X1)*( (L3 - L2)*tan(fBeamBend) + (X2 - X3)*cos(fBeamBend) )+ (L3 - L2)*L1 ;

  double denomTerm1, denomTerm2, denom; 
  denomTerm1 = sqrt( L1*L1 + (a + X1)*(a + X1) );
  denomTerm2 = sqrt( 
               (L3 - L2)*(L3 - L2) 
             + TMath::Power( ( (L3 - L2)*tan(fBeamBend) 
                             + (X2 - X3)*cos(fBeamBend) ), 2 ) );
                            
  denom = denomTerm1*denomTerm2;
  cosTheta = cosTheta / denom;

  //
 
  double a = (X2*L3 - X3*L2) / (L3 - L2);
   
  double cosTheta = (a - X1)*(X3 - X2)*cos(fBeamBend) + (L3 - L2 )*( L1 + (a - X1)*tan(fBeamBend) );
    
  double denomTerm1, denomTerm2, denom; 
  denomTerm1 = sqrt( L1*L1 + (a - X1)*(a - X1) );
  denomTerm2 = sqrt( 
               (L3 - L2)*(L3 - L2) 
             + TMath::Power( ( (L3 - L2)*tan(fBeamBend) 
                             + (X3 - X2)*cos(fBeamBend) ), 2 ) );
  
  denom = denomTerm1*denomTerm2;
  cosTheta = cosTheta / denom;
 */

    //double a = L2*tan(fBeamBend) + X2*cos(fBeamBend) - ( (L3 - L2)*tan(fBeamBend) + (X3 - X2)*cos(fBeamBend) )*( L2 - X2*sin(fBeamBend) );
    double a =  ( (L3 - L2)*tan(fBeamBend) + (X3 - X2)*cos(fBeamBend) )*( L2 - X2*sin(fBeamBend) );
    a = a / (L3 - X3*sin(fBeamBend) - L2 + X2*sin(fBeamBend) );
    a = L2*tan(fBeamBend) + X2*cos(fBeamBend) - a;
  
    double numTerm = (a - X1)*(L3 - L2)*tan(fBeamBend) + (a - X1)*(X3 - X2)*cos(fBeamBend) + L1*( (L3 - L2) - (X3 - X2)*sin(fBeamBend) );
  
    double denomTerm1, denomTerm2, denom;
    denomTerm1 = sqrt( L1*L1 + (a - X1)*(a - X1) );
    denomTerm2 = sqrt( TMath::Power( ( (L3 - L2)*tan(fBeamBend) + (X3 - X2)*cos(fBeamBend) ),2)
                     + TMath::Power( ( (L3 - L2)                - (X3 - X2)*sin(fBeamBend) ),2) );
    denom = denomTerm1 * denomTerm2; 

    double cosTheta = numTerm/denom;
  
  return cosTheta;
}

//Gets info from FBMs matching in time
void proto::BeamEvent::GetPairedFBMInfo(beam::ProtoDUNEBeamEvent beamevt, double Time){

  //This method goes through the paired FBMs and gets 2D hits in their profile that match in time
  //to the input time and fills 2D histograms. 
  //If there are any missing triggers in a pair, then it skips filling the respective hist
  //
  //Need to figure out what to do with multiple active fibers

  std::map<std::string, std::vector<size_t> > goodTriggers;
  for(size_t ip = 0; ip < fPairedDevices.size(); ++ip){
    std::string name = fPairedDevices[ip].first;
    std::vector<size_t> triggers;
//    std::cout << ip << " " << name << " " << fPairedDevices[ip].second << std::endl;
    for(size_t itN = 0; itN < beamevt.GetNFBMTriggers(name); ++itN){
      if ( ( (beamevt.DecodeFiberTime(name, itN, fOffsetTAI ) - Time) < fTolerance ) && ( (beamevt.DecodeFiberTime(name, itN, fOffsetTAI ) - Time)     >= 0 ) ){
//        std::cout << "Found Good Time " << name << " " << beamevt.DecodeFiberTime(name, itN ) << std::endl;
        triggers.push_back(itN);
      }
    }
    goodTriggers[name] = triggers;

    name = fPairedDevices[ip].second;
    triggers.clear();
//    std::cout << ip << " " << name << " " << fPairedDevices[ip].second << std::endl;
    for(size_t itN = 0; itN < beamevt.GetNFBMTriggers(name); ++itN){
      if ( ( (beamevt.DecodeFiberTime(name, itN, fOffsetTAI ) - Time) < fTolerance ) && ( (beamevt.DecodeFiberTime(name, itN, fOffsetTAI ) - Time)     >= 0 ) ){
//        std::cout << "Found Good Time " << name << " " << beamevt.DecodeFiberTime(name, itN ) << std::endl;
        triggers.push_back(itN);
      }
    }
    goodTriggers[name] = triggers;
  }
   
  for(size_t ip = 0; ip < fPairedDevices.size(); ++ip){
    std::string nameOne = fPairedDevices[ip].first;
    std::string nameTwo = fPairedDevices[ip].second;

    size_t triggerSizeOne = goodTriggers[nameOne].size();
    size_t triggerSizeTwo = goodTriggers[nameTwo].size();

    if(triggerSizeOne < 1){
//      std::cout << "Missing trigger for " << nameOne << std::endl;
    }
    if(triggerSizeTwo < 1){
//      std::cout << "Missing trigger for " << nameTwo << std::endl;
    }

    if(triggerSizeOne > 0 && triggerSizeTwo > 0){
      //std::cout << "Found good triggers for " << nameOne << " " << nameTwo << std::endl;
      size_t triggerOne = goodTriggers[nameOne][0];
      size_t triggerTwo = goodTriggers[nameTwo][0];
      fBeamProf2D[ip]->Fill(beamevt.GetActiveFibers(nameOne,triggerOne)[0], beamevt.GetActiveFibers(nameTwo,triggerTwo)[0]);
    }
  }

}

void proto::BeamEvent::GetPairedStraightFBMInfo(beam::ProtoDUNEBeamEvent beamevt, double Time){

  //This method goes through the paired FBMs and gets 2D hits in their profile that match in time
  //to the input time and fills 2D histograms. 
  //If there are any missing triggers in a pair, then it skips filling the respective hist
  //
  //Need to figure out what to do with multiple active fibers

  std::map<std::string, std::vector<size_t> > goodTriggers;
  for(size_t ip = 0; ip < fPairedStraightDevices.size(); ++ip){
    std::string name = fPairedStraightDevices[ip].first;
    std::vector<size_t> triggers;
    //std::cout << ip << " " << name << " " << fPairedStraightDevices[ip].second << std::endl;
    for(size_t itN = 0; itN < beamevt.GetNFBMTriggers(name); ++itN){
      if ( ( (beamevt.DecodeFiberTime(name, itN, fOffsetTAI ) - Time) < fTolerance ) && ( (beamevt.DecodeFiberTime(name, itN, fOffsetTAI ) - Time)     >= 0 ) ){
        //std::cout << "Found Good Time " << name << " " << beamevt.DecodeFiberTime(name, itN ) << std::endl;
        triggers.push_back(itN);
      }
    }
    goodTriggers[name] = triggers;

    name = fPairedStraightDevices[ip].second;
    triggers.clear();
    //std::cout << ip << " " << name << " " << fPairedStraightDevices[ip].second << std::endl;
    for(size_t itN = 0; itN < beamevt.GetNFBMTriggers(name); ++itN){
      if ( ( (beamevt.DecodeFiberTime(name, itN, fOffsetTAI ) - Time) < fTolerance ) && ( (beamevt.DecodeFiberTime(name, itN, fOffsetTAI ) - Time)     >= 0 ) ){
        //std::cout << "Found Good Time " << name << " " << beamevt.DecodeFiberTime(name, itN ) << std::endl;
        triggers.push_back(itN);
      }
    }
    goodTriggers[name] = triggers;
  }
   
  for(size_t ip = 0; ip < fPairedStraightDevices.size(); ++ip){
    std::string nameOne = fPairedStraightDevices[ip].first;
    std::string nameTwo = fPairedStraightDevices[ip].second;

    size_t triggerSizeOne = goodTriggers[nameOne].size();
    size_t triggerSizeTwo = goodTriggers[nameTwo].size();

    if(triggerSizeOne < 1){
      //std::cout << "Missing trigger for " << nameOne << std::endl;
    }
    if(triggerSizeTwo < 1){
      //std::cout << "Missing trigger for " << nameTwo << std::endl;
    }

    if(triggerSizeOne > 0 && triggerSizeTwo > 0){
      //std::cout << "Found good triggers for " << nameOne << " " << nameTwo << std::endl;
      size_t triggerOne = goodTriggers[nameOne][0];
      size_t triggerTwo = goodTriggers[nameTwo][0];
      fBeamProf2D[ip]->Fill(beamevt.GetActiveFibers(nameOne,triggerOne)[0], beamevt.GetActiveFibers(nameTwo,triggerTwo)[0]);
    }
  }

}

void proto::BeamEvent::GetUnpairedFBMInfo(beam::ProtoDUNEBeamEvent beamevt, double Time){
  //This method goes through the unpaired devices, as well as the paired devices individually
  //and gets the info that matches to the inputted time
  //
  //Currently allows for multiple fibers in an event


  //Going through paired devices individually
  for(size_t ip = 0; ip < fPairedDevices.size(); ++ip){
    std::string name = fPairedDevices[ip].first;
    //std::cout << ip << " " << name << " " << fPairedDevices[ip].second << std::endl;
    for(size_t itN = 0; itN < beamevt.GetNFBMTriggers(name); ++itN){
      if ( ( (beamevt.DecodeFiberTime(name, itN, fOffsetTAI ) - Time) < fTolerance ) && ( (beamevt.DecodeFiberTime(name, itN, fOffsetTAI ) - Time)     >= 0 ) ){
        //std::cout << "Found Good Time " << name << " " << beamevt.DecodeFiberTime(name, itN ) << std::endl;
        for(size_t iFiber = 0; iFiber < beamevt.GetActiveFibers(name,itN).size(); ++iFiber){
	  fBeamProf1D[name]->Fill(beamevt.GetActiveFibers(name, itN)[iFiber]);
	} 
      }
    }

    name = fPairedDevices[ip].second;
    //std::cout << ip << " " << name << " " << fPairedDevices[ip].second << std::endl;
    for(size_t itN = 0; itN < beamevt.GetNFBMTriggers(name); ++itN){
      if ( ( (beamevt.DecodeFiberTime(name, itN, fOffsetTAI ) - Time) < fTolerance ) && ( (beamevt.DecodeFiberTime(name, itN, fOffsetTAI ) - Time)     >= 0 ) ){
        //std::cout << "Found Good Time " << name << " " << beamevt.DecodeFiberTime(name, itN ) << std::endl;
        for(size_t iFiber = 0; iFiber < beamevt.GetActiveFibers(name,itN).size(); ++iFiber){
	  fBeamProf1D[name]->Fill(beamevt.GetActiveFibers(name, itN)[iFiber]);
	} 
      }
    }
  }

  //Now going through unpaired
  for(size_t id = 0; id < fDevices.size(); ++id){
    std::string name = fDevices[id];
    //std::cout << id << " " << name << " " << fDevices[id] << std::endl;
    for(size_t itN = 0; itN < beamevt.GetNFBMTriggers(name); ++itN){
      if ( ( (beamevt.DecodeFiberTime(name, itN, fOffsetTAI ) - Time) < fTolerance ) && ( (beamevt.DecodeFiberTime(name, itN, fOffsetTAI ) - Time)     >= 0 ) ){
        //std::cout << "Found Good Time " << name << " " << beamevt.DecodeFiberTime(name, itN ) << std::endl;
        for(size_t iFiber = 0; iFiber < beamevt.GetActiveFibers(name,itN).size(); ++iFiber){
	  fBeamProf1D[name]->Fill(beamevt.GetActiveFibers(name, itN)[iFiber]);
	} 
      }
    }
  }

}

TVector3 proto::BeamEvent::ProjectToTPC(TVector3 firstPoint, TVector3 secondPoint){
  TVector3 dR = (secondPoint - firstPoint);
  
  double deltaZ = -1.*secondPoint.Z();
  double deltaX = deltaZ * (dR.X() / dR.Z());
  double deltaY = deltaZ * (dR.Y() / dR.Z());

  TVector3 lastPoint = secondPoint + TVector3(deltaX, deltaY, deltaZ);
  return lastPoint;
}

double proto::BeamEvent::GetPosition(std::string deviceName, int fiberIdx){
  //NEEDS WORK
  if(fiberIdx > 192){ std::cout << "Please select fiber in range [0,191]" << std::endl; return -1.;}
  double size = fFiberDimension[deviceName];
  //double size = 1.;
  
  //Define 0th fiber as farthest positive. Last fiber is farthest negative. Center is between 96 and 97 
  double pos = size*(96 - fiberIdx) - size/2.;
  return pos;
}

TVector3 proto::BeamEvent::TranslateDeviceToDetector(TVector3 globalDeviceCoords){
  //fGlobalDetCoords given by fcl
  //Translate position of device w.r.t. Detector
  TVector3 inDetCoords = globalDeviceCoords - fGlobalDetCoords;
   
  //Rotate into detector coordinates
  inDetCoords.RotateX(fDetRotation[0]);
  inDetCoords.RotateY(fDetRotation[1]);
  inDetCoords.RotateZ(fDetRotation[2]);
  return inDetCoords;
}


DEFINE_ART_MODULE(proto::BeamEvent)
