#include "dune/Protodune/Analysis/ProtoDUNEDataUtils.h"

#include "lardataobj/RawData/RDTimeStamp.h"
#include "lardataobj/RawData/RawDigit.h"

#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"

#include "dune-raw-data/Services/ChannelMap/PdspChannelMapService.h"
#include "dune/DuneObj/ProtoDUNEBeamEvent.h"

protoana::ProtoDUNEDataUtils::ProtoDUNEDataUtils(fhicl::ParameterSet const& p){
  this->reconfigure(p);
}

protoana::ProtoDUNEDataUtils::~ProtoDUNEDataUtils(){

}

void protoana::ProtoDUNEDataUtils::reconfigure(fhicl::ParameterSet const& p){
  fTimingTag = p.get<art::InputTag>("TimingTag");
  fBeamEventTag = p.get<art::InputTag>("BeamEventTag");

  fStrictTOF = p.get<bool>("StrictTOF");
  fStrictCherenkov = p.get<bool>("StrictCherenkov");

  fTOFDistance = p.get<float>("TOFDistance"); // meters
  fMomentumScaleFactor = p.get<float>("MomentumScaleFactor"); 
  fMomentumOffset = p.get<float>("MomentumOffset"); // GeV/c
  fTOFScaleFactor = p.get<float>("TOFScaleFactor");
  fTOFOffset = p.get<float>("TOFOffset"); // ns

  fTOFElectronCuts = p.get<std::vector<float> >("TOFElectronCuts"); // ns
  fTOFMuonCuts = p.get<std::vector<float> >("TOFMuonCuts"); // ns
  fTOFPionCuts = p.get<std::vector<float> >("TOFPionCuts"); // ns
  fTOFKaonCuts = p.get<std::vector<float> >("TOFKaonCuts"); // ns
  fTOFProtonCuts = p.get<std::vector<float> >("TOFProtonCuts"); // ns
  if(fTOFElectronCuts.size() != 2)
  {
    throw cet::exception("ProtoDUNEDataUtils")<< "TOFElectronCuts parameter needs to be a vector of size 2 not "<<fTOFElectronCuts.size();
  }
  if(fTOFMuonCuts.size() != 2)
  {
    throw cet::exception("ProtoDUNEDataUtils")<< "TOFMuonCuts parameter needs to be a vector of size 2 not "<<fTOFMuonCuts.size();
  }
  if(fTOFPionCuts.size() != 2)
  {
    throw cet::exception("ProtoDUNEDataUtils")<< "TOFPionCuts parameter needs to be a vector of size 2 not "<<fTOFPionCuts.size();
  }
  if(fTOFKaonCuts.size() != 2)
  {
    throw cet::exception("ProtoDUNEDataUtils")<< "TOFKaonCuts parameter needs to be a vector of size 2 not "<<fTOFKaonCuts.size();
  }
  if(fTOFProtonCuts.size() != 2)
  {
    throw cet::exception("ProtoDUNEDataUtils")<< "TOFProtonCuts parameter needs to be a vector of size 2 not "<<fTOFProtonCuts.size();
  }
}

// Access the trigger information to see if this is a beam trigger
bool protoana::ProtoDUNEDataUtils::IsBeamTrigger(art::Event const & evt) const{

  bool isBeam = false;

  // Accessing the trigger information as done in DataPrepModule
  // The information is stored in the time stamps
  art::Handle<std::vector<raw::RDTimeStamp>> timeStamps;
  evt.getByLabel(fTimingTag,timeStamps);

  // Return false if we have no time stamps
  if(!timeStamps.isValid()) return isBeam;
  // We should only have one RDTimeStamp
  if(timeStamps->size() > 1) return isBeam;
  
  // Access the trigger information. Beam trigger flag = 0xc
  const raw::RDTimeStamp& timeStamp = timeStamps->at(0);
  isBeam = (timeStamp.GetFlags() == 0xc);
  
  return isBeam;
}

// ----------------------------------------------------------------------------
int protoana::ProtoDUNEDataUtils::GetNActiveFembsForAPA(art::Event const & evt, int apa) const {


// Get pd channel map
  art::ServiceHandle<dune::PdspChannelMapService> channelMap;

  // set only saves unique elements
  std::set<int> apaset;

  // Get raw digits time stamps
  art::Handle< std::vector<raw::RDTimeStamp> > RawdigitTSListHandle;
  std::vector<art::Ptr<raw::RDTimeStamp> > digitTSlist;
  

   // Get raw digits
  art::Handle< std::vector<raw::RawDigit> > RawdigitListHandle;
  std::vector<art::Ptr<raw::RawDigit> > digitlist;

  if (evt.getByLabel("tpcrawdecoder", "daq", RawdigitListHandle)){
  	
    art::fill_ptr_vector(digitlist, RawdigitListHandle);  

    for(auto const & dptr : digitlist) {
    const raw::RawDigit& digit = *dptr;
    // Get the channel number for this digit
    uint32_t chan = digit.Channel();
    
    int iapa = channelMap->APAFromOfflineChannel(chan);
    if(iapa != apa) continue;
    // Get the channel FEMB and WIB
    int WIB = channelMap->WIBFromOfflineChannel(chan); // 0-4
    int FEMB = channelMap->FEMBFromOfflineChannel(chan); // 1-4
    //int FEMBchan = channelMap->FEMBChannelFromOfflineChannel(chan);
    int iFEMB = ((WIB*4)+(FEMB-1)); //index of the FEMB 0-19

    apaset.insert(iFEMB);
  }
}

else{ // if raw digits have been dropped use RDTimeStamps instead
	evt.getByLabel("tpcrawdecoder", "daq", RawdigitTSListHandle);
	
	art::fill_ptr_vector(digitTSlist, RawdigitTSListHandle);  
	
  for(auto const & dptr : digitTSlist) {
  	
    const raw::RDTimeStamp & digit = *dptr;
    
    // Get the channel number for this digit
    uint16_t chan = digit.GetFlags();
    
    int iapa = channelMap->APAFromOfflineChannel(chan);
    if(iapa != apa) continue;
    // Get the channel FEMB and WIB
    int WIB = channelMap->WIBFromOfflineChannel(chan); // 0-4
    int FEMB = channelMap->FEMBFromOfflineChannel(chan); // 1-4
    //int FEMBchan = channelMap->FEMBChannelFromOfflineChannel(chan);
    int iFEMB = ((WIB*4)+(FEMB-1)); //index of the FEMB 0-19

    apaset.insert(iFEMB);
  }
}

  return (apaset.size());

}
// ----------------------------------------------------------------------------
bool protoana::ProtoDUNEDataUtils::IsGoodBeamlineTrigger(art::Event const & evt) const{
  return true;
}

std::vector<double> protoana::ProtoDUNEDataUtils::GetBeamlineMass(art::Event const & evt) const{
  std::vector<double> tofs;
  std::vector<double> momenta;
  std::vector<double> result;

  std::vector<art::Ptr<beam::ProtoDUNEBeamEvent>> beamVec;
  auto beamHand = evt.getValidHandle<std::vector<beam::ProtoDUNEBeamEvent>>(fBeamEventTag);
  if(beamHand.isValid())
  {
    art::fill_ptr_vector(beamVec, beamHand);
  }

  for(size_t iBeamEvent=0; iBeamEvent < beamVec.size(); iBeamEvent++)
  {
    const beam::ProtoDUNEBeamEvent& beamEvent = *(beamVec.at(iBeamEvent));
    tofs.push_back(beamEvent.GetTOF());

    const std::vector<double> & beamMomenta = beamEvent.GetRecoBeamMomenta();
    for(size_t iMom=0; iMom < beamMomenta.size(); iMom++)
    {
      momenta.push_back(beamEvent.GetRecoBeamMomentum(iMom));
    }
  }
  for(const auto & tof : tofs)
  {
    for(const auto & momentum : momenta)
    {
        double massSquared = std::pow(momentum*fMomentumScaleFactor-fMomentumOffset,2) 
                    * ((tof*fTOFScaleFactor-fTOFOffset)/(fTOFDistance/2.99e8*1e9) - 1.);
        double mass = std::sqrt(massSquared);
        result.push_back(mass);
    }
  }
  return result;
}

protoana::PossibleParticleCands protoana::ProtoDUNEDataUtils::GetCherenkovParticleID(art::Event const & evt, const float beamEnergyGeV) const{
  std::vector<art::Ptr<beam::ProtoDUNEBeamEvent>> beamVec;
  auto beamHand = evt.getValidHandle<std::vector<beam::ProtoDUNEBeamEvent>>(fBeamEventTag);
  if(beamHand.isValid())
  {
    art::fill_ptr_vector(beamVec, beamHand);
  }
  int CKov0Status = -5;
  int CKov1Status = -5;
  for(size_t iBeamEvent=0; iBeamEvent < beamVec.size(); iBeamEvent++)
  {
    const beam::ProtoDUNEBeamEvent& beamEvent = *(beamVec.at(iBeamEvent));
    CKov0Status = beamEvent.GetCKov0Status();
    CKov1Status = beamEvent.GetCKov1Status();
  }
  if (CKov0Status < 0 || CKov1Status < 0)
  {
    if (fStrictCherenkov)
    {
      protoana::PossibleParticleCands result = {false,false,false,false,false};
      return result;
    }
    else
    {
      protoana::PossibleParticleCands result = {true,true,true,true,true};
      return result;
    }
  }
  if (beamEnergyGeV < 2.5)
  {
    if(CKov1Status == 1)
    {
      protoana::PossibleParticleCands result = {true,false,false,false,false};
      return result;
    }
    else
    {
      protoana::PossibleParticleCands result = {false,true,true,true,true};
      return result;
    }
  }
  else if (beamEnergyGeV < 3.5)
  {
    if(CKov0Status == 1 && CKov1Status == 1) // electron
    {
      protoana::PossibleParticleCands result = {true,false,false,false,false};
      return result;
    }
    else if(CKov1Status == 1) // pi/mu
    {
      protoana::PossibleParticleCands result = {false,true,true,false,false};
      return result;
    }
    else // kaon/proton
    {
      protoana::PossibleParticleCands result = {false,false,false,true,true};
      return result;
    }
  }
  else // 6 and 7 GeV
  {
    if(CKov0Status == 1 && CKov1Status == 1) // electron/muon/pion
    {
      protoana::PossibleParticleCands result = {true,true,true,false,false};
      return result;
    }
    else if(CKov1Status == 1) // kaon
    {
      protoana::PossibleParticleCands result = {false,false,false,true,false};
      return result;
    }
    else // proton
    {
      protoana::PossibleParticleCands result = {false,false,false,false,true};
      return result;
    }
  }
}

protoana::PossibleParticleCands protoana::ProtoDUNEDataUtils::GetTOFParticleID(art::Event const & evt, const float beamEnergyGeV) const{
  protoana::PossibleParticleCands result = {false,false,false,false};
  std::vector<art::Ptr<beam::ProtoDUNEBeamEvent>> beamVec;
  auto beamHand = evt.getValidHandle<std::vector<beam::ProtoDUNEBeamEvent>>(fBeamEventTag);
  if(beamHand.isValid())
  {
    art::fill_ptr_vector(beamVec, beamHand);
  }
  float TOF = -1e6;
  for(size_t iBeamEvent=0; iBeamEvent < beamVec.size(); iBeamEvent++)
  {
    const beam::ProtoDUNEBeamEvent& beamEvent = *(beamVec.at(iBeamEvent));
    if (!beamEvent.CheckIsMatched())
    {
      if(fStrictTOF)
      {
        result = {false,false,false,false,false};
        return result;
      }
      else
      {
        result = {true,true,true,true,true};
        return result;
      }
    }
    TOF = beamEvent.GetTOF();
    break;
  }
  if(TOF >= fTOFElectronCuts.at(0) && TOF <= fTOFElectronCuts.at(1))
  {
    result.electron = true;
  }
  if(TOF >= fTOFMuonCuts.at(0) && TOF <= fTOFMuonCuts.at(1))
  {
    result.muon = true;
  }
  if(TOF >= fTOFPionCuts.at(0) && TOF <= fTOFPionCuts.at(1))
  {
    result.pion = true;
  }
  if(TOF >= fTOFKaonCuts.at(0) && TOF <= fTOFKaonCuts.at(1))
  {
    result.kaon = true;
  }
  if(TOF >= fTOFProtonCuts.at(0) && TOF <= fTOFProtonCuts.at(1))
  {
    result.proton = true;
  }
  
  return result;
}

protoana::PossibleParticleCands protoana::ProtoDUNEDataUtils::GetBeamlineParticleID(art::Event const & evt, const float beamEnergyGeV) const {
  const auto tof = GetTOFParticleID(evt,beamEnergyGeV);
  const auto chkov = GetCherenkovParticleID(evt,beamEnergyGeV);
  auto result = tof && chkov;
  result.electron = chkov.electron; // disregard tof for electrons
  return result;
}

const std::tuple<double,double,int,int> protoana::ProtoDUNEDataUtils::GetBeamlineVars(art::Event const & evt) const {

  double momentum = -99999.;
  double tof = -99999.;
  int ckov0 = -99999.;
  int ckov1 = -99999.;

  std::vector<art::Ptr<beam::ProtoDUNEBeamEvent>> beamVec;
  auto beamHand = evt.getValidHandle<std::vector<beam::ProtoDUNEBeamEvent>>(fBeamEventTag);
  if(beamHand.isValid())
  {
    art::fill_ptr_vector(beamVec, beamHand);
  }
  for(size_t iBeamEvent=0; iBeamEvent < beamVec.size(); iBeamEvent++)
  {
    const beam::ProtoDUNEBeamEvent& beamEvent = *(beamVec.at(iBeamEvent));
    tof = beamEvent.GetTOF();
    ckov0 = beamEvent.GetCKov0Status();
    ckov1 = beamEvent.GetCKov1Status();

    const std::vector<double> & beamMomenta = beamEvent.GetRecoBeamMomenta();
    if (beamMomenta.size() > 0)
    {
        momentum = beamEvent.GetRecoBeamMomentum(0);
    }
  }
  return std::make_tuple(momentum, tof, ckov0, ckov1);
}

const std::tuple<double,double,int,int,int,double,double,int,int,bool> protoana::ProtoDUNEDataUtils::GetBeamlineVarsAndStatus(art::Event const & evt) const {

  double momentum = -99999.;
  double tof = -99999.;
  int tofChannel = -99999.;
  int ckov0 = -99999.;
  int ckov1 = -99999.;
  double ckov0Pressure = -99999.;
  double ckov1Pressure = -99999.;
  int timingTrigger = -99999.;
  int BITrigger = -99999.;
  bool areBIAndTimingMatched = false;

  std::vector<art::Ptr<beam::ProtoDUNEBeamEvent>> beamVec;
  auto beamHand = evt.getValidHandle<std::vector<beam::ProtoDUNEBeamEvent>>(fBeamEventTag);
  if(beamHand.isValid())
  {
    art::fill_ptr_vector(beamVec, beamHand);
  }
  for(size_t iBeamEvent=0; iBeamEvent < beamVec.size(); iBeamEvent++)
  {
    const beam::ProtoDUNEBeamEvent& beamEvent = *(beamVec.at(iBeamEvent));
    tof = beamEvent.GetTOF();
    tofChannel = beamEvent.GetTOFChan();
    ckov0 = beamEvent.GetCKov0Status();
    ckov1 = beamEvent.GetCKov1Status();
    ckov0Pressure = beamEvent.GetCKov0Pressure();
    ckov1Pressure = beamEvent.GetCKov1Pressure();
    timingTrigger = beamEvent.GetTimingTrigger();
    BITrigger = beamEvent.GetBITrigger();
    areBIAndTimingMatched = beamEvent.CheckIsMatched();

    const std::vector<double> & beamMomenta = beamEvent.GetRecoBeamMomenta();
    if (beamMomenta.size() > 0)
    {
        momentum = beamEvent.GetRecoBeamMomentum(0);
    }
    break;
  }
  return std::make_tuple(momentum, tof, tofChannel, ckov0, ckov1, ckov0Pressure, ckov1Pressure, timingTrigger, BITrigger, areBIAndTimingMatched);
}

