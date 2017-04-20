/////////////////////////////////////////////////////////////////////////
// Class:       OverlayRawDataDetailDUNE35t
// Module Type: producer
// File:        OverlayRawDataDetailDUNE35t_module.cc
//
// This borrows a lot from the microboone mixing module:
//      OverlayRawDataMicroBooNE_module.cc
////////////////////////////////////////////////////////////////////////

/*

   TODO: Implement mixing for OpDetWaveforms in 35t (I have no idea how). What about AuxDet stuff???

 */


#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Modules/MixFilter.h"
#include "art/Framework/IO/ProductMix/MixHelper.h"
#include "art/Framework/Core/PtrRemapper.h"
#include "art/Persistency/Common/CollectionUtilities.h"
#include "art/Framework/Services/System/FileCatalogMetadata.h"
#include "canvas/Utilities/InputTag.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "IFDH_service.h"

#include <memory>
#include <string>
#include <vector>
#include <unordered_map>
#include <exception>
#include <sstream>
#include <unistd.h>

#include "larcore/Geometry/Geometry.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"

#include "larevt/CalibrationDBI/Interface/ChannelStatusService.h"
#include "larevt/CalibrationDBI/Interface/ChannelStatusProvider.h"

#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/MCParticle.h"

#include "lardataobj/Simulation/SimChannel.h"
#include "lardataobj/Simulation/SimPhotons.h"
#include "lardataobj/Simulation/AuxDetSimChannel.h"

#include "lardataobj/MCBase/MCTrack.h"
#include "lardataobj/MCBase/MCShower.h"

#include "DataOverlay/RawDigitMixer.h"
#include "lardataobj/RawData/RawDigit.h"
#include "lardataobj/RawData/ExternalTrigger.h"
#include "larcoreobj/SimpleTypesAndConstants/RawTypes.h"

#include "TFile.h"
#include "TTreeReader.h"

//#include "DataOverlay/OpDetWaveformMixer.h"
//#include "lardataobj/RawData/OpDetWaveform.h"

#include "DataOverlayProducts/EventMixingSummary.h"

namespace mix {
  class OverlayRawDataDetailDUNE35t;
  typedef art::MixFilter<OverlayRawDataDetailDUNE35t> OverlayRawDataDUNE35t;
}

class mix::OverlayRawDataDetailDUNE35t : public boost::noncopyable {
public:

  OverlayRawDataDetailDUNE35t(fhicl::ParameterSet const& p,
                              art::MixHelper &helper);
  ~OverlayRawDataDetailDUNE35t();

  void startEvent(const art::Event&);  //called at the start of every event
  void finalizeEvent(art::Event &);    //called at the end of every event

  size_t nSecondaries() {
    return fEventsToMix;
  }

  void processEventIDs(art::EventIDSequence const& seq); //bookkepping for event IDs

  // Mixing Functions

  // For now, allow exactly one  input. Assume MC inputs have been merged
  // previously and one detsim output created if needed. This could be changed
  // but would require mixing functions for MC here.

  //a lot of MC collections are just simple copies of the collections...
  template<typename T>
  bool MixSimpleCopy( std::vector< std::vector<T> const*> const& inputs,
                      std::vector< T > & output,
                      art::PtrRemapper const &);

  bool MixRawDigits( std::vector< std::vector<raw::RawDigit> const* > const& inputs,
                     std::vector<raw::RawDigit> & output,
                     art::PtrRemapper const &);

  //bool MixOpDetWaveforms_HighGain( std::vector< std::vector<raw::OpDetWaveform> const* > const& inputs,
  //				   std::vector<raw::OpDetWaveform> & output,
  //				   art::PtrRemapper const &);
  //bool MixOpDetWaveforms_LowGain( std::vector< std::vector<raw::OpDetWaveform> const* > const& inputs,
  //				  std::vector<raw::OpDetWaveform> & output,
  //				  art::PtrRemapper const &);

  // Choose mix file.

  std::string getMixFile();


private:

  // Declare member data here.
  RawDigitMixer fRDMixer;
  //OpDetWaveformMixer         fODMixer;

  fhicl::ParameterSet fpset;
  short fDefaultRawDigitSatPoint;
  //short                fDefaultOpDetSatPoint;
  //size_t               fOpDetMinSampleSize;
  bool fInputFileIsData;

  std::string fRawDigitDataModuleLabel;
  //std::string          fOpDetDataModuleLabel;
  std::string fRawDigitMCModuleLabel;
  //std::string          fOpDetMCModuleLabel;

  std::string fRawDigitInputSourceModuleLabel;
  //std::string          fOpDetInputSourceModuleLabel;
  std::string fRawDigitMixerSourceModuleLabel;
  //std::string          fOpDetMixerSourceModuleLabel;

  std::string fG4InputModuleLabel;
  std::string fGeneratorInputModuleLabel;
  std::string fTriggerInputModuleLabel;

  bool fDoMCReco;
  std::string fMCRecoInputModuleLabel;

  size_t fEventsToMix;
  float fDefaultMCRawDigitScale;
  //float                fDefaultMCOpDetScale;

  bool fForceStuckBitRetention;
  size_t fDataMixStartTick;
  size_t fDataMixEndTick;
  size_t fMCMixStartTick;
  size_t fMCMixEndTick;

  std::string fSamDefname;
  std::string fSamProject;
  std::string fSamStation;
  std::string fSamAppFamily;
  std::string fSamAppName;
  std::string fSamAppVersion;
  std::string fSamUser;
  std::string fSamDescription;
  int fSamFileLimit;
  std::string fSamSchema;

  std::string fSamProjectURI;
  std::string fSamProcessID;
  std::string fSamCurrentFileURI;
  std::string fSamCurrentFileName;

  std::string fChannelGainFile;

  art::Handle< std::vector<raw::RawDigit> > inputDigitHandle;
  //art::Handle< std::vector<raw::OpDetWaveform> > inputOpDetHandle_HighGain;
  //art::Handle< std::vector<raw::OpDetWaveform> > inputOpDetHandle_LowGain;

  void GenerateMCRawDigitScaleMap(std::vector<raw::RawDigit> const&);
  std::unordered_map<raw::ChannelID_t,float> fMCRawDigitScaleMap;

  //void GenerateMCOpDetHighGainScaleMap(std::vector<raw::OpDetWaveform> const&);
  //std::unordered_map<raw::Channel_t,float> fMCOpDetHighGainScaleMap;

  //void GenerateMCOpDetLowGainScaleMap(std::vector<raw::OpDetWaveform> const&);
  //std::unordered_map<raw::Channel_t,float> fMCOpDetLowGainScaleMap;

  std::unique_ptr< std::vector<dunemix::EventMixingSummary> > fEventMixingSummary;

};


mix::OverlayRawDataDetailDUNE35t::OverlayRawDataDetailDUNE35t(fhicl::ParameterSet const& p,
                                                              art::MixHelper &helper)
  :
    fRDMixer(false), //print warnings turned off
    //fODMixer(false), //print warnings turned off
    fpset(p.get<fhicl::ParameterSet>("detail")),
    fDefaultRawDigitSatPoint(fpset.get<short>("DefaultRawDigitSaturationPoint",4096)),
    //fDefaultOpDetSatPoint(fpset.get<short>("DefaultOpDetSaturationPoint",4096)),
    //fOpDetMinSampleSize(fpset.get<size_t>("OpDetMinSampleSize",100)),
    fInputFileIsData(fpset.get<bool>("InputFileIsData")),
    fRawDigitDataModuleLabel(fpset.get<std::string>("RawDigitDataModuleLabel")),
    //fOpDetDataModuleLabel(fpset.get<std::string>("OpDetMCModuleLabel")),
    fRawDigitMCModuleLabel(fpset.get<std::string>("RawDigitMCModuleLabel")),
    //fOpDetMCModuleLabel(fpset.get<std::string>("OpDetMCModuleLabel")),
    fEventsToMix(fpset.get<size_t>("EventsToMix",1)),
    fDefaultMCRawDigitScale(fpset.get<float>("DefaultMCRawDigitScale",1)),
    //fDefaultMCOpDetScale(fpset.get<float>("DefaultMCOpDetScale",1)),

    fForceStuckBitRetention(fpset.get<bool>("ForceStuckBitRetention",false)),
    fDataMixStartTick(fpset.get<size_t>("DataMixStartTick",9800)),
    fDataMixEndTick(fpset.get<size_t>("DataMixEndTick",15000)),
    fMCMixStartTick(fpset.get<size_t>("MCMixStartTick",0)),
    fMCMixEndTick(fpset.get<size_t>("MCMixEndTick",5200)),

// Get sam related parameters.
// These parameters should normally be set by the work flow.
// Usually, the only ones that should need to be set are "SamDefname" and "SamProject."

    fSamDefname(fpset.get<std::string>("SamDefname", "")),
    fSamProject(fpset.get<std::string>("SamProject", "")),
    fSamStation(fpset.get<std::string>("SamStation", "")),
    fSamAppFamily(fpset.get<std::string>("SamAppFamily", "art")),
    fSamAppName(fpset.get<std::string>("SamAppName", "mix")),
    fSamAppVersion(fpset.get<std::string>("SamAppVersion", "1")),
    fSamUser(fpset.get<std::string>("SamUser", "")),
    fSamDescription(fpset.get<std::string>("SamDescription", "")),
    fSamFileLimit(fpset.get<int>("SamFileLimit", 100)),
    fSamSchema(fpset.get<std::string>("SamSchema", "root")),  // xrootd by default.

    // if desired, give a path to a root file containing a ttree described below in the function
    fChannelGainFile(fpset.get<std::string>("ChannelGainFile","")),

    fEventMixingSummary(nullptr)
{

  if(fEventsToMix!=1) {
      std::stringstream err_str;
      err_str << "ERROR! Really sorry, but we can only do mixing for one collection right now! ";
      err_str << "\nYep. We're gonna throw an exception now. You should change your fcl to set 'EventsToMix' to 1";
      throw cet::exception("OverlayRawDataDUNE35t") << err_str.str() << std::endl;
    }

  if(fInputFileIsData) {
      fRawDigitInputSourceModuleLabel = fRawDigitDataModuleLabel;
      //fOpDetInputSourceModuleLabel    = fOpDetDataModuleLabel;
      fRawDigitMixerSourceModuleLabel = fRawDigitMCModuleLabel;
      //fOpDetMixerSourceModuleLabel    = fOpDetMCModuleLabel;
    }
  else if(!fInputFileIsData) {
      fRawDigitInputSourceModuleLabel = fRawDigitMCModuleLabel;
      //fOpDetInputSourceModuleLabel    = fOpDetMCModuleLabel;
      fRawDigitMixerSourceModuleLabel = fRawDigitDataModuleLabel;
      //fOpDetMixerSourceModuleLabel    = fOpDetDataModuleLabel;
    }

  if(fInputFileIsData) {
      fDoMCReco = fpset.get_if_present<std::string>("MCRecoInputModuleLabel",fMCRecoInputModuleLabel);
      fG4InputModuleLabel = fpset.get<std::string>("G4InputModuleLabel");
      fGeneratorInputModuleLabel = fpset.get<std::string>("GeneratorInputModuleLabel");
      fTriggerInputModuleLabel = fpset.get<std::string>("TriggerInputModuleLabel");

      std::string instance = "MC";

      //MC generator info is a simple copy
      helper.declareMixOp( art::InputTag(fGeneratorInputModuleLabel),
                           instance,
                           &OverlayRawDataDetailDUNE35t::MixSimpleCopy<simb::MCTruth>,
                           *this );

      //Simple copies of G4 SimPhotons, MCParticles, SimChannels, and SimAuxDetChannel
      helper.declareMixOp( art::InputTag(fG4InputModuleLabel),
                           instance,
                           &OverlayRawDataDetailDUNE35t::MixSimpleCopy<simb::MCParticle>,
                           *this );
      //helper.declareMixOp( art::InputTag(fG4InputModuleLabel),
      //                   instance,
      //			 &OverlayRawDataDetailDUNE35t::MixSimpleCopy<sim::SimPhotons>,
      //			 *this );
      helper.declareMixOp( art::InputTag(fG4InputModuleLabel),
                           instance,
                           &OverlayRawDataDetailDUNE35t::MixSimpleCopy<sim::SimChannel>,
                           *this );
      helper.declareMixOp( art::InputTag(fG4InputModuleLabel),
                           instance,
                           &OverlayRawDataDetailDUNE35t::MixSimpleCopy<sim::AuxDetSimChannel>,
                           *this );

      helper.declareMixOp( art::InputTag(fTriggerInputModuleLabel),
                           instance,
                           &OverlayRawDataDetailDUNE35t::MixSimpleCopy<raw::ExternalTrigger>,
                           *this );

      helper.declareMixOp( art::InputTag(fRawDigitMCModuleLabel),
                           instance,
                           &OverlayRawDataDetailDUNE35t::MixSimpleCopy<raw::RawDigit>,
                           *this );

      /*
         //Associations of MCParticles to MCTruth...hopefully a simple copy is enough
         helper.declareMixOp( art::InputTag(fG4InputModuleLabel),
         &OverlayRawDataDetailDUNE35t::MixSimpleCopy
                         < art::Assns<simb::MCTruth,simb::MCParticle,void> >,
         *this );
       */

    } //end if file is input data

  if (!fInputFileIsData)
    {
      helper.declareMixOp( art::InputTag(fRawDigitDataModuleLabel),
                           "DATA",
                           &OverlayRawDataDetailDUNE35t::MixSimpleCopy<raw::RawDigit>,
                           *this );
    }

  helper.declareMixOp( art::InputTag(fRawDigitMixerSourceModuleLabel),
                       "MIXED",
                       &OverlayRawDataDetailDUNE35t::MixRawDigits,
                       *this );

  //helper.declareMixOp( art::InputTag(fOpDetMixerSourceModuleLabel,"OpdetBeamHighGain"),
  //		       &OverlayRawDataDetailDUNE35t::MixOpDetWaveforms_HighGain,
  //		       *this );
  //helper.declareMixOp( art::InputTag(fOpDetMixerSourceModuleLabel,"OpdetBeamLowGain"),
  //		       &OverlayRawDataDetailDUNE35t::MixOpDetWaveforms_LowGain,
  //		       *this );

  //If it produces something on its own, declare it here
  helper.produces< std::vector<dunemix::EventMixingSummary> >();

  // Following block handles case of mix input from sam.

  if(!fSamDefname.empty()) {

      // Register getMixFile method with MixHelper.

      helper.registerSecondaryFileNameProvider(std::bind(&mix::OverlayRawDataDetailDUNE35t::getMixFile, this));

      // Get IFDH art service.

      art::ServiceHandle<ifdh_ns::IFDH> ifdh;

      // Get sam station.
      // If the station was not specified by a fcl parameter, use environment variable
      // $SAM_STATION, or else use a default value of "dune."

      if(fSamStation.empty()) {
          const char* c = getenv("SAM_STATION");
          if(c == 0 || *c == 0)
            c = "dune";
          fSamStation = c;
          std::cout << "Mix SAM: Station = " << fSamStation << std::endl;
        }

      // Find project uri.

      fSamProjectURI = ifdh->findProject(fSamProject, fSamStation);
      std::cout << "Mix SAM: project uri = " << fSamProjectURI << std::endl;
      if(fSamProjectURI.empty())
        throw cet::exception("OverlayRawDataDUNE35t") << "Failed to find project uri.";

      // Get hostname.

      char hostname[256];
      gethostname(hostname, sizeof hostname);

      // Get user.
      // If the user was not specified by a fcl parameter, use environment variable
      // $SAM_USER (this should work on grid), or else use environment variable $LOGNAME.

      if(fSamUser.empty()) {
          const char* c = getenv("SAM_USER");
          if(c == 0 || *c == 0)
            c = getenv("LOGNAME");
          if(c != 0 && *c != 0)
            fSamUser = c;
          std::cout << "Mix SAM: User = " << fSamUser << std::endl;
        }

      // Join project.

      std::cout << "Mix SAM: fSamProjectURI = " << fSamProjectURI << "\n"
                << "Mix SAM: fSamAppName = " << fSamAppName << "\n"
                << "Mix SAM: fSamAppVersion = " << fSamAppVersion << "\n"
                << "Mix SAM: hostname = " << hostname << "\n"
                << "Mix SAM: fSamUser = " << fSamUser << "\n"
                << "Mix SAM: fSamAppFamily = " << fSamAppFamily << "\n"
                << "Mix SAM: fSamDescription = " << fSamDescription << "\n"
                << "Mix SAM: fSamFileLimit = " << fSamFileLimit << "\n"
                << "Mix SAM: fSamSchema = " << fSamSchema << "\n";

      fSamProcessID = ifdh->establishProcess(fSamProjectURI,
                                             fSamAppName,
                                             fSamAppVersion,
                                             hostname,
                                             fSamUser,
                                             fSamAppFamily,
                                             fSamDescription,
                                             fSamFileLimit,
                                             fSamSchema);
      mf::LogInfo("OverlayRawDigitDUNE35t") << "Overlay sam definition: " << fSamDefname << "\n"
                                            << "Overlay sam project: " << fSamProject;

      std::cout << "Mix SAM: process id = " << fSamProcessID << std::endl;
      if(fSamProcessID.empty())
        throw cet::exception("OverlayRawDataDUNE35t") << "Failed to start sam process.";
    }
}

// Destructor.
mix::OverlayRawDataDetailDUNE35t::~OverlayRawDataDetailDUNE35t()
{
  if(!fSamProcessID.empty()) {

      // Get IFDH art service.

      art::ServiceHandle<ifdh_ns::IFDH> ifdh;

      // Mark current file as consumed.

      if(!fSamCurrentFileName.empty()) {
          ifdh->updateFileStatus(fSamProjectURI,
                                 fSamProcessID,
                                 fSamCurrentFileName,
                                 "consumed");
          std::cout << "Mix SAM: File " << fSamCurrentFileName << " status changed to consumed." << std::endl;
        }

      // Stop process.

      ifdh->endProcess(fSamProjectURI, fSamProcessID);
      std::cout << "Mix SAM: End process." << std::endl;
    }
}

//Initialize for each event
void mix::OverlayRawDataDetailDUNE35t::startEvent(const art::Event& event) {

  if(!( (event.isRealData() && fInputFileIsData) || (!event.isRealData() && !fInputFileIsData)))
    throw cet::exception("OverlayRawDataDUNE35t") << "Input file claimed to be data/not data, but it's not." << std::endl; ;


  event.getByLabel(fRawDigitInputSourceModuleLabel,inputDigitHandle);
  if(!inputDigitHandle.isValid())
    throw cet::exception("OverlayRawDataDUNE35t") << "Bad input digit handle." << std::endl; ;
  //fRDMixer.SetSaturationPoint(fDefaultRawDigitSatPoint);

  fRDMixer.SetStuckBitRetentionMethod(fForceStuckBitRetention);
  fRDMixer.SetDataMixTicks(fDataMixStartTick,fDataMixEndTick);
  fRDMixer.SetMCMixTicks(fMCMixStartTick,fMCMixEndTick);

  //event.getByLabel(fOpDetInputSourceModuleLabel,"OpdetBeamLowGain",inputOpDetHandle_LowGain);
  //if(!inputOpDetHandle_LowGain.isValid())
  //  throw cet::exception("OverlayRawDataDUNE35t") << "Bad input opdet lowgain handle." << std::endl;;

  //event.getByLabel(fOpDetInputSourceModuleLabel,"OpdetBeamHighGain",inputOpDetHandle_HighGain);
  //if(!inputOpDetHandle_HighGain.isValid())
  //  throw cet::exception("OverlayRawDataDUNE35t") << "Bad input opdet highgain handle." << std::endl;;

  //fODMixer.SetSaturationPoint(fDefaultOpDetSatPoint);
  //fODMixer.SetMinSampleSize(fOpDetMinSampleSize);


  fEventMixingSummary.reset(new std::vector<dunemix::EventMixingSummary>);
}

//For each of the mixed in events...bookkepping for event IDs
void mix::OverlayRawDataDetailDUNE35t::processEventIDs(art::EventIDSequence const& seq){
  for (auto const& id : seq)
    fEventMixingSummary->emplace_back(id.event(),id.subRun(),id.run());
}


//End each event
void mix::OverlayRawDataDetailDUNE35t::finalizeEvent(art::Event& event) {
  event.put(std::move(fEventMixingSummary));
}

template<typename T>
bool mix::OverlayRawDataDetailDUNE35t::MixSimpleCopy( std::vector< std::vector<T> const*> const& inputs,
                                                      std::vector< T > & output,
                                                      art::PtrRemapper const &){
  art::flattenCollections(inputs,output);
  return true;
}

void mix::OverlayRawDataDetailDUNE35t::GenerateMCRawDigitScaleMap(std::vector<raw::RawDigit> const& dataDigitVector){

  // read in data file which contains a TTree with all of the channel gain corrections,
  // calculated using normal (unscaled) MC and real data
  // NOTE: implementing this channel gain factors must be done in a two-step process:
  //    1. Run full data overlay and reconstruction (without additional gain corrections)
  //    2. Calculate channel gain correction factors, then re-run the entire dataoverlay
  //       process (and necessarily, the rest of the reco/ana as well) with this method
  //       implementing the new gains
  std::map<raw::ChannelID_t,double> channelgains;
  art::ServiceHandle<geo::Geometry> fGeom;
  try
    {
      TFile * gainsfile = TFile::Open(fChannelGainFile.c_str(),"READ");
      if (gainsfile && !gainsfile->IsZombie())
        {
          TTreeReader reader("channelgains",gainsfile);
          TTreeReaderValue<Double_t> changain(reader,"changain");
          TTreeReaderValue<Int_t> channel(reader,"channel");
          TTreeReaderValue<Double_t> datampverr(reader,"datampverr");
          TTreeReaderValue<Double_t> datawidtherr(reader,"datawidtherr");
          TTreeReaderValue<Double_t> simmpverr(reader,"simmpverr");
          TTreeReaderValue<Double_t> simwidtherr(reader,"simwidtherr");
          while (reader.Next())
            {
              raw::ChannelID_t chanid = *channel;
              double gain;
              if (*datampverr < 200 && *datawidtherr < 400 && *simmpverr < 200 && *simwidtherr < 400)
                {
                  gain = *changain;
                }
              else
                {
                  gain = 1.0;
                }
              channelgains[chanid] = gain;
            }
        }
    }
  catch (...)
    {
      std::cout << "Can't open gains file, setting channelgains all equal to 1" << std::endl;
      for (raw::ChannelID_t chan = 0; chan < fGeom->Nchannels(); ++chan)
        {
          channelgains[chan] = 1.0;
        }
    }

  //right now, assume the number of channels is the number in the collection
  //and, loop through the channels one by one to get the right channel number
  //note: we will put here access to the channel database to determine dead channels
  fMCRawDigitScaleMap.clear();

  const lariov::ChannelStatusProvider& chanStatus = art::ServiceHandle<lariov::ChannelStatusService>()->GetProvider();

  for(auto const& d : dataDigitVector) {
      if(chanStatus.IsBad(d.Channel()))
        {
          if (fGeom->SignalType(d.Channel()) == geo::kCollection) std::cout << "Channel " << d.Channel() << " IN BAD CHANNEL LIST. Scale=0.0" << std::endl;
          fMCRawDigitScaleMap[d.Channel()] = 0.0;
        }
      else if (channelgains.find(d.Channel()) == channelgains.end())
        {
          if (fGeom->SignalType(d.Channel()) == geo::kCollection) std::cout << "Channel " << d.Channel() << " NOT IN channelgains. Scale=" << fDefaultMCRawDigitScale << std::endl;
          fMCRawDigitScaleMap[d.Channel()] = fDefaultMCRawDigitScale;
        }
      else
        {
          if (fGeom->SignalType(d.Channel()) == geo::kCollection) std::cout << "Channel " << d.Channel() << " GOOD. Scale=" << fDefaultMCRawDigitScale * channelgains[d.Channel()] << std::endl;
          fMCRawDigitScaleMap[d.Channel()] = fDefaultMCRawDigitScale * channelgains[d.Channel()];
        }
    }
}

bool mix::OverlayRawDataDetailDUNE35t::MixRawDigits( std::vector< std::vector<raw::RawDigit> const* > const& inputs,
                                                     std::vector<raw::RawDigit> & output,
                                                     art::PtrRemapper const & remap) {

  //make sure we only have two collections for now
  if(inputs.size()!=fEventsToMix || (inputs.size()!=1 && !fInputFileIsData)) {
      std::stringstream err_str;
      err_str << "ERROR! We have the wrong number of collections of raw digits we are adding! " << inputs.size();
      throw std::runtime_error(err_str.str());
    }

  if(fInputFileIsData) {
      GenerateMCRawDigitScaleMap(*inputDigitHandle);
      fRDMixer.DeclareData(*inputDigitHandle);
      for(auto const& icol : inputs)
        fRDMixer.Mix(*icol,fMCRawDigitScaleMap);
    }
  else if(!fInputFileIsData) {
      GenerateMCRawDigitScaleMap(*(inputs[0]));
      fRDMixer.DeclareData(*(inputs[0]));
      fRDMixer.Mix(*inputDigitHandle,fMCRawDigitScaleMap);
    }

  fRDMixer.FillRawDigitOutput(output);

  return true;
}

/*
   void mix::OverlayRawDataDetailDUNE35t::GenerateMCOpDetHighGainScaleMap(std::vector<raw::OpDetWaveform> const& dataVector){
   //right now, assume the number of channels is the number in the collection
   //and, loop through the channels one by one to get the right channel number
   //note: we will put here access to the channel database to determine dead channels
   fMCOpDetHighGainScaleMap.clear();
   for(auto const& d : dataVector)
    fMCOpDetHighGainScaleMap[d.ChannelNumber()] = fDefaultMCOpDetScale;
   }

   void mix::OverlayRawDataDetailDUNE35t::GenerateMCOpDetLowGainScaleMap(std::vector<raw::OpDetWaveform> const& dataVector){
   //right now, assume the number of channels is the number in the collection
   //and, loop through the channels one by one to get the right channel number
   //note: we will put here access to the channel database to determine dead channels
   fMCOpDetLowGainScaleMap.clear();
   for(auto const& d : dataVector)
    fMCOpDetLowGainScaleMap[d.ChannelNumber()] = fDefaultMCOpDetScale;
   }

   bool mix::OverlayRawDataDetailDUNE35t::MixOpDetWaveforms_HighGain( std::vector< std::vector<raw::OpDetWaveform> const* > const& inputs,
                                                                      std::vector<raw::OpDetWaveform> & output,
                                                                      art::PtrRemapper const & remap) {

   //make sure we only have two collections for now
   if(inputs.size()!=fEventsToMix || (inputs.size()!=1 && !fInputFileIsData)){
    std::stringstream err_str;
    err_str << "ERROR! We have the wrong number of collections of raw digits we are adding! " << inputs.size();
    throw std::runtime_error(err_str.str());
   }


   if(fInputFileIsData){
    GenerateMCOpDetHighGainScaleMap(*inputOpDetHandle_HighGain);
    fODMixer.DeclareData(*inputOpDetHandle_HighGain,output);
    for(auto const& icol : inputs)
      fODMixer.Mix(*icol,fMCOpDetHighGainScaleMap,output);
   }
   else if(!fInputFileIsData){
    GenerateMCOpDetHighGainScaleMap(*(inputs[0]));
    fODMixer.DeclareData(*(inputs[0]),output);
    fODMixer.Mix(*inputOpDetHandle_HighGain,fMCOpDetHighGainScaleMap,output);
   }

   return true;
   }

   bool mix::OverlayRawDataDetailDUNE35t::MixOpDetWaveforms_LowGain( std::vector< std::vector<raw::OpDetWaveform> const* > const& inputs,
                                                                     std::vector<raw::OpDetWaveform> & output,
                                                                     art::PtrRemapper const & remap) {

   //make sure we only have two collections for now
   if(inputs.size()!=fEventsToMix || (inputs.size()!=1 && !fInputFileIsData)){
    std::stringstream err_str;
    err_str << "ERROR! We have the wrong number of collections of raw digits we are adding! " << inputs.size();
    throw std::runtime_error(err_str.str());
   }


   if(fInputFileIsData){
    GenerateMCOpDetLowGainScaleMap(*inputOpDetHandle_LowGain);
    fODMixer.DeclareData(*inputOpDetHandle_LowGain,output);
    for(auto const& icol : inputs)
      fODMixer.Mix(*icol,fMCOpDetLowGainScaleMap,output);
   }
   else if(!fInputFileIsData){
    GenerateMCOpDetLowGainScaleMap(*(inputs[0]));
    fODMixer.DeclareData(*(inputs[0]),output);
    fODMixer.Mix(*inputOpDetHandle_LowGain,fMCOpDetLowGainScaleMap,output);
   }

   return true;
   }
 */

// Return next file to mix.

std::string mix::OverlayRawDataDetailDUNE35t::getMixFile()
{
  std::string result;

  if(!fSamProcessID.empty()) {

      // Get IFDH art service.

      art::ServiceHandle<ifdh_ns::IFDH> ifdh;

      // Update status of current file, if any, to "consumed."

      if(!fSamCurrentFileName.empty()) {
          ifdh->updateFileStatus(fSamProjectURI,
                                 fSamProcessID,
                                 fSamCurrentFileName,
                                 "consumed");

          std::cout << "Mix SAM: File " << fSamCurrentFileName << " status changed to consumed." << std::endl;
          fSamCurrentFileURI = std::string();
          fSamCurrentFileName = std::string();
        }

      // Get next file uri.

      fSamCurrentFileURI = fSamCurrentFileURI = ifdh->getNextFile(fSamProjectURI, fSamProcessID);
      unsigned int n = fSamCurrentFileURI.find_last_of('/') + 1;
      fSamCurrentFileName = fSamCurrentFileURI.substr(n);
      std::cout << "Mix SAM: Next file uri = " << fSamCurrentFileURI << std::endl;
      std::cout << "Mix SAM: Next file name = " << fSamCurrentFileName << std::endl;
      mf::LogInfo("OverlayRawDigitDUNE35t") << "Next mix file uri: " << fSamCurrentFileURI << "\n"
                                            << "Next mix file name: " << fSamCurrentFileName;

      // Throw an exception if we didn't get a next file.

      if(fSamCurrentFileURI.empty() || fSamCurrentFileName.empty())
        throw cet::exception("OverlayRawDataDUNE35t") << "Failed to get next mix file.";

      // Here is where we would copy the file to the local node, if that were necessary.
      // Since we are using schema "root" (i.e. xrootd) to stream files, copying the
      // file is not necessary.
      // Note further that we should not update the file status to "transferred" for
      // streaming files, since that can in principle allow the file to be deleted from
      // disk cache.

      // Update metadata.

      art::ServiceHandle<art::FileCatalogMetadata> md;
      md->addMetadataString("mixparent", fSamCurrentFileName);

      // Done.

      result = fSamCurrentFileURI;
    }

  return result;
}


DEFINE_ART_MODULE(mix::OverlayRawDataDUNE35t)
