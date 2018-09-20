// DataPrepModule_module.cc

// David Adams
// July 2016
// November 2016 - Add option to group channels.
//
// Module that reads RawData and writes Wire and their associations.
// It uses RawDigitPrepService to build the wires.
//
// It is possible to also write intermediate states. List the states in the fcl
// vector IntermediateStates. Supported values are:
//   extracted - After pedestal subtraction
//   mitigated - After mitigation (e.g. stuck bit interpolation)
//   noiseRemoved - After noise removal
// The states are written to the containers with the same name. To save space
// and time, e.g. in standard production, this vector should be empty.
//
// If the flag is set set, channels are processed in groups specified by
// the ChannelGroupService. This can save memory because the transient ADC
// channel data is deleted after each group is processed.
//
// Configuration parameters:
//             LogLevel - Usual logging level.
//           DigitLabel - Full label for the input digit container, e.g. daq
//             WireName - Name for the output wire container.
//   IntermediateStates - Names of intermediate states to record.
//             DoGroups - Process channels in groups.

#include "art/Framework/Core/ModuleMacros.h" 
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Principal/Event.h" 
#include "art/Framework/Services/Registry/ServiceHandle.h" 
#include "dune/Protodune/singlephase/RawDecoding/data/RDStatus.h"
#include "lardataobj/RawData/RawDigit.h"
#include "lardataobj/RecoBase/Wire.h"
#include "lardata/Utilities/AssociationUtil.h"
#include "dune/DuneInterface/RawDigitPrepService.h"
#include "dune/DuneInterface/ChannelGroupService.h"
#include "dune/DuneInterface/Tool/IndexMapTool.h"
#include "dune/DuneCommon/DuneTimeConverter.h"
#include "dune/ArtSupport/DuneToolManager.h"
#include "TTimeStamp.h"
#include "lardataobj/RawData/RDTimeStamp.h"

using std::cout;
using std::endl;
using std::string;
using std::vector;
using std::move;
using art::ServiceHandle;
using art::Timestamp;
using raw::RDStatus;
using recob::Wire;

//**********************************************************************

class DataPrepModule : public art::EDProducer {

public:
    
  // Ctor.
  explicit DataPrepModule(fhicl::ParameterSet const& pset); 

  // Dtor.
  ~DataPrepModule();
    
  // Producer methods.
  void produce(art::Event& evt);
  void beginJob();
  void endJob();
  void reconfigure(fhicl::ParameterSet const& p);
    
private:
    
  // Configuration parameters.
  int m_LogLevel;
  std::string m_DigitLabel;  ///< Full label for the input digit container, e.g. daq:
  std::string m_WireName;    ///< Second field in full label for the output wire container.
  std::vector<std::string> m_IntermediateStates;
  bool m_DoAssns = false;
  bool m_DoGroups = false;
  AdcChannel m_KeepChannelBegin =0;
  AdcChannel m_KeepChannelEnd =0;
  AdcChannelVector m_SkipChannels;
  AdcChannelVector m_KeepFembs;

  // Split label into producer and name: PRODUCER or PRODUCER:NAME
  std::string m_DigitProducer;
  std::string m_DigitName;

  // Accessed services.
  RawDigitPrepService* m_pRawDigitPrepService = nullptr;
  ChannelGroupService* m_pChannelGroupService = nullptr;

  // Tools.
  std::string m_OnlineChannelMapTool;
  std::unique_ptr<IndexMapTool> m_onlineChannelMapTool;

  // Processed event count.
  unsigned int m_nproc =0;

  // Skipped event count.
  unsigned int m_nskip =0;

};

DEFINE_ART_MODULE(DataPrepModule)
  
//**********************************************************************

DataPrepModule::DataPrepModule(fhicl::ParameterSet const& pset) {
  this->reconfigure(pset);
  produces<std::vector<recob::Wire>>(m_WireName);
  if ( m_DoAssns ) {
    produces<art::Assns<raw::RawDigit, recob::Wire>>(m_WireName);
  }
  for ( string sname : m_IntermediateStates ) {
    produces<std::vector<recob::Wire>>(sname);
  }
}
  
//**********************************************************************

DataPrepModule::~DataPrepModule() { }

//**********************************************************************

void DataPrepModule::reconfigure(fhicl::ParameterSet const& pset) {
  const string myname = "DataPrepModule::reconfigure: ";
  m_LogLevel   = pset.get<int>("LogLevel");
  m_DigitLabel = pset.get<std::string>("DigitLabel", "daq");
  m_WireName   = pset.get<std::string>("WireName", "");
  m_DoAssns    = pset.get<bool>("DoAssns");
  m_DoGroups   = pset.get<bool>("DoGroups");
  m_IntermediateStates = pset.get<vector<string>>("IntermediateStates");
  pset.get_if_present<AdcChannel>("KeepChannelBegin", m_KeepChannelBegin);
  pset.get_if_present<AdcChannel>("KeepChannelEnd", m_KeepChannelEnd);
  pset.get_if_present<AdcChannelVector>("SkipChannels", m_SkipChannels);
  pset.get_if_present<AdcChannelVector>("KeepFembs", m_KeepFembs);
  pset.get_if_present<std::string>("OnlineChannelMapTool", m_OnlineChannelMapTool);

  size_t ipos = m_DigitLabel.find(":");
  if ( ipos == std::string::npos ) {
    m_DigitProducer = m_DigitLabel;
  } else {
    m_DigitProducer = m_DigitLabel.substr(0, ipos);
    m_DigitName = m_DigitLabel.substr(ipos + 1);
  }

  m_pRawDigitPrepService = &*ServiceHandle<RawDigitPrepService>();
  if ( m_DoGroups ) m_pChannelGroupService = &*ServiceHandle<ChannelGroupService>();

  if ( m_OnlineChannelMapTool.size() ) {
    DuneToolManager* ptm = DuneToolManager::instance();
    m_onlineChannelMapTool = ptm->getPrivate<IndexMapTool>(m_OnlineChannelMapTool);
  }

  if ( m_LogLevel >= 1 ) {
    cout << myname << "             LogLevel: " << m_LogLevel << endl;
    cout << myname << "           DigitLabel: " << m_DigitLabel << " (" << m_DigitProducer
                   << ", " << m_DigitName << ")" << endl;
    cout << myname << "             WireName: " << m_WireName << endl;
    cout << myname << "              DoAssns: " << m_DoAssns << endl;
    cout << myname << "             DoGroups: " << m_DoGroups << endl;
    cout << myname << "   IntermediateStates: [";
    int count = 0;
    for ( string sname : m_IntermediateStates ) cout << (count++ == 0 ? "" : " ") << sname;
    cout << "]" << endl;
    cout << myname << "  OnlineChannelMapTool: " << m_OnlineChannelMapTool << endl;
    cout << myname << "      KeepChannelBegin: " << m_KeepChannelBegin << endl;
    cout << myname << "        KeepChannelEnd: " << m_KeepChannelEnd << endl;
    cout << myname << "          SkipChannels: " << "{";
    bool first = true;
    for ( AdcChannel ich : m_SkipChannels ) {
      if ( first ) first = false;
      else cout << ", ";
      cout << ich;
    }
    cout << "}" << endl;
    cout << myname << "             KeepFembs: " << "{";
    first = true;
    for ( AdcChannel ifmb : m_KeepFembs ) {
      if ( first ) first = false;
      else cout << ", ";
      cout << ifmb;
    }
    cout << "}" << endl;
  }
}

//**********************************************************************

void DataPrepModule::beginJob() {
  const string myname = "DataPrepModule::beginJob: ";
  if ( m_LogLevel >= 2 ) cout << myname << "Starting job." << endl;
  m_nproc = 0;
  m_nskip = 0;
}

//**********************************************************************

void DataPrepModule::endJob() {
  const string myname = "DataPrepModule::endJob: ";
  if ( m_LogLevel >= 1 ) {
    cout << myname << "# events processed: " << m_nproc << endl;
    cout << myname << "  # events skipped: " << m_nskip << endl;
  }
}
  
//**********************************************************************

void DataPrepModule::produce(art::Event& evt) {      
  const string myname = "DataPrepModule::produce: ";

  // Control flags.
  bool skipAllEvents = false;
  bool skipEventsWithCorruptDataDropped = false;

  // Fetch the event time.
  Timestamp beginTime = evt.time();

  // Fetch the timing clock.
  string m_TimingProducer = "timingrawdecoder";
  AdcLongIndex timingClock = 0;
  if ( true ) {
    art::Handle<std::vector<raw::RDTimeStamp>> htims;
    //evt.getByLabel(m_DigitProducer, m_DigitName, htims);
    evt.getByLabel("timingrawdecoder", "daq", htims);
    if ( ! htims.isValid() ) {
      cout << myname << "WARNING: Timing clocks product not found." << endl;
    } else if (  htims->size() != 1 ) {
      cout << myname << "WARNING: Unexpected timing clocks size: " << htims->size() << endl;
      for ( unsigned int itim=0; itim<htims->size() && itim<50; ++itim ) {
        cout << myname << "  " << htims->at(itim).GetTimeStamp() << endl;
      }
    } else {
      const raw::RDTimeStamp& tim = htims->at(0);
      cout << myname << "Timing clock: " << tim.GetTimeStamp() << endl;
      timingClock = tim.GetTimeStamp();
    }
  }

  // Read the raw digit status.
  art::Handle<std::vector<raw::RDStatus>> hrdstats;
  evt.getByLabel(m_DigitProducer, m_DigitName, hrdstats);
  string srdstat;
  bool skipEvent = skipAllEvents;
  if ( ! hrdstats.isValid() ) {
    cout << myname << "WARNING: Raw data status product not found." << endl;
  } else {
    if ( hrdstats->size() != 1 ) {
      cout << myname << "WARNING: Unexpected raw data status size: " << hrdstats->size() << endl;
    }
    const RDStatus rdstat = hrdstats->at(0);
    if ( false ) {
      cout << myname << "Raw data status: " << rdstat.GetStatWord();
      if ( rdstat.GetCorruptDataDroppedFlag() ) cout << " (Corrupt data was dropped.)";
      if ( rdstat.GetCorruptDataKeptFlag() ) cout << " (Corrupt data was retained.)";
      cout << endl;
    }
    srdstat = "rdstat=" + std::to_string(rdstat.GetStatWord());
    skipEvent |= skipEventsWithCorruptDataDropped && rdstat.GetCorruptDataDroppedFlag();
  }

  // Read in the digits. 
  if ( m_LogLevel >= 2 ) {
    cout << myname << "Run " << evt.run();
    if ( evt.subRun() ) cout << "-" << evt.subRun();
    cout << ", event " << evt.event();
    if ( srdstat.size() ) cout << ", " << srdstat;
    cout << ", nproc=" << m_nproc;
    if ( m_nskip ) cout << ", nskip=" << m_nskip;
    cout << endl;
    if ( m_LogLevel >= 3 ) cout << myname << "Reading raw digits for producer, name: " << m_DigitProducer << ", " << m_DigitName << endl;
    // July 2018. ProtoDUNE real data has zero in high field and unix time in low field.
    if ( beginTime.timeLow() == 0 ) {
      cout << myname << "Sim data event time: " << DuneTimeConverter::toString(beginTime) << endl;
    } else {
      unsigned int itim = beginTime.timeHigh();
      unsigned int itimrem = 0;
      if ( itim == 0 ) {
        itimrem = itim;
        itim = beginTime.timeLow();
      }
      TTimeStamp rtim(itim, itimrem);
      string stim = string(rtim.AsString("s")) + " UTC";
      cout << myname << "Real data event time: " << itim << " (" << stim << ")" << endl;
    }
  }
  if ( m_LogLevel >= 3 ) {
    cout << myname << "Event time high, low: " << beginTime.timeHigh() << ", " << beginTime.timeLow() << endl;
  }

  // Read in the digits. 
  art::Handle<std::vector<raw::RawDigit>> hdigits;
  evt.getByLabel(m_DigitProducer, m_DigitName, hdigits);
  if ( m_LogLevel >= 3 ) {
    cout << myname << "# digits read: " << hdigits->size() << endl;
  }
  if ( hdigits->size() == 0 ) mf::LogWarning("DataPrepModule") << "Input digit container is empty";

  // Create the container to hold the output wires.
  std::unique_ptr<std::vector<recob::Wire>> pwires(new std::vector<recob::Wire>);
  pwires->reserve(hdigits->size());

  // Create the association container.
  std::unique_ptr<art::Assns<raw::RawDigit,recob::Wire>> passns(new art::Assns<raw::RawDigit,recob::Wire>);

  // If status was bad, skip this event.
  // We store empty results to avoid exception.
  // We have to have read digits to store those results (yech).
  if ( skipEvent ) {
    cout << myname << "Skipping event with " << srdstat << endl;
    evt.put(std::move(pwires), m_WireName);
    if ( m_DoAssns ) evt.put(std::move(passns), m_WireName);
    ++m_nskip;
    return;
  }
  
  // Prepare the intermediate state cache.
  // Note that transient data is retained between groups and so most of the memory saving
  // of groups is lost if  intermediate states are recorded.
  WiredAdcChannelDataMap* pintStates = nullptr;
  if ( m_IntermediateStates.size() ) {
    pintStates = new WiredAdcChannelDataMap(m_IntermediateStates, hdigits->size());
  }

  // Create the transient data map and copy the digits there.
  AdcChannelDataMap fulldatamap;
  bool checkKeep = m_KeepChannelEnd > m_KeepChannelBegin;
  unsigned int nkeep = 0;
  unsigned int nskip = 0;
  unsigned int ndigi = hdigits->size();
  for ( unsigned int idig=0; idig<ndigi; ++idig ) {
    const raw::RawDigit& dig = (*hdigits)[idig];
    AdcChannel chan = dig.Channel();
    if ( checkKeep ) {
      if ( chan < m_KeepChannelBegin || chan >= m_KeepChannelEnd ) {
        ++nskip;
        continue;
      }
    }
    if ( std::find(m_SkipChannels.begin(), m_SkipChannels.end(), chan) != m_SkipChannels.end() ) {
      ++nskip;
      continue;
    }
    if ( fulldatamap.find(chan) != fulldatamap.end() ) {
      cout << myname << "WARNING: Skipping duplicate channel " << chan << "." << endl;
      ++nskip;
      continue;
    }
    // Fetch the online ID.
    bool haveFemb = false;
    AdcChannel fembID = -1;
    AdcChannel fembChannel = -1;
    if ( m_onlineChannelMapTool ) {
      unsigned int ichOn = m_onlineChannelMapTool->get(chan);
      if ( ichOn != IndexMapTool::badIndex() ) {
        fembID = ichOn/128;
        fembChannel = ichOn % 128;
        haveFemb = true;
      }
    }
    if ( m_KeepFembs.size() ) {
      if ( find(m_KeepFembs.begin(), m_KeepFembs.end(), fembID) == m_KeepFembs.end() ) {
        continue;
        ++nskip;
      }
    }
    // Build the channel data.
    AdcChannelData& acd = fulldatamap[chan];
    acd.run = evt.run();
    acd.subRun = evt.subRun();
    acd.event = evt.event();
    acd.channel = chan;
    acd.digitIndex = idig;
    acd.digit = &dig;
    if ( haveFemb ) {
      acd.fembID = fembID;
      acd.fembChannel = fembChannel;
    }
    acd.triggerClock = timingClock;
    acd.metadata["ndigi"] = ndigi;
    ++nkeep;
  }

  // Create a vector of data maps with an entry for each group.
  unsigned int nproc = 0;
  vector<AdcChannelDataMap> datamaps;
  if ( m_DoGroups ) {
    if ( m_pChannelGroupService == nullptr ) {
      mf::LogError("DataPrepModule") << "Channel group service not found." << endl;
      return;
    }
    unsigned int ngrp = m_pChannelGroupService->size();
    for ( unsigned int igrp=0; igrp<ngrp; ++igrp ) {
      datamaps.emplace_back();
      AdcChannelDataMap& datamap = datamaps.back();
      for ( AdcChannel chan : m_pChannelGroupService->channels(igrp) ) {
        if ( fulldatamap.find(chan) != fulldatamap.end() ) {
          datamap.emplace(chan, move(fulldatamap[chan]));
          ++nproc;
        }
      }
    }
  } else {
    datamaps.emplace_back(move(fulldatamap));
    nproc = datamaps.front().size();
  }
  if ( m_LogLevel >= 2 ) {
    cout << myname << "              # input digits: " << ndigi << endl;
    cout << myname << "         # channels selected: " << nkeep << endl;
    cout << myname << "          # channels skipped: " << nskip << endl;
    cout << myname << "  # channels to be processed: " << nproc << endl;
  }

  for ( AdcChannelDataMap& datamap : datamaps ) {

    if ( datamap.size() == 0 ) continue;

    // Use the data preparation service to build the wires and intermediate states.
    int rstat = m_pRawDigitPrepService->prepare(datamap, pwires.get(), pintStates);
    if ( rstat != 0 ) mf::LogWarning("DataPrepModule") << "Data preparation service returned error " << rstat;

    // Build associations between wires and digits.
    if ( m_DoAssns ) {
      for ( const AdcChannelDataMap::value_type& iacd : datamap ) {
        const AdcChannelData& acd = iacd.second;
        AdcIndex idig = acd.digitIndex;
        if ( idig == AdcChannelData::badIndex )
          throw art::Exception(art::errors::ProductRegistrationFailure) << "Digit index is not set.";
        AdcIndex iwir = acd.wireIndex;
        if ( iwir == AdcChannelData::badIndex ) continue;
        art::Ptr<raw::RawDigit> pdig(hdigits, idig);
        bool success = util::CreateAssn(*this, evt, *pwires, pdig, *passns, m_WireName, iwir);
        if ( !success ) throw art::Exception(art::errors::ProductRegistrationFailure)
                              << "Can't associate wire " << iwir << " with raw digit " << idig;
      }
    }

    // Delete the entries from the current channel map.
    // This is an easy way to clear the transient data.
    datamap.erase(datamap.begin(), datamap.end());

  }  // end loop over groups

  if ( m_LogLevel >= 2 ) {
    cout << myname << "Created wire count: " << pwires->size() << endl;
  }
  if ( pwires->size() == 0 ) mf::LogWarning("DataPrepModule") << "No wires made for this event.";

  // Record wires and associations in the event.
  evt.put(std::move(pwires), m_WireName);
  if ( m_DoAssns ) {
    evt.put(std::move(passns), m_WireName);
  }

  // Record intermediate state wires.
  for ( string sname : m_IntermediateStates ) {
    vector<Wire>* pintWires = pintStates->wires[sname];
    if ( pintWires == nullptr ) {
      cout << myname << "WARNING: Wires not found for state " << sname << "." << endl;
      continue;
    }
    if ( m_LogLevel >=2 ) {
      cout << myname << "Recording intermediate state " << sname << "  with "
                     << pintWires->size() << " channels." << endl;
    }
    if ( m_LogLevel >=3 ) {
      for ( const Wire& wire : *pintWires ) {
        cout << myname << "   Channel " << wire.Channel() << " has "
                       << wire.SignalROI().n_ranges() << " ROI and "
                       << wire.SignalROI().count() << "/" << wire.SignalROI().size()
                       << " ticks." << endl;
      }
    }
    std::unique_ptr<std::vector<recob::Wire>> pintWiresWrapped(pintWires);
    evt.put(std::move(pintWiresWrapped), sname);
  }

  ++m_nproc;
  return;
}
