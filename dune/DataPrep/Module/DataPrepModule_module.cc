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
//             DoGroups - Process channels in groups obtained from ChannelGroupService
//                        if ChannelRanges is empty.
//        ChannelRanges - Process channels in groups corresponding to these range names.
//                        The range for each name is obtained from the tool channelRanges.
//       BeamEventLabel - Label for the BeamEvent data product. If blank, it is not used.

#include "art/Framework/Core/ModuleMacros.h" 
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Principal/Event.h" 
#include "art/Framework/Services/Registry/ServiceHandle.h" 
#include "dune/DuneObj/PDSPTPCDataInterfaceParent.h"
#include "lardataobj/RawData/RawDigit.h"
#include "lardataobj/RecoBase/Wire.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardata/Utilities/AssociationUtil.h"
#include "larevt/CalibrationDBI/Interface/ChannelStatusService.h"
#include "larevt/CalibrationDBI/Interface/ChannelStatusProvider.h"
#include "dune/DuneInterface/Service/RawDigitPrepService.h"
#include "dune/DuneInterface/Service/ChannelGroupService.h"
#include "dune/DuneInterface/Tool/IndexMapTool.h"
#include "dune/DuneInterface/Tool/IndexRangeTool.h"
#include "dune/DuneCommon/Utility/DuneTimeConverter.h"
#include "dune/ArtSupport/DuneToolManager.h"
#include "TTimeStamp.h"
#include "lardataobj/RawData/RDTimeStamp.h"
#include "dune/DuneObj/ProtoDUNEBeamEvent.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include <iomanip>

using std::cout;
using std::endl;
using std::string;
using std::vector;
using std::move;
using std::setw;
using art::ServiceHandle;
using art::Timestamp;
using raw::RDStatus;
using recob::Wire;

//**********************************************************************

class DataPrepModule : public art::EDProducer {

public:
    
  using Index = unsigned int;
  using Name = std::string;
  using NameVector = std::vector<Name>;

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
  Name m_DecoderTool; // Name for the decoder tool
  Name m_DigitLabel;  ///< Full label for the input digit container, e.g. daq:
  Name m_TimeStampName;    // Label for the output RDTimeStamp conainer
  Name m_OutputDigitName;  // Label for the output raw::RawDigit conainer
  Name m_WireName;    ///< Second field in full label for the output wire container.
  std::vector<std::string> m_IntermediateStates;
  bool m_DoAssns = false;
  bool m_DoGroups = false;
  NameVector m_ChannelRanges;
  std::string m_BeamEventLabel;
  AdcChannel m_KeepChannelBegin =0;
  AdcChannel m_KeepChannelEnd =0;
  AdcChannelVector m_SkipChannels;
  AdcChannelVector m_KeepFembs;

  // Split label into producer and name: PRODUCER or PRODUCER:NAME
  std::string m_DigitProducer;
  std::string m_DigitName;

  // Accessed services.
  const lariov::ChannelStatusProvider* m_pChannelStatusProvider;
  RawDigitPrepService* m_pRawDigitPrepService = nullptr;
  ChannelGroupService* m_pChannelGroupService = nullptr;

  // Tools.
  std::string m_OnlineChannelMapTool;

  std::unique_ptr<PDSPTPCDataInterfaceParent> m_pDecoderTool;
  std::unique_ptr<IndexMapTool> m_onlineChannelMapTool;

  // Processed event count.
  unsigned int m_nproc =0;

  // Skipped event count.
  unsigned int m_nskip =0;

};

DEFINE_ART_MODULE(DataPrepModule)
  
//**********************************************************************

DataPrepModule::DataPrepModule(fhicl::ParameterSet const& pset) : EDProducer{pset} {
  const Name myname = "DataPrepModule::ctor: ";
  this->reconfigure(pset);
  produces<std::vector<recob::Wire>>(m_WireName);
  if ( m_DoAssns ) {
    produces<art::Assns<raw::RawDigit, recob::Wire>>(m_WireName);
  }
  for ( string sname : m_IntermediateStates ) {
    if ( m_LogLevel > 0 ) cout << myname << "Module will produce intermediate Wires with name " << sname << endl;
    produces<std::vector<recob::Wire>>(sname);
  }
  //produces<std::vector<raw::RawDigit>>("dataprep");
  if ( m_DecoderTool.size() ) {
    if ( m_TimeStampName.size() ) {
      if ( m_LogLevel > 0 ) {
        cout << myname << "Module will produce RDTimeStamps with name " << m_TimeStampName << endl;
      }
      produces<std::vector<raw::RDTimeStamp>>(m_TimeStampName);
    }
    if ( m_OutputDigitName.size() ) {
      if ( m_LogLevel > 0 ) {
        cout << myname << "Module will produce digits with name " << m_OutputDigitName << endl;
      }
      produces<std::vector<raw::RawDigit>>(m_OutputDigitName);
    }
    //produces<std::vector<raw::RDStatus>>("dataprep");
  }
}
  
//**********************************************************************

DataPrepModule::~DataPrepModule() { }

//**********************************************************************

void DataPrepModule::reconfigure(fhicl::ParameterSet const& pset) {
  const string myname = "DataPrepModule::reconfigure: ";
  m_LogLevel           = pset.get<int>("LogLevel");
  m_DecoderTool        = pset.get<Name>("DecoderTool");
  m_DigitLabel         = pset.get<Name>("DigitLabel", "daq");
  m_TimeStampName      = pset.get<Name>("TimeStampName");
  m_OutputDigitName    = pset.get<Name>("OutputDigitName");
  m_WireName           = pset.get<Name>("WireName", "");
  m_DoAssns            = pset.get<bool>("DoAssns");
  m_DoGroups           = pset.get<bool>("DoGroups");
  m_ChannelRanges      = pset.get<NameVector>("ChannelRanges");
  m_BeamEventLabel     = pset.get<string>("BeamEventLabel");
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

  m_pChannelStatusProvider = &art::ServiceHandle<lariov::ChannelStatusService>()->GetProvider();
  if ( m_pChannelStatusProvider == nullptr ) {
    cout << myname << "WARNING: Channel status provider not found." << endl;
  }

  m_pRawDigitPrepService = &*ServiceHandle<RawDigitPrepService>();
  if ( m_DoGroups ) m_pChannelGroupService = &*ServiceHandle<ChannelGroupService>();

  if ( m_DecoderTool.size() ) {
    DuneToolManager* ptm = DuneToolManager::instance();
    m_pDecoderTool = ptm->getPrivate<PDSPTPCDataInterfaceParent>(m_DecoderTool);
    if ( ! m_pDecoderTool ) {
      cout << myname << "ERROR: Decoder tool not found: " << m_DecoderTool << endl;
    }
  }
  if ( m_pDecoderTool ) {
    cout << myname << "Raw digits will be obtained with decoder tool " << m_DecoderTool << endl;
  } else {
    cout << myname << "Raw digits will be obtained from event data store." << endl;
  }

  if ( m_OnlineChannelMapTool.size() ) {
    DuneToolManager* ptm = DuneToolManager::instance();
    m_onlineChannelMapTool = ptm->getPrivate<IndexMapTool>(m_OnlineChannelMapTool);
  }

  if ( m_LogLevel >= 1 ) {
    cout << myname << "             LogLevel: " << m_LogLevel << endl;
    cout << myname << "          DecoderTool: " << m_DecoderTool << endl;
    cout << myname << "           DigitLabel: " << m_DigitLabel << " (" << m_DigitProducer
                   << ", " << m_DigitName << ")" << endl;
    cout << myname << "        TimeStampName: " << m_TimeStampName << endl;
    cout << myname << "      OutputDigitName: " << m_OutputDigitName << endl;
    cout << myname << "             WireName: " << m_WireName << endl;
    cout << myname << "              DoAssns: " << m_DoAssns << endl;
    cout << myname << "             DoGroups: " << m_DoGroups << endl;
    cout << myname << "        ChannelRanges: [";
    bool first = true;
    for ( Name rnam : m_ChannelRanges ) {
      if ( first ) first = false;
      else cout << ", ";
      cout << rnam;
    }
    cout << "]" << endl;
    cout << myname << "       BeamEventLabel: " << m_BeamEventLabel << endl;
    cout << myname << "   IntermediateStates: [";
    int count = 0;
    for ( string sname : m_IntermediateStates ) cout << (count++ == 0 ? "" : " ") << sname;
    cout << "]" << endl;
    cout << myname << "  OnlineChannelMapTool: " << m_OnlineChannelMapTool << endl;
    cout << myname << "      KeepChannelBegin: " << m_KeepChannelBegin << endl;
    cout << myname << "        KeepChannelEnd: " << m_KeepChannelEnd << endl;
    cout << myname << "          SkipChannels: " << "[";
    first = true;
    for ( AdcChannel ich : m_SkipChannels ) {
      if ( first ) first = false;
      else cout << ", ";
      cout << ich;
    }
    cout << "]" << endl;
    cout << myname << "             KeepFembs: " << "[";
    first = true;
    for ( AdcChannel ifmb : m_KeepFembs ) {
      if ( first ) first = false;
      else cout << ", ";
      cout << ifmb;
    }
    cout << "]" << endl;
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

  // Flag indicating that non-verbose info level messages should be logged.
  bool logInfo = m_LogLevel >= 2;

  // Factor to convert event clock to ticks.
  unsigned long triggerPerTick = 25;

  // Control flags.
  bool skipAllEvents = false;
  bool skipEventsWithCorruptDataDropped = false;

  // Decode the event time.
  Timestamp beginTime = evt.time();
  time_t itim = beginTime.timeHigh();
  int itimrem = beginTime.timeLow();
  // Older protoDUNE data has time in low field.
  if ( itim == 0 && itimrem != 0 ) {
    itimrem = itim;
    itim = beginTime.timeLow();
  }

  // Log event processing header.
  if ( logInfo ) {
    cout << myname << "Run " << evt.run();
    if ( evt.subRun() ) cout << "-" << evt.subRun();
    cout << ", event " << evt.event();
    cout << ", nproc=" << m_nproc;
    if ( m_nskip ) cout << ", nskip=" << m_nskip;
    cout << endl;
    if ( m_LogLevel >= 3 ) cout << myname << "Reading raw digits for producer, name: " << m_DigitProducer << ", " << m_DigitName << endl;
    if ( evt.isRealData() ) {
      TTimeStamp rtim(itim, itimrem);
      string stim = string(rtim.AsString("s")) + " UTC";
      cout << myname << "Real data event time: " << itim << " (" << stim << ")" << endl;
    } else {
      cout << myname << "Sim data event time: " << DuneTimeConverter::toString(beginTime) << endl;
    }
  }
  if ( m_LogLevel >= 3 ) {
    cout << myname << "Event time high, low: " << beginTime.timeHigh() << ", " << beginTime.timeLow() << endl;
  }

  // Fetch the event trigger and timing clock.
  AdcIndex trigFlag = 0;
  AdcLongIndex timingClock = 0;
  using TimeVector  = std::vector<raw::RDTimeStamp>;
  const TimeVector* ptims = nullptr;
  if ( true ) {
    art::InputTag itag1("timingrawdecoder", "daq");
    auto htims = evt.getHandle<TimeVector>(itag1);
    if ( htims ) {
      ptims = &*htims;
      if ( ptims->size() != 1 ) {
        cout << myname << "WARNING: Unexpected timing clocks size: " << ptims->size() << endl;
        for ( unsigned int itim=0; itim<ptims->size() && itim<50; ++itim ) {
          cout << myname << "  " << ptims->at(itim).GetTimeStamp() << endl;
        }
      } else {
        const raw::RDTimeStamp& tim = ptims->at(0);
        if ( logInfo ) cout << myname << "Timing clock: " << tim.GetTimeStamp() << endl;
        timingClock = tim.GetTimeStamp();
        // See https://twiki.cern.ch/twiki/bin/view/CENF/TimingSystemAdvancedOp#Reference_info
        trigFlag = tim.GetFlags();
        if ( m_LogLevel >= 2 ) cout << myname << "Trigger flag: " << trigFlag << " (";
        bool isBeam = trigFlag == 0xc;
        bool isCrt = trigFlag == 13;
        bool isFake = trigFlag >= 0x8 && trigFlag <= 0xb;
        if ( logInfo ) {
          if ( isBeam ) cout << "Beam";
          else if ( isCrt ) cout << "CRT";
          else if ( isFake ) cout << "Fake";
          else cout << "Unexpected";
          cout << ")" << endl;
        }
      }
    } else {
      if ( logInfo ) {
        cout << myname << "WARNING: Event timing clocks product not found." << endl;
      }
    }
  }

  // If the decoder tool is used, use it to retrive the raw digits, their status
  // and the channelclocks.
  bool useDecoderTool = bool(m_pDecoderTool);
  using StatVector  = std::vector<raw::RDStatus>;
  using DigitVector = std::vector<raw::RawDigit>;
  std::unique_ptr<TimeVector>  ptimsFromTool;
  std::unique_ptr<StatVector>  pstatsFromTool;
  std::unique_ptr<DigitVector> pdigitsFromTool;
  if ( useDecoderTool ) {
    if ( useDecoderTool ) {
      ptimsFromTool.reset(new TimeVector);
      pstatsFromTool.reset(new StatVector);
      pdigitsFromTool.reset(new DigitVector);
    }
    std::vector<int> apas = {-1};
    if ( logInfo ) cout << myname << "Fetching digits and clocks with decoder tool." << endl;
    int decodeStat = m_pDecoderTool->
      retrieveDataForSpecifiedAPAs(evt, *pdigitsFromTool.get(), *ptimsFromTool.get(),
                                   *pstatsFromTool.get(), apas);
    if ( m_LogLevel >= 3 ) {    // Decoder tool can return any value for success 
      cout << myname << "WARNING: Decoder tool returned " << decodeStat << endl;
    }
    if ( logInfo ) {
      cout << myname << "  Digit count from tool: " << pdigitsFromTool->size() << endl;
      cout << myname << "  Stats count from tool: " << pstatsFromTool->size() << endl;
      cout << myname << "  Clock count from tool: " << ptimsFromTool->size() << endl;
    }
  }
  
  // Fetch and check the channel clocks from the tool.
  vector<AdcLongIndex> channelClocks;
  vector<ULong64_t> tzeroClockCandidates;   // Candidates for t0 = 0.
  float trigTickOffset = -500.5;
  if ( ptimsFromTool ) {
    using ClockCounter = std::map<ULong64_t, AdcIndex>;
    ClockCounter clockCounts;
    ULong64_t chClock = 0;
    long chClockDiff = 0;
    float tickdiff = 0.0;
    for ( raw::RDTimeStamp chts : *ptimsFromTool ) {
      chClock = chts.GetTimeStamp();
      channelClocks.push_back(chClock);
      chClockDiff = chClock > timingClock ?  (chClock - timingClock)
                                               : -(timingClock - chClock);
      bool nearTrigger = fabs(tickdiff - trigTickOffset) < 1.0;
      tickdiff = chClockDiff/triggerPerTick;
      if ( clockCounts.find(chClock) == clockCounts.end() ) {
        clockCounts[chClock] = 1;
        if ( nearTrigger ) tzeroClockCandidates.push_back(chClock);
      } else {
        ++clockCounts[chClock];
      }
      if ( ! nearTrigger ) {
        if ( m_LogLevel >= 3 ) {
          cout << myname << "WARNING: Channel timing difference: " << chClockDiff
               << " (" << tickdiff << " ticks)." << endl;
        }
      }
    }
    if ( clockCounts.size() > 1 ) {
      if ( logInfo ) {
        cout << myname << "WARNING: Channel clocks are not consistent." << endl;
        cout << myname << "WARNING:     Clock     ticks   count" << endl;
      }
      for ( ClockCounter::value_type iclk : clockCounts ) {
        ULong64_t chClock = iclk.first;
        AdcIndex count = iclk.second;
        long chClockDiff = chClock > timingClock ?  (chClock - timingClock)
                                                 : -(timingClock - chClock);
        float tickdiff = chClockDiff/triggerPerTick;
        if ( logInfo ) {
          cout << myname << "WARNING:" << setw(10) << chClockDiff << setw(10) << tickdiff
               << setw(8) << count << endl;
        }
      }
    } else {
      if ( logInfo ) cout << myname << "Channel counts are consistent with an offset of "
                          << tickdiff << " ticks." << endl;
    }
  } else {
    if ( logInfo ) cout << myname << "Channel clocks not checked for data retrieval from stroe." << endl;
  }

  // Read the raw digit status.
  const std::vector<raw::RDStatus>* pstats = nullptr;
  if ( useDecoderTool ) {
    pstats = pstatsFromTool.get();
  } else {
    art::InputTag itag2(m_DigitProducer, m_DigitName);
    auto hstats = evt.getHandle<std::vector<raw::RDStatus>>(itag2);
    if ( hstats ) {
      pstats = &*hstats;
    } else {
      if ( evt.isRealData() ) {
        cout << myname << "WARNING: Raw data status product not found." << endl;
      }
    }
  }

  // Check read status.
  string srdstat;
  bool skipEvent = skipAllEvents;
  if ( pstats != nullptr ) {
    if ( pstats->size() != 1 ) {
      cout << myname << "WARNING: Unexpected raw data status size: " << pstats->size() << endl;
    }
    const RDStatus rdstat = pstats->at(0);
    if ( false ) { 
      cout << myname << "Raw data status: " << rdstat.GetStatWord();
      if ( rdstat.GetCorruptDataDroppedFlag() ) cout << " (Corrupt data was dropped.)";
      if ( rdstat.GetCorruptDataKeptFlag() ) cout << " (Corrupt data was retained.)";
      cout << endl;
    }
    cout << myname << "Raw data read status: " << std::to_string(rdstat.GetStatWord()) << endl;
    srdstat = "rdstat=" + std::to_string(rdstat.GetStatWord());
    skipEvent |= skipEventsWithCorruptDataDropped && rdstat.GetCorruptDataDroppedFlag();
  }

  // Fetch beam information
  float beamTof = 0.0;    // Time of flight.
  if ( m_BeamEventLabel.size() ) {
    std::vector< art::Ptr<beam::ProtoDUNEBeamEvent> > beaminfo;
    auto pdbeamHandle = evt.getHandle< std::vector<beam::ProtoDUNEBeamEvent> >(m_BeamEventLabel);
    if ( pdbeamHandle ) {
      art::fill_ptr_vector(beaminfo, pdbeamHandle);
      if ( beaminfo.size() == 0 ) {
        cout << myname << "WARNING: Beam event vector is empty." << endl;
      } else {
        if ( beaminfo.size() > 1 ) {
          cout << myname << "WARNING: Beam event vector has size " << beaminfo.size() << endl;
        }
        AdcIndex beamTrigFlag = beaminfo[0]->GetTimingTrigger();
        if ( beamTrigFlag != trigFlag ) {
          if ( logInfo ) cout << myname << "Beam event and timing trigger flags differ: " << beamTrigFlag << " != " << trigFlag << endl;
        } else if ( beamTrigFlag != 12 ) {
          if ( logInfo ) cout << myname << "Beam event trigger is not beam: it is " << beamTrigFlag << endl;
        //} else if ( ! beaminfo[0]->CheckIsMatched() ) {
        //  cout << myname << "Beam event is not matched." << endl;
        } else if ( beaminfo[0]->GetTOFChan() == -1 ) {
          if ( logInfo ) cout << myname << "Beam event index does not indicate match." << endl;
        } else {
          int beamChan = beaminfo[0]->GetTOFChan();
          beamTof = beaminfo[0]->GetTOF();
          if ( logInfo ) cout << myname << "Beam event TOF[" << beamChan << "]: " << beamTof << endl;
        }
      }
    } else {
      if ( logInfo ) cout << myname << "Beam event data product not found: " << m_BeamEventLabel << endl;
    }
  }
            
  // Read in the digits. 
  const DigitVector* pdigits = nullptr;
  art::Handle<DigitVector> hdigits;
  if ( useDecoderTool ) {
    pdigits = pdigitsFromTool.get();
  } else {
    art::InputTag itag3(m_DigitProducer, m_DigitName);
    hdigits = evt.getHandle<DigitVector>(itag3);
    pdigits = &*hdigits;
  }
  if ( m_LogLevel >= 3 ) {
    cout << myname << "# digits read: " << pdigits->size() << endl;
  }
  if ( pdigits->size() == 0 ) mf::LogWarning("DataPrepModule") << "Input digit container is empty";

  // Create the container to hold the output wires.
  std::unique_ptr<std::vector<recob::Wire>> pwires(new std::vector<recob::Wire>);
  pwires->reserve(pdigits->size());

  // Create the association container.
  std::unique_ptr<art::Assns<raw::RawDigit,recob::Wire>> passns(new art::Assns<raw::RawDigit,recob::Wire>);

  // If status was bad, skip this event.
  // We store empty results to avoid exception.
  // We have to have read digits to store those results (yech).
  if ( skipEvent ) {
    if ( logInfo ) cout << myname << "Skipping event with " << srdstat << endl;
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
    pintStates = new WiredAdcChannelDataMap(m_IntermediateStates, pdigits->size());
  }

  // Create the transient data map and copy the digits there.
  AdcChannelDataMap fulldatamap;
  bool checkKeep = m_KeepChannelEnd > m_KeepChannelBegin;
  unsigned int nkeep = 0;
  unsigned int nskip = 0;
  unsigned int ndigi = pdigits->size();
  using EventInfo = AdcChannelData::EventInfo;
  EventInfo* pevt = new EventInfo;
  pevt->run = evt.run();
  pevt->subRun = evt.subRun();
  pevt->time = itim;
  pevt->timerem = itimrem;
  pevt->triggerClock = timingClock;
  pevt->trigger = trigFlag;
  AdcChannelData::EventInfoPtr pevtShared(pevt);
  for ( unsigned int idig=0; idig<ndigi; ++idig ) {
    const raw::RawDigit& dig = (*pdigits)[idig];
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
    // Fetch the channel status.
    Index chanStat = AdcChannelStatusGood;
    if ( m_pChannelStatusProvider != nullptr ) {
      if ( m_pChannelStatusProvider->IsNoisy(chan) ) chanStat = AdcChannelStatusNoisy;
      if ( m_pChannelStatusProvider->IsBad(chan)   ) chanStat = AdcChannelStatusBad;
    }
    // Fetch the online ID.
    AdcChannel fembID = AdcChannelData::badIndex();
    AdcChannel fembChannel = AdcChannelData::badIndex();
    if ( m_onlineChannelMapTool ) {
      unsigned int ichOn = m_onlineChannelMapTool->get(chan);
      if ( ichOn != IndexMapTool::badIndex() ) {
        fembID = ichOn/128;
        fembChannel = ichOn % 128;
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
    acd.setEventInfo(pevtShared);
    acd.setChannelInfo(chan, fembID, fembChannel, chanStat);
    acd.digitIndex = idig;
    acd.digit = &dig;
    if ( channelClocks.size() > idig ) acd.channelClock = channelClocks[idig];
    acd.metadata["ndigi"] = ndigi;
    if ( m_BeamEventLabel.size() ) {
      acd.metadata["beamTof"] = beamTof;
    }
    ++nkeep;
  }

  // Create a vector of data maps with an entry for each group.
  unsigned int nproc = 0;
  vector<AdcChannelDataMap> datamaps;
  if ( m_ChannelRanges.size() ) {
    const IndexRangeTool* pcrt = nullptr;
    string errmsg;
    DuneToolManager* ptm = DuneToolManager::instance();
    if ( ptm == nullptr ) {
      errmsg = "Tool manager not found.";
    } else {
      pcrt = ptm->getShared<IndexRangeTool>("channelRanges");
      if ( pcrt == nullptr ) errmsg = "Unable to find IndexRangeTool with name channelRanges.";
    }
    if ( pcrt == nullptr ) {
      cout << myname << "ERROR: IndexRangeTool not found: channelRanges" << endl;
    } else {
      for ( Name crn : m_ChannelRanges ) {
        IndexRange ran = pcrt->get(crn);
        if ( ran.isValid() ) {
          datamaps.emplace_back();
          AdcChannelDataMap& datamap = datamaps.back();
          for ( Index icha=ran.begin; icha<ran.end; ++icha ) {
            if ( fulldatamap.find(icha) != fulldatamap.end() ) {
              datamap.emplace(icha, move(fulldatamap[icha]));
            }
            ++nproc;
          }
        } else {
          cout << myname << "WARNING: Channel range not found: " << crn << endl;
        }
      }
    }
  } else if ( m_DoGroups ) {
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
  if ( logInfo ) {
    cout << myname << "              # input digits: " << ndigi << endl;
    cout << myname << "         # channels selected: " << nkeep << endl;
    cout << myname << "          # channels skipped: " << nskip << endl;
    cout << myname << "  # channels to be processed: " << nproc << endl;
  }

  DuneEventInfo devt;
  devt.run = evt.run();
  devt.subRun = evt.subRun();
  devt.event = evt.event();
  devt.time = itim;
  devt.timerem = itimrem;
  devt.triggerClock = timingClock;
  devt.triggerTick0 = timingClock/triggerPerTick;
  int bstat = m_pRawDigitPrepService->beginEvent(devt);
  if ( bstat ) cout << myname << "WARNING: Event initialization failed." << endl;

  auto const clockData = art::ServiceHandle<detinfo::DetectorClocksService const>()->DataFor(evt);
  for ( AdcChannelDataMap& datamap : datamaps ) {

    Index nacd = datamap.size();
    if ( nacd == 0 ) {
      if ( m_LogLevel >= 3 ) {
        if ( logInfo ) cout << myname << "Skipping empty data map." << endl;
      }
      continue;
    }

    // Use the data preparation service to build the wires and intermediate states.
    if ( m_LogLevel >= 3 ) {
      cout << myname << "Preparing " << nacd << " channel" << (nacd == 1 ? "" : "s") << "." << endl;
    }
    int rstat = m_pRawDigitPrepService->prepare(clockData, datamap, pwires.get(), pintStates);
    if ( rstat != 0 ) mf::LogWarning("DataPrepModule") << "Data preparation service returned error " << rstat;

    // Build associations between wires and digits.
    if ( m_DoAssns && !useDecoderTool ) {
      for ( const AdcChannelDataMap::value_type& iacd : datamap ) {
        const AdcChannelData& acd = iacd.second;
        AdcIndex idig = acd.digitIndex;
        if ( idig == AdcChannelData::badIndex() )
          throw art::Exception(art::errors::ProductRegistrationFailure) << "Digit index is not set.";
        AdcIndex iwir = acd.wireIndex;
        if ( iwir == AdcChannelData::badIndex() ) continue;
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

  int estat = m_pRawDigitPrepService->endEvent(devt);
  if ( estat ) cout << myname << "WARNING: Event finalization failed." << endl;

  if ( logInfo ) cout << myname << "Created wire count: " << pwires->size() << endl;
  if ( pwires->size() == 0 ) mf::LogWarning("DataPrepModule") << "No wires made for this event.";

  // Record wires and associations in the event.
  evt.put(std::move(pwires), m_WireName);
  if ( m_DoAssns ) {
    evt.put(std::move(passns), m_WireName);
  }

  // Record decoder containers.
  if ( useDecoderTool ) {
    if ( m_TimeStampName.size() ) evt.put(std::move(ptimsFromTool), m_TimeStampName);
    if ( m_OutputDigitName.size() ) evt.put(std::move(pdigitsFromTool), m_OutputDigitName);
    //evt.put(std::move(pstatsFromTool), "dataprep");
  }

  // Record intermediate state wires.
  for ( string sname : m_IntermediateStates ) {
    vector<Wire>* pintWires = pintStates->wires[sname];
    if ( pintWires == nullptr ) {
      cout << myname << "WARNING: Wires not found for state " << sname << "." << endl;
      continue;
    }
    if ( logInfo ) {
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
