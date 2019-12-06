
// Class:       SSPRawDecoder
// Module Type: producer
// File:        SSPRawDecoder_module.cc
//
// Generated at Thu Jul  6 18:31:48 2017 by Antonino Sergi,32 1-A14,+41227678738, using artmod
// from cetpkgsupport v1_11_00.
//
// Most recent additions by Bryan Ramson of FNAL (bjrams87@fnal.gov), Thursday September 20, 2018   
////////////////////////////////////////////////////////////////////////

// art includes
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art_root_io/TFileService.h"
#include "dune-raw-data/Services/ChannelMap/PdspChannelMapService.h"

// artdaq and dune-raw-data includes
#include "dune-raw-data/Overlays/SSPFragment.hh"
#include "artdaq-core/Data/ContainerFragment.hh"
// larsoft includes
#include "lardataobj/RawData/OpDetWaveform.h"
#include "lardataobj/RecoBase/OpHit.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "larcoreobj/SimpleTypesAndConstants/RawTypes.h"

// ROOT includes
#include "TH1.h"
#include "TH2.h"
#include "TGraph.h"

// C++ Includes
#include <memory>
#include <map>

namespace dune {
  class SSPRawDecoder;
}

class dune::SSPRawDecoder : public art::EDProducer {
public:
  explicit SSPRawDecoder(fhicl::ParameterSet const & p);
  // The destructor generated by the compiler is fine for classes
  // without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  SSPRawDecoder(SSPRawDecoder const &) = delete;
  SSPRawDecoder(SSPRawDecoder &&) = delete;
  SSPRawDecoder & operator = (SSPRawDecoder const &) = delete;
  SSPRawDecoder & operator = (SSPRawDecoder &&) = delete;

  // Required functions.
  void produce(art::Event & e) override;
  void reconfigure(const fhicl::ParameterSet &pset);
  void printParameterSet();
  
  struct trig_variables {
    unsigned int header;
    unsigned short length;
    unsigned short type;
    unsigned short status_flags;
    unsigned short header_type;
    unsigned short trig_id;
    unsigned short module_id;
    unsigned short channel_id;
    unsigned long timestamp_sync_delay;
    unsigned long timestamp_sync_count;
    uint64_t timestamp_nova; //make sure the data type actually contains the data...
    uint32_t peaksum;
    unsigned short peaktime;
    unsigned int prerise;
    unsigned int intsum;
    unsigned long baseline;
    unsigned long baselinesum;
    unsigned long cfd_interpol[4];
    unsigned long internal_interpol;
    uint64_t internal_timestamp; //internal timestamps necessary for 150 Mhz sample matching, if desired.
  };
  void readHeader(const SSPDAQ::EventHeader* daqHeader, struct trig_variables* tv);

  void getFragments(art::Event &evt,std::vector<artdaq::Fragment>* fragments);
  void beginJob() override;
  void endJob() override;
  void beginEvent(art::EventNumber_t eventNumber);
  void endEvent  (art::EventNumber_t eventNumber);

  void setRootObjects();

  recob::OpHit ConstructOpHit(trig_variables &trig, unsigned int channel);

private:

  // Declare member data here.
  std::string fRawDataLabel;
  bool        fSplitTriggers;
  std::string fOutputDataLabel;
  std::string fExtTrigOutputLabel;
  std::string fIntTrigOutputLabel;
  bool fUseChannelMap;
  bool fDebug;
  raw::Compress_t        fCompression;      ///< compression type to use
  unsigned int           fZeroThreshold;    ///< Zero suppression threshold

  uint32_t n_adc_counter_;  //counter of total number of ALL adc values in an event
  uint64_t adc_cumulative_; //cumulative total of ALL adc values in an event
  
  //debug and scaling variables scale all times to sample 0, per SSP
  uint64_t intreftime_[24];
  uint64_t int_ireftime_[24];
  uint64_t extreftime_[24];
  uint64_t ext_ireftime_[24];
  uint64_t allreftime;
  //int diff_time = 5; //~33 ns window for coincidences

  uint32_t         verb_adcs_;
  bool             verb_meta_;
  bool             timed_;

  TH1D * n_event_packets_; //diagnostic histos
  TH1D * frag_sizes_;
  TH1D * trig_ref_time_;
  //TH1D * trig_abs_time_;
  //TH2D * trig_adc_time_;

  // m1, m2, i1, i2
  double m1,m2,i1,i2;

  //double NOvAClockFrequency;
  double SPESize;

  //const uint64_t preread_ = 75000; //~.5 msecs before the external trigger (will be changed to 25 usecs soon)
  const uint32_t ext_trig_samp_time = 375000; //~2.5 msecs after the external trigger 
  //const uint32_t spillsamptime_ = 720000000; //~4.8 sec beam spill time
  
  //global vectors for smuggling relevant variables out of the producer method/function
  std::vector<unsigned short> coin_module_id, coin_channel_id;
  std::vector<double> coin_adc_peak;
  std::vector<uint64_t> coin_ext_time;
  std::vector<int32_t> coin_int_time;
                                                                                                                                        

  int number_of_packets = 12;  // 12 channels per SSP
  
  //mapping for SSPs to simple array
  std::map<int,int> ssp_map_ =
    { {11,0},
      {12,1},
      {13,2},
      {14,3},
      {21,4},
      {22,5},
      {23,6},
      {24,7},
      {31,8},
      {32,9},
      {33,10},
      {34,11},
      {41,12},
      {42,13},
      {43,14},
      {44,15},
      {51,16},
      {52,17},
      {53,18},
      {54,19},
      {61,20},
      {62,21},
      {63,22},
      {64,23} };

  std::map<size_t,TH1D*> trigger_type_; // internal vs. external  (16 internal, 48 external)

  //int smooth; // unused
  
  std::vector<raw::OpDetWaveform> waveforms; 
  std::vector<raw::OpDetWaveform> ext_waveforms;
  std::vector<raw::OpDetWaveform> int_waveforms;
  std::vector<recob::OpHit> hits;
  std::vector<recob::OpHit> ext_hits;
  std::vector<recob::OpHit> int_hits;

};


dune::SSPRawDecoder::SSPRawDecoder(fhicl::ParameterSet const & pset)
  : EDProducer{pset}
{
  reconfigure(pset);
  if (!fSplitTriggers) {
    produces< std::vector<raw::OpDetWaveform> > (fOutputDataLabel);
    produces< std::vector<recob::OpHit> > (fOutputDataLabel);
  }
  else{
    produces< std::vector<raw::OpDetWaveform> > (fExtTrigOutputLabel);
    produces< std::vector<raw::OpDetWaveform> > (fIntTrigOutputLabel);
    produces< std::vector<recob::OpHit> > (fExtTrigOutputLabel);
    produces< std::vector<recob::OpHit> > (fIntTrigOutputLabel);
  }
}

void dune::SSPRawDecoder::reconfigure(fhicl::ParameterSet const& pset) {

  fRawDataLabel = pset.get<std::string>("RawDataLabel");
  fSplitTriggers = pset.get<bool>("SplitTriggers");
  fOutputDataLabel = pset.get<std::string>("OutputDataLabel");
  fExtTrigOutputLabel = pset.get<std::string>("ExtTrigOutputLabel");
  fIntTrigOutputLabel = pset.get<std::string>("IntTrigOutputLabel");
  fUseChannelMap = pset.get<bool>("UseChannelMap");
  number_of_packets=pset.get<int>("number_of_packets");
  fDebug = pset.get<bool>("Debug");
  fZeroThreshold=0;
  fCompression=raw::kNone;

  if(fDebug) printParameterSet();

  verb_adcs_=pset.get<uint32_t>        ("verbose_adcs", 10000); 
  verb_meta_=pset.get<bool>            ("verbose_metadata", false); 
  n_adc_counter_=0; 
  adc_cumulative_=0; 
  
  // m1, i1, i2
  m1=pset.get<int>("SSP_m1"); 
  m2=pset.get<int>("SSP_m2");
  i1=pset.get<int>("SSP_i1"); 
  i2=pset.get<int>("SSP_i2"); 
  //NOvAClockFrequency=pset.get<double>("NOvAClockFrequency"); // in MHz
  SPESize=pset.get<double>("SPESize");
                                                       
  std::cout << "Parameters from the fcl file" << std::endl;
  std::cout << "m1: " << m1 << std::endl;
  std::cout << "m2: " << m2 << std::endl;
  std::cout << "i1: " << i1 << std::endl;
  std::cout << "i2: " << i2 << std::endl;
  //std::cout << "NOvAClockFrequency: " << NOvAClockFrequency << std::endl; 
  std::cout << "SPESize: " << SPESize << std::endl;
  std::cout << std::endl;

  //startTime=0;

}

void dune::SSPRawDecoder::printParameterSet(){

  for(int i=0;i<20;i++) std::cout << "=";
  std::cout << std::endl;
  std::cout << "Parameter Set" << std::endl;
  for(int i=0;i<20;i++) std::cout << "=";
  std::cout << std::endl;

  std::cout << "fRawDataLabel: " << fRawDataLabel << std::endl;
  if (!fSplitTriggers) {
    std::cout << "Not splitting triggers" << std::endl;
    std::cout << "fOutputDataLabel: " << fOutputDataLabel << std::endl;
  }
  else{
    std::cout << "Splitting triggers" << std::endl;
    std::cout << "fExtTrigOutputLabel: " << fExtTrigOutputLabel << std::endl;
    std::cout << "fIntTrigOutputLabel: " << fIntTrigOutputLabel << std::endl;
  }    
  std::cout << "fDebug: ";
  if(fDebug) std::cout << "true" << std::endl;
  else std::cout << "false" << std::endl;

  for(int i=0;i<20;i++) std::cout << "=";
  std::cout << std::endl;    
}

void dune::SSPRawDecoder::setRootObjects(){
  art::ServiceHandle<art::TFileService> tFileService;

  n_event_packets_ = tFileService->make<TH1D>("ssp_n_event_packets","SSP: n_event_packets",960,-0.5,959.5);  
  frag_sizes_ = tFileService->make<TH1D>("ssp_frag_sizes","SSP: frag_sizes",960,0,2e6);  
  trig_ref_time_ = tFileService->make<TH1D>("trig_ref_time_","trig_ref_time_",3750,0,ext_trig_samp_time);
  //trig_abs_time_ = tFileService->make<TH1D>("trig_abs_time_","trig_abs_time_",1000000,0,23000000000.0);
  //trig_adc_time_ = tFileService->make<TH2D>("trig_adc_time_","trig_adc_time_",10000,0.0,23000000000.0,1000,0.0,4000.0);
}

void dune::SSPRawDecoder::readHeader(const SSPDAQ::EventHeader* daqHeader, struct trig_variables* tv){

  
  tv->header = daqHeader->header;                          // the 'start of header word' (should always be 0xAAAAAAAA)
  tv->length = daqHeader->length;                          // the length of the packet in unsigned ints (including header)
  tv->type = ((daqHeader->group1 & 0xFF00) >> 8);          // packet type
  tv->status_flags = ((daqHeader->group1 & 0x00F0) >> 4);  // status flags
  tv->header_type = ((daqHeader->group1 & 0x000F) >> 0);   // header type
  tv->trig_id = daqHeader->triggerID;                      // the packet ID
  tv->module_id = ((daqHeader->group2 & 0xFFF0) >> 4);     // module ID
  tv->channel_id = ((daqHeader->group2 & 0x000F) >> 0);    // channel ID

                                                           // external timestamp sync delay (FP mode)
  tv->timestamp_sync_delay = ((unsigned int)(daqHeader->timestamp[1]) << 16) + (unsigned int)(daqHeader->timestamp[0]);
  // external timestamp sync count (FP mode)
  tv->timestamp_sync_count = ((unsigned int)(daqHeader->timestamp[3]) << 16) + (unsigned int)(daqHeader->timestamp[2]);
  // get the external timestamp (NOvA mode)
  tv->timestamp_nova = ((unsigned long)daqHeader->timestamp[3] << 48) + ((unsigned long)daqHeader->timestamp[2] << 32)
    + ((unsigned long)daqHeader->timestamp[1] << 16) + ((unsigned long)daqHeader->timestamp[0] << 0);
 

  tv->peaksum = ((daqHeader->group3 & 0x00FF) >> 16) + daqHeader->peakSumLow;  // peak sum
  if(tv->peaksum & 0x00800000) {
    tv->peaksum |= 0xFF000000;
  }

  tv->peaktime = ((daqHeader->group3 & 0xFF00) >> 8);                                  // peak time
  tv->prerise = ((daqHeader->group4 & 0x00FF) << 16) + daqHeader->preriseLow;          // prerise
  tv->intsum = ((unsigned int)(daqHeader->intSumHigh) << 8) + (((unsigned int)(daqHeader->group4) & 0xFF00) >> 8);  // integrated sum
  tv->baseline = daqHeader->baseline;                                                  // baseline
  tv->baselinesum = ((daqHeader->group4 & 0x00FF) << 16) + daqHeader->preriseLow;      // baselinesum
  for(unsigned int i_cfdi = 0; i_cfdi < 4; i_cfdi++)                                   // CFD timestamp interpolation points
    tv->cfd_interpol[i_cfdi] = daqHeader->cfdPoint[i_cfdi];

  tv->internal_interpol = daqHeader->intTimestamp[0];                                  // internal interpolation point
                                                                                       // internal timestamp
  tv->internal_timestamp = ((uint64_t)((uint64_t)daqHeader->intTimestamp[3] << 32)) + ((uint64_t)((uint64_t)daqHeader->intTimestamp[2]) << 16) + ((uint64_t)((uint64_t)daqHeader->intTimestamp[1]));
  
}

void dune::SSPRawDecoder::getFragments(art::Event &evt, std::vector<artdaq::Fragment> *fragments){

  art::EventNumber_t eventNumber = evt.event();

  art::Handle<artdaq::Fragments> rawFragments;
  art::Handle<artdaq::Fragments> containerFragments;

  bool have_data = true;

  /// look for Container Fragments:
  evt.getByLabel(fRawDataLabel, "ContainerPHOTON", containerFragments);
  // Check if there is SSP data in this event
  // Don't crash code if not present, just don't save anything    
  try { containerFragments->size(); }
  catch(std::exception e)  {
    //std::cout << "WARNING: Container SSP data not found in event " << eventNumber << std::endl;
    have_data = false;
  }

  if (have_data)
    {
      //Check that the data are valid
      if(!containerFragments.isValid()){
        MF_LOG_ERROR("SSPRawDecoder") << "Run: " << evt.run()
                                   << ", SubRun: " << evt.subRun()
                                   << ", Event: " << eventNumber
                                   << " Container Fragments found but NOT VALID";
        return;
      }

      for (auto cont : *containerFragments)
        {
          //std::cout << "container fragment type: " << (unsigned)cont.type() << std::endl;
          artdaq::ContainerFragment contf(cont);
          for (size_t ii = 0; ii < contf.block_count(); ++ii)
            {
              size_t fragSize = contf.fragSize(ii);
              frag_sizes_->Fill(fragSize);
              //artdaq::Fragment thisfrag;
              //thisfrag.resizeBytes(fragSize);
            
              //memcpy(thisfrag.headerAddress(), contf.at(ii), fragSize);
              fragments->emplace_back(*contf[ii]);
            }
        }
    }

  /// Look for non-container Raw Fragments:

  bool have_data2=true;

  evt.getByLabel(fRawDataLabel, "PHOTON", rawFragments);
    
  // Check if there is SSP data in this event
  // Don't crash code if not present, just don't save anything
  try { rawFragments->size(); }
  catch(std::exception e) {
    //std::cout << "WARNING: Raw SSP data not found in event " << eventNumber << std::endl;
    have_data2=false;
  }

  if (have_data2)
    {
      //Check that the data is valid
      if(!rawFragments.isValid()){

        MF_LOG_ERROR("SSPRawDecoder") << "Run: " << evt.run()
                                   << ", SubRun: " << evt.subRun()
                                   << ", Event: " << eventNumber
                                   << " Non-Container Fragments found but NOT VALID";
        return;
      }
      for(auto const& rawfrag: *rawFragments){
        fragments->emplace_back( rawfrag );
      }
    }
  
}

void dune::SSPRawDecoder::beginJob(){
  //intializing normalizing time references for use in the producer method
  setRootObjects();
  allreftime=0;
  for(int i=0;i<24;i++){ 
    int_ireftime_[i]=0;
    ext_ireftime_[i]=0;
  }
}

void dune::SSPRawDecoder::beginEvent(art::EventNumber_t /*eventNumber*/)
{
  //intializing adc counters and internal references
  n_adc_counter_  = 0;
  adc_cumulative_ = 0;
  for(int i=0;i<24;i++) {
  intreftime_[i]=0;
  // extreftime_[i]=0;
  }
  timed_ = false;
  
}

void dune::SSPRawDecoder::endEvent(art::EventNumber_t eventNumber)
{
  
  //These should only work with "internal" trig.trigger_type=16 triggers.
  if(coin_ext_time.size() > 0){
    for(size_t i=0;i<coin_ext_time.size();i++){
      trig_ref_time_->Fill(coin_int_time[i]);                                
      //trig_abs_time_->Fill(coin_ext_time[i]-allreftime);
      //trig_adc_time_->Fill(coin_ext_time[i]-allreftime,coin_adc_peak[i]);
        
    }
  }
  coin_module_id.clear();
  coin_channel_id.clear();
  coin_int_time.clear();
  coin_ext_time.clear();
  coin_adc_peak.clear();
}
void dune::SSPRawDecoder::endJob(){

}

void dune::SSPRawDecoder::produce(art::Event & evt){

  art::ServiceHandle<art::TFileService> tFileService;
  art::ServiceHandle<dune::PdspChannelMapService> channelMap;

  //MF_LOG_INFO("SSPRawDecoder") << "-------------------- SSP RawDecoder -------------------";
  // Implementation of required member function here.

  art::EventNumber_t eventNumber = evt.event();  
  
  /// Get the fragments (Container or Raw)
  std::vector<artdaq::Fragment> fragments;
  getFragments(evt,&fragments);

  unsigned int allPacketsProcessed = 0;
  unsigned int waveform_counter = 0;
  uint64_t ssptrigtime = 0;
  std::map<int, int> packets_per_fragment;

  // just to make sure -- the std::move from the previous event should clear them out, but this is
  // not guaranteed by the standard.

  waveforms.clear();
  int_waveforms.clear();
  ext_waveforms.clear();
  hits.clear();
  int_hits.clear();
  ext_hits.clear();
  
  /// Process all packets:
  
  for(auto const& frag: fragments){
    if((unsigned)frag.type() != 3) continue;
 
    ///> Create a SSPFragment from the generic artdaq fragment
    dune::SSPFragment sspf(frag);
    
    const SSPDAQ::MillisliceHeader* meta=0;
    
     ///> get the information from the header
    if(frag.hasMetadata()) meta = &(frag.metadata<SSPFragment::Metadata>()->sliceHeader); ///> get the metadata
    else std::cout << "SSP fragment has no metadata associated with it." << std::endl;
    
    ///> get a pointer to the first packet in the millislice
    const unsigned int* dataPointer = sspf.dataBegin();
    
    ///> loop over the packets in the millislice
    unsigned int packetsProcessed=0;
    while(( meta==0 || packetsProcessed<meta->nTriggers) && dataPointer<sspf.dataEnd() ){
      
      ///> get the packet header
      const SSPDAQ::EventHeader* daqHeader=reinterpret_cast<const SSPDAQ::EventHeader*>(dataPointer);
      
      /// read the header to provide the trigger variables structure        
      struct trig_variables trig;
      readHeader(daqHeader, &trig);
      
      auto const* ts = lar::providerFrom<detinfo::DetectorClocksService>();
      
      /// time
      ////long double time = trig.timestamp_nova*ts->OpticalClock().TickPeriod(); //in experiment microseconds
      // DO NOT USE ts->OpticalClock().TickPeriod()!!!! It is not precise enough
      // use OpticalClock().Frequency, and do the division yourself with high precission.
      double time = double(trig.timestamp_nova % 1000000000 ) / double(ts->OpticalClock().Frequency());
      //true time truncated by 10 digits in order to make sure the math works correctly
      //std::cout << time << std::endl;
      unsigned int channel = ((trunc(frag.fragmentID()/10) -1 )*4 + frag.fragmentID()%10 -1 )*number_of_packets + trig.channel_id;

      //set external reference time to the first time stamp of the run, if a lower time stamp is found, adjust
      if(allreftime==0 || (allreftime > trig.timestamp_nova)) allreftime=trig.timestamp_nova;
      
      //internal and external reference times on external (beam/cosmic window triggers). Can be used to time external timestamp down to the internal timesample. Might want to try this at some point in production...
      if(trig.type==48) {
        if (ssptrigtime==0) {
          ssptrigtime=trig.timestamp_nova;
          if(verb_meta_) std::cout << "SSP: " << ssptrigtime << std::endl; 
        }
        
        intreftime_[ssp_map_[trig.module_id]] = trig.internal_timestamp;  
        extreftime_[ssp_map_[trig.module_id]] = trig.timestamp_nova; 
        if(int_ireftime_[ssp_map_[trig.module_id]] == 0) int_ireftime_[ssp_map_[trig.module_id]] = trig.internal_timestamp;
        if(ext_ireftime_[ssp_map_[trig.module_id]] == 0) ext_ireftime_[ssp_map_[trig.module_id]] = trig.timestamp_nova;  
      }
     
      // Trigger type histogram
      if (trigger_type_.find(channel) == trigger_type_.end())
        {
          TH1D* tth = tFileService->make<TH1D>(Form("trigger_type_channel_%d",channel),Form("trigger_type_channel_%d",channel),4,0,3);
          tth->SetTitle(Form("Trigger type - Channel %d",channel));
          tth->GetXaxis()->SetTitle("Trigger type");
          tth->GetXaxis()->SetBinLabel(2,"Internal (16)");
          tth->GetXaxis()->SetBinLabel(3,"External (48)");
          trigger_type_[channel] = tth;
        }
      
      if ( trig.type == 16 ) trigger_type_[channel]->Fill(1);
      if ( trig.type == 48 ) trigger_type_[channel]->Fill(2);
      
      ///> increment the data pointer past the packet header
      dataPointer+=sizeof(SSPDAQ::EventHeader)/sizeof(unsigned int);
      
      //>get the information from the data
      bool verb_values = false;
      
      ///> get the number of ADC values in the packet
      unsigned int nADC=(trig.length-sizeof(SSPDAQ::EventHeader)/sizeof(unsigned int))*2;
      
      ///> get a pointer to the first ADC value
      const unsigned short* adcPointer=reinterpret_cast<const unsigned short*>(dataPointer);
      
      char histname[100];
      sprintf(histname,"evt%i_frag%d_wav%d",eventNumber, frag.fragmentID(), packetsProcessed);
      
      TH1D* hist=new TH1D("hist","hist",nADC,0,nADC);
      
      // map the channel number to offline if requested
      
      unsigned int mappedchannel = channel;
      if (fUseChannelMap) mappedchannel = channelMap->SSPOfflineChannelFromOnlineChannel(channel);
      
      // Get information from the header, //added by Jingbo
      
      unsigned short     OpChannel   =  (unsigned short) mappedchannel;   ///< Derived Optical channel
      raw::OpDetWaveform Waveform(time, OpChannel, nADC);
      
      //calculating relevant values in decoder because what comes out of the trigger header seems incorrect-Bryan Ramson
      unsigned long calbasesum = 0;
      unsigned short maxadc = 0;
      unsigned long calintsum = 0;
      unsigned short  calpeaktime = 0;
      unsigned int calpeaksum =0;
      ///> copy the waveforms
      for(size_t idata = 0; idata < nADC; idata++) {
        ///> get the 'idata-th' ADC value
        const unsigned short* adc = adcPointer + idata;
        if(idata < i1) calbasesum +=  static_cast<unsigned long>(*adc); //added by Bryan Ramson
        if(idata > i1+m1 && idata <= i2+i1+m1) calintsum += static_cast<unsigned long>(*adc); //added by Bryan Ramson
        if(idata >= i1+m1+m2 && idata <= i1+2*m1+m2) calpeaksum += static_cast<unsigned int>(*adc); //added by Bryan Ramson
        
        maxadc = std::max(maxadc,*adc); //added by Bryan Ramson
        if(maxadc == *adc) calpeaktime = idata; //added by Bryan Ramson
        Waveform.push_back(*adc); //added by Jingbo
        waveform_counter++;
        
        n_adc_counter_++;
        adc_cumulative_ += (uint64_t)(*adc);
        
        ///> Waveform 
        hist->SetBinContent(idata+1,*adc);
        
        if (idata >= verb_adcs_) verb_values = false;
        
        verb_values = false; //don't print adc. Added by J.Wang
        if(verb_values) {
          if(idata == 0&&verb_adcs_>0) std::cout << "Printing the " << nADC << " ADC values saved with the packet:" << std::endl;
          std::cout << *adc << " ";
        }
      }// idata

      // pedestal, area and peak (according to the Register table, the  SSP User Manual has i1 and i2 inverted)
      double pedestal = calbasesum / ((double)i1);    
      double area = calintsum-(pedestal*i2);
      if(area<0) area=0; //On external triggers area over "peak" less pedestal could be negative which is nonsense.
      double peak = maxadc;
     
      trig.baselinesum = calbasesum;
      trig.intsum = calintsum;
      trig.peaktime = calpeaktime;
      trig.peaksum = calpeaksum;
      
      if(verb_meta_) {
        std::cout
          << "Channel:                            " << channel                   << std::endl
          << "Header:                             " << trig.header               << std::endl
          << "Length:                             " << trig.length               << std::endl
          << "Trigger type:                       " << trig.type                 << std::endl
          << "Status flags:                       " << trig.status_flags         << std::endl
          << "Header type:                        " << trig.header_type          << std::endl
          << "Trigger ID:                         " << trig.trig_id              << std::endl
          << "Module ID:                          " << trig.module_id            << std::endl
          << "Channel ID:                         " << trig.channel_id           << std::endl
          << "External timestamp (F mode):        "                              << std::endl
          << "  Sync delay:                       " << trig.timestamp_sync_delay << std::endl
          << "  Sync count:                       " << trig.timestamp_sync_count << std::endl
          << "External timestamp (NOvA mode):     " << trig.timestamp_nova       << std::endl
          << "Peak sum:                           " << trig.peaksum              << std::endl
          << "Peak time:                          " << trig.peaktime             << std::endl
          << "Prerise:                            " << trig.prerise              << std::endl
          << "Integrated sum:                     " << trig.intsum               << std::endl
          << "Baseline sum:                       " << trig.baselinesum          << std::endl
          << "CFD Timestamp interpolation points: " << trig.cfd_interpol[0]      << " " << trig.cfd_interpol[1] << " " << trig.cfd_interpol[2]   << " " << trig.cfd_interpol[3] << std::endl
          << "Internal interpolation point:       " << trig.internal_interpol    << std::endl
          << "Internal timestamp:                 " << trig.internal_timestamp   << std::endl
          << std::endl
          << "Pedestal                            " << pedestal                  << std::endl
          << "Area                                " << area                      << std::endl
          << "Peak heigth                         " << peak                      << std::endl
          << std::endl;
      }
            
      //Fill diagnostic vectors for use at the end of the event
      coin_module_id.push_back(trig.module_id);
      coin_channel_id.push_back(trig.channel_id);
      coin_int_time.push_back(trig.internal_timestamp-intreftime_[ssp_map_[trig.module_id]]);
      coin_ext_time.push_back(trig.timestamp_nova);
      coin_adc_peak.push_back(peak);
     
      ///> increment the data pointer to the end of the current packet (to the start of the next packet header, if available)
      dataPointer+=nADC/2;
      
      // Put waveform and ophit into collections
      // Split into internal and external triggers if that has been set.
      if (!fSplitTriggers) {
        waveforms.emplace_back( Waveform );
        hits.emplace_back( ConstructOpHit(trig, mappedchannel) );
      }
      else{
        if (trig.type == 48 ) {
          ext_waveforms.emplace_back( Waveform );
          ext_hits.emplace_back( ConstructOpHit(trig, mappedchannel) );
        }
        else if (trig.type == 16) {
          int_waveforms.emplace_back( Waveform );
          int_hits.emplace_back( ConstructOpHit(trig, mappedchannel) );
        }
        else {
          std::cerr << "Unknown trigger type " << trig.type << ", cannot assign to appropriate data product with SplitTriggers enabled." << std::endl;
        }
      }

      hist->Delete();
    
      ++packetsProcessed; // packets
    }
    
    packets_per_fragment[frag.fragmentID()] = packetsProcessed;
    allPacketsProcessed += packetsProcessed;
  }//frag: fragments
  
  n_event_packets_->Fill(allPacketsProcessed);
  
  endEvent(eventNumber);


  if (!fSplitTriggers) {
    evt.put(std::make_unique<decltype(waveforms)>(std::move(waveforms)), fOutputDataLabel);
    evt.put(std::make_unique<decltype(hits)>(     std::move(hits)),      fOutputDataLabel);
  }
  else {
    evt.put(std::make_unique<decltype(ext_waveforms)>(std::move(ext_waveforms)), fExtTrigOutputLabel);
    evt.put(std::make_unique<decltype(ext_hits)>(     std::move(ext_hits)),      fExtTrigOutputLabel);
    evt.put(std::make_unique<decltype(int_waveforms)>(std::move(int_waveforms)), fIntTrigOutputLabel);
    evt.put(std::make_unique<decltype(int_hits)>(     std::move(int_hits)),      fIntTrigOutputLabel);
  }
}



recob::OpHit dune::SSPRawDecoder::ConstructOpHit(trig_variables &trig, unsigned int channel)
{
  // Get basic information from the header
  unsigned short     OpChannel   = channel;         ///< Derived Optical channel
  unsigned long      FirstSample = trig.timestamp_nova;
  double             TimeStamp   = ((double)FirstSample); ///< first sample experiment time in microseconds

  auto const* ts = lar::providerFrom<detinfo::DetectorClocksService>();
  double peakTime = ((double) trig.peaktime) * ts->OpticalClock().TickPeriod(); // microseconds
  double width = ((double)i1) * ts->OpticalClock().TickPeriod(); // microseconds
  double pedestal = ( (double) trig.baselinesum ) / ( (double) i1 );
  double area =     ( (double) trig.intsum      ) - pedestal * ( (double) i2 );
  double peak =     ( (double) trig.peaksum     ) / ( (double) m1 ) - pedestal;
  
  
  recob::OpHit ophit(OpChannel,
                     TimeStamp+peakTime,  // Relative Time
                     TimeStamp+peakTime,  // Absolute time
                     0,          // Frame, not used by DUNE
                     width,
                     area,
                     peak,
                     area / SPESize, // PE
                     0.);
  
  return ophit;
}


DEFINE_ART_MODULE(dune::SSPRawDecoder)
