//////////////////////////////////////////////////////////////////////////
// Class:       ImportSingle311Event
// Module Type: producer
// File:        test_code.cc
//
// Generated at Fri Mar 17 13:56:31 2017 by Kevin Fusshoeller
//////////////////////////////////////////////////////////////////////////

#include <vector>
#include <string>
#include <sstream>
#include <iostream>

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"

#include "fhiclcpp/ParameterSet.h"
#include "dune/3x1x1dp/DataImport/Services/EventDecoder.h"

#include "lardataobj/RawData/raw.h"
#include "lardataobj/RawData/RawDigit.h"
#include "larcore/Geometry/Geometry.h"

using std::ostringstream;
using std::string;

using namespace dlardaq;

namespace EventGen{

  class ImportSingle311Event : public art::EDProducer {
  public:
    explicit ImportSingle311Event(fhicl::ParameterSet const & p);
    
    void produce(art::Event &evt);
    void beginJob();
    void endJob();

  private: 

    void SplitAdc(const std::vector<adc16_t> adc, size_t channel, uint32_t num_samples, std::vector<short> &adclist, size_t LArchan);
    size_t Get311Chan(size_t LAr_chan);

    double GetPedMean(size_t RawChan, std::vector< std::pair<double, double> > *fPedMap){ return fPedMap->at(RawChan).first; }
    double GetPedRMS(size_t RawChan, std::vector< std::pair<double, double> > *fPedMap){ return fPedMap->at(RawChan).second; }

    EventDecoder DataDecode{1280, 1};

    // fcl parameters
    std::string fFilename;
    size_t fEvt_num;
    std::string fPedestalFile;

    std::vector< std::pair<double, double> > fPedMap; 
  }; // Class ImportSingle311Event    


  //////////////////////////////////////////////////////////////////////////
  // function to read in pedestals from file
  /////////////////////////////////////////////////////////////////////////
  void ReadPedestalFile(std::string PedestalFileName, std::vector < std::pair<double, double> > &PedMap)
  {
    //initialize the channel-ped value map
    std::ifstream file;
    file.open(PedestalFileName);
    if( !file.is_open() )
    {
      throw art::Exception( art::errors::FileReadError )
      		<< "failed to open input file " << PedestalFileName << "\n";
    }
    
    while(!file.eof())
    {
      size_t ch, cryo, crate, rawch;
      double mean, rms;
      file >> rawch >> cryo >> crate >> ch >> mean >> rms;
      PedMap.emplace_back(mean, rms);
    }
    
    file.close();
    return;
  } // ReadPedestalFile()

  
  /////////////////////////////////////////////////////////////////////////////////////////////////////
  //
  /////////////////////////////////////////////////////////////////////////////////////////////////////

  DEFINE_ART_MODULE(EventGen::ImportSingle311Event)  

  ImportSingle311Event::ImportSingle311Event(fhicl::ParameterSet const & p)
  {
    fFilename = p.get<std::string>("Filename");
    fEvt_num = p.get<size_t>("Evt_num");
    fPedestalFile = p.get<std::string>("PedestalFile");
    produces< std::vector<raw::RawDigit> >();
  }


  void ImportSingle311Event::beginJob(){ 
    DataDecode.Open(fFilename);
  }


  void ImportSingle311Event::produce(art::Event &evt){
      ReadPedestalFile(fPedestalFile, fPedMap);

      //ImportSingle311Event(fhicl::ParameterSet const & p);
      //beginJob();

      uint32_t nsample;
      size_t nch = DataDecode.GetNCh();
      dlardaq::evheader_t EveHead;
      std::vector<adc16_t> adcvec;
      
      ssize_t blabla = DataDecode.GetEvent(fEvt_num, EveHead, adcvec);
      std::cout << "Fetching event " << blabla << "\n";
      nsample = (uint32_t) EveHead.ev_size/(nch*1.5);
      std::cout << "The size of the event is " << nsample << "\n";
      auto max_adcvec = std::max_element(std::begin(adcvec), std::end(adcvec));
      std::cout << "Maximum element of adcvec is: " << *max_adcvec << "\n";
      std::cout << "The size of adcvec is: " << adcvec.size() << "\n";
      // std::cout << "The size in bytes of adcvec is: " << sizeof(adcvec) << "\n";

      std::vector<short> adclist;
      std::unique_ptr<std::vector<raw::RawDigit>> digcol(new std::vector<raw::RawDigit>);
      for(size_t LAr_chan = 0; LAr_chan < nch; LAr_chan++){
        adclist.clear();
	size_t Chan311 = Get311Chan(LAr_chan);
        SplitAdc(adcvec, Chan311, nsample, adclist, LAr_chan);

        short unsigned int nTickReadout = nsample;
        raw::ChannelID_t channel = LAr_chan;
        raw::Compress_t comp = raw::kNone;
	
        raw::RawDigit rd(channel, nTickReadout, adclist, comp);
        double pedval = ImportSingle311Event::GetPedMean(Chan311, &fPedMap);
        double pedrms = ImportSingle311Event::GetPedRMS(Chan311, &fPedMap);
        rd.SetPedestal(pedval, pedrms);

        digcol->push_back(rd);
	}

      evt.put(std::move(digcol));
      //endJob();
  }  


  void ImportSingle311Event::SplitAdc(const std::vector<adc16_t> adc, size_t channel, uint32_t num_samples, std::vector<short> &adclist, size_t LArchan){
    adclist.resize(num_samples);
    for (uint32_t i = 0; i<num_samples; i++){
      adclist[i] = static_cast<short>(adc[channel*num_samples + i]);
    }
  }
  
  //-------------------------------
  size_t ImportSingle311Event::Get311Chan(size_t LAr_chan){

    size_t crate = LAr_chan / 320;
    size_t Chan311;

    LAr_chan = 8*(LAr_chan/8+1)-LAr_chan%8 -1;

    if(crate == 0)
      {
	LAr_chan = 32*(LAr_chan/32+1)-LAr_chan%32 -1;
        size_t card = 4 - ((LAr_chan / 32) % 5);
        if(LAr_chan > 159)
          {
            size_t shift = 31 - (LAr_chan % 32);
            Chan311 = (2*card)*32 + shift;
          }
        else
          {
            size_t shift = 31 - (LAr_chan % 32);
            Chan311 = (2*card + 1)*32 + shift;
          }
      }
    else
      {
        size_t new_LAr_chan = LAr_chan - crate*320;
        size_t card = ((new_LAr_chan / 32) % 5);
        if(new_LAr_chan > 159)
          {
            size_t shift = new_LAr_chan % 32;
            Chan311 = (2*card)*32 + shift;
          }
        else
          {
            size_t shift = new_LAr_chan % 32;
            Chan311 = (2*card + 1)*32 + shift;
          }
        Chan311 = Chan311 + crate*320;
      } // end of if/else statementi

    return Chan311;
  } // Get311Chan

  void ImportSingle311Event::endJob(){
    DataDecode.Close();
  } 
} // EventGen namespace


