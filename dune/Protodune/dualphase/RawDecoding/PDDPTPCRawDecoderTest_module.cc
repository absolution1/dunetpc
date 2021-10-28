////////////////////////////////////////////////////////////////////////
// Class:       PDDPTPCRawDecoderTest
// Plugin Type: analyzer (art v3_02_05)
// File:        PDDPTPCRawDecoderTest_module.cc
//
// Generated at Tue Jun 11 16:37:45 2019 by Vyacheslav Galymov using cetskelgen
// from cetlib version v3_07_02.
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "art_root_io/TFileService.h" 
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "lardataobj/RawData/raw.h"
#include "lardataobj/RawData/RawDigit.h"

//#include "TH1.h"
#include "TH2.h"

#include <string>

class PDDPTPCRawDecoderTest;


class PDDPTPCRawDecoderTest : public art::EDAnalyzer {
public:
  explicit PDDPTPCRawDecoderTest(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  PDDPTPCRawDecoderTest(PDDPTPCRawDecoderTest const&) = delete;
  PDDPTPCRawDecoderTest(PDDPTPCRawDecoderTest&&) = delete;
  PDDPTPCRawDecoderTest& operator=(PDDPTPCRawDecoderTest const&) = delete;
  PDDPTPCRawDecoderTest& operator=(PDDPTPCRawDecoderTest&&) = delete;

  // Required functions.
  void analyze(art::Event const& e) override;

  // Selected optional functions.
  void beginJob() override;
  void endJob() override;

private:
  // Declare member data here.
  std::string __RawDigitLabel;
  unsigned    __Ticks;
  unsigned    __Chans;
    
  TH2I* __AdcData;
};


PDDPTPCRawDecoderTest::PDDPTPCRawDecoderTest(fhicl::ParameterSet const& p)
  : EDAnalyzer{p}  // ,
  // More initializers here.
{
  // Call appropriate consumes<>() for any products to be retrieved by this module.
  __RawDigitLabel  =  p.get< std::string >("RawDigitLabel");
  __Ticks          =  p.get< unsigned >("Ticks");
  __Chans          =  p.get< unsigned >("Chans");
  
}

void PDDPTPCRawDecoderTest::analyze(art::Event const& e)
{
  auto Raw = e.getHandle< std::vector<raw::RawDigit> >(__RawDigitLabel);
  std::vector< art::Ptr<raw::RawDigit> >  Digits;
  art::fill_ptr_vector(Digits, Raw);

  //loop through all RawDigits (over entire channels)
  mf::LogInfo("") << "Total number of channels "<<Digits.size();
  
  for(auto &digit : Digits)
    {
      auto chan    = digit->Channel();
      auto samples = digit->Samples();
      //mf::LogInfo("") << "vector size " << chan <<" "<< samples;

      raw::RawDigit::ADCvector_t data(samples);
      // note: only works with uncompressed data for now
      raw::Uncompress(digit->ADCs(), data, digit->Compression());
      
      if( chan >= __Chans ) break; //continue;
      unsigned tick = 0;
      for( auto &adc : data )
	{
	  if( tick >= __Ticks ) break;
	  __AdcData->Fill( chan, tick++, adc );
	}
    }
}

void PDDPTPCRawDecoderTest::beginJob()
{
  // Implementation of optional member function here.
  art::ServiceHandle<art::TFileService> tfs;
  
  // put all channels in single histogram
  __AdcData = tfs->make<TH2I>( "amc_data","AMC Data",
			       __Chans, 0, __Chans, __Ticks, 0, __Ticks );
}

void PDDPTPCRawDecoderTest::endJob()
{
  // Implementation of optional member function here.
}

DEFINE_ART_MODULE(PDDPTPCRawDecoderTest)
