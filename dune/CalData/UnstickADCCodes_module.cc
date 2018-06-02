////////////////////////////////////////////////////////////////////////
//
// UnstickADCCodes class
//
// jti3@fnal.gov
//
// Module to interpolate over stuck ADC codes (0x00 and 0x3F) 
// in 35t raw digits' 6 LSBs
//
////////////////////////////////////////////////////////////////////////

// C/C++ standard libraries
#include <string>
#include <vector>
#include <utility> // std::move()
#include <memory> // std::unique_ptr<>

// ROOT libraries
#include "TComplex.h"

// framework libraries
#include "cetlib_except/exception.h"
#include "cetlib/search_path.h"
#include "fhiclcpp/ParameterSet.h" 
#include "messagefacility/MessageLogger/MessageLogger.h" 
#include "art/Framework/Core/ModuleMacros.h" 
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Principal/Event.h" 
#include "art/Framework/Principal/Handle.h" 
#include "canvas/Persistency/Common/Ptr.h" 
#include "art/Framework/Services/Registry/ServiceHandle.h" 

// LArSoft libraries
#include "larcoreobj/SimpleTypesAndConstants/RawTypes.h" // raw::ChannelID_t
#include "larcore/Geometry/Geometry.h"
#include "lardataobj/RawData/RawDigit.h"
#include "lardataobj/RawData/raw.h"
#include "lardata/Utilities/AssociationUtil.h"

///creation of calibrated signals on wires
namespace unstick {

  class UnstickADCCodes : public art::EDProducer {

  public:
    
    // create calibrated signals on wires. this class runs 
    // an fft to remove the electronics shaping.     
    explicit UnstickADCCodes(fhicl::ParameterSet const& pset); 
    virtual ~UnstickADCCodes();
    
    void produce(art::Event& evt);
    void beginJob();
    void endJob();
    void reconfigure(fhicl::ParameterSet const& p);
    
  private:
    
    std::string  fDigitModuleLabel;  ///< module that made digits
                                                       ///< constants
    std::string  fSpillName;  ///< nominal spill is an empty string
                              ///< it is set by the DigitModuleLabel
                              ///< ex.:  "daq:preSpill" for prespill data

    raw::Compress_t        fCompression;      ///< compression type to use

    const unsigned int onemask = 0x003f; ///< Unsigned int ending in 111111 used to select 6 LSBs with bitwise AND
    unsigned int           fStickyADCCodesLimit; ///< Number of ADC codes to check for stickiness at 0x00 or 0x3f before stopping

  protected: 
    
  }; // class UnstickADCCodes

  DEFINE_ART_MODULE(UnstickADCCodes)
  
  //-------------------------------------------------
  UnstickADCCodes::UnstickADCCodes(fhicl::ParameterSet const& pset)
  {
    this->reconfigure(pset);
    produces< std::vector<raw::RawDigit>   >(fSpillName);

  }
  
  //-------------------------------------------------
  UnstickADCCodes::~UnstickADCCodes()
  {
  }

  //////////////////////////////////////////////////////
  void UnstickADCCodes::reconfigure(fhicl::ParameterSet const& p)
  {
 
    fDigitModuleLabel = p.get< std::string >("DigitModuleLabel", "daq");
    fStickyADCCodesLimit = p.get< unsigned int    >("StickyADCCodesLimit");
    fSpillName.clear();
    
    size_t pos = fDigitModuleLabel.find(":");
    if( pos!=std::string::npos ) {
      fSpillName = fDigitModuleLabel.substr( pos+1 );
      fDigitModuleLabel = fDigitModuleLabel.substr( 0, pos );
    }
    
  }

  //-------------------------------------------------
  void UnstickADCCodes::beginJob()
  {  
  }

  //////////////////////////////////////////////////////
  void UnstickADCCodes::endJob()
  {  
  }
  
  //////////////////////////////////////////////////////
  void UnstickADCCodes::produce(art::Event& evt)
  {      

    // make an unique_ptr of sim::SimDigits that allows ownership of the produced
    // digits to be transferred to the art::Event after the put statement below
    std::unique_ptr< std::vector<raw::RawDigit>   >  digcol(new std::vector<raw::RawDigit>);

    // Read in the digit List object(s). 
    art::Handle< std::vector<raw::RawDigit> > digitVecHandle;
    evt.getByLabel(fDigitModuleLabel, fSpillName, digitVecHandle);

    if (!digitVecHandle->size())  return;
    mf::LogInfo("UnstickADCCodes") << "UnstickADCCodes:: digitVecHandle size is " << digitVecHandle->size();

    // Use the handle to get a particular (0th) element of collection.
    art::Ptr<raw::RawDigit> digitVec0(digitVecHandle, 0);

    unsigned int dataSize = digitVec0->Samples(); //size of raw data vectors

  
    raw::ChannelID_t channel = raw::InvalidChannelID; // channel number

    std::vector<short> rawadc(dataSize);  // vector holding uncompressed adc values
    
    short *rawadc_a=0;

    // loop over all raw digits  
    digcol->reserve(digitVecHandle->size());
    for(size_t rdIter = 0; rdIter < digitVecHandle->size(); ++rdIter){ // ++ move
      
      // get the reference to the current raw::RawDigit
      art::Ptr<raw::RawDigit> digitVec(digitVecHandle, rdIter);
      channel = digitVec->Channel();
      fCompression = digitVec->Compression();
      float pedestal = digitVec->GetPedestal();
      // uncompress the data

      int pedestal_value = (int) digitVec->GetPedestal();
      raw::Uncompress(digitVec->ADCs(), rawadc, pedestal_value, fCompression);
      dataSize = rawadc.size();
      rawadc_a = rawadc.data();

      // loop over raw ADC vector and interpolate over stuck values
      for(size_t i = 0; i < dataSize; ++i){

	unsigned int sixlsbs = rawadc_a[i] & onemask;

	if(sixlsbs==onemask || sixlsbs==0){ //ADC code is stuck at 0x3f or 0x00
	  

	  if(i==0){ //if first ADC code is stuck, set it equal to pedestal value
	    rawadc_a[i] = (short) pedestal;
	    continue;
	  }

	  //find nearest preceding unstuck ADC code
	  //which should be immediately preceding ADC code since all earlier ADC codes have been interpolated to non-stuck values
	  size_t last_unstuck = i > 0 ? i - 1 : 0;
	  

	  //find nearest following unstuck ADC code
	  size_t next_unstuck = i;
	  unsigned short next_unstuck_sixlsbs;
	  unsigned short sixlsbs_stuck_in_a_row = 0;

	  do{

	    if(next_unstuck < dataSize - 1 && sixlsbs_stuck_in_a_row < fStickyADCCodesLimit)
	      ++next_unstuck;
	    else{ 
	      rawadc_a[next_unstuck] = (short) pedestal;
	      break;
	    }
	    next_unstuck_sixlsbs = rawadc_a[next_unstuck] & onemask;
	    if(next_unstuck_sixlsbs==0 || next_unstuck_sixlsbs==onemask)
	      ++sixlsbs_stuck_in_a_row;
	    else
	      sixlsbs_stuck_in_a_row = 0;

	  }
	  while(next_unstuck_sixlsbs==onemask || next_unstuck_sixlsbs==0);

	  // With last and next unstuck ADC codes found, interpolate linearly and replace stuck ADC code


	  float interpolated_unstuck_value = 1.0*(rawadc_a[next_unstuck] - rawadc_a[last_unstuck])/(next_unstuck-last_unstuck)*(i-last_unstuck)+rawadc_a[last_unstuck];
	  rawadc_a[i] = (short) interpolated_unstuck_value;

	}

	
      }

      raw::RawDigit rd(channel, dataSize, rawadc, raw::kNone );
      rd.SetPedestal( digitVec->GetPedestal(), digitVec->GetSigma() );
      digcol->push_back(rd);
 
    }

    evt.put(std::move(digcol), fSpillName);
 
    return;
  }
  


} // end namespace unstick


