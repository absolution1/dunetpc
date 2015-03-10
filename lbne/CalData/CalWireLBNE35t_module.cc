////////////////////////////////////////////////////////////////////////
//
// CalWireLBNE35t class
//
// brebel@fnal.gov
//
// 11-3-09 Pulled all FFT code out and put into Utilitiess/LArFFT
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
#include "cetlib/exception.h"
#include "cetlib/search_path.h"
#include "fhiclcpp/ParameterSet.h" 
#include "messagefacility/MessageLogger/MessageLogger.h" 
#include "art/Framework/Core/ModuleMacros.h" 
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Principal/Event.h" 
#include "art/Framework/Principal/Handle.h" 
#include "art/Persistency/Common/Ptr.h" 
#include "art/Framework/Services/Registry/ServiceHandle.h" 

// LArSoft libraries
#include "SimpleTypesAndConstants/RawTypes.h" // raw::ChannelID_t
#include "Geometry/Geometry.h"
#include "Filters/ChannelFilter.h"
#include "RawData/RawDigit.h"
#include "RawData/raw.h"
#include "RecoBase/Wire.h"
#include "RecoBaseArt/WireCreator.h"
#include "Utilities/LArFFT.h"
#include "Utilities/AssociationUtil.h"
#include "lbne/Utilities/SignalShapingServiceLBNE35t.h"


///creation of calibrated signals on wires
namespace caldata {

  class CalWireLBNE35t : public art::EDProducer {

  public:
    
    // create calibrated signals on wires. this class runs 
    // an fft to remove the electronics shaping.     
    explicit CalWireLBNE35t(fhicl::ParameterSet const& pset); 
    virtual ~CalWireLBNE35t();
    
    void produce(art::Event& evt);
    void beginJob();
    void endJob();
    void reconfigure(fhicl::ParameterSet const& p);
    
  private:
    
    int          fDataSize;          ///< size of raw data on one wire
    int          fPostsample;        ///< number of postsample bins
    std::string  fDigitModuleLabel;  ///< module that made digits
                                                       ///< constants
    std::string  fSpillName;  ///< nominal spill is an empty string
                              ///< it is set by the DigitModuleLabel
                              ///< ex.:  "daq:preSpill" for prespill data

  protected: 
    
  }; // class CalWireLBNE35t

  DEFINE_ART_MODULE(CalWireLBNE35t)
  
  //-------------------------------------------------
  CalWireLBNE35t::CalWireLBNE35t(fhicl::ParameterSet const& pset)
  {
    this->reconfigure(pset);

    produces< std::vector<recob::Wire> >(fSpillName);
    produces<art::Assns<raw::RawDigit, recob::Wire>>(fSpillName);
  }
  
  //-------------------------------------------------
  CalWireLBNE35t::~CalWireLBNE35t()
  {
  }

  //////////////////////////////////////////////////////
  void CalWireLBNE35t::reconfigure(fhicl::ParameterSet const& p)
  {
    fDigitModuleLabel = p.get< std::string >("DigitModuleLabel", "daq");
    fPostsample       = p.get< int >        ("PostsampleBins");
    
    fSpillName.clear();
    
    size_t pos = fDigitModuleLabel.find(":");
    if( pos!=std::string::npos ) {
      fSpillName = fDigitModuleLabel.substr( pos+1 );
      fDigitModuleLabel = fDigitModuleLabel.substr( 0, pos );
    }
    
  }

  //-------------------------------------------------
  void CalWireLBNE35t::beginJob()
  {  
  }

  //////////////////////////////////////////////////////
  void CalWireLBNE35t::endJob()
  {  
  }
  
  //////////////////////////////////////////////////////
  void CalWireLBNE35t::produce(art::Event& evt)
  {      
    // get the geometry
    art::ServiceHandle<geo::Geometry> geom;

    // get the FFT service to have access to the FFT size
    art::ServiceHandle<util::LArFFT> fFFT;
    int transformSize = fFFT->FFTSize();

    // Get signal shaping service.
    art::ServiceHandle<util::SignalShapingServiceLBNE35t> sss;

    // make a collection of Wires
    std::unique_ptr<std::vector<recob::Wire> > wirecol(new std::vector<recob::Wire>);
    // ... and an association set
    std::unique_ptr<art::Assns<raw::RawDigit,recob::Wire> > WireDigitAssn
      (new art::Assns<raw::RawDigit,recob::Wire>);
    
    // Read in the digit List object(s). 
    art::Handle< std::vector<raw::RawDigit> > digitVecHandle;
    evt.getByLabel(fDigitModuleLabel, fSpillName, digitVecHandle);

    if (!digitVecHandle->size())  return;
    mf::LogInfo("CalWireLBNE35t") << "CalWireLBNE35t:: digitVecHandle size is " << digitVecHandle->size();

    // Use the handle to get a particular (0th) element of collection.
    art::Ptr<raw::RawDigit> digitVec0(digitVecHandle, 0);
        
    //  if (digitVec0->Compression() != raw::kZeroSuppression) {
    //   throw art::Exception(art::errors::UnimplementedFeature)
    //	<< "CalGausHFLBNE only supports zero-suppressed raw digit input!";
    //} // if



    unsigned int dataSize = digitVec0->Samples(); //size of raw data vectors
    //std::cout << "Xin " << dataSize << std::endl;

    raw::ChannelID_t channel = raw::InvalidChannelID; // channel number
    unsigned int bin(0);     // time bin loop variable
    
    filter::ChannelFilter *chanFilt = new filter::ChannelFilter();  

    std::vector<float> holder;                // holds signal data
    std::vector<short> rawadc(transformSize);  // vector holding uncompressed adc values
    std::vector<TComplex> freqHolder(transformSize+1); // temporary frequency data
    
    // loop over all wires    
    wirecol->reserve(digitVecHandle->size());
    for(size_t rdIter = 0; rdIter < digitVecHandle->size(); ++rdIter){ // ++ move
      holder.clear();
      
      // get the reference to the current raw::RawDigit
      art::Ptr<raw::RawDigit> digitVec(digitVecHandle, rdIter);
      channel = digitVec->Channel();

      // skip bad channels
      if(!chanFilt->BadChannel(channel)) {
	holder.resize(transformSize);
	
	// uncompress the data
	raw::Uncompress(digitVec->ADCs(), rawadc, digitVec->Compression());
	
	// loop over all adc values and subtract the pedestal
	for(bin = 0; bin < dataSize; ++bin) 
	  holder[bin]=(rawadc[bin]-digitVec->GetPedestal());

	// Do deconvolution.
	sss->Deconvolute(channel, holder);

      } // end if not a bad channel 
      
      holder.resize(dataSize,1e-5);

      //This restores the DC component to signal removed by the deconvolution.
      if(fPostsample) {
        double average=0.0;
	for(bin=0; bin < (unsigned int)fPostsample; ++bin) 
	  average+=holder[holder.size()-1-bin]/(double)fPostsample;
        for(bin = 0; bin < holder.size(); ++bin) holder[bin]-=average;
      }  
      // Make a single ROI that spans the entire data size
      wirecol->push_back(recob::WireCreator(holder,*digitVec).move());
      // add an association between the last object in wirecol
      // (that we just inserted) and digitVec
      if (!util::CreateAssn(*this, evt, *wirecol, digitVec, *WireDigitAssn, fSpillName)) {
        throw art::Exception(art::errors::InsertFailure)
          << "Can't associate wire #" << (wirecol->size() - 1)
          << " with raw digit #" << digitVec.key();
      } // if failed to add association
    }
    
    if(wirecol->size() == 0)
      mf::LogWarning("CalWireLBNE35t") << "No wires made for this event.";

    evt.put(std::move(wirecol), fSpillName);
    evt.put(std::move(WireDigitAssn), fSpillName);


    delete chanFilt;
    return;
  }
  
} // end namespace caldata
