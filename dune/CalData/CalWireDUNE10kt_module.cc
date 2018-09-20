////////////////////////////////////////////////////////////////////////
//
// CalWireDUNE10kt class
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
#include <iostream>
#include <iomanip>

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
#include "larevt/Filters/ChannelFilter.h"
#include "lardataobj/RawData/RawDigit.h"
#include "lardataobj/RawData/raw.h"
#include "lardataobj/RecoBase/Wire.h"
#include "lardata/ArtDataHelper/WireCreator.h"
#include "lardata/Utilities/LArFFT.h"
#include "lardata/Utilities/AssociationUtil.h"
#include "larevt/CalibrationDBI/Interface/ChannelStatusService.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"

// Dune includes.
#include "dune/Utilities/SignalShapingServiceDUNE.h"
#include "dune/DuneInterface/ChannelMappingService.h"

using std::string;
using std::cout;
using std::endl;
using std::setw;
using lariov::ChannelStatusService;
using lariov::ChannelStatusProvider;

///creation of calibrated signals on wires
namespace caldata {

class CalWireDUNE10kt : public art::EDProducer {

public:
    
  // create calibrated signals on wires. this class runs 
  // an fft to remove the electronics shaping.     
  explicit CalWireDUNE10kt(fhicl::ParameterSet const& pset); 
  virtual ~CalWireDUNE10kt();
  
  void produce(art::Event& evt); 
  void beginJob();
  void endJob();
  void reconfigure(fhicl::ParameterSet const& p);
 
private:
    
  //int            fDataSize;         ///< size of raw data on one wire // unused
  int            fPostsampleBins;   ///< number of postsample bins
  int            fDoBaselineSub;    ///< number of postsample bins
  std::string    fDigitModuleLabel; ///< module that made digits
                                                     ///< constants
  std::string    fSpillName;        ///< nominal spill is an empty string
                                    ///< it is set by the DigitModuleLabel
                                    ///< ex.:  "daq:preSpill" for prespill data
  unsigned short fPreROIPad;        ///< ROI padding
  unsigned short fPostROIPad;       ///< ROI padding
  double         fSigThrFact;       ///< Signal shreshold factor 
  bool           fSkipBadChannels;  ///< Skip bad channels.
  int            fLogLevel;         ///< Log level: 0=none, 1=init
  void SubtractBaseline(std::vector<float>& holder);
  
}; // class CalWireDUNE10kt

DEFINE_ART_MODULE(CalWireDUNE10kt)
  
//////////////////////////////////////////////////////

CalWireDUNE10kt::CalWireDUNE10kt(fhicl::ParameterSet const& pset) {
  this->reconfigure(pset);

  produces< std::vector<recob::Wire> >(fSpillName);
  produces<art::Assns<raw::RawDigit, recob::Wire>>(fSpillName);
}
  
//////////////////////////////////////////////////////

CalWireDUNE10kt::~CalWireDUNE10kt() { }

//////////////////////////////////////////////////////

void CalWireDUNE10kt::reconfigure(fhicl::ParameterSet const& pset) {
  const string myname = "CalWireDUNE10kt::reconfigure: ";
  std::vector<unsigned short> uin;    std::vector<unsigned short> vin;
  std::vector<unsigned short> zin;

  fDigitModuleLabel = pset.get<std::string>("DigitModuleLabel", "daq");
  fPostsampleBins   = pset.get<int>("PostsampleBins");
  fDoBaselineSub    = pset.get<bool>("DoBaselineSub");
  uin               = pset.get<std::vector<unsigned short>>("PlaneROIPad");
  fSigThrFact       = pset.get<double>("SigThrFact", 3.0);
  fSkipBadChannels = false;
  pset.get_if_present<int>("SkipBadChannels", fLogLevel);
  fLogLevel         = 1;
  pset.get_if_present<int>("LogLevel", fLogLevel);
  // put the ROI pad sizes into more convenient vectors
  fPreROIPad  = uin[0];
  fPostROIPad = uin[1];
  fSpillName.clear();
  
  size_t pos = fDigitModuleLabel.find(":");
  if( pos!=std::string::npos ) {
    fSpillName = fDigitModuleLabel.substr( pos+1 );
    fDigitModuleLabel = fDigitModuleLabel.substr( 0, pos );
  }

  if ( fLogLevel >= 1 ) {
    cout << myname << "  DigitModuleLabel: " << fDigitModuleLabel << endl;
    cout << myname << "    PostsampleBins: " << fPostsampleBins << endl;
    cout << myname << "     DoBaselineSub: " << fDoBaselineSub << endl;
    cout << myname << "    PlaneROIPad[0]: " << fPreROIPad << endl;
    cout << myname << "    PlaneROIPad[1]: " << fPostROIPad << endl;
    cout << myname << "        SigThrFact: " << fSigThrFact << endl;
    cout << myname << "   SkipBadChannels: " << fSkipBadChannels << endl;
    cout << myname << "          LogLevel: " << fLogLevel << endl;
  }
  
}

//////////////////////////////////////////////////////

void CalWireDUNE10kt::beginJob() { }

//////////////////////////////////////////////////////

void CalWireDUNE10kt::endJob() { }
  
//////////////////////////////////////////////////////

void CalWireDUNE10kt::produce(art::Event& evt) {      
  const string myname = "CalWireDUNE10kt::produce: ";

  // get the geometry
  art::ServiceHandle<geo::Geometry> geom;

  // get the FFT service to have access to the FFT size
  art::ServiceHandle<util::LArFFT> fFFT;
  int transformSize = fFFT->FFTSize();
  if ( fLogLevel >= 2 ) cout << myname << "FFT size: " << transformSize << endl;

  // Get signal shaping service.
  art::ServiceHandle<util::SignalShapingServiceDUNE> sss;
  double DeconNorm = sss->GetDeconNorm();
 
  // Fetch the channel mapping and channel status services.
  const ChannelMappingService* pchanmap = nullptr;
  const ChannelStatusProvider* pcsp = nullptr;
  if ( fSkipBadChannels ) {
    art::ServiceHandle<ChannelMappingService> hchanmap;
    pchanmap = &*hchanmap;
    art::ServiceHandle<ChannelStatusService> cssHandle;
    pcsp = cssHandle->GetProviderPtr();
    if ( pcsp == nullptr ) {
      cout << myname << "ERROR: Channel status service not found." << endl;
      abort();
    }
  }

  // make a collection of Wires
  std::unique_ptr<std::vector<recob::Wire> > wirecol(new std::vector<recob::Wire>);
  // ... and an association set
  std::unique_ptr<art::Assns<raw::RawDigit,recob::Wire> > WireDigitAssn
    (new art::Assns<raw::RawDigit,recob::Wire>);
  
  // Read in the digit List object(s). 
  art::Handle< std::vector<raw::RawDigit> > digitVecHandle;
  evt.getByLabel(fDigitModuleLabel, fSpillName, digitVecHandle);

  if (!digitVecHandle->size())  return;
  mf::LogInfo("CalWireDUNE10kt") << "CalWireDUNE10kt:: digitVecHandle size is " << digitVecHandle->size();

  // Use the handle to get a particular (0th) element of collection.
  art::Ptr<raw::RawDigit> digitVec0(digitVecHandle, 0);
      
  unsigned int dataSize = digitVec0->Samples(); //size of raw data vectors
  if ( fLogLevel >= 2 ) cout << myname << "Expected raw data size: " << dataSize << endl;
  int readoutwindowsize = art::ServiceHandle<detinfo::DetectorPropertiesService>()->provider()->ReadOutWindowSize();
  if (int(dataSize) != readoutwindowsize){
    throw art::Exception(art::errors::Configuration)
      << "ReadOutWindowSize "<<readoutwindowsize<<" does not match data size "<<dataSize<<". Please set services.DetectorPropertiesService.NumberTimeSamples and services.DetectorPropertiesService.ReadOutWindowSize in fcl file to "<<dataSize;
  }
  
  raw::ChannelID_t channel = raw::InvalidChannelID; // channel number
  unsigned int bin(0);     // time bin loop variable
  
  std::vector<float> holder;                // holds signal data
  std::vector<short> rawadc(transformSize);  // vector holding uncompressed adc values
  std::vector<TComplex> freqHolder(transformSize+1); // temporary frequency data
  
  // loop over all wires    
  wirecol->reserve(digitVecHandle->size());
  int maxdump = 1;
  int ndump = 0;
  for ( size_t rdIter = 0; rdIter < digitVecHandle->size(); ++rdIter ) { // ++ move
    holder.clear();
    
    // get the reference to the current raw::RawDigit
    art::Ptr<raw::RawDigit> digitVec(digitVecHandle, rdIter);
    channel = digitVec->Channel();
    unsigned int onlineChannel = channel;
    if ( pchanmap != nullptr ) onlineChannel = pchanmap->online(channel);

    // Skip bad channels
    if ( fSkipBadChannels && pcsp->IsBad(onlineChannel) ) {
      if ( fLogLevel >= 2 ) cout << myname << "Skipping bad channel " << channel << endl;
    } else {

      holder.resize(transformSize);

      // uncompress the data
      raw::Uncompress(digitVec->ADCs(), rawadc, digitVec->GetPedestal(), digitVec->Compression());
      for ( unsigned int ibin=0; ibin<rawadc.size(); ++ibin ) rawadc[ibin] -= digitVec->GetPedestal();
      if ( fLogLevel >= 2 && ndump<maxdump ) {
        cout << myname << "Processing channel " << channel << endl;
        if ( fLogLevel == 2 ) cout << myname << "  Uncompressed raw data size: " << rawadc.size() << endl;
        if ( fLogLevel >= 3 ) {
          cout << "  Raw data vector has " << rawadc.size() << " entries." << endl;
          for ( unsigned int ibin=0; ibin<rawadc.size(); ++ibin ) {
            cout << setw(8) << ibin << ": " << rawadc[ibin] << endl;
          }
        }
        ++ndump;
      }

      for ( bin=0; bin<dataSize; ++bin ) holder[bin] = rawadc[bin];
      //Xin fill the remaining bin with data
      for ( bin=dataSize; bin<holder.size(); ++bin) holder[bin] = holder[bin-dataSize];

      // Do deconvolution.
      sss->Deconvolute(channel, holder);
      for ( bin = 0; bin < holder.size(); ++bin) holder[bin]=holder[bin]/DeconNorm;
    } // end if not a bad channel 
      
    holder.resize(dataSize,1e-5);

    // Restore the DC component to signal removed by the deconvolution.
    if ( fPostsampleBins ) {
      double average=0.0;
      for ( bin=0; bin < (unsigned int)fPostsampleBins; ++bin)
        average+=holder[holder.size()-1-bin]/(double)fPostsampleBins;
      for ( bin = 0; bin < holder.size(); ++bin ) holder[bin]-=average;
    }  

    // Adaptive baseline subtraction
    if ( fDoBaselineSub ) SubtractBaseline(holder);
      
    // Work out the ROI 
    recob::Wire::RegionsOfInterest_t ROIVec;
    std::vector<std::pair<unsigned int, unsigned int>> holderInfo;
    std::vector<std::pair<unsigned int, unsigned int>> rois;
      
    double max = 0;
    double deconNoise = sss->GetDeconNoise(channel);
    // find out all ROI
    unsigned int roiStart = 0;
    for ( bin = 0; bin < dataSize; ++bin ) {
      double SigVal = holder[bin];
      if (SigVal > max) max = SigVal;
      if(roiStart == 0) {
        if (SigVal > fSigThrFact*deconNoise) roiStart = bin; // n sigma above noise
      }else{
        if (SigVal < deconNoise){
          rois.push_back(std::make_pair(roiStart, bin));
          roiStart = 0;
        }
      }
    }
    if (roiStart!=0){
      rois.push_back(std::make_pair(roiStart, dataSize-1));
      roiStart = 0;
    }
      
    // pad them
    // if (channel==512){
    // 	for (bin = 0; bin< holder.size();++bin){
    // 	  if (fabs(holder[bin]) > 2)
    // 	      std::cout << "Xin1: " << holder[bin] << std::endl;
    // 	}
    // }
    //std::cout << "Xin: "  << max << " "<< channel << " " << deconNoise << " " << rois.size() << std::endl;

    if ( rois.size() == 0 ) continue;
    holderInfo.clear();
    for( unsigned int ii = 0; ii<rois.size(); ++ii ) {
      // low ROI end
      int low = rois[ii].first - fPreROIPad;
      if(low < 0) low = 0;
      rois[ii].first = low;
      // high ROI end
      unsigned int high = rois[ii].second + fPostROIPad;
      if(high >= dataSize) high = dataSize-1;
      rois[ii].second = high;
    }
    // merge them
    if ( rois.size() >= 1 ) {
      // temporary vector for merged ROIs

      for ( unsigned int ii = 0; ii<rois.size(); ++ii ) {
        unsigned int roiStart = rois[ii].first;
        unsigned int roiEnd = rois[ii].second;
        int flag1 = 1;
        unsigned int jj=ii+1;
        while ( flag1 ) {	
          if ( jj<rois.size() ) {
            if( rois[jj].first <= roiEnd ) {
              roiEnd = rois[jj].second;
	      ii = jj;
              jj = ii+1;
            } else {
	      flag1 = 0;
            }
          } else {
            flag1 = 0;
          }
        }
        std::vector<float> sigTemp;
        for(unsigned int kk = roiStart; kk < roiEnd; ++kk) {
          sigTemp.push_back(holder[kk]);
        } // jj
        //	  std::cout << "Xin: " << roiStart << std::endl;
        ROIVec.add_range(roiStart, std::move(sigTemp));
        //trois.push_back(std::make_pair(roiStart,roiEnd));	    
      }
    }// else{
    // 	unsigned int roiStart = rois[0].first;
    // 	unsigned int roiEnd = rois[0].second;
    // 	std::vector<float> sigTemp;
    // 	  for(unsigned int kk = roiStart; kk < roiEnd; ++kk) {
    // 	    sigTemp.push_back(holder[kk]);
    // 	  } // jj
    // 	  //	  std::cout << "Xin: " << roiStart << std::endl;
    // 	  ROIVec.add_range(roiStart, std::move(sigTemp));
    // }
    
    // save them
    wirecol->push_back(recob::WireCreator(std::move(ROIVec),*digitVec).move());

    // Make a single ROI that spans the entire data size
    //wirecol->push_back(recob::WireCreator(holder,*digitVec).move());
    // add an association between the last object in wirecol
    // (that we just inserted) and digitVec
    if ( ! util::CreateAssn(*this, evt, *wirecol, digitVec, *WireDigitAssn, fSpillName) ) {
      throw art::Exception(art::errors::ProductRegistrationFailure)
        << "Can't associate wire #" << (wirecol->size() - 1)
        << " with raw digit #" << digitVec.key();
    } // if failed to add association
  }
    
  if ( wirecol->size() == 0 )
    mf::LogWarning("CalWireDUNE10kt") << "No wires made for this event.";

  evt.put(std::move(wirecol), fSpillName);
  evt.put(std::move(WireDigitAssn), fSpillName);

  return;
}
  
//////////////////////////////////////////////////////

void CalWireDUNE10kt::SubtractBaseline(std::vector<float>& holder) {

  float min = 0,max=0;
  for (unsigned int bin = 0; bin < holder.size(); bin++){
    if (holder[bin] > max) max = holder[bin];
    if (holder[bin] < min) min = holder[bin];
  }
  int nbin = max - min;
  if (nbin!=0){
    TH1F *h1 = new TH1F("h1","h1",nbin,min,max);
    for (unsigned int bin = 0; bin < holder.size(); bin++){
      h1->Fill(holder[bin]);
    }
    float ped = h1->GetMaximum();
    float ave=0,ncount = 0;
    for (unsigned int bin = 0; bin < holder.size(); bin++){
      if (fabs(holder[bin]-ped)<2){
        ave +=holder[bin];
        ncount ++;
      }
    }
    if (ncount==0) ncount=1;
    ave = ave/ncount;
    for (unsigned int bin = 0; bin < holder.size(); bin++){
      holder[bin] -= ave;
    }
    h1->Delete();
  }
}

//////////////////////////////////////////////////////

} // end namespace caldata
