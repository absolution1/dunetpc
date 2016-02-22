#ifndef FILTERWF_H
#define FILTERWF_H
////////////////////////////////////////////////////////////////////////
//
// FilterWF_module
//
// bkirby@bnl.gov
//
// Updated by m.wallbank@sheffield.ac.uk
//
////////////////////////////////////////////////////////////////////////

#include <string>
#include <vector>
#include <iostream>
#include <TProfile.h>
#include <TTree.h>
#include <TMath.h>
#include <TFile.h>
//#include <TFileDirectory.h>

#include "art/Framework/Principal/Handle.h" 
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "messagefacility/MessageLogger/MessageLogger.h" 
#include "art/Framework/Services/Optional/detail/TH1AddDirectorySentry.h"
//#include "art/Framework/Services/Optional/detail/RootDirectorySentry.h"

#include "lardata/RawData/RawDigit.h"
//#include "DetectorInfoServices/DetectorPropertiesService.h"

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "art/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include <memory>
#include <iostream>
#include <fstream>
#include <sstream>

namespace lbne {
  class FilterWF;
}

class lbne::FilterWF : public art::EDProducer {

public:
  explicit FilterWF(fhicl::ParameterSet const& pset);
  virtual ~FilterWF();

  virtual void produce(art::Event& evt);
  void reconfigure(fhicl::ParameterSet const& pset);
  void beginRun(art::Run const& run);
  void endRun(art::Run const& run);

  void buildChannelMap(std::string const& mapFile);
    
private:

  //******************************
  //Variables Taken from FHICL File
  std::string       fRawDigitModuleLabel, fRawDigitModuleInstance, fChannelMapFile;   //label for rawdigit module
  //std::string fOutputModuleLabel;

  std::map<int, int> fOfflineToOnlineChannel, fOnlineToOfflineChannel;

}; //end class GetWF


//-------------------------------------------------------------------
lbne::FilterWF::FilterWF(fhicl::ParameterSet const& pset) {
  this->reconfigure(pset);
  // Only need to build channel map once -- do it in constructor
  this->buildChannelMap(fChannelMapFile);
  produces< std::vector<raw::RawDigit> >();
}

//-------------------------------------------------------------------
lbne::FilterWF::~FilterWF(){}

//-------------------------------------------------------------------
void lbne::FilterWF::reconfigure(fhicl::ParameterSet const& pset){
  fRawDigitModuleLabel = pset.get<std::string>("RawDigitModuleLabel");
  fRawDigitModuleInstance = pset.get<std::string>("RawDigitModuleInstance");
  fChannelMapFile = pset.get<std::string>("ChannelMapFile");
  //fOutputModuleLabel = "filterwf";
}

//-------------------------------------------------------------------
void lbne::FilterWF::beginRun(art::Run const& run){
  return;
}

//-------------------------------------------------------------------
void lbne::FilterWF::endRun(art::Run const& run){
  return;
}
  
//-------------------------------------------------------------------
void lbne::FilterWF::produce(art::Event& evt){

  art::Handle< std::vector<raw::RawDigit> > rawDigitHandle;
  evt.getByLabel(fRawDigitModuleLabel, fRawDigitModuleInstance,rawDigitHandle);
  std::vector<raw::RawDigit> const& rawDigitVector(*rawDigitHandle);

  std::cout << "Event # " << evt.event() << std::endl;

  //define temporary vector to hold filtered waveforms
  const unsigned int n_channels = rawDigitVector.size();    
  std::vector< std::vector<short> > filterWf(n_channels); // MW: I had to make this a vector of shorts to make it compile (the RawDigit constructor does not accept ADC as float)

  //loop over channels - get pedestal mean (no stuck ADC codes)
  int maxNumBins = -1;
  std::vector<float> meanVec;
  for(unsigned int ich=0; ich<n_channels; ich++){
    const size_t n_samp = rawDigitVector.at(ich).NADC();
    if (maxNumBins<0) maxNumBins = n_samp;
    const int offlineChannel = rawDigitVector.at(ich).Channel();
    const int onlineChannel = fOfflineToOnlineChannel.at(offlineChannel);
    if ((int)ich != onlineChannel){
      std::cout << "WARNING! Problem in first channel loop. Online channel " << onlineChannel << " (offline channel " << offlineChannel << ") is not equal to loop index, " << ich << std::endl;
      abort();
    }
    double mean = 0;
    int count = 0;
    for(unsigned int s = 0 ; s < n_samp ; s++){
      short adc = rawDigitVector.at(ich).ADC(s);
      if( adc < 10 ) 
	continue;
      if( (adc & 0x3F ) == 0x0 || (adc & 0x3F ) == 0x3F )
	continue;
      mean += adc;
      count++;
    }
    if( count > 0 )
      mean = mean / (double) count;
    meanVec.push_back(mean);
  }

  //define set of induction and collection channels in each regulator group
  std::vector<int> indCh;
  std::vector<int> colCh;
  for(int i = 0 ; i < 64 ; i++){
    if( i < 32 || i == 32 || i == 47 || i == 48 || i == 63 )
      indCh.push_back(i);
    else
      colCh.push_back(i);
  }

  //derive correction factors - require raw adc waveform and pedestal for each channel
  std::vector<Double_t> corrVals;
  //loop through time slices
  //int maxNumBins = detprop->ReadOutWindowSize();
  for(int s = 0 ; s < maxNumBins ; s++ ){
    //loop through regulator groups
    for(int g = 0 ; g < 16*2 ; g++){
      int baseCh = g*64;
      //get induction plane correction
      if(1){
	corrVals.clear();
	for(unsigned int c = 0 ; c < indCh.size() ; c++){
	  int ch = baseCh + indCh.at(c);
	  int adc = rawDigitVector.at(ch).ADC(s); // NB/ MW: possible bug; changed .at(c) to .at(ch)
	  if( (adc & 0x3F ) == 0x0 || (adc & 0x3F ) == 0x3F )
	    continue;
	  if( adc < 10  ) //skip "sample dropping" problem
	    continue;
	  double mean = meanVec.at(ch);
	  if( mean < 10 )
	    continue;
	  corrVals.push_back(adc - mean);
	}
			
	unsigned int corrValSize = corrVals.size();
	sort(corrVals.begin(),corrVals.end());
	double correction = 0;
	if(corrValSize < 2)
	  correction = 0.0;
	else if((corrValSize % 2) == 0)
	  correction = (corrVals[corrValSize/2] + corrVals[(corrValSize/2)-1])/2.0;
	else
	  correction = corrVals[(corrValSize-1)/2];
	
	for(unsigned int c = 0 ; c < indCh.size() ; c++){
	  int ch = baseCh + indCh.at(c);
	  int adc = rawDigitVector.at(ch).ADC(s); //again c->ch
	  double newAdc = adc - correction;
//	  if( (adc & 0x3F ) == 0x0 || (adc & 0x3F ) == 0x3F )
//	    newAdc = 0;
//	  if( adc < 10  ) //skip "sample dropping" problem
//	    newAdc = 0;
	  filterWf.at(ch).push_back(newAdc);
	}
      }//end induction plane correction
      //get collection plane correction
      if(1){
	corrVals.clear();
	for(unsigned int c = 0 ; c < colCh.size() ; c++){
	  int ch = baseCh + colCh.at(c);
	  int adc = rawDigitVector.at(ch).ADC(s); // MW: fix 'bug'(?) again
	  if( (adc & 0x3F ) == 0x0 || (adc & 0x3F ) == 0x3F )
	    continue;
	  if( adc < 10  ) //skip "sample dropping" problem
	    continue;
	  double mean = meanVec.at(ch);
	  if( mean < 10 )
	    continue;
	  corrVals.push_back(adc - mean);
	}
			
	unsigned int corrValSize = corrVals.size();
	sort(corrVals.begin(),corrVals.end());
	double correction = 0;
	if(corrValSize < 2)
	  correction = 0.0;
	else if((corrValSize % 2) == 0)
	  correction = (corrVals[corrValSize/2] + corrVals[(corrValSize/2)-1])/2.0;
	else
	  correction = corrVals[(corrValSize-1)/2];

	for(unsigned int c = 0 ; c < colCh.size() ; c++){
	  int ch = baseCh + colCh.at(c);
	  int adc = rawDigitVector.at(ch).ADC(s); // MW: and again
	  double newAdc = adc - correction;
//	  if( (adc & 0x3F ) == 0x0 || (adc & 0x3F ) == 0x3F )
//	    newAdc = 0;
//	  if( adc < 10  ) //skip "sample dropping" problem
//	    newAdc = 0;
	  filterWf.at(ch).push_back(newAdc);
	}
      }//end collection plane correction
    }//end loop over regulator groups
  }//end loop over samples

  //loop over channels - save filtered waveforms into digits
  std::vector<raw::RawDigit> filterRawDigitVector;
  for(unsigned int ich=0; ich<n_channels; ich++){
    short fChan = rawDigitVector.at(ich).Channel();
    raw::RawDigit theRawDigit(fChan, filterWf.at(ich).size(), filterWf.at(ich));//, raw::kNone );
    filterRawDigitVector.push_back(theRawDigit);            // add this digit to the collection
  }

  //save filtered waveforms to event
  evt.put(std::make_unique<decltype(filterRawDigitVector)>(std::move(filterRawDigitVector)));//, fOutputModuleLabel);
}

void lbne::FilterWF::buildChannelMap(std::string const& channelMapFile) {

  /// Builds TPC channel map from the map txt file

  fOnlineToOfflineChannel.clear();
  fOfflineToOnlineChannel.clear();

  int onlineChannel;
  int offlineChannel;
    
  std::string fullname;
  cet::search_path sp("FW_SEARCH_PATH");
  sp.find_file(channelMapFile, fullname);
    
  if (fullname.empty())
    mf::LogWarning("DAQToOffline") << "Input TPC channel map file " << channelMapFile << " not found in FW_SEARCH_PATH.  Using online channel numbers!" << std::endl;

  else {
    mf::LogVerbatim("DAQToOffline") << "Build TPC Online->Offline channel Map from " << fullname;
    std::ifstream infile(fullname);
    while (infile.good()) {
      infile >> onlineChannel >> offlineChannel;
      fOnlineToOfflineChannel.insert(std::make_pair(onlineChannel,offlineChannel));
      fOfflineToOnlineChannel.insert(std::make_pair(offlineChannel,onlineChannel));
      mf::LogVerbatim("DAQToOffline") << "   " << onlineChannel << " -> " << offlineChannel;
    }
  }

}

DEFINE_ART_MODULE(lbne::FilterWF)

#endif //FILTERWF_H
