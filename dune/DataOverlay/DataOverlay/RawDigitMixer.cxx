#ifndef OVERLAY_DATAOVERLAY_RAWDIGITMIXER_CXX
#define OVERLAY_DATAOVERLAY_RAWDIGITMIXER_CXX

#include "RawDigitMixer.h"
#include <limits>
#include <iostream>
#include <stdexcept>
#include <algorithm>

void mix::RawDigitMixer::DeclareData(std::vector<raw::RawDigit> const& dataVector){

  fChannelIndexMap.clear();
  fOutputWaveforms.clear();
  fOutputWaveforms.resize(dataVector.size());

  for(size_t i_rd=0; i_rd<dataVector.size(); i_rd++){
    
    std::vector<short> waveform;
    if(dataVector[i_rd].Compression()!=raw::kNone)
      {
	raw::Uncompress(dataVector[i_rd].ADCs(),
			waveform,
			dataVector[i_rd].GetPedestal(),
			dataVector[i_rd].Compression());
      }
    else 
      {
	waveform = dataVector[i_rd].ADCs();
      }    

    fChannelIndexMap[dataVector[i_rd].Channel()] = i_rd;
    
    //initialize adc vector for output
    fOutputWaveforms[i_rd].waveform.resize(dataVector[i_rd].Samples(),0);
    fOutputWaveforms[i_rd].channel = dataVector[i_rd].Channel();
    fOutputWaveforms[i_rd].ped     = dataVector[i_rd].GetPedestal();
    fOutputWaveforms[i_rd].sigma   = dataVector[i_rd].GetSigma();
    
    
    //do the steps for filling that output vector
    //This is data, so set scales to one.
    fRDAdderAlg.SetScaleInputs(1.0,1.0);
    fRDAdderAlg.SetPedestalInputs(0.0,0.0);
    fRDAdderAlg.AddRawDigits(waveform,fOutputWaveforms[i_rd].waveform);
  }
  
}

void mix::RawDigitMixer::Mix(std::vector<raw::RawDigit> const& mcVector,
			     std::unordered_map<raw::ChannelID_t,float> const& scale_map){


  for( auto const& rd : mcVector){

    if(rd.Compression()!=raw::kNone)
      {
	throw std::runtime_error("Error in RawDigitMixer::Mix : Compressed MC waveforms not supported. Turn it off and try again."); 
      }

    auto it_ch = fChannelIndexMap.find(rd.Channel());

    //if this channel is not in the data, skip this channel!
    if(it_ch==fChannelIndexMap.end())
      continue;

    size_t i_output = it_ch->second;

    fRDAdderAlg.SetPedestalInputs(rd.GetPedestal(),0.0);
    fRDAdderAlg.SetScaleInputs(scale_map.at(rd.Channel()),1.0);
    fRDAdderAlg.SetStuckBitRetentionMethod(_stuckRetention);
    
    // hard-code 15000 here because the production split/stitched data files should have 15000 ticks exactly
    if (fOutputWaveforms[i_output].waveform.size() == 15000 && fOutputWaveforms[i_output].waveform.size() > rd.ADCs().size()+1)
      {
	size_t offset = fOutputWaveforms[i_output].waveform.size()-rd.Samples()-1;
	std::vector<short> const data_trimmed(fOutputWaveforms[i_output].waveform.begin()+offset,
					fOutputWaveforms[i_output].waveform.begin()+offset+rd.Samples());
	//std::vector<short> const& mc(rd.ADCs());
	fRDAdderAlg.AddRawDigits(rd.ADCs(),data_trimmed,fOutputWaveforms[i_output].waveform);
      }
    else
      {
	throw std::runtime_error("Error in RawDigitMixer::Mix : The waveform sizes aren't what they're expected to be. I need real data to be 15000 ticks and MC to be 5200 ticks. If this has changed, then recompile, or come up with a more clever solution than throwing an exception.");
      }
    
  }

}

void mix::RawDigitMixer::FillRawDigitOutput(std::vector<raw::RawDigit> & output){

  for(auto const& rd_i : fOutputWaveforms){
    
    //now emplace back onto output collection...
    output.emplace_back(rd_i.channel,
			rd_i.waveform.size(),
			rd_i.waveform);
    
    //set pedestal and rms to be same as data
    output.back().SetPedestal(rd_i.ped,rd_i.sigma);
  }
  
}


#endif
