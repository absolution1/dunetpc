// art 
#include "art/Framework/Core/FileBlock.h"
#include "art/Framework/Core/InputSourceMacros.h"
#include "art/Framework/Core/ProductRegistryHelper.h"
#include "art/Framework/IO/Root/rootNames.h"
#include "art/Framework/IO/Sources/Source.h"
#include "art/Framework/IO/Sources/SourceHelper.h"
#include "art/Framework/IO/Sources/SourceTraits.h"
#include "art/Framework/IO/Sources/put_product_in_principal.h"
#include "art/Persistency/Provenance/EventID.h"
#include "art/Persistency/Provenance/MasterProductRegistry.h"
#include "art/Persistency/Provenance/RunID.h"
#include "art/Persistency/Provenance/SubRunID.h"
#include "art/Utilities/InputTag.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "art/Persistency/Provenance/EventAuxiliary.h"
#include "SimulationBase/MCParticle.h"
#include "Simulation/SimChannel.h"

//Pedestal stuff...
#include "CalibrationDBI/Interface/IDetPedestalService.h"
#include "CalibrationDBI/Interface/IDetPedestalProvider.h"
#include "CalibrationDBI/Interface/IChannelStatusService.h"
#include "CalibrationDBI/Interface/IChannelStatusProvider.h"
#include "dune/RunHistory/DetPedestalDUNE.h"
#include "cetlib/getenv.h"

// artdaq 
#include "artdaq-core/Data/Fragments.hh"

// lardata
#include "RawData/RawDigit.h"
#include "RawData/ExternalTrigger.h"
#include "Utilities/TimeService.h"

// dune
#include "tpcFragmentToRawDigits.h"
#include "SSPReformatterAlgs.h"
#include "PennToOffline.h"

// C++ 
#include <functional>
#include <memory>
#include <regex>
#include <string>
#include <vector>

// ROOT
#include "TTree.h"
#include "TFile.h"

using raw::RawDigit;
using raw::OpDetWaveform;
using raw::ExternalTrigger;
using std::vector;
using std::string;

// TODO:
//  Handle cases of missing data products.  What if SSP or Penn or RCE data are just not there?
//         If SSP or Penn not there it is fine, but if no RCE then it won't know how when to stop making the event...
//  Deal with ZS data -- currently this assumes non-ZS data
//===================================================================================================================

namespace {
  
  struct LoadedWaveforms {
    
    LoadedWaveforms() : waveforms() {}
    vector<OpDetWaveform> waveforms;
    
    void load( vector<OpDetWaveform> const & v, int fDebugLevel ) {
      waveforms = v;
    }
    
    // not really used 
    bool empty() const { 
      if (waveforms.size() == 0) return true;
      return false;
    }

    // Clear waveforms vector.
    void clear(int fDebugLevel) {
      if (fDebugLevel > 3) std::cout << "Clearing LoadedWaveforms." << std::endl;
      waveforms.clear();
      empty();
    }

    std::vector<OpDetWaveform> TakeAll() {
      return waveforms;
    }
    
    void findinrange(std::vector<OpDetWaveform> &wbo, 
                     lbne::TpcNanoSlice::Header::nova_timestamp_t first_timestamp,
                     lbne::TpcNanoSlice::Header::nova_timestamp_t last_timestamp,
		     lbne::TpcNanoSlice::Header::nova_timestamp_t event_timestamp,
                     double NovaTicksPerSSPTick,
		     int fDebugLevel) {
      int hh = 0;
      for (auto wf : waveforms) { // see if any waveforms have pieces inside this TPC boundary
	lbne::TpcNanoSlice::Header::nova_timestamp_t WaveformTimestamp = wf.TimeStamp() * NovaTicksPerSSPTick;
        if ( hh < 5 && fDebugLevel > 2) { // Just want to write out the first few waveforms
	  std::cout << "Looking at waveform number " << hh << ". It was on channel " << wf.ChannelNumber() << " at time " << WaveformTimestamp
		    << ". The times I passed were " << first_timestamp << " and " << last_timestamp
		    << std::endl;
	}
	if (WaveformTimestamp <= last_timestamp && WaveformTimestamp >= first_timestamp) {
	  lbne::TpcNanoSlice::Header::nova_timestamp_t NewNovaTime = WaveformTimestamp - event_timestamp;
	  double NewTime = NewNovaTime / NovaTicksPerSSPTick;
	  raw::OpDetWaveform odw( NewTime, wf.ChannelNumber(), wf.size() );
          if (fDebugLevel > 3)
	    std::cout << "Pushing back waveform " << hh << " on channel " << odw.ChannelNumber() << " at time " << odw.TimeStamp() << " ("<< odw.TimeStamp()*NovaTicksPerSSPTick <<")" << std::endl;
	  for (size_t WaveSize = 0; WaveSize < wf.size(); ++WaveSize ) {
	    odw.emplace_back(wf[WaveSize]);
          }
	  wbo.emplace_back(std::move(odw));
        }
        ++hh;
      } // auto waveforms
      if (fDebugLevel > 1) std::cout << "At the end of Waveform findinrange, wbo has size " << wbo.size() << std::endl;
    } // findinrange
    
    //=======================================================================================
    bool PhotonTrigger(lbne::TpcNanoSlice::Header::nova_timestamp_t this_timestamp, double fWaveformADCThreshold, int fWaveformADCsOverThreshold,
                       double fWaveformADCWidth, double NovaTicksPerSSPTick, int fDebugLevel ) { // Triggering on photon detectors
      int HighADCWaveforms = 0;

      double SumWaveforms = 0;
      int ADCs = 0;
      int More1550 = 0;
      int More1750 = 0;
      int More2000 = 0;
      double Lowest = 1e7;
      double Biggest = 0;
      for (auto wf : waveforms) { // see if any waveforms have pieces inside this TPC boundary
	lbne::TpcNanoSlice::Header::nova_timestamp_t WaveformTimestamp = wf.TimeStamp() * NovaTicksPerSSPTick;
	lbne::TpcNanoSlice::Header::nova_timestamp_t StartTimestamp = WaveformTimestamp - ( fWaveformADCWidth * NovaTicksPerSSPTick * 0.5);
	lbne::TpcNanoSlice::Header::nova_timestamp_t EndTimestamp   = WaveformTimestamp + ( fWaveformADCWidth * NovaTicksPerSSPTick * 0.5);
	if (this_timestamp <= EndTimestamp && this_timestamp >= StartTimestamp) {
          for (int ii=0;ii<(int)wf.size();++ii) {
            if ((int)wf.Waveform()[ii] >  fWaveformADCThreshold) {
              ++HighADCWaveforms;
              
              SumWaveforms += (int)wf.Waveform()[ii];
              if ( wf.Waveform()[ii] < Lowest ) Lowest = wf.Waveform()[ii];
              if ( wf.Waveform()[ii] > Biggest) Biggest= wf.Waveform()[ii];
              if ( wf.Waveform()[ii] > 1550 ) ++More1550;
              if ( wf.Waveform()[ii] > 1750 ) ++More1750;
              if ( wf.Waveform()[ii] > 2000 ) ++More2000;
              ++ADCs;
            } // If ADC count more than threshold (1550)
          } // Loop over wf ADC's
        } // If times match
      } // Loop over waveforms
      if (fDebugLevel > 3) {
	std::cout << "Lowest " << Lowest << ", Biggest " << Biggest << ", Average " << SumWaveforms / ADCs
		  << ", Total " << ADCs << ", More1550 " << More1550 << ", More1750 " << More1750 << ", More2000 " << More2000
		  << std::endl;
      }
      if ( HighADCWaveforms > fWaveformADCsOverThreshold ) return true;
      else return false;
    } // Photon Trigger
  }; // Loaded Waveforms
  
  //=============================== loaded OpHits ========================================================
  struct LoadedOpHits {
    
    LoadedOpHits() : ophits() {}
    vector<recob::OpHit> ophits;
    
    void load( vector<recob::OpHit> const & v, int fDebugLevel ) {
      ophits = v;
    }
    
    // not really used 
    bool empty() const { 
      if (ophits.size() == 0) return true;
      return false;
    }

    // Clear waveforms vector.
    void clear(int fDebugLevel) {
      if (fDebugLevel > 3) std::cout << "Clearing LoadedOphits." << std::endl;
      ophits.clear();
      empty();
    }

    std::vector<recob::OpHit> TakeAll() {
      return ophits;
    }
    
    void findinrange(std::vector<recob::OpHit> &obo, 
                     lbne::TpcNanoSlice::Header::nova_timestamp_t first_timestamp,
                     lbne::TpcNanoSlice::Header::nova_timestamp_t last_timestamp,
		     lbne::TpcNanoSlice::Header::nova_timestamp_t event_timestamp,
                     double NovaTicksPerSSPTick,
		     int fDebugLevel) {
      int hh = 0;
      for (auto wf : ophits) { // see if any ophits have pieces inside this TPC boundary
	lbne::TpcNanoSlice::Header::nova_timestamp_t HitTimestamp = wf.PeakTime() * NovaTicksPerSSPTick;
        if ( hh < 5 && fDebugLevel > 2) { // Just want to write out the first few ophits
	  std::cout << "Looking at waveform number " << hh << ". It was on channel " << wf.OpChannel() << " at time " << HitTimestamp
		    << ". The times I passed were " << first_timestamp << " and " << last_timestamp
		    << std::endl;
	}
	if (HitTimestamp <= last_timestamp && HitTimestamp >= first_timestamp) {
	  double NewTime    = wf.PeakTime()    - (event_timestamp/NovaTicksPerSSPTick);
	  double NewAbsTime = wf.PeakTimeAbs() - (event_timestamp/NovaTicksPerSSPTick);
	  recob::OpHit newHit( wf.OpChannel(), NewTime, NewAbsTime, wf.Frame(), wf.Width(), wf.Area(), wf.Amplitude(), wf.PE(), wf.FastToTotal() );
	  obo.emplace_back(std::move(newHit));
	  
          if (fDebugLevel > 3)
	    std::cout << "Pushing back waveform " << hh << " on channel " << wf.OpChannel() << " at corrected time " << NewTime << "("<<HitTimestamp 
		      << "). The times I passed were " << first_timestamp << " and " << last_timestamp << ", event_timestamp " << event_timestamp <<  std::endl;
	}
        ++hh;
      } // auto ophits
      if (fDebugLevel > 1) std::cout << "At the end of Waveform findinrange, wbo has size " << obo.size() << std::endl;
    } // findinrange
    
    //=======================================================================================
    bool OpHitTrigger(lbne::TpcNanoSlice::Header::nova_timestamp_t this_timestamp, double fOpHitADCThreshold, int fOpHitADCsOverThreshold,
                       double fOpHitADCWidth, double NovaTicksPerSSPTick, int fDebugLevel ) { // Triggering on photon detectors
      int hh = 0;
      int HighADCHits = 0;
      for (auto wf : ophits) { // see if any waveforms have pieces inside this TPC boundary
	std::cout << "Looping through ophits, looking at element " << hh << ". It is channel " << wf.OpChannel() << " at time " << wf.PeakTime()
		  << ", it has Amplitude " << wf.Amplitude() << " and " << wf.PE() << " PE's." << std::endl;
	++hh;
	if ( wf.Amplitude() > fOpHitADCThreshold ) ++HighADCHits;
      } // Loop over waveforms
      if ( HighADCHits > fOpHitADCsOverThreshold ) return true;
      else return false;
    } // Photon Trigger
  }; // Loaded OpHits
  
  //=============================== loaded Digits ========================================================
  struct LoadedDigits {
    
    LoadedDigits() : digits(), index(), firstTimestamp(0) {}
    
    vector<RawDigit> digits;
    size_t index;
    lbne::TpcNanoSlice::Header::nova_timestamp_t firstTimestamp; // timestamp of the first nanoslice in the digits vector
    
    void loadTimestamp(lbne::TpcNanoSlice::Header::nova_timestamp_t ts) {
      firstTimestamp = ts;
    }
    
    // Clear digits vector.
    void clear(int fDebugLevel) {
      if (fDebugLevel > 3) std::cout << "Clearing LoadedDigits." << std::endl;
      digits.clear();
      empty(fDebugLevel);
    }

    // copy rawdigits to LoadedDigits and uncompress if necessary.
    void load( vector<RawDigit> const & v, int fDebugLevel ) {
      if (v.size() == 0 || v.back().Compression() == raw::kNone) {
        digits = v;
	if (fDebugLevel > 1) std::cout << "RCE information is not compressed." << std::endl;
      }
      else {
        if (fDebugLevel > 1) std::cout << "RCE information is compressed." << std::endl;
	// make a new raw::RawDigit object for each compressed one
        // to think about optimization -- two copies of the uncompressed raw digits here.
        digits = std::vector<RawDigit>();
        for (auto idigit = v.begin(); idigit != v.end(); ++idigit) {
          std::vector<short> uncompressed;
          int pedestal = (int)idigit->GetPedestal();
          raw::Uncompress(idigit->ADCs(),
                          uncompressed,
                          pedestal,
                          idigit->Compression()
                          );
          raw::RawDigit digit(idigit->Channel(),
                              idigit->Samples(),
                              uncompressed,
                              raw::kNone);
          digits.push_back(digit);
        }
      }
      bool GoodDig = true;
      for (unsigned int j=0; j<digits.size(); j=j+128) {
        if ( digits[0].NADC() - digits[j].NADC() ) GoodDig = false;
        if (fDebugLevel > 2)
	  std::cout << "digits[0] has " << digits[0].NADC() << " ADC's, whilst digits["<<j<<"] has " << digits[j].NADC() << " ADCs. Still a good digit? " << GoodDig << std::endl;
      }
      if (!GoodDig) {
        if (fDebugLevel) std::cout << "Got a bad digit, so want to clear it...\n" << std::endl;
        clear(fDebugLevel);
      }
      
      index = 0ul;
    } // load

    bool empty(int fDebugLevel) const {
      if (digits.size() == 0) {
	if (fDebugLevel > 3) std::cout << "digits.size is 0" << std::endl;
	return true;
      }
      if (digits[0].Samples() == 0) {
	if (fDebugLevel > 3) std::cout << "digits[0].samples is 0" << std::endl;
	return true;
      }
      if (index >= digits[0].Samples()) {
	if (fDebugLevel > 3) std::cout << "The index is more than the number of samples" << std::endl;
	return true;
      }
      return false;
    } // empty
    
    RawDigit::ADCvector_t next() {
      ++index;
      RawDigit::ADCvector_t digitsatindex;
      for (size_t ichan=0; ichan<digits.size(); ichan++) {
        digitsatindex.push_back(digits[ichan].ADC(index-1));
      }
      return digitsatindex;
    }
    
    lbne::TpcNanoSlice::Header::nova_timestamp_t getTimeStampAtIndex(size_t index_local, double novatickspertpctick)
    {
      return firstTimestamp + index_local*novatickspertpctick;
    }
    
  }; // loadedDigits
     //=============================== loaded Counters ========================================================
  struct LoadedCounters {
    
    LoadedCounters() : counters(), index() {}
    
    vector<ExternalTrigger> counters;
    size_t index;
    
    void load( vector<ExternalTrigger> const & v, int fDebugLevel ) {
      
      if (v.size() == 0 )
        counters = v;
      else {
        // make a new raw::ExternalTrigger object for each compressed one
        // to think about optimization -- two copies of the uncompressed raw ExternalTrigger here.
        counters = std::vector<ExternalTrigger>();
        int ii = 0;
        for (auto icounter = v.begin(); icounter != v.end(); ++icounter) {
          std::vector<short> uncompressed;
          raw::ExternalTrigger counter(icounter->GetTrigID(),
                                       icounter->GetTrigTime() );
          counters.push_back(counter);
          ++ii;
        }
      }
      index = 0ul;
    } // load
    
    // not really used 
    bool empty() const { 
      if (counters.size() == 0) return true;
      return false;
    }   

    lbne::TpcNanoSlice::Header::nova_timestamp_t ConvCounterTick(lbne::TpcNanoSlice::Header::nova_timestamp_t TrigTime, double novatickspercounttick) {
      return TrigTime*novatickspercounttick;
    }

    std::vector<ExternalTrigger> TakeAll() {
      return counters;
    }

    void clear(int fDebugLevel) {
      if (fDebugLevel > 3) std::cout << "Clearing LoadedCounters." << std::endl;
      counters.clear();
      empty();
    }

    void findinrange(std::vector<ExternalTrigger> &cbo,
                     lbne::TpcNanoSlice::Header::nova_timestamp_t first_timestamp,
                     lbne::TpcNanoSlice::Header::nova_timestamp_t last_timestamp,
		     lbne::TpcNanoSlice::Header::nova_timestamp_t event_timestamp,
                     unsigned int novatickspercounttick,
		     int fDebugLevel) {
      int hh = 0;
      for (auto count : counters) { // Loop through all the counters
        lbne::TpcNanoSlice::Header::nova_timestamp_t TimeStamp = ConvCounterTick(count.GetTrigTime(), novatickspercounttick);
        if ( hh < 5 && fDebugLevel > 2) { // Just want to look at first few counters.
          std::cout << "Looking at muon counter " << hh << " It has ID " << count.GetTrigID() << " and time " << TimeStamp << "."
		    << " The times I passed were " << first_timestamp << " and " << last_timestamp << std::endl;
        }
        if (TimeStamp <= last_timestamp && TimeStamp >= first_timestamp) {
	  lbne::TpcNanoSlice::Header::nova_timestamp_t NewTime = TimeStamp - (event_timestamp/novatickspercounttick);
	  raw::ExternalTrigger ET(count.GetTrigID(), NewTime );
          cbo.emplace_back(std::move(ET));
	  if (fDebugLevel > 3 )
	    std::cout << "Pushing back counter, channel " << ET.GetTrigID() << "(" << count.GetTrigID() << ") and corrected time " << ET.GetTrigTime() << " (" << count.GetTrigTime() << ")" << std::endl;
 
        } // If within range
	++hh;
      } // auto waveforms
      if (fDebugLevel > 1) std::cout << "At the end of Counter findinrange, cbo has size " << cbo.size() << std::endl;
    } // findinrange

    //=======================================================================================
    bool PTBTrigger(lbne::TpcNanoSlice::Header::nova_timestamp_t this_timestamp, double NovaTicksPerCountTick, double NovaTicksPerTPCTick, int fDebugLevel, std::vector<unsigned int> SpecialChan ) { // Triggering on muon counters
      for ( auto count : counters ) {
        for ( size_t ChanLoop = 0; ChanLoop < SpecialChan.size(); ++ChanLoop) {
	  if ( SpecialChan[ChanLoop] == count.GetTrigID() ) {
	    lbne::TpcNanoSlice::Header::nova_timestamp_t EffecTimeStamp = ConvCounterTick( count.GetTrigTime(), NovaTicksPerCountTick );
	    if ( EffecTimeStamp > this_timestamp
		 && EffecTimeStamp < (this_timestamp+NovaTicksPerTPCTick)
		 ) { // If timestamps match!
	      if (fDebugLevel)
		std::cout << "Triggering on Counter " << count.GetTrigID() << ", " << EffecTimeStamp << ", TPC tick " << this_timestamp << std::endl;
	      return true;
	    } // If timestamps match
	  } // If a special TrigID
	} // Loop through special counters.
      } // Loop through counters
      return false;
    } // CounterTrigger
  
  }; // LoadedCounters

  //=============================== Monte Carlo Stuff ========================================================
  struct MonteCarlo {
    MonteCarlo() : MCParts() {} //, SimChans(), DetSimChans(), SimPhots() {}

    vector<simb::MCParticle> MCParts;
    vector<sim::SimChannel> SimChans;
    //vector<sim::AuxDetSimChannels> DetSimChans;
    //vector<sim::SimPhotonsLits> SimPhots;

    void loadMCParts( vector<simb::MCParticle> const & v ) {
      MCParts = v;
    }

    void loadSimChans( vector<sim::SimChannel> const & b ) {
      SimChans = b;
    }
    
    vector<simb::MCParticle> TakeMCParts(lbne::TpcNanoSlice::Header::nova_timestamp_t start_timestamp,
					 lbne::TpcNanoSlice::Header::nova_timestamp_t end_timestamp,
					 lbne::TpcNanoSlice::Header::nova_timestamp_t event_timestamp,
					 double fNanoSecondsPerNovaTick, int fDebugLevel) {
      vector<simb::MCParticle> retParts;
      if (fDebugLevel > 2) std::cout << "In TakeMCParts....MCParts has size " << MCParts.size() << " I gave it timestamps " << start_timestamp << " " << end_timestamp << std::endl;
      int hh=0;
      for (auto part: MCParts) {
	if ( part.T() > (start_timestamp*fNanoSecondsPerNovaTick) && part.T() < (end_timestamp*fNanoSecondsPerNovaTick) ) {
	  simb::MCParticle newPart = simb::MCParticle(part.TrackId(), part.PdgCode(), part.Process(), part.Mother(), part.Mass(), part.StatusCode());
	  for (size_t qq=0; qq < part.NumberTrajectoryPoints(); ++qq) {
	    const TLorentzVector pos = TLorentzVector( part.Vx(qq), part.Vy(qq), part.Vz(qq), part.T(qq) - (event_timestamp*fNanoSecondsPerNovaTick) );
	    const TLorentzVector mom = TLorentzVector( part.Px(qq), part.Py(qq), part.Pz(qq), part.T(qq) - (event_timestamp*fNanoSecondsPerNovaTick) );
	    newPart.AddTrajectoryPoint( pos, mom );
	  }
	  newPart.SetGvtx( part.Gvx(), part.Gvy(), part.Gvz(), part.Gvt() );
	  newPart.SetEndProcess( part.EndProcess() );
	  newPart.SetPolarization ( part.Polarization() );
	  newPart.SetRescatter( part.Rescatter() );
	  retParts.emplace_back(newPart);
	  if (fDebugLevel > 2) std::cout << "Made a new MCParticle, it has " << newPart.NumberTrajectoryPoints() << " traj points, and time " << newPart.T() << std::endl;
	}
	++hh;
      }
      if (fDebugLevel > 1) std::cout << "At the end of TakeMCParts I am returning a vector of MCParticles with size " << retParts.size() << std::endl;
      return retParts;
    }
    
    vector<sim::SimChannel> TakeSimChans(lbne::TpcNanoSlice::Header::nova_timestamp_t start_timestamp,
					 lbne::TpcNanoSlice::Header::nova_timestamp_t end_timestamp,
					 lbne::TpcNanoSlice::Header::nova_timestamp_t event_timestamp,
					 double fNovaTicksPerTPCTick, int fDebugLevel) {
      vector<sim::SimChannel> retSimChans;
      if (fDebugLevel > 2) std::cout << "In TakeSimChans....SimChans has size " << SimChans.size() << " I gave it timestamps " << start_timestamp << " " << end_timestamp << std::endl;
      int qq = 0;
      for (auto LoopSimChan: SimChans) {
	sim::SimChannel newSimChan = sim::SimChannel( LoopSimChan.Channel() );
	for (std::map<unsigned short,std::vector<sim::IDE> >::const_iterator ideMap = LoopSimChan.TDCIDEMap().begin(); ideMap != LoopSimChan.TDCIDEMap().end(); ++ideMap ) {
	  if ( ideMap->first > ( start_timestamp / fNovaTicksPerTPCTick ) && ideMap->first < ( end_timestamp / fNovaTicksPerTPCTick ) ) {
	    const std::vector<sim::IDE> OldSimIDE = ideMap->second;
	    unsigned short NewTime = ideMap->first - ( event_timestamp / fNovaTicksPerTPCTick );
	    for (size_t zz = 0; zz<OldSimIDE.size(); ++zz) {
	      double IDEPos[3] = { OldSimIDE[zz].x, OldSimIDE[zz].y, OldSimIDE[zz].z };
	      newSimChan.AddIonizationElectrons(OldSimIDE[zz].trackID, NewTime, OldSimIDE[zz].numElectrons, IDEPos, OldSimIDE[zz].energy );
	    }
	    if (fDebugLevel > 2 )
	      std::cout << "The original key value is inbetween my two times..." << start_timestamp <<"("<<start_timestamp / 32<<"), and " << end_timestamp <<"("<<end_timestamp / 32<<")"
			<< ". It has value " << ideMap->first << ", vector has size " << ideMap->second.size()
			<< ". My new map has start time " << NewTime << " and size " << newSimChan.TDCIDEMap().size() << std::endl;
	  }
	}
	retSimChans.emplace_back(newSimChan);
	++qq;
      }
      if (fDebugLevel > 1) std::cout << "At the end of TakeSimChans I am returning a vector of SimChannels with size " << retSimChans.size() << std::endl;
      return retSimChans;
    }
  }; // MonteCarlo
  //===============================================================================================  
  // Retrieves branch name (a la art convention) where object resides
  const char* getBranchName( art::InputTag const & tag, const string inputDataProduct ) {
    std::ostringstream pat_s;
    pat_s << inputDataProduct << "s" 
          << '_'
          << tag.label()
          << '_'
          << tag.instance()
          << '_'
          << tag.process()
          << ".obj";
    
    return pat_s.str().data();
  }
  
  artdaq::Fragments*
  getFragments( TBranch* br, unsigned entry ) {
    br->GetEntry( entry );
    return reinterpret_cast<artdaq::Fragments*>( br->GetAddress() );
  }

  vector<raw::RawDigit>*
  getRawDigits( TBranch* br, unsigned entry ) {
    br->GetEntry( entry );
    return reinterpret_cast<vector<raw::RawDigit>*>( br->GetAddress() );
  }
  
  vector<raw::OpDetWaveform>*
  getSSPWaveforms( TBranch* br, unsigned entry ) {
    br->GetEntry( entry );
    return reinterpret_cast<vector<raw::OpDetWaveform>*>( br->GetAddress() );
  }

  vector<recob::OpHit>*
  getOpHitWaveforms( TBranch* br, unsigned entry ) {
    br->GetEntry( entry );
    return reinterpret_cast<vector<recob::OpHit>*>( br->GetAddress() );
  }
  
  vector<raw::ExternalTrigger>*
  getRawExternalTriggers( TBranch* br, unsigned entry ) {
    br->GetEntry( entry );
    return reinterpret_cast<vector<raw::ExternalTrigger>*>( br->GetAddress() );
  }

  vector<simb::MCParticle>* getMCParticles(TBranch* br, unsigned entry) {
    br->GetEntry( entry );
    return reinterpret_cast<vector<simb::MCParticle>*>(br->GetAddress() );
  }

  vector<sim::SimChannel>* getMCSimChans(TBranch* br, unsigned entry) {
    br->GetEntry( entry );
    return reinterpret_cast<vector<sim::SimChannel>*>(br->GetAddress() );
  }
}

//==========================================================================
namespace raw {
  // Enable 'pset.get<raw::Compress_t>("compression")'
  void decode (boost::any const & a, Compress_t & result ){
    unsigned tmp;
    fhicl::detail::decode(a,tmp);
    result = static_cast<Compress_t>(tmp);
  }
}
//==========================================================================

namespace DAQToOffline {
  // The class Splitter is to be used as the template parameter for art::Source<T>. It understands how to read art's ROOT data files,
  // to extract artdaq::Fragments from them, how to convert the artdaq::Fragments into vectors of raw::RawDigit objects, and
  // finally how to re-combine those raw::RawDigit objects, raw::OpDetWaveform objects,  and external trigger objects to present
  // the user of the art::Source<Splitter> with a sequence of art::Events that have the desired event structure, different from
  // the event structure of the data file(s) being read.
  
  class Splitter {
  public:
    Splitter(fhicl::ParameterSet const& ps,
             art::ProductRegistryHelper& prh,
             art::SourceHelper& sh);
    
    // See art/Framework/IO/Sources/Source.h for a description of each of the public member functions of Splitter.
    bool readFile(string const& filename, art::FileBlock*& fb);
    
    bool readNext(art::RunPrincipal* const& inR,
                  art::SubRunPrincipal* const& inSR,
                  art::RunPrincipal*& outR,
                  art::SubRunPrincipal*& outSR,
                  art::EventPrincipal*& outE);
    
    void closeCurrentFile();
    
  private:
    
    using rawDigits_t = vector<RawDigit>;
    using SSPWaveforms_t = vector<OpDetWaveform>;
    using OpHits_t = vector<recob::OpHit>;
    using PennCounters_t = vector<ExternalTrigger>;

    using MCPart_t = vector<simb::MCParticle>;
    using MCSimChan_t = vector<sim::SimChannel>;
    
    string                 sourceName_;
    string                 lastFileName_;
    std::unique_ptr<TFile> file_;
    bool                   doneWithFiles_;
    art::InputTag          TPCinputTag_;
    art::InputTag          SSPinputTag_;
    art::InputTag          OpHitinputTag_;
    art::InputTag          PenninputTag_;
    art::InputTag          MCPartinputTag_;
    art::InputTag          MCSimChaninputTag_;
    string                 TPCinputDataProduct_;
    string                 SSPinputDataProduct_;
    string                 OpHitinputDataProduct_;
    string                 PenninputDataProduct_;
    string                 MCPartinputDataProduct_;
    string                 MCSimChaninputDataProduct_;
    SSPReformatterAlgs     sspReform;
    string                 fTPCChannelMapFile;
    string                 fPTBMapFile;
    string                 fPTBMapDir;
    art::SourceHelper      sh_;
    TBranch*               TPCinputBranch_;
    TBranch*               SSPinputBranch_;
    TBranch*               OpHitinputBranch_;
    TBranch*               PenninputBranch_;
    TBranch*               MCPartinputBranch_;
    TBranch*               MCSimChaninputBranch_;
    TBranch*               EventAuxBranch_;
    LoadedDigits           loadedDigits_;
    LoadedWaveforms        loadedWaveforms_;
    LoadedOpHits           loadedOpHits_;
    LoadedCounters         loadedCounters_;
    MonteCarlo             MonteCarlo_;
    size_t                 nInputEvts_;
    size_t                 treeIndex_;
    art::RunNumber_t       runNumber_;
    art::SubRunNumber_t    subRunNumber_;
    art::EventNumber_t     eventNumber_;
    art::RunNumber_t       cachedRunNumber_;
    art::SubRunNumber_t    cachedSubRunNumber_;
    art::RunNumber_t       inputRunNumber_;
    art::SubRunNumber_t    inputSubRunNumber_;
    art::EventNumber_t     inputEventNumber_;
    art::Timestamp         inputEventTime_; 
    size_t                 ticksPerEvent_;
    rawDigits_t            bufferedDigits_;
    std::vector<RawDigit::ADCvector_t>  dbuf_;
    SSPWaveforms_t         wbuf_;
    OpHits_t               hbuf_;
    PennCounters_t         cbuf_;
    MCPart_t               mcbuf_;
    MCSimChan_t            simchanbuf_;
    unsigned short         fTicksAccumulated;

    bool                   fTrigger = false; 
    size_t                 fLastTriggerIndex = 0;
    size_t                 fLastTreeIndex = 0;
    int                    fDiffFromLastTrig = 0;
    lbne::TpcNanoSlice::Header::nova_timestamp_t fLastTimeStamp = 0;
    lbne::TpcNanoSlice::Header::nova_timestamp_t first_timestamp=0;
    lbne::TpcNanoSlice::Header::nova_timestamp_t Event_timestamp=0;
    lbne::TpcNanoSlice::Header::nova_timestamp_t last_timestamp=0;
    lbne::TpcNanoSlice::Header::nova_timestamp_t this_timestamp=0;
    lbne::TpcNanoSlice::Header::nova_timestamp_t prev_timestamp=0;

    std::map<int,int>      TPCChannelMap;

    std::map<uint64_t,size_t> EventTreeMap;

    std::function<rawDigits_t(artdaq::Fragments const&, lbne::TpcNanoSlice::Header::nova_timestamp_t&, std::map<int,int> const&)> fragmentsToDigits_;

    bool eventIsFull_(rawDigits_t const & v);

    bool loadEvents_( size_t &InputTree );
    bool LoadPTBInformation( size_t LoadTree );
    void LoadSSPInformation( size_t LoadTree );
    void LoadOpHitInformation( size_t LoadTree );
    void LoadRCEInformation( size_t LoadTree );

    void makeEventAndPutDigits_( art::EventPrincipal*& outE, art::Timestamp art_timestamp=0);

    void Reset();

    void CheckTimestamps(bool &JumpEvent, size_t &JumpNADC );

    bool NoRCEsCase(art::RunPrincipal*& outR, art::SubRunPrincipal*& outSR, art::EventPrincipal*& outE);

    void Triggering(std::map<int,int> &PrevChanADC, std::vector<short> ADCdigits, bool NewTree);
    
    void CheckTrigger();
    
    bool TicklerTrigger( std::map<int,int> &PrevChanADC, std::vector<short> ADCdigits );

    art::EventAuxiliary    evAux_;
    art::EventAuxiliary*   pevaux_;

    bool         fRequireRCE;
    bool         fRequireSSP;
    bool         fRequireOpHit;
    bool         fRequirePTB;
    bool         fMonteCarlo;
    size_t       fPostTriggerTicks;
    size_t       fPreTriggerTicks;
    double       fNovaTicksPerTPCTick;
    double       fNovaTicksPerSSPTick;
    double       fNovaTicksPerCountTick;
    double       fNanoSecondsPerNovaTick;
    int          fDebugLevel;
    double       fTimeStampThreshold;
    int          fMCTrigLevel;
    std::vector<unsigned int> fWhichTrigger;
    std::vector<unsigned int> fPTBTrigs;
    int          fTrigSeparation;
    double       fWaveformADCWidth;
    double       fWaveformADCThreshold;
    int          fWaveformADCsOverThreshold;
    double       fOpHitADCWidth;
    double       fOpHitADCThreshold;
    int          fOpHitADCsOverThreshold;
    int          fADCdiffThreshold;
    int          fADCsOverThreshold;
    bool         fUsePedestalDefault;
    bool         fUsePedestalFile;
    bool         fUsePedestalFileSearchPath;
    std::string  fPedestalFile;
    std::string  fPedestalFileSearchPath;
    int          fSkipNInputEvents;
    int          fSkipNOutputEvents;

    bool RCEsNotPresent = false;
    int  gSkippedOuputEvents = 0;
    std::vector<std::pair<size_t,size_t> > GoodEvents;

    std::pair <std::pair<lbne::PennMicroSlice::Payload_Header::short_nova_timestamp_t, std::bitset<TypeSizes::CounterWordSize> >,
               std::pair<lbne::PennMicroSlice::Payload_Header::short_nova_timestamp_t, std::bitset<TypeSizes::TriggerWordSize> > > PrevTimeStampWords;

    std::map<uint16_t, std::map <size_t, std::pair<float,float> > > AllPedMap;
    std::map<int,int> fPTBMap;
  };
}

//=======================================================================================
DAQToOffline::Splitter::Splitter(fhicl::ParameterSet const& ps,
                                 art::ProductRegistryHelper& prh,
                                 art::SourceHelper& sh) :
  sourceName_("SplitterInput"),
  lastFileName_(ps.get<vector<string>>("fileNames",{}).back()),
  file_(),
  doneWithFiles_(false),
  //  TPCinputTag_("daq:TPC:DAQ"), // "moduleLabel:instance:processName"
  TPCinputTag_         (ps.get<string>("TPCInputTag")),
  SSPinputTag_         (ps.get<string>("SSPInputTag")),
  OpHitinputTag_       (ps.get<string>("OpHitInputTag")),
  PenninputTag_        (ps.get<string>("PennInputTag")),
  MCPartinputTag_      (ps.get<string>("MCPartInputTag")),
  MCSimChaninputTag_   (ps.get<string>("MCSimChanInputTag")),
  TPCinputDataProduct_ (ps.get<string>("TPCInputDataProduct")),
  SSPinputDataProduct_ (ps.get<string>("SSPInputDataProduct")),
  OpHitinputDataProduct_(ps.get<string>("OpHitInputDataProduct")),
  PenninputDataProduct_(ps.get<string>("PennInputDataProduct")),
  MCPartinputDataProduct_(ps.get<string>("MCPartInputDataProduct")),
  MCSimChaninputDataProduct_(ps.get<string>("MCSimChanInputDataProduct")),
  sspReform            (ps.get<fhicl::ParameterSet>("SSPReformatter")),
  fTPCChannelMapFile   (ps.get<string>("TPCChannelMapFile")),
  fPTBMapFile          (ps.get<std::string>("PTBMapFile")),
  fPTBMapDir           (ps.get<std::string>("PTBMapDir")),
  sh_(sh),
  TPCinputBranch_(nullptr),
  SSPinputBranch_(nullptr),
  EventAuxBranch_(nullptr),
  nInputEvts_(),
  runNumber_(1),      // Defaults 
  subRunNumber_(0),   
  eventNumber_(),
  cachedRunNumber_(-1),
  cachedSubRunNumber_(-1),
  inputRunNumber_(1),      
  inputSubRunNumber_(1),   
  inputEventNumber_(1),
  bufferedDigits_(),
  dbuf_(),
  wbuf_(),
  hbuf_(),
  cbuf_(),
  fTicksAccumulated(0),
  fragmentsToDigits_( std::bind( DAQToOffline::tpcFragmentToRawDigits,
                                 std::placeholders::_1, // artdaq::Fragments
                                 std::placeholders::_2, // lbne::TpcNanoSlice::Header::nova_timestamp_t& firstTimestamp
                                 std::placeholders::_3, // the channel map
                                 ps.get<bool>("debug",false),
                                 ps.get<raw::Compress_t>("compression",raw::kNone),
                                 ps.get<unsigned>("zeroThreshold",0) ) ),
  fRequireRCE            (ps.get<bool>  ("RequireRCE")),
  fRequireSSP            (ps.get<bool>  ("RequireSSP")),
  fRequireOpHit          (ps.get<bool>  ("RequireOpHit")),
  fRequirePTB            (ps.get<bool>  ("RequirePTB")),
  fMonteCarlo            (ps.get<bool>  ("MonteCarlo")),
  fPostTriggerTicks      (ps.get<size_t>("PostTriggerTicks")),
  fPreTriggerTicks       (ps.get<size_t>("PreTriggerTicks")),
  fNovaTicksPerTPCTick   (ps.get<double>("NovaTicksPerTPCTick")),
  fNovaTicksPerSSPTick   (ps.get<double>("NovaTicksPerSSPTick")),
  fNovaTicksPerCountTick (ps.get<double>("NovaTicksPerCountTick")),
  fNanoSecondsPerNovaTick(ps.get<double>("NanoSecondsPerNovaTick")),
  fDebugLevel            (ps.get<int>   ("DebugLevel")),
  fTimeStampThreshold    (ps.get<double>("TimeStampThreshold")),
  fMCTrigLevel           (ps.get<int>   ("MCTrigLevel")),
  fWhichTrigger          (ps.get<std::vector<unsigned int> >("WhichTrigger")),
  fPTBTrigs              (ps.get<std::vector<unsigned int> >("PTBTrigs")),
  fTrigSeparation        (ps.get<int>   ("TrigSeparation")),
  fWaveformADCWidth      (ps.get<double>("WaveformADCWidth")),
  fWaveformADCThreshold  (ps.get<double>("WaveformADCThreshold")),
  fWaveformADCsOverThreshold(ps.get<double>("WaveformADCsOverThreshold")),
  fOpHitADCWidth         (ps.get<double>("OpHitADCWidth")),
  fOpHitADCThreshold     (ps.get<double>("OpHitADCThreshold")),
  fOpHitADCsOverThreshold(ps.get<double>("OpHitADCsOverThreshold")),
  fADCdiffThreshold      (ps.get<int>   ("ADCdiffThreshold")),
  fADCsOverThreshold     (ps.get<int>   ("ADCsOverThreshold")),
  fUsePedestalDefault    (ps.get<bool>  ("UsePedestalDefault")),
  fUsePedestalFile       (ps.get<bool>  ("UsePedestalFile")),
  fUsePedestalFileSearchPath(ps.get<bool>    ("UsePedestalFileSearchPath",false)),
  fPedestalFile          (ps.get<std::string>("PedestalFile")),
  fPedestalFileSearchPath(ps.get<std::string>("PedestalFileSearchPath", "")),
  fSkipNInputEvents      (ps.get<int>   ("SkipNInputEvents")),
  fSkipNOutputEvents     (ps.get<int>   ("SkipNOutputEvents"))
{
  ticksPerEvent_ = fPostTriggerTicks + fPreTriggerTicks;
  // Will use same instance names for the outgoing products as for the
  // incoming ones.
  prh.reconstitutes<rawDigits_t,art::InEvent>( sourceName_, TPCinputTag_.instance() );
  prh.reconstitutes<SSPWaveforms_t,art::InEvent>( sourceName_, SSPinputTag_.instance() );
  prh.reconstitutes<OpHits_t,art::InEvent>( sourceName_, OpHitinputTag_.instance() );
  prh.reconstitutes<PennCounters_t,art::InEvent>( sourceName_, PenninputTag_.instance() );
  // If looking at Monte Carlo, also want to copy the truth information to the split event.
  if (fMonteCarlo) {
    prh.reconstitutes<MCPart_t,art::InEvent>( sourceName_, MCPartinputTag_.instance() );
    prh.reconstitutes<MCSimChan_t,art::InEvent>( sourceName_, MCSimChaninputTag_.instance() );
  }

  BuildTPCChannelMap(fTPCChannelMapFile, TPCChannelMap);
  //std::cout << "Built TPC Channel Map" << std::endl;
  BuildPTBChannelMap(fPTBMapDir, fPTBMapFile, fPTBMap);
}

//=======================================================================================
bool DAQToOffline::Splitter::readFile(string const& filename, art::FileBlock*& fb) {
  //std::cout << "At the top of readFile" << std::endl;
  // Get fragments branches
  file_.reset( new TFile(filename.data()) );
  TTree* evtree    = reinterpret_cast<TTree*>(file_->Get(art::rootNames::eventTreeName().c_str()));
  
  TPCinputBranch_ = evtree->GetBranch( getBranchName(TPCinputTag_, TPCinputDataProduct_ ) ); // get branch for TPC input tag
  SSPinputBranch_ = evtree->GetBranch( getBranchName(SSPinputTag_, SSPinputDataProduct_ ) ); // get branch for SSP input tag
  OpHitinputBranch_ = evtree->GetBranch( getBranchName(OpHitinputTag_, OpHitinputDataProduct_ ) ); // get branch for OpHit input tag
  PenninputBranch_ = evtree->GetBranch( getBranchName(PenninputTag_, PenninputDataProduct_ ) ); // get branch for Penn Board input tag

  MCPartinputBranch_ = evtree->GetBranch( getBranchName(MCPartinputTag_, MCPartinputDataProduct_ ) );
  MCSimChaninputBranch_ = evtree->GetBranch( getBranchName(MCSimChaninputTag_, MCSimChaninputDataProduct_ ) );

  if (TPCinputBranch_) nInputEvts_      = static_cast<size_t>( TPCinputBranch_->GetEntries() );
  size_t nevt_ssp  = 0;
  if (SSPinputBranch_) nevt_ssp = static_cast<size_t>( SSPinputBranch_->GetEntries() );
  size_t nevt_ophit = 0;
  if (OpHitinputBranch_) nevt_ophit = static_cast<size_t>( OpHitinputBranch_->GetEntries() );
  size_t nevt_penn  = 0;
  if (PenninputBranch_) nevt_penn  = static_cast<size_t>( PenninputBranch_->GetEntries());
  
  if (nevt_ssp != nInputEvts_&& nevt_ssp) throw cet::exception("35-ton SplitterInput: Different numbers of RCE and SSP input events in file");
  if (nevt_ophit !=  nInputEvts_&& nevt_ophit) throw cet::exception("35-ton SplitterInput: Different numbers of RCE and OpHit input events in file");
  if (nevt_penn != nInputEvts_&& nevt_penn) throw cet::exception("35-ton SplitterInput: Different numbers of RCE and Penn input events in file");
  treeIndex_ = 0ul;
  
  EventAuxBranch_ = evtree->GetBranch( "EventAuxiliary" );
  pevaux_ = &evAux_;
  EventAuxBranch_->SetAddress(&pevaux_);

  // ------------ Make my sorted event branch --------------------
  EventTreeMap.clear();
  art::RunNumber_t TreeRunNumber; uint16_t IntRunNumber;
  art::SubRunNumber_t TreeSubRunNumber; uint16_t IntSubRunNumber;
  art::EventNumber_t TreeEventNumber; uint32_t IntEventNumber;
  uint64_t CombinedInt;

  std::cout << "\nfDebugLevel is " << fDebugLevel << ", PedestalDefault is " << fUsePedestalDefault << " fUsePedestalFile is " << fUsePedestalFile << std::endl;

  art::RunNumber_t PrevRunNumber = -1;
  for (size_t Tree=1; Tree < nInputEvts_; ++Tree) {
    EventAuxBranch_->GetEntry(Tree); TreeRunNumber = evAux_.run(); //Get the run number
    IntRunNumber    = (int)TreeRunNumber; TreeSubRunNumber = evAux_.subRun(); IntSubRunNumber = (int)TreeSubRunNumber; //Get the subrun number
    TreeEventNumber = evAux_.event();   IntEventNumber  = (int)TreeEventNumber; //Get the event number
    CombinedInt = (uint64_t) IntRunNumber << 16 | IntSubRunNumber << 16 | IntEventNumber; // Combine them as a 64 bit int.
    if (fDebugLevel > 4 )
      std::cout << "Looking at Tree " << Tree << ", RunNumber " << IntRunNumber << ", SubRunNumber " << IntSubRunNumber << ", EventNumber " << IntEventNumber << ", CrazyNumber " << CombinedInt << std::endl;
    EventTreeMap[CombinedInt] = Tree; // Add that to a tree - use the fact that this will sort them by Run, Subrun, Event.
  
    art::RunNumber_t ThisNumber = evAux_.run();;
    if (ThisNumber != PrevRunNumber ) {
      if ( evAux_.isRealData() ) { // If real data subtract pedestal conditions
	dune::DetPedestalDUNE pedestals("dune35t");
	pedestals.SetDetName("dune35t");
	pedestals.SetUseDefaults(fUsePedestalDefault);
	if ( fUsePedestalFile ) {
          std::string fullname;
          if( fUsePedestalFileSearchPath ){
            if( fPedestalFileSearchPath != "" ){
              std::cout << "SplitterStitcher: fPedestalFileSearchPath: " << fPedestalFileSearchPath << " fPedestalFile: " << fPedestalFile << std::endl;
              cet::search_path sp(fPedestalFileSearchPath);
              if(fPedestalFile != "" ) sp.find_file(fPedestalFile, fullname);
              else fullname = cet::getenv(fPedestalFileSearchPath);
            }//fPedestalFileSearchPath != ""
            else{
              std::cout << "SplitterStitcher: fPedestalFileSearchPath: " << fPedestalFileSearchPath << " fPedestalFile: " << fPedestalFile << std::endl;
              fullname = fPedestalFile;
            }
          }//fUsePedestalFileSearchPath
          else fullname = fPedestalFile;

	  std::cout << "SplitterStitcher: Setting CSVFileName to " << fullname << std::endl;
	  pedestals.SetCSVFileName(fullname);
	} else {
	  pedestals.SetUseDB(true);
	}
	pedestals.Update(ThisNumber);
	for (size_t ichan=0;ichan<2048;ichan++) {
	  AllPedMap[ThisNumber][ichan].first  = pedestals.PedMean(ichan);
	  AllPedMap[ThisNumber][ichan].second = pedestals.PedMeanErr(ichan);
	  if (fDebugLevel > 2) {
	    std::cout << "AllPedMap["<<ThisNumber<<"]["<<ichan<<"] has mean " << pedestals.PedMean(ichan) << " (" << AllPedMap[ThisNumber][ichan].first << ")"
		      << "and error " << pedestals.PedMeanErr(ichan) << " (" << AllPedMap[ThisNumber][ichan].second << "). " << std::endl;
	  }
	}
	PrevRunNumber = ThisNumber;
      }
    }
  }

  // New fileblock
  fb = new art::FileBlock(art::FileFormatVersion(),filename);
  if ( fb == nullptr ) {
    throw art::Exception(art::errors::FileOpenError)
      << "Unable to open file " << filename << ".\n";
  }

  //std::cout << "At the end of ReadFile" << std::endl;

  return true;
}

//=======================================================================================
bool DAQToOffline::Splitter::readNext(art::RunPrincipal*    const& inR,
                                 art::SubRunPrincipal* const& inSR,
                                 art::RunPrincipal*    & outR,
                                 art::SubRunPrincipal* & outSR,
                                 art::EventPrincipal*  & outE) {
  //std::cout << "At the start of readNext..." << std::endl;
  if (!fRequireRCE) {
    std::cout << "Entering NoRCEsCase" << std::endl;
    bool Return = NoRCEsCase(outR, outSR, outE);
    std::cout << "Left NoRCEsCase" << std::endl;
    return Return;
  }

  if ( doneWithFiles_ ) {
    return false;
  }
  first_timestamp = Event_timestamp = last_timestamp = this_timestamp = prev_timestamp = 0;
  double FirstDigIndex;
  size_t FirstDigTree;
  bool   first_tick = true; // The earliest time in this new split event, so want to calculate time of this for use with first_timestamp variable only!
  bool   NewTree;
  bool   JumpEvent = false;
  size_t JumpNADC  = 0;

  std::map<int,int> PrevChanADC;
  if (fDebugLevel > 3 ) {
    std::cout << "\nAt the top of readNext....what do I increment here? " << fTicksAccumulated << " " << ticksPerEvent_ << " " << loadedDigits_.empty(fDebugLevel)
	      << " " << wbuf_.size() << " " << cbuf_.size() << " " << hbuf_.size() << std::endl;
  }
  while ( fTicksAccumulated < ticksPerEvent_ ) {  
    ++fDiffFromLastTrig;
    NewTree = false;
    prev_timestamp = this_timestamp; // set prev_timestamp to old timestamp

    // ************* Check if loadedDigits is empty... ************************
    while (loadedDigits_.empty(fDebugLevel)) {
      if (fDebugLevel > 3 ) std::cout << "\nLoaded digits is empty..." << std::endl;
      if ( fTrigger ) { // Want to load wbuf with end of last event, before loading new data.
        loadedWaveforms_.findinrange(wbuf_, first_timestamp, last_timestamp, Event_timestamp, fNovaTicksPerSSPTick, fDebugLevel);
	loadedOpHits_.findinrange   (hbuf_, first_timestamp, last_timestamp, Event_timestamp, fNovaTicksPerSSPTick, fDebugLevel);
        loadedCounters_.findinrange (cbuf_, first_timestamp, last_timestamp, Event_timestamp, fNovaTicksPerCountTick, fDebugLevel);
	if (fDebugLevel > 2 ) {
	  std::cout << "Loaded digits was empty, will be refilled..."
		    << "\nwbuf_ has size " << wbuf_.size() << " at " << first_timestamp << " " << last_timestamp << " " << fNovaTicksPerSSPTick
		    << "\nhbuf_ has size " << hbuf_.size() << " at " << first_timestamp << " " << last_timestamp << " " << fNovaTicksPerSSPTick
		    << "\ncbuf_ has size " << cbuf_.size() << " at " << first_timestamp << " " << last_timestamp << " " << fNovaTicksPerCountTick
		    << std::endl;
	}
      }
      bool rc = loadEvents_(treeIndex_);
      if (RCEsNotPresent) {
	if (fDebugLevel) std::cout << "The RCEs aren't present, so switching to the don't require RCEs case...." << std::endl;
	fRequireRCE = false;
	bool Return = NoRCEsCase(outR, outSR, outE);
	return Return;
      } // RCEsNotPresent
      if (fDebugLevel > 2) std::cout << "There are a total of " << loadedDigits_.digits[0].NADC() << " ADC's " << std::endl;
      if (!rc) {
        doneWithFiles_ = (file_->GetName() == lastFileName_);
        return false;
      }
      NewTree = true;
      if (treeIndex_ == fLastTreeIndex) 
        loadedDigits_.index = fLastTriggerIndex; // If have to re-load an old tree, need to jump back to previous position in that tree.
      first_tick = true; // Just loaded in new event, so want to reset first_timestamp...doesn't effect previous loaded event.
      first_timestamp=0, last_timestamp=0; // Want to reset the timestamps.
      // ******* Check that the time stamps lead on from one another!! ***********
      if ( loadedDigits_.digits.size() != 0 && loadedDigits_.digits[0].NADC() ) {
        bool NewEvent = true;
        for (unsigned int GoodEvSize = 0; GoodEvSize < GoodEvents.size(); ++GoodEvSize ) {
          if (GoodEvents[GoodEvSize].first == treeIndex_-1 ) NewEvent = false;
        }
        if (NewEvent) {
          if (fDebugLevel > 2) std::cout << "Adding a new event to goodEvents" << std::endl;
          GoodEvents.push_back( std::make_pair(treeIndex_-1,loadedDigits_.digits[0].NADC()) );
        }
        if (fTrigger ) {
	  CheckTimestamps( JumpEvent, JumpNADC ); // Check that the time stamps lead on from one another!!
	}
      }
    } // loadedDigits_.empty()
    
    if (NewTree) {
      if (JumpEvent) {
        loadedDigits_.index = fPreTriggerTicks - JumpNADC;
        JumpEvent = false;
        JumpNADC  = 0;
      }
      if (fDebugLevel) {
	std::cout << "\nLooking at treeIndex " << treeIndex_-1 << ", index " << loadedDigits_.index 
	  //<< ". I have missed " << fDiffFromLastTrig << " ticks since my last trigger at treeIndex " << fLastTreeIndex << ", tick " << fLastTriggerIndex 
		  << std::endl;
      }
    }
    
    std::vector<short> nextdigits = loadedDigits_.next();
    this_timestamp = loadedDigits_.getTimeStampAtIndex(loadedDigits_.index, fNovaTicksPerTPCTick);
    if (fTrigger && fDebugLevel > 4) {
      std::cout << "Index " << loadedDigits_.index << " " << fTicksAccumulated << " " << prev_timestamp << " " << this_timestamp << std::endl;
    }    
    // ******* See if can trigger on this tick...only want to do this if haven't already triggered. We also don't want two triggers too close together!... *****************
    if ( fTicksAccumulated == 0 && fDiffFromLastTrig >= fTrigSeparation ) {
      Triggering (PrevChanADC, nextdigits, NewTree);
    } // if TickAccumulated == 0
    // ******* See if can trigger on this tick...only want to do this if haven't already triggered.... *****************
    if ( fTrigger ) { // Have triggered!
      if (fTicksAccumulated == 0 ) {
        FirstDigIndex = loadedDigits_.index;
        FirstDigTree  = treeIndex_ - 1;
        Event_timestamp = this_timestamp;
        if (fDebugLevel) {
	  std::cout << "\nThe trigger is good so triggering on, treeIndex " << fLastTreeIndex
		    << ", loadedDigits_.index() " << fLastTriggerIndex << ", with timestamp " << fLastTimeStamp
		    << "\nThe first tick in this event is in tree index " << inputEventNumber_ << ", loadedDigits index " << loadedDigits_.index
		    << ". It has timestamp " << this_timestamp << "\n"
		    << std::endl;
	}
      }
      // ************* Work out first and last SSP Timestamp for wbuf_, hbuf_ and cbuf_ ************************
      if (first_tick) { // First tick in split event and/or first tick in newly loaded event.
        first_timestamp = this_timestamp;
	first_tick = false;
      }
      last_timestamp = this_timestamp;
      
      // ************* Now want to load the RCE information into dbuf_ ************************
      if (dbuf_.size() == 0) {
        RawDigit::ADCvector_t emptyvector;
        for (size_t ichan=0;ichan<nextdigits.size();ichan++) dbuf_.push_back(emptyvector);
      }
      for (size_t ichan=0;ichan<nextdigits.size();ichan++) {
        dbuf_[ichan].push_back(nextdigits[ichan]);
      }
      fTicksAccumulated ++;  
    } // If triggered on this tick!
    // ************* Now Done for this tick ************************
    // If I want to skip N Output events, I want to delete the event now that I have made it...
    if (fTicksAccumulated == ticksPerEvent_ && gSkippedOuputEvents < fSkipNOutputEvents) {
      if (fDebugLevel) std::cout << "I have skipped " << gSkippedOuputEvents << " events of a desired " << fSkipNOutputEvents << ". I am now voiding another...\n" << std::endl;
      ++gSkippedOuputEvents;
      Reset();
    }
  } // while ticks accumulated < ticksperEvent
  
  if (fDebugLevel > 1)
    std::cout << "Got enough ticks now start building the events...Start/End time for Waveforms/Counters are " << first_timestamp << " " << last_timestamp <<  std::endl;
  // ************* Fill wbuf_ with the SSP information within time range ************************
  if (fDebugLevel > 1)
    std::cout << "Loading the Waveforms...wbuf_ has size " << wbuf_.size() << " " << fNovaTicksPerSSPTick << std::endl;
  loadedWaveforms_.findinrange(wbuf_, first_timestamp, last_timestamp, Event_timestamp, fNovaTicksPerSSPTick, fDebugLevel);
  if (fDebugLevel > 1)
    std::cout << "wbuf_ now has size " << wbuf_.size() << std::endl;

  // ************* Fill hbuf_ with the OpHit information within time range ************************
  if (fDebugLevel > 1)
    std::cout << "Loading the Waveforms...hbuf_ has size " << hbuf_.size() << " " << fNovaTicksPerSSPTick << std::endl;
  loadedOpHits_.findinrange(hbuf_, first_timestamp, last_timestamp, Event_timestamp, fNovaTicksPerSSPTick, fDebugLevel);
  if (fDebugLevel > 1)
    std::cout << "hbuf_ now has size " << hbuf_.size() << std::endl;
  
  // ************* Fill cbuf_ with the PTB information within time range ************************
  if (fDebugLevel > 1)
    std::cout << "Loading the Counters! cbuf size " << cbuf_.size() << " counterticks " << fNovaTicksPerCountTick << std::endl;
  loadedCounters_.findinrange(cbuf_, first_timestamp, last_timestamp, Event_timestamp, fNovaTicksPerCountTick, fDebugLevel);
  if (fDebugLevel > 1)
    std::cout << "Now cbuf has size " << cbuf_.size() << std::endl;

  // ******** Now Build the event *********
  runNumber_ = inputRunNumber_;
  subRunNumber_ = inputSubRunNumber_;
  //art::Timestamp ts; // LBNE should decide how to initialize this -- use first_timestamp converted into an art::Timestamp
  //FIXME - This is a first attempt at interpreting the novatimestamp from the tpc data to create an art event timestamp
  art::Timestamp this_art_event_timestamp = DAQToOffline::make_art_timestamp_from_nova_timestamp(Event_timestamp);
  if ( runNumber_ != cachedRunNumber_ ) {
    outR = sh_.makeRunPrincipal(runNumber_,this_art_event_timestamp);
    cachedRunNumber_ = runNumber_;
    eventNumber_ = 0ul;
  }
  if ( subRunNumber_ != cachedSubRunNumber_ ) {
    outSR = sh_.makeSubRunPrincipal(runNumber_,subRunNumber_,this_art_event_timestamp);
    cachedSubRunNumber_ = subRunNumber_;
    eventNumber_ = 0ul;
  }

  // ************* Now fill dbuf_ with TPC information collected ************************
  if (fDebugLevel > 1) std::cout << "Just about to fill d " << fTicksAccumulated << " " << dbuf_.size() << std::endl;
  for (size_t ichan=0;ichan<dbuf_.size();ichan++) {
    RawDigit d(loadedDigits_.digits[ichan].Channel(),
               fTicksAccumulated,
               dbuf_[ichan]
               //,loadedDigits_.digits[ichan].Compression()
               );
    if (evAux_.isRealData() ) //If looking at real data!
      d.SetPedestal(AllPedMap[runNumber_][loadedDigits_.digits[ichan].Channel()].first,
                    AllPedMap[runNumber_][loadedDigits_.digits[ichan].Channel()].second );
    else //If looking at Truth
      d.SetPedestal(loadedDigits_.digits[ichan].GetPedestal(),
                    loadedDigits_.digits[ichan].GetSigma() );
    if (fDebugLevel > 3) {
      std::cout << "digit[0] corresponding to channel " << d.Channel() << " ("<<ichan<<") has ADC value " << d.ADC(0)
		<< ", pedestal "<<d.GetPedestal()<<" ["<<AllPedMap[runNumber_][loadedDigits_.digits[ichan].Channel()].first <<"],"
		<< " and sigma "<<d.GetSigma()   <<" ["<<AllPedMap[runNumber_][loadedDigits_.digits[ichan].Channel()].second<<"]."
		<< std::endl;
    }
    bufferedDigits_.emplace_back(d);
  }

  // If looking at Monte Carlo, now want to add the Truth stuff
  if (fMonteCarlo) {
    mcbuf_ = MonteCarlo_.TakeMCParts   (first_timestamp, last_timestamp, Event_timestamp, fNanoSecondsPerNovaTick, fDebugLevel);
    if (fDebugLevel) std::cout << "Have now returned from TakeMCParts, it has size " << mcbuf_.size() << std::endl;
    simchanbuf_ = MonteCarlo_.TakeSimChans(first_timestamp, last_timestamp, Event_timestamp, fNovaTicksPerTPCTick, fDebugLevel);
    if (fDebugLevel) std::cout << "Have now returned from TakeMCSimChans, it has size " << simchanbuf_.size() << std::endl;
  } 


  //inputEventTime_ is the art::Timestamp() of the online art::Event() used to create the offline art::Event()
  makeEventAndPutDigits_( outE, inputEventTime_ );
  
  // ******** Reset loadedDigits_.index and TreeIndex_ to where the trigger was *********
  if (fDebugLevel) {
    std::cout << "\nMaking an event which triggered on Tree Index " << fLastTreeIndex << ", tick " << fLastTriggerIndex << ".\n"
	      << "It went from Tree index " << FirstDigTree  << ", tick " << FirstDigIndex << " to Tree index " << treeIndex_-1 << ", tick " << loadedDigits_.index << ".\n" 
	      << "I want to reset the tick value for sure, but do I need to reload the digits because treeIndex is different?" << std::endl;
  }
  if ( treeIndex_-1 != fLastTreeIndex ) {
    if (fDebugLevel) std::cout << "Yes I do! Changing treeIndex_ to fLastTreeIndex...Also want to clear loadedDigits." << std::endl;
    treeIndex_ = fLastTreeIndex;
    loadedDigits_.index = 0;
    loadedDigits_.clear(fDebugLevel);
    loadEvents_(treeIndex_);
  } else {
    if (fDebugLevel) std::cout << "No, I'm still looking at the same tree!\n" << std::endl;
  }
  loadedDigits_.index = fLastTriggerIndex;
  this_timestamp      = fLastTimeStamp;
  
  if (fDebugLevel) {
    std::cout << "This is a list of the good events seen so far!" << std::endl;
    for( size_t el=0; el<GoodEvents.size(); ++el) {
      std::cout << "Tree index " << GoodEvents[el].first << " was a good event, which had " << GoodEvents[el].second << " ADC values." << std::endl;
    }
    std::cout << std::endl;
  }
  return true;
} // read next

//=======================================================================================
void DAQToOffline::Splitter::closeCurrentFile() {
  file_.reset(nullptr);
}

//=======================================================================================
bool DAQToOffline::Splitter::eventIsFull_( vector<RawDigit> const & v ) {
  return v.size() == ticksPerEvent_;
}

//=======================================================================================
bool DAQToOffline::Splitter::loadEvents_( size_t &InputTree ) {
  //std::cout << "At the start of loadEvents..." << std::endl;
  if ( InputTree != nInputEvts_ ) {
    if ( !loadedDigits_.empty(fDebugLevel) && fRequireRCE ) return false;
    
    // I want to look through my map to find correct tree for this event, whilst ensuring I skip the neccessary number of events....
    bool GoodTree = false;
    size_t LoadTree = 0;
    while (!GoodTree) {
      int LookingAtIndex = 0;
      for (std::map<uint64_t,size_t>::iterator it=EventTreeMap.begin(); it!=EventTreeMap.end(); ++it ) {
	++LookingAtIndex;
	if ( LookingAtIndex == (int)InputTree ) {
	  // Is this index high enough?
	  if (LookingAtIndex > fSkipNInputEvents ) {
	    LoadTree = it->second;
	    GoodTree = true;
	  }
	  break;
	}
      } // For Loop
      if (fDebugLevel > 2)
	std::cout << "Looking for event " << InputTree << ", it was found at TreeIndex " << LookingAtIndex << " but I want to skip " << fSkipNInputEvents << ". Do I load this tree? " << GoodTree << std::endl;
      if (!GoodTree) ++InputTree;
    } // While Loop
    // I want to look through my map to find correct tree for this event!
      
    EventAuxBranch_->GetEntry(LoadTree);
    inputRunNumber_ = evAux_.run();
    inputSubRunNumber_ = evAux_.subRun();
    inputEventNumber_ = evAux_.event();
    inputEventTime_ = evAux_.time();
    if (fDebugLevel > 1)
      std::cout << "\nLoading event " << inputEventNumber_ << " on Tree " << InputTree << std::endl;

    bool PTBTrigPresent = false;
    if (fRequirePTB) {
      PTBTrigPresent= LoadPTBInformation( LoadTree );
      if (fDebugLevel > 1) std::cout << "Is there a PTB trigger present? " << PTBTrigPresent << std::endl;
    }
    
    if ( fWhichTrigger.size() == 1 && fWhichTrigger[0] == 3 ) { // If looking for only looking for triggers from the PTB
      if ( PTBTrigPresent || fTrigger ) { // Only load RCEs and SSPs if a trigger is present.
	if (fRequireSSP) LoadSSPInformation( LoadTree );
	if (fRequireOpHit) LoadOpHitInformation( LoadTree);
	if (fRequireRCE) LoadRCEInformation( LoadTree );
      } else {                // If no PTB trigger is present make sure to clear Digits and Waveforms.
	if (fRequireSSP) loadedWaveforms_.clear(fDebugLevel);
	if (fRequireOpHit) loadedOpHits_.clear(fDebugLevel);
	if (fRequireRCE) loadedDigits_.clear(fDebugLevel);
      }
    } else { // If either not looking for PTB triggers at all, or looking for additional triggers too, then always load RCE/SSP.
      if (fRequireSSP) LoadSSPInformation( LoadTree );
      if (fRequireOpHit) LoadOpHitInformation( LoadTree);
      if (fRequireRCE) LoadRCEInformation( LoadTree );
    }
    // If Looking at MonteCarlo, also want to load in some other stuff....
    if (fMonteCarlo) {
      auto *particles = getMCParticles(MCPartinputBranch_, LoadTree);
      MonteCarlo_.loadMCParts(*particles);
      if (fDebugLevel) std::cout << "Particles has size " << particles->size() << std::endl;
      auto *simchans = getMCSimChans(MCSimChaninputBranch_, LoadTree);
      MonteCarlo_.loadSimChans(*simchans);
      if (fDebugLevel) std::cout << "SimChans has size " << simchans->size() << std::endl;
    }
    
    InputTree++;
    return true;
  }
  else return false;
} // load digits
//=======================================================================================
bool DAQToOffline::Splitter::LoadPTBInformation( size_t LoadTree ) {
  bool TrigPresent = false;
  if (PenninputDataProduct_.find("Fragment") != std::string::npos) {
    auto* PennFragments = getFragments ( PenninputBranch_, LoadTree );
    std::vector<raw::ExternalTrigger> counters = DAQToOffline::PennFragmentToExternalTrigger( *PennFragments, fPTBMap );
    loadedCounters_.load( counters, fDebugLevel );
    if (fDebugLevel > 1) std::cout << "Counters has size " << counters.size() << std::endl;
    
    for (size_t CountLoop = 0; CountLoop < counters.size(); ++CountLoop) {
      if (fDebugLevel > 3 )
	std::cout << "Looking at counters[" << CountLoop << "] has CounterID " << counters[CountLoop].GetTrigID() << " and Timestamp " << counters[CountLoop].GetTrigTime() << std::endl;
      for (size_t PTB = 0; PTB < fPTBTrigs.size(); ++PTB) {
	if ( counters[CountLoop].GetTrigID() == fPTBTrigs[PTB] ) {
	  if (fDebugLevel) std::cout << "Looking at event " << inputEventNumber_ << ", there is a trigger here on channel " << counters[CountLoop].GetTrigID() << std::endl;
	  TrigPresent = true;
	}
      }
    }
  } else {
    auto *counters = getRawExternalTriggers(PenninputBranch_, LoadTree );
    loadedCounters_.load( *counters, fDebugLevel );
    if (fDebugLevel > 1) {
      std::cout << "Counters has size" << counters->size() << std::endl;
    }

    for (auto count: *counters) {
      if (fDebugLevel > 3 )
	std::cout << "Looking at a counter which has CounterID " << count.GetTrigID() << " and Timestamp " << count.GetTrigTime() << std::endl;
      for (size_t PTB = 0; PTB < fPTBTrigs.size(); ++PTB) {
	if ( count.GetTrigID() == fPTBTrigs[PTB] ) {
	  if (fDebugLevel) std::cout << "Looking at event " << inputEventNumber_ << ", there is a trigger here on channel " << count.GetTrigID() << std::endl;
	  TrigPresent = true;
	}
      }
    }
  }
  return TrigPresent;
}
//=======================================================================================
void DAQToOffline::Splitter::LoadSSPInformation( size_t LoadTree ) {
  if (SSPinputDataProduct_.find("Fragment") != std::string::npos) {
    auto* SSPfragments = getFragments( SSPinputBranch_, LoadTree );
    std::vector<raw::OpDetWaveform> waveforms = sspReform.SSPFragmentToOpDetWaveform(*SSPfragments);
    for ( size_t WaveLoop=0; WaveLoop < waveforms.size(); ++WaveLoop ) {
      int64_t SSPTime = waveforms[WaveLoop].TimeStamp()*fNovaTicksPerSSPTick;
      if (fDebugLevel > 3 )
	std::cout << "Looking at waveform[" << WaveLoop << "] it has channel number " << waveforms[WaveLoop].ChannelNumber()
		  << " and timestamp " << SSPTime << ", and size " << waveforms[WaveLoop].size() << std::endl;
    }
    
    if (fDebugLevel > 1) std::cout << "Loaded waveforms has size " << waveforms.size() << std::endl;
    loadedWaveforms_.load( waveforms, fDebugLevel );
  }
  else {
    auto* waveforms = getSSPWaveforms(SSPinputBranch_, LoadTree );
    if (fDebugLevel > 1) std::cout << "Loaded waveforms has size " << waveforms->size() << std::endl;
    loadedWaveforms_.load( *waveforms, fDebugLevel );
  }
  return;
}
//=======================================================================================
void DAQToOffline::Splitter::LoadOpHitInformation( size_t LoadTree ) {
  if (OpHitinputDataProduct_.find("Fragment") != std::string::npos) {
    auto* OpHitfragments = getFragments( OpHitinputBranch_, LoadTree );
    std::vector<recob::OpHit> OpHits = sspReform.SSPHeaderToOpHit(*OpHitfragments);
    for ( size_t HitLoop=0; HitLoop < OpHits.size(); ++HitLoop ) {
      int64_t OpHitTime = OpHits[HitLoop].PeakTime()*fNovaTicksPerSSPTick;
      if (fDebugLevel > 3 )
	std::cout << "Looking at waveform[" << HitLoop << "] it is on channel number " << OpHits[HitLoop].OpChannel() << " at timestamp " << OpHitTime << std::endl;
    }
    
    if (fDebugLevel > 1) std::cout << "Loaded OpHits has size " << OpHits.size() << std::endl;
    loadedOpHits_.load( OpHits, fDebugLevel );
  }
  else {
    auto* OpHits = getOpHitWaveforms(OpHitinputBranch_, LoadTree );
    if (fDebugLevel > 1) std::cout << "Loaded waveforms has size " << OpHits->size() << std::endl;
    loadedOpHits_.load( *OpHits, fDebugLevel );
  }
  return;
}
//=======================================================================================
void DAQToOffline::Splitter::LoadRCEInformation( size_t LoadTree ) {
  if (TPCinputDataProduct_.find("Fragment") != std::string::npos) {
    lbne::TpcNanoSlice::Header::nova_timestamp_t firstTimestamp = 0;
    auto* fragments = getFragments( TPCinputBranch_, LoadTree );
    rawDigits_t const digits = fragmentsToDigits_( *fragments, firstTimestamp, TPCChannelMap );
    if (!digits.size() ) {
      RCEsNotPresent = true;
    }
    loadedDigits_.load( digits, fDebugLevel );
    loadedDigits_.loadTimestamp( firstTimestamp );
    if (fDebugLevel > 1) std::cout << "Loaded RCE information with timestamp " << firstTimestamp << std::endl;
  }
  else {
    auto* digits = getRawDigits(TPCinputBranch_, LoadTree );
    loadedDigits_.load( *digits, fDebugLevel);
    loadedDigits_.loadTimestamp(0); // MC timestamp is zero (? assume?)
    if (fDebugLevel > 1) std::cout << "Loaded MC RCE's with timestamp 0" << std::endl;
  }
}
//=======================================================================================
void DAQToOffline::Splitter::makeEventAndPutDigits_(art::EventPrincipal*& outE, art::Timestamp art_timestamp) {
  // just keep incrementing the event number as we split along
  ++eventNumber_;
  
  for (size_t TrigSize = 0; TrigSize < fWhichTrigger.size(); ++TrigSize) {
    if ( fWhichTrigger[TrigSize] == 0 && fDebugLevel) {
      std::cout << "\n\n\nI hope you know that you are triggering on a random number of ticks and not any sort of data! Check that fwhichTrigger is set correctly.\n\n\n" << std::endl;
    }
  }
  std::cout << "Making an event with RunNumber " << runNumber_ << ", subRunNumber " << subRunNumber_ << ", EventNumber " << eventNumber_ << " and art_timestamp " << art_timestamp.value() << std::endl;

  outE = sh_.makeEventPrincipal( runNumber_, subRunNumber_, eventNumber_, art_timestamp );
  art::put_product_in_principal( std::make_unique<rawDigits_t>(bufferedDigits_),
                                 *outE,
                                 sourceName_,
                                 TPCinputTag_.instance() );
  art::put_product_in_principal( std::make_unique<SSPWaveforms_t>(wbuf_),
                                 *outE,
                                 sourceName_,
                                 SSPinputTag_.instance() );
  art::put_product_in_principal( std::make_unique<OpHits_t>(hbuf_),
                                 *outE,
                                 sourceName_,
                                 OpHitinputTag_.instance() );
  art::put_product_in_principal( std::make_unique<PennCounters_t>(cbuf_),
                                 *outE,
                                 sourceName_,
                                 PenninputTag_.instance() );

  if (fMonteCarlo) {
    art::put_product_in_principal( std::make_unique<MCPart_t>(mcbuf_),
				   *outE,
				   sourceName_,
				   MCPartinputTag_.instance() );
    art::put_product_in_principal( std::make_unique<MCSimChan_t>(simchanbuf_),
				   *outE,
				   sourceName_,
				   MCSimChaninputTag_.instance() );
  }
  

  mf::LogDebug("SplitterFunc") << "Producing event: " << outE->id() << " with " << bufferedDigits_.size() << " RCE digits and " <<
    wbuf_.size() << " SSP waveforms, " << hbuf_.size() << " OpHits and " << cbuf_.size() << " External Triggers (muon counters)";
  Reset();
}
//=======================================================================================
void DAQToOffline::Splitter::Reset() {
  bufferedDigits_.clear();
  for (size_t ichan=0;ichan<dbuf_.size();ichan++) { dbuf_[ichan].clear(); }
  dbuf_.clear();
  wbuf_.clear();
  hbuf_.clear();
  cbuf_.clear();
  Event_timestamp = 0;
  fTicksAccumulated = 0; // No longer have any RCE data...
  fTrigger = false;      // Need to re-decide where to trigger
  fDiffFromLastTrig = 0; // Reset trigger counter.
  if (fDebugLevel > 1) std::cout << "Resetting everything (dbuf, cbuf, wbuf, hbuf, Trigger, etc)" << std::endl;
}
//=======================================================================================
void DAQToOffline::Splitter::CheckTimestamps(bool &JumpEvent, size_t &JumpNADC ) {
  int StampDiff = loadedDigits_.getTimeStampAtIndex(loadedDigits_.index, fNovaTicksPerTPCTick) - prev_timestamp;
  if (fDebugLevel) std::cout << "\n" << StampDiff << " = " << loadedDigits_.getTimeStampAtIndex(loadedDigits_.index, fNovaTicksPerTPCTick) << " - " << prev_timestamp << std::endl;
  if ( fabs(StampDiff) > fTimeStampThreshold ) { // Timestamps of old and new file too far apart. So want to clear all previously loaded event.
    if (fDebugLevel) std::cout << "\nThe absolute gap between timestamps is " << fabs(StampDiff) << " which is more than the threshold " << fTimeStampThreshold << std::endl;
    bool fixed = false;
    if ( StampDiff < 0 ) { // got overlapping timestamps.
      if (fDebugLevel) std::cout << "Stamp diff is negative...need to figure out how to try and fix..." << std::endl;
      //fixed = true;
    } else {
      if (fDebugLevel) std::cout << "Stamp diff is positive...need to figure out how to try and fix..." << std::endl;
      //fixed=true;
    }
    if (!fixed) {
      if (fDebugLevel) std::cout << "\nCan't reconcile the timestamps, so voiding this trigger :( \n" << std::endl;
      Reset();
      loadedDigits_.index = fPreTriggerTicks;
      if (fPreTriggerTicks > loadedDigits_.digits[0].NADC() ) {
	JumpEvent = true;
	JumpNADC  = loadedDigits_.digits[0].NADC();
      }
      fDiffFromLastTrig   = fPreTriggerTicks;
      this_timestamp      = loadedDigits_.getTimeStampAtIndex(loadedDigits_.index, fNovaTicksPerTPCTick);
      fLastTriggerIndex = 0;
    } else {
      if (fDebugLevel) std::cout << "\nRectified the timestamps, carry on building event :D\n" << std::endl;
    }
  } else {
    if (fDebugLevel) std::cout << "\nTimestamps lead on from each other, carry on :)\n" << std::endl;
  }
} // Check Timestamps
//=======================================================================================
bool DAQToOffline::Splitter::NoRCEsCase(art::RunPrincipal*& outR, art::SubRunPrincipal*& outSR, art::EventPrincipal*& outE) {
  while (!fTrigger) {
    bool NewTree = false;
    // Whilst LoadedWaveforms and LoadedCounters are empty, load a new event...
    while ( loadedWaveforms_.empty() || loadedCounters_.empty() ) {
      bool rc = loadEvents_(treeIndex_);
      if (!rc) {
	doneWithFiles_ = (file_->GetName() == lastFileName_);
	return false;
      }
      NewTree = true; 
    } // while empty
    std::map<int,int> PrevChanADC;
    std::vector<short> ADCdigits;
    Triggering(PrevChanADC, ADCdigits, NewTree);
  }
  wbuf_ = loadedWaveforms_.TakeAll();
  hbuf_ = loadedOpHits_.TakeAll();
  cbuf_ = loadedCounters_.TakeAll();
  if (fDebugLevel > 2) std::cout << "After looking at treeIndex_ " << treeIndex_-1 << " fTrigger is " << fTrigger << " and wbuf and cbuf have sizes " << wbuf_.size() << " and " << cbuf_.size() << std::endl;
  
  // ******** Now Build the event *********
  runNumber_ = inputRunNumber_;
  subRunNumber_ = inputSubRunNumber_;
  //art::Timestamp ts; // LBNE should decide how to initialize this -- use first_timestamp converted into an art::Timestamp
  //FIXME - This is a first attempt at interpreting the novatimestamp from the tpc data to create an art event timestamp
  art::Timestamp this_art_event_timestamp = DAQToOffline::make_art_timestamp_from_nova_timestamp(Event_timestamp);
  if ( runNumber_ != cachedRunNumber_ ) {
    outR = sh_.makeRunPrincipal(runNumber_,this_art_event_timestamp);
    cachedRunNumber_ = runNumber_;
    eventNumber_ = 0ul;
  }
  if ( subRunNumber_ != cachedSubRunNumber_ ) {
    outSR = sh_.makeSubRunPrincipal(runNumber_,subRunNumber_,this_art_event_timestamp);
    cachedSubRunNumber_ = subRunNumber_;
    eventNumber_ = 0ul;
  }
  //inputEventTime_ is the art::Timestamp() of the online art::Event() used to create the offline art::Event()
  makeEventAndPutDigits_( outE, inputEventTime_ );
  Reset();
  loadedWaveforms_.clear(fDebugLevel);
  loadedCounters_.clear(fDebugLevel);
  std::cout << "After reset() is Waveforms empty? " << loadedWaveforms_.empty() << ", what about Counters? " << loadedCounters_.empty() << std::endl;
  return true;
}
//=======================================================================================
void DAQToOffline::Splitter::CheckTrigger() {
  double TempTriggerIndex = loadedDigits_.index;
  size_t TempTreeIndex    = treeIndex_ -1;
  lbne::TpcNanoSlice::Header::nova_timestamp_t TempTimeStamp = this_timestamp;
  size_t TempNADCs        = loadedDigits_.digits[0].NADC();
  art::EventNumber_t TempEventNumber = inputEventNumber_;
  if (fDebugLevel) std::cout << "\nTrying to Trigger on timestamp " << this_timestamp << ", last trigger was on " << fLastTimeStamp << "...." << this_timestamp - fLastTimeStamp << std::endl;
  
  //******** Now to sort out the prebuffer!!! ***********
  int BufferResidual =  loadedDigits_.index - fPreTriggerTicks;
  if ( BufferResidual > 0 ) {
    if (fDebugLevel) {
      std::cout << "I have enough previous digits in this event for the prebuffer (" << BufferResidual << " = " << loadedDigits_.index << " - " << fPreTriggerTicks 
		<< ")! Moving loadedDigits_.index to " << BufferResidual << std::endl;
    }
    loadedDigits_.index = BufferResidual;
  }
  else { // Don't have enough ticks in the event for the prebuffer :( so need to load previous events!
    BufferResidual = -BufferResidual;
    lbne::TpcNanoSlice::Header::nova_timestamp_t TrigEvStart = loadedDigits_.getTimeStampAtIndex(loadedDigits_.index, fNovaTicksPerTPCTick);
    fTrigger = false;
    if (fDebugLevel) {
      std::cout << "I don't have enough previous digits :(, I need an extra " << BufferResidual << " ticks from previous events. TrigEvStart = " << TrigEvStart << std::endl;
    }
    
    size_t LoadEv = 0, LoadInd = 0;
    for ( size_t el=2; el<GoodEvents.size()+1; ++el) {
      if (fDebugLevel > 3) {
	std::cout << "Going backwards...Tree index " << GoodEvents[GoodEvents.size()-el].first << " was a good event, which had "
		  << GoodEvents[GoodEvents.size()-el].second << " ADC values."
		  << "\nThis means TempBuffer has gone from " << BufferResidual << " to " << BufferResidual-(int)GoodEvents[GoodEvents.size()-el].second
		  << std::endl;
      }
      if (BufferResidual-(int)GoodEvents[GoodEvents.size()-el].second < 0) {
        if (fDebugLevel) std::cout << "I can satisfy the prebuffer on treeIndex " << GoodEvents[GoodEvents.size()-el].first << " at index " << (int)GoodEvents[GoodEvents.size()-el].second - BufferResidual << std::endl;
        LoadEv = GoodEvents[GoodEvents.size()-el].first;
        LoadInd= (int)GoodEvents[GoodEvents.size()-el].second - BufferResidual;
        fTrigger = true;
        break;
      } else {
        BufferResidual = BufferResidual-(int)GoodEvents[GoodEvents.size()-el].second;
        if (fDebugLevel > 1) std::cout << "Can't quite trigger on treeIndex " << GoodEvents[GoodEvents.size()-el].first << ", need another " << BufferResidual << "ticks." << std::endl;
      }
    }
    if (fTrigger) {
      loadedDigits_.index = 0;
      loadedDigits_.clear(fDebugLevel);
      treeIndex_ = LoadEv;
      loadEvents_(treeIndex_);
      loadedDigits_.index = LoadInd;
    }
  }
  this_timestamp = loadedDigits_.getTimeStampAtIndex(loadedDigits_.index, fNovaTicksPerTPCTick); // Set this_timestamp to whatever it now is....
  
  if ( fTrigger ) { // If trigger is still good!
    fLastTriggerIndex = TempTriggerIndex;
    fLastTreeIndex    = TempTreeIndex;
    fLastTimeStamp    = TempTimeStamp;
  }
  else {
    loadedDigits_.index = TempTriggerIndex + BufferResidual; // Jump to where the trigger was plus buffer residual
    if (fDebugLevel) {
      std::cout << "Trigger isn't good so I'm going back to where I triggered..."
		<< "Attempted trigger was in event " << TempEventNumber << " at index " << TempTriggerIndex
		<< " at timestamp " << TempTimeStamp << ", it had " << TempNADCs << " adc's"
		<< "\nI'm now at event " << inputEventNumber_ << " index " << loadedDigits_.index
		<< " and timestamp " << loadedDigits_.getTimeStampAtIndex(loadedDigits_.index, fNovaTicksPerTPCTick) 
		<< " and " << loadedDigits_.digits[0].NADC() << " adcs, loadedDigits empty? " << loadedDigits_.empty(fDebugLevel) << "\n"
		<< std::endl;
    }
  }
}
//===================================================================================================================================
void DAQToOffline::Splitter::Triggering(std::map<int,int> &PrevChanADC, std::vector<short> ADCdigits, bool NewTree) {
  if ( treeIndex_-1 != fLastTreeIndex ) fLastTimeStamp = 0; // No longer looking at same treeIndex as previous trigger, so reset lastTimeStamp
    
  for (size_t TrigSize = 0; TrigSize < fWhichTrigger.size(); ++TrigSize) {
    // Trigger on Monte Carlo whichTrigger == 0
    if ( fWhichTrigger[TrigSize] == 0 ) {
      if ( fDiffFromLastTrig > fMCTrigLevel ) fTrigger = true;
    }
    // Trigger on new files
    else if ( fWhichTrigger[TrigSize] == 1 && NewTree ) {
      fTrigger = true;
    }
    // Trigger on Photon Detectors
    else if ( fWhichTrigger[TrigSize] == 2 ) {
      fTrigger = loadedWaveforms_.PhotonTrigger( prev_timestamp, fWaveformADCThreshold, fWaveformADCsOverThreshold, fWaveformADCWidth, fNovaTicksPerSSPTick, fDebugLevel );
      }
    // Trigger using the user defined PTB Trigger vector
    else if ( fWhichTrigger[TrigSize] == 3 ) {
      fTrigger = loadedCounters_.PTBTrigger( this_timestamp, fNovaTicksPerCountTick, fNovaTicksPerTPCTick, fDebugLevel, fPTBTrigs);
    }
    // Trigger on "Tickler" / TPC information
    else if ( fWhichTrigger[TrigSize] == 4 ) {
      fTrigger = TicklerTrigger( PrevChanADC, ADCdigits);
    }
    
    if (fTrigger && fRequireRCE) CheckTrigger();
    if (fTrigger) break;
  }
}
//===================================================================================================================================
bool DAQToOffline::Splitter::TicklerTrigger( std::map<int,int> &PrevChanADC, std::vector<short> ADCdigits ) {
  int HitsOverThreshold = 0;
  if (PrevChanADC.size() != 0) {
    for (unsigned int achan=0; achan<ADCdigits.size(); ++achan)
      if ( fabs( ADCdigits[achan] - PrevChanADC[achan] ) > fADCdiffThreshold ) {
        ++HitsOverThreshold;
        if (fDebugLevel > 3) {
	  std::cout << "Looking at index " << loadedDigits_.index << " channel " << achan << "..."  << ADCdigits[achan] << " - "
		    << PrevChanADC[achan] << " = " << fabs( ADCdigits[achan] - PrevChanADC[achan] ) << " > " << fADCdiffThreshold
		    << std::endl;
	}
      }
    if ( HitsOverThreshold != 0 )
      if (fDebugLevel > 2) std::cout << " after looking through all the channels ("<<ADCdigits.size()<<") I had " << HitsOverThreshold << " ticks with diff more than " << fADCdiffThreshold << std::endl;
    if ( HitsOverThreshold > fADCsOverThreshold ) {
      if (fDebugLevel > 2) std::cout << "Looking at index " << loadedDigits_.index << ", which had " << HitsOverThreshold << " hits over diff threshold. Trigger threshold is " << fADCsOverThreshold << std::endl;
      return true;
    }
  } // if PrevChanADC not empty.
  for (unsigned int bchan=0; bchan<ADCdigits.size(); ++bchan)
    PrevChanADC[bchan] = ADCdigits[bchan];
  return false;
}
//===================================================================================================================================
DEFINE_ART_INPUT_SOURCE(art::Source<DAQToOffline::Splitter>)
//===================================================================================================================================
