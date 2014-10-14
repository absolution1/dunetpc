////////////////////////////////////////////////////////////////////////
// Class:       SSPDump
// Module Type: analyzer
// File:        SSPDump_module.cc
// Description: Prints out information about each event.
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Optional/TFileService.h"

#include "art/Utilities/Exception.h"
//#include "lbne-artdaq/Overlays/SSPFragment.hh"
#include "lbne-raw-data/Overlays/SSPFragment.hh"
#include "artdaq-core/Data/Fragments.hh"
//#include "lbne-artdaq/Generators/anlBoard/anlTypes.h"

#include "TH1.h"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <vector>
#include <iostream>
#include <limits>

namespace lbne {
  class SSPDump;
}

class lbne::SSPDump : public art::EDAnalyzer {
public:
  explicit SSPDump(fhicl::ParameterSet const & pset);
  virtual ~SSPDump();

  virtual void analyze(art::Event const & evt);

private:
  void beginJob() override;
  void endJob  () override;
  void beginEvent(art::EventNumber_t eventNumber);
  void endEvent  (art::EventNumber_t eventNumber);

  std::string raw_data_label_;
  std::string frag_type_;

  //sets the verbosity of std::cout printing
  //std::vector<int> verb_microslice_ids_;
  //std::vector<int> verb_nanoslice_ids_;
  uint32_t         verb_adcs_;
  bool             verb_meta_;

  //histograms, counters, etc
  TH1D * adc_values_;
  uint32_t n_adc_counter_;  //counter of total number of ALL adc values in an event
  uint64_t adc_cumulative_; //cumulative total of ALL adc values in an event
};


lbne::SSPDump::SSPDump(fhicl::ParameterSet const & pset)
    : EDAnalyzer(pset),
      raw_data_label_(pset.get<std::string>("raw_data_label")),
      frag_type_(pset.get<std::string>("frag_type")),
      //verb_microslice_ids_(pset.get<std::vector<int>>("verbose_microslice_ids", std::vector<int>(1,0))),
      //verb_nanoslice_ids_ (pset.get<std::vector<int>>("verbose_nanoslice_ids",  std::vector<int>(1,0))),
      verb_adcs_(pset.get<uint32_t>        ("verbose_adcs", 15)),
      verb_meta_(pset.get<bool>            ("verbose_metadata", true)),
      adc_values_(nullptr),
      n_adc_counter_(0),
      adc_cumulative_(0)
{
}

void lbne::SSPDump::beginJob()
{
  art::ServiceHandle<art::TFileService> tfs;
  adc_values_ = tfs->make<TH1D>("adc_values","adc_values",4096,-0.5,4095.5);  
}

void lbne::SSPDump::beginEvent(art::EventNumber_t /*eventNumber*/)
{
  //reset ADC histogram
  adc_values_->Reset();
  //reset counters
  n_adc_counter_  = 0;
  adc_cumulative_ = 0;
}

void lbne::SSPDump::endEvent(art::EventNumber_t eventNumber)
{
  //write the ADC histogram for the given event
  if(n_adc_counter_)
    adc_values_->Write(Form("adc_values:event_%d", eventNumber));
}

void lbne::SSPDump::endJob()
{
  delete adc_values_;
}

lbne::SSPDump::~SSPDump()
{
}

void lbne::SSPDump::analyze(art::Event const & evt)
{

  art::EventNumber_t eventNumber = evt.event();
  beginEvent(eventNumber);

  // ***********************
  // *** SSP Fragments ***
  // ***********************

  // look for raw SSP data

  art::Handle<artdaq::Fragments> raw;
  evt.getByLabel(raw_data_label_, frag_type_, raw);

  if (raw.isValid()) {
    std::cout << "######################################################################" << std::endl;
    std::cout << std::endl;
    std::cout << "Run " << evt.run() << ", subrun " << evt.subRun()
              << ", event " << eventNumber << " has " << raw->size()
              << " fragment(s) of type " << frag_type_ << std::endl;

    for (size_t idx = 0; idx < raw->size(); ++idx) {
      const auto& frag((*raw)[idx]);

      SSPFragment sspf(frag);

      std::cout << std::endl;
      std::cout << "SSP fragment "     << frag.fragmentID() 
		<< " has total size: " << sspf.hdr_event_size()
		<< " and run number: " << sspf.hdr_run_number()
		<< " with " << sspf.total_adc_values() << " total ADC values"
		<< std::endl;
      std::cout << std::endl;

      //get the information from the header
      if(frag.hasMetadata())
      {
	const SSPFragment::Metadata * meta = frag.metadata<SSPFragment::Metadata>();
	// peaksum arrives as a signed 24 bit number.  It must be sign extended when stored as a 32 bit int 
	uint32_t peaksum = ((meta->daqHeader.group3 & 0x00FF) >> 16) + meta->daqHeader.peakSumLow;
	if(peaksum & 0x00800000) {
	  peaksum |= 0xFF000000;
	}
	if(verb_meta_) {
	  std::cout
	    << "Header:                             " << meta->daqHeader.header   << std::endl
	    << "Length:                             " << meta->daqHeader.length   << std::endl
	    << "Trigger type:                       " << ((meta->daqHeader.group1 & 0xFF00) >> 8) << std::endl
	    << "Status flags:                       " << ((meta->daqHeader.group1 & 0x00F0) >> 4) << std::endl
	    << "Header type:                        " << ((meta->daqHeader.group1 & 0x000F) >> 0) << std::endl
	    << "Trigger ID:                         " << meta->daqHeader.triggerID << std::endl
	    << "Module ID:                          " << ((meta->daqHeader.group2 & 0xFFF0) >> 4) << std::endl
	    << "Channel ID:                         " << ((meta->daqHeader.group2 & 0x000F) >> 0) << std::endl
	    << "Sync delay:                         " << ((unsigned int)(meta->daqHeader.timestamp[1]) << 16) + (unsigned int)(meta->daqHeader.timestamp[0]) << std::endl
	    << "Sync count:                         " << ((unsigned int)(meta->daqHeader.timestamp[3]) << 16) + (unsigned int)(meta->daqHeader.timestamp[2]) << std::endl
	    << "Peak sum:                           " << peaksum << std::endl
	    << "Peak time:                          " << ((meta->daqHeader.group3 & 0xFF00) >> 8) << std::endl
	    << "Prerise:                            " << ((meta->daqHeader.group4 & 0x00FF) << 16) + meta->daqHeader.preriseLow << std::endl
	    << "Integrated sum:                     " << ((unsigned int)(meta->daqHeader.intSumHigh) << 8) + (((unsigned int)(meta->daqHeader.group4) & 0xFF00) >> 8) << std::endl
	    << "Baseline:                           " << meta->daqHeader.baseline << std::endl
	    << "CFD Timestamp interpolation points: " << meta->daqHeader.cfdPoint[0] << " " << meta->daqHeader.cfdPoint[1] << " " << meta->daqHeader.cfdPoint[2] << " " << meta->daqHeader.cfdPoint[3] << std::endl
	    << "Internal interpolation point:       " << meta->daqHeader.intTimestamp[0] << std::endl
	    << "Internal timestamp:                 " << ((uint64_t)((uint64_t)meta->daqHeader.intTimestamp[3] << 32)) + ((uint64_t)((uint64_t)meta->daqHeader.intTimestamp[2] << 16)) + ((uint64_t)((uint64_t)meta->daqHeader.intTimestamp[1] << 0)) << std::endl
	    << std::endl;
	}
      }
      else
      {
	std::cout << "SSP fragment has no metadata associated with it." << std::endl;
      }

      //get the information from the data
      bool verb_values = true;
      for(size_t idata = 0; idata < sspf.total_adc_values(); idata++) {
	if(idata >= verb_adcs_)
	  verb_values = false;
	else if(idata == 0)
	  std::cout << "Data values: ";

	//const SSPFragment::adc_t * adc = sspf.dataBegin() + idata;
	const uint16_t * adc = (const uint16_t*)sspf.dataBegin() + idata;
	adc_values_->Fill(*adc);
	n_adc_counter_++;
	adc_cumulative_ += (uint64_t)adc;

	//std::cout << n_adc_counter_ << "\t" << adc_cumulative_ << std::endl;

	if(verb_values)
	  std::cout << *adc << " ";
      }
      std::cout << "..." << std::endl;

    }//idx

    std::cout << std::endl
	      << "Event ADC average is (from counter):   " << (double)adc_cumulative_/(double)n_adc_counter_
	      << std::endl
	      << "Event ADC average is (from histogram): " << adc_values_->GetMean()
	      << std::endl;
  }//raw.IsValid()?
  else {
    std::cout << "Run " << evt.run() << ", subrun " << evt.subRun()
              << ", event " << eventNumber << " has zero"
              << " SSP fragments." << std::endl;
  }
  std::cout << std::endl;

  beginEvent(eventNumber);
}

DEFINE_ART_MODULE(lbne::SSPDump)
