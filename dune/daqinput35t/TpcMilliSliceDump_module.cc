////////////////////////////////////////////////////////////////////////
// Class:       TpcMilliSliceDump
// Module Type: analyzer
// File:        TpcMilliSliceDump_module.cc
// Description: Prints out information about each event.
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Optional/TFileService.h"

#include "canvas/Utilities/Exception.h"
#include "lbne-raw-data/Overlays/TpcMilliSliceFragment.hh"
#include "artdaq-core/Data/Fragment.hh"

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
class TpcMilliSliceDump;
}

class lbne::TpcMilliSliceDump : public art::EDAnalyzer {
public:
	explicit TpcMilliSliceDump(fhicl::ParameterSet const & pset);
	virtual ~TpcMilliSliceDump();

	virtual void analyze(art::Event const & evt);

private:
	void beginJob();

	std::string raw_data_label_;
	std::string frag_type_;
	std::vector<int> verb_microslice_ids_;
	std::vector<int> verb_nanoslice_ids_;
	uint32_t         verb_nanoslice_adcs_;
	TH1D * adc_values_;
};


lbne::TpcMilliSliceDump::TpcMilliSliceDump(fhicl::ParameterSet const & pset)
: EDAnalyzer(pset),
  raw_data_label_(pset.get<std::string>("raw_data_label")),
  frag_type_(pset.get<std::string>("frag_type")),
  verb_microslice_ids_(pset.get<std::vector<int>>("verbose_microslice_ids", std::vector<int>(1,0))),
  verb_nanoslice_ids_ (pset.get<std::vector<int>>("verbose_nanoslice_ids",  std::vector<int>(1,0))),
  verb_nanoslice_adcs_(pset.get<uint32_t>        ("verbose_nanoslice_adcs", 6)),
  adc_values_(nullptr)
{
}

void lbne::TpcMilliSliceDump::beginJob()
{
	art::ServiceHandle<art::TFileService> tfs;
	adc_values_ = tfs->make<TH1D>("adc_values","adc_values",4096,-0.5,4095.5);
}

lbne::TpcMilliSliceDump::~TpcMilliSliceDump()
{
}

void lbne::TpcMilliSliceDump::analyze(art::Event const & evt)
{
	//counter of total number of ALL adc values in the event
	uint32_t n_adc_counter(0);
	//cumulative total of ALL adc values in the event
	uint64_t adc_cumulative(0);

	art::EventNumber_t eventNumber = evt.event();

	// ***********************
	// *** TpcMilliSlice Fragments ***
	// ***********************

	// look for raw TpcMilliSlice data

	art::Handle<artdaq::Fragments> raw;
	evt.getByLabel(raw_data_label_, frag_type_, raw);

	if (raw.isValid()) {
		std::cout << "######################################################################" << std::endl;
		std::cout << std::endl;
		std::cout << "Run " << evt.run() << ", subrun " << evt.subRun()
            		  << ", event " << eventNumber << " has " << raw->size()
            		  << " fragment(s) of type " << frag_type_ << std::endl;

		for (size_t idx = 0; idx < raw->size(); ++idx)
		{
			const auto& frag((*raw)[idx]);

			TpcMilliSliceFragment msf(frag);

			std::cout << std::endl;
			std::cout << "TpcMilliSlice fragment " << frag.fragmentID() << " consists of: " << msf.size() << " bytes containing "
					<< msf.microSliceCount() << " microslices" << std::endl;
			std::cout << std::endl;

			for (uint32_t i_ms = 0; i_ms < msf.microSliceCount(); ++i_ms)
			{
				bool verb_microslice = (std::find(verb_microslice_ids_.begin(), verb_microslice_ids_.end(), i_ms) != verb_microslice_ids_.end()) ? true : false;

				std::unique_ptr<const TpcMicroSlice> microslice = msf.microSlice(i_ms);

				if (!microslice) {
					throw cet::exception("Error in TpcMilliSliceDump module: unable to find requested microslice");
				}

				if(verb_microslice) {
					std::cout << "TpcMilliSlice fragment " << frag.fragmentID()
			    		<< ", microslice " << i_ms
			    		<< " has sequence ID " << microslice->sequence_id()
			    		<< " size " << microslice->size()
			    		<< " and consists of: "<< microslice->nanoSliceCount() << " nanoslices" << std::endl;
				}

				if (microslice->nanoSliceCount() == 0) continue;

				// First pass through nanoslices listed for verbose output. This is faster than iterating through many nanoslices
				// and searching for the index in the verbose list

				for (std::vector<int>::iterator i_nano = verb_nanoslice_ids_.begin(); i_nano != verb_nanoslice_ids_.end(); i_nano++)
				{
					bool verb_nanoslice_header = verb_microslice;
					if (verb_nanoslice_header && (std::find(verb_nanoslice_ids_.begin(), verb_nanoslice_ids_.end(), *i_nano) == verb_nanoslice_ids_.end()))
					{
						verb_nanoslice_header = false;
					}

					std::unique_ptr<const TpcNanoSlice> nanoslice = microslice->nanoSlice(*i_nano);
					if (! nanoslice )
					{
						throw cet::exception("Error in TpcMilliSliceDump module: Unable to locate requested nanoslice");
					}
					if (verb_nanoslice_header)
					{

						std::cout << "  NanoSlice " << *i_nano << " size " << nanoslice->size() << " bytes" << std::endl << std::endl;

						const uint64_t* nanoslice_raw = nanoslice->raw();

						std::cout << "    Raw Header : 0x"
								<< std::hex << std::setw(16) << std::setfill('0') << nanoslice_raw[0] <<std::dec << std::endl;

						std::cout << "    Raw Data Payload:" << std::endl << "    ";
						uint16_t num_nanoslice_words = (nanoslice->size() - sizeof(lbne::TpcNanoSlice::Header))/sizeof(uint64_t);
						size_t word_offset = sizeof(lbne::TpcNanoSlice::Header)/sizeof(uint64_t);

						for (int i_adc_word = 0; i_adc_word < num_nanoslice_words; i_adc_word++)
						{
							std::cout << "0x" << std::hex << std::setfill('0') << std::setw(16) << nanoslice_raw[word_offset+i_adc_word] << std::dec << " ";
							if ((i_adc_word % 4) == 3) std::cout << std::endl << "    ";
						}
						std::cout << std::endl;

						//std::cout << "    Decoded Header:" << std::endl;
						std::cout << "    NOvA timestamp   : 0x" << std::hex << std::setw(14) << std::setfill('0')
							      << nanoslice->nova_timestamp() << std::dec << std::endl;


						std::cout << "    Decoded ADC Values:" << std::endl << "    ";
						for (uint32_t i_samp = 0; i_samp < nanoslice->getNChannels(); ++i_samp)
						{
							uint16_t val = std::numeric_limits<uint16_t>::max();
							nanoslice->sampleValue(i_samp, val);
							std::cout << std::setfill(' ') << std::setw(4) << i_samp << " : " << std::setw(4) << val << " ";
							if ((i_samp % 8) == 7) std::cout << std::endl << "    ";

						} // i_samp
						std::cout << std::endl;
					}
				} // i_nano

				std::cout << std::endl;

				// Make a second pass through all nanoslices to build some statistics and fill histograms
				for (uint32_t i_nano = 0; i_nano < microslice->nanoSliceCount(); i_nano++)
				{
					std::unique_ptr<const TpcNanoSlice> nanoslice = microslice->nanoSlice(i_nano);

					if (!nanoslice)
					{
						throw cet::exception("Error in TpcMilliSliceDump module: Unable to find requested nanoslice");
					}

					for (uint32_t i_samp = 0; i_samp < nanoslice->sampleCount(); i_samp++)
					{
						uint16_t val = std::numeric_limits<uint16_t>::max();
						bool success = nanoslice->sampleValue(i_samp, val);

						if (!success)
						{
							throw cet::exception("Error in TpcMilliSliceDump module: Unable to find the requested sample value in nanoslice");
						}
						n_adc_counter++;
						adc_cumulative += (uint64_t)val;
						adc_values_->Fill(val);
					} // i_samp
				} // i_nano
			} // i_ms

		} //idx

		std::cout << std::endl
				<< "Event ADC average is (from counter):   " << (double)adc_cumulative/(double)n_adc_counter
				<< std::endl
				<< "Event ADC average is (from histogram): " << adc_values_->GetMean()
				<< std::endl;

	}//raw.IsValid()?
	else {
		std::cout << "Run " << evt.run() << ", subrun " << evt.subRun()
				  << ", event " << eventNumber << " has zero"
				  << " TpcMilliSlice fragments." << std::endl;
	}
	std::cout << std::endl;

}

DEFINE_ART_MODULE(lbne::TpcMilliSliceDump)
