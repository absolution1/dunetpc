////////////////////////////////////////////////////////////////////////
// Class:       TriggerPrimitiveFinderPass1
// Plugin Type: service (art v2_10_03)
// File:        TriggerPrimitiveFinderPass1_service.cc
//
// Generated at Tue Jun  5 07:51:38 2018 by Philip Rodrigues using cetskelgen
// from cetlib version v3_02_00.
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Services/Registry/ActivityRegistry.h"
#include "art/Framework/Services/Registry/ServiceMacros.h"
#include "fhiclcpp/ParameterSet.h"

#include "dune/DAQSimAna/TriggerPrimitiveFinderService.h"

class TriggerPrimitiveFinderPass1 : public TriggerPrimitiveFinderService {
public:
  explicit TriggerPrimitiveFinderPass1(fhicl::ParameterSet const & p, art::ActivityRegistry & areg);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

    virtual std::vector<TriggerPrimitiveFinderService::Hit>
    findHits(const std::vector<unsigned int>& channel_numbers, 
             const std::vector<std::vector<short>>& collection_samples);
    

private:
    unsigned int m_threshold;
};


TriggerPrimitiveFinderPass1::TriggerPrimitiveFinderPass1(fhicl::ParameterSet const & p, art::ActivityRegistry&)
    : m_threshold(p.get<unsigned int>("Threshold", 10))
// Initialize member data here.
{
    std::cout << "Threshold is " << m_threshold << std::endl;
}

std::vector<TriggerPrimitiveFinderService::Hit>
TriggerPrimitiveFinderPass1::findHits(const std::vector<unsigned int>& channel_numbers, 
                                      const std::vector<std::vector<short>>& collection_samples)
{
    auto hits=std::vector<TriggerPrimitiveFinderService::Hit>();
    std::cout << "findHits called with " << collection_samples.size()
              << " channels. First chan has " << collection_samples[0].size() << " samples" << std::endl;
    // std::cout << "First few samples: ";
    // for(int i=0; i<10; ++i) std::cout << collection_samples[0][i] << " ";
    // std::cout << std::endl;

    for(size_t ich=0; ich<collection_samples.size(); ++ich){
        // std::cout << ich << std::endl;
        const std::vector<short>& waveform=collection_samples[ich];
        bool is_hit=false;
        bool was_hit=false;
        TriggerPrimitiveFinderService::Hit hit(channel_numbers[ich], 0, 0, 0);
        for(size_t isample=0; isample<waveform.size(); ++isample){
            // if(ich>11510) std::cout << isample << " " << std::flush;
            short adc=waveform[isample];
            is_hit=adc>500+(short)m_threshold;
            if(is_hit && !was_hit){
                // We just started a hit. Set the start time
                hit.startTime=isample;
                hit.charge=adc;
                hit.timeOverThreshold=1;
            }
            if(is_hit && was_hit){
                hit.charge+=adc;
                hit.timeOverThreshold++;
            }
            if(!is_hit && was_hit){
                // The hit is over. Push it to the output vector
                hits.push_back(hit);
            }
            was_hit=is_hit;
        }
        // std::cout << std::endl;
    }
    std::cout << "Returning " << hits.size() << " hits" << std::endl;
    return hits;
}

DECLARE_ART_SERVICE_INTERFACE_IMPL( TriggerPrimitiveFinderPass1, TriggerPrimitiveFinderService, LEGACY)
DEFINE_ART_SERVICE_INTERFACE_IMPL(  TriggerPrimitiveFinderPass1, TriggerPrimitiveFinderService)
